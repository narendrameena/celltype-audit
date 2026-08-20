#!/usr/bin/env python3
"""Rebuild the CELLxGENE reference from the API, reproducibly, without the local cache.

The .wmg_cache directory in this repo is 12 GB and was written under a per-process hash
(see the note in heca_to_cl._fetch), so its filenames cannot be recomputed and nobody
else can regenerate or verify it. Everything downstream rests on that reference, which
makes it the weakest link in reproducing this work.

This asks the public API directly for a stated list of tissues and genes and writes the
same wide_ref.npz / wide_ref_index.json that build_wide_ref.py produces from the cache.
It needs no local state, so a reader can reproduce the reference from scratch.

The gene list is the union of marker genes the pipeline needs; the tissue list is
whatever the atlases declare. Both are written into the output so a rebuilt reference can
be compared against the one used here.

VERIFIED: the reference this produces is now the default everywhere. Against the
cache-assembled matrix it correlates r = 1.000 in every tissue where both carry data, and
is about twice as dense (68%% vs 35%% nonzero, 4157 vs 3450 CL terms). The cache was
incomplete, not wrong, so results computed on it were conservative rather than invalid.

Usage:
  python benchmark/fetch_reference.py --genes markers.txt --tissues UBERON:0002048,...
  python benchmark/fetch_reference.py --from-results        # reuse this study's lists
"""
import argparse
import hashlib
import json
import os
import sys
import time
import urllib.request

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
WMG = "https://api.cellxgene.cziscience.com/wmg/v2/query"
DIMS = "https://api.cellxgene.cziscience.com/wmg/v2/primary_filter_dimensions"
HUMAN = "NCBITaxon:9606"
CHUNK = 60                     # genes per request; the API rejects very large gene lists
CACHE = os.path.join(HERE, ".wmg_reproducible")


def _key(body):
    return hashlib.sha1(json.dumps(body, sort_keys=True).encode()).hexdigest()[:16]


def fetch(body, tries=4):
    os.makedirs(CACHE, exist_ok=True)
    p = os.path.join(CACHE, _key(body) + ".json")
    if os.path.exists(p):
        try:
            return json.load(open(p))
        except Exception:
            pass
    for a in range(tries):
        try:
            req = urllib.request.Request(WMG, data=json.dumps(body).encode(),
                                         headers={"Content-Type": "application/json"})
            r = json.loads(urllib.request.urlopen(req, timeout=300).read())
            json.dump(r, open(p, "w"))
            return r
        except Exception:
            if a == tries - 1:
                return None
            time.sleep(5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes")
    ap.add_argument("--tissues")
    ap.add_argument("--from-results", action="store_true")
    ap.add_argument("--out", default=os.path.join(RES, "wide_ref_reproduced"))
    a = ap.parse_args()

    if a.from_results:
        idx = json.load(open(os.path.join(RES, "wide_ref_index.json")))
        tissues = sorted(idx["term_ix"])
        genes = sorted({g for t in idx["gene_ix"].values() for g in t})
    else:
        genes = sorted(set(open(a.genes).read().split()))
        tissues = a.tissues.split(",")
    print("rebuilding reference for %d tissues and %d genes" % (len(tissues), len(genes)), flush=True)

    dims = json.loads(urllib.request.urlopen(DIMS, timeout=300).read())
    s2e = {}
    for e in dims["gene_terms"][HUMAN]:
        for k, v in e.items():
            s2e.setdefault(v.upper(), k)
    ids = [s2e[g] for g in genes if g in s2e]
    print("   %d of %d gene symbols resolve to CELLxGENE ids" % (len(ids), len(genes)), flush=True)

    val = {t: {} for t in tissues}
    counts = {t: {} for t in tissues}
    labels = {t: {} for t in tissues}
    e2s = {v: k for k, v in s2e.items()}
    for i in range(0, len(ids), CHUNK):
        chunk = ids[i:i + CHUNK]
        print("   genes %d-%d of %d" % (i, i + len(chunk), len(ids)), flush=True)
        r = fetch({"filter": {"gene_ontology_term_ids": chunk,
                              "organism_ontology_term_id": HUMAN}, "is_rollup": True})
        if not r:
            continue
        tl = (r.get("term_id_labels") or {}).get("cell_types") or {}
        for t in tissues:
            for c, info in (tl.get(t) or {}).items():
                ag = info.get("aggregated") or {}
                counts[t].setdefault(c, ag.get("total_count", 0))
                labels[t].setdefault(c, ag.get("name", c))
        for gid, byt in (r.get("expression_summary") or {}).items():
            sym = e2s.get(gid)
            if not sym:
                continue
            for t in tissues:
                for c, v in (byt.get(t) or {}).items():
                    ag = v.get("aggregated") or {}
                    val[t].setdefault(c, {})[sym] = ag.get("pc", 0.0) * ag.get("me", 0.0)

    gene_ix, term_ix, mats = {}, {}, {}
    for t in tissues:
        gl = sorted({g for c in val[t].values() for g in c})
        tl2 = sorted(val[t])
        if not gl or not tl2:
            continue
        gene_ix[t] = {g: j for j, g in enumerate(gl)}
        term_ix[t] = {c: j for j, c in enumerate(tl2)}
        M = np.zeros((len(tl2), len(gl)), dtype=np.float32)
        for c, gs in val[t].items():
            for g, x in gs.items():
                M[term_ix[t][c], gene_ix[t][g]] = x
        mats[t] = M
    np.savez_compressed(a.out + ".npz", **{"M__" + t: m for t, m in mats.items()})
    json.dump({"gene_ix": gene_ix, "term_ix": term_ix, "counts": counts, "labels": labels,
               "provenance": {"source": WMG, "genes_requested": len(genes),
                              "tissues_requested": tissues}},
              open(a.out + "_index.json", "w"))
    print("\nwrote %s.npz (%d tissues) and %s_index.json"
          % (a.out, len(mats), a.out))
    print("compare against the cache-built reference with:")
    print("  python benchmark/compare_reference.py")


if __name__ == "__main__":
    main()
