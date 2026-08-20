#!/usr/bin/env python3
"""Reference matrix: every CL term in a tissue x the marker genes we need to score.

The contradiction sweep compares lineage ANCHOR sets, so it is blind to an error that
stays inside one lineage. To see those, we need the RANK of the term the label asserts
among all candidate terms in that tissue -- which needs every term scored, not the five
that heca_to_cl.py kept.

Genes come from BOTH atlases: hECA's marker sets and Tabula Sapiens', because the same
matrix is used to audit both and a gene missing from the index cannot contribute. Built
in ONE walk over the WMG cache (its filenames were written under a per-process
hash and cannot be recomputed, so it is read content-addressed), into a float32 matrix
per tissue: terms x genes, a few MB each. Saved to results/wide_ref.npz.

Usage: python benchmark/build_wide_ref.py
"""
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
CACHE = os.path.join(HERE, ".wmg_cache")
OUT = os.path.join(RES, "wide_ref.npz")
MINC = 500


def main():
    import heca_to_cl as H
    e2s = {}
    for e in H.dims["gene_terms"][H.HUMAN]:
        for k, v in e.items():
            e2s[k] = v.upper()

    # pass A -- which tissue, and which marker genes, per organ
    organs = {}
    for p in sorted(glob.glob(os.path.join(RES, "heca_markers_deep_*.json"))):
        organ = os.path.basename(p)[len("heca_markers_deep_"):-len(".json")]
        mp = os.path.join(RES, "heca_to_cl_%s.json" % organ)
        if not os.path.exists(mp):
            continue
        ub = json.load(open(mp))["uberon"]
        genes = set()
        for t, v in json.load(open(p))["types"].items():
            if v["n_cells"] >= MINC:
                genes |= {m["gene"] for m in v["markers"]}
        if genes:
            organs[organ] = {"uberon": ub, "genes": sorted(genes)}
    # Tabula Sapiens marker genes, keyed by the tissue TS itself declares
    import h5py
    from audit_ts import ROLLUP          # same part_of rollup the audit uses
    TSDIR = os.path.join(os.path.dirname(HERE), "..", "ts_data", "TS-%s.h5ad")
    ts_genes = {}
    for p in sorted(glob.glob(os.path.join(RES, "ts_markers_*.json"))):
        o = os.path.basename(p)[len("ts_markers_"):-len(".json")]
        hp = os.path.abspath(TSDIR % o)
        if not os.path.exists(hp):
            continue
        try:
            f = h5py.File(hp, "r")
            gg = f["obs/tissue_ontology_term_id"]
            cats = [c.decode() if isinstance(c, bytes) else c for c in gg["categories"][:]]
            codes = gg["codes"][:]
            vals, cnt = np.unique(codes[codes >= 0], return_counts=True)
            ub = cats[int(vals[np.argmax(cnt)])]
            f.close()
        except Exception:
            continue
        # a TS file can hold several tissues; register its genes under EVERY tissue it
        # declares (rolled up where CELLxGENE lacks the sub-tissue), not just the
        # commonest one, or the clusters in the others have no genes to score on
        try:
            f = h5py.File(hp, "r")
            gg = f["obs/tissue_ontology_term_id"]
            allt = {ROLLUP.get(str(c.decode() if isinstance(c, bytes) else c),
                               str(c.decode() if isinstance(c, bytes) else c))
                    for c in gg["categories"][:]}
            f.close()
        except Exception:
            allt = {ub}
        gs = {m["gene"] for v in json.load(open(p))["types"].values()
              if v["n_cells"] >= MINC for m in v["markers"]}
        if gs:
            for u in allt:
                ts_genes.setdefault(u, set()).update(gs)
    print("   Tabula Sapiens adds %d genes across %d tissues"
          % (len({g for v in ts_genes.values() for g in v}), len(ts_genes)), flush=True)

    tissues = {o["uberon"] for o in organs.values()} | set(ts_genes)
    print("pass A: %d organs, %d tissues, %d distinct marker genes"
          % (len(organs), len(tissues), len({g for o in organs.values() for g in o["genes"]})), flush=True)

    # every CL term present in each tissue, from the term labels any cached response carries
    terms, counts = {}, {}
    for p in sorted(glob.glob(os.path.join(CACHE, "*.json")), key=os.path.getsize, reverse=True):
        try:
            d = json.load(open(p))
        except Exception:
            continue
        tl = (d.get("term_id_labels") or {}).get("cell_types") or {}
        for ub in tissues:
            ct = tl.get(ub)
            if ct and len(ct) > len(terms.get(ub, {})):
                terms[ub] = {c: ((v.get("aggregated") or {}).get("name", c)) for c, v in ct.items()}
                counts[ub] = {c: ((v.get("aggregated") or {}).get("total_count", 0)) for c, v in ct.items()}
        if len(terms) == len(tissues) and min(len(v) for v in terms.values()) > 50:
            break
    print("   CL terms per tissue: %s" % ", ".join("%s=%d" % (k.split(":")[-1], len(v))
                                                   for k, v in list(terms.items())[:6]), flush=True)

    gene_ix, term_ix, mats = {}, {}, {}
    for ub in tissues:
        gl = sorted({g for o in organs.values() if o["uberon"] == ub for g in o["genes"]}
                    | ts_genes.get(ub, set()))
        tl = sorted(terms.get(ub, {}))
        if not gl or not tl:
            continue
        gene_ix[ub] = {g: i for i, g in enumerate(gl)}
        term_ix[ub] = {c: i for i, c in enumerate(tl)}
        mats[ub] = np.zeros((len(tl), len(gl)), dtype=np.float32)
    print("   allocated %d tissue matrices, %.1f MB"
          % (len(mats), sum(m.nbytes for m in mats.values()) / 1e6), flush=True)

    files = sorted(glob.glob(os.path.join(CACHE, "*.json")))
    for i, fp in enumerate(files):
        if i % 500 == 0:
            print("   cache pass %d/%d" % (i, len(files)), flush=True)
        try:
            d = json.load(open(fp))
        except Exception:
            continue
        for gid, byt in (d.get("expression_summary") or {}).items():
            sym = e2s.get(gid)
            if not sym:
                continue
            for ub, M in mats.items():
                gi = gene_ix[ub].get(sym)
                if gi is None:
                    continue
                for c, v in (byt.get(ub) or {}).items():
                    ti = term_ix[ub].get(c)
                    if ti is None:
                        continue
                    a = v.get("aggregated") or {}
                    M[ti, gi] = a.get("pc", 0.0) * a.get("me", 0.0)

    np.savez_compressed(OUT, **{"M__" + ub: m for ub, m in mats.items()})
    json.dump({"organs": organs,
               "gene_ix": {k: v for k, v in gene_ix.items()},
               "term_ix": {k: v for k, v in term_ix.items()},
               "labels": terms, "counts": counts},
              open(os.path.join(RES, "wide_ref_index.json"), "w"))
    print("\nwrote results/wide_ref.npz (%.1f MB) and wide_ref_index.json"
          % (os.path.getsize(OUT) / 1e6))


if __name__ == "__main__":
    main()
