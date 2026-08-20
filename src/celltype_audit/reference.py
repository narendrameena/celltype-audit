"""The expression reference: CELLxGENE WMG, per tissue, per CL term.

Scoring a cluster needs something to score it against. This fetches that from the public
CELLxGENE API for exactly the tissues and genes required, and caches it. It deliberately
does NOT depend on any local corpus: a reference a reader cannot regenerate is the weakest
link in reproducing an audit.

The statistic is the one the WMG API exposes: detection rate x mean expression (pc * me)
per gene per CL term per tissue. A cluster scores against a term as the mean of that over
the cluster's marker genes.
"""
import hashlib
import json
import os
import time
import urllib.request

import numpy as np

from .ontology import default_cache

#: CELLxGENE's expression vocabulary is coarser than UBERON. A dataset annotated to a
#: sub-structure has no reference of its own, so it is rolled up to the part that does.
#: These are part_of relations, not guesses, and the rollup is reported in the output.
ROLLUP = {
    "UBERON:0000362": "UBERON:0002113",   # renal medulla          -> kidney
    "UBERON:0001224": "UBERON:0002113",   # renal pelvis           -> kidney
    "UBERON:0001225": "UBERON:0002113",   # cortex of kidney       -> kidney
    "UBERON:0003517": "UBERON:0002113",   # kidney blood vessel    -> kidney
    "UBERON:0001348": "UBERON:0001013",   # brown adipose          -> adipose tissue
    "UBERON:0002190": "UBERON:0001013",   # subcutaneous adipose   -> adipose tissue
    "UBERON:0000017": "UBERON:0001264",   # exocrine pancreas      -> pancreas
    "UBERON:0000016": "UBERON:0001264",   # endocrine pancreas     -> pancreas
    "UBERON:0001416": "UBERON:0002097",   # skin of abdomen        -> skin of body
    "UBERON:0001868": "UBERON:0002097",   # skin of chest          -> skin of body
}

WMG = "https://api.cellxgene.cziscience.com/wmg/v2/query"
DIMS = "https://api.cellxgene.cziscience.com/wmg/v2/primary_filter_dimensions"
HUMAN = "NCBITaxon:9606"
CHUNK = 60                    # genes per request; the API rejects very long gene lists


def _key(body):
    return hashlib.sha1(json.dumps(body, sort_keys=True).encode()).hexdigest()[:16]


def _fetch(body, cache, tries=4):
    """Cached POST. The key is a STABLE hash -- an earlier version of this code keyed on
    hash(), which Python salts per process, so every run missed the cache and refetched."""
    os.makedirs(cache, exist_ok=True)
    p = os.path.join(cache, _key(body) + ".json")
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


class Reference:
    """Per-tissue matrices of CL term x gene, plus each term's reference cell count."""

    def __init__(self, mats, gene_ix, term_ix, counts, labels):
        self.mats, self.gene_ix, self.term_ix = mats, gene_ix, term_ix
        self.counts, self.labels = counts, labels

    # ------------------------------------------------------------------ build
    @classmethod
    def fetch(cls, tissues, genes, cache=None, verbose=True):
        cache = os.path.join(cache or default_cache(), "wmg")
        dims = json.loads(urllib.request.urlopen(DIMS, timeout=300).read())
        s2e = {}
        for e in dims["gene_terms"][HUMAN]:
            for k, v in e.items():
                s2e.setdefault(v.upper(), k)
        e2s = {v: k for k, v in s2e.items()}
        ids = [s2e[g] for g in genes if g in s2e]
        val = {t: {} for t in tissues}
        counts, labels = {t: {} for t in tissues}, {t: {} for t in tissues}
        for i in range(0, len(ids), CHUNK):
            if verbose:
                print("   reference: genes %d-%d of %d" % (i, min(i + CHUNK, len(ids)), len(ids)),
                      flush=True)
            r = _fetch({"filter": {"gene_ontology_term_ids": ids[i:i + CHUNK],
                                   "organism_ontology_term_id": HUMAN},
                        "is_rollup": True}, cache)
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
        mats, gix, tix = {}, {}, {}
        for t in tissues:
            gl = sorted({g for c in val[t].values() for g in c})
            tl2 = sorted(val[t])
            if not gl or not tl2:
                continue
            gix[t] = {g: j for j, g in enumerate(gl)}
            tix[t] = {c: j for j, c in enumerate(tl2)}
            M = np.zeros((len(tl2), len(gl)), dtype=np.float32)
            for c, gs in val[t].items():
                for g, x in gs.items():
                    M[tix[t][c], gix[t][g]] = x
            mats[t] = M
        return cls(mats, gix, tix, counts, labels)

    # ------------------------------------------------------------------ io
    def save(self, stem):
        np.savez_compressed(stem + ".npz", **{"M__" + t: m for t, m in self.mats.items()})
        json.dump({"gene_ix": self.gene_ix, "term_ix": self.term_ix,
                   "counts": self.counts, "labels": self.labels},
                  open(stem + "_index.json", "w"))

    @classmethod
    def load(cls, stem):
        npz = np.load(stem + ".npz")
        idx = json.load(open(stem + "_index.json"))
        return cls({k[3:]: npz[k] for k in npz.files if k.startswith("M__")},
                   idx["gene_ix"], idx["term_ix"], idx["counts"], idx["labels"])

    # ------------------------------------------------------------------ score
    def score(self, tissue, markers, min_ref=100):
        """-> {curie: score} for every sufficiently supported term in the tissue."""
        M, gix, tix = self.mats.get(tissue), self.gene_ix.get(tissue), self.term_ix.get(tissue)
        if M is None or not gix or not tix:
            return {}
        cols = [gix[g] for g in markers if g in gix]
        if len(cols) < 3:
            return {}
        cnt = self.counts.get(tissue, {})
        inv = {v: k for k, v in tix.items()}
        s = M[:, cols].mean(axis=1)
        return {inv[i]: float(s[i]) for i in range(len(s))
                if cnt.get(inv[i], 0) >= min_ref}

    def support(self, tissue, curie):
        return self.counts.get(tissue, {}).get(curie, 0)

    def tissues(self):
        return sorted(self.mats)
