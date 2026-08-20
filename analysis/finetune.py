#!/usr/bin/env python3
"""SingleR-style ITERATIVE FINE-TUNING on top of the subspace scorer.

We sit at ~59% top-1 but ~85% top-5: the right CL term is usually in the shortlist, just not
first. SingleR (Aran et al., Nat Immunol 2019) solves exactly this by re-ranking within the
shortlist using only the genes that separate the SHORTLIST members from each other:

    "in the fine-tuning step, SingleR reruns the correlation analysis, but only for the top cell
     types from the previous step ... only on variable genes between these cell types. The lowest
     value cell type is removed ... and then this step is repeated until only two cell types
     remain."

Here the same idea, adapted: initial ranking by cosine in the union-of-markers subspace, then
repeatedly (a) take the genes most variable ACROSS THE SURVIVING CANDIDATES in the reference,
(b) re-score the query against just those, (c) drop the weakest candidate. The elimination order
gives the final ranking.

This needs a DENSE query profile (the query's detection rate for arbitrary genes), because the
discriminating genes for a shortlist need not be among the query's own top markers. So we build
a (gene x cell type) profile straight from the h5ad, handling CSC / CSR / dense layouts.

Usage: python benchmark/finetune.py [Organ ...]   (default: the four gold organs)
"""
import glob
import json
import os
import sys
import urllib.request
from collections import defaultdict

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from scoring_variants import ok                                    # noqa: E402
from heca_markers import marker_table, NOT_A_CELL_TYPE             # noqa: E402

RES = os.path.join(HERE, "results")
TOPK_MARK = int(os.environ.get("TOPK", 20))     # markers per type -> defines the subspace
SHORTLIST = int(os.environ.get("SHORTLIST", 5))  # candidates entering fine-tuning
NVAR = int(os.environ.get("NVAR", 60))           # genes kept per fine-tuning round

ORG = {"Pancreas": ("UBERON:0001264", "../heca_data/RNA-Pancreas.h5ad"),
       "Liver": ("UBERON:0002107", "../heca_data/RNA-Liver.h5ad"),
       "Blood": ("UBERON:0000178", "../heca_data/RNA-Blood.h5ad"),
       "Bone_marrow": ("UBERON:0002371", "../heca_data/RNA-Bone_marrow.h5ad")}


def load_reference(tissues):
    """One pass over the WMG cache -> {uberon: {CL term: {gene: pc*me}}}."""
    ref = {u: defaultdict(dict) for u in tissues}
    for p in glob.glob(os.path.join(HERE, ".wmg_cache", "*.json")):
        try:
            d = json.load(open(p))
        except Exception:
            continue
        es = d.get("expression_summary")
        if not es:
            continue
        for gene, by_t in es.items():
            for u in tissues:
                ct = by_t.get(u)
                if not ct:
                    continue
                for c, v in ct.items():
                    a = v.get("aggregated") or {}
                    if a.get("pc") is not None:
                        ref[u][c][gene] = a["pc"] * a.get("me", 0.0)
    return {u: {c: v for c, v in r.items() if len(v) >= 100 and c != "CL:0000000"}
            for u, r in ref.items()}


def profile_matrix(h5, want_symbols):
    """(len(want) x n_types) detection-rate matrix, for CSC / CSR / dense h5ad alike."""
    f = h5py.File(h5, "r")
    gs = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var/Gene_symbol"][:]])
    gi = {s: i for i, s in enumerate(gs)}
    ctg = f["obs/cell_type"]
    cats = np.array([c.decode() if isinstance(c, bytes) else c for c in ctg["categories"][:]])
    codes = ctg["codes"][:].astype(np.int64)
    nt = len(cats)
    npt = np.bincount(codes[codes >= 0], minlength=nt).astype(float)
    cols = [gi.get(s, -1) for s in want_symbols]
    P = np.zeros((len(want_symbols), nt))
    X = f["X"]
    enc = dict(X.attrs).get("encoding-type")
    if enc == "csc_matrix":
        indptr = f["X/indptr"][:]
        D, I = f["X/data"], f["X/indices"]
        for k, ci in enumerate(cols):
            if ci < 0:
                continue
            lo, hi = int(indptr[ci]), int(indptr[ci + 1])
            if hi <= lo:
                continue
            rows = I[lo:hi]
            nz = np.bincount(codes[rows], weights=(D[lo:hi] > 0).astype(float), minlength=nt)
            P[k] = np.where(npt > 0, nz / np.maximum(npt, 1), 0.0)
    elif enc == "csr_matrix":
        indptr = f["X/indptr"][:]
        D, I = f["X/data"], f["X/indices"]
        keep = {ci: k for k, ci in enumerate(cols) if ci >= 0}
        n_cells = len(codes)
        for s in range(0, n_cells, 20000):
            e = min(s + 20000, n_cells)
            lo, hi = int(indptr[s]), int(indptr[e])
            if hi <= lo:
                continue
            gcol = I[lo:hi]
            vals = D[lo:hi]
            cnt = np.diff(indptr[s:e + 1]).astype(np.int64)
            tt = np.repeat(codes[s:e], cnt)
            m = np.isin(gcol, list(keep.keys())) & (vals > 0)
            if not m.any():
                continue
            kk = np.array([keep[int(c)] for c in gcol[m]])
            np.add.at(P, (kk, tt[m]), 1.0)
        P = np.where(npt > 0, P / np.maximum(npt, 1), 0.0)
    else:                                                   # dense
        n_cells = X.shape[0]
        good = [(k, ci) for k, ci in enumerate(cols) if ci >= 0]
        idx = np.array([ci for _, ci in good])
        blk = max(1, int(2e8 // max(len(idx), 1)))
        for s in range(0, n_cells, blk):
            e = min(s + blk, n_cells)
            B = X[s:e, :][:, idx]
            cb = codes[s:e]
            for t in np.unique(cb):
                if t < 0:
                    continue
                P[[k for k, _ in good], t] += (B[cb == t] > 0).sum(axis=0)
        P = np.where(npt > 0, P / np.maximum(npt, 1), 0.0)
    f.close()
    return cats, npt, P


def cosine(R, q):
    return R.dot(q) / (np.linalg.norm(R, axis=1) * np.linalg.norm(q) + 1e-9)


def finetune(R, RT, q_dense, shortlist_idx, nvar=NVAR, delta=0.05):
    """SingleR-style elimination. Returns candidates best-first."""
    cand = list(shortlist_idx)
    order = []
    while len(cand) > 1:
        sub = R[cand]                                   # candidates x genes
        spread = sub.max(axis=0) - sub.min(axis=0)      # genes that SEPARATE the survivors
        js = np.argsort(-spread)[:nvar]
        s = cosine(sub[:, js], q_dense[js])
        worst = int(np.argmin(s))
        best = s.max()
        drop = [i for i, v in enumerate(s) if v < best - delta]
        drop = drop or [worst]
        for i in sorted(drop, reverse=True):
            order.append(cand.pop(i))
        if len(cand) == 1:
            break
    return [RT[i] for i in cand] + [RT[i] for i in reversed(order)]


if __name__ == "__main__":
    organs = sys.argv[1:] or list(ORG)
    ref = load_reference({ORG[o][0] for o in organs})
    dims = json.loads(urllib.request.urlopen(
        "https://api.cellxgene.cziscience.com/wmg/v2/primary_filter_dimensions", timeout=180).read())
    E2S = {}
    for e in dims["gene_terms"]["NCBITaxon:9606"]:
        for k, v in e.items():
            E2S[k] = v.upper()

    print("%-13s %5s | %-17s | %-17s" % ("organ", "n", "SUBSPACE", "SUBSPACE+FINETUNE"))
    print("%-13s %5s | %8s %8s | %8s %8s" % ("", "", "top-1", "top-5", "top-1", "top-5"))
    print("-" * 64)
    tot = [0, 0, 0, 0, 0]
    for organ in organs:
        ub, h5 = ORG[organ]
        r = ref.get(ub) or {}
        if not r:
            continue
        GEN = sorted({g for v in r.values() for g in v})
        RT = sorted(r)
        R = np.array([[r[c].get(g, 0.0) for g in GEN] for c in RT])
        sym = [E2S.get(x, "") for x in GEN]
        keep = [i for i, s in enumerate(sym) if s]
        R, sym = R[:, keep], [sym[i] for i in keep]
        s2j = {s: j for j, s in enumerate(sym)}
        GOLD = {k: v for k, v in json.load(open(os.path.join(HERE, "%s_gold.json" % organ.lower()))).items()
                if not k.startswith("_") and v}
        PROD = json.load(open(os.path.join(RES, "heca_to_cl_%s.json" % organ)))["types"]
        tbl, _ = marker_table(h5, topk=TOPK_MARK, min_cells=50, verbose=False)
        uni = [g for g in sorted({m["gene"] for v in tbl.values() for m in v["markers"]}) if g in s2j]
        js_sub = [s2j[g] for g in uni]
        cats, npt, P = profile_matrix(h5, sym)          # dense query over the reference genes
        a, n = [0, 0, 0, 0], 0
        for t, gd in GOLD.items():
            if t not in PROD or PROD[t]["n_cells"] < 500 or t not in tbl:
                continue
            w = np.where(cats == t)[0]
            if not len(w):
                continue
            n += 1
            mine = {m["gene"]: m["pc_in"] for m in tbl[t]["markers"]}
            q_sub = np.array([mine.get(g, 0.0) for g in uni])
            s0 = cosine(R[:, js_sub], q_sub)
            rank0 = list(np.argsort(-s0))
            top0 = [RT[i] for i in rank0[:5]]
            a[0] += ok(top0[0], gd)
            a[1] += any(ok(x, gd) for x in top0)
            ft = finetune(R, RT, P[:, int(w[0])], rank0[:SHORTLIST])
            a[2] += ok(ft[0], gd)
            a[3] += any(ok(x, gd) for x in ft[:5])
        if not n:
            continue
        print("%-13s %5d | %7.1f%% %7.1f%% | %7.1f%% %7.1f%%"
              % (organ, n, 100 * a[0] / n, 100 * a[1] / n, 100 * a[2] / n, 100 * a[3] / n), flush=True)
        for i in range(4):
            tot[i] += a[i]
        tot[4] += n
    n = tot[4]
    if n:
        print("-" * 64)
        print("%-13s %5d | %7.1f%% %7.1f%% | %7.1f%% %7.1f%%"
              % ("POOLED", n, 100 * tot[0] / n, 100 * tot[1] / n, 100 * tot[2] / n, 100 * tot[3] / n))
