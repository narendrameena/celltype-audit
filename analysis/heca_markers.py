#!/usr/bin/env python3
"""Compute marker panels for every uHAF cell type in an hECA v2.0 organ, FROM THE DATA.

No curated marker table is trusted (uHAF's own is row-misaligned; see PUBLICATION_PLAN.md).
Markers are derived from the expression matrix itself, then used to map each cell type to the
Cell Ontology by expression via the CELLxGENE WMG API.

Per gene g and cell type t:
    pc[t,g] = fraction of t's cells with non-zero g        (detection rate)
    me[t,g] = mean expression of g among those cells
    binary_score(t,g) = pc[t,g] * (1 - max_{t'!=t} pc[t',g])    (NS-Forest-style: on in t, off elsewhere)

Reads .h5ad directly with h5py (CSC, gene-major) so scanpy/anndata are not required.

Usage:  python benchmark/heca_markers.py <organ.h5ad> [--topk 5] [--min-cells 50]
Output: results/heca_markers_<organ>.json
"""
import json
import os
import sys

import re

import h5py
import numpy as np

# Genes that are never useful cell-type markers: mitochondrial transcripts and MT/ribosomal
# pseudogenes. Low-quality clusters otherwise "mark" on these (observed for hECA PP cell/Neuron).
JUNK = re.compile(r"^(MT-|MT[A-Z]{2}\d|MTRNR|MTND|MTCO|MTATP|MTCY)|^RP[LS]\d+.*P\d+$|^RNA5|^RNU\d")
# Sex-linked genes track the DONOR, not the cell type; they surface as spurious markers
# whenever a cluster is dominated by one donor (observed: RPS4Y1 for adrenal Schwann precursors).
SEX = {"XIST", "TSIX", "RPS4Y1", "RPS4Y2", "DDX3Y", "UTY", "USP9Y", "EIF1AY", "KDM5D",
       "NLGN4Y", "ZFY", "TXLNGY", "TMSB4Y", "PRKY", "UTX", "NLGN4X"}
# Categories that are not cell types. A large catch-all ("Unclassified") is especially harmful:
# it contains cells of every type, so as a competitor it suppresses every real marker.
NOT_A_CELL_TYPE = {"unclassified", "unknown", "other", "others", "mixed", "doublet", "doublets",
                   "proliferating cell", "proliferative cell", "cycling cell", "na", "n/a",
                   "undetermined", "unassigned", "ambiguous", "low quality", "filtered"}

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)


def _strs(a):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def marker_table(path, topk=5, min_cells=50, block=2000, verbose=True):
    f = h5py.File(path, "r")
    X = f["X"]
    enc = dict(X.attrs).get("encoding-type")
    # hECA ships organs in BOTH layouts: sparse CSC (gene-major) and a plain dense array.
    if enc in ("csc_matrix", "csr_matrix"):
        n_cells, n_genes = (int(x) for x in X.attrs["shape"])
    elif isinstance(X, h5py.Dataset):
        n_cells, n_genes = int(X.shape[0]), int(X.shape[1])
    else:
        raise SystemExit("unsupported X encoding: %s" % enc)
    genes = _strs(f["var/Gene_symbol"][:])
    ctg = f["obs/cell_type"]
    cats = _strs(ctg["categories"][:])
    codes = ctg["codes"][:].astype(np.int32)
    n_t = len(cats)
    n_per_t = np.bincount(codes[codes >= 0], minlength=n_t).astype(np.float64)

    nnz_ct = np.zeros((n_genes, n_t), dtype=np.float64)   # cells of type t expressing gene g
    sum_ct = np.zeros((n_genes, n_t), dtype=np.float64)   # summed expression

    if enc == "csc_matrix":
        indptr = f["X/indptr"][:]
        data_ds, idx_ds = f["X/data"], f["X/indices"]
        for s in range(0, n_genes, block):
            e = min(s + block, n_genes)
            lo, hi = int(indptr[s]), int(indptr[e])
            if hi <= lo:
                continue
            vals = data_ds[lo:hi]
            rows = idx_ds[lo:hi]
            counts = np.diff(indptr[s:e + 1]).astype(np.int64)
            gene_local = np.repeat(np.arange(e - s, dtype=np.int64), counts)
            key = gene_local * n_t + codes[rows]
            m = (e - s) * n_t
            nnz_ct[s:e] = np.bincount(key, minlength=m)[:m].reshape(e - s, n_t)
            sum_ct[s:e] = np.bincount(key, weights=vals, minlength=m)[:m].reshape(e - s, n_t)
            if verbose and (s // block) % 5 == 0:
                print("    genes %d/%d" % (e, n_genes), flush=True)
    elif enc == "csr_matrix":
        # cell-major: indptr runs over cells, indices are gene ids
        indptr = f["X/indptr"][:]
        data_ds, idx_ds = f["X/data"], f["X/indices"]
        cblock = 20000
        m = n_genes * n_t
        for s in range(0, n_cells, cblock):
            e = min(s + cblock, n_cells)
            lo, hi = int(indptr[s]), int(indptr[e])
            if hi <= lo:
                continue
            vals = data_ds[lo:hi]
            cols = idx_ds[lo:hi].astype(np.int64)
            counts = np.diff(indptr[s:e + 1]).astype(np.int64)
            t_per_entry = np.repeat(codes[s:e].astype(np.int64), counts)
            good = t_per_entry >= 0
            key = cols[good] * n_t + t_per_entry[good]
            nnz_ct += np.bincount(key, minlength=m)[:m].reshape(n_genes, n_t)
            sum_ct += np.bincount(key, weights=vals[good], minlength=m)[:m].reshape(n_genes, n_t)
            if verbose:
                print("    cells %d/%d" % (e, n_cells), flush=True)
    else:
        cblock = max(1, int(2e8 // max(n_genes, 1)))       # ~1.6GB per read
        for s in range(0, n_cells, cblock):
            e = min(s + cblock, n_cells)
            B = X[s:e, :]
            cb = codes[s:e]
            for t in np.unique(cb):
                if t < 0:
                    continue
                sub = B[cb == t]
                nnz_ct[:, t] += (sub > 0).sum(axis=0)
                sum_ct[:, t] += sub.sum(axis=0)
            if verbose:
                print("    cells %d/%d" % (e, n_cells), flush=True)

    with np.errstate(divide="ignore", invalid="ignore"):
        pc = np.where(n_per_t > 0, nnz_ct / np.maximum(n_per_t, 1), 0.0)      # (genes, types)
        me = np.where(nnz_ct > 0, sum_ct / np.maximum(nnz_ct, 1), 0.0)

    ok_gene = np.array([(not JUNK.match(g)) and (g.upper() not in SEX) for g in genes])
    out = {}
    dropped = [str(cats[i]) for i in range(n_t)
               if str(cats[i]).strip().lower() in NOT_A_CELL_TYPE and n_per_t[i] > 0]
    keep = [i for i in range(n_t)
            if n_per_t[i] >= min_cells and str(cats[i]).strip().lower() not in NOT_A_CELL_TYPE]
    for i in keep:
        others = [j for j in keep if j != i]
        if others:
            sub = pc[:, others]
            order = np.sort(sub, axis=1)
            max_out = order[:, -1]
            # robust: ignore the single closest competitor, so a DUPLICATE category (e.g. uHAF
            # lists both "Gamma cell" and "PP cell") cannot cancel a real marker
            max_out2 = order[:, -2] if len(others) > 1 else np.zeros(n_genes)
        else:
            max_out = max_out2 = np.zeros(n_genes)
        score = pc[:, i] * (1.0 - max_out) * ok_gene
        score2 = pc[:, i] * (1.0 - max_out2) * ok_gene
        top = np.argsort(-score)[:topk]
        top2 = np.argsort(-score2)[:topk]
        rival = None
        if others and score2[top2[0]] > 3 * max(score[top[0]], 1e-9):
            j = others[int(np.argmax(pc[top2[0], others]))]
            rival = str(cats[j])          # this label's markers are masked by that label
        out[str(cats[i])] = {
            "n_cells": int(n_per_t[i]),
            "markers": [{"gene": str(genes[g]), "binary_score": round(float(score[g]), 4),
                         "pc_in": round(float(pc[g, i]), 4),
                         "pc_out_max": round(float(max_out[g]), 4),
                         "mean_expr": round(float(me[g, i]), 3)} for g in top],
            "markers_robust": [{"gene": str(genes[g]), "binary_score": round(float(score2[g]), 4),
                                "pc_in": round(float(pc[g, i]), 4)} for g in top2],
            "masked_by": rival,
        }
    f.close()
    return out, {"n_cells": n_cells, "n_genes": n_genes, "n_types": n_t,
                 "types_kept": len(keep), "min_cells": min_cells, "topk": topk,
                 "dropped_not_cell_types": dropped,
                 "dropped_cells": int(sum(n_per_t[i] for i in range(n_t)
                                          if str(cats[i]).strip().lower() in NOT_A_CELL_TYPE))}


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    path = sys.argv[1]
    topk = int(sys.argv[sys.argv.index("--topk") + 1]) if "--topk" in sys.argv else 5
    minc = int(sys.argv[sys.argv.index("--min-cells") + 1]) if "--min-cells" in sys.argv else 50
    organ = os.path.basename(path).replace("RNA-", "").replace(".h5ad", "")
    print("computing markers for %s ..." % organ, flush=True)
    tbl, meta = marker_table(path, topk=topk, min_cells=minc)
    meta["organ"] = organ
    json.dump({"meta": meta, "types": tbl},
              open(os.path.join(RES, "heca_markers_%s.json" % organ), "w"), indent=1)
    print("\n%s: %d cells, %d types (%d kept with >=%d cells)"
          % (organ, meta["n_cells"], meta["n_types"], meta["types_kept"], minc))
    if meta["dropped_not_cell_types"]:
        print("   excluded as NOT cell types: %s  (%d cells)"
              % (", ".join(meta["dropped_not_cell_types"]), meta["dropped_cells"]))
    for t, v in sorted(tbl.items(), key=lambda x: -x[1]["n_cells"])[:40]:
        tag = ("  [masked by: %s]" % v["masked_by"]) if v.get("masked_by") else ""
        print("   %-38s n=%-7d %s%s" % (t[:38], v["n_cells"],
                                        ", ".join(m["gene"] for m in v["markers"]), tag))
    dup = {t: v["masked_by"] for t, v in tbl.items() if v.get("masked_by")}
    if dup:
        print("\n  REDUNDANT/OVERLAPPING annotation categories detected from data:")
        for t, r in dup.items():
            print("     %-40s masked by  %s" % (t[:40], r))
            print("       -> robust markers: %s"
                  % ", ".join(m["gene"] for m in tbl[t]["markers_robust"]))
    print("\nwrote results/heca_markers_%s.json" % organ)
