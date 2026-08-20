#!/usr/bin/env python3
"""Markers for Tabula Sapiens cell types, computed exactly as for hECA.

A flag says a label's own markers point at a different lineage. The strongest check is
whether the SAME label in an INDEPENDENT atlas carries the markers it should: if Tabula
Sapiens' adipose fibroblasts look like fibroblasts and hECA's look like macrophages, the
disagreement is in hECA's annotation, not in our method.

Same statistic as heca_markers.py -- the NS-Forest-style binary score pc_in * (1 - max
pc_out) over detection rates -- so the two atlases' marker sets are directly comparable.
Tabula Sapiens is CELLxGENE-standardised, so genes live in var/feature_name and the CL
term is carried in obs/cell_type_ontology_term_id.

Usage: python benchmark/ts_markers.py [Organ ...]
"""
import glob
import json
import os
import re
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
TS = os.path.join(os.path.dirname(HERE), "..", "ts_data", "TS-%s.h5ad")
from heca_markers import JUNK, SEX, NOT_A_CELL_TYPE                   # noqa: E402

MINC = 50
TOPK = 20


def markers(path, topk=TOPK, minc=MINC, block=20000):
    f = h5py.File(path, "r")
    gv = f["var/feature_name"]
    if "categories" in gv:
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in gv["categories"][:]])
        genes = cats[gv["codes"][:]]
    else:
        genes = np.array([c.decode() if isinstance(c, bytes) else c for c in gv[:]])
    ctg = f["obs/cell_type"]
    tcats = np.array([c.decode() if isinstance(c, bytes) else c for c in ctg["categories"][:]])
    codes = ctg["codes"][:].astype(np.int64)
    clg = f["obs/cell_type_ontology_term_id"]
    ccats = np.array([c.decode() if isinstance(c, bytes) else c for c in clg["categories"][:]])
    ccodes = clg["codes"][:].astype(np.int64)
    nt = len(tcats)
    npt = np.bincount(codes[codes >= 0], minlength=nt).astype(float)

    n_cells, n_genes = f["X"].attrs["shape"] if "shape" in f["X"].attrs else (len(codes), len(genes))
    det = np.zeros((nt, int(n_genes)), dtype=np.float32)
    indptr = f["X/indptr"][:]
    D, I = f["X/data"], f["X/indices"]
    for s in range(0, int(n_cells), block):
        e = min(s + block, int(n_cells))
        lo, hi = int(indptr[s]), int(indptr[e])
        if hi <= lo:
            continue
        gcol, vals = I[lo:hi], D[lo:hi]
        cnt = np.diff(indptr[s:e + 1]).astype(np.int64)
        tt = np.repeat(codes[s:e], cnt)
        m = (vals > 0) & (tt >= 0)
        np.add.at(det, (tt[m], gcol[m]), 1.0)
    det = det / np.maximum(npt, 1)[:, None]

    keep = np.array([bool(g) and not JUNK.search(g) and g not in SEX for g in genes])
    out = {}
    for k in range(nt):
        if npt[k] < minc or tcats[k] in NOT_A_CELL_TYPE:
            continue
        inn = det[k]
        oth = np.delete(det, k, axis=0)
        score = inn * (1.0 - oth.max(axis=0))
        score = np.where(keep, score, -1.0)
        top = np.argsort(-score)[:topk]
        cl = ccats[ccodes[codes == k][0]] if (codes == k).any() else ""
        out[tcats[k]] = {"n_cells": int(npt[k]), "cl": cl,
                         "markers": [{"gene": str(genes[j]), "score": round(float(score[j]), 4),
                                      "pc_in": round(float(inn[j]), 3)} for j in top if score[j] > 0]}
    f.close()
    return out


if __name__ == "__main__":
    organs = sys.argv[1:] or [os.path.basename(p)[3:-5] for p in sorted(glob.glob(TS % "*"))]
    for o in organs:
        p = os.path.abspath(TS % o)
        out = os.path.join(RES, "ts_markers_%s.json" % o)
        if os.path.exists(out):
            print("skip %s" % o, flush=True)
            continue
        if not os.path.exists(p):
            print("missing %s" % p, flush=True)
            continue
        print("TS %s ..." % o, flush=True)
        m = markers(p)
        json.dump({"organ": o, "types": m}, open(out, "w"), indent=1)
        print("   %d types" % len(m), flush=True)
    print("TS MARKERS DONE")
