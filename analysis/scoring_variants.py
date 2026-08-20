#!/usr/bin/env python3
"""Compare scoring strategies for mapping an atlas cell type to a CL term by expression.

Baseline (what heca_to_cl.py ships) scores a CL term as the mean over 5 data-derived markers of
(pc x me). That throws away almost all the reference signal — but naively using the FULL profile
is worse at rank 1, because the profile is dominated by abundance/housekeeping rather than
identity. So the interesting axis is WEIGHTING, not gene count.

Variants
  M5     5 data-derived markers, mean(pc*me)                      [current baseline]
  M20    same but 20 markers                                       [does more help?]
  COS    cosine on z-scored full profiles                          [full profile, unweighted]
  IDF    cosine on specificity-weighted profiles                   [down-weight genes that are
         w(g) = 1/log(1 + #CL terms expressing g)                   "on" in many CL terms]
  SPEC   like IDF but weights come from the QUERY side too          [both-sided specificity]
  HIER   IDF, then walk to the most specific CL term whose score
         is within `tol` of the best                               [targets the top1/top5 gap]

Usage: python benchmark/scoring_variants.py [Organ ...]
"""
import glob
import json
import os
import sys
import urllib.request
from collections import defaultdict, deque

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(HERE, "results")
SHORT = lambda i: i.rsplit("/", 1)[-1].replace("_", ":")

g = json.load(open(os.path.join(ROOT, "cl-full.json")))["graphs"][0]
LAB = {SHORT(n["id"]): n["lbl"] for n in g["nodes"] if n.get("lbl")}
PAR, CH = defaultdict(set), defaultdict(set)
for e in g["edges"]:
    if e.get("pred") == "is_a":
        PAR[SHORT(e["sub"])].add(SHORT(e["obj"]))
        CH[SHORT(e["obj"])].add(SHORT(e["sub"]))


def within(c, adj, d=3):
    out, fr = set(), {c}
    for _ in range(d):
        nx = set()
        for x in fr:
            nx |= adj.get(x, set())
        nx -= out
        out |= nx
        fr = nx
    return out


def ok(p, gd):
    return p == gd or p in within(gd, PAR) or p in within(gd, CH)


def depth(c):
    d, fr = 0, {c}
    while fr and d < 25:
        nx = set()
        for x in fr:
            nx |= PAR.get(x, set())
        if not nx:
            break
        d += 1
        fr = nx
    return d


_dims = None


def ens2sym():
    global _dims
    if _dims is None:
        _dims = json.loads(urllib.request.urlopen(
            "https://api.cellxgene.cziscience.com/wmg/v2/primary_filter_dimensions",
            timeout=180).read())
    m = {}
    for e in _dims["gene_terms"]["NCBITaxon:9606"]:
        for k, v in e.items():
            m[k] = v.upper()
    return m


def reference(uberon, min_genes=100):
    """CL term -> {ensembl gene: pc*me} for one tissue, from every cached WMG response."""
    ref = defaultdict(dict)
    for p in glob.glob(os.path.join(HERE, ".wmg_cache", "*.json")):
        try:
            d = json.load(open(p))
        except Exception:
            continue
        es = d.get("expression_summary")
        if not es:
            continue
        for gene, by_t in es.items():
            ct = by_t.get(uberon)
            if not ct:
                continue
            for c, v in ct.items():
                a = v.get("aggregated") or {}
                if a.get("pc") is not None:
                    ref[c][gene] = a["pc"] * a.get("me", 0.0)
    return {c: v for c, v in ref.items() if len(v) >= min_genes and c != "CL:0000000"}


def query(h5, genes_sym):
    """uHAF cell type -> profile over the given gene symbols (pc*me), from the h5ad."""
    f = h5py.File(h5, "r")
    gs = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var/Gene_symbol"][:]])
    gi = {s: i for i, s in enumerate(gs)}
    ctg = f["obs/cell_type"]
    cats = np.array([c.decode() if isinstance(c, bytes) else c for c in ctg["categories"][:]])
    codes = ctg["codes"][:]
    nt = len(cats)
    npt = np.bincount(codes[codes >= 0], minlength=nt).astype(float)
    Q = np.zeros((len(genes_sym), nt))
    enc = dict(f["X"].attrs).get("encoding-type")
    if enc != "csc_matrix":
        f.close()
        return cats, npt, None                      # only CSC supported in this experiment
    indptr = f["X/indptr"][:]
    D, I = f["X/data"], f["X/indices"]
    for j, s in enumerate(genes_sym):
        ci = gi.get(s, -1)
        if ci < 0:
            continue
        lo, hi = int(indptr[ci]), int(indptr[ci + 1])
        if hi <= lo:
            continue
        rows, vals = I[lo:hi], D[lo:hi]
        tt = codes[rows]
        nz = np.bincount(tt, weights=(vals > 0).astype(float), minlength=nt)
        sm = np.bincount(tt, weights=vals, minlength=nt)
        pc = np.where(npt > 0, nz / np.maximum(npt, 1), 0.0)
        me = np.where(nz > 0, sm / np.maximum(nz, 1), 0.0)
        Q[j] = pc * me
    f.close()
    return cats, npt, Q


def _cos(R, q):
    num = R.dot(q)
    den = np.linalg.norm(R, axis=1) * np.linalg.norm(q) + 1e-9
    return num / den


def rank_variants(R, RT, q, marker_scores, tol=0.02):
    """Return {variant: ranked CL list}."""
    out = {}
    zr = (R - R.mean(0)) / (R.std(0) + 1e-9)
    zq = (q - q.mean()) / (q.std() + 1e-9)
    out["COS"] = [RT[i] for i in np.argsort(-_cos(zr, zq))]
    present = (R > 0).sum(axis=0)
    w = 1.0 / np.log1p(1 + present)
    out["IDF"] = [RT[i] for i in np.argsort(-_cos(zr * w, zq * w))]
    qw = w * (q / (q.max() + 1e-9))
    out["SPEC"] = [RT[i] for i in np.argsort(-_cos(zr * w, zq * qw / (qw.max() + 1e-9)))]
    for k, v in marker_scores.items():
        out[k] = v
    best = out["IDF"]
    s = _cos(zr * w, zq * w)
    smax = s.max()
    cand = [RT[i] for i in np.argsort(-s) if s[i] >= smax - tol]
    out["HIER"] = sorted(cand, key=lambda c: -depth(c)) + [c for c in best if c not in cand]
    return out
