#!/usr/bin/env python3
"""Compare the cache-built reference against one refetched from the API.

The cache-built reference (build_wide_ref.py) was assembled from 3,888 WMG responses that
had been requested for OTHER purposes, each carrying only the genes of its own query. A
gene therefore has data only where some earlier query happened to ask for it, so the
matrix is sparse in a way that has nothing to do with biology.

fetch_reference.py asks the API for exactly the genes and tissues needed. This quantifies
the difference, because if the two disagree materially then results computed on the cache
need recomputing on the refetched reference -- which is the honest reason to check rather
than assume.

Usage: python benchmark/compare_reference.py [--a wide_ref] [--b wide_ref_repro]
"""
import argparse
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")


def load(stem):
    npz = np.load(os.path.join(RES, stem + ".npz"))
    idx = json.load(open(os.path.join(RES, stem + "_index.json")))
    return {k[3:]: npz[k] for k in npz.files if k.startswith("M__")}, idx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", default="wide_ref")
    ap.add_argument("--b", default="wide_ref_repro")
    x = ap.parse_args()
    A, ia = load(x.a)
    B, ib = load(x.b)
    common = sorted(set(A) & set(B))
    print("REFERENCE COMPARISON  %s (cache)  vs  %s (API)\n" % (x.a, x.b))
    print("  tissues: %d vs %d, %d in common\n" % (len(A), len(B), len(common)))
    print("  %-18s %9s %9s %9s %9s %10s" % ("tissue", "terms A", "terms B",
                                            "nonzero A", "nonzero B", "corr"))
    print("  " + "-" * 70)
    corrs, denser = [], 0
    for ub in common[:24]:
        ga, gb = ia["gene_ix"][ub], ib["gene_ix"][ub]
        ta, tb = ia["term_ix"][ub], ib["term_ix"][ub]
        genes = sorted(set(ga) & set(gb))
        terms = sorted(set(ta) & set(tb))
        if len(genes) < 5 or len(terms) < 5:
            continue
        Ma = A[ub][np.ix_([ta[t] for t in terms], [ga[g] for g in genes])]
        Mb = B[ub][np.ix_([tb[t] for t in terms], [gb[g] for g in genes])]
        nza, nzb = (Ma > 0).mean(), (Mb > 0).mean()
        both = (Ma > 0) & (Mb > 0)
        c = np.corrcoef(Ma[both], Mb[both])[0, 1] if both.sum() > 20 else float("nan")
        corrs.append(c)
        denser += nzb > nza
        print("  %-18s %9d %9d %8.0f%% %8.0f%% %10.3f"
              % (ub, len(ta), len(tb), 100 * nza, 100 * nzb, c))
    if corrs:
        print("\n  median correlation where both have data : %.3f" % np.nanmedian(corrs))
        print("  tissues where the refetched matrix is denser : %d of %d" % (denser, len(corrs)))
        print("\n  A high correlation means the cache was not WRONG where it had data; a")
        print("  large density gap means it was INCOMPLETE. Those have different remedies:")
        print("  the first would invalidate results, the second only costs sensitivity.")


if __name__ == "__main__":
    main()
