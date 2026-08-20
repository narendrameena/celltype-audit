#!/usr/bin/env python3
"""Cross-atlas confirmation: ask an independent atlas what a flagged cluster looks like.

Comparing the two atlases' top marker lists directly does NOT work, and it is worth
saying why: an NS-Forest-style marker is chosen to discriminate a cell type from the
OTHER types present in that atlas, so the same cell type gets different top markers in
different atlases. hECA and Tabula Sapiens share no top markers even for macrophage.

So the query is turned around. Take the markers hECA derived for a flagged cluster and
ask which Tabula Sapiens cell type expresses them: detection rate for those genes, per TS
cell type, in the same organ. If hECA's "fibroblast" markers light up TS macrophages, the
label is wrong in hECA, and an atlas that never saw our pipeline says so.

Usage: python benchmark/ts_query.py
"""
import json
import os
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
TS = os.path.join(os.path.dirname(HERE), "..", "ts_data", "TS-%s.h5ad")


def detection(organ, genes, block=20000):
    """-> (type names, n cells, detection-rate matrix types x genes) in Tabula Sapiens."""
    p = os.path.abspath(TS % organ)
    if not os.path.exists(p):
        return None
    f = h5py.File(p, "r")
    gv = f["var/feature_name"]
    if "categories" in gv:
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in gv["categories"][:]])
        gnames = cats[gv["codes"][:]]
    else:
        gnames = np.array([c.decode() if isinstance(c, bytes) else c for c in gv[:]])
    gi = {}
    for j, g in enumerate(gnames):
        gi.setdefault(str(g).upper(), j)
    cols = [gi.get(g.upper(), -1) for g in genes]
    ctg = f["obs/cell_type"]
    tcats = np.array([c.decode() if isinstance(c, bytes) else c for c in ctg["categories"][:]])
    codes = ctg["codes"][:].astype(np.int64)
    nt = len(tcats)
    npt = np.bincount(codes[codes >= 0], minlength=nt).astype(float)
    want = {c: k for k, c in enumerate(cols) if c >= 0}
    det = np.zeros((nt, len(genes)))
    n_cells = len(codes)
    indptr = f["X/indptr"][:]
    D, I = f["X/data"], f["X/indices"]
    for s in range(0, n_cells, block):
        e = min(s + block, n_cells)
        lo, hi = int(indptr[s]), int(indptr[e])
        if hi <= lo:
            continue
        gcol, vals = I[lo:hi], D[lo:hi]
        cnt = np.diff(indptr[s:e + 1]).astype(np.int64)
        tt = np.repeat(codes[s:e], cnt)
        m = (vals > 0) & (tt >= 0) & np.isin(gcol, list(want))
        if not m.any():
            continue
        kk = np.array([want[int(c)] for c in gcol[m]])
        np.add.at(det, (tt[m], kk), 1.0)
    f.close()
    return tcats, npt, det / np.maximum(npt, 1)[:, None]


if __name__ == "__main__":
    flags = json.load(open(os.path.join(RES, "misannotation_flags.json")))
    heca, results = {}, []
    print("CROSS-ATLAS CONFIRMATION — hECA flags queried against Tabula Sapiens\n")
    done = 0
    for fl in flags:
        organ = fl["organ"]
        if not os.path.exists(os.path.abspath(TS % organ)):
            continue
        hm = json.load(open(os.path.join(RES, "heca_markers_%s.json" % organ)))["types"]
        v = hm.get(fl["label"])
        if not v:
            continue
        genes = [m["gene"] for m in v["markers"]]
        r = detection(organ, genes)
        if r is None:
            continue
        tcats, npt, det = r
        mean = det.mean(axis=1)
        ok = npt >= 50
        order = np.argsort(-np.where(ok, mean, -1))
        best = [(tcats[i], float(mean[i]), int(npt[i])) for i in order[:3]]
        done += 1
        # verdict: the winning TS type must be a real winner (its detection clearly above
        # the field) and its lineage must be disjoint from what the label asserts
        # Verdict by CONSENSUS across the TS types the markers actually light up, not by
        # a single winner: when the top hits are fibroblast and adipose MSC, they agree
        # with each other, and a margin rule would wrongly call that ambiguous.
        # resolve TS type names with the stage-1 resolver, not raw lexical lookup:
        # "macrophage" matches CL:0000235 and the Drosophila synonym on plasmatocyte, and
        # a raw lookup drops it as ambiguous exactly where it matters most.
        from cl_lineage import anchor_set
        from cl_resolve import resolve as resolve2
        A = set(fl.get("asserted_anchors", []))
        panel = []
        for name, mval, nn in best:
            cur, _how = resolve2(name, organ=organ)
            a = anchor_set(cur) if cur else set()
            if a and mval >= 0.15:
                panel.append((name, mval, a))
        tsa = sorted(set().union(*[a for _, _, a in panel])) if panel else []
        # Weigh the best OFF-lineage match against the best ON-lineage one. Requiring every
        # match to be off-lineage is too strict -- one weak on-lineage hit far down the
        # list should not veto a dominant off-lineage one.
        off = max([m for _, m, a in panel if not (A & a)], default=0.0)
        on = max([m for _, m, a in panel if (A & a)], default=0.0)
        if not panel or best[0][1] < 0.20 or not A:
            verdict = "uninformative"
        elif off >= 0.20 and off >= 1.15 * max(on, 1e-9):
            verdict = "CONFIRMED"
        elif on >= 0.20 and on >= 1.15 * max(off, 1e-9):
            verdict = "refuted"
        else:
            verdict = "uninformative"
        results_extra = {"best_off_lineage": round(off, 3), "best_on_lineage": round(on, 3)}
        results.append({"organ": organ, "label": fl["label"], "n_cells": v["n_cells"],
                        "markers": genes, "asserted_anchors": sorted(A),
                        "ts_top": [{"type": b[0], "detection": round(b[1], 3), "n": b[2]} for b in best],
                        "ts_top_anchors": sorted(tsa), "verdict": verdict,
                        **results_extra})
        print("  %s  %s  (n=%d)" % (organ, fl["label"], v["n_cells"]))
        print("     hECA markers      : %s" % ", ".join(genes))
        print("     label asserts     : %s" % "/".join(fl.get("asserted_anchors", [])))
        print("     top TS cell types : %s" % "; ".join("%s (%.2f, n=%d)" % b for b in best))
        print("     VERDICT           : %s" % results[-1]["verdict"])
        print()
    from collections import Counter
    c = Counter(r["verdict"] for r in results)
    print("queried %d flags against Tabula Sapiens" % done)
    print("   CONFIRMED %d | refuted %d | uninformative %d"
          % (c["CONFIRMED"], c["refuted"], c["uninformative"]))
    json.dump(results, open(os.path.join(RES, "cross_atlas_confirmation.json"), "w"), indent=1)
    print("wrote results/cross_atlas_confirmation.json")
