#!/usr/bin/env python3
"""Score CellTypist and this method against the same gold, on the same cell types.

CellTypist is the right comparator rather than an arbitrary one: Osumi-Sutherland et al.
(Nat. Cell Biol. 23, 1129-1135, 2021) record it as an automated annotator that "maps all
cell types to CL", so it and this method answer the same question -- which CL term names
this cluster -- by different routes, label transfer against a reference versus marker
evidence against an expression reference.

The comparison is PAIRED: only cell types where both produce a call and a hand-curated
gold term exists, scored by the same rule used everywhere else (correct if the gold term
or within three is_a steps of it, either direction). Reporting either method's accuracy on
its own coverage would compare different populations.

Writes results/celltypist_head_to_head.json.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")

import score_organ as SO                                    # noqa: E402
from cl_lineage import ancestors, load                      # noqa: E402


def depth(curie):
    d, frontier, seen = 0, {curie}, set()
    while frontier and d <= 40:
        nxt = set()
        for x in frontier:
            if x in seen:
                continue
            seen.add(x)
            nxt |= set(ancestors(x)) - seen
        if not nxt:
            break
        d += 1
        frontier = nxt
    return d


def main():
    g = load()
    L = g["label"]
    CT = json.load(open(os.path.join(RES, "baseline_celltypist.json")))
    rows = []
    for organ, blob in CT.items():
        gp = os.path.join(HERE, "%s_gold.json" % organ.lower())
        mp = os.path.join(RES, "heca_to_cl_%s.json" % organ)
        if not (os.path.exists(gp) and os.path.exists(mp)):
            continue
        GOLD = {k: v for k, v in json.load(open(gp)).items()
                if not k.startswith("_") and v}
        M = json.load(open(mp))["types"]
        for label, rec in (blob.get("types") or {}).items():
            gold, mine = GOLD.get(label), M.get(label)
            if not gold or not mine or not mine.get("cl"):
                continue
            if mine.get("n_cells", 0) < 500:
                continue
            ct = (rec.get("cl") or [None])[0]
            if not ct:
                continue
            ours = mine["cl"][0]["curie"]
            rows.append({
                "organ": organ, "label": label, "n_cells": mine["n_cells"],
                "gold": gold, "gold_name": L.get(gold, gold),
                "celltypist": ct, "celltypist_name": L.get(ct, ct),
                "ours": ours, "ours_name": L.get(ours, ours),
                "celltypist_correct": SO.verdict(ct, gold) != "wrong",
                "ours_correct": SO.verdict(ours, gold) != "wrong",
                "ours_top5": min([SO.verdict(c["curie"], gold) for c in mine["cl"]],
                                 key=SO.ORDER.index) != "wrong",
            })
    n = len(rows)
    if not n:
        raise SystemExit("no paired cell types; run the mapper and the CellTypist baseline first")
    ct_ok = sum(r["celltypist_correct"] for r in rows)
    our_ok = sum(r["ours_correct"] for r in rows)
    top5 = sum(r["ours_top5"] for r in rows)
    both = sum(1 for r in rows if r["celltypist_correct"] and r["ours_correct"])
    ct_only = sum(1 for r in rows if r["celltypist_correct"] and not r["ours_correct"])
    our_only = sum(1 for r in rows if r["ours_correct"] and not r["celltypist_correct"])
    neither = n - both - ct_only - our_only
    try:
        from scipy.stats import binomtest
        p = binomtest(our_only, ct_only + our_only, 0.5).pvalue if ct_only + our_only else 1.0
    except Exception:
        p = float("nan")
    out = {"n": n, "organs": sorted({r["organ"] for r in rows}),
           "celltypist_top1": round(100.0 * ct_ok / n, 1),
           "ours_top1": round(100.0 * our_ok / n, 1),
           "ours_top5": round(100.0 * top5 / n, 1),
           "both_correct": both, "celltypist_only": ct_only, "ours_only": our_only,
           "neither": neither, "mcnemar_p": round(float(p), 3),
           "median_depth": {"celltypist": sorted(depth(r["celltypist"]) for r in rows)[n // 2],
                            "ours": sorted(depth(r["ours"]) for r in rows)[n // 2],
                            "gold": sorted(depth(r["gold"]) for r in rows)[n // 2]},
           "rows": rows}
    json.dump(out, open(os.path.join(RES, "celltypist_head_to_head.json"), "w"), indent=1)
    print("paired against the hand-curated gold, %d cell types, %d organs"
          % (n, len(out["organs"])))
    print("   CellTypist        top-1 %5.1f%%  (%d/%d)" % (out["celltypist_top1"], ct_ok, n))
    print("   celltype-audit    top-1 %5.1f%%  (%d/%d)   top-5 %5.1f%%"
          % (out["ours_top1"], our_ok, n, out["ours_top5"]))
    print("   both %d | CellTypist only %d | audit only %d | neither %d  (McNemar p = %.3f)"
          % (both, ct_only, our_only, neither, p))
    print("wrote results/celltypist_head_to_head.json")


if __name__ == "__main__":
    main()
