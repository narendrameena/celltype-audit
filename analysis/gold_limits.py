#!/usr/bin/env python3
"""Why a gold standard cannot be built automatically from expert marker databases.

Validation needs a LABEL-INDEPENDENT gold, and the obvious way to scale one is to match
each cluster's data-derived markers against expert marker sets (ASCT+B, CellMarker) and
take the best hit. This measures whether that works, against the 108 hand-curated types.

It does not. Agreement sits near 70% and does not improve when the thresholds are
tightened or when both expert resources are required to concur independently -- it gets
worse, because tightening shrinks the overlap with the hand-curated set faster than it
removes errors. The mechanism is SIZE BIAS, measured in the third output. Expert marker sets vary from 2
to over 1300 genes, and a bigger set is more likely to share three genes with anything, so
best-overlap matching drifts towards whichever term happens to be described at length. The
set chosen as "best" has a median of 29 markers against a typical candidate's 10 -- and
when the automated call disagrees with the hand-curated answer the winner is bigger still.
(An earlier guess, that expert sets for related cell types overlap heavily, was tested and
refuted: median Jaccard is 0.000 for related and unrelated pairs alike.)

Usage: python benchmark/gold_limits.py
"""
import glob
import json
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")

from scoring_variants import ok                                       # noqa: E402
from cl_lineage import load, ancestors                                # noqa: E402
from marker_gold import (asctb_by_organ, cellmarker_by_organ,          # noqa: E402
                         ORGAN2ASCTB, ORGAN2CELLMARKER)

from gold_organs import curated
GOLDO = tuple(curated())
SETTINGS = [(3, 2.0), (3, 1.5), (4, 1.5), (4, 2.0), (5, 2.0), (5, 3.0), (6, 3.0)]


def hand_curated():
    hc = {}
    for o in GOLDO:
        for k, v in json.load(open(os.path.join(HERE, "%s_gold.json" % o.lower()))).items():
            if not k.startswith("_") and v:
                hc[(o, k)] = v if isinstance(v, list) else [v]
    return hc


def best(cands, mine, minov, margin):
    sc = sorted(((len(mine & gs), c) for c, (gs, lab) in cands.items()), reverse=True)
    if not sc or sc[0][0] < minov:
        return None
    sec = sc[1][0] if len(sc) > 1 else 0
    if sec and sc[0][0] < margin * sec:
        return None
    return sc[0][1]


def sweep():
    AB, CM = asctb_by_organ(), cellmarker_by_organ()
    HC = hand_curated()
    deep = sorted(glob.glob(os.path.join(RES, "heca_markers_deep_*.json")))
    out = []
    for mode in ("union", "consensus"):
        for minov, margin in SETTINGS:
            gold, bad = {}, []
            for p in deep:
                organ = os.path.basename(p)[len("heca_markers_deep_"):-len(".json")]
                a = AB.get(ORGAN2ASCTB.get(organ, ""), {})
                c = CM.get(ORGAN2CELLMARKER.get(organ, ""), {})
                if mode == "consensus" and not (a and c):
                    continue
                pool = dict(c)
                pool.update(a)
                if not pool:
                    continue
                for t, v in json.load(open(p))["types"].items():
                    if v["n_cells"] < 500:
                        continue
                    mine = {m["gene"] for m in v["markers"]}
                    if mode == "union":
                        b = best(pool, mine, minov, margin)
                    else:
                        ba, bc = best(a, mine, minov, margin), best(c, mine, minov, margin)
                        b = ba if (ba and bc and (ok(ba, bc) or ok(bc, ba))) else None
                    if b:
                        gold[(organ, t)] = b
            both = [k for k in gold if k in HC]
            agree = 0
            for k in both:
                if any(ok(gold[k], g) for g in HC[k]):
                    agree += 1
                else:
                    bad.append({"organ": k[0], "label": k[1], "auto": gold[k],
                                "hand": HC[k][0]})
            out.append({"mode": mode, "min_overlap": minov, "margin": margin,
                        "n_gold": len(gold), "n_checked": len(both), "n_agree": agree,
                        "agreement": round(100 * agree / max(len(both), 1), 1),
                        "disagreements": bad[:12]})
    return out


def size_bias():
    """Is the winner just the biggest expert set? Sizes of pool, winners, right, wrong."""
    AB, CM = asctb_by_organ(), cellmarker_by_organ()
    HC = hand_curated()
    pool_sizes, chosen, right, wrong, best_ov = [], [], [], [], []
    for p in sorted(glob.glob(os.path.join(RES, "heca_markers_deep_*.json"))):
        organ = os.path.basename(p)[len("heca_markers_deep_"):-len(".json")]
        pool = dict(CM.get(ORGAN2CELLMARKER.get(organ, ""), {}))
        pool.update(AB.get(ORGAN2ASCTB.get(organ, ""), {}))
        if not pool:
            continue
        pool_sizes += [len(gs) for gs, _ in pool.values()]
        for t, v in json.load(open(p))["types"].items():
            if v["n_cells"] < 500:
                continue
            mine = {m["gene"] for m in v["markers"]}
            sc = sorted(((len(mine & gs), c) for c, (gs, lab) in pool.items()), reverse=True)
            if not sc or sc[0][0] < 3:
                continue
            sz = len(pool[sc[0][1]][0])
            chosen.append(sz)
            best_ov.append(sc[0][0])
            k = (organ, t)
            if k in HC:
                (right if any(ok(sc[0][1], g) for g in HC[k]) else wrong).append(sz)
    return {"pool": pool_sizes, "chosen": chosen, "right": right, "wrong": wrong,
            "best_overlap": best_ov}


if __name__ == "__main__":
    sw = sweep()
    print("AUTOMATED GOLD vs 108 HAND-CURATED TYPES\n")
    print("%-11s %6s %7s %8s %9s %11s" % ("pool", "minov", "margin", "n_gold", "checked", "agreement"))
    print("-" * 58)
    for r in sw:
        print("%-11s %6d %7.1f %8d %9d %10.0f%%"
              % (r["mode"], r["min_overlap"], r["margin"], r["n_gold"], r["n_checked"], r["agreement"]))
    sb = size_bias()
    import statistics as st
    print("\nSIZE BIAS — expert marker-set sizes")
    for k in ("pool", "chosen", "right", "wrong"):
        v = sb[k]
        print("   %-8s n=%-6d median %5.1f" % (k, len(v), st.median(v) if v else -1))
    json.dump({"sweep": sw, "size_bias": sb}, open(os.path.join(RES, "gold_limits.json"), "w"))
    print("\nwrote results/gold_limits.json")
