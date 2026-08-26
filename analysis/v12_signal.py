#!/usr/bin/env python3
"""Measure the confidence signals V12 declines to threshold on, so the claim is derived.

The proposals page asserts that a pass/fail verdict cannot rest on the scorer's confidence
and quotes two AUCs for it. Those were measured once and then carried as prose, which the
page's own admissibility rule forbids: a claim counts only if it is grounded, derived or
measured. This derives them, on the same 297 hand-curated cell types the mapper scores, and
writes them where the page and the V12 docstring can quote a computed number.

Definitions, stated because the value depends on them:
  correct  a top-1 verdict that is not "wrong" under score_organ.verdict (gold term, or
           within three is_a steps of it, either direction)
  margin   top-1 score minus the score of the best candidate NOT is_a-related to top-1;
           0.0 when the shortlist holds no unrelated candidate
  depth    number of is_a levels above the predicted term, counted to the roots
  AUC      Mann-Whitney: P(signal higher on a correct call than on a wrong one), ties half
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")

import gold_organs                                          # noqa: E402
import score_organ as SO                                    # noqa: E402
from cl_lineage import ancestors, load                      # noqa: E402


def related(a, b):
    return a == b or b in ancestors(a) or a in ancestors(b)


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


def auc(rows, i):
    pos = [r[i] for r in rows if r[0]]
    neg = [r[i] for r in rows if not r[0]]
    if not pos or not neg:
        return float("nan")
    wins = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg)
    return wins / (len(pos) * len(neg))


def main():
    load()
    rows = []
    for organ in gold_organs.curated():
        gp = os.path.join(HERE, "%s_gold.json" % organ.lower())
        mp = os.path.join(RES, "heca_to_cl_%s.json" % organ)
        if not (os.path.exists(gp) and os.path.exists(mp)):
            continue
        GOLD = {k: v for k, v in json.load(open(gp)).items()
                if not k.startswith("_") and v}
        for label, v in json.load(open(mp))["types"].items():
            gold = GOLD.get(label)
            cl = v.get("cl") or []
            if not gold or not cl or v.get("n_cells", 0) < 500:
                continue
            top = cl[0]
            runner = next((c for c in cl[1:] if not related(c["curie"], top["curie"])), None)
            rows.append((SO.verdict(top["curie"], gold) != "wrong",
                         top["score"] - (runner["score"] if runner else 0.0),
                         depth(top["curie"]),
                         top["score"]))
    right = sorted(r[3] for r in rows if r[0])
    p10 = right[max(int(0.10 * len(right)) - 1, 0)] if right else 0.0
    wrong = [r for r in rows if not r[0]]
    over = sum(1 for r in wrong if r[3] >= p10) / max(len(wrong), 1)
    out = {"n": len(rows),
           "n_correct": sum(1 for r in rows if r[0]),
           "n_wrong": len(wrong),
           "auc_margin_over_unrelated_runner_up": round(auc(rows, 1), 3),
           "auc_predicted_term_depth": round(auc(rows, 2), 3),
           "wrong_calls_above_10th_pct_of_right": round(over, 3)}
    json.dump(out, open(os.path.join(RES, "v12_signal.json"), "w"), indent=1)
    print("  n = %d (%d correct, %d wrong)" % (out["n"], out["n_correct"], out["n_wrong"]))
    print("  margin over best unrelated runner-up : AUC %.3f" % out["auc_margin_over_unrelated_runner_up"])
    print("  predicted term is_a depth            : AUC %.3f" % out["auc_predicted_term_depth"])
    print("  wrong calls above the 10th percentile of right ones : %.0f%%"
          % (100 * out["wrong_calls_above_10th_pct_of_right"]))
    print("wrote results/v12_signal.json")


if __name__ == "__main__":
    main()
