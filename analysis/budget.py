#!/usr/bin/env python3
"""What should a triage result be reported AS?

Early drafts led with recall at a fixed review budget: a curator reviewing the top 33
candidates finds 64% of the known errors. That number then fell to 50%, then 47%, without
the method changing at all. Curating three more organs added errors to the denominator
while the budget stayed pinned at 33, so recall dropped. The slide was a property of the
gold standard's size, not of the ranking -- and a headline statistic that moves when the
thing it measures does not is the wrong statistic.

This replays the ranking over the bases as they actually grew and asks which summary
survives. The queue is orderable per type: within_lineage.queue sorts on
(related_to_best, ratio), both computed inside a single cluster's own tissue, so
restricting it to an organ subset reproduces exactly the order the pipeline would have
produced on that subset. Nothing is re-ranked here -- re-implementing the ranking in an
analysis script has silently disagreed with the pipeline before.

The answer is that precision at a fixed budget is invariant and recall is not, and the
invariance is not the trivial kind: over the growth path the top 33 turns over by 13 rows
and three errors from the newly added organs displace three older ones, yet the count of
errors in the window stays at 9. Precision is also the quantity a curator actually
experiences -- how often opening a candidate is worth the time -- and it is a rate, so it
composes across bases the way recall cannot.

The three organs curated last are worth reading separately. The queue takes no fitted
parameters, so they are a prospective held-out set: the ranking was fixed before those
cell types existed in the gold.

Usage:
    python benchmark/budget.py
"""
import json
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")

BUDGET = 33
# the bases in the order they were curated, so the replay is the real history
FIRST7 = ["Blood", "Bone_marrow", "Heart", "Kidney", "Liver", "Lung", "Pancreas"]
ADDED = ["Skin", "Spleen", "Muscle"]
GROWTH = [("7 organs", FIRST7),
          ("8 (+Skin)", FIRST7 + ADDED[:1]),
          ("9 (+Spleen)", FIRST7 + ADDED[:2]),
          ("10 (+Muscle)", FIRST7 + ADDED)]


def queue():
    return json.load(open(os.path.join(RES, "within_lineage.json")))["queue"]


def stats(q, k=BUDGET):
    n = len(q)
    e = sum(1 for r in q if r["error"])
    f = sum(1 for r in q[:k] if r["error"])
    base = e / n
    return {"n": n, "errors": e, "found": f,
            "recall": f / e, "precision": f / k, "enrichment": (f / k) / base,
            "base_rate": base}


def reviews_for(q, target):
    e = sum(1 for r in q if r["error"])
    seen = 0
    for i, r in enumerate(q, 1):
        seen += bool(r["error"])
        if seen >= target * e:
            return i
    return None


def main():
    q = queue()
    print("Which summary is invariant to the size of the gold?\n")
    print("  %-13s %5s %5s | %8s %9s %10s" %
          ("base", "N", "err", "recall@33", "prec@33", "enrich@33"))
    print("  " + "-" * 60)
    series = {"recall": [], "precision": [], "enrichment": []}
    for name, orgs in GROWTH:
        s = stats([r for r in q if r["organ"] in orgs])
        for k in series:
            series[k].append(s[k])
        print("  %-13s %5d %5d | %7.0f%% %8.1f%% %9.2fx"
              % (name, s["n"], s["errors"], 100 * s["recall"],
                 100 * s["precision"], s["enrichment"]))
    print()
    for k, v in series.items():
        print("  %-11s spread %4.0f%% of mean   [%s]"
              % (k, 100 * (max(v) - min(v)) / statistics.mean(v),
                 " ".join("%.2f" % x for x in v)))

    # the invariance is by substitution, not by the window standing still
    top_now = q[:BUDGET]
    turnover = sum(1 for r in top_now if r["organ"] in ADDED)
    print("\n  %d of the top %d rows now come from the three organs added last, and their"
          % (turnover, BUDGET))
    print("  errors rank %s -- the window's composition changed; its yield did not."
          % ", ".join(str(i) for i, r in enumerate(q, 1)
                      if r["error"] and r["organ"] in ADDED and i <= BUDGET))

    held = [r for r in q if r["organ"] in ADDED]
    ranks = [i for i, r in enumerate(q, 1) if r["error"] and r["organ"] in ADDED]
    print("\n  prospective held-out (the queue takes no fitted parameters, so these three")
    print("  organs were ranked before they were curated): %d errors in %d types,"
          % (len(ranks), len(held)))
    print("  ranked %s of %d; %d in the top 20."
          % (", ".join(map(str, ranks)), len(q), sum(1 for r in ranks if r <= 20)))

    full = stats(q)
    print("\n  reported form: a curator reviewing the top %d opens %d real errors --"
          % (BUDGET, full["found"]))
    print("  %.0f%% of reviews land on one, against a %.1f%% base rate (%.1fx)."
          % (100 * full["precision"], 100 * full["base_rate"], full["enrichment"]))
    for t in (0.5, 0.8):
        r = reviews_for(q, t)
        print("  reaching %2.0f%% of all known errors takes %d reviews (%.0f%% of types)."
              % (100 * t, r, 100 * r / full["n"]))

    curve = [dict(stats(q, k), k=k) for k in range(1, len(q) + 1)]
    out = {"budget": BUDGET, "growth": [{"base": n, "organs": o,
                                         **stats([r for r in q if r["organ"] in o])}
                                        for n, o in GROWTH],
           "spread_pct_of_mean": {k: round(100 * (max(v) - min(v))
                                           / statistics.mean(v), 1)
                                  for k, v in series.items()},
           "held_out": {"organs": ADDED, "n_types": len(held), "errors": len(ranks),
                        "ranks": ranks, "in_top_20": sum(1 for r in ranks if r <= 20)},
           "top_window_turnover": turnover,
           "reviews_for_50pct": reviews_for(q, 0.5),
           "reviews_for_80pct": reviews_for(q, 0.8),
           "headline": {k: full[k] for k in
                        ("n", "errors", "found", "precision", "recall",
                         "enrichment", "base_rate")},
           "curve": curve}
    p = os.path.join(RES, "budget.json")
    json.dump(out, open(p, "w"), indent=1)
    print("\nwrote %s" % p)


if __name__ == "__main__":
    main()
