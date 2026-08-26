#!/usr/bin/env python3
"""Run the contradiction sweep across the WHOLE atlas, and adjudicate what it raises.

The manuscript reported "14 flags from 611 resolved cell types ... 5 where the label is
wrong, 4 where our mapping is, and 5 on a boundary". That split was a hand review of an
earlier run, and nothing recomputes it, so it went stale the moment the scorer changed.
This replaces it with a derived version.

Two things make the adjudication computable now that it was not before. The gold covers
ten organs and 388 curated cell types, and for a flag falling in one of them the gold
decides which side of the disagreement is wrong:

  the gold agrees with the LABEL's term      -> our mapping is wrong, the label was right
  the gold agrees with the WINNER            -> the label is wrong
  the gold is a third term, or an abstention -> a boundary, where neither reading wins

Flags outside the curated organs stay unadjudicated and are reported as such. They are not
counted as errors, which is the mistake the original passage exists to warn against: a
method counting every flag as an atlas error overstates its yield.

Thresholds: the sweep's (margin, floor) are fitted leave-one-organ-out inside the gold. An
organ with no gold is held out by construction, so it takes the fit over all ten. Curated
organs keep their own leave-one-out setting, read from within_lineage.json, so no organ is
ever scored under a threshold fitted on itself.
"""
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")

import within_lineage as wl                                 # noqa: E402
from cl_lineage import load                                 # noqa: E402
from cl_resolve import resolve as resolve2                  # noqa: E402
from scoring_variants import ok                             # noqa: E402

SUBK = 20


def fit_over_all(rows):
    """The (margin, floor) the sweep would use for an organ outside the gold entirely."""
    floors = sorted({round(float(f), 4) for f in np.quantile(
        [r["best_score"] for r in rows], [0.0, 0.1, 0.2, 0.3, 0.4])})
    best, par = -1.0, (0.8, floors[0])
    for m in np.arange(0.5, 1.001, 0.05):
        for fl in floors:
            tp = sum(1 for r in rows if wl._call(r, m, fl) and r["error"])
            fp = sum(1 for r in rows if wl._call(r, m, fl) and not r["error"])
            fn = sum(1 for r in rows if not wl._call(r, m, fl) and r["error"])
            prec, rec = tp / max(tp + fp, 1), tp / max(tp + fn, 1)
            f1 = 2 * prec * rec / max(prec + rec, 1e-9)
            if prec >= 0.5 and f1 > best:
                best, par = f1, (float(m), float(fl))
    return par


def main():
    g = load()["label"]
    idx, mats = wl.load_ref()
    gold_rows = wl.build_rows(idx, mats)
    fitted = wl.evaluate(gold_rows)
    per_organ = {r["organ"]: (r["margin_thr"], r["floor"]) for r in fitted}
    default = fit_over_all(gold_rows)
    print("thresholds: %d curated organs keep their leave-one-out fit; every other organ "
          "uses margin %.2f, floor %.4f fitted over all ten" % (len(per_organ), *default))

    GOLDS = {}
    for o in per_organ:
        p = os.path.join(HERE, "%s_gold.json" % o.lower())
        if os.path.exists(p):
            GOLDS[o] = {k: v for k, v in json.load(open(p)).items() if not k.startswith("_")}

    resolved, flags = 0, []
    for mp in sorted(glob.glob(os.path.join(RES, "heca_to_cl_*.json"))):
        organ = os.path.basename(mp)[len("heca_to_cl_"):-len(".json")]
        dp = os.path.join(RES, "heca_markers_deep_%s.json" % organ)
        if not os.path.exists(dp):
            continue
        D = json.load(open(dp))["types"]
        M = json.load(open(mp))["types"]
        ctx = {c["curie"] for v in M.values() for c in v.get("cl", [])}
        sub = sorted({m["gene"] for v in D.values() for m in v["markers"][:SUBK]})
        margin, floor = per_organ.get(organ, default)
        for t, v in D.items():
            if v["n_cells"] < wl.MINC:
                continue
            cur, _how = resolve2(t, ctx, organ=organ)
            if not cur:
                continue
            resolved += 1
            rates = {m["gene"]: float(m.get("pc_in", 0.0)) for m in v["markers"][:SUBK]}
            f = wl.features(idx, mats, organ, [m["gene"] for m in v["markers"]], cur,
                            subspace=sub, rates=rates)
            if not f:
                continue
            pct, ratio, r, n, best, bscore, related = f
            row = {"ratio": ratio, "related_to_best": related, "best_score": float(bscore)}
            if not wl._call(row, margin, floor):
                continue
            gold = (GOLDS.get(organ) or {}).get(t)
            if organ not in GOLDS:
                verdict = "no gold for this organ"
            elif t not in (GOLDS.get(organ) or {}):
                verdict = "cluster not curated"
            elif not gold:
                verdict = "boundary (curator abstained)"
            elif ok(cur, gold):
                verdict = "our mapping is wrong"
            elif ok(best, gold):
                verdict = "the label is wrong"
            else:
                verdict = "boundary (gold is a third term)"
            flags.append({"organ": organ, "label": t, "n_cells": v["n_cells"],
                          "asserted": cur, "asserted_name": g.get(cur, cur),
                          "best": best, "best_name": g.get(best, best),
                          "ratio": round(ratio, 4), "gold": gold,
                          "gold_name": g.get(gold, gold) if gold else None,
                          "verdict": verdict})

    from collections import Counter
    c = Counter(f["verdict"] for f in flags)
    out = {"resolved_cell_types": resolved, "flags": len(flags),
           "adjudication": dict(c), "default_threshold": default, "rows": flags}
    json.dump(out, open(os.path.join(RES, "atlas_sweep.json"), "w"), indent=1)
    print("\n  %d flags from %d resolved cell types across %d organs"
          % (len(flags), resolved, len({f["organ"] for f in flags})))
    print("\n  adjudicated by the hand-curated gold:")
    for k in ("the label is wrong", "our mapping is wrong",
              "boundary (gold is a third term)", "boundary (curator abstained)"):
        if c.get(k):
            print("     %-34s %d" % (k, c[k]))
    un = c.get("no gold for this organ", 0) + c.get("cluster not curated", 0)
    print("     %-34s %d" % ("unadjudicated (no gold)", un))
    print("\n  the flags:")
    for f in sorted(flags, key=lambda x: -x["n_cells"]):
        print("     %-13s %-26s %7d  %-24s -> %-24s [%s]"
              % (f["organ"], f["label"][:26], f["n_cells"], f["asserted_name"][:24],
                 f["best_name"][:24], f["verdict"]))
    print("\nwrote results/atlas_sweep.json")


if __name__ == "__main__":
    main()
