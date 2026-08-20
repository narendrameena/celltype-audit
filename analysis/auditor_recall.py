#!/usr/bin/env python3
"""How many real annotation errors does the contradiction sweep actually catch?

Precision was measurable from the start -- flags can be inspected. RECALL was not, because
it needs a set of KNOWN errors, and that only exists now that seven organs are hand-curated
from markers. An error here is a cell type whose label resolves (stage 1) to a CL term that
the marker-based gold rejects.

The decomposition matters more than the headline. The sweep tests whether the label's
lineage ANCHOR set is disjoint from the evidence's, so it can only ever see errors that
CROSS a lineage boundary. An error like "neutrophilic granulocyte" whose markers say
classical monocyte is invisible to it by construction -- both are haematopoietic -- and no
amount of tuning changes that. Reporting a single recall number would hide this.

Usage: python benchmark/auditor_recall.py
"""
import json
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")

from scoring_variants import ok                                       # noqa: E402
from cl_lineage import load, anchor_set                               # noqa: E402
from cl_resolve import resolve as resolve2                            # noqa: E402

ORGANS = ["Pancreas", "Liver", "Blood", "Bone_marrow", "Lung", "Kidney", "Heart"]
MINC, K = 500, 5


def main():
    g = load()["label"]
    rows = []
    for o in ORGANS:
        gp = os.path.join(HERE, "%s_gold.json" % o.lower())
        mp = os.path.join(RES, "heca_to_cl_%s.json" % o)
        if not (os.path.exists(gp) and os.path.exists(mp)):
            continue
        G = {k: v for k, v in json.load(open(gp)).items() if not k.startswith("_") and v}
        M = json.load(open(mp))["types"]
        ctx = {c["curie"] for v in M.values() for c in v.get("cl", [])}
        for t, gd in G.items():
            v = M.get(t)
            if not v or v["n_cells"] < MINC or len(v.get("cl", [])) < K:
                continue
            cur, how = resolve2(t, ctx, organ=o)
            if not cur:
                continue
            A = anchor_set(cur)
            GA = anchor_set(gd)
            anc = [anchor_set(c["curie"]) for c in v["cl"][:K]]
            flagged = bool(A) and any(anc) and all(B and not (A & B) for B in anc)
            err = not ok(cur, gd)
            # can a lineage test see this error at all?
            crosses = bool(A and GA and not (A & GA))
            rows.append({"organ": o, "label": t, "n_cells": v["n_cells"],
                         "label_cl": cur, "label_cl_name": g.get(cur, cur),
                         "gold_cl": gd, "gold_cl_name": g.get(gd, gd),
                         "error": err, "crosses_lineage": crosses, "flagged": flagged,
                         "markers": v["markers"][:4]})
    n = len(rows)
    err = [r for r in rows if r["error"]]
    cross = [r for r in err if r["crosses_lineage"]]
    within = [r for r in err if not r["crosses_lineage"]]
    fl = [r for r in rows if r["flagged"]]
    tp = [r for r in fl if r["error"]]

    print("AUDITOR RECALL against %d hand-curated cell types in %d organs\n" % (n, len(ORGANS)))
    print("  label errors found by curation      : %d (%.0f%% of resolved types)"
          % (len(err), 100 * len(err) / max(n, 1)))
    print("     of which CROSS a lineage boundary: %d   <- the sweep can see these" % len(cross))
    print("     of which stay WITHIN one lineage : %d   <- invisible to it by construction" % len(within))
    print("\n  flags raised                        : %d" % len(fl))
    print("  precision (flags that are real)     : %d/%d = %.0f%%"
          % (len(tp), len(fl), 100 * len(tp) / max(len(fl), 1)))
    print("  recall over ALL errors              : %d/%d = %.0f%%"
          % (len(tp), len(err), 100 * len(tp) / max(len(err), 1)))
    print("  recall over CATCHABLE errors        : %d/%d = %.0f%%"
          % (len([r for r in tp if r["crosses_lineage"]]), len(cross),
             100 * len([r for r in tp if r["crosses_lineage"]]) / max(len(cross), 1)))

    print("\n  errors the sweep CANNOT see (same lineage on both sides):")
    for r in sorted(within, key=lambda x: -x["n_cells"])[:8]:
        print("     %-8s %-28s %7d  label=%-22s markers=%s"
              % (r["organ"][:8], r["label"][:28], r["n_cells"],
                 r["label_cl_name"][:22], r["gold_cl_name"][:24]))
    json.dump(rows, open(os.path.join(RES, "auditor_recall.json"), "w"), indent=1)
    print("\nwrote results/auditor_recall.json")


if __name__ == "__main__":
    main()
