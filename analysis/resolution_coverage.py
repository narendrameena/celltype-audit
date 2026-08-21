#!/usr/bin/env python3
"""How much of the atlas does lexical resolution actually reach?

The manuscript quotes a label-coverage and a cell-coverage figure for the resolver. Both
were previously typed into a docstring and went stale: the label figure was 80.4% and had
not been updated after the alias table and organ-qualified lookup were added, which between
them resolve 133 of the 602 well-powered labels. Computing it here keeps the number honest.

Usage: python benchmark/resolution_coverage.py
"""
import collections
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")

from cl_resolve import resolve as resolve2                            # noqa: E402

MINC = 500


def main():
    how = collections.Counter()
    n_t = n_r = n_c = n_rc = n_all = n_rall = 0
    per = []
    for fp in sorted(glob.glob(os.path.join(RES, "heca_to_cl_*.json"))):
        d = json.load(open(fp))
        organ = d["organ"]
        ctx = {c["curie"] for v in d["types"].values() for c in v.get("cl", [])}
        t = r = c = rc = 0
        for label, v in d["types"].items():
            cells = v.get("n_cells", 0)
            cur, status = resolve2(label, ctx, organ=organ)
            n_all += 1
            n_rall += bool(cur)
            if cells < MINC:
                continue
            t += 1
            c += cells
            if cur:
                r += 1
                rc += cells
                how[str(status)] += 1
        per.append((organ, t, r, c, rc))
        n_t += t; n_r += r; n_c += c; n_rc += rc
    gen = sum(v for k, v in how.items() if k.startswith("general"))
    print("resolution coverage, %d organs, size floor %d cells\n" % (len(per), MINC))
    print("  well-powered labels resolving      : %d of %d = %.1f%%" % (n_r, n_t, 100 * n_r / n_t))
    print("    of which an EXACT CL identity    : %d of %d = %.1f%%"
          % (n_r - gen, n_t, 100 * (n_r - gen) / n_t))
    print("    generalised (superclass only)    : %d" % gen)
    print("  cells under a resolved label       : %s of %s = %.2f%%"
          % (format(n_rc, ","), format(n_c, ","), 100 * n_rc / n_c))
    print("  every label, ignoring the floor    : %d of %d = %.1f%%"
          % (n_rall, n_all, 100 * n_rall / n_all))
    print("\n  by resolution route: %s"
          % ", ".join("%s %d" % (k, v) for k, v in how.most_common()))
    json.dump({"organs": len(per), "min_cells": MINC,
               "labels_total": n_t, "labels_resolved": n_r, "labels_exact": n_r - gen,
               "labels_generalised": gen, "cells_total": n_c, "cells_resolved": n_rc,
               "all_labels_total": n_all, "all_labels_resolved": n_rall,
               "by_route": dict(how)},
              open(os.path.join(RES, "resolution_coverage.json"), "w"), indent=1)
    print("\nwrote results/resolution_coverage.json")


if __name__ == "__main__":
    main()
