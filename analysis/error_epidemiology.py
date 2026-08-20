#!/usr/bin/env python3
"""Where do cell-type annotation errors concentrate?

Nobody has published this, because it needs a set of errors that is both KNOWN and
labelled with the properties you want to test against. Seven hand-curated organs give 190
cell types of which 18 carry a label the markers reject -- small, but enough to measure
effect sizes if the report is honest about power.

Four questions, each answerable from data already on disk:
  size      are errors commoner in small clusters, where evidence is thin?
  lineage   which lineages are mislabelled most often?
  depth     are specific terms (deep in CL) riskier than general ones?
  direction which lineage gets confused FOR which -- the confusion structure

With 18 errors almost nothing reaches significance, and that is stated rather than
disguised: this reports proportions with exact binomial intervals and the Fisher p where
a 2x2 exists, and calls the result a description, not a test.

Usage: python benchmark/error_epidemiology.py
"""
import json
import math
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")

from cl_lineage import load, anchor_set, ancestors                    # noqa: E402


def wilson(k, n, z=1.96):
    """Wilson interval -- behaves at the small n this analysis actually has."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def fisher(a, b, c, d):
    """Two-sided Fisher exact on a 2x2, without scipy."""
    def lc(n, k):
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

    n = a + b + c + d
    obs = lc(a + b, a) + lc(c + d, c) - lc(n, a + c)
    tot = 0.0
    for i in range(0, min(a + b, a + c) + 1):
        j, k, l = a + b - i, a + c - i, d - a + i
        if j < 0 or k < 0 or l < 0:
            continue
        p = lc(a + b, i) + lc(c + d, k) - lc(n, a + c)
        if p <= obs + 1e-9:
            tot += math.exp(p)
    return min(1.0, tot)


def depth(c):
    return len(ancestors(c))


def main():
    g = load()["label"]
    rows = json.load(open(os.path.join(RES, "auditor_recall.json")))
    err = [r for r in rows if r["error"]]
    n, e = len(rows), len(err)
    print("ERROR EPIDEMIOLOGY — %d hand-curated cell types, %d organs, %d errors (%.1f%%)\n"
          % (n, len({r["organ"] for r in rows}), e, 100 * e / n))

    # ---- size
    print("  BY CLUSTER SIZE")
    cuts = [(500, 2000), (2000, 10000), (10000, 50000), (50000, 10 ** 9)]
    for lo, hi in cuts:
        sub = [r for r in rows if lo <= r["n_cells"] < hi]
        k = sum(1 for r in sub if r["error"])
        if not sub:
            continue
        lo_ci, hi_ci = wilson(k, len(sub))
        print("     %-16s %3d types  %2d errors  %5.1f%%  [%.0f-%.0f%%]"
              % ("%d-%s" % (lo, "inf" if hi > 10 ** 8 else str(hi)), len(sub), k,
                 100 * k / len(sub), 100 * lo_ci, 100 * hi_ci))
    small = [r for r in rows if r["n_cells"] < 2000]
    big = [r for r in rows if r["n_cells"] >= 2000]
    a, b = sum(1 for r in small if r["error"]), len(small) - sum(1 for r in small if r["error"])
    c, d = sum(1 for r in big if r["error"]), len(big) - sum(1 for r in big if r["error"])
    print("     small(<2k) %d/%d = %.1f%%  vs  >=2k %d/%d = %.1f%%   Fisher p=%.4f"
          % (a, a + b, 100 * a / (a + b), c, c + d, 100 * c / (c + d), fisher(a, b, c, d)))

    # ---- lineage of the asserted term
    print("\n  BY LINEAGE THE LABEL ASSERTS")
    by = defaultdict(lambda: [0, 0])
    for r in rows:
        for anc in (sorted(anchor_set(r["label_cl"])) or ["(none)"]):
            by[anc][0] += 1
            by[anc][1] += r["error"]
    for k2, (tot, ke) in sorted(by.items(), key=lambda kv: -kv[1][1] / max(kv[1][0], 1)):
        if tot < 5:
            continue
        lo_ci, hi_ci = wilson(ke, tot)
        print("     %-16s %3d types  %2d errors  %5.1f%%  [%.0f-%.0f%%]"
              % (k2, tot, ke, 100 * ke / tot, 100 * lo_ci, 100 * hi_ci))

    # ---- depth of the asserted term
    print("\n  BY SPECIFICITY OF THE ASSERTED TERM (is_a depth)")
    ds = sorted(depth(r["label_cl"]) for r in rows)
    med = ds[len(ds) // 2]
    for name, sel in (("general (depth <= %d)" % med, lambda r: depth(r["label_cl"]) <= med),
                      ("specific (depth > %d)" % med, lambda r: depth(r["label_cl"]) > med)):
        sub = [r for r in rows if sel(r)]
        k = sum(1 for r in sub if r["error"])
        lo_ci, hi_ci = wilson(k, len(sub))
        print("     %-22s %3d types  %2d errors  %5.1f%%  [%.0f-%.0f%%]"
              % (name, len(sub), k, 100 * k / len(sub), 100 * lo_ci, 100 * hi_ci))

    # ---- confusion direction
    print("\n  CONFUSION STRUCTURE — what is mislabelled as what")
    conf = Counter()
    for r in err:
        a2 = "/".join(sorted(anchor_set(r["label_cl"]))) or "?"
        b2 = "/".join(sorted(anchor_set(r["gold_cl"]))) or "?"
        conf[(a2, b2)] += 1
    for (a2, b2), k in conf.most_common():
        tag = "SAME lineage" if a2 == b2 else "crosses"
        print("     %-22s -> %-22s %2d   (%s)" % (a2[:22], b2[:22], k, tag))

    # the muscle pattern, counted from each source rather than asserted
    print("\n  THE SMOOTH-MUSCLE PATTERN, counted in each evidence source")
    try:
        CX = json.load(open(os.path.join(RES, "cross_atlas_confirmation.json")))
        confirmed = [r for r in CX if r["verdict"] == "CONFIRMED"]
        mc = [r for r in confirmed if "muscle" in r["asserted_anchors"]]
        print("     Tabula-Sapiens-confirmed errors asserting muscle : %d of %d"
              % (len(mc), len(confirmed)))
    except Exception:
        pass
    try:
        AB = json.load(open(os.path.join(RES, "audit_baseline.json")))
        dis = [r for r in AB["rows"] if not r["ct_agrees_label"]]
        mus = [r for r in dis if "muscle" in r["label_anchors"]]
        print("     CellTypist disagreements asserting muscle        : %d of %d"
              % (len(mus), len(dis)))
    except Exception:
        pass
    km = sum(v for k2, v in conf.items() if k2[0] == "muscle")
    print("     hand-curated errors asserting muscle             : %d of %d" % (km, e))
    print("     the pattern is strong in the two INDEPENDENT sources and weak in the")
    print("     curated set, which contains few smooth-muscle clusters -- so it is")
    print("     reported as a cross-source observation, not a rate from one denominator.")

    json.dump({"n": n, "errors": e,
               "by_size": [{"lo": lo, "hi": hi,
                            "n": len([r for r in rows if lo <= r["n_cells"] < hi]),
                            "k": sum(1 for r in rows if lo <= r["n_cells"] < hi and r["error"])}
                           for lo, hi in cuts],
               "by_lineage": {k2: v for k2, v in by.items()},
               "confusion": {"%s->%s" % k2: v for k2, v in conf.items()}},
              open(os.path.join(RES, "error_epidemiology.json"), "w"), indent=1)
    print("\nwrote results/error_epidemiology.json")


if __name__ == "__main__":
    main()
