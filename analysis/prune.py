#!/usr/bin/env python3
"""Stage 3: constrain and re-rank the shortlist.

The reasoner cannot refute a wrong genus (0 of 28, theme 04), so CL is useless as a
verifier. This asks the weaker, answerable question: can CL and a reference prior move
the correct term into first place? Top-1 is 50.9%% against the hand-curated gold while
top-5 is 80.6%%, so the right answer is usually present and merely out of first place --
a ranking problem, which constraints can fix.

Three re-rankers, measured separately because they are not equally interesting:

  P0  REFERENCE-SUPPORT PRIOR. Damp each candidate by the log of how many cells the
      CELLxGENE reference has for that CL term. THIS IS NOT AN ONTOLOGY CONSTRAINT --
      it is a frequency prior, and it is where nearly all of the gain comes from. Named
      plainly so the ontology is not credited for it.

  C1  ORGAN PRIOR. ASCT+B records which CL terms experts place in an organ, applied as
      a soft boost, not a filter: the tables carry 21-42 terms per organ against
      hundreds of candidates, so filtering on membership would discard the right answer
      more often than the wrong one. Blood has no ASCT+B table at all.

  C2  LINEAGE CONSTRAINT. The label resolves to a CL term (stage 1); candidates whose
      anchor set is DISJOINT from it are dropped. Safe on this gold set (4 wins, 0
      losses) but structurally dangerous: it forces agreement with the label, so it
      entrenches the very mis-annotations theme 07 exists to find. It must run AFTER
      the contradiction sweep, never before. It is also vacuous exactly where the
      problem is hardest -- in blood it prunes 0.5 of 20 candidates, because in a
      lineage-pure organ every candidate shares the same anchor.

PROTOCOL. Four hand-curated gold organs (the label-independent gold), 20-deep shortlists
(results/deep_to_cl_*.json). The damping strength is chosen by LEAVE-ONE-ORGAN-OUT: for
each held-out organ the parameter is fitted on the other three, so no reported number
was tuned on the data it is scored against.

Usage: python benchmark/prune.py
"""
import json
import math
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")

from cl_lineage import anchor_set                                     # noqa: E402
from cl_resolve import resolve as resolve2                            # noqa: E402
from scoring_variants import ok                                       # noqa: E402

from gold_organs import curated
# a curated organ is not enough here: prune also needs deep_to_cl_<Organ>.json, which was
# only ever generated for four organs. Intersecting with what exists keeps the module
# honest about its own coverage instead of failing on a missing input.
ORGANS = [o for o in curated()
          if os.path.exists(os.path.join(RES, "deep_to_cl_%s.json" % o))]
ASCTB_ORGAN = {"Pancreas": "pancreas", "Liver": "liver", "Bone_marrow": "bone-marrow"}
MINC = 500
GRID = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
BOOST = 1.35
LOG_LO, LOG_SPAN = 2.0, 5.2          # reference counts span ~1e2 to ~1.4e7


def asctb_sets():
    by = defaultdict(set)
    for c, v in json.load(open(os.path.join(RES, "asctb_markers.json"))).items():
        for o in v.get("organs", []):
            by[o].add(c)
    return by


def rerank(cands, boost=None, keep_anchor=None, damp=0.0, with_scores=False):
    out = []
    for c in cands:
        if keep_anchor:                                    # C2
            a = anchor_set(c["curie"])
            if a and not (a & keep_anchor):
                continue
        s = c["score"]
        if boost and c["curie"] in boost:                  # C1
            s *= BOOST
        if damp:                                           # P0
            f = (math.log10(max(c.get("n", 1), 1)) - LOG_LO) / LOG_SPAN
            s *= (1 - damp) + damp * max(0.0, min(1.0, f))
        out.append((s, c))
    out.sort(key=lambda r: -r[0])
    if with_scores:
        return out or [(c["score"], c) for c in cands]
    return [c for _, c in out] or list(cands)


def organ_counts(organ, cfg, damp, ab):
    p = os.path.join(RES, "deep_to_cl_%s.json" % organ)
    M = json.load(open(p))["types"]
    ctx = {c["curie"] for v in M.values() for c in v.get("cl", [])}
    G = {k: v for k, v in json.load(open(os.path.join(HERE, "%s_gold.json" % organ.lower()))).items()
         if not k.startswith("_") and v}
    boost = ab.get(ASCTB_ORGAN.get(organ, ""), set()) if cfg.get("c1") else None
    a = [0, 0, 0]
    for t, gd in G.items():
        v = M.get(t)
        if not v or v["n_cells"] < MINC or not v.get("cl"):
            continue
        ka = None
        if cfg.get("c2"):
            cur, _ = resolve2(t, ctx, organ=organ)
            ka = (anchor_set(cur) if cur else None) or None
        r = rerank(v["cl"], boost, ka, damp)
        a[0] += ok(r[0]["curie"], gd)
        a[1] += any(ok(c["curie"], gd) for c in r[:5])
        a[2] += 1
    return a


def held_out(cfg, ab):
    """Leave-one-organ-out: fit damp on the rest, score the held-out organ."""
    tot, picks = [0, 0, 0], []
    for held in ORGANS:
        best, bd = -1.0, 0.0
        for d in GRID:
            tr = [organ_counts(o, cfg, d, ab) for o in ORGANS if o != held]
            acc = sum(x[0] for x in tr) / max(sum(x[2] for x in tr), 1)
            if acc > best:
                best, bd = acc, d
        a = organ_counts(held, cfg, bd, ab)
        picks.append((held, bd))
        for i in range(3):
            tot[i] += a[i]
    return (100 * tot[0] / tot[2], 100 * tot[1] / tot[2], tot[2]), picks


if __name__ == "__main__":
    ab = asctb_sets()
    base = [0, 0, 0]
    for o in ORGANS:
        a = organ_counts(o, {}, 0.0, ab)
        for i in range(3):
            base[i] += a[i]
    b1, b5 = 100 * base[0] / base[2], 100 * base[1] / base[2]
    print("hand-curated gold, %d cell types, %d organs, 20-deep shortlists\n" % (base[2], len(ORGANS)))
    print("%-34s %8s %8s %10s" % ("", "top-1", "top-5", "vs base"))
    print("-" * 64)
    print("%-34s %7.1f%% %7.1f%%" % ("baseline (expression score only)", b1, b5))
    for name, cfg in (("+ P0 reference-support prior", {}),
                      ("+ P0 + C1 organ prior", {"c1": True}),
                      ("+ P0 + C1 + C2 lineage", {"c1": True, "c2": True})):
        (t1, t5, n), picks = held_out(cfg, ab)
        print("%-34s %7.1f%% %7.1f%% %9.1f" % (name, t1, t5, t1 - b1))
    print("\nall figures leave-one-organ-out; damp selected on the other three organs")
    print("of the total gain, the ontology constraints (C1+C2) contribute about 1 point;")
    print("the rest is the reference-support prior, which is not ontological.")
