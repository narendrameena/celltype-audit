#!/usr/bin/env python3
"""Figure 6 — refining atlas annotation with the Cell Ontology, under gates.

a  the funnel: 602 well-powered cell types down to 101 confident refinements, and what
   each gate removes
b  the binding limit is EVIDENCE, not CL: most terms have subtypes in the ontology that
   the reference cannot support
c  validation — proposals land inside the hand-curated gold term's subtree 43 times in
   44, and the margin separates confident splits from mixed clusters

Refinements are taken from the term the LABEL resolves to (stage 1), never from the
expression top-1 used by the earlier version of this figure: top-1 is right about half
the time, so refining from it put roughly half the proposals on a base already wrong.

Outputs figure/fig6_refinement.{svg,pdf,png} + sourceData/*.tsv
"""
import glob
import json
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "figures", "_shared"))
sys.path.insert(0, os.path.join(ROOT, "cellscribe_tool", "benchmark"))

import figstyle as S                                            # noqa: E402
S.apply()
import matplotlib.pyplot as plt                                 # noqa: E402

from cl_lineage import load, ancestors                          # noqa: E402
from cl_resolve import resolve as resolve2                      # noqa: E402
from refine import children_map, descendants, score_terms, MIN_REF, MIN_MARGIN, MINC, EXACT_OK  # noqa: E402

RES = os.path.join(ROOT, "cellscribe_tool", "benchmark", "results")
G = load()
CH = children_map()
REFC = json.load(open(os.path.join(RES, "refine_reference.json")))
R = json.load(open(os.path.join(RES, "refinements.json")))
F = R["funnel"]

# ---- panel b: for every type that got past the gates, how many CL subtypes exist and
# how many of those the reference can actually support
pairs = []
for p in sorted(glob.glob(os.path.join(RES, "heca_to_cl_*.json"))):
    d = json.load(open(p))
    ctx = {c["curie"] for v in d["types"].values() for c in v.get("cl", [])}
    ru = REFC.get(d["uberon"])
    if not ru:
        continue
    for t, v in d["types"].items():
        if v["n_cells"] < MINC or not v.get("cl"):
            continue
        cur, how = resolve2(t, ctx)
        if not cur or how not in EXACT_OK:
            continue
        kids = descendants(cur, CH)
        if not kids:
            continue
        ev = [k for k in kids if ru["cnt"].get(k, [0, ""])[0] >= MIN_REF]
        pairs.append((len(kids), len(ev)))

# ---- panel c: margin, and gold-subtree compatibility where a gold standard exists
GOLDO = ["Pancreas", "Liver", "Blood", "Bone_marrow"]
gold = {}
for o in GOLDO:
    fp = os.path.join(ROOT, "cellscribe_tool", "benchmark", "%s_gold.json" % o.lower())
    gold[o] = {k: ([v] if isinstance(v, str) else list(v))
               for k, v in json.load(open(fp)).items() if not k.startswith("_") and v}
mrows, sole = [], []
for r in R["proposals"]:
    compat = None
    if r["organ"] in gold and r["label"] in gold[r["organ"]]:
        anc = ancestors(r["to_curie"])
        compat = any(gd in anc for gd in gold[r["organ"]][r["label"]])
    # a proposal with one evidenced subtype has no runner-up, so it has no margin; it is
    # counted separately rather than parked at the axis maximum as a fake spike
    (sole if r["n_evidenced_subtypes"] < 2 else mrows).append((r["margin"], compat, r["confident"], r))

fig = plt.figure(figsize=(S.W2, 3.55))
gs = fig.add_gridspec(1, 3, wspace=0.52, left=0.075, right=0.985, top=0.845, bottom=0.20,
                      width_ratios=[1.18, 0.92, 1.10])

# ------------------------------------------------------------------ a  funnel
ax = fig.add_subplot(gs[0, 0])
steps = [("well-powered cell types", F["well_powered"], None),
         ("G1  exact CL identity", F["pass_G1"], F["G1_no_exact_identity"]),
         ("G2  not contradicted", F["pass_G2"], F["G2_contradicted"]),
         ("G3  evidenced, better-scoring subtype", F["PROPOSED"],
          F["G3_no_cl_subtype"] + F["G3_no_evidenced_subtype"] + F["G3_not_discriminated"]),
         ("G4  not a mixture", F["pass_G4"], F["G4_mixture"])]
y = np.arange(len(steps))[::-1]
for (lab, val, loss), yy in zip(steps, y):
    ax.barh(yy, val, height=0.52, color=S.GREEN if lab == "not a mixture" else S.BLUE, lw=0)
    ax.text(val + 14, yy, "%d" % val, va="center", fontsize=6, color=S.INK, fontweight="bold")
    ax.text(-16, yy, lab, va="center", ha="right", fontsize=5.6, color=S.INK)
    if loss:
        ax.text(val + 100, yy, "\u2212%d" % loss, va="center", fontsize=5.5, color=S.VERM)
ax.set_xlim(0, 830)
ax.set_ylim(-0.62, len(steps) - 0.38)
ax.set_yticks([])
ax.set_xticks([0, 200, 400, 600])
ax.set_xlabel("atlas cell types")
for s in ("left",):
    ax.spines[s].set_visible(False)
ax.text(0.0, -0.235, "losses shown in orange; gates run in this order, because G2 must\n"
        "clear before a disputed label is refined", transform=ax.transAxes,
        fontsize=5.2, color=S.INK2, ha="left", va="top", linespacing=1.35)
ax.set_title("Every gate a refinement must pass", loc="left", pad=5)
S.panel(ax, "a", dx=-0.52, dy=1.14)

# ------------------------------------------------------------------ b  evidence limit
ax = fig.add_subplot(gs[0, 1])
fr = [100.0 * b / a for a, b in pairs if a > 0]
ax.hist(fr, bins=np.arange(0, 105, 5), color=S.BLUE, lw=0)
ax.set_xlabel("%% of a term's CL subtypes that the\nreference can support (>=%d cells)" % MIN_REF,
              linespacing=1.35)
ax.set_ylabel("atlas cell types")
zero = sum(1 for f in fr if f == 0)
ax.annotate("%d of %d (%.0f%%) have\nno supported subtype" % (zero, len(fr), 100.0 * zero / len(fr)),
            xy=(4, zero * 0.92), xytext=(31, zero * 1.24), fontsize=5.5, color=S.VERM,
            linespacing=1.3, arrowprops=dict(arrowstyle="-", lw=0.5, color=S.RULE))
ma, me = int(np.median([a for a, _ in pairs])), int(np.median([b for _, b in pairs]))
ax.text(0.97, 0.40, "median: %d subtypes exist,\n%d supported" % (ma, me),
        transform=ax.transAxes, ha="right", fontsize=5.4, color=S.INK2, linespacing=1.3)
ax.set_title("Evidence, not CL, is the limit", loc="left", pad=5)
S.panel(ax, "b", dx=-0.34, dy=1.14)

# ------------------------------------------------------------------ c  validation
ax = fig.add_subplot(gs[0, 2])
rng = np.random.default_rng(0)
for m, compat, conf, r in mrows:
    jit = rng.uniform(-0.30, 0.30)
    if compat is None:
        ax.scatter(m, jit, s=7, color=S.RULE, lw=0, zorder=2)
    else:
        ax.scatter(m, jit, s=17, color=S.GREEN if compat else S.VERM, lw=0, zorder=3)
ax.axvline(MIN_MARGIN, color=S.INK2, lw=0.7, ls="--", zorder=1)
nc = sum(1 for m, c, cf, r in mrows if cf) + len(sole)
allrows = mrows + sole
ok = [c for m, c, cf, r in allrows if c is not None]
ax.text(MIN_MARGIN * 1.14, 0.62, "margin gate %.2f" % MIN_MARGIN,
        fontsize=5.4, color=S.INK2, va="top")
ax.text(1.02, -0.94, "%d rejected as mixtures" % sum(1 for m, c, cf, r in mrows if not cf),
        fontsize=5.4, color=S.VERM, ha="left", va="bottom")
ax.text(0.985, 0.985, "+%d had a single evidenced\nsubtype (no runner-up)" % len(sole),
        transform=ax.transAxes, fontsize=5.2, color=S.INK2, ha="right", va="top",
        linespacing=1.3)
ax.set_xscale("log")
ax.set_xlim(0.96, 700)
ax.set_ylim(-1.02, 0.78)
ax.set_yticks([])
ax.set_xticks([1, 2, 5, 20, 100, 400])
ax.set_xticklabels(["1", "2", "5", "20", "100", "400"])
ax.set_xlabel("margin: best subtype ÷ runner-up")
for s in ("left",):
    ax.spines[s].set_visible(False)
ax.scatter([], [], s=17, color=S.GREEN, lw=0,
           label="inside the gold subtree (%d of %d)" % (sum(1 for c in ok if c), len(ok)))
ax.scatter([], [], s=17, color=S.VERM, lw=0,
           label="leaves it (%d)" % sum(1 for c in ok if not c))
ax.scatter([], [], s=7, color=S.RULE, lw=0, label="no gold standard (%d)" % (len(allrows) - len(ok)))
ax.legend(loc="upper left", bbox_to_anchor=(-0.03, -0.26), ncol=1, handletextpad=0.4,
          borderpad=0, labelspacing=0.35, fontsize=5.3)
ax.set_title("Validated against hand-curated gold", loc="left", pad=5)
S.panel(ax, "c", dx=-0.16, dy=1.14)

tsv = ["panel\tkey\tvalue"]
for lab, val, loss in steps:
    tsv.append("a\t%s\t%d" % (lab.replace(" ", "_"), val))
for k, v in sorted(F.items()):
    tsv.append("a\tfunnel_%s\t%d" % (k, v))
tsv.append("b\tmedian_subtypes_existing\t%d" % ma)
tsv.append("b\tmedian_subtypes_supported\t%d" % me)
tsv.append("b\ttypes_with_no_supported_subtype\t%d" % zero)
tsv.append("b\ttypes_assessed\t%d" % len(fr))
tsv.append("c\tproposals_total\t%d" % len(allrows))
tsv.append("c\tsole_evidenced_subtype_no_margin\t%d" % len(sole))
tsv.append("c\tconfident_margin_ge_%.2f\t%d" % (MIN_MARGIN, nc))
tsv.append("c\tgold_checked\t%d" % len(ok))
tsv.append("c\tgold_compatible\t%d" % sum(1 for c in ok if c))
for m, c, cf, r in sorted(allrows, key=lambda z: -z[0]):
    tsv.append("c\t%s|%s -> %s\tmargin=%.3f;confident=%s;gold_compatible=%s"
               % (r["organ"], r["label"], r["to_label"], r["margin"], cf, c))
out = S.save(fig, HERE, "fig6_refinement", "\n".join(tsv) + "\n")
print("funnel:", {k: F[k] for k in ("well_powered", "pass_G1", "pass_G2", "PROPOSED", "pass_G4")})
print("gold-checked %d, compatible %d" % (len(ok), sum(1 for c in ok if c)))
print("\n".join(out))
