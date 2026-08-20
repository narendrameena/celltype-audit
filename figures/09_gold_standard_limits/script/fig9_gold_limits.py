#!/usr/bin/env python3
"""Figure 9 — a gold standard cannot be built automatically from marker databases.

Every claim in this work rests on a LABEL-INDEPENDENT gold, and the hand-curated one is
108 cell types across 4 organs. The obvious way to scale it is to match each cluster's
data-derived markers against expert marker sets (ASCT+B, CellMarker 2.0) and take the best
hit. This figure is the record that the obvious way does not work, so that nobody -- us
included -- reports numbers against a gold built that way.

a  agreement with the hand-curated types never approaches usable, at any threshold, and
   tightening makes it worse rather than better
b  requiring two independent expert resources to concur does not rescue it either; it
   just shrinks the set
c  the mechanism is SIZE BIAS: a longer expert marker list is more likely to share three
   genes with anything, so best-overlap matching drifts to whichever cell type happens to
   be described at length -- and those are the calls that turn out wrong

Outputs figure/fig9_gold_limits.{svg,pdf,png} + sourceData/*.tsv
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "figures", "_shared"))
import figstyle as S                                            # noqa: E402
S.apply()
import matplotlib.pyplot as plt                                 # noqa: E402

RES = os.path.join(ROOT, "cellscribe_tool", "benchmark", "results")
D = json.load(open(os.path.join(RES, "gold_limits.json")))
SW, SB = D["sweep"], D["size_bias"]

fig = plt.figure(figsize=(S.W2, 3.15))
gs = fig.add_gridspec(1, 3, width_ratios=[1.12, 1.0, 1.0], wspace=0.46,
                      left=0.075, right=0.985, top=0.83, bottom=0.255)

# ------------------------------------------------------------------ a agreement
ax = fig.add_subplot(gs[0, 0])
S.panel(ax, "a", dx=-0.26, dy=1.15)
ax.set_title("Automated gold never becomes usable", loc="left", pad=5)
un = [r for r in SW if r["mode"] == "union" and r["n_checked"] >= 5]
co = [r for r in SW if r["mode"] == "consensus" and r["n_checked"] >= 5]
for rows, col, lab in ((un, S.BLUE, "single pool"), (co, S.VERM, "two-resource consensus")):
    x = [r["n_gold"] for r in rows]
    y = [r["agreement"] for r in rows]
    ax.plot(x, y, "-o", ms=3.4, lw=0.8, color=col, label=lab, mew=0)
ax.axhline(100, color=S.GREY, lw=0.7, ls=":")
ax.text(160, 100, "hand-curated\nreference", fontsize=5.3, color=S.GREY, va="center",
        ha="right", linespacing=1.25)
ax.set_xlabel("cell types the automated gold covers")
ax.set_ylabel("agreement with hand-curated (%)")
ax.set_ylim(0, 116)
ax.set_xlim(0, 175)
ax.legend(loc="lower right", fontsize=5.4, handlelength=1.3, borderpad=0.1)
ax.text(0.03, 0.20, "tightening the thresholds\nmoves LEFT, not up",
        transform=ax.transAxes, fontsize=5.3, color=S.INK2, linespacing=1.3)

# ------------------------------------------------------------------ b coverage cost
ax = fig.add_subplot(gs[0, 1])
S.panel(ax, "b", dx=-0.30, dy=1.15)
ax.set_title("What each setting actually buys", loc="left", pad=5)
# only settings with enough hand-curated overlap to give a meaningful percentage;
# the tightest consensus settings check 1-2 types and would print a hollow "100%"
allr = [r for r in SW if r["n_checked"] >= 5]
lab = ["%d / %.1f" % (r["min_overlap"], r["margin"]) for r in allr]
x = np.arange(len(allr))
ax.bar(x, [r["n_gold"] for r in allr], width=0.62,
       color=[S.BLUE if r["mode"] == "union" else S.VERM for r in allr], lw=0)
for i, r in enumerate(allr):
    ax.annotate("%d%%" % r["agreement"], (i, r["n_gold"]), xytext=(0, 4), fontsize=5.2,
                textcoords="offset points", ha="center", color=S.INK2)
ax.set_xticks(x)
ax.set_xticklabels(lab, fontsize=5.0, rotation=90)
ax.set_ylabel("cell types covered")
ax.set_ylim(0, max(r["n_gold"] for r in allr) * 1.30)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor=S.BLUE, label="single pool"),
                   Patch(facecolor=S.VERM, label="two-resource consensus")],
          loc="upper right", fontsize=5.2, handlelength=1.0, borderpad=0.1,
          labelspacing=0.3)
ax.text(0.5, -0.40, "min shared markers / margin;  label above bar = agreement\n"
        "(settings checking fewer than 5 hand-curated types are omitted)",
        transform=ax.transAxes, fontsize=5.2, color=S.INK2, ha="center", linespacing=1.35)

# ------------------------------------------------------------------ c size bias
ax = fig.add_subplot(gs[0, 2])
S.panel(ax, "c", dx=-0.30, dy=1.15)
ax.set_title("The winner is just the longest list", loc="left", pad=5)
groups = [("any candidate", SB["pool"], S.GREY),
          ("chosen as best", SB["chosen"], S.BLUE),
          ("chosen, and right", SB["right"], S.GREEN),
          ("chosen, and wrong", SB["wrong"], S.VERM)]
pos = np.arange(len(groups))
for i, (name, vals, col) in enumerate(groups):
    v = np.array(vals, dtype=float)
    v = v[v > 0]
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    ax.plot([q1, q3], [i, i], lw=2.2, color=col, solid_capstyle="round", alpha=0.5)
    ax.scatter([med], [i], s=22, color=col, zorder=3, lw=0)
    ax.annotate("%d" % med, (med, i), xytext=(0, 6), textcoords="offset points",
                fontsize=5.4, ha="center", color=col, fontweight="bold")
ax.set_yticks(pos)
ax.set_yticklabels(["%s\n(n=%d)" % (n, len([x for x in v if x > 0])) for n, v, _ in groups],
                   fontsize=5.3, linespacing=1.25)
ax.set_xscale("log")
ax.set_xlim(1.5, 700)
ax.set_ylim(-0.62, len(groups) - 0.30)
ax.set_xlabel("genes in the expert marker set (log)")
ax.invert_yaxis()
for s in ("left",):
    ax.spines[s].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.text(0.5, -0.30, "bars span the interquartile range; point is the median",
        transform=ax.transAxes, fontsize=5.2, color=S.INK2, ha="center")

tsv = ["panel\tkey\tvalue"]
for r in SW:
    tsv.append("a/b\t%s|minov=%d|margin=%.1f\tn_gold=%d;checked=%d;agree=%d;agreement=%.1f"
               % (r["mode"], r["min_overlap"], r["margin"], r["n_gold"], r["n_checked"],
                  r["n_agree"], r["agreement"]))
import statistics as st
for name, vals, _ in groups:
    v = [x for x in vals if x > 0]
    tsv.append("c\t%s\tn=%d;median=%.1f;q1=%.1f;q3=%.1f"
               % (name.replace(" ", "_"), len(v), st.median(v),
                  np.percentile(v, 25), np.percentile(v, 75)))
for r in SW:
    for d in r["disagreements"][:4]:
        tsv.append("a\tdisagreement|%s|%s\tauto=%s;hand=%s" % (d["organ"], d["label"], d["auto"], d["hand"]))
out = S.save(fig, HERE, "fig9_gold_limits", "\n".join(tsv) + "\n")
print("best agreement anywhere: %.0f%%" % max(r["agreement"] for r in SW if r["n_checked"] >= 5))
print("\n".join(out))
