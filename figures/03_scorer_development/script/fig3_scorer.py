#!/usr/bin/env python3
"""Figure 3 — developing the scorer: four ideas failed, one worked.

Form follows the data's job, so nothing has to be decoded twice:
  a  comparison against a BASELINE   -> delta lollipop; 0 = production scorer.
                                        Filled dot top-1, open dot top-5.
  b  a trend that is not there       -> plain line
  c  PAIRED CHANGE per organ         -> dumbbell; the arrow IS the improvement

Outputs fig3_scorer.{svg,pdf,png} + fig3_scorer_source_data.tsv
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "_shared")))
import figstyle as S                                            # noqa: E402
S.apply()
import matplotlib.pyplot as plt                                 # noqa: E402

SHOW_TOP5 = os.environ.get("TOP5", "1") != "0"
NAME = "fig3_scorer" if SHOW_TOP5 else "fig3_scorer_top1only"

RES = os.path.abspath(os.path.join(HERE, "..", "..", "..",
                                   "cellscribe_tool", "benchmark", "results"))

BASE1, BASE5 = 81.8, 95.5                     # production scorer, pancreas
VARIANTS = [("full profile cosine", 68.2, 86.4),
            ("specificity weighted (IDF)", 59.1, 86.4),
            ("two-sided weighting", 68.2, 90.9),
            ("most-specific above threshold", 59.1, 86.4),
            ("discriminative SUBSPACE", 86.4, 95.5)]
SWEEP = [(5, 63.6, 86.4), (10, 59.1, 86.4), (20, 59.1, 81.8), (50, 68.2, 90.9)]

sub = json.load(open(os.path.join(RES, "subspace_test.json")))["rows"]
ORD = [o for o in ["Pancreas", "Liver", "Blood", "Bone_marrow", "POOLED"] if o in sub]

fig = plt.figure(figsize=(S.W2, 3.3))
gs = fig.add_gridspec(1, 3, wspace=0.95, left=0.20, right=0.92, top=0.78, bottom=0.22,
                      width_ratios=[1.15, 0.82, 1.05])

# ---------------------------------------------------------------- a  delta lollipop
ax = fig.add_subplot(gs[0, 0])
v = sorted(VARIANTS, key=lambda r: r[1] - BASE1)
y = np.arange(len(v))
ax.axvline(0, color=S.INK2, lw=0.8, zorder=0)
for i, r in enumerate(v):
    d1, d5 = r[1] - BASE1, r[2] - BASE5
    c = S.GREEN if d1 > 0 else S.GREY
    ax.plot([0, d1], [i, i], color=c, lw=1.5, solid_capstyle="round", zorder=1)
    ax.plot([d1], [i], "o", ms=5.5, color=c, zorder=3)
    if SHOW_TOP5:
        ax.plot([d5], [i], "o", ms=4.2, mfc="white", mec=S.INK2, mew=0.9, zorder=2)
    ax.annotate("%+.1f" % d1, xy=(d1, i), xytext=(0, 6), textcoords="offset points",
                ha="center", fontsize=5.8, color=c)
ax.set_yticks(y)
ax.set_yticklabels([r[0] for r in v], fontsize=5.8)
ax.set_ylim(-1.1, len(v) - 0.3)
ax.set_xlim(-30, 15)
ax.set_xlabel("change vs production scorer (percentage points)\npancreas, n=22 gold cell types")
if SHOW_TOP5:
    ax.plot([], [], "o", ms=5.5, color=S.GREY, label="top-1")
    ax.plot([], [], "o", ms=4.2, mfc="white", mec=S.INK2, mew=0.9, label="top-5")
    ax.legend(loc="lower left", handlelength=0.8, fontsize=5.4, borderpad=0.2,
              bbox_to_anchor=(0.0, -0.03))
ax.set_title("Four ideas failed, one worked", loc="left", pad=4)
S.panel(ax, "a", dx=-0.66, dy=1.13)

# ---------------------------------------------------------------- b  no-trend line
ax = fig.add_subplot(gs[0, 1])
ax.plot([s[0] for s in SWEEP], [s[2] for s in SWEEP], marker="o", ms=4.0, lw=1.2,
        mfc="white", mec=S.INK2, color=S.INK2, label="top-5")
ax.plot([s[0] for s in SWEEP], [s[1] for s in SWEEP], marker="o", ms=4.5, lw=1.3,
        color=S.VERM, label="top-1")
for k, v1, v5 in SWEEP:
    ax.annotate("%.1f" % v1, xy=(k, v1), xytext=(0, -9), textcoords="offset points",
                ha="center", fontsize=5.4, color=S.VERM)
    ax.annotate("%.1f" % v5, xy=(k, v5), xytext=(0, 6), textcoords="offset points",
                ha="center", fontsize=5.4, color=S.INK2)
ax.legend(loc="lower left", handlelength=1.0, fontsize=5.4, borderpad=0.2,
          bbox_to_anchor=(0.0, -0.03))
ax.set_xscale("log")
ax.set_xticks([s[0] for s in SWEEP])
ax.set_xticklabels([str(s[0]) for s in SWEEP])
ax.set_xlabel("markers per cell type")
ax.set_ylabel("accuracy (%)")
ax.set_ylim(45, 98)
ax.set_title("Marker count is not the driver", loc="left", pad=4)
S.panel(ax, "b", dx=-0.40, dy=1.13)

# ------------------------------------------------- c  the shortlist gap, and it closing
# Each organ gets two spans running from top-1 (left dot) to top-5 (right dot): the span IS the
# shortlist gap. The subspace scorer barely moves top-5 but lifts top-1, so the span SHRINKS —
# it converts shortlist hits into rank-1 hits, which is the actual claim.
ax = fig.add_subplot(gs[0, 2])
y = np.arange(len(ORD))[::-1].astype(float)
off = 0.17
for i, o in zip(y, ORD):
    p1, p5 = sub[o]["prod_top1"], sub[o]["prod_top5"]
    s1, s5 = sub[o]["sub_top1"], sub[o]["sub_top5"]
    ax.plot([p1, p5], [i + off] * 2, color=S.GREY, lw=1.3, solid_capstyle="round", zorder=1)
    ax.plot([p1], [i + off], "o", ms=4.2, color=S.GREY, zorder=3)
    ax.plot([p5], [i + off], "o", ms=4.2, mfc="white", mec=S.GREY, mew=1.0, zorder=3)
    ax.plot([s1, s5], [i - off] * 2, color=S.GREEN, lw=1.3, solid_capstyle="round", zorder=1)
    ax.plot([s1], [i - off], "o", ms=4.2, color=S.GREEN, zorder=3)
    ax.plot([s5], [i - off], "o", ms=4.2, mfc="white", mec=S.GREEN, mew=1.0, zorder=3)
    ax.annotate("", xy=(s1, i - off + 0.02), xytext=(p1, i + off - 0.02),
                arrowprops=dict(arrowstyle="-", lw=0.5, color=S.RULE), zorder=0)
    g0, g1 = p5 - p1, s5 - s1
    ax.text(103, i, "gap %.0f→%.0f" % (g0, g1), va="center", ha="left",
            fontsize=5.2, color=S.GREEN if g1 < g0 else S.VERM)
ax.set_yticks(y)
ax.set_yticklabels([o.replace("_", " ").replace("POOLED", "pooled") for o in ORD], fontsize=5.8)
ax.set_ylim(-0.95, len(ORD) - 0.35)
ax.set_xlim(28, 102)
ax.set_xlabel("accuracy (%)   ● top-1   ○ top-5")
ax.plot([], [], "o", ms=4.2, color=S.GREY, label="production")
ax.plot([], [], "o", ms=4.2, color=S.GREEN, label="subspace")
ax.legend(loc="lower left", handlelength=0.8, fontsize=5.4, borderpad=0.2,
          bbox_to_anchor=(0.0, -0.03))
ax.set_title("Shortlist hits become rank-1 hits", loc="left", pad=4)
ax.text(0.99, -0.155, "gap narrows in 4 of 5; blood is the exception",
        transform=ax.transAxes, ha="right", va="top", fontsize=5.2, color=S.INK2)
S.panel(ax, "c", dx=-0.40, dy=1.13)

tsv = ["panel\tlabel\tmetric\tvalue"]
tsv.append("a\tproduction (baseline)\ttop1\t%.1f" % BASE1)
tsv.append("a\tproduction (baseline)\ttop5\t%.1f" % BASE5)
for r in VARIANTS:
    tsv.append("a\t%s\ttop1\t%.1f" % (r[0], r[1]))
    tsv.append("a\t%s\ttop5\t%.1f" % (r[0], r[2]))
for k, v1, v5 in SWEEP:
    tsv.append("b\tmarkers_%d\ttop1\t%.1f" % (k, v1))
    tsv.append("b\tmarkers_%d\ttop5\t%.1f" % (k, v5))
for o in ORD:
    for m in ("prod_top1", "sub_top1", "prod_top5", "sub_top5"):
        tsv.append("c\t%s\t%s\t%.1f" % (o, m, sub[o][m]))
out = S.save(fig, HERE, NAME, "\n".join(tsv) + "\n")
print("\n".join(out))
