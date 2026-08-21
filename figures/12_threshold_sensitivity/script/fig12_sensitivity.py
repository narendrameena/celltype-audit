#!/usr/bin/env python3
"""Supplementary figure — do the conclusions survive the thresholds?

The pipeline carries five tuned constants. Each was chosen for a stated reason, but a
reader is entitled to ask whether the headline claims move when they move. Each constant
is swept one at a time, everything else held at its default, and four claims re-evaluated
through the SAME code path the main figures use (within_lineage.build_rows / queue), not a
re-implementation of it.

a  C1  hECA shows more cross-lineage disagreement than Tabula Sapiens
b  C2  the anchor sweep's flags are mostly real
c  C3  a fixed review budget of 33 cell types finds most known errors
d  C4  small clusters are mislabelled more often

Outputs figure/fig12_sensitivity.{svg,pdf,png} + sourceData/*.tsv
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
from matplotlib.lines import Line2D                             # noqa: E402

RES = os.path.join(ROOT, "cellscribe_tool", "benchmark", "results")
D = json.load(open(os.path.join(RES, "sensitivity.json")))
SW, DEF = D["sweeps"], D["default"]

# the swept constants, in the order the pipeline applies them
GROUPS = [("support", "reference-support floor", ["0.0", "0.05", "0.1", "0.2", "0.35"]),
          ("min_ref", "min. reference cells", ["50", "100", "250", "500"]),
          ("size_floor", "min. cluster cells", ["200", "500", "1000", "2000"]),
          ("topk", "shortlist depth", ["3", "5", "8"])]
DEFAULT_OF = {"support": "0.1", "min_ref": "100", "size_floor": "500", "topk": "5"}
SHORT = {"1000": "1k", "2000": "2k"}          # four-digit ticks collide at 5 pt

xs, labels, isdef, rows, bounds = [], [], [], [], []
x = 0.0
for key, _title, vals in GROUPS:
    lo = len(xs)                                 # index into xs, not an x coordinate
    for v in vals:
        xs.append(x)
        labels.append(SHORT.get(v, v))
        isdef.append(v == DEFAULT_OF[key])
        rows.append(SW[key][v])
        x += 1.35
    bounds.append((lo, len(xs) - 1))
    x += 1.0                                     # a gap between parameter groups
xs = np.array(xs)
N = len(xs)


def mark_na(ax, xi, text, y):
    """A setting where the claim is UNDEFINED, not false. The label sits on the axis floor
    with a guide line, so it cannot be misread as annotating the nearest data point."""
    ax.axvline(xi, lw=0.5, ls=":", color=S.RULE, zorder=0)
    ax.text(xi, y, text, fontsize=5.0, color=S.INK2, ha="center", va="bottom",
            linespacing=1.2, style="italic")


def frame(ax, letter, title, ticks=False):
    S.panel(ax, letter, dx=-0.115, dy=1.20)
    ax.set_title(title, pad=5, loc="left", linespacing=1.25)
    for lo, hi in bounds[:-1]:                   # separators sit in the gaps
        ax.axvline(xs[hi] + 1.18, lw=0.5, color=S.FAINT, zorder=0)
    ax.set_xlim(xs[0] - 1.0, xs[-1] + 1.0)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels if ticks else [""] * N, fontsize=5.4)
    if ticks:
        for t, d in zip(ax.get_xticklabels(), isdef):
            if d:
                t.set_fontweight("bold")
                t.set_color(S.INK)


fig = plt.figure(figsize=(S.W2, 4.5))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.24,
                      left=0.062, right=0.988, top=0.905, bottom=0.155)

# ------------------------------------------------------- a  C1 discrimination
ax = fig.add_subplot(gs[0, 0])
h = np.array([r["heca_cross"] for r in rows])
t = np.array([r["ts_cross"] for r in rows])
frame(ax, "a", "C1  cross-lineage disagreement, atlas vs atlas")
for xi, hi, ti in zip(xs, h, t):
    ax.plot([xi, xi], [ti, hi], lw=0.7, color=S.RULE, zorder=1)
ax.scatter(xs, h, s=14, color=S.VERM, zorder=3, lw=0, label="hECA v2.0")
ax.scatter(xs, t, s=14, color=S.GREY, zorder=3, lw=0, label="Tabula Sapiens")
ax.set_ylabel("cell types crossing a\nlineage boundary (%)", linespacing=1.3)
ax.set_ylim(-0.5, 9.2)
ax.legend(loc="upper right", handletextpad=0.25, borderpad=0.1, labelspacing=0.3,
          fontsize=5.6, scatterpoints=1)
ax.text(0.0, -0.115, "hECA is above Tabula Sapiens in all %d settings" % N,
        transform=ax.transAxes, fontsize=5.4, color=S.INK2)

# ------------------------------------------------------- b  C2 flag precision
ax = fig.add_subplot(gs[0, 1])
pr = np.array([r["precision"] for r in rows], dtype=float)
fl = np.array([r["flags"] for r in rows])
frame(ax, "b", "C2  precision of the anchor sweep's flags")
ok = ~np.isnan(pr)
ax.scatter(xs[ok], pr[ok], s=14, color=S.BLUE, zorder=3, lw=0)
ax.axhline(DEF["precision"], lw=0.6, ls="--", color=S.RULE, zorder=0)
ax.text(xs[-1] + 0.9, DEF["precision"] + 2.0, "default %.0f%%" % DEF["precision"],
        fontsize=5.3, color=S.INK2, ha="right")
for xi, p_, f_ in zip(xs, pr, fl):
    if np.isnan(p_):
        mark_na(ax, xi, "no flags\nfired", 39.5)
    elif f_ != DEF["flags"]:
        ax.text(xi, p_ - 4.5, "%d flags" % f_, fontsize=5.0, color=S.INK2, ha="center",
                va="top")
ax.set_ylabel("flags that are real errors (%)")
ax.set_ylim(38, 100)
ax.text(0.0, -0.115, "precision is stable except where the sweep stops firing "
        "— at 2–3 flags the denominator, not the method, is the story",
        transform=ax.transAxes, fontsize=5.4, color=S.INK2)

# ------------------------------------------------------- c  C3 recall at budget
ax = fig.add_subplot(gs[1, 0])
tt = np.array([r["two_tier"] for r in rows])
ne = np.array([r["n_err"] for r in rows])
nt = np.array([r["n_types"] for r in rows])
frame(ax, "c", "C3  errors found at a fixed review budget (33 cell types)", ticks=True)
rand = 100 * 33 / nt                                   # reviewing 33 at random
ax.plot(xs, rand, lw=0.8, ls="--", color=S.GREY, zorder=1)
ax.scatter(xs, tt, s=14, color=S.BLUE, zorder=3, lw=0)
ax.text(xs[0] - 0.8, rand[0] - 4.0, "reviewing at random", fontsize=5.3, color=S.INK2,
        va="top")
ax.set_ylabel("known errors found (%)")
ax.set_ylim(0, 78)
ax.text(0.0, -0.30, "between %.0f%% and %.0f%% of known errors, against %.0f–%.0f%% "
        "for the same effort spent at random" % (tt.min(), tt.max(), rand.min(), rand.max()),
        transform=ax.transAxes, fontsize=5.4, color=S.INK2)

# ------------------------------------------------------- d  C4 size enrichment
ax = fig.add_subplot(gs[1, 1])
fo = np.array([r["fold"] for r in rows])
pv = np.array([r["p"] for r in rows])
frame(ax, "d", "C4  mislabelling enrichment in small clusters", ticks=True)
ax.axhline(1.0, lw=0.6, ls="--", color=S.GREY, zorder=0)
ax.text(xs[-1] + 0.9, 1.06, "no enrichment", fontsize=5.3, color=S.INK2, ha="right")
live = fo > 0
ax.scatter(xs[live], fo[live], s=14, color=S.BLUE, zorder=3, lw=0)
for xi in xs[~live]:
    mark_na(ax, xi, "small bin\nis empty", 0.46)
ax.set_ylabel("error rate, <2,000 cells\nrelative to ≥2,000", linespacing=1.3)
ax.set_ylim(0.4, 3.15)
ax.text(0.0, -0.30, "the DIRECTION holds in every setting (%.1f–%.1f×); the "
        "p-value does not (%.2f–%.2f). The claim is the enrichment, not significance."
        % (fo[live].min(), fo[live].max(), pv[live].min(), pv[live].max()),
        transform=ax.transAxes, fontsize=5.4, color=S.INK2)

# group headers under the bottom row
for (key, title, vals), (lo, hi) in zip(GROUPS, bounds):
    mid = (xs[lo] + xs[hi]) / 2
    for col in (0, 1):
        a = fig.axes[2 + col]
        a.text(mid, -0.175, title, transform=a.get_xaxis_transform(), fontsize=5.3,
               color=S.INK, ha="center", va="top")

tsv = ["panel\tparameter\tvalue\tis_default\theca_cross_pct\tts_cross_pct\t"
       "precision_pct\tflags\ttwo_tier_recall_pct\tqueue_only_recall_pct\t"
       "n_types\tn_errors\tfold\tfisher_p"]
i = 0
for key, _t, vals in GROUPS:
    for v in vals:
        r = rows[i]
        tsv.append("a-d\t%s\t%s\t%s\t%.2f\t%.2f\t%s\t%d\t%.1f\t%.1f\t%d\t%d\t%s\t%.4f"
                   % (key, v, isdef[i], r["heca_cross"], r["ts_cross"],
                      "NA" if np.isnan(r["precision"]) else "%.1f" % r["precision"],
                      r["flags"], r["two_tier"], r["queue_only"], r["n_types"],
                      r["n_err"], "NA" if r["fold"] == 0 else "%.2f" % r["fold"], r["p"]))
        i += 1
S.save(fig, HERE, "fig12_sensitivity", "\n".join(tsv) + "\n")
print("C1 holds in %d/%d   C3 range %.0f-%.0f%%   C4 range %.1f-%.1fx  p %.3f-%.3f"
      % (int((h > t).sum()), N, tt.min(), tt.max(),
         fo[live].min(), fo[live].max(), pv[live].min(), pv[live].max()))
