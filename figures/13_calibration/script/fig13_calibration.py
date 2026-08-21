#!/usr/bin/env python3
"""Supplementary figure — a confidence score cannot gate assignment.

The expression mapper is right about 57% of the time and its top-5 contains the answer 86%
of the time, so the tempting fix is selective prediction: score how likely the rank-1 call
is to be right, accept above a threshold, abstain below. It does not work, and the way it
fails is the point — the score ranks better than chance, but there is no operating point
at which accepting a call unread would be defensible.

a  precision against coverage: 80% costs all but a fifth of the coverage, 90% is never met
b  ROC: the signal is real (AUC 0.65-0.71) and that is not enough
c  the most confident errors are biologically adjacent, which no score shape can detect

Outputs figure/fig13_calibration.{svg,pdf,png} + sourceData/*.tsv
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
D = json.load(open(os.path.join(RES, "calibration.json")))


def curve(conf, ok):
    """Precision and coverage when the top-k most confident calls are accepted."""
    o = np.argsort(-np.asarray(conf, float))
    y = np.asarray(ok, float)[o]
    k = np.arange(1, len(y) + 1)
    return k / len(y), np.cumsum(y) / k


def roc(conf, ok):
    o = np.argsort(-np.asarray(conf, float))
    y = np.asarray(ok, float)[o]
    tp = np.cumsum(y); fp = np.cumsum(1 - y)
    return np.r_[0, fp / max(fp[-1], 1)], np.r_[0, tp / max(tp[-1], 1)]


CX = [(r["conf"], r["correct"]) for r in D["per_type"]]
GD = [(r["conf"], r["correct"]) for r in D["gold_per_type"]]
SETS = [("HRA crosswalks (n=%d)" % len(CX), CX, S.VERM, D["auc_crosswalks"]),
        ("hand-curated gold (n=%d)" % len(GD), GD, S.BLUE, D["auc_gold_organs"])]

fig = plt.figure(figsize=(S.W2, 3.35))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.74, 1.30], wspace=0.86,
                      left=0.068, right=0.988, top=0.86, bottom=0.40)

# ------------------------------------------------------------ a precision vs coverage
ax = fig.add_subplot(gs[0, 0])
S.panel(ax, "a", dx=-0.28, dy=1.18)
ax.set_title("Useful precision costs almost all the coverage", loc="left", pad=5)
stat = {}
for name, dat, col, _a in SETS:
    cov, pre = curve(*zip(*dat))
    keep = cov >= 0.03                       # the first few points are pure noise
    ax.plot(100 * cov[keep], 100 * pre[keep], lw=1.4, color=col, label=name)
    i = int(np.argmax(pre[keep]))
    stat[name] = {"max": 100 * pre[keep][i], "at": 100 * cov[keep][i],
                  "cov80": 100 * cov[pre >= 0.80].max() if (pre >= 0.80).any() else None,
                  "n": len(dat)}
    ax.scatter([100 * cov[keep][i]], [100 * pre[keep][i]], s=16, color=col, zorder=5, lw=0)
for tgt, lab in ((80, "80%"), (90, "90%"), (95, "95%")):
    ax.axhline(tgt, lw=0.5, ls=":", color=S.RULE, zorder=0)
    ax.text(101, tgt, lab, fontsize=5.0, color=S.INK2, va="center")
g = stat["hand-curated gold (n=%d)" % len(GD)]
c = stat["HRA crosswalks (n=%d)" % len(CX)]
ax.annotate("%.1f%% at %.0f%% coverage" % (g["max"], g["at"]), (g["at"], g["max"]),
            xytext=(9, 5), textcoords="offset points", fontsize=5.2, color=S.BLUE)
ax.annotate("%.1f%% at %.0f%% coverage" % (c["max"], c["at"]), (c["at"], c["max"]),
            xytext=(14, -34), textcoords="offset points", fontsize=5.2, color=S.VERM,
            arrowprops=dict(arrowstyle="-", lw=0.4, color=S.RULE))
ax.set_xlabel("cell types accepted (% of the set)")
ax.set_ylabel("precision among accepted (%)")
ax.set_xlim(0, 100); ax.set_ylim(35, 100)
ax.legend(loc="lower left", fontsize=5.3, handlelength=1.4, borderpad=0.2,
          labelspacing=0.35)
ax.text(0, -0.40, "most confident calls accepted first. 80%% precision is\n"
        "reachable on the gold, but only over the top %.0f%% of\n"
        "cell types; 90%% is reached on neither set."
        % g["cov80"],
        transform=ax.transAxes, fontsize=5.3, color=S.INK2, linespacing=1.4)

# ------------------------------------------------------------ b ROC
ax = fig.add_subplot(gs[0, 1])
S.panel(ax, "b", dx=-0.30, dy=1.18)
ax.set_title("The signal is real", loc="left", pad=5)
ax.plot([0, 1], [0, 1], lw=0.7, ls="--", color=S.RULE, zorder=0)
for name, dat, col, a in SETS:
    x, y = roc(*zip(*dat))
    ax.plot(x, y, lw=1.4, color=col, label="%s  %.2f" % (name.split(" (")[0], a))
ax.text(0.97, 0.05, "score shape alone: %.2f" % D["auc_crosswalks_score_only"],
        fontsize=5.2, color=S.INK2, ha="right")
ax.set_xlabel("false-positive rate"); ax.set_ylabel("true-positive rate")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_aspect("equal", adjustable="box")
ax.legend(loc="upper left", bbox_to_anchor=(-0.02, 1.02), fontsize=5.2, handlelength=1.3,
          borderpad=0.15, labelspacing=0.3, title="AUC", title_fontsize=5.2)
ax.text(0, -0.52, "leave-one-organ-out. Ranking above\nchance is necessary for a gate,\n"
        "and nowhere near sufficient.",
        transform=ax.transAxes, fontsize=5.3, color=S.INK2, linespacing=1.4)

# ------------------------------------------------------------ c why it fails
ax = fig.add_subplot(gs[0, 2])
S.panel(ax, "c", dx=-0.20, dy=1.18)
ax.set_title("The confident errors are near-misses", loc="left", pad=5)
wrong = sorted([r for r in D["per_type"] if not r["correct"]],
               key=lambda r: -r["conf"])[:7]
y = np.arange(len(wrong))[::-1]
for r, yy in zip(wrong, y):
    ax.barh(yy, r["conf"], height=0.5, color=S.VERM, lw=0)
    short = lambda x: (x.replace(" cell", "").replace("progenitor", "prog.")
                       .replace("Bronchial smooth muscle", "Bronchial sm. muscle"))
    ax.text(-0.02, yy, "%s → %s" % (short(r["label"]), short(r["pred_label"])),
            ha="right", va="center", fontsize=5.0, color=S.INK)
    ax.text(r["conf"] + 0.015, yy, "%.2f%s" % (r["conf"], "  ✓top-5" if r["top5_ok"] else ""),
            ha="left", va="center", fontsize=5.2, color=S.INK2)
ax.set_xlim(0, 1.52); ax.set_ylim(-0.8, len(wrong) - 0.2)
ax.set_yticks([]); ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
ax.set_xlabel("model confidence in a call that is wrong")
for sp in ("left",):
    ax.spines[sp].set_visible(False)
n5 = sum(1 for r in wrong if r["top5_ok"])
ax.text(-0.26, -0.44, "the seven most confident wrong calls. %d of %d have the correct\n"
        "term inside their own top-5, so these are RANKING failures between\n"
        "adjacent types, not retrieval failures — and no score shape separates\n"
        "them from a hit." % (n5, len(wrong)), transform=ax.transAxes, fontsize=5.3, color=S.INK2,
        linespacing=1.4)

tsv = ["panel\tkey\tvalue"]
for name, dat, _c, a in SETS:
    cov, pre = curve(*zip(*dat))
    keep = cov >= 0.03
    tsv.append("a\t%s\tn=%d;max_precision=%.1f;at_coverage=%.0f;coverage_at_80pct_precision=%s"
               % (name, len(dat), stat[name]["max"], stat[name]["at"],
                  ("%.0f" % stat[name]["cov80"]) if stat[name]["cov80"] else "unreachable"))
    tsv.append("b\t%s\tAUC=%.4f" % (name, a))
tsv.append("b\tHRA crosswalks, score-shape features only\tAUC=%.4f"
           % D["auc_crosswalks_score_only"])
for r in wrong:
    tsv.append("c\t%s %s -> %s\tconf=%.4f;n_cells=%d;correct_in_top5=%s"
               % (r["organ"], r["label"], r["pred_label"], r["conf"], r["n_cells"],
                  bool(r["top5_ok"])))
S.save(fig, HERE, "fig13_calibration", "\n".join(tsv) + "\n")
print("gold max %.1f%%  crosswalk max %.1f%%  AUC gold %.3f  crosswalks %.3f (score-only %.3f)  top-5 among the 7: %d"
      % (g["max"], c["max"], D["auc_gold_organs"], D["auc_crosswalks"],
         D["auc_crosswalks_score_only"], n5))
