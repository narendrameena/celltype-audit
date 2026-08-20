#!/usr/bin/env python3
"""Figure 10 — where cell-type annotation errors concentrate.

Which annotation decisions are unreliable has not been measured, because it needs errors
that are both KNOWN and labelled with the properties to test against. Seven hand-curated
organs give 190 cell types of which 18 carry a label the markers reject.

a  cluster size. Clusters of 500-2,000 cells are mislabelled about three times as often
   as larger ones. The DIRECTION is robust -- sensitivity.py sweeps five pipeline
   thresholds and finds 1.7-2.5x enrichment in every setting -- but the p-value is not:
   it is 0.011 over every cluster whose label resolves and 0.138 over the subset that can
   also be scored against the reference. The defensible claim is the enrichment, not the
   significance, and the panel says so.
b  the confusion structure -- which lineage is mislabelled as which. The diagonal is what
   an anchor-set test cannot see, and it is where most of the mass sits.
c  the smooth-muscle pattern, counted separately in each evidence source rather than
   pooled, because the denominators are not comparable.

With 18 errors almost nothing else reaches significance. Panels b and c are described as
observations; only panel a carries a test. Intervals are Wilson, which behaves at this n.

Outputs figure/fig10_epidemiology.{svg,pdf,png} + sourceData/*.tsv
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "figures", "_shared"))
sys.path.insert(0, os.path.join(ROOT, "cellscribe_tool", "benchmark"))
import figstyle as S                                            # noqa: E402
S.apply()
import matplotlib.pyplot as plt                                 # noqa: E402
from matplotlib.patches import Patch                            # noqa: E402

from cl_lineage import anchor_set                               # noqa: E402
from error_epidemiology import wilson, fisher                   # noqa: E402

RES = os.path.join(ROOT, "cellscribe_tool", "benchmark", "results")
rows = json.load(open(os.path.join(RES, "auditor_recall.json")))
err = [r for r in rows if r["error"]]

fig = plt.figure(figsize=(S.W2, 3.05))
gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.25, 0.85], wspace=0.50,
                      left=0.075, right=0.985, top=0.825, bottom=0.315)

# ------------------------------------------------------------------ a size
ax = fig.add_subplot(gs[0, 0])
S.panel(ax, "a", dx=-0.28, dy=1.16)
_sm = [r for r in rows if r["n_cells"] < 2000]
_bg = [r for r in rows if r["n_cells"] >= 2000]
_rs = (sum(1 for r in _sm if r["error"]) / max(len(_sm), 1))
_rb = (sum(1 for r in _bg if r["error"]) / max(len(_bg), 1))
ax.set_title("Small clusters are mislabelled\n%.1f times as often" % (_rs / max(_rb, 1e-9)),
             loc="left", pad=5, linespacing=1.25)
CUTS = [(500, 2000, "0.5-2k"), (2000, 10000, "2-10k"),
        (10000, 50000, "10-50k"), (50000, 10 ** 9, ">50k")]
xs, rates, los, his, ns = [], [], [], [], []
for i, (lo, hi, lab) in enumerate(CUTS):
    sub = [r for r in rows if lo <= r["n_cells"] < hi]
    k = sum(1 for r in sub if r["error"])
    a, b = wilson(k, len(sub))
    xs.append(i); rates.append(100 * k / max(len(sub), 1))
    los.append(100 * a); his.append(100 * b); ns.append(len(sub))
ax.errorbar(xs, rates, yerr=[np.array(rates) - los, np.array(his) - np.array(rates)],
            fmt="none", elinewidth=0.9, capsize=2.2, ecolor=S.RULE, zorder=3)
for i, r in enumerate(rates):
    ax.scatter([i], [r], s=28, color=S.VERM if i == 0 else S.BLUE, zorder=4, lw=0)
ax.set_xticks(xs)
ax.set_xticklabels(["%s\n(n=%d)" % (c[2], n) for c, n in zip(CUTS, ns)], fontsize=5.6,
                   linespacing=1.3)
ax.set_ylabel("cell types mislabelled (%)")
ax.set_xlabel("cluster size")
ax.set_ylim(0, max(his) * 1.12)
ax.set_xlim(-0.5, len(CUTS) - 0.5)
sm = [r for r in rows if r["n_cells"] < 2000]
bg = [r for r in rows if r["n_cells"] >= 2000]
a1 = sum(1 for r in sm if r["error"]); b1 = len(sm) - a1
c1 = sum(1 for r in bg if r["error"]); d1 = len(bg) - c1
ax.text(0.97, 0.95, "<2k: %.1f%%     \u2265 2k: %.1f%%\nFisher p = %.4f"
        % (100 * a1 / len(sm), 100 * c1 / len(bg), fisher(a1, b1, c1, d1)),
        transform=ax.transAxes, ha="right", va="top", fontsize=5.5, color=S.INK2,
        linespacing=1.35)
ax.text(0, -0.50, "bars are Wilson 95%% intervals; %d errors in %d cell types.\n"
        "Direction holds across five swept thresholds (1.7-2.5\u00d7); the p-value does\n"
        "not (0.011 here, 0.138 on the reference-scoreable subset)."
        % (len(err), len(rows)), transform=ax.transAxes, fontsize=5.2, color=S.INK2,
        linespacing=1.4)

# ------------------------------------------------------------------ b confusion
ax = fig.add_subplot(gs[0, 1])
S.panel(ax, "b", dx=-0.22, dy=1.16)
ax.set_title("What is mislabelled as what", loc="left", pad=5)
LIN = ["haematopoietic", "epithelial", "connective", "muscle", "endothelial", "neural"]
SHORT = {"haematopoietic": "haemato.", "epithelial": "epithel.", "connective": "connect.",
         "muscle": "muscle", "endothelial": "endothel.", "neural": "neural"}
M = np.zeros((len(LIN), len(LIN)))
for r in err:
    a2 = sorted(anchor_set(r["label_cl"]))
    b2 = sorted(anchor_set(r["gold_cl"]))
    if not a2 or not b2:
        continue
    if a2[0] in LIN and b2[0] in LIN:
        M[LIN.index(a2[0]), LIN.index(b2[0])] += 1
ax.imshow(np.ma.masked_where(M == 0, M), cmap="Oranges", vmin=0, vmax=M.max(),
          aspect="auto")
for i in range(len(LIN)):
    for j in range(len(LIN)):
        if M[i, j]:
            ax.text(j, i, "%d" % M[i, j], ha="center", va="center", fontsize=6,
                    color="white" if M[i, j] > M.max() * 0.6 else S.INK, fontweight="bold")
    ax.add_patch(plt.Rectangle((i - 0.5, i - 0.5), 1, 1, fill=False, lw=0.9,
                               edgecolor=S.BLUE, zorder=5))
ax.set_xticks(range(len(LIN)))
ax.set_xticklabels([SHORT[l] for l in LIN], fontsize=5.4, rotation=40, ha="right",
                   rotation_mode="anchor")
ax.set_yticks(range(len(LIN)))
ax.set_yticklabels([SHORT[l] for l in LIN], fontsize=5.4)
ax.set_xlabel("what the markers say")
ax.set_ylabel("what the label asserts")
diag = int(np.trace(M))
ax.text(0, -0.52, "boxed diagonal = label and evidence share a lineage.\n"
        "%d of %d placed errors sit there, invisible to an anchor-set test."
        % (diag, int(M.sum())), transform=ax.transAxes, fontsize=5.2, color=S.INK2,
        linespacing=1.4)

# ------------------------------------------------------------------ c muscle pattern
ax = fig.add_subplot(gs[0, 2])
S.panel(ax, "c", dx=-0.34, dy=1.16)
ax.set_title("A smooth-muscle label is\nthe commonest single trap", loc="left", pad=5,
             linespacing=1.25)
CX = json.load(open(os.path.join(RES, "cross_atlas_confirmation.json")))
AB = json.load(open(os.path.join(RES, "audit_baseline.json")))
conf = [r for r in CX if r["verdict"] == "CONFIRMED"]
dis = [r for r in AB["rows"] if not r["ct_agrees_label"]]
src = [("Tabula Sapiens\nconfirmed", sum(1 for r in conf if "muscle" in r["asserted_anchors"]), len(conf)),
       ("CellTypist\ndisagreements", sum(1 for r in dis if "muscle" in r["label_anchors"]), len(dis)),
       ("hand-curated\nerrors", sum(1 for r in err if "muscle" in anchor_set(r["label_cl"])), len(err))]
y = np.arange(len(src))[::-1]
for (name, k, tot), yy in zip(src, y):
    ax.barh(yy, 100 * k / tot, height=0.5, color=S.VERM, lw=0)
    ax.barh(yy, 100, height=0.5, color=S.FAINT, zorder=0, lw=0)
    ax.text(101, yy, "%d/%d" % (k, tot), va="center", fontsize=5.6, color=S.INK2)
    ax.text(-3, yy, name, ha="right", va="center", fontsize=5.5, linespacing=1.25)
ax.set_xlim(0, 128)
ax.set_ylim(-0.65, len(src) - 0.35)
ax.set_yticks([])
ax.set_xticks([0, 50, 100])
ax.set_xlabel("% of that source's errors")
for sp in ("left",):
    ax.spines[sp].set_visible(False)
ax.text(-0.62, -0.56, "counted per source, not pooled \u2014 the\n"
        "denominators are different populations; the\n"
        "curated organs hold few smooth-muscle clusters.",
        transform=ax.transAxes, fontsize=5.2, color=S.INK2, ha="left", va="top",
        linespacing=1.4)

tsv = ["panel\tkey\tvalue"]
for (lo, hi, lab), n, r, l, h in zip(CUTS, ns, rates, los, his):
    tsv.append("a\tsize_%s\tn=%d;error_pct=%.1f;ci=%.1f-%.1f" % (lab, n, r, l, h))
tsv.append("a\tfisher_lt2k_vs_ge2k\t%.5f" % fisher(a1, b1, c1, d1))
for i, li in enumerate(LIN):
    for j, lj in enumerate(LIN):
        if M[i, j]:
            tsv.append("b\t%s->%s\t%d" % (li, lj, int(M[i, j])))
tsv.append("b\tdiagonal_within_lineage\t%d of %d" % (diag, int(M.sum())))
for name, k, tot in src:
    tsv.append("c\t%s\t%d of %d" % (name.replace("\n", " "), k, tot))
out = S.save(fig, HERE, "fig10_epidemiology", "\n".join(tsv) + "\n")
print("size p=%.4f | confusion diagonal %d/%d" % (fisher(a1, b1, c1, d1), diag, int(M.sum())))
print("\n".join(out))
