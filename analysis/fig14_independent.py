#!/usr/bin/env python3
"""Figure 14 - the audit on an atlas the expression reference has never seen.

Every other number in this paper is a lower bound for one reason: the atlases audited are
part of the corpus the reference is built from, so a cluster they mislabel contributes to
the wrong term's profile and helps hide its own error. Testing without that circularity
needs data the reference has not seen, and inside CELLxGENE that is nearly unobtainable,
because the reference IS CELLxGENE. Of 2,216 indexed datasets exactly ONE is human,
non-diseased, non-fetal, carries a tissue the reference covers, was not already audited
here, and postdates the reference snapshot: HLiCA's liver endothelial atlas, 86,530 cells
in 8 CL-annotated types. n = 1 is stated rather than hidden; the scarcity is the finding.

a  What each tool returns per cluster. The atlas distinguishes eight endothelial subtypes.
   The audit resolves all eight and scores five; CellTypist's liver vocabulary holds
   seventeen labels of which one is endothelial, so the eight collapse to a single call.
b  The cluster the audit flags. 22,202 cells labelled "endothelial cell of pericentral
   hepatic sinusoid" whose eight strongest markers are hepatocyte genes. No other cluster
   carries any of them, so this is not blanket ambient signal.
c  Why it was missed until now. The sweep needs the asserted term's anchors disjoint from
   EVERY candidate's, and a term with no anchor set - which asserts nothing - was vetoing
   the flag. Taking unanimity over the candidates that do assert a lineage nearly doubles
   the sweep's reach at unchanged precision.
d  A second consortium, where the clusters carry no name at all. HuBMAP publishes Leiden
   clusters without a cell-type assertion, so there is nothing to audit - but naming them
   is what the shortlist is for. On 8,278 kidney cells the top-ranked term is a kidney
   anatomical type for 6 of 7 clusters and appears somewhere in the top five for all 7.
   This is coherence, not accuracy: no gold exists for these clusters, so it says the
   shortlist lands in the right organ, not that it picks the right segment.

Outputs figure/fig14_independent.{svg,pdf,png} + sourceData/fig14_independent_source_data.tsv
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "figures", "_shared"))
DATA = os.environ.get("FIG14_DATA", os.path.join(ROOT, "results_independent"))

import figstyle as S                                            # noqa: E402
import matplotlib.pyplot as plt                                 # noqa: E402

S.apply()

A = json.load(open(os.path.join(DATA, "hlica_audit.json")))
C = json.load(open(os.path.join(DATA, "hlica_celltypist.json")))
V = json.load(open(os.path.join(DATA, "anchor_veto.json")))

recs = {r["atlas_label"]: r for r in A["annotations"]}
order = sorted(C, key=lambda t: -C[t]["n"])
short = {t: (t.replace("endothelial cell of ", "")
             .replace(" hepatic sinusoid", " sinusoid")) for t in order}

fig = plt.figure(figsize=(180 * S.MM, 205 * S.MM))
gs = fig.add_gridspec(3, 2, height_ratios=[0.86, 0.72, 0.90], width_ratios=[1.22, 1.0],
                      hspace=0.98, wspace=0.52, left=0.245, right=0.985,
                      top=0.905, bottom=0.070)

# ------------------------------------------------------------------ a
ax = fig.add_subplot(gs[0, :])
S.panel(ax, "a", dx=-0.118, dy=1.30)
ax.set_title("One atlas, two tools: what each returns for the same eight clusters",
             loc="left", pad=16)
ax.text(0, 1.075, "HLiCA liver endothelium, 86,530 cells, published after the reference "
                  "snapshot.  CellTypist returns “Endothelial cells” for all eight.",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=5.8, color=S.INK2)
y = np.arange(len(order))[::-1]
for t, yy in zip(order, y):
    n = C[t]["n"]
    r = recs.get(t)
    a = (r or {}).get("audit", {})
    flagged = bool(a.get("contradicted"))
    ax.barh(yy, n, height=0.58, color=S.VERM if flagged else S.FAINT,
            edgecolor=S.RULE, lw=0.4)
    ax.text(-800, yy, short[t][:32], ha="right", va="center", fontsize=6,
            color=S.INK if flagged else S.INK2)
    call = ("not audited, under the 500-cell floor" if not r else
            ("no reference profile" if a.get("ratio") is None
             else "audit → " + (a.get("best_term_label") or "-")
             .replace("endothelial cell of ", "endothelial, ")
             .replace(" hepatic sinusoid", " sinusoid")))
    ax.text(n + 600, yy, "%s    %s" % (format(n, ","), call[:44]),
            va="center", fontsize=5.6, color=S.INK2)
ax.set_xlim(0, 56000)
ax.set_ylim(-0.7, len(order) - 0.3)
ax.set_yticks([])
ax.set_xlabel("cells")
for sp in ("top", "right", "left"):
    ax.spines[sp].set_visible(False)
ax.scatter([], [], marker="s", s=20, color=S.VERM, label="flagged by the audit")
ax.legend(loc="lower right", frameon=False, fontsize=5.8, handletextpad=.5,
          bbox_to_anchor=(1.0, -0.02))

# ------------------------------------------------------------------ b
ax = fig.add_subplot(gs[1, 0])
S.panel(ax, "b", dx=-0.20, dy=1.34)
flag = [r for r in A["annotations"] if r["audit"].get("contradicted")][0]
mk = (flag.get("evidence") or {}).get("markers") or []
genes = [m["gene"] for m in mk[:8]]
pcs = [100 * m["pc_in"] for m in mk[:8]]
ax.bar(np.arange(len(genes)), pcs, width=0.62, color=S.VERM, lw=0)
ax.set_xticks(np.arange(len(genes)))
ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=5.8)
ax.set_ylabel("detected in % of the cluster's cells")
ax.set_ylim(0, 50)
ax.set_title("The flagged cluster's markers are hepatocyte genes",
             loc="left", pad=17, fontsize=7)
ax.text(0, 1.035, "22,202 cells labelled “endothelial cell of pericentral hepatic "
                  "sinusoid”.\nNo other cluster in the atlas carries any of them.",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=5.6, color=S.INK2,
        linespacing=1.35)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)

# ------------------------------------------------------------------ c
ax = fig.add_subplot(gs[1, 1])
S.panel(ax, "c", dx=-0.26, dy=1.34)
fires = [V["atlas_fires_strict"], V["atlas_fires_loose"]]
ax.bar([0, 1], fires, width=0.5, color=[S.GREY, S.BLUE], lw=0)
for i, v in enumerate(fires):
    ax.text(i, v + 0.5, str(v), ha="center", fontsize=7, color=S.INK)
ax.set_xticks([0, 1])
ax.set_xticklabels(["unanimity over\nall candidates",
                    "over candidates that\nassert a lineage"], fontsize=5.8)
ax.set_ylabel("flags raised across the atlas")
ax.set_ylim(0, max(fires) * 1.30)
ax.set_title("Terms that assert nothing vetoed flags", loc="left", pad=17, fontsize=7)
ax.text(0, 1.035, "On the ten curated organs %d → %d flags,\nat %d%% precision either way."
        % (V["gold_strict"], V["gold_loose"], V["gold_precision"]),
        transform=ax.transAxes, ha="left", va="bottom", fontsize=5.6, color=S.INK2,
        linespacing=1.35)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)

# ------------------------------------------------------------------ d
H = json.load(open(os.path.join(DATA, "hubmap_coherence.json")))
ax = fig.add_subplot(gs[2, :])
S.panel(ax, "d", dx=-0.118, dy=1.24)
ax.set_title("A second consortium, whose clusters carry no name to audit", loc="left", pad=16)
ax.text(0, 1.055, "HuBMAP kidney, 8,278 cells in 7 Leiden clusters with no cell-type "
        "assertion.  Naming them is what the shortlist is for.",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=5.8, color=S.INK2)
hc = H["clusters"]
yy = np.arange(len(hc))[::-1]
for r, y2 in zip(hc, yy):
    ok = r["top1_kidney"]
    ax.barh(y2, r["n"], height=0.58, color=S.GREEN if ok else S.FAINT,
            edgecolor=S.RULE, lw=0.4)
    ax.text(-40, y2, "cluster %s" % r["cluster"], ha="right", va="center", fontsize=6,
            color=S.INK if ok else S.INK2)
    ax.text(r["n"] + 30, y2, "%s    %s" % (format(r["n"], ","), r["top1"][:44]),
            va="center", fontsize=5.6, color=S.INK2)
ax.set_xlim(0, 2650)
ax.set_ylim(-0.7, len(hc) - 0.3)
ax.set_yticks([])
ax.set_xlabel("cells")
for sp in ("top", "right", "left"):
    ax.spines[sp].set_visible(False)
ax.scatter([], [], marker="s", s=20, color=S.GREEN,
           label="top-ranked term is a kidney anatomical type: %d of %d clusters, "
                 "and all %d carry one somewhere in the top five"
                 % (H["top1_kidney_specific"], H["n_clusters"], H["any5_kidney_specific"]))
ax.legend(loc="upper right", frameon=False, fontsize=5.6, handletextpad=.5,
          bbox_to_anchor=(1.005, 1.115))

tsv = ["panel\tkey\tvalue"]
for t in order:
    r = recs.get(t)
    a = (r or {}).get("audit", {})
    tsv.append("a\t%s\tcells=%d;audited=%s;scoreable=%s;flagged=%s;audit_best=%s;celltypist=%s"
               % (t, C[t]["n"], bool(r), (a.get("ratio") is not None) if r else False,
                  bool(a.get("contradicted")) if r else False,
                  (a.get("best_term_label") or "") if r else "", C[t]["celltypist"]))
for gname, pc in zip(genes, pcs):
    tsv.append("b\t%s\t%.1f%% of the flagged cluster's cells" % (gname, pc))
tsv.append("c\tatlas_flags_unanimity_all\t%d" % V["atlas_fires_strict"])
tsv.append("c\tatlas_flags_lineage_asserting_only\t%d" % V["atlas_fires_loose"])
tsv.append("c\tgold_flags\t%d -> %d at %d%% precision"
           % (V["gold_strict"], V["gold_loose"], V["gold_precision"]))
for r in hc:
    tsv.append("d\tcluster %s\t%d cells; top-1 %s; kidney-specific top-1 %s; in top-5 %s"
               % (r["cluster"], r["n"], r["top1"], r["top1_kidney"], r["any5_kidney"]))
out = S.save(fig, HERE, "fig14_independent", "\n".join(tsv) + "\n")
print(out if isinstance(out, str) else "\n".join(map(str, out)))
