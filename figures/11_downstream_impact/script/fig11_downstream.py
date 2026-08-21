#!/usr/bin/env python3
"""Figure — one correction changes a published statistic.

Finding an error only matters if it propagates. The largest confirmed error is a
56,394-cell lung cluster labelled "Neutrophilic granulocyte" whose markers are a classical
monocyte's. This shows what that does to a composition statistic anyone might report, and
checks the corrected value against an atlas that was never involved.

a  the marker evidence, gene by gene: what a neutrophil requires and does not have
b  the neutrophil:monocyte ratio as published, as curated, and in an independent atlas
c  how much of the atlas carries a label the evidence contradicts

Outputs figure/fig11_downstream.{svg,pdf,png} + sourceData/*.tsv
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

RES = os.path.join(ROOT, "cellscribe_tool", "benchmark", "results")
D = json.load(open(os.path.join(RES, "downstream_impact.json")))
CX = json.load(open(os.path.join(RES, "cross_atlas_confirmation.json")))
REC = json.load(open(os.path.join(RES, "auditor_recall.json")))
RT = json.load(open(os.path.join(RES, "referee_tests.json")))["test_b"]

fig = plt.figure(figsize=(S.W2, 5.9))
gs = fig.add_gridspec(2, 2, width_ratios=[1.06, 1.0], height_ratios=[1.0, 1.0],
                      wspace=0.42, hspace=1.05,
                      left=0.085, right=0.985, top=0.90, bottom=0.145)

# ------------------------------------------------------------------ a marker evidence
ax = fig.add_subplot(gs[0, 0])
S.panel(ax, "a", dx=-0.30, dy=1.18)
ax.set_title("The cluster's own markers", loc="left", pad=5)
present = [("FCN1", 1), ("CD14", 1), ("S100A12", 1), ("MAFB", 1)]
absent = [("FCGR3B", 0), ("ELANE", 0)]
rows = present + absent
y = np.arange(len(rows))[::-1]
for (g, p), yy in zip(rows, y):
    ax.barh(yy, 1 if p else 0.0, height=0.55, color=S.VERM if p else S.FAINT, lw=0)
    if not p:
        ax.barh(yy, 1, height=0.55, color="none", lw=0.6, edgecolor=S.RULE, linestyle=":")
    ax.text(-0.06, yy, g, ha="right", va="center", fontsize=6,
            color=S.INK if p else S.INK2, style="normal" if p else "italic")
    ax.text(1.06, yy, "monocyte" if p else "required for neutrophil\nand ABSENT",
            ha="left", va="center", fontsize=5.3, color=S.VERM if p else S.INK2,
            linespacing=1.25)
ax.set_xlim(-0.55, 2.35)
ax.set_ylim(-0.75, len(rows) - 0.25)
ax.set_yticks([])
ax.set_xticks([])
for sp in ("left", "bottom"):
    ax.spines[sp].set_visible(False)
ax.text(-0.30, -0.30,
        "lung cluster labelled “Neutrophilic granulocyte”, 56,394 cells. A neutrophil\n"
        "without FCGR3B or ELANE is the anomaly that needs explaining. It is not\n"
        "degradation: the cluster detects a median %d genes per cell, MORE than the\n"
        "lung’s own monocytes (%d) and above the median cell type (%d). The three\n"
        "genuinely sparse clusters here run 278–420."
        % (RT["median_genes_flagged"], RT["median_genes_monocyte"],
           RT["median_across_types"]),
        transform=ax.transAxes, fontsize=5.3, color=S.INK2, ha="left", va="top",
        linespacing=1.4)

# ------------------------------------------------------- b what the sources called them
ax = fig.add_subplot(gs[0, 1])
S.panel(ax, "b", dx=-0.30, dy=1.18)
ax.set_title("What the source datasets called these cells", loc="left", pad=5)
nprov = RT["n_cells"]
prov = [("monocyte, macrophage,\nmyeloid or dendritic", RT["mononuclear_phagocyte"], S.VERM),
        ("other or unspecified", nprov - RT["mononuclear_phagocyte"] - RT["neutrophil"], S.GREY),
        ("neutrophil or\ngranulocyte", RT["neutrophil"], S.BLUE)]
yp = np.arange(len(prov))[::-1]
for (name, v, col), yy in zip(prov, yp):
    ax.barh(yy, 100 * v / nprov, height=0.5, color=col, lw=0)
    ax.text(100 * v / nprov + 1.2, yy, "%.1f%%   %s cells" % (100 * v / nprov, format(v, ",")),
            va="center", fontsize=5.5, color=S.INK)
    ax.text(-1.6, yy, name, ha="right", va="center", fontsize=5.4, linespacing=1.25)
ax.set_xlim(0, 148)
ax.set_ylim(-0.8, len(prov) - 0.2)
ax.set_yticks([])
ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xlabel("% of the 56,394 cells")
for sp in ("left",):
    ax.spines[sp].set_visible(False)
_top = sorted(RT["original_names"].items(), key=lambda x: -x[1])[:4]
ax.text(-0.30, -0.30,
        "hECA records each cell’s ORIGINAL name from the dataset it came from:\n"
        "%s,\n%s.\nThe harmonised label disagrees with the atlas’s own provenance. %s cells\n"
        "carry a mangled original name (“Monocytenocyte…”, %d distinct strings) — the\n"
        "signature of a broken string substitution in the harmonisation step."
        % (", ".join("%s %s" % (k, format(v, ",")) for k, v in _top[:2]),
           ", ".join("%s %s" % (k, format(v, ",")) for k, v in _top[2:]),
           format(RT["mangled_cells"], ","), RT["mangled_strings"]),
        transform=ax.transAxes, fontsize=5.3, color=S.INK2, ha="left", va="top",
        linespacing=1.4)

# ------------------------------------------------------------------ c the ratio
ax = fig.add_subplot(gs[1, 0])
S.panel(ax, "c", dx=-0.24, dy=1.18)
ax.set_title("Neutrophil : monocyte in lung", loc="left", pad=5)
b, a, t = D["before"], D["after"], D["tabula_sapiens"]
sets = [("as\npublished", b, S.VERM),
        ("one cluster\ncorrected", a, S.BLUE),
        ("independent\natlas", t, S.GREY)]
x = np.arange(len(sets))
for i, (name, d, col) in enumerate(sets):
    r = d["neutrophil"] / max(d["monocyte"], 1)
    ax.bar(i, max(r, 0.004), width=0.55, color=col, lw=0)
    ax.annotate("%.2f : 1" % r, (i, max(r, 0.004)), xytext=(0, 5), fontsize=6,
                textcoords="offset points", ha="center", color=S.INK, fontweight="bold")
ax.set_yscale("log")
ax.set_ylim(0.003, 6)
ax.set_xticks(x)
ax.set_xticklabels([s[0] for s in sets], fontsize=5.5, linespacing=1.3)
ax.set_ylabel("ratio (log)")
ax.axhline(1.0, lw=0.6, ls="--", color=S.RULE)
ax.text(2.42, 1.05, "parity", fontsize=5.2, color=S.INK2, ha="right")
ax.set_xlim(-0.62, 2.62)
ax.text(-0.34, -0.40, "As published, neutrophils OUTNUMBER\n"
        "monocytes in human lung; the independent\n"
        "atlas gives the reverse at ~14:1. One\n"
        "cluster accounts for an %d-fold gap."
        % round(D["fold_overstatement"]),
        transform=ax.transAxes, fontsize=5.3, color=S.INK2, ha="left", va="top",
        linespacing=1.4)

# ------------------------------------------------------------------ d affected cells
ax = fig.add_subplot(gs[1, 1])
S.panel(ax, "d", dx=-0.30, dy=1.18)
ax.set_title("Cells under a disputed label", loc="left", pad=5)
conf = [r for r in CX if r["verdict"] == "CONFIRMED"]
u = {(r["organ"], r["label"]): r["n_cells"] for r in conf}
for r in REC:
    if r["error"]:
        u[(r["organ"], r["label"])] = r["n_cells"]
# the CL-sibling artifact is excluded and said so
art = [k for k in u if k[1].startswith("Ventricle cardiomyocyte")]
n_art = sum(u[k] for k in art)
tot_all = sum(u.values())
tot = tot_all - n_art
atlas = 10071582
vals = [("confirmed by an\nindependent atlas", sum(r["n_cells"] for r in conf), S.VERM),
        ("plus hand-curated\nerrors", tot, S.BLUE)]
y = np.arange(len(vals))[::-1]
for (name, v, col), yy in zip(vals, y):
    ax.barh(yy, 100 * v / atlas, height=0.5, color=col, lw=0)
    ax.text(100 * v / atlas + 0.08, yy, "%s cells  (%.2f%%)" % (format(v, ","), 100 * v / atlas),
            va="center", fontsize=5.5, color=S.INK)
    ax.text(-0.12, yy, name, ha="right", va="center", fontsize=5.5, linespacing=1.25)
ax.set_xlim(0, 5.6)
ax.set_ylim(-0.85, len(vals) - 0.25)
ax.set_yticks([])
ax.set_xlabel("% of the atlas (10.07M cells)")
for sp in ("left",):
    ax.spines[sp].set_visible(False)
ax.text(-0.62, -0.40, "A further %s cells sit under a CL naming artifact\n"
        "(two sibling ventricular-cardiomyocyte terms with no is_a\n"
        "path between them) and are excluded as a labelling\nquestion, not an error."
        % format(n_art, ","), transform=ax.transAxes, fontsize=5.3, color=S.INK2,
        ha="left", va="top", linespacing=1.4)

tsv = ["panel\tkey\tvalue"]
for g, p in rows:
    tsv.append("a\t%s\t%s" % (g, "present" if p else "absent"))
for name, v, _ in prov:
    tsv.append("b\t%s\t%d cells (%.2f%% of %d)"
               % (name.replace("\n", " "), v, 100 * v / nprov, nprov))
tsv.append("b\tmangled_original_names\t%d cells in %d strings"
           % (RT["mangled_cells"], RT["mangled_strings"]))
tsv.append("a\tmedian_genes_detected\tflagged=%d;lung_monocyte=%d;median_cell_type=%d"
           % (RT["median_genes_flagged"], RT["median_genes_monocyte"], RT["median_across_types"]))
for name, d, _ in sets:
    tsv.append("c\t%s\tneutrophil=%d;monocyte=%d;ratio=%.4f"
               % (name.replace("\n", " "), d["neutrophil"], d["monocyte"],
                  d["neutrophil"] / max(d["monocyte"], 1)))
tsv.append("c\tfold_overstatement\t%.1f" % D["fold_overstatement"])
for name, v, _ in vals:
    tsv.append("d\t%s\t%d cells (%.3f%% of %d)" % (name.replace("\n", " "), v,
                                                   100 * v / atlas, atlas))
tsv.append("d\texcluded_cl_artifact\t%d" % n_art)
out = S.save(fig, HERE, "fig11_downstream", "\n".join(tsv) + "\n")
print("fold=%.0f  confirmed cells=%d  all-errors cells=%d (%.2f%% of atlas)"
      % (D["fold_overstatement"], sum(r["n_cells"] for r in conf), tot, 100 * tot / atlas))
print("provenance: %.1f%% mononuclear phagocyte, %.1f%% neutrophil; %d genes/cell vs %d in monocytes"
      % (100 * RT["mononuclear_phagocyte"] / nprov, 100 * RT["neutrophil"] / nprov,
         RT["median_genes_flagged"], RT["median_genes_monocyte"]))
print("\n".join(out))
