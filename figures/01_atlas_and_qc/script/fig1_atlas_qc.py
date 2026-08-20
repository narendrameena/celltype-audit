#!/usr/bin/env python3
"""Figure 1 — hECA v2.0 atlas composition and QC for CL mapping.

a  cells per organ
b  cell types per organ: annotated / retained (>=50 cells) / mapped to CL
c  cluster-size ECDF, with the >=500-cell analysis threshold
d  CELLxGENE reference depth per tissue (CL terms with >=100 cells)
e  what is excluded, and why
f  mapping coverage per organ

Outputs fig1_atlas_qc.{svg,pdf,png} + fig1_atlas_qc_source_data.tsv
"""
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "_shared")))
import figstyle as S                                            # noqa: E402
S.apply()
import matplotlib.pyplot as plt                                 # noqa: E402

BENCH = os.path.abspath(os.path.join(HERE, "..", "..", "..", "cellscribe_tool", "benchmark"))
RES = os.path.join(BENCH, "results")

# ---------------------------------------------------------------- gather
rows, sizes = [], []
for p in sorted(glob.glob(os.path.join(RES, "heca_markers_*.json"))):
    organ = os.path.basename(p)[len("heca_markers_"):-len(".json")]
    M = json.load(open(p))
    meta, types = M["meta"], M["types"]
    cp = os.path.join(RES, "heca_to_cl_%s.json" % organ)
    mapped = len(json.load(open(cp))["types"]) if os.path.exists(cp) else 0
    n500 = sum(1 for v in types.values() if v["n_cells"] >= 500)
    sizes += [v["n_cells"] for v in types.values()]
    rows.append({"organ": organ, "cells": meta["n_cells"], "types": meta["n_types"],
                 "kept": meta["types_kept"], "mapped": mapped, "n500": n500,
                 "dropped_cells": meta.get("dropped_cells", 0),
                 "dropped": len(meta.get("dropped_not_cell_types", []))})
rows.sort(key=lambda r: -r["cells"])

summ = json.load(open(os.path.join(RES, "all_organs_summary.json"))) \
    if os.path.exists(os.path.join(RES, "all_organs_summary.json")) else {"organs": []}
refdepth = {r["organ"]: (r.get("ref_terms", 0), r.get("tissue", "")) for r in summ.get("organs", [])}

fig = plt.figure(figsize=(S.W2, 6.6))
gs = fig.add_gridspec(3, 2, hspace=0.55, wspace=0.34,
                      left=0.13, right=0.98, top=0.96, bottom=0.06)

# --- a  cells per organ
ax = fig.add_subplot(gs[0, 0])
top = rows[:14]
y = np.arange(len(top))
ax.barh(y, [r["cells"] for r in top], height=0.7, color=S.BLUE, linewidth=0)
ax.set_yticks(y)
ax.set_yticklabels([r["organ"].replace("_", " ") for r in top])
ax.invert_yaxis()
ax.set_xscale("log")
ax.set_xlabel("cells (log scale)")
ax.set_title("Atlas composition — 14 largest organs", loc="left", pad=3)
tot = sum(r["cells"] for r in rows)
ax.text(0.98, 0.04, "%d organs · %.1f M cells" % (len(rows), tot / 1e6),
        transform=ax.transAxes, ha="right", fontsize=6, color=S.INK2)
S.panel(ax, "a", dx=-0.30)

# --- b  types annotated / kept / mapped
ax = fig.add_subplot(gs[0, 1])
tb = rows[:14]
y = np.arange(len(tb))
h = 0.26
ax.barh(y - h, [r["types"] for r in tb], height=h, color=S.GREY, linewidth=0, label="annotated")
ax.barh(y, [r["kept"] for r in tb], height=h, color=S.BLUE, linewidth=0, label="retained (≥50 cells)")
ax.barh(y + h, [r["mapped"] for r in tb], height=h, color=S.GREEN, linewidth=0, label="mapped to CL")
ax.set_yticks(y)
ax.set_yticklabels([r["organ"].replace("_", " ") for r in tb])
ax.invert_yaxis()
ax.set_xlabel("cell types")
ax.set_title("Cell types retained and mapped", loc="left", pad=3)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=3,
          handlelength=1.1, borderpad=0.2, columnspacing=1.2)
S.panel(ax, "b", dx=-0.30)

# --- c  cluster-size ECDF
ax = fig.add_subplot(gs[1, 0])
s = np.sort(np.array(sizes))
ax.step(s, np.arange(1, len(s) + 1) / len(s), where="post", color=S.BLUE, lw=1.2)
ax.set_xscale("log")
ax.axvline(500, color=S.VERM, lw=0.8, ls="--")
frac = (s >= 500).mean()
ax.text(560, 0.16, "≥500 cells\n%.0f%% of types" % (100 * frac), fontsize=6, color=S.VERM)
ax.axvline(50, color=S.GREY, lw=0.8, ls=":")
ax.text(52, 0.86, "≥50 retained", fontsize=6, color=S.INK2)
ax.set_xlabel("cells per annotated cell type")
ax.set_ylabel("cumulative fraction")
ax.set_title("Cluster size governs whether markers are computable", loc="left", pad=3)
S.panel(ax, "c")

# --- d  reference depth
ax = fig.add_subplot(gs[1, 1])
rd = sorted([(o, v[0]) for o, v in refdepth.items() if v[0]], key=lambda x: -x[1])[:14]
if rd:
    y = np.arange(len(rd))
    ax.barh(y, [v for _, v in rd], height=0.7, color=S.PINK, linewidth=0)
    ax.set_yticks(y)
    ax.set_yticklabels([o.replace("_", " ") for o, _ in rd])
    ax.invert_yaxis()
ax.set_xlabel("CL terms with ≥100 cells in CELLxGENE")
ax.set_title("Reference depth bounds what can be recovered", loc="left", pad=3)
S.panel(ax, "d", dx=-0.30)

# --- e  exclusions
ax = fig.add_subplot(gs[2, 0])
drop_cells = sum(r["dropped_cells"] for r in rows)
lab = ["retained", "catch-all categories\n(Unclassified, etc.)"]
val = [tot - drop_cells, drop_cells]
ax.barh([0], [val[0]], height=0.5, color=S.BLUE, linewidth=0)
ax.barh([0], [val[1]], left=[val[0]], height=0.5, color=S.GREY, linewidth=0)
ax.set_yticks([])
ax.set_xlabel("cells")
ax.set_xlim(0, tot)
ax.text(val[0] / 2, 0, "%s\n%.2f M (%.1f%%)" % (lab[0], val[0] / 1e6, 100 * val[0] / tot),
        ha="center", va="center", fontsize=5.8, color="white")
ax.annotate("%s\n%.2f M (%.1f%%)" % (lab[1], val[1] / 1e6, 100 * val[1] / tot),
            xy=(val[0] + val[1] / 2, 0.28), xytext=(val[0] * 0.86, 0.62),
            fontsize=5.8, color=S.INK2, ha="center",
            arrowprops=dict(arrowstyle="-", lw=0.5, color=S.RULE))
ax.set_ylim(-0.5, 1.0)
ax.set_title("Excluded before marker computation", loc="left", pad=3)
S.panel(ax, "e")

# --- f  mapping coverage
ax = fig.add_subplot(gs[2, 1])
mapped_o = [r for r in rows if r["mapped"]]
unmapped_o = [r for r in rows if not r["mapped"]]
ax.bar([0], [len(mapped_o)], width=0.55, color=S.GREEN, linewidth=0)
ax.bar([1], [len(unmapped_o)], width=0.55, color=S.GREY, linewidth=0)
ax.set_xticks([0, 1])
ax.set_xticklabels(["mapped", "no CELLxGENE\ntissue equivalent"])
ax.set_ylabel("organs")
for i, v in enumerate([len(mapped_o), len(unmapped_o)]):
    ax.text(i, v + 0.5, str(v), ha="center", fontsize=6.5, color=S.INK)
if unmapped_o:
    ax.text(0.97, 0.97, "no equivalent:\n" + "\n".join(
                r["organ"].replace("_", " ") for r in unmapped_o[:6]),
            transform=ax.transAxes, ha="right", va="top",
            fontsize=5.6, color=S.INK2, linespacing=1.35)
ax.set_title("Mapping coverage across organs", loc="left", pad=3)
S.panel(ax, "f")

tsv = ["organ\tcells\ttypes_annotated\ttypes_retained\ttypes_mapped\ttypes_ge500\t"
       "cells_dropped_catchall\tref_CL_terms_ge100"]
for r in rows:
    tsv.append("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d" % (
        r["organ"], r["cells"], r["types"], r["kept"], r["mapped"], r["n500"],
        r["dropped_cells"], refdepth.get(r["organ"], (0, ""))[0]))
out = S.save(fig, HERE, "fig1_atlas_qc", "\n".join(tsv) + "\n")
print("\n".join(out))
