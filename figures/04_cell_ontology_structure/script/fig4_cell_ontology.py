#!/usr/bin/env python3
"""Figure 4 — what the Cell Ontology can and cannot support for automated curation.

a  CL's formal content: how much of each axis is actually populated
b  lineage disjointness: what can be asserted, and what real polyhierarchy blocks
c  refutation power: the upper bound, and what a reasoner actually returns
d  the CL graph mixes lineage with location, so graph distance is not similarity

Outputs fig4_cell_ontology.{svg,pdf,png} + fig4_cell_ontology_source_data.tsv
"""
import json
import math
import os
import sys
from collections import defaultdict, deque

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "_shared")))
import figstyle as S                                            # noqa: E402
S.apply()
import matplotlib.pyplot as plt                                 # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", "cellscribe_tool"))
RES = os.path.join(ROOT, "benchmark", "results")
SHORT = lambda i: i.rsplit("/", 1)[-1].replace("_", ":")

g = json.load(open(os.path.join(ROOT, "cl-full.json")))["graphs"][0]
LAB = {SHORT(n["id"]): n["lbl"] for n in g["nodes"] if n.get("lbl")}
PAR, CH = defaultdict(set), defaultdict(set)
for e in g["edges"]:
    if e.get("pred") == "is_a":
        PAR[SHORT(e["sub"])].add(SHORT(e["obj"]))
        CH[SHORT(e["obj"])].add(SHORT(e["sub"]))
CLT = [c for c in LAB if c.startswith("CL:")]
HASDEF = {SHORT(n["id"]) for n in g["nodes"]
          if (((n.get("meta") or {}).get("definition") or {}).get("val"))
          and SHORT(n["id"]).startswith("CL:")}
ld = defaultdict(lambda: {"genus": set(), "part_of": set(), "hpmp": set(), "go": set()})
for a in g.get("logicalDefinitionAxioms", []):
    cur = SHORT(a["definedClassId"])
    if not cur.startswith("CL:"):
        continue
    for gid in a.get("genusIds", []):
        ld[cur]["genus"].add(SHORT(gid))
    for r in a.get("restrictions") or []:
        if not r or not r.get("fillerId"):
            continue
        p, f = r.get("propertyId", ""), SHORT(r["fillerId"])
        if p.endswith("BFO_0000050") and f.startswith("UBERON:"):
            ld[cur]["part_of"].add(f)
        elif p.endswith("RO_0002104") and f.startswith("PR:"):
            ld[cur]["hpmp"].add(f)
        elif p.endswith("RO_0002215") and f.startswith("GO:"):
            ld[cur]["go"].add(f)

D = json.load(open(os.path.join(RES, "disjointness_experiment.json")))
N_CL = len(CLT)

fig = plt.figure(figsize=(S.W2, 3.5))
gs = fig.add_gridspec(1, 3, wspace=0.85, left=0.14, right=0.96, top=0.80, bottom=0.20,
                      width_ratios=[1.0, 0.95, 1.05])

# ---------------------------------------------------------------- a  formal content
ax = fig.add_subplot(gs[0, 0])
axes_ = [("is_a parent", sum(1 for c in CLT if PAR.get(c))),
         ("textual definition", len(HASDEF)),
         ("logical genus", sum(1 for c in CLT if ld.get(c, {}).get("genus"))),
         ("part_of (Uberon)", sum(1 for c in CLT if ld.get(c, {}).get("part_of"))),
         ("capable_of (GO)", sum(1 for c in CLT if ld.get(c, {}).get("go"))),
         ("surface marker (PRO)", sum(1 for c in CLT if ld.get(c, {}).get("hpmp"))),
         ("disjointness", 35)]
axes_.sort(key=lambda r: r[1])
y = np.arange(len(axes_))
frac = [100.0 * v / N_CL for _, v in axes_]
cols = [S.VERM if f < 15 else (S.BLUE if f < 60 else S.GREEN) for f in frac]
ax.barh(y, frac, height=0.66, color=cols, linewidth=0)
for i, ((lab, v), f) in enumerate(zip(axes_, frac)):
    ax.text(f + 1.5, i, "%.0f%%  (%d)" % (f, v), va="center", fontsize=5.4, color=S.INK2)
ax.set_yticks(y)
ax.set_yticklabels([lab for lab, _ in axes_], fontsize=5.8)
ax.set_xlim(0, 128)
ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xlabel("%% of the %d labelled CL terms" % N_CL)
ax.set_title("What CL actually encodes", loc="left", pad=4)
ax.text(0.98, 0.28, "markers and disjointness\nare nearly empty", transform=ax.transAxes,
        ha="right", va="top", fontsize=5.4, color=S.VERM, linespacing=1.3)
S.panel(ax, "a", dx=-0.52, dy=1.13)

# ---------------------------------------------------------------- b  disjointness
ax = fig.add_subplot(gs[0, 1])
safe, blocked = D["D1"]["safe_pairs"], D["D1"]["blocked_pairs"]
ax.bar([0], [safe], width=0.55, color=S.BLUE, linewidth=0)
ax.bar([0], [blocked], bottom=[safe], width=0.55, color=S.GREY, linewidth=0)
ax.bar([1], [35], width=0.55, color=S.VERM, linewidth=0)
ax.set_xticks([0, 1])
ax.set_xticklabels(["lineage pairs\namong 23 lineages", "asserted in CL\n(all cell types)"],
                   fontsize=5.6)
ax.set_ylabel("number of pairs")
ax.text(0.33, safe / 2, "%d assertable" % safe, ha="left", va="center",
        fontsize=5.4, color=S.BLUE)
ax.text(0.33, safe + blocked / 2, "%d blocked by\nreal biological overlap" % blocked,
        ha="left", va="center", fontsize=5.4, color=S.INK2, linespacing=1.25)
ax.text(1, 45, "35", ha="center", fontsize=6, color=S.VERM)
ax.set_ylim(0, 385)
ax.text(0.5, 0.90, "microglia are both leukocyte and glia;\nependymal cells both epithelial and glia",
        transform=ax.transAxes, ha="center", va="top", fontsize=5.2, color=S.INK2,
        linespacing=1.3)
ax.set_title("Cell identity is not a partition", loc="left", pad=4)
S.panel(ax, "b", dx=-0.34, dy=1.13)

# ---------------------------------------------------------------- c  refutation power
ax = fig.add_subplot(gs[0, 2])
lab = ["errors a reasoner\ncould catch in principle", "errors actually\nrefuted by ELK"]
val = [100 * D["D2"]["A1"]["rate"], 0.0]
ax.barh([1, 0], val, height=0.5, color=[S.BLUE, S.VERM], linewidth=0)
ax.set_yticks([1, 0])
ax.set_yticklabels(lab, fontsize=5.6)
ax.set_xlim(0, 100)
ax.set_xlabel("% of observed genus errors")
ax.text(val[0] + 2, 1, "%.0f%%  upper bound" % val[0], va="center", fontsize=5.6, color=S.BLUE)
ax.text(2, 0, "0 of 28\n(all %d axioms injected)" % safe,
        va="center", fontsize=5.4, color=S.VERM, linespacing=1.3)
ax.set_ylim(-0.6, 1.7)
ax.set_title("A reasoner cannot refute a wrong genus", loc="left", pad=4)
ax.text(0.98, 0.03, "the exclusions exist —\nnothing entails the second class",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=5.2,
        color=S.INK2, linespacing=1.3)
S.panel(ax, "c", dx=-0.42, dy=1.13)

tsv = ["panel\tlabel\tvalue\tdenominator"]
for lab_, v in axes_:
    tsv.append("a\t%s\t%d\t%d" % (lab_, v, N_CL))
tsv.append("b\tassertable_lineage_pairs\t%d\t%d" % (safe, safe + blocked))
tsv.append("b\tblocked_by_overlap\t%d\t%d" % (blocked, safe + blocked))
tsv.append("b\tdisjointness_asserted_in_CL\t35\t%d" % N_CL)
tsv.append("c\tupper_bound_catchable_pct\t%.1f\t100" % val[0])
tsv.append("c\tactually_refuted\t0\t%d" % D["D3"]["candidates"])
out = S.save(fig, HERE, "fig4_cell_ontology", "\n".join(tsv) + "\n")
print("\n".join(out))
