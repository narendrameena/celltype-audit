#!/usr/bin/env python3
"""Figure 5 — the atlas–ontology granularity gap.

a  what the 3,537 CL cell-type terms are, and how few an atlas reaches
b  coverage measured against the IN-SCOPE ontology, not against all of CL
c  per organ: cell types the atlas annotates vs CL terms available in that tissue
d  reached terms are overwhelmingly organ-specific; only six are ubiquitous

Outputs figure/fig5_granularity_gap.{svg,pdf,png} + sourceData/*.tsv
"""
import glob
import json
import os
import re
import sys
from collections import Counter, defaultdict, deque

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
PAR, CH, TAX = defaultdict(set), defaultdict(set), defaultdict(set)
for e in g["edges"]:
    p, a, b = e.get("pred"), SHORT(e["sub"]), SHORT(e["obj"])
    if p == "is_a":
        PAR[a].add(b)
        CH[b].add(a)
    elif p and p.endswith("RO_0002162"):
        TAX[a].add(b)


def desc(c):
    s, q = set(), deque([c])
    while q:
        for k in CH.get(q.popleft(), ()):
            if k not in s:
                s.add(k)
                q.append(k)
    return s


CLT = {c for c in LAB if c.startswith("CL:")}
reached, organs_of = set(), defaultdict(set)
for p in glob.glob(os.path.join(RES, "heca_to_cl_*.json")):
    o = os.path.basename(p)[len("heca_to_cl_"):-len(".json")]
    for t, v in json.load(open(p))["types"].items():
        c = v["cl"][0]["curie"]
        reached.add(c)
        organs_of[c].add(o)

IN_VITRO = desc("CL:0001034") | {"CL:0001034"}
ABNORMAL = desc("CL:0001061") | {"CL:0001061"}
DEVSET = desc("CL:0011115") | desc("CL:0002321") | desc("CL:0000034")
NONHUM = re.compile(r"\(sensu|\(Mmus\)|\(mouse|\(rat\b|\(Xenopus|\(zebrafish|\(Drosophila", re.I)
DEVW = re.compile(r"progenitor|precursor|blast\b|embryonic|fetal|primordial|immature|neonat", re.I)
cat = Counter()
for c in CLT - reached:
    l = LAB.get(c, "")
    taxa = {SHORT(t) for t in TAX.get(c, ())}
    if NONHUM.search(l) or (taxa and "NCBITaxon:9606" not in taxa):
        cat["non-human /\nspecies-specific"] += 1
    elif c in IN_VITRO or c in ABNORMAL:
        cat["in vitro /\nabnormal"] += 1
    elif c in DEVSET or DEVW.search(l):
        cat["developmental /\nprecursor"] += 1
    else:
        cat["human, mature, in vivo\n— NOT reached"] += 1
inscope = cat["human, mature, in vivo\n— NOT reached"] + len(reached)

fig = plt.figure(figsize=(S.W2, 3.6))
gs = fig.add_gridspec(1, 3, wspace=0.60, left=0.155, right=0.985, top=0.80, bottom=0.20,
                      width_ratios=[1.0, 1.2, 1.0])

# ------------------------------------------------- a  what CL contains vs what is reached
ax = fig.add_subplot(gs[0, 0])
segs = [("reached by the atlas", len(reached), S.GREEN),
        ("human, mature, in vivo\n— NOT reached", cat["human, mature, in vivo\n— NOT reached"], S.VERM),
        ("developmental / precursor", cat["developmental /\nprecursor"], S.GREY),
        ("non-human / species-specific", cat["non-human /\nspecies-specific"], S.RULE),
        ("in vitro / abnormal", cat["in vitro /\nabnormal"], S.RULE)]
segs_sorted = sorted(segs, key=lambda r: r[1])
y = np.arange(len(segs_sorted))
ax.barh(y, [v for _, v, _ in segs_sorted], height=0.62,
        color=[c for _, _, c in segs_sorted], linewidth=0)
for i, (lab, v, c) in enumerate(segs_sorted):
    ax.text(v + 60, i, "%d  (%.0f%%)" % (v, 100.0 * v / len(CLT)), va="center",
            fontsize=5.6, color=S.INK2)
ax.set_yticks(y)
ax.set_yticklabels([lab for lab, _, _ in segs_sorted], fontsize=5.6)
ax.set_xlim(0, 3050)
ax.set_xlabel("Cell Ontology cell-type terms")
ax.set_title("CL contains %d cell-type terms" % len(CLT), loc="left", pad=4)
ax.text(0.97, 0.22, "in scope for a human atlas: %d\natlas reaches %.1f%% of them"
        % (inscope, 100.0 * len(reached) / inscope),
        transform=ax.transAxes, ha="right", va="bottom", fontsize=5.4,
        color=S.VERM, linespacing=1.3)
S.panel(ax, "a", dx=-0.62, dy=1.15)

# ------------------------------------------------- b  per-organ granularity gap
ax = fig.add_subplot(gs[0, 1])
summ = json.load(open(os.path.join(RES, "all_organs_summary.json")))["organs"]
pts = [(r["organ"], r["mapped"], r["ref_terms"]) for r in summ if r.get("ref_terms")]
xs = [p[1] for p in pts]
ys = [p[2] for p in pts]
ax.scatter(xs, ys, s=13, color=S.BLUE, linewidths=0, zorder=3)
lim = max(max(xs), max(ys)) * 1.12
ax.plot([0, lim], [0, lim], color=S.RULE, lw=0.7, ls="--", zorder=1)
ax.text(lim * 0.66, lim * 0.72, "equal granularity", fontsize=5.2, color=S.INK2, rotation=38)
for o, x, y in pts:
    if o in ("Lung", "Brain", "Blood", "Kidney", "Pancreas", "Bone_marrow"):
        ax.annotate(o.replace("_", " "), xy=(x, y), xytext=(3, 4),
                    textcoords="offset points", fontsize=5.2, color=S.INK2)
ax.set_xlim(0, lim)
ax.set_ylim(0, lim)
ax.set_xlabel("cell types the atlas annotates")
ax.set_ylabel("CL terms available in that tissue")
ax.set_title("The ontology is finer than the annotation", loc="left", pad=4)
ax.text(0.97, 0.06, "every organ sits above the line",
        transform=ax.transAxes, ha="right", fontsize=5.4, color=S.INK2)
S.panel(ax, "b", dx=-0.20, dy=1.15)

# ------------------------------------------------- c  organ specificity of reached terms
ax = fig.add_subplot(gs[0, 2])
cnt = Counter(len(v) for v in organs_of.values())
ks = sorted(cnt)
ax.bar(ks, [cnt[k] for k in ks], width=0.82, color=S.GREEN, linewidth=0)
ax.set_xlabel("organs a CL term is used in")
ax.set_ylabel("number of CL terms")
one = cnt[1]
ax.annotate("%d terms (%.0f%%)\nused in one organ only" % (one, 100.0 * one / len(reached)),
            xy=(1, one), xytext=(3.4, one * 0.82), fontsize=5.4, color=S.INK2,
            linespacing=1.3, arrowprops=dict(arrowstyle="-", lw=0.5, color=S.RULE))
ubiq = sorted([(c, len(v)) for c, v in organs_of.items() if len(v) >= 10], key=lambda z: -z[1])
SHORT_NAME = {"endothelial cell of lymphatic vessel": "lymphatic endothelial cell",
              "vascular associated smooth muscle cell": "vascular smooth muscle cell"}
ax.text(0.97, 0.58, "ubiquitous (≥10 organs):\n" + "\n".join(
            SHORT_NAME.get(LAB.get(c, c), LAB.get(c, c)) for c, _ in ubiq[:6]),
        transform=ax.transAxes, ha="right", va="top", fontsize=5.0,
        color=S.INK2, linespacing=1.35)
ax.set_title("Most reached terms are organ-specific", loc="left", pad=4)
S.panel(ax, "c", dx=-0.30, dy=1.15)

tsv = ["panel\tlabel\tvalue"]
for lab, v, _ in segs:
    tsv.append("a\t%s\t%d" % (lab.replace("\n", " "), v))
tsv.append("a\tin_scope_total\t%d" % inscope)
tsv.append("a\tin_scope_coverage_pct\t%.1f" % (100.0 * len(reached) / inscope))
for o, x, y in pts:
    tsv.append("b\t%s\tatlas_types=%d;cl_terms_in_tissue=%d" % (o, x, y))
for k in ks:
    tsv.append("c\torgans_%d\t%d" % (k, cnt[k]))
out = S.save(fig, HERE, "fig5_granularity_gap", "\n".join(tsv) + "\n")
print("\n".join(out))
print("in-scope %d, coverage %.1f%%" % (inscope, 100.0 * len(reached) / inscope))
