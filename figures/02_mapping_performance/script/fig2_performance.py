#!/usr/bin/env python3
"""Figure 2 — expression-based mapping of atlas cell types to Cell Ontology terms.

a  per-organ agreement with the HuBMAP HRA crosswalks
b  seven hand-curated gold organs, against the achievable ceiling
c  the top-1 / top-5 gap — why the tool ships ranked candidates, not one answer

WHAT PANEL a IS NOT. The HRA crosswalk is itself a lexical label -> CL mapping, so it
trusts the label: it records which CL term an expert judged the label STRING to name,
not what the cells are. Panel a is therefore agreement between two independent routes
to a CL term -- lexical and expression -- and not accuracy against ground truth. The
lexical resolver in cellscribe_tool/benchmark/cl_resolve.py reproduces the crosswalk on
312 of these 315 types, which is what a circular comparison looks like when both sides
read the same string.

Panel b is the measure that does not trust the label: a human read the marker genes and
decided. It is the lower and the honest number, and the two should be read together --
a disagreement in panel a can mean the mapper is wrong OR the annotation is.

Outputs fig2_performance.{svg,pdf,png} + fig2_performance_source_data.tsv
"""
import glob
import json
import os
import sys
from collections import defaultdict, deque

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "_shared")))
import figstyle as S                                            # noqa: E402
S.apply()
import matplotlib.pyplot as plt                                 # noqa: E402

BENCH = os.path.abspath(os.path.join(HERE, "..", "..", "..", "cellscribe_tool", "benchmark"))
RES = os.path.join(BENCH, "results")
CLJ = os.path.join(os.path.dirname(BENCH), "cl-full.json")
SHORT = lambda i: i.rsplit("/", 1)[-1].replace("_", ":")

g = json.load(open(CLJ))["graphs"][0]
PAR, CH = defaultdict(set), defaultdict(set)
for e in g["edges"]:
    if e.get("pred") == "is_a":
        PAR[SHORT(e["sub"])].add(SHORT(e["obj"]))
        CH[SHORT(e["obj"])].add(SHORT(e["sub"]))


def within(c, adj, d=3):
    out, fr = set(), {c}
    for _ in range(d):
        nx = set()
        for x in fr:
            nx |= adj.get(x, set())
        nx -= out
        out |= nx
        fr = nx
    return out


def ok(p, G):
    return any(p == gd or p in within(gd, PAR) or p in within(gd, CH) for gd in G)


XW = defaultdict(set)
for r in json.load(open(os.path.join(RES, "hra_crosswalks.json"))):
    XW[r["label"].lower()].add(r["cl"])

per_organ, pooled, bysize = [], [0, 0, 0], []
for p in sorted(glob.glob(os.path.join(RES, "heca_to_cl_*.json"))):
    organ = os.path.basename(p)[len("heca_to_cl_"):-len(".json")]
    M = json.load(open(p))["types"]
    r1 = r5 = n = 0
    for t, v in M.items():
        G = XW.get(t.lower())
        if not G or v["n_cells"] < 500:
            continue
        n += 1
        a = ok(v["cl"][0]["curie"], G)
        b = any(ok(c["curie"], G) for c in v["cl"])
        r1 += a
        r5 += b
        bysize.append((v["n_cells"], a))
    if n >= 6:
        per_organ.append({"organ": organ, "n": n, "t1": 100 * r1 / n, "t5": 100 * r5 / n})
    pooled[0] += r1
    pooled[1] += r5
    pooled[2] += n
per_organ.sort(key=lambda r: -r["t1"])

gold = {s["organ"]: s for s in json.load(open(os.path.join(RES, "organ_scores.json")))}

fig = plt.figure(figsize=(S.W2, 5.35))
gs = fig.add_gridspec(2, 2, hspace=0.95, wspace=0.30,
                      left=0.10, right=0.98, top=0.90, bottom=0.135,
                      height_ratios=[1.25, 1.0])

# --- a per-organ vs HRA crosswalks
ax = fig.add_subplot(gs[0, :])
x = np.arange(len(per_organ))
w = 0.38
ax.bar(x - w / 2, [r["t5"] for r in per_organ], width=w, color=S.BLUE, linewidth=0, label="top-5")
ax.bar(x + w / 2, [r["t1"] for r in per_organ], width=w, color=S.VERM, linewidth=0, label="top-1")
ax.set_xticks(x)
ax.set_xticklabels(["%s (n=%d)" % (r["organ"].replace("_", " "), r["n"]) for r in per_organ],
                   rotation=90, ha="center", va="top", fontsize=5.2)
ax.set_ylabel("agreement with crosswalk (%)")
ax.set_ylim(0, 105)
pt1 = 100 * pooled[0] / pooled[2]
pt5 = 100 * pooled[1] / pooled[2]
ax.axhline(pt1, color=S.VERM, lw=0.6, ls="--")
ax.axhline(pt5, color=S.BLUE, lw=0.6, ls="--")
ax.set_xlim(-0.8, len(per_organ) + 2.6)
ax.text(len(per_organ) + 0.1, pt5, "pooled\ntop-5 %.1f%%" % pt5, fontsize=5.6,
        color=S.BLUE, ha="left", va="center", linespacing=1.25)
ax.text(len(per_organ) + 0.1, pt1, "pooled\ntop-1 %.1f%%" % pt1, fontsize=5.6,
        color=S.VERM, ha="left", va="center", linespacing=1.25)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.19), ncol=2,
          handlelength=1.1, columnspacing=1.4)
ax.set_title("Agreement with HuBMAP HRA crosswalks — a label\u2192CL translation, not ground "
             "truth (%d cell types, %d organs)" % (pooled[2], len(per_organ)), loc="left", pad=3)
S.panel(ax, "a", dx=-0.055)

# --- b gold organs vs ceiling
ax = fig.add_subplot(gs[1, 0])
order = ["Pancreas", "Liver", "Blood", "Bone_marrow", "Lung", "Kidney", "Heart"]
HELD_OUT = {"Lung", "Kidney", "Heart"}   # curated after the method was fixed; never tuned on
og = [o for o in order if o in gold]
x = np.arange(len(og))
key = "n>=500 & achievable"
t1 = [gold[o].get(key, {}).get("top1", 0) for o in og]
t5 = [gold[o].get(key, {}).get("top5", 0) for o in og]
ceil = [gold[o].get("ceiling", 0) for o in og]
ax.bar(x - 0.2, t5, width=0.38, color=S.BLUE, linewidth=0, label="top-5")
ax.bar(x + 0.2, t1, width=0.38, color=S.VERM, linewidth=0, label="top-1")
ax.plot(x, ceil, marker="_", ms=13, mew=1.4, ls="none", color=S.GREY, label="ceiling")
ax.set_xticks(x)
ax.set_xticklabels([o.replace("_", " ") for o in og], fontsize=5.4, rotation=30,
                   ha="right", rotation_mode="anchor")
for i, o in enumerate(og):                       # mark the held-out organs
    if o in HELD_OUT:
        ax.plot([i], [-6], marker="^", ms=3.0, color=S.INK2, clip_on=False, mew=0)
ax.set_ylabel("accuracy (%)")
ax.set_ylim(0, 124)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.40), ncol=3,
          handlelength=1.0, columnspacing=1.2, fontsize=5.4)
hc_n = sum(gold[o].get("n>=500", {}).get("n", 0) for o in og)
hc_ok = sum(gold[o].get("n>=500", {}).get("n", 0) * gold[o].get("n>=500", {}).get("top1", 0) / 100.0
            for o in og)
hc_t1 = 100 * hc_ok / max(hc_n, 1)
ax.axhline(hc_t1, color=S.INK2, lw=0.6, ls=":")
ax.text(len(og) - 0.45, 116, "\u22ef pooled top-1 %.1f%% (n=%d, %d organs)"
        % (hc_t1, hc_n, len(og)), fontsize=5.3, color=S.INK2, ha="right", va="center")
ax.set_title("Hand-curated gold — the measure that does not trust the label",
             loc="left", pad=3)
ax.text(0.5, -0.60, "\u25b2 curated after the method was fixed; never tuned on",
        transform=ax.transAxes, fontsize=5.2, color=S.INK2, ha="center")
S.panel(ax, "b")

# --- c the top-1/top-5 gap
ax = fig.add_subplot(gs[1, 1])
gaps = sorted([(r["organ"], r["t5"] - r["t1"]) for r in per_organ], key=lambda z: z[1])
y = np.arange(len(gaps))
ax.barh(y, [v for _, v in gaps], height=0.68, color=S.GREEN, linewidth=0)
ax.set_yticks(y)
ax.set_yticklabels([o.replace("_", " ") for o, _ in gaps], fontsize=5.2)
ax.set_xlabel("top-5 minus top-1 (percentage points)")
med = np.median([v for _, v in gaps])
ax.axvline(med, color=S.INK2, lw=0.6, ls="--")
ax.text(med + 1, 0.2, "median %.0f pts" % med, fontsize=5.6, color=S.INK2)
ax.set_title("The correct term is usually in the shortlist", loc="left", pad=3)
S.panel(ax, "c", dx=-0.30)

tsv = ["# panel a = agreement between two label->CL routes (lexical crosswalk vs expression);",
       "# it is NOT accuracy against ground truth. panel b = hand-curated gold, label-independent.",
       "panel\torgan\tn_types\ttop1_pct\ttop5_pct\tceiling_pct"]
for r in per_organ:
    tsv.append("a\t%s\t%d\t%.1f\t%.1f\t" % (r["organ"], r["n"], r["t1"], r["t5"]))
tsv.append("a\tPOOLED\t%d\t%.1f\t%.1f\t" % (pooled[2], pt1, pt5))
tsv.append("b\tPOOLED_handcurated_n>=500\t%d\t%.1f\t\t" % (hc_n, hc_t1))
for o in og:
    s = gold[o].get(key, {})
    tsv.append("b\t%s\t%d\t%.1f\t%.1f\t%.1f" % (o, s.get("n", 0), s.get("top1", 0),
                                                s.get("top5", 0), gold[o].get("ceiling", 0)))
out = S.save(fig, HERE, "fig2_performance", "\n".join(tsv) + "\n")
print("\n".join(out))
