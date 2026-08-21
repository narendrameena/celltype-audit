#!/usr/bin/env python3
"""Figure 8 — the flags hold up against evidence that never saw this pipeline.

Internal validation cannot establish that a flagged annotation is actually wrong; only
an independent source can. Two are used here.

a  Tabula Sapiens. hECA's markers for each flagged cluster are scored across TS cell
   types of the same organ, and the best OFF-lineage match is set against the best
   ON-lineage one. Comparing the two atlases' marker LISTS does not work -- an NS-Forest
   marker discriminates against whatever else is in that atlas, and the two share no top
   markers even for macrophage -- so the query is turned around instead.
b  CellTypist, run on the same cells and translated to CL through the HRA crosswalk. It
   agrees with the atlas label on 179 of 187 well-powered types; the 8 disagreements are
   dominated by one systematic pattern.
c  Where all three lines of evidence converge.
d  Tabula Sapiens audited as a SUBJECT rather than a witness, by the same pipeline. Its
   disagreements are overwhelmingly about granularity; hECA's cross lineage boundaries
   three times as often. Read as a lower bound for TS -- it is one of the datasets behind
   the CELLxGENE reference it is scored against, so a cluster it mislabels helps the wrong
   term fit and hides its own error.

Outputs figure/fig8_external.{svg,pdf,png} + sourceData/*.tsv
"""
import json
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "figures", "_shared"))
import figstyle as S                                            # noqa: E402
S.apply()
import matplotlib.pyplot as plt                                 # noqa: E402

sys.path.insert(0, os.path.join(ROOT, "cellscribe_tool", "benchmark"))
RES = os.path.join(ROOT, "cellscribe_tool", "benchmark", "results")


def unauditable_ts():
    """TS clusters passing the size floor whose tissue has no reference at all.

    Computed, not asserted: this sentence used to carry a hardcoded 31 over three named
    tissues, which missed a single skin-of-chest cluster and was off by one.
    """
    import glob
    from audit_ts import tissue_per_type
    from cl_lineage import load as _load
    have = set(json.load(open(os.path.join(RES, "wide_ref_repro_index.json")))["term_ix"])
    n, tis = 0, Counter()
    for fp in sorted(glob.glob(os.path.join(RES, "ts_markers_*.json"))):
        organ = os.path.basename(fp)[len("ts_markers_"):-len(".json")]
        t2ub = tissue_per_type(organ)
        for t, v in json.load(open(fp))["types"].items():
            if v["n_cells"] < 500 or not v.get("cl"):
                continue
            ub = t2ub.get(t)
            if not ub or ub not in have:
                n += 1
                tis[organ] += 1
    return n, tis
CX = json.load(open(os.path.join(RES, "cross_atlas_confirmation.json")))
AB = json.load(open(os.path.join(RES, "audit_baseline.json")))

fig = plt.figure(figsize=(S.W2, 6.30))
gs = fig.add_gridspec(2, 2, width_ratios=[1.30, 1.05], height_ratios=[1.0, 0.90],
                      wspace=0.40, hspace=0.52,
                      left=0.005, right=0.995, top=0.925, bottom=0.055)

# ------------------------------------------------------------------ a cross-atlas
ax = fig.add_subplot(gs[0, 0])
S.panel(ax, "a", dx=-0.505, dy=1.145)
ax.set_title("An independent atlas is asked what\nthe flagged cells are", pad=6,
             loc="left", x=0.03, linespacing=1.25)
rows = sorted(CX, key=lambda r: (r["verdict"] != "CONFIRMED", -r["best_off_lineage"]))
y = np.arange(len(rows))[::-1]
for r, yy in zip(rows, y):
    conf = r["verdict"] == "CONFIRMED"
    on, off = r["best_on_lineage"], r["best_off_lineage"]
    ax.plot([on, off], [yy, yy], lw=0.7, color=S.RULE, zorder=1, solid_capstyle="butt")
    ax.scatter([on], [yy], s=15, color=S.GREY, zorder=3, lw=0)
    ax.scatter([off], [yy], s=17, color=S.VERM if conf else S.INK2, zorder=3, lw=0,
               alpha=1.0 if conf else 0.45)
    ax.text(-0.025, yy, "%s  %s" % (r["organ"].replace("_", " "), r["label"]),
            ha="right", va="center", fontsize=5.6, color=S.INK if conf else S.INK2)
    ax.text(max(off, on) + 0.022, yy, r["ts_top"][0]["type"][:24],
            ha="left", va="center", fontsize=5.2,
            color=S.VERM if conf else S.INK2, style="italic")
ax.set_xlim(-0.60, 1.12)
ax.set_ylim(-1.15, len(rows) - 0.35)
ax.set_yticks([])
ax.set_xticks([0, 0.2, 0.4, 0.6])
ax.set_xlabel("detection of hECA's markers in Tabula Sapiens")
for s in ("left",):
    ax.spines[s].set_visible(False)
ax.scatter([], [], s=17, color=S.VERM, lw=0,
           label="best OFF-lineage TS match (%d confirmed)"
                 % sum(1 for r in rows if r["verdict"] == "CONFIRMED"))
ax.scatter([], [], s=15, color=S.GREY, lw=0, label="best ON-lineage TS match")
ax.legend(loc="upper left", bbox_to_anchor=(-0.52, -0.105), ncol=1, handletextpad=0.4,
          borderpad=0, labelspacing=0.35, fontsize=5.4)

# ------------------------------------------------------------------ b celltypist
ax = fig.add_subplot(gs[0, 1])
S.panel(ax, "b", dx=-0.30, dy=1.145)
ax.set_title("CellTypist, run on the same cells\nand mapped to CL", pad=6, loc="left",
             linespacing=1.25)
dis = sorted([r for r in AB["rows"] if not r["ct_agrees_label"]], key=lambda r: -r["n_cells"])
yy = np.arange(len(dis))[::-1]
for r, k in zip(dis, yy):
    col = S.VERM if r["flagged"] else S.INK2
    ax.scatter([0], [k], s=14, color=S.GREY, lw=0, zorder=3)
    ax.scatter([1], [k], s=14, color=col, lw=0, zorder=3, alpha=1.0 if r["flagged"] else 0.5)
    ax.plot([0, 1], [k, k], lw=0.6, color=S.RULE, zorder=1, solid_capstyle="butt")
    ax.text(-0.06, k, "%s  %s" % (r["organ"].replace("_", " "), r["label"][:25]),
            ha="right", va="center", fontsize=5.3,
            color=S.INK if r["flagged"] else S.INK2)
    ax.text(1.06, k, r["celltypist"][:19], ha="left", va="center", fontsize=5.2, color=col)
ax.set_xlim(-1.62, 2.05)
ax.set_ylim(-1.35, len(dis) + 0.45)
ax.set_yticks([])
ax.set_xticks([])
for s in ("left", "bottom"):
    ax.spines[s].set_visible(False)
ax.tick_params(length=0)
ax.text(0.40, len(dis) - 0.55, "atlas label", ha="right", va="bottom", fontsize=5.5, color=S.INK2)
ax.text(0.60, len(dis) - 0.55, "CellTypist says", ha="left", va="bottom", fontsize=5.5, color=S.INK2)
ax.text(-1.60, -0.55, "agrees with the label on %d of %d well-powered types;\n"
        "%d disagreements, and %d of them are a smooth-muscle\n"
        "label that CellTypist reads as pericyte or fibroblast"
        % (AB["agree"], AB["n"], len(dis),
           sum(1 for r in dis if "muscle" in r["label_anchors"])),
        fontsize=5.3, color=S.INK2, ha="left", va="top", linespacing=1.4)
ax.scatter([], [], s=14, color=S.VERM, lw=0, label="also flagged by the auditor (%d)"
           % sum(1 for r in dis if r["flagged"]))
ax.legend(loc="upper left", bbox_to_anchor=(-0.30, -0.145), ncol=1, handletextpad=0.4,
          borderpad=0, labelspacing=0.3, fontsize=5.4)

# ------------------------------------------------------------------ c convergence
ax = fig.add_subplot(gs[1, 0])
S.panel(ax, "c", dx=-0.10, dy=1.145)
ax.set_title("Three independent lines of evidence,\none conclusion", pad=6, loc="left",
             linespacing=1.25)
cases = []
for r in CX:
    hit = [x for x in AB["rows"] if x["organ"] == r["organ"] and x["label"] == r["label"]]
    if r["verdict"] == "CONFIRMED" and hit and hit[0]["flagged"]:
        cases.append((r, hit[0]))
cases.sort(key=lambda z: -z[0]["n_cells"])
ax.axis("off")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
top = 0.97
for r, h in cases:
    ax.text(0.0, top, "%s  %s   (n=%d)" % (r["organ"], r["label"], r["n_cells"]),
            fontsize=6.0, color=S.INK, fontweight="bold", va="top")
    ax.text(0.02, top - 0.085, "label asserts", fontsize=5.3, color=S.INK2, va="top")
    ax.text(0.40, top - 0.085, "%s" % "/".join(r["asserted_anchors"]), fontsize=5.3,
            color=S.INK, va="top")
    lines = [("its own markers", ", ".join(r["markers"][:4])),
             ("Tabula Sapiens", "%s  (%.2f vs %.2f)" % (r["ts_top"][0]["type"][:22],
                                                        r["best_off_lineage"], r["best_on_lineage"])),
             ("CellTypist", h["celltypist"][:26])]
    dy = 0.085
    for i, (k, v) in enumerate(lines):
        ax.text(0.02, top - 0.085 - dy * (i + 1), k, fontsize=5.3, color=S.INK2, va="top")
        ax.text(0.40, top - 0.085 - dy * (i + 1), v, fontsize=5.3, color=S.VERM, va="top")
    ax.plot([0.0, 0.98], [top - 0.40, top - 0.40], lw=0.5, color=S.RULE)
    top -= 0.47
ax.text(0.0, top + 0.02, "%d of %d flags queryable against Tabula Sapiens are\n"
        "confirmed; none is refuted. Two are confirmed twice over."
        % (sum(1 for r in CX if r["verdict"] == "CONFIRMED"), len(CX)),
        fontsize=5.4, color=S.INK2, va="top", linespacing=1.4)

# ------------------------------------------------------------------ d atlas comparison
from collections import Counter
sys.path.insert(0, os.path.join(ROOT, "cellscribe_tool", "benchmark"))
from cl_lineage import anchor_set                                   # noqa: E402
TSA = json.load(open(os.path.join(RES, "audit_ts.json")))
WL = json.load(open(os.path.join(RES, "within_lineage.json")))["queue"]
IDX = json.load(open(os.path.join(RES, "wide_ref_index.json")))


def _kind(asserted, best, related, sup_a, sup_b):
    if sup_b < 0.10 * max(sup_a, 1):
        return "thin"
    if related:
        return "agree"
    a, b = anchor_set(asserted), anchor_set(best)
    return "cross" if (a and b and not (a & b)) else "within"


cts = Counter(_kind(r["asserted"], r["best"], r["related"],
                    r["support_asserted"], r["support_best"]) for r in TSA["rows"])
che = Counter()
for r in WL:
    o = IDX["organs"][r["organ"]]
    c = IDX["counts"].get(o["uberon"], {})
    che[_kind(r["asserted"], r.get("best_curie", ""), r["related_to_best"],
              c.get(r["asserted"], 0), c.get(r.get("best_curie", ""), 0))] += 1

ax = fig.add_subplot(gs[1, 1])
S.panel(ax, "d", dx=-0.20, dy=1.13)
ax.set_title("The same audit, applied to each atlas", loc="left", pad=5)
CATS = [("agree", "best term IS the label, or an is_a relative", S.GREY),
        ("within", "a different term in the SAME lineage", S.BLUE),
        ("cross", "a term in a DIFFERENT lineage", S.VERM)]
sets = [("Tabula Sapiens", cts), ("hECA", che)]   # drawn bottom-up, so hECA reads first
for i, (name, c) in enumerate(sets):
    tot = sum(v for k, v in c.items() if k != "thin")
    left = 0.0
    for key, _lab, col in CATS:
        w = 100.0 * c[key] / max(tot, 1)
        ax.barh(i, w, left=left, height=0.46, color=col, lw=0)
        if w > 6:
            ax.text(left + w / 2, i, "%d" % c[key], ha="center", va="center",
                    fontsize=5.8, color="white", fontweight="bold")
        left += w
    ax.text(-2, i, "%s\n(n=%d)" % (name, tot), ha="right", va="center",
            fontsize=5.9, linespacing=1.25, color=S.INK)
    ax.text(101, i, "%.1f%% cross-lineage" % (100.0 * c["cross"] / max(tot, 1)),
            ha="left", va="center", fontsize=5.5, color=S.VERM)
ax.set_yticks([])
ax.set_xlim(-30, 132)
ax.set_ylim(-0.72, len(sets) - 0.10)
ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xlabel("% of assessed cell types")
for sp in ("left",):
    ax.spines[sp].set_visible(False)
from matplotlib.patches import Patch                              # noqa: E402
ax.legend(handles=[Patch(facecolor=c, label=l) for _k, l, c in CATS],
          loc="upper left", bbox_to_anchor=(-0.30, -0.20), ncol=1, fontsize=5.4,
          handlelength=1.1, borderpad=0, labelspacing=0.35)
N_UNAUD, UNAUD_TIS = unauditable_ts()
ax.text(-0.30, -0.60,
        "same pipeline, same support guard, same reference. Tabula Sapiens is a LOWER\n"
        "bound: it is one of the datasets behind that reference, so an error of its own\n"
        "helps the wrong term fit. A further %d TS clusters sit in tissues CELLxGENE has\n"
        "no reference for at all (%s) and cannot be audited."
        % (N_UNAUD, ", ".join("%s %d" % (o.lower().replace("_", " "), c)
                              for o, c in UNAUD_TIS.most_common())),
        transform=ax.transAxes, fontsize=5.3, color=S.INK2, ha="left", va="top",
        linespacing=1.42)

tsv = ["panel\torgan\tlabel\tn_cells\tasserted\tts_best_match\toff_lineage\ton_lineage\tverdict\tcelltypist"]
for name, c in sets:
    tsv.append("d\t-\t%s\t%d\t-\t-\t-\t-\tagree=%d;within=%d;cross=%d;thin_discarded=%d\t-"
               % (name, sum(v for k, v in c.items() if k != "thin"),
                  c["agree"], c["within"], c["cross"], c["thin"]))
ctmap = {(r["organ"], r["label"]): r for r in AB["rows"]}
for r in rows:
    h = ctmap.get((r["organ"], r["label"]), {})
    tsv.append("a\t%s\t%s\t%d\t%s\t%s\t%.3f\t%.3f\t%s\t%s" % (
        r["organ"], r["label"], r["n_cells"], "/".join(r["asserted_anchors"]),
        r["ts_top"][0]["type"], r["best_off_lineage"], r["best_on_lineage"], r["verdict"],
        h.get("celltypist", "")))
for r in dis:
    tsv.append("b\t%s\t%s\t%d\t%s\t\t\t\tflagged=%s\t%s" % (
        r["organ"], r["label"], r["n_cells"], "/".join(r["label_anchors"]),
        r["flagged"], r["celltypist"]))
tsv.append("b\t-\tSUMMARY\t%d\t\t\t\t\tagree=%d;disagree=%d\t-" % (AB["n"], AB["agree"], len(dis)))
out = S.save(fig, HERE, "fig8_external", "\n".join(tsv) + "\n")
print("cross-atlas confirmed %d/%d | celltypist agree %d/%d | convergent cases %d"
      % (sum(1 for r in CX if r["verdict"] == "CONFIRMED"), len(CX), AB["agree"], AB["n"], len(cases)))
print("\n".join(out))
