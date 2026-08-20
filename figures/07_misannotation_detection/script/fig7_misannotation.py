#!/usr/bin/env python3
"""Figure 7 - Detecting mis-annotated cell types with the Cell Ontology.

The detector asks one question of every atlas cell type: does the CL term its label
asserts share a lineage with the CL terms its OWN expression evidence points to?
When the two anchor sets are disjoint across the whole shortlist, the label and the
data disagree about what the cells are.

Panel a is the methodological point. The detector is only as good as its lineage
oracle, and a keyword classifier over label strings is not one: it is blind to terms
whose names carry no lineage word ('type B pancreatic cell') and wrong on terms whose
names mislead ('pericyte' -> muscle, when CL places it under connective tissue). The
oracle used here instead reads lineage off CL's is_a graph, as the set of anchor
terms a class descends from -- which handles polyhierarchy natively, since a
microglial cell is both haematopoietic and neural and must not conflict with either.

Panel d is the recall measurement, which only became possible once seven organs were
hand-curated from markers. It reports what the sweep CANNOT see as well as what it misses:
the test compares lineage ANCHOR sets, so an error that stays inside one lineage --
"neutrophilic granulocyte" whose markers say classical monocyte, both haematopoietic -- is
invisible to it by construction, and no tuning changes that.

Run:  python figures/07_misannotation_detection/script/fig7_misannotation.py
"""
import glob
import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "figures", "_shared"))
sys.path.insert(0, os.path.join(ROOT, "cellscribe_tool", "benchmark"))

import figstyle as F                                                  # noqa: E402
from cl_lineage import resolve, anchor_set                            # noqa: E402
from cl_resolve import resolve as resolve2                            # noqa: E402
from lineage import lineage_of                                        # noqa: E402

RES = os.path.join(ROOT, "cellscribe_tool", "benchmark", "results")
K, MINC = 5, 200

# 'stromal' in the keyword vocabulary and 'connective' in CL name the same lineage.
# Normalising them keeps panel a a test of correctness, not of vocabulary.
SYN = {"stromal": "connective"}

# Full probe set for the panel-a summary count; the listed rows are the ones drawn.
PROBE = ["type B pancreatic cell", "fat cell", "myofibroblast cell", "astrocyte",
         "tuft cell", "endothelial cell", "fibroblast", "pericyte", "mesothelial cell",
         "Kupffer cell", "juxtaglomerular cell", "enteroendocrine cell", "Schwann cell",
         "hepatocyte", "podocyte", "chondrocyte", "osteoblast", "microglial cell",
         "Purkinje cell", "mast cell", "platelet", "keratinocyte", "melanocyte",
         "adipocyte", "natural killer cell", "plasma cell"]
SHOW = ["type B pancreatic cell", "fat cell", "juxtaglomerular cell",
        "myofibroblast cell", "pericyte", "mesothelial cell",
        "microglial cell", "fibroblast", "Schwann cell"]

SHORT = {"haematopoietic": "haemato.", "connective": "connect.", "endothelial": "endothel.",
         "epithelial": "epithel.", "muscle": "muscle", "neural": "neural",
         "germ": "germ", "germline": "germline"}

# Manual review of each flag's marker genes places it in one of three classes. A flag is
# a disagreement, and a disagreement has more than one possible culprit.
CONFIRMED = {                      # markers are unambiguously of the EVIDENCE lineage
    ("Adipose", "Fibroblast"), ("Lung", "Tuft cell"),
    ("Adipose", "Lymphatic endothelial cell"),
    ("Kidney", "Juxtaglomerular cell"), ("Kidney", "Neuroglial cell")}
MAPPING_FAIL = {                   # markers SUPPORT the label; the CL candidates are wrong
    ("Skin", "Conventional dendritic cell 1"), ("Testis", "Smooth muscle cell"),
    ("Adrenal_gland", "Endothelial cell"), ("Adipose", "Committed preadipocyte")}


def verdict_of(x):
    k = (x["organ"], x["label"])
    return "confirmed" if k in CONFIRMED else ("mapping" if k in MAPPING_FAIL else "boundary")


def verdict(term):
    """Compare the two oracles on one CL term name."""
    cs = resolve(term)
    g = sorted(anchor_set(next(iter(cs)))) if len(cs) == 1 else []
    r = lineage_of(term)
    rn = SYN.get(r, r)
    if not g:
        return r, g, "silent"
    if r == "other":
        return r, g, "blind"
    return (r, g, "agree") if rn in g else (r, g, "wrong")


def organ_context():
    """CL terms CELLxGENE observes in each organ, used to break resolver ties."""
    ctx = {}
    for p in sorted(glob.glob(os.path.join(RES, "heca_to_cl_*.json"))):
        d = json.load(open(p))
        ctx[d["organ"]] = {c["curie"] for v in d["types"].values() for c in v.get("cl", [])}
    return ctx


def detect(oracle, ctx=None):
    """oracle: 'graph' (CL is_a anchor sets) or 'keyword' (string classifier)."""
    flags, tested = [], 0
    for p in sorted(glob.glob(os.path.join(RES, "heca_to_cl_*.json"))):
        d = json.load(open(p))
        for t, v in d["types"].items():
            if v["n_cells"] < MINC or len(v.get("cl", [])) < K:
                continue
            cand = v["cl"][:K]
            if oracle == "graph":
                cur, _ = resolve2(t, (ctx or {}).get(d["organ"]))
                if not cur:
                    continue
                A = anchor_set(cur)
                anc = [anchor_set(c["curie"]) for c in cand]
            else:
                A = {SYN.get(lineage_of(t), lineage_of(t))} - {"other"}
                anc = [{SYN.get(lineage_of(c["label"]), lineage_of(c["label"]))} - {"other"}
                       for c in cand]
            if not A or not any(anc):
                continue
            tested += 1
            if all(B and not (A & B) for B in anc):
                ev = sorted(set().union(*[B for B in anc if B]))
                flags.append(dict(organ=d["organ"], label=t, n=v["n_cells"],
                                  markers=v["markers"][:4], asserted=sorted(A),
                                  evidence=ev, cands=[c["label"] for c in cand]))
    return tested, flags


def main():
    F.apply()
    ctx = organ_context()
    t_g, f_g = detect("graph", ctx)
    t_k, f_k = detect("keyword")
    G = {(x["organ"], x["label"]) for x in f_g}
    Kk = {(x["organ"], x["label"]) for x in f_k}
    artifacts = [x for x in f_k if (x["organ"], x["label"]) not in G]
    shared = len(G & Kk)

    allv = {t: verdict(t) for t in PROBE}
    n_wrong = sum(1 for v in allv.values() if v[2] == "wrong")
    n_blind = sum(1 for v in allv.values() if v[2] == "blind")
    n_sil = sum(1 for v in allv.values() if v[2] == "silent")

    fig = plt.figure(figsize=(F.W2, 6.45))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.52, 1.0], height_ratios=[1.0, 1.32],
                          wspace=0.34, hspace=0.50,
                          left=0.005, right=0.995, top=0.925, bottom=0.055)

    # ---------------------------------------------------------------- panel a
    ax = fig.add_subplot(gs[0, 0])
    F.panel(ax, "a", dx=-0.012, dy=1.150)
    ax.set_title("A keyword oracle mis-reads the lineage\nCL already states",
                 pad=6, loc="left", x=0.052, linespacing=1.25)
    XK, XG = 0.66, 1.46
    y = np.arange(len(SHOW))[::-1]
    for term, yy in zip(SHOW, y):
        r, g, v = allv[term]
        ax.text(0.20, yy, term, ha="right", va="center", fontsize=5.9, color=F.INK)
        if v == "silent":
            ax.text(XK, yy, SHORT.get(r, r), ha="center", va="center", fontsize=5.9,
                    color=F.INK2)
            ax.text(XG, yy, "no anchor", ha="center", va="center", fontsize=5.9,
                    color=F.RULE, style="italic")
        else:
            ck = {"wrong": F.VERM, "blind": F.VERM, "agree": F.INK2}[v]
            ax.text(XK, yy, "—" if r == "other" else SHORT.get(r, r), ha="center",
                    va="center", fontsize=5.9, color=ck,
                    fontweight="bold" if v != "agree" else "normal")
            ax.text(XG, yy, "/".join(SHORT.get(x, x) for x in g), ha="center",
                    va="center", fontsize=5.9, color=F.BLUE,
                    fontweight="bold" if v != "agree" else "normal")
    ax.text(XK, len(SHOW) - 0.28, "keyword\nclassifier", ha="center", va="bottom",
            fontsize=5.9, color=F.INK2, linespacing=1.25)
    ax.text(XG, len(SHOW) - 0.28, "CL is_a\nanchor set", ha="center", va="bottom",
            fontsize=5.9, color=F.BLUE, fontweight="bold", linespacing=1.25)
    for yy, tag in ((y[0] + 0.5, None),):
        pass
    ax.plot([-1.06, 1.90], [len(SHOW) - 0.42] * 2, lw=0.5, color=F.RULE)
    ax.text(-1.13, y[1], "blind", fontsize=5.6, color=F.VERM, ha="center",
            va="center", rotation=90)
    ax.text(-1.13, y[4], "wrong", fontsize=5.6, color=F.VERM, ha="center",
            va="center", rotation=90)
    ax.plot([-1.06, -1.06], [y[2] - 0.34, y[0] + 0.34], lw=1.1, color=F.VERM,
            solid_capstyle="butt")
    ax.plot([-1.06, -1.06], [y[5] - 0.34, y[3] + 0.34], lw=1.1, color=F.VERM,
            solid_capstyle="butt")
    ax.plot([-1.06, -1.06], [y[8] - 0.34, y[6] + 0.34], lw=1.1, color=F.RULE,
            solid_capstyle="butt")
    ax.set_xlim(-1.30, 1.94)
    ax.set_ylim(-1.30, len(SHOW) + 0.95)
    ax.axis("off")
    ax.text(-1.26, -0.72, "over %d common CL terms: keyword oracle wrong on %d, blind on %d;\n"
            "graph oracle silent on %d. Polyhierarchy is kept — a microglial cell is\n"
            "both haematopoietic and neural, so it conflicts with neither."
            % (len(PROBE), n_wrong, n_blind, n_sil),
            fontsize=5.4, color=F.INK2, ha="left", va="top", linespacing=1.4)

    # ---------------------------------------------------------------- panel b
    ax = fig.add_subplot(gs[0, 1])
    F.panel(ax, "b", dx=-0.52, dy=1.150)
    ax.set_title("Which oracle\nraised the flag", pad=6, loc="left", linespacing=1.25)
    vals = [len(artifacts), shared, len(f_g) - shared]
    ax.bar(range(3), vals, width=0.60, color=[F.VERM, F.GREY, F.BLUE], lw=0)
    for i, v in enumerate(vals):
        ax.annotate("%d" % v, (i, v), xytext=(0, 4), fontsize=6.5, fontweight="bold",
                    textcoords="offset points", ha="center", color=F.INK)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["keyword\nonly", "both", "CL graph\nonly"], fontsize=5.9,
                       linespacing=1.3)
    ax.set_ylabel("atlas cell types flagged")
    ax.set_yticks([0, 2, 4, 6])
    ax.set_ylim(0, max(vals) * 1.90)
    ax.set_xlim(-0.70, 2.70)
    ax.text(0, max(vals) * 1.62, "all %d refuted by\ntheir own markers" % len(artifacts),
            fontsize=5.5, color=F.VERM, ha="center", va="center", linespacing=1.3)
    ax.text(-0.62, -0.215, "keyword-only flags include Schwann-cell precursors\n"
            "called by S100B and a uterine fibroblast called by\nEMX2 and CXCL12",
            transform=ax.transAxes, fontsize=5.4, color=F.INK2, ha="left", va="top",
            linespacing=1.4)

    # ---------------------------------------------------------------- panel c
    ax = fig.add_subplot(gs[1, 0])
    F.panel(ax, "c", dx=-0.335, dy=1.150)
    ax.set_title("Every flag, and which side of it is wrong\n"
                 "(%d flags from %d resolved cell types)" % (len(f_g), t_g),
                 pad=6, loc="left", linespacing=1.25)
    fl = sorted(f_g, key=lambda x: x["n"])
    COL = {"confirmed": F.VERM, "mapping": F.BLUE, "boundary": F.INK2}
    for i, x in enumerate(fl):
        vd = verdict_of(x)
        col = COL[vd]
        ax.plot([0, 1], [i, i], lw=0.6, color=F.RULE, zorder=1, solid_capstyle="butt")
        ax.scatter([0], [i], s=15, color=F.GREY, zorder=3, lw=0)
        ax.scatter([1], [i], s=15, color=col, zorder=3, lw=0,
                   alpha=1.0 if vd != "boundary" else 0.5)
        ax.text(-0.06, i, "%s  %s" % (x["organ"].replace("_", " "), x["label"]),
                ha="right", va="center", fontsize=5.8,
                color=F.INK if vd != "boundary" else F.INK2)
        ax.text(1.06, i, "/".join(SHORT.get(e, e) for e in x["evidence"]),
                ha="left", va="center", fontsize=5.8, color=col)
        ax.text(0.5, i + 0.26, ", ".join(x["markers"][:3]), ha="center", va="bottom",
                fontsize=5.0, color=F.INK2, style="italic")
    ax.set_xlim(-1.62, 2.42)
    ax.set_ylim(-0.95, len(fl) - 0.05)
    ax.set_yticks([])
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["asserted\nby the label", "implied by its\nown markers"],
                       fontsize=5.9, linespacing=1.3)
    for s in ("left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)
    nc = sum(1 for x in fl if verdict_of(x) == "confirmed")
    nm = sum(1 for x in fl if verdict_of(x) == "mapping")
    nb = len(fl) - nc - nm
    ax.scatter([], [], s=15, color=F.VERM, lw=0,
               label="markers of the evidence lineage \u2014 label is wrong (%d)" % nc)
    ax.scatter([], [], s=15, color=F.BLUE, lw=0,
               label="markers of the label's lineage \u2014 our mapping is wrong (%d)" % nm)
    ax.scatter([], [], s=15, color=F.INK2, alpha=0.5, lw=0,
               label="boundary between adjacent lineages (%d)" % nb)
    ax.legend(loc="upper left", bbox_to_anchor=(-0.46, -0.095), ncol=1,
              handletextpad=0.4, borderpad=0, labelspacing=0.42, fontsize=5.4)

    # ---------------------------------------------------------------- panel d
    WL = json.load(open(os.path.join(RES, "within_lineage.json")))
    qrows = WL["queue"]
    REC = {(r["organ"], r["label"]): r
           for r in json.load(open(os.path.join(RES, "auditor_recall.json")))}
    # one common denominator: the cell types BOTH tests can score
    keys = [(r["organ"], r["label"]) for r in qrows]
    E = [k for k, r in zip(keys, qrows) if r["error"]]
    anchor = [k for k in keys if REC.get(k, {}).get("flagged")]
    nE, nT = len(E), len(keys)

    ax = fig.add_subplot(gs[1, 1])
    F.panel(ax, "d", dx=-0.26, dy=1.10)
    ax.set_title("Errors found for review effort spent\n"
                 "(%d cell types, %d known errors)" % (nT, nE), pad=6, loc="left",
                 linespacing=1.25)
    # queue-only curve
    xs, ys = [0], [0]
    found = 0
    for i, r in enumerate(qrows, 1):
        found += bool(r["error"])
        xs.append(i)
        ys.append(found)
    ax.plot(xs, ys, lw=1.4, color=F.BLUE, label="marker queue, ranked")
    # two-tier: anchor flags reviewed first, then the queue
    seen, xs2, ys2, f2 = set(), [0], [0], 0
    for k in anchor:
        seen.add(k)
        f2 += k in E
        xs2.append(len(seen))
        ys2.append(f2)
    for k, r in zip(keys, qrows):
        if k in seen:
            continue
        seen.add(k)
        f2 += r["error"]
        xs2.append(len(seen))
        ys2.append(f2)
    ax.plot(xs2, ys2, lw=1.4, color=F.VERM, label="anchor sweep, then queue")
    ax.plot([0, nT], [0, nE], lw=0.7, ls="--", color=F.GREY, label="reviewing at random")
    K = 33
    yk = next((y for x, y in zip(xs2, ys2) if x >= K), ys2[-1])
    ax.scatter([K], [yk], s=26, color=F.VERM, zorder=5, lw=0)
    ax.annotate("review %d of %d (%.0f%%)\nfind %d of %d errors (%.0f%%)"
                % (K, nT, 100 * K / nT, yk, nE, 100 * yk / max(nE, 1)),
                (K, yk), xytext=(16, -20), textcoords="offset points", fontsize=5.4,
                color=F.VERM, linespacing=1.3,
                arrowprops=dict(arrowstyle="-", lw=0.5, color=F.RULE))
    ax.set_xlabel("cell types a curator reviews")
    ax.set_ylabel("known errors found")
    ax.set_xlim(0, nT)
    ax.set_ylim(0, nE + 0.8)
    ax.legend(loc="lower right", fontsize=5.3, handlelength=1.6, borderpad=0.2,
              labelspacing=0.35)
    ax.text(-0.155, -0.275,
            "the anchor sweep alone stops at %d errors \u2014 it compares lineage anchor\n"
            "sets, so an error inside one lineage is invisible to it. The marker queue\n"
            "asks instead which CL term best explains the cluster, and catches lung\n"
            "\u201cneutrophilic granulocyte\u201d (56,394 cells) whose best term is classical monocyte"
            % sum(1 for k in anchor if k in E),
            transform=ax.transAxes, fontsize=5.3, color=F.INK2, ha="left", va="top",
            linespacing=1.42)

    tsv = ["panel\torgan\tlabel\tn_cells\tmarkers\tasserted_lineage\tevidence_lineage\tcl_candidates\tnote"]
    for i, r in enumerate(qrows[:40], 1):
        tsv.append("d\t%s\t%s\t%d\t-\t%s\t%s\t-\tqueue_rank=%d;ratio=%.3f;error=%s;anchor_flagged=%s"
                   % (r["organ"], r["label"], r["n_cells"], r["asserted_name"],
                      r["best_term"], i, r["ratio"], r["error"],
                      REC.get((r["organ"], r["label"]), {}).get("flagged", False)))
    tsv.append("d\t-\tSUMMARY\t%d\t-\t-\t-\t-\t"
               "errors=%d;anchor_flags=%d;two_tier_reviewed=%d;two_tier_found=%d"
               % (nT, nE, len(anchor), K, yk))

    for x in fl:
        tsv.append("c\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s" % (
            x["organ"], x["label"], x["n"], ";".join(x["markers"]), "/".join(x["asserted"]),
            "/".join(x["evidence"]), ";".join(x["cands"]), verdict_of(x)))
    for t in PROBE:
        r, g, v = allv[t]
        tsv.append("a\t-\t%s\t-\t-\t%s\t%s\t-\t%s" % (t, r, "/".join(g) or "-", v))
    tsv.append("b\t-\tkeyword_only\t%d\t-\t-\t-\t-\ttested=%d" % (len(artifacts), t_k))
    tsv.append("b\t-\tboth_oracles\t%d\t-\t-\t-\t-\t-" % shared)
    tsv.append("b\t-\tgraph_only\t%d\t-\t-\t-\t-\ttested=%d" % (len(f_g) - shared, t_g))
    for x in artifacts:
        tsv.append("b\t%s\t%s\t%d\t%s\t%s\t%s\t-\tkeyword_only_refuted" % (
            x["organ"], x["label"], x["n"], ";".join(x["markers"]),
            "/".join(x["asserted"]), "/".join(x["evidence"])))

    out = F.save(fig, HERE, "fig7_misannotation", "\n".join(tsv) + "\n")
    print("oracle probe: wrong=%d blind=%d silent=%d of %d" % (n_wrong, n_blind, n_sil, len(PROBE)))
    print("tested graph=%d keyword=%d | keyword-only %d, both %d, graph-only %d"
          % (t_g, t_k, len(artifacts), shared, len(f_g) - shared))
    for p in out:
        print("  ", p)


if __name__ == "__main__":
    main()
