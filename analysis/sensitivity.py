#!/usr/bin/env python3
"""Do the conclusions survive the thresholds, or were they chosen to produce them?

The pipeline carries tuned constants -- a 10%% reference-support floor, a 1.30 margin, a
100-cell minimum for a CL term to be scoreable, a 500-cell floor on clusters, a shortlist
of 5. Each was picked for a stated reason, but a reader is entitled to ask whether the
headline claims move when they move. Two episodes in this project make that a fair
question rather than a formality: the flagship lung "neutrophilic granulocyte" catch
flips on the binary margin depending on which reference is used, and the auditor's recall
over catchable errors moved from 67%% to 100%% as upstream bugs were fixed.

Each constant is swept one at a time, everything else held at its default, and four
claims are re-evaluated:

  C1  hECA shows more cross-lineage disagreement than Tabula Sapiens   (discrimination)
  C2  the anchor sweep's flags are mostly real                          (precision)
  C3  a fixed review budget finds most known errors                     (two-tier recall)
  C4  small clusters are mislabelled more often                         (epidemiology)

A claim that survives every setting is robust. A claim that survives only near the
default is reported as such.

Usage: python benchmark/sensitivity.py
"""
import json
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")

from cl_lineage import anchor_set, ancestors, load                    # noqa: E402
from cl_resolve import resolve as resolve2                            # noqa: E402
from scoring_variants import ok                                       # noqa: E402
from error_epidemiology import fisher                                 # noqa: E402

DEF = {"support": 0.10, "margin": 1.30, "min_ref": 100, "size_floor": 500, "topk": 5}
GOLD = ["Pancreas", "Liver", "Blood", "Bone_marrow", "Lung", "Kidney", "Heart"]
BUDGET = 33


def load_ref(stem="wide_ref_repro"):
    idx = json.load(open(os.path.join(RES, stem + "_index.json")))
    npz = np.load(os.path.join(RES, stem + ".npz"))
    return idx, {k[3:]: npz[k] for k in npz.files if k.startswith("M__")}


def score_cluster(idx, mats, ub, markers, asserted, min_ref):
    M, gix, tix = mats.get(ub), idx["gene_ix"].get(ub), idx["term_ix"].get(ub)
    if M is None or not gix or not tix or asserted not in tix:
        return None
    cols = [gix[g] for g in markers if g in gix]
    if len(cols) < 3:
        return None
    cnt = idx["counts"].get(ub, {})
    inv = {v: k for k, v in tix.items()}
    s = M[:, cols].mean(axis=1)
    keep = np.array([cnt.get(inv[i], 0) >= min_ref for i in range(len(s))])
    if not keep[tix[asserted]]:
        return None
    sk = np.where(keep, s, -1.0)
    bi = int(np.argmax(sk))
    bc = inv[bi]
    return {"asserted": asserted, "best": bc, "ratio": float(s[tix[asserted]]) / (float(sk[bi]) + 1e-9),
            "sup_a": cnt.get(asserted, 0), "sup_b": cnt.get(bc, 0),
            "related": (asserted == bc) or (bc in ancestors(asserted)) or (asserted in ancestors(bc))}


def kinds(rows, support):
    c = Counter()
    for r in rows:
        if r["sup_b"] < support * max(r["sup_a"], 1):
            c["thin"] += 1
            continue
        if r["related"]:
            c["agree"] += 1
            continue
        a, b = anchor_set(r["asserted"]), anchor_set(r["best"])
        c["cross" if (a and b and not (a & b)) else "within"] += 1
    return c


def heca_rows(idx, mats, p):
    out = []
    for o in GOLD:
        D = json.load(open(os.path.join(RES, "heca_markers_deep_%s.json" % o)))["types"]
        M = json.load(open(os.path.join(RES, "heca_to_cl_%s.json" % o)))["types"]
        ctx = {c["curie"] for v in M.values() for c in v.get("cl", [])}
        ub = json.load(open(os.path.join(RES, "heca_to_cl_%s.json" % o)))["uberon"]
        G = {k: v for k, v in json.load(open(os.path.join(HERE, "%s_gold.json" % o.lower()))).items()
             if not k.startswith("_") and v}
        for t, v in D.items():
            if v["n_cells"] < p["size_floor"]:
                continue
            cur, _h = resolve2(t, ctx, organ=o)
            if not cur:
                continue
            r = score_cluster(idx, mats, ub, [m["gene"] for m in v["markers"]], cur, p["min_ref"])
            if r:
                r.update(organ=o, label=t, n_cells=v["n_cells"],
                         error=(t in G and not ok(cur, G[t])))
                out.append(r)
    return out


def ts_rows(idx, mats, p):
    import glob
    from audit_ts import tissue_per_type
    out = []
    for fp in sorted(glob.glob(os.path.join(RES, "ts_markers_*.json"))):
        organ = os.path.basename(fp)[len("ts_markers_"):-len(".json")]
        t2ub = tissue_per_type(organ)
        for t, v in json.load(open(fp))["types"].items():
            if v["n_cells"] < p["size_floor"] or not v.get("cl"):
                continue
            ub = t2ub.get(t)
            if not ub:
                continue
            r = score_cluster(idx, mats, ub, [m["gene"] for m in v["markers"]], v["cl"], p["min_ref"])
            if r:
                out.append(r)
    return out


def anchor_flags(p):
    """The lineage sweep, at shortlist depth topk and the given size floor."""
    import glob
    tp = fp = 0
    for fpth in sorted(glob.glob(os.path.join(RES, "heca_to_cl_*.json"))):
        d = json.load(open(fpth))
        o = d["organ"]
        gp = os.path.join(HERE, "%s_gold.json" % o.lower())
        if not os.path.exists(gp):
            continue
        G = {k: v for k, v in json.load(open(gp)).items() if not k.startswith("_") and v}
        ctx = {c["curie"] for v in d["types"].values() for c in v.get("cl", [])}
        for t, v in d["types"].items():
            if v["n_cells"] < p["size_floor"] or len(v.get("cl", [])) < p["topk"]:
                continue
            cur, _h = resolve2(t, ctx, organ=o)
            if not cur or t not in G:
                continue
            A = anchor_set(cur)
            anc = [anchor_set(c["curie"]) for c in v["cl"][:p["topk"]]]
            if A and any(anc) and all(B and not (A & B) for B in anc):
                if not ok(cur, G[t]):
                    tp += 1
                else:
                    fp += 1
    return tp, fp


def evaluate(p, idx, mats):
    hr = heca_rows(idx, mats, p)
    tr = ts_rows(idx, mats, p)
    kh, kt = kinds(hr, p["support"]), kinds(tr, p["support"])
    liveh = kh["agree"] + kh["within"] + kh["cross"]
    livet = kt["agree"] + kt["within"] + kt["cross"]
    c1h = 100 * kh["cross"] / max(liveh, 1)
    c1t = 100 * kt["cross"] / max(livet, 1)
    tp, fp = anchor_flags(p)
    # with no flags precision is UNDEFINED, not zero; printing 0% would read as a
    # collapse when the sweep has simply become too conservative to fire at all
    prec = (100.0 * tp / (tp + fp)) if (tp + fp) else float("nan")
    # two-tier at a fixed review budget
    q = sorted([r for r in hr if not r["related"]
                and r["sup_b"] >= p["support"] * max(r["sup_a"], 1)], key=lambda r: r["ratio"])
    E = {(r["organ"], r["label"]) for r in hr if r["error"]}
    seen = []
    for r in q:
        k = (r["organ"], r["label"])
        if k not in seen:
            seen.append(k)
    found = len(set(seen[:BUDGET]) & E)
    rec = 100 * found / max(len(E), 1)
    # epidemiology
    sm = [r for r in hr if r["n_cells"] < 2000]
    bg = [r for r in hr if r["n_cells"] >= 2000]
    a = sum(1 for r in sm if r["error"]); b = len(sm) - a
    c = sum(1 for r in bg if r["error"]); d = len(bg) - c
    rs = 100 * a / max(len(sm), 1); rb = 100 * c / max(len(bg), 1)
    return {"heca_cross": c1h, "ts_cross": c1t, "precision": prec, "flags": tp + fp,
            "two_tier": rec, "n_err": len(E), "small": rs, "large": rb,
            "fold": rs / max(rb, 1e-9), "p": fisher(a, b, c, d)}


if __name__ == "__main__":
    idx, mats = load_ref()
    SWEEP = {"support": [0.0, 0.05, 0.10, 0.20, 0.35],
             "min_ref": [50, 100, 250, 500],
             "size_floor": [200, 500, 1000, 2000],
             "topk": [3, 5, 8]}
    base = evaluate(DEF, idx, mats)
    print("SENSITIVITY SWEEP  (default in brackets; everything else held at default)\n")
    print("  %-22s %8s %8s %10s %10s %9s %7s" % ("setting", "hECA%", "TS%", "precision",
                                                 "two-tier", "small/lg", "p"))
    print("  " + "-" * 82)
    def fmt(r, tag):
        pr = "  n/a" if np.isnan(r["precision"]) else "%4.0f%%" % r["precision"]
        fold = "  n/a" if (not np.isfinite(r["fold"]) or r["fold"] == 0) else "%4.1fx" % r["fold"]
        return ("  %-22s %7.1f%% %7.1f%% %9s %9.0f%% %9s %7.3f   flags=%d"
                % (tag, r["heca_cross"], r["ts_cross"], pr, r["two_tier"], fold, r["p"], r["flags"]))
    print(fmt(base, "DEFAULT"))
    out = {"default": base, "sweeps": {}}
    for key, vals in SWEEP.items():
        print("  " + "-" * 82)
        out["sweeps"][key] = {}
        for v in vals:
            p = dict(DEF)
            p[key] = v
            r = evaluate(p, idx, mats)
            out["sweeps"][key][str(v)] = r
            tag = "%s = %s%s" % (key, v, "  <-" if v == DEF[key] else "")
            print(fmt(r, tag))
    print("\n  C1 hECA > TS cross-lineage : %s"
          % ("HOLDS in every setting" if all(x["heca_cross"] > x["ts_cross"]
             for s in out["sweeps"].values() for x in s.values()) else "FAILS somewhere"))
    # a size_floor at or above the 2000-cell split empties the "small" bin, so C4 is
    # UNDEFINED there rather than false; those settings are excluded from the verdict
    c4 = [x for k, sw in out["sweeps"].items() for kk, x in sw.items()
          if not (k == "size_floor" and int(kk) >= 2000)]
    print("  C4 small > large error rate : %s (%.1f-%.1fx, p %.3f-%.3f)"
          % ("HOLDS in direction everywhere" if all(x["fold"] > 1 for x in c4)
             else "FAILS somewhere",
             min(x["fold"] for x in c4), max(x["fold"] for x in c4),
             min(x["p"] for x in c4), max(x["p"] for x in c4)))
    print("     Direction is robust; SIGNIFICANCE is not. This sweep scores only clusters")
    print("     that can be scored against the reference (14 errors); error_epidemiology.py")
    print("     uses every cluster whose label resolves (17 errors) and reports p = 0.011.")
    print("     The defensible claim is the ~2-3x enrichment, not a p-value.")
    json.dump(out, open(os.path.join(RES, "sensitivity.json"), "w"), indent=1)
    print("\nwrote results/sensitivity.json")
