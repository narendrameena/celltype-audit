#!/usr/bin/env python3
"""Stage 2: calibrated confidence and abstention.

Top-1 accuracy is 56.8% and top-5 is 85.7%, so assigning the rank-1 term everywhere is
wrong about four times in ten, while the right term is nearly always somewhere in the
shortlist. The fix is not a better point estimate but SELECTIVE PREDICTION: score how
likely the rank-1 call is to be right, auto-accept above a threshold, and abstain
below it so a curator sees a shortlist instead of a wrong label.

Confidence is fitted on the HuBMAP HRA expert crosswalks -- an independent gold
standard, not our own gold organs -- with LEAVE-ONE-ORGAN-OUT cross-validation, so the
reported precision/coverage is out-of-fold and the threshold is never chosen on the
data it is scored against. A threshold picked in-sample would look far better and mean
nothing.

Correctness follows the established benchmark rule: a prediction counts as right if it
is a gold term or within 3 is_a hops of one.

RESULT: THIS DOES NOT WORK, AND THE REASON MATTERS.

The score RANKS somewhat better than chance but cannot be turned into a gate.
Leave-one-organ-out AUC is 0.653 against the HRA crosswalks (n=315) on score-shape
features alone and 0.684 once ontology shape is added; against the hand-curated gold
(n=200, 7 organs) it is 0.707. That is real signal, and it is still useless for
assignment: precision plateaus at 73.9%% and NO threshold reaches 80%% at any coverage,
90%% or 95%% at any coverage on either set. The most confident single call is wrong, and
the top decile runs at 58%% precision, so the ranking is not even monotone where it
matters. The high-confidence errors are systematically BIOLOGICALLY ADJACENT
calls (Megakaryocyte -> platelet, Mast cell -> basophil), and nothing about the shape of
a score distribution can distinguish a near-miss from a hit.

The second finding reframes the benchmark. The HRA crosswalk is itself a lexical
label -> CL mapping, so scoring a lexical resolver against it is circular: cl_resolve
reproduces it on 312 of 315 types. That means the crosswalk measures label -> CL
TRANSLATION, not annotation correctness -- it trusts the label. The 56.8%% "accuracy" of
the expression mapper against it is therefore a DISAGREEMENT RATE between two routes to
a CL term, not an accuracy against truth. Against the hand-curated gold, where a human
read the markers, top-1 is 55.5%%.

So the expression mapper should not be the assigner and does not need a confidence gate
to become one. Lexical resolution assigns (stage 1, 98.8%% of cells); expression audits,
and its disagreements feed the contradiction sweep; abstention applies where lexical
fails (1.2%% of cells), and there the expression top-5 gives a curator a shortlist.

Usage: python benchmark/calibrate.py
"""
import glob
import json
import os
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
CLJ = os.path.join(os.path.dirname(HERE), "cl-full.json")
MINC = 500
SHORT = lambda i: i.rsplit("/", 1)[-1].replace("_", ":")     # noqa: E731

FEATS = ["top_score", "margin", "rel_margin", "spread", "entropy",
         "log_cells", "n_markers", "log_ref_n",
         # ontology shape of the PREDICTED term. Non-circular: these describe the
         # candidate, never its relation to the atlas label (which the gold encodes).
         "pred_depth", "pred_n_children", "pred_is_leaf", "anchor_agree", "n_anchors"]


def load_graph():
    g = json.load(open(CLJ))["graphs"][0]
    par, ch = defaultdict(set), defaultdict(set)
    for e in g["edges"]:
        if e.get("pred") == "is_a":
            par[SHORT(e["sub"])].add(SHORT(e["obj"]))
            ch[SHORT(e["obj"])].add(SHORT(e["sub"]))
    return par, ch


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


def features(v, ont=None):
    s = np.array([c["score"] for c in v["cl"]], dtype=float)
    s1 = float(s[0])
    s2 = float(s[1]) if len(s) > 1 else 0.0
    s5 = float(s[min(len(s), 5) - 1])
    p = np.clip(s[:5], 1e-9, None)
    p = p / p.sum()
    ent = float(-(p * np.log(p)).sum() / np.log(len(p))) if len(p) > 1 else 0.0
    f = [s1, s1 - s2, (s1 - s2) / (s1 + 1e-9), s1 - s5, ent,
         np.log10(max(v["n_cells"], 1)), len(v.get("markers", [])),
         np.log10(max(v["cl"][0].get("n", 1), 1))]
    if ont is not None:
        depth, nch, anch = ont
        top = v["cl"][0]["curie"]
        a0 = anch(top)
        agree = np.mean([1.0 if anch(c["curie"]) & a0 else 0.0 for c in v["cl"][:5]]) if a0 else 0.0
        distinct = len({frozenset(anch(c["curie"])) for c in v["cl"][:5]})
        f += [depth.get(top, 0), len(nch.get(top, ())),
              1.0 if not nch.get(top) else 0.0, float(agree), float(distinct)]
    return f


def build(use_ontology=True):
    par, ch = load_graph()
    ont = None
    if use_ontology:
        from cl_lineage import anchor_set
        depth = {}
        for c in set(par) | set(ch):
            d, fr, seen = 0, {c}, {c}
            while fr:
                nx = set()
                for x in fr:
                    nx |= par.get(x, set())
                nx -= seen
                if not nx:
                    break
                seen |= nx
                fr = nx
                d += 1
            depth[c] = d
        ont = (depth, ch, anchor_set)
    xw = defaultdict(set)
    for r in json.load(open(os.path.join(RES, "hra_crosswalks.json"))):
        xw[r["label"].lower()].add(r["cl"])

    def ok(p, G):
        return any(p == gd or p in within(gd, par) or p in within(gd, ch) for gd in G)

    X, y, organ, meta = [], [], [], []
    for p in sorted(glob.glob(os.path.join(RES, "heca_to_cl_*.json"))):
        d = json.load(open(p))
        for t, v in d["types"].items():
            G = xw.get(t.lower())
            if not G or v["n_cells"] < MINC or not v.get("cl"):
                continue
            X.append(features(v, ont))
            y.append(int(ok(v["cl"][0]["curie"], G)))
            organ.append(d["organ"])
            meta.append({"organ": d["organ"], "label": t, "n_cells": v["n_cells"],
                         "pred": v["cl"][0]["curie"], "pred_label": v["cl"][0]["label"],
                         "top5_ok": int(any(ok(c["curie"], G) for c in v["cl"]))})
    return np.array(X, float), np.array(y, int), np.array(organ), meta


def fit(X, y, l2=1.0, iters=4000, lr=0.12):
    """Plain logistic regression, gradient descent with L2. No sklearn in this env."""
    Xb = np.hstack([np.ones((len(X), 1)), X])
    w = np.zeros(Xb.shape[1])
    for _ in range(iters):
        z = Xb @ w
        pr = 1.0 / (1.0 + np.exp(-np.clip(z, -30, 30)))
        gr = Xb.T @ (pr - y) / len(y)
        gr[1:] += l2 * w[1:] / len(y)
        w -= lr * gr
    return w


def predict(w, X):
    Xb = np.hstack([np.ones((len(X), 1)), X])
    return 1.0 / (1.0 + np.exp(-np.clip(Xb @ w, -30, 30)))


def oof(X, y, organ):
    """Leave-one-organ-out out-of-fold confidence."""
    conf = np.zeros(len(y))
    mu, sd = X.mean(0), X.std(0) + 1e-9
    for o in np.unique(organ):
        te = organ == o
        tr = ~te
        if tr.sum() < 20 or len(np.unique(y[tr])) < 2:
            conf[te] = y[tr].mean() if tr.sum() else 0.5
            continue
        w = fit((X[tr] - mu) / sd, y[tr])
        conf[te] = predict(w, (X[te] - mu) / sd)
    return conf


def curve(conf, y):
    """Precision and coverage as the accept-threshold sweeps."""
    o = np.argsort(-conf)
    yy = y[o]
    k = np.arange(1, len(yy) + 1)
    return conf[o], np.cumsum(yy) / k, k / len(yy)


def at_precision(conf, y, target):
    thr, prec, cov = curve(conf, y)
    good = np.where(prec >= target)[0]
    if not len(good):
        return None
    i = good[-1]                       # deepest point still meeting the target
    return dict(threshold=float(thr[i]), precision=float(prec[i]),
                coverage=float(cov[i]), n_accepted=int(i + 1))


def gold_organ_set():
    """The non-circular calibration set: organs where a human read the markers."""
    from scoring_variants import ok as ok_gold
    par, ch = load_graph()
    from cl_lineage import anchor_set
    depth = {}
    for c in set(par) | set(ch):
        d, fr, seen = 0, {c}, {c}
        while fr:
            nx = set()
            for x in fr:
                nx |= par.get(x, set())
            nx -= seen
            if not nx:
                break
            seen |= nx
            fr = nx
            d += 1
        depth[c] = d
    ont = (depth, ch, anchor_set)
    X, y, organ = [], [], []
    # every curated organ. This list used to stop at the four that existed when the
    # calibration was first run, so lung, kidney and heart -- curated later, and the three
    # never tuned on -- were silently excluded from the one non-circular check in the study.
    for o in ("Pancreas", "Liver", "Blood", "Bone_marrow", "Lung", "Kidney", "Heart"):
        G = {k: v for k, v in json.load(open(os.path.join(HERE, "%s_gold.json" % o.lower()))).items()
             if not k.startswith("_") and v}
        M = json.load(open(os.path.join(RES, "heca_to_cl_%s.json" % o)))["types"]
        for t, gd in G.items():
            v = M.get(t)
            if not v or v["n_cells"] < MINC or not v.get("cl"):
                continue
            X.append(features(v, ont))
            y.append(int(ok_gold(v["cl"][0]["curie"], gd)))
            organ.append(o)
    return np.array(X, float), np.array(y, int), np.array(organ)


def auc(s, t):
    o = np.argsort(s)
    r = np.empty(len(s))
    r[o] = np.arange(1, len(s) + 1)
    n1, n0 = t.sum(), len(t) - t.sum()
    return float((r[t == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)) if n1 and n0 else float("nan")


if __name__ == "__main__":
    X, y, organ, meta = build()
    print("calibration set: %d cell types, %d organs, base top-1 accuracy %.1f%%"
          % (len(y), len(set(organ)), 100 * y.mean()))
    print("                 top-5 accuracy %.1f%%\n"
          % (100 * np.mean([m["top5_ok"] for m in meta])))

    conf = oof(X, y, organ)
    print("out-of-fold (leave-one-organ-out)")
    for tgt in (0.95, 0.90, 0.85, 0.80):
        r = at_precision(conf, y, tgt)
        if r:
            print("  precision >= %.0f%%  ->  coverage %5.1f%%  (%3d of %d types,"
                  " threshold %.3f, realised precision %.1f%%)"
                  % (100 * tgt, 100 * r["coverage"], r["n_accepted"], len(y),
                     r["threshold"], 100 * r["precision"]))
        else:
            print("  precision >= %.0f%%  ->  not reachable at any threshold" % (100 * tgt))

    print("\nbaselines (single feature, same out-of-fold protocol)")
    for name, col in (("top score", 0), ("margin", 1), ("rel. margin", 2)):
        c1 = oof(X[:, [col]], y, organ)
        r = at_precision(c1, y, 0.90)
        print("  %-12s precision >= 90%% -> coverage %s"
              % (name, "%5.1f%%" % (100 * r["coverage"]) if r else "unreachable"))

    mu, sd = X.mean(0), X.std(0) + 1e-9
    w = fit((X - mu) / sd, y)
    print("\nfitted weights (standardised, full data — for inspection only)")
    for n, wi in sorted(zip(FEATS, w[1:]), key=lambda r: -abs(r[1])):
        print("  %-12s %+.3f" % (n, wi))

    # the docstring quotes an AUC for score-shape features alone; compute it rather than
    # assert it, since the ontology-shape lift is the whole point of reporting both
    n_score = len([f for f in FEATS if not f.startswith(("pred_", "anchor_", "n_anchors"))])
    auc_score_only = auc(oof(X[:, :n_score], y, organ), y)
    print("\nAUC against the crosswalks: %.3f on the %d score-shape features alone, "
          "%.3f with ontology shape" % (auc_score_only, n_score, auc(conf, y)))

    Xg, yg, og = gold_organ_set()
    cg = oof(Xg, yg, og)
    print("\nnon-circular check — hand-curated gold organs (%d types, %d organs)" % (len(yg), len(set(og))))
    print("  top-1 %.1f%%, leave-one-organ-out AUC %.3f" % (100 * yg.mean(), auc(cg, yg)))
    for tgt in (0.95, 0.90):
        r = at_precision(cg, yg, tgt)
        print("  precision >= %.0f%% -> %s" % (100 * tgt,
              "coverage %.1f%% (n=%d) — not deployable" % (100 * r["coverage"], r["n_accepted"])
              if r else "unreachable"))

    json.dump({"verdict": "confidence not learnable from mapper output; see module docstring",
               "gold_per_type": [{"organ": str(o), "conf": float(c), "correct": int(v)}
                                 for o, c, v in zip(og, cg, yg)],
               "auc_crosswalks": auc(conf, y), "auc_crosswalks_score_only": auc_score_only,
               "auc_gold_organs": auc(cg, yg), "n_gold": len(yg),
               "n_gold_organs": len(set(og)),
               "features": FEATS, "n": len(y), "base_top1": float(y.mean()),
               "operating_points": {str(t): at_precision(conf, y, t) for t in (0.95, 0.9, 0.85, 0.8)},
               "per_type": [dict(m, conf=float(c), correct=int(v))
                            for m, c, v in zip(meta, conf, y)]},
              open(os.path.join(RES, "calibration.json"), "w"), indent=1)
    print("\nwrote results/calibration.json")
