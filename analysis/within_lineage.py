#!/usr/bin/env python3
"""Catch annotation errors that stay INSIDE one lineage.

The contradiction sweep asks whether the label's lineage anchor set is disjoint from the
evidence's. That question cannot be asked of "neutrophilic granulocyte" whose markers say
classical monocyte -- both are haematopoietic -- so 12 of the 18 known errors in the
hand-curated set are invisible to it, capping recall at 22%.

This asks a different question, of the same data: HOW WELL DOES THE TERM THE LABEL ASSERTS
ACTUALLY EXPLAIN THIS CLUSTER? Every CL term in the tissue is scored on the cluster's own
markers, and the asserted term is located in that ranking. A correctly labelled cluster
puts its own term at or near the top; a mislabelled one does not, whatever the lineages
involved. No ontology structure is used, which is precisely why it sees what the anchor
test cannot.

Ranking the asserted term turned out to be the WRONG question. "Neutrophilic granulocyte"
whose markers are FCN1/CD14/S100A12 still puts neutrophil at rank 7 of 247, because
neutrophils and monocytes share a great deal -- yet the top-scoring term is classical
monocyte, which is the curated answer. So the test asks who WINS, and whether the winner
is genuinely a different cell type:

  margin      score(asserted) / score(winner); low means the winner explains it better
  related     is the winner is_a-related to the asserted term? then the label is merely
              coarse or fine, not wrong, and nothing is flagged
  informative is the winner's absolute score above a floor? when every term scores badly
              the evidence is uninformative and the honest output is an abstention, not a
              flag -- this is what stops a liver neutrophil cluster being called a Schwann
              cell on noise

RESULT: THIS IS A TRIAGE QUEUE, NOT AN AUTOMATIC FLAG, and it is reported as one.

As a binary test it does not reach usable precision: leave-one-organ-out it tops out near
33% precision at 33% recall, because only about 12% of the candidate rows are real errors
and the score distributions overlap (errors median 0.42, non-errors 0.67). Shipping that
as an auto-flag would bury a curator in false positives.

Ranked, it works. Sorting every cell type by suspicion puts 3 of the top 5 and 8 of the
top 30 on real errors -- 8.7x and 3.9x the 6.9% base rate. Paired with the high-precision
anchor sweep it gives a two-tier system: the sweep auto-flags what it is sure of, the
queue is what a curator reads next. Together they surface 8 of 12 known errors while
reviewing 33 of 174 cell types.

The flagship catch is the one the anchor test structurally cannot make: lung
"neutrophilic granulocyte", 56,394 cells, whose top-scoring term is classical monocyte.

Usage: python benchmark/within_lineage.py
"""
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")

from scoring_variants import ok                                       # noqa: E402
from cl_lineage import load, anchor_set                               # noqa: E402
from cl_resolve import resolve as resolve2                            # noqa: E402

GOLD = ["Pancreas", "Liver", "Blood", "Bone_marrow", "Lung", "Kidney", "Heart"]
MINC, MIN_REF = 500, 100


def load_ref(stem=None):
    """stem: 'wide_ref' is the cache-assembled matrix; 'wide_ref_repro' is refetched from
    the API for exactly these tissues and genes and is ~2x denser (68% vs 35% nonzero,
    4157 vs 3450 CL terms). compare_reference.py shows the two agree perfectly where both
    have data (r=1.000), so the cache was incomplete rather than wrong."""
    stem = stem or os.environ.get("WIDE_REF", "wide_ref_repro")
    idx = json.load(open(os.path.join(RES, stem + "_index.json")))
    npz = np.load(os.path.join(RES, stem + ".npz"))
    mats = {k[3:]: npz[k] for k in npz.files if k.startswith("M__")}
    return idx, mats


def _related(a, b, hops=3):
    """is_a-related within a few hops, in either direction."""
    from cl_lineage import ancestors
    if a == b:
        return True
    return b in ancestors(a) or a in ancestors(b)


def features(idx, mats, organ, markers, asserted, min_ref=None):
    """-> (percentile, ratio, rank, n_terms, best_term, best_score, related) or None.

    min_ref overrides the module default; sensitivity.py sweeps it, and must reach this
    code path rather than re-implementing it, or the sweep measures a different pipeline.
    """
    min_ref = MIN_REF if min_ref is None else min_ref
    o = idx["organs"].get(organ)
    if not o:
        return None
    ub = o["uberon"]
    M, gix, tix = mats.get(ub), idx["gene_ix"].get(ub), idx["term_ix"].get(ub)
    if M is None or not gix or not tix or asserted not in tix:
        return None
    cols = [gix[g] for g in markers if g in gix]
    if not cols:
        return None
    s = M[:, cols].mean(axis=1)
    cnt = idx["counts"].get(ub, {})
    inv = {v: k for k, v in tix.items()}
    keep = np.array([cnt.get(inv[i], 0) >= min_ref for i in range(len(s))])
    if not keep.any() or not keep[tix[asserted]]:
        return None
    sk = np.where(keep, s, -1.0)
    order = np.argsort(-sk)
    ranks = np.empty(len(sk), dtype=int)
    ranks[order] = np.arange(len(sk))
    n = int(keep.sum())
    r = int(ranks[tix[asserted]])
    best = float(sk.max())
    bc = inv[int(order[0])]
    return (1.0 - r / max(n - 1, 1), float(s[tix[asserted]]) / (best + 1e-9),
            r, n, bc, best, _related(asserted, bc))


def build_rows(idx, mats, minc=None, min_ref=None):
    """One row per gold-annotated cluster that can be scored against the reference."""
    minc = MINC if minc is None else minc
    g = load()["label"]
    rows = []
    for o in GOLD:
        gp = os.path.join(HERE, "%s_gold.json" % o.lower())
        dp = os.path.join(RES, "heca_markers_deep_%s.json" % o)
        mp = os.path.join(RES, "heca_to_cl_%s.json" % o)
        if not all(os.path.exists(x) for x in (gp, dp, mp)):
            continue
        G = {k: v for k, v in json.load(open(gp)).items() if not k.startswith("_") and v}
        D = json.load(open(dp))["types"]
        M = json.load(open(mp))["types"]
        ctx = {c["curie"] for v in M.values() for c in v.get("cl", [])}
        for t, gd in G.items():
            v = D.get(t)
            if not v or v["n_cells"] < minc:
                continue
            cur, how = resolve2(t, ctx, organ=o)
            if not cur:
                continue
            f = features(idx, mats, o, [m["gene"] for m in v["markers"]], cur,
                         min_ref=min_ref)
            if not f:
                continue
            pct, ratio, r, n, best, bscore, related = f
            A, GA = anchor_set(cur), anchor_set(gd)
            rows.append({"organ": o, "label": t, "n_cells": v["n_cells"],
                         "asserted": cur, "asserted_name": g.get(cur, cur),
                         "gold": gd, "gold_name": g.get(gd, gd),
                         "percentile": round(pct, 4), "ratio": round(ratio, 4),
                         "rank": r, "n_terms": n, "best_term": g.get(best, best),
                         "best_curie": best, "best_score": round(float(bscore), 4),
                         "related_to_best": bool(related),
                         "error": not ok(cur, gd),
                         "crosses_lineage": bool(A and GA and not (A & GA))})
    return rows


def _call(r, margin_thr, floor):
    """A flag needs all three: a better winner, a DIFFERENT cell type, real evidence."""
    return (r["ratio"] < margin_thr and not r["related_to_best"]
            and r["best_score"] >= floor)


def evaluate(rows):
    """Leave-one-organ-out grid over (margin, evidence floor)."""
    floors = sorted({round(float(f), 4) for f in np.quantile(
        [r["best_score"] for r in rows], [0.0, 0.1, 0.2, 0.3, 0.4])})
    out = []
    for held in sorted({r["organ"] for r in rows}):
        tr = [r for r in rows if r["organ"] != held]
        bestf1, par = -1.0, (0.8, floors[0])
        for m in np.arange(0.5, 1.001, 0.05):
            for fl in floors:
                tp = sum(1 for r in tr if _call(r, m, fl) and r["error"])
                fp = sum(1 for r in tr if _call(r, m, fl) and not r["error"])
                fn = sum(1 for r in tr if not _call(r, m, fl) and r["error"])
                prec, rec = tp / max(tp + fp, 1), tp / max(tp + fn, 1)
                f1 = 2 * prec * rec / max(prec + rec, 1e-9)
                if prec >= 0.5 and f1 > bestf1:
                    bestf1, par = f1, (float(m), float(fl))
        for r in [x for x in rows if x["organ"] == held]:
            out.append(dict(r, flagged=_call(r, *par), margin_thr=par[0], floor=par[1]))
    return out


def queue(rows):
    """Most-suspicious first: an is_a-unrelated winner outranks a related one, then by
    how much better that winner explains the cluster."""
    return sorted(rows, key=lambda r: (r["related_to_best"], r["ratio"]))


if __name__ == "__main__":
    idx, mats = load_ref()
    rows = build_rows(idx, mats)
    ev = evaluate(rows)
    err = [r for r in ev if r["error"]]
    fl = [r for r in ev if r["flagged"]]
    tp = [r for r in fl if r["error"]]
    within = [r for r in err if not r["crosses_lineage"]]
    wcaught = [r for r in within if r["flagged"]]
    print("WITHIN-LINEAGE MARKER TEST — %d hand-curated types, %d organs\n" % (len(ev), len({r['organ'] for r in ev})))
    print("  known label errors        : %d (%d within-lineage)" % (len(err), len(within)))
    print("  flags raised              : %d" % len(fl))
    print("  precision                 : %d/%d = %.0f%%" % (len(tp), len(fl), 100 * len(tp) / max(len(fl), 1)))
    print("  recall over ALL errors    : %d/%d = %.0f%%" % (len(tp), len(err), 100 * len(tp) / max(len(err), 1)))
    print("  recall, within-lineage    : %d/%d = %.0f%%" % (len(wcaught), len(within), 100 * len(wcaught) / max(len(within), 1)))
    print("  (thresholds chosen leave-one-organ-out; a flag needs a better winner, an\n"
          "   is_a-UNRELATED winner, and evidence above the floor)")
    print("\n  errors caught that the anchor sweep cannot see:")
    for r in sorted(wcaught, key=lambda x: -x["n_cells"])[:10]:
        print("     %-8s %-27s %7d  asserted %-24s rank %d/%d  best=%s"
              % (r["organ"][:8], r["label"][:27], r["n_cells"], r["asserted_name"][:24],
                 r["rank"], r["n_terms"], r["best_term"][:24]))
    q = queue(rows)
    E = [r for r in rows if r["error"]]
    base = len(E) / max(len(rows), 1)
    print("\n  AS A RANKED REVIEW QUEUE (base rate %.1f%%)" % (100 * base))
    print("     %5s %8s %12s %11s" % ("depth", "errors", "precision@k", "enrichment"))
    for k in (5, 10, 20, 30):
        tp = sum(1 for r in q[:k] if r["error"])
        print("     %5d %8d %11.0f%% %10.1fx" % (k, tp, 100 * tp / k, (tp / k) / base))
    try:
        REC = {(r["organ"], r["label"]): r
               for r in json.load(open(os.path.join(RES, "auditor_recall.json")))}
        # A curator has a REVIEW BUDGET, so the budget is what is held fixed: anchor
        # flags are read first, then the queue, up to BUDGET items in total. Fixing the
        # queue depth instead and adding anchors on top would quietly vary the budget
        # with the number of flags, and would not match how the figure reports it.
        BUDGET = 33
        anchor = [k for k in [(r["organ"], r["label"]) for r in q]
                  if REC.get(k, {}).get("flagged")]
        seen = []
        for k in anchor + [(r["organ"], r["label"]) for r in q]:
            if k not in seen:
                seen.append(k)
        reviewed = set(seen[:BUDGET])
        keys = {(r["organ"], r["label"]) for r in E}
        tp = len(reviewed & keys)
        print("\n  TWO-TIER (anchor sweep read first, then the queue, %d reviewed in total)" % BUDGET)
        print("     reviewed %d of %d cell types (%.0f%%), found %d of %d errors = %.0f%% recall"
              % (len(reviewed), len(rows), 100 * len(reviewed) / len(rows), tp, len(E),
                 100 * tp / max(len(E), 1)))
    except Exception:
        pass
    json.dump({"rows": ev, "queue": [{"rank": i + 1, **r} for i, r in enumerate(q)]},
              open(os.path.join(RES, "within_lineage.json"), "w"), indent=1)
    print("\nwrote results/within_lineage.json")
