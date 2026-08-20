#!/usr/bin/env python3
"""Score expression-based CL mapping for one organ against a hand-curated gold standard.

Reports top-1 / top-5 accuracy, and — crucially — the ACHIEVABLE CEILING: a gold term that is
absent from the CELLxGENE reference for that tissue can never be returned, so accuracy must be
read against that ceiling rather than against 100%.

Verdicts: exact | ancestor (gold's superclass, <=3 is_a steps) | descendant (<=3 steps) | wrong.
"accepted" = anything but wrong.

Usage: python benchmark/score_organ.py <Organ> [<Organ> ...]
       expects benchmark/<organ>_gold.json and results/heca_to_cl_<Organ>.json
"""
import glob
import json
import os
import sys
from collections import defaultdict, deque

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(HERE, "results")
SHORT = lambda i: i.rsplit("/", 1)[-1].replace("_", ":")
ORDER = ["exact", "descendant", "ancestor", "wrong"]

g = json.load(open(os.environ.get("CL_JSON", os.path.join(ROOT, "cl-full.json"))))["graphs"][0]
LAB = {SHORT(n["id"]): n["lbl"] for n in g["nodes"] if n.get("lbl")}
PAR, CH = defaultdict(set), defaultdict(set)
for e in g["edges"]:
    if e.get("pred") == "is_a":
        PAR[SHORT(e["sub"])].add(SHORT(e["obj"]))
        CH[SHORT(e["obj"])].add(SHORT(e["sub"]))


def within(c, adj, depth=3):
    out, frontier = set(), {c}
    for _ in range(depth):
        nxt = set()
        for x in frontier:
            nxt |= adj.get(x, set())
        nxt -= out
        out |= nxt
        frontier = nxt
    return out


def verdict(pred, gold):
    if pred == gold:
        return "exact"
    if pred in within(gold, PAR):
        return "ancestor"
    if pred in within(gold, CH):
        return "descendant"
    return "wrong"


_REFCOUNTS = None


def reference_counts(uberon):
    """total_count per CL term in the CELLxGENE reference for this tissue.

    Built in ONE pass over the ~12 GB WMG cache for every tissue at once and memoised.
    This used to glob the whole cache per organ, which is fine for one organ and
    quadratic misery for seven -- the same trap that made summarise_all.py and refine.py
    slow before they were fixed.
    """
    global _REFCOUNTS
    if _REFCOUNTS is None:
        best, _REFCOUNTS = {}, {}
        for p in glob.glob(os.path.join(HERE, ".wmg_cache", "*.json")):
            try:
                d = json.load(open(p))
            except Exception:
                continue
            for ub, ct in ((d.get("term_id_labels") or {}).get("cell_types") or {}).items():
                if ct and len(ct) > best.get(ub, 0):
                    best[ub] = len(ct)
                    _REFCOUNTS[ub] = {k: ((v.get("aggregated") or {}).get("total_count", 0))
                                      for k, v in ct.items()}
    return _REFCOUNTS.get(uberon, {})


UH = json.load(open(os.path.join(RES, "uhaf_parsed.json")))
UH_CHILD = defaultdict(set)
for a, b in UH["edges"]:
    UH_CHILD[a].add(b)


def has_annotated_child(label, present):
    """True if a MORE SPECIFIC uHAF label is separately annotated in the same organ.

    The exclusive binary score cannot find markers for such a parent: every discriminative gene
    is claimed by one of its children, leaving only residual/generic transcripts.
    """
    seen, q = set(), deque(UH_CHILD.get(label, ()))
    while q:
        c = q.popleft()
        if c in seen:
            continue
        seen.add(c)
        if c in present:
            return True
        q.extend(UH_CHILD.get(c, ()))
    return False


def score(organ, min_cells_strata=(0, 500)):
    gp = os.path.join(HERE, "%s_gold.json" % organ.lower())
    mp = os.path.join(RES, "heca_to_cl_%s.json" % organ)
    if not (os.path.exists(gp) and os.path.exists(mp)):
        print("  %s: missing gold or mapping" % organ)
        return None
    GOLD = {k: v for k, v in json.load(open(gp)).items() if not k.startswith("_")}
    M = json.load(open(mp))
    ref = reference_counts(M["uberon"])
    rows = []
    for t, v in M["types"].items():
        gold = GOLD.get(t)
        if not gold:
            continue
        top = v["cl"][0]["curie"]
        v1 = verdict(top, gold)
        v5 = min([verdict(c["curie"], gold) for c in v["cl"]], key=ORDER.index)
        rows.append({"label": t, "n": v["n_cells"], "pred": top,
                     "pred_label": v["cl"][0]["label"], "gold": gold,
                     "gold_label": LAB.get(gold, "?"), "v1": v1, "v5": v5,
                     "achievable": ref.get(gold, 0) >= 100,
                     "parent_of_annotated": has_annotated_child(t, set(M["types"]))})
    if not rows:
        print("  %s: no scorable types" % organ)
        return None
    unmatched = [t for t in GOLD if t not in M["types"] and GOLD[t]]
    ok = lambda r, k: r[k] != "wrong"
    print("\n=== %s === (tissue %s; %d gold labels, %d present in the data & mapped)"
          % (organ, M["tissue"], len(GOLD), len(rows)))
    print("%-40s %-7s %-12s %-30s %s" % ("uHAF label", "n", "predicted", "CL label", "verdict"))
    print("-" * 104)
    for r in sorted(rows, key=lambda x: -x["n"]):
        flag = "" if r["achievable"] else "  (gold ABSENT from reference)"
        print("%-40s %-7d %-12s %-30s %s%s"
              % (r["label"][:40], r["n"], r["pred"], r["pred_label"][:30],
                 r["v1"] if r["v1"] == r["v5"] else "%s->top5:%s" % (r["v1"], r["v5"]), flag))
    out = {"organ": organ, "n": len(rows)}
    print()
    npar = sum(r["parent_of_annotated"] for r in rows)
    for mc in min_cells_strata:
        for ach, leaf in ((False, False), (True, False), (True, True)):
            sub = [r for r in rows if r["n"] >= mc and (r["achievable"] if ach else True)
                   and (not r["parent_of_annotated"] if leaf else True)]
            if not sub:
                continue
            k = "n>=%d%s%s" % (mc, " & achievable" if ach else "", " & leaf" if leaf else "")
            t1 = 100.0 * sum(ok(r, "v1") for r in sub) / len(sub)
            t5 = 100.0 * sum(ok(r, "v5") for r in sub) / len(sub)
            out[k] = {"n": len(sub), "top1": round(t1, 1), "top5": round(t5, 1)}
            print("  %-26s n=%-4d top-1 %5.1f%%   top-5 %5.1f%%" % (k, len(sub), t1, t5))
    nach = sum(r["achievable"] for r in rows)
    print("  ceiling: %d/%d gold terms present in the CELLxGENE reference = %.0f%%"
          % (nach, len(rows), 100.0 * nach / len(rows)))
    print("  parent categories whose children are separately annotated: %d/%d" % (npar, len(rows)))
    out["ceiling"] = round(100.0 * nach / len(rows), 1)
    if unmatched:
        print("  gold labels not present in this organ's data: %d" % len(unmatched))
    return out


if __name__ == "__main__":
    organs = sys.argv[1:] or ["Pancreas"]
    allout = [o for o in (score(x) for x in organs) if o]
    if len(allout) > 1:
        print("\n=== SUMMARY ===")
        print("%-16s %5s %8s %8s %9s" % ("organ", "n", "top-1", "top-5", "ceiling"))
        for o in allout:
            k = "n>=500 & achievable & leaf"
            s = o.get(k) or o.get("n>=0 & achievable") or {}
            print("%-16s %5d %7.1f%% %7.1f%% %8.1f%%"
                  % (o["organ"], s.get("n", 0), s.get("top1", 0), s.get("top5", 0), o["ceiling"]))
    json.dump(allout, open(os.path.join(RES, "organ_scores.json"), "w"), indent=1)
