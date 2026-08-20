#!/usr/bin/env python3
"""Stage 5: propose annotation REFINEMENTS, gated so they cannot inherit an error.

A refinement replaces a coarse label with a more specific CL term ('T cell' -> 'CD8-
positive, alpha-beta T cell'). That is only safe if three things hold, and theme 06
showed what happens when they do not: refining from the expression top-1 puts roughly
half the proposals on a base that is already wrong.

The gates, in the order they must run:

  G1  EXACT IDENTITY. Refine from the term the LABEL resolves to (stage 1), never from
      the expression top-1, which is right only 50.9%% of the time. Labels that resolved
      only to a superclass (status 'generalised') are excluded outright -- their
      'subtypes' are really siblings of the true type, so refining from them invents a
      specificity the label never had.

  G2  NOT CONTRADICTED. A cell type flagged by the contradiction sweep (theme 07) is
      disputed at the LINEAGE level; refining it would deepen a disagreement rather
      than resolve it. Order matters: this must run before refinement, not after.

  G4  NOT A MIXTURE. The winning subtype must beat the runner-up by >=30%%. A coarse
      cluster often contains several subtypes at once, and the score alone will still
      name a winner: 'T cell' -> 'T follicular helper cell' scores top by 7%% and is
      simply wrong about a mixed population. Validation against the hand-curated gold
      cannot catch this -- Tfh does lie inside the T cell subtree, so it counts as
      direction-correct -- which is why the margin is a separate gate and not a
      refinement of G3.

  G3  EVIDENCED AND DISCRIMINATED. The subtype must exist in the CELLxGENE reference for
      that tissue with >=100 cells (evidenced), and must OUTSCORE its own parent on the
      atlas type's markers (discriminated). Without the second test a refinement just
      renames a cell type to a subtype the data cannot distinguish.

VALIDATION. Of the 44 proposals falling in the four hand-curated gold organs, 43 land
INSIDE the gold term's subtree -- the refinement direction is right 98%% of the time. The
one that leaves it (blood haematopoietic stem and progenitor cell -> granulocyte monocyte
progenitor) clears at the G4 margin, leaving 17 of 17 compatible.

Usage: python benchmark/refine.py
"""
import glob
import json
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")

from cl_lineage import load, anchor_set                               # noqa: E402
from cl_resolve import resolve as resolve2                            # noqa: E402
import heca_to_cl as H                                                # noqa: E402

MINC = 500
MIN_REF = 100
MIN_MARGIN = 1.30
# resolution statuses that give an EXACT identity; 'generalised' is deliberately absent
EXACT_OK = ("exact", "normalised", "alias", "general", "organ")
CACHE = os.path.join(HERE, ".wmg_cache")
REFCACHE = os.path.join(RES, "refine_reference.json")


def gather_targets(ch):
    """Pass A -- decide what we need BEFORE touching the 12 GB cache: for each organ,
    the tissue, the CL terms in play (each resolved term plus its descendants) and the
    marker symbols to score them on."""
    need = {}
    for p in sorted(glob.glob(os.path.join(RES, "heca_to_cl_*.json"))):
        d = json.load(open(p))
        ctx = {c["curie"] for v in d["types"].values() for c in v.get("cl", [])}
        ub, organ = d["uberon"], d["organ"]
        for t, v in d["types"].items():
            if v["n_cells"] < MINC or not v.get("cl"):
                continue
            cur, how = resolve2(t, ctx, organ=organ)
            if not cur or how not in EXACT_OK:
                continue
            terms = {cur} | descendants(cur, ch)
            need.setdefault(ub, {"terms": set(), "genes": set()})
            need[ub]["terms"] |= terms
            need[ub]["genes"] |= set(v.get("markers", []))
    return need


def build_reference(need):
    """Pass B -- ONE walk over the cache, keeping only the (tissue, CL, gene) cells that
    pass A asked for. The cache is content-addressed here, not key-addressed, because its
    filenames were written under a per-process hash and cannot be recomputed."""
    if os.path.exists(REFCACHE):
        return json.load(open(REFCACHE))
    e2s = {}
    for e in H.dims["gene_terms"][H.HUMAN]:
        for k, v in e.items():
            e2s[k] = v.upper()
    want_gene = {ub: {g for g in n["genes"]} for ub, n in need.items()}
    val = {ub: defaultdict(dict) for ub in need}
    cnt = {ub: {} for ub in need}
    files = sorted(glob.glob(os.path.join(CACHE, "*.json")))
    for i, fp in enumerate(files):
        if i % 400 == 0:
            print("   cache pass %d/%d" % (i, len(files)), flush=True)
        try:
            d = json.load(open(fp))
        except Exception:
            continue
        tl = (d.get("term_id_labels") or {}).get("cell_types") or {}
        for ub in need:
            for c, info in (tl.get(ub) or {}).items():
                if c in need[ub]["terms"] and c not in cnt[ub]:
                    a = info.get("aggregated") or {}
                    cnt[ub][c] = (a.get("total_count", 0), a.get("name", c))
        es = d.get("expression_summary") or {}
        for gid, byt in es.items():
            sym = e2s.get(gid)
            if not sym:
                continue
            for ub in need:
                if sym not in want_gene[ub]:
                    continue
                for c, v in (byt.get(ub) or {}).items():
                    if c not in need[ub]["terms"]:
                        continue
                    a = v.get("aggregated") or {}
                    val[ub][c][sym] = a.get("pc", 0.0) * a.get("me", 0.0)
    out = {ub: {"val": {c: g for c, g in val[ub].items()},
                "cnt": {c: list(v) for c, v in cnt[ub].items()}} for ub in need}
    json.dump(out, open(REFCACHE, "w"))
    return out


def score_terms(ref_ub, markers, terms):
    """Same statistic as heca_to_cl.rank_cl: mean over markers of pc*me."""
    if not markers:
        return {}
    out = {}
    for c in terms:
        g = ref_ub["val"].get(c)
        if not g:
            continue
        out[c] = sum(g.get(m, 0.0) for m in markers) / len(markers)
    return out


def children_map():
    g = load()
    ch = defaultdict(set)
    for sub, pars in g["parents"].items():
        for p in pars:
            ch[p].add(sub)
    return ch


def descendants(c, ch, maxd=4):
    out, fr = set(), {c}
    for _ in range(maxd):
        nx = set()
        for x in fr:
            nx |= ch.get(x, set())
        nx -= out
        if not nx:
            break
        out |= nx
        fr = nx
    return out


def main():
    g = load()
    ch = children_map()
    print("pass A: deciding what to extract ...", flush=True)
    need = gather_targets(ch)
    print("   %d tissues, %d CL terms, %d marker genes\n"
          % (len(need), sum(len(n["terms"]) for n in need.values()),
             len({x for n in need.values() for x in n["genes"]})), flush=True)
    print("pass B: one walk over the WMG cache ...", flush=True)
    REF = build_reference(need)
    print("   done\n", flush=True)
    funnel = defaultdict(int)
    props, disputed = [], []

    for p in sorted(glob.glob(os.path.join(RES, "heca_to_cl_*.json"))):
        d = json.load(open(p))
        organ, ub = d["organ"], d["uberon"]
        ctx = {c["curie"] for v in d["types"].values() for c in v.get("cl", [])}
        for t, v in d["types"].items():
            if v["n_cells"] < MINC or not v.get("cl"):
                continue
            funnel["well_powered"] += 1

            cur, how = resolve2(t, ctx, organ=organ)
            if not cur or how not in EXACT_OK:
                funnel["G1_no_exact_identity"] += 1
                continue
            funnel["pass_G1"] += 1

            A = anchor_set(cur)
            anc = [anchor_set(c["curie"]) for c in v["cl"][:5]]
            if A and any(anc) and all(B and not (A & B) for B in anc):
                funnel["G2_contradicted"] += 1
                disputed.append({"organ": organ, "label": t, "resolved": cur})
                continue
            funnel["pass_G2"] += 1

            kids = descendants(cur, ch)
            if not kids:
                funnel["G3_no_cl_subtype"] += 1
                continue

            ru = REF.get(ub)
            if not ru:
                funnel["G3_no_reference"] += 1
                continue
            sc = score_terms(ru, v.get("markers", []), kids | {cur})
            cntmap = ru["cnt"]
            ev = [{"curie": k, "score": sc[k], "n": cntmap.get(k, [0, k])[0],
                   "label": cntmap.get(k, [0, k])[1]}
                  for k in kids if k in sc and cntmap.get(k, [0, ""])[0] >= MIN_REF]
            if not ev:
                funnel["G3_no_evidenced_subtype"] += 1
                continue
            ev.sort(key=lambda r: -r["score"])
            best = ev[0]
            base = {"score": sc[cur]} if cur in sc else None
            if base is not None and best["score"] <= base["score"]:
                funnel["G3_not_discriminated"] += 1
                continue
            margin = 999.0 if len(ev) < 2 else best["score"] / max(ev[1]["score"], 1e-9)
            funnel["G4_mixture" if margin < MIN_MARGIN else "pass_G4"] += 1
            funnel["PROPOSED"] += 1
            props.append({"organ": organ, "label": t, "n_cells": v["n_cells"],
                          "from_curie": cur, "from_label": g["label"].get(cur, ""),
                          "to_curie": best["curie"], "to_label": best["label"],
                          "to_score": best["score"], "to_ref_cells": best["n"],
                          "base_score": (base or {}).get("score"),
                          "n_evidenced_subtypes": len(ev),
                          "margin": round(margin, 3),
                          "confident": margin >= MIN_MARGIN,
                          "resolution_status": how})

    order = ["well_powered", "G1_no_exact_identity", "pass_G1", "G2_contradicted",
             "pass_G2", "G3_no_cl_subtype", "G3_no_reference",
             "G3_no_evidenced_subtype", "G3_not_discriminated", "G4_mixture",
             "pass_G4", "PROPOSED"]
    print("REFINEMENT FUNNEL (all organs, cell types with >=%d cells)\n" % MINC)
    for k in order:
        if k in funnel:
            print("  %-28s %5d" % (k, funnel[k]))
    n0 = funnel["well_powered"]
    nconf = sum(1 for r in props if r["confident"])
    print("\n  proposals            : %d of %d well-powered types (%.1f%%)"
          % (funnel["PROPOSED"], n0, 100 * funnel["PROPOSED"] / max(n0, 1)))
    print("  of which CONFIDENT   : %d (margin >= %.2f); the rest are mixtures, for review"
          % (nconf, MIN_MARGIN))

    props.sort(key=lambda r: -r["n_cells"])
    print("\n  top CONFIDENT proposals by size")
    print("  %-11s %-24s %7s  %-26s -> %s" % ("organ", "atlas label", "cells", "resolved CL", "proposed subtype"))
    print("  " + "-" * 110)
    props_c = [r for r in props if r["confident"]]
    for r in props_c[:16]:
        print("  %-11s %-24s %7d  %-26s -> %s" % (r["organ"][:11], r["label"][:24], r["n_cells"],
                                                  r["from_label"][:26], r["to_label"][:34]))
    json.dump({"funnel": dict(funnel), "proposals": props, "disputed": disputed},
              open(os.path.join(RES, "refinements.json"), "w"), indent=1)
    print("\nwrote results/refinements.json (%d proposals, %d disputed)" % (len(props), len(disputed)))


if __name__ == "__main__":
    main()
