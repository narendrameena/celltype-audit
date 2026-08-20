#!/usr/bin/env python3
"""Stage 6: emit the annotation as a record, not a string.

A cell type annotated as the bare text 'Fibroblast' cannot be checked, merged across
atlases, or corrected later -- there is nothing to disagree with. This writes, for every
well-powered cell type, what was assigned, how confident that is, what evidence supports
it, what disputes it, and what could not be decided at all.

THE ASSIGNER IS LEXICAL, NOT EXPRESSION. Stage 2 established that the expression mapper
cannot be calibrated into an assigner (leave-one-organ-out AUC 0.563 on hand-curated
gold; precision caps near 73%%). Its role here is the AUDITOR: it supplies evidence and
raises contradictions, and never sets `assignment.curie`. Where the label cannot be
resolved the record ABSTAINS and carries the ranked expression shortlist for a curator,
because an abstention a human can act on is worth more than a guess they cannot see.

Fields that matter downstream:
  assignment.source     lexical | lexical-generalised | abstained
  assignment.exact      false when only a superclass was recoverable -- such records are
                        valid for lineage checks but MUST NOT be used for refinement
  audit.contradicted    the label's lineage is disjoint from all expression evidence
  refinement.confident  the split is discriminated and the cluster is not a mixture
  provenance            CL release and CELLxGENE snapshot, so a record can be re-derived

Usage: python benchmark/writeback.py
"""
import glob
import json
import math
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")

from cl_lineage import load, anchor_set                               # noqa: E402
from cl_resolve import resolve as resolve2                            # noqa: E402
from prune import asctb_sets, ASCTB_ORGAN, rerank                     # noqa: E402

MINC = 500
SCHEMA = "cellscribe-annotation/1"
EXACT = ("exact", "normalised", "alias", "general", "organ")


def provenance():
    g = json.load(open(os.path.join(os.path.dirname(HERE), "cl-full.json")))["graphs"][0]
    ver = [x.get("val") for x in (g.get("meta") or {}).get("basicPropertyValues", [])
           if "version" in str(x.get("pred", "")).lower()]
    snap = None
    for p in sorted(glob.glob(os.path.join(HERE, ".wmg_cache", "*.json")),
                    key=os.path.getsize, reverse=True)[:1]:
        snap = json.load(open(p)).get("snapshot_id")
    return {"cell_ontology_release": (ver or [None])[0],
            "cellxgene_wmg_snapshot": snap,
            "assigner": "cl_resolve.py (lexical)",
            "auditor": "heca_to_cl.py + prune.py (expression)",
            "refiner": "refine.py",
            "note": "expression evidence never sets assignment.curie; see stage 2"}


def main():
    g = load()
    ab = asctb_sets()
    ref = json.load(open(os.path.join(RES, "refinements.json")))
    byref = {(r["organ"], r["label"]): r for r in ref["proposals"]}

    out, counts = [], Counter()
    for p in sorted(glob.glob(os.path.join(RES, "heca_to_cl_*.json"))):
        d = json.load(open(p))
        organ, ub = d["organ"], d["uberon"]
        ctx = {c["curie"] for v in d["types"].values() for c in v.get("cl", [])}
        for t, v in sorted(d["types"].items(), key=lambda z: -z[1]["n_cells"]):
            if v["n_cells"] < MINC or not v.get("cl"):
                continue
            cur, how = resolve2(t, ctx, organ=organ)
            exact = how in EXACT
            A = anchor_set(cur) if cur else set()
            anc = [anchor_set(c["curie"]) for c in v["cl"][:5]]
            contradicted = bool(A) and any(anc) and all(B and not (A & B) for B in anc)

            # stage 3 re-rank: reference-support prior + organ prior; the lineage
            # constraint is withheld from contradicted types, which it would entrench
            boost = ab.get(ASCTB_ORGAN.get(organ, ""), set())
            ka = A if (A and not contradicted) else None
            scored = rerank(v["cl"], boost, ka, damp=1.0, with_scores=True)
            ranked = [c for _, c in scored]

            if not cur:
                src, counts["abstained"] = "abstained", counts["abstained"] + 1
            elif exact:
                src, counts["lexical"] = "lexical", counts["lexical"] + 1
            else:
                src, counts["lexical-generalised"] = ("lexical-generalised",
                                                      counts["lexical-generalised"] + 1)
            counts["contradicted"] += contradicted
            rp = byref.get((organ, t))

            out.append({
                "organ": organ, "uberon": ub, "atlas_label": t, "n_cells": v["n_cells"],
                "assignment": {
                    "curie": cur, "label": g["label"].get(cur) if cur else None,
                    "source": src, "resolution_status": how, "exact": bool(cur and exact),
                },
                "abstained": cur is None,
                "abstain_reason": None if cur else how,
                "evidence": {
                    "markers": v.get("markers", []),
                    # `score` is the raw expression statistic; `adjusted_score` is what
                    # the ranking uses (stage 3: reference-support prior + organ prior,
                    # and the lineage constraint where it was not withheld)
                    "expression_candidates": [
                        {"rank": i + 1, "curie": c["curie"], "label": c["label"],
                         "score": c["score"], "adjusted_score": round(adj, 4),
                         "reference_cells": c.get("n")}
                        for i, (adj, c) in enumerate(scored[:5])],
                    "expression_agrees_with_label": bool(
                        cur and ranked and (anchor_set(ranked[0]["curie"]) & A)),
                },
                "audit": {
                    "contradicted": contradicted,
                    "asserted_lineage": sorted(A),
                    "evidence_lineage": sorted(set().union(*[B for B in anc if B]) or set()),
                },
                "refinement": None if not rp else {
                    "curie": rp["to_curie"], "label": rp["to_label"],
                    "margin": rp["margin"], "confident": rp["confident"],
                    "n_evidenced_subtypes": rp["n_evidenced_subtypes"],
                },
            })

    doc = {"schema": SCHEMA, "provenance": provenance(),
           "counts": {k: int(v) for k, v in counts.items()},
           "n_records": len(out), "annotations": out}
    json.dump(doc, open(os.path.join(RES, "annotations.json"), "w"), indent=1)

    tsv = os.path.join(RES, "annotations.tsv")
    with open(tsv, "w") as fh:
        fh.write("organ\tatlas_label\tn_cells\tcl_curie\tcl_label\tsource\texact\t"
                 "abstained\tcontradicted\ttop_expression_candidate\trefine_to\trefine_confident\n")
        for r in out:
            a, e = r["assignment"], r["evidence"]["expression_candidates"]
            rf = r["refinement"] or {}
            fh.write("\t".join(str(x) for x in [
                r["organ"], r["atlas_label"], r["n_cells"], a["curie"] or "", a["label"] or "",
                a["source"], a["exact"], r["abstained"], r["audit"]["contradicted"],
                e[0]["label"] if e else "", rf.get("label", ""), rf.get("confident", "")]) + "\n")

    n = len(out)
    print("STAGE 6 — annotation records written\n")
    print("  records (cell types >= %d cells) : %d" % (MINC, n))
    for k in ("lexical", "lexical-generalised", "abstained"):
        print("    %-22s %4d (%.1f%%)" % (k, counts[k], 100 * counts[k] / n))
    print("    %-22s %4d" % ("contradicted (audit)", counts["contradicted"]))
    nr = sum(1 for r in out if r["refinement"])
    nc = sum(1 for r in out if (r["refinement"] or {}).get("confident"))
    print("    %-22s %4d (%d confident)" % ("with a refinement", nr, nc))
    print("\n  provenance: CL %s, CELLxGENE WMG snapshot %s"
          % (doc["provenance"]["cell_ontology_release"], doc["provenance"]["cellxgene_wmg_snapshot"]))
    print("  wrote results/annotations.json and results/annotations.tsv")


if __name__ == "__main__":
    main()
