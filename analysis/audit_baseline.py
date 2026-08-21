#!/usr/bin/env python3
"""Does the auditor catch errors that an established annotator makes -- and does an
established annotator independently corroborate the auditor's flags?

Two questions, one table:

  (1) CellTypist is run on the same hECA cells (baselines.py) and its calls are translated
      to CL through the HRA CTann crosswalk. Where CellTypist's lineage is disjoint from
      the atlas label's, two independent annotators disagree about what the cells are.

  (2) For every cell type the contradiction sweep flagged, CellTypist casts an independent
      vote. If CellTypist's lineage matches the EVIDENCE lineage rather than the label's,
      the flag is CORROBORATED by a tool that never saw our pipeline. That is external
      confirmation of a mis-annotation, which no amount of internal validation can supply.

Usage: python benchmark/audit_baseline.py
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

MINC = 500
K = 5


def anchors_of(curies):
    a = set()
    for c in curies:
        a |= anchor_set(c)
    return a


def main():
    g = load()
    B = json.load(open(os.path.join(RES, "baseline_celltypist.json")))
    rows, flags = [], {}
    for p in sorted(glob.glob(os.path.join(RES, "heca_to_cl_*.json"))):
        d = json.load(open(p))
        organ = d["organ"]
        if organ not in B:
            continue
        ctx = {c["curie"] for v in d["types"].values() for c in v.get("cl", [])}
        for t, v in d["types"].items():
            if v["n_cells"] < MINC or len(v.get("cl", [])) < K:
                continue
            bt = B[organ]["types"].get(t)
            if not bt or not bt["cl"]:
                continue
            cur, how = resolve2(t, ctx, organ=organ)
            A = anchor_set(cur) if cur else set()
            E = anchors_of([c["curie"] for c in v["cl"][:K]])       # our expression evidence
            C = anchors_of(bt["cl"])                                 # CellTypist
            if not A or not C:
                continue
            contradicted = bool(E) and all(anchor_set(c["curie"]) and
                                           not (A & anchor_set(c["curie"])) for c in v["cl"][:K])
            rows.append({"organ": organ, "label": t, "n_cells": v["n_cells"],
                         "markers": v["markers"][:4],
                         "label_cl": cur, "label_anchors": sorted(A),
                         "celltypist": bt["celltypist"], "ct_cl": bt["cl"],
                         "ct_anchors": sorted(C), "evidence_anchors": sorted(E),
                         "ct_agrees_label": bool(A & C),
                         "ct_agrees_evidence": bool(E & C),
                         "flagged": contradicted})
    n = len(rows)
    agree = sum(r["ct_agrees_label"] for r in rows)
    fl = [r for r in rows if r["flagged"]]
    corro = [r for r in fl if r["ct_agrees_evidence"] and not r["ct_agrees_label"]]
    ct_dis = [r for r in rows if not r["ct_agrees_label"]]
    caught = [r for r in ct_dis if r["flagged"]]

    print("AUDITOR vs AN ESTABLISHED ANNOTATOR  (CellTypist, %d organs)\n" % len(B))
    print("  cell types with both a resolved label and a CellTypist call : %d" % n)
    print("  CellTypist lineage agrees with the atlas label              : %d (%.0f%%)"
          % (agree, 100 * agree / max(n, 1)))
    print("  CellTypist disagrees with the atlas label                   : %d" % len(ct_dis))
    print("\n  our contradiction flags in these organs                     : %d" % len(fl))
    print("  of those, CORROBORATED by CellTypist (it sides with the\n"
          "  expression evidence, not the label)                          : %d" % len(corro))

    if corro:
        print("\n  EXTERNALLY CONFIRMED MIS-ANNOTATIONS")
        print("  %-9s %-24s %7s %-26s %-22s %s" % ("organ", "atlas label", "cells",
                                                   "its markers", "CellTypist says", "lineage"))
        print("  " + "-" * 118)
        for r in sorted(corro, key=lambda x: -x["n_cells"]):
            print("  %-9s %-24s %7d %-26s %-22s %s→%s"
                  % (r["organ"][:9], r["label"][:24], r["n_cells"], ",".join(r["markers"])[:26],
                     r["celltypist"][:22], "/".join(r["label_anchors"])[:9],
                     "/".join(r["ct_anchors"])[:12]))

    # CellTypist AS A COMPETING DETECTOR. The obvious baseline for "find the mislabelled
    # clusters" is to run an established annotator and flag where it disagrees with the
    # label. Scoring that the same way the sweep is scored -- against the hand-curated gold,
    # in the organs where a gold exists -- is the comparison a referee will construct if we
    # do not, so compute it here rather than leaving it implicit.
    from scoring_variants import ok as ok_gold
    GOLD = {}
    for o in sorted({r["organ"] for r in rows}):
        gp = os.path.join(HERE, "%s_gold.json" % o.lower())
        if os.path.exists(gp):
            GOLD[o] = {k: v for k, v in json.load(open(gp)).items()
                       if not k.startswith("_") and v}
    def is_error(r):
        g = GOLD.get(r["organ"], {}).get(r["label"])
        return None if not g else (not ok_gold(r["label_cl"], g))
    ct_scored = [(r, is_error(r)) for r in ct_dis]
    ct_tp = sum(1 for _r, v in ct_scored if v is True)
    ct_fp = sum(1 for _r, v in ct_scored if v is False)
    ct_unk = sum(1 for _r, v in ct_scored if v is None)
    base = {"flags": len(ct_dis), "scoreable": ct_tp + ct_fp, "true": ct_tp,
            "false": ct_fp, "uncurated": ct_unk,
            "precision": (100.0 * ct_tp / (ct_tp + ct_fp)) if (ct_tp + ct_fp) else None}
    print("\n  CELLTYPIST AS A COMPETING DETECTOR (flag where it disagrees with the label)")
    print("     flags raised                     : %d of %d cell types" % (len(ct_dis), n))
    print("     falling in a curated organ       : %d (%d not curated, unscoreable)"
          % (base["scoreable"], ct_unk))
    print("     of those, real errors            : %d  -> precision %s"
          % (ct_tp, "%.0f%%" % base["precision"] if base["precision"] is not None else "n/a"))
    for r, v in sorted(ct_scored, key=lambda x: -x[0]["n_cells"]):
        print("       %-9s %-26s %-24s %s"
              % (r["organ"][:9], r["label"][:26], r["celltypist"][:24],
                 {True: "REAL ERROR", False: "false flag", None: "not curated"}[v]))

    print("\n  CellTypist calls the auditor disputes but we did NOT flag   : %d"
          % (len(ct_dis) - len(caught)))
    for r in sorted([x for x in ct_dis if not x["flagged"]], key=lambda x: -x["n_cells"])[:8]:
        print("     %-9s %-24s label=%-16s CellTypist=%s"
              % (r["organ"][:9], r["label"][:24], "/".join(r["label_anchors"])[:16], r["celltypist"][:24]))

    json.dump({"n": n, "agree": agree, "flags": len(fl), "corroborated": len(corro),
               "celltypist_as_detector": base, "rows": rows},
              open(os.path.join(RES, "audit_baseline.json"), "w"), indent=1)
    print("\nwrote results/audit_baseline.json")


if __name__ == "__main__":
    main()
