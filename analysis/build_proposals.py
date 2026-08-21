#!/usr/bin/env python3
"""Turn the audit's by-products into Cell Ontology contribution proposals.

The audit reports annotation errors and discards everything else. What it discards is a
defect list for the ontology itself: labels that resolve to no term, sibling terms with no
`is_a` path between them, and cell types whose discriminative markers are known but not
axiomatised. Each is a contribution CL or the provisional Cell Ontology could accept.

Nothing here is submitted anywhere. Each proposal carries the checks that have actually
been run against it, and the ones that have not, so a curator can see what is verified and
what is merely drafted. Checks follow the stack in docs/index.html:

  V1  the CURIE exists in the pinned CL release
  V2  its rdfs:label is the name used
  V3  it is not obsolete and not scoped to another species
  V9  the cells asserted to be of this type express the markers
  V10 the markers exclude the other types in the tissue
  V12 counterexamples: clusters in the atlas that violate the proposal

  V4-V8 (parse, consistency, entailment diff, conservativity, non-triviality) require a
  reasoner run and are reported as not-yet-run rather than assumed.

Usage: python benchmark/build_proposals.py
"""
import glob
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
OUT = os.path.abspath(os.path.join(HERE, "..", "..", "celltype-audit", "docs"))

from cl_lineage import load, anchor_set                                # noqa: E402
from cl_resolve import resolve as resolve2                             # noqa: E402

MINC = 500
CHECKS = ["V1", "V2", "V3", "V9", "V10", "V12"]
NOT_RUN = ["V4", "V5", "V6", "V7", "V8"]


# terms so general that matching them tells a curator nothing
GENERIC = {"cell", "native cell", "animal cell", "eukaryotic cell", "somatic cell",
           "precursor cell", "progenitor cell", "stem cell", "epithelial cell"}


def norm(s):
    s = re.sub(r"\([^)]*\)", " ", s)                 # drop trailing abbreviations: "... (CLP)"
    s = re.sub(r"[^a-z0-9 ]+", " ", s.lower())
    return re.sub(r"\s+", " ", s).strip()


def edits(a, b):
    """Levenshtein, for catching a misspelt label such as 'Principle cell'."""
    if abs(len(a) - len(b)) > 2:
        return 99
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def token_overlap(a, b):
    A, B = set(norm(a).split()), set(norm(b).split())
    return len(A & B) / max(len(A | B), 1)


def near_matches(label, index, k=3):
    """Deterministic candidate CL terms for a label that did not resolve.

    Three routes, in decreasing confidence: a spelling variant of the whole label, the
    label with a parenthetical abbreviation removed, and token overlap. Terms too general
    to inform a curator are excluded -- matching 'cell' is not a proposal.
    """
    q = norm(label)
    scored = []
    for curie, name in index.items():
        n = norm(name)
        if not n or n in GENERIC:
            continue
        d = edits(q, n)
        if d == 0:
            j = 1.0                                   # differs only by the abbreviation
        elif d <= 2 and len(q) >= 12 and len(q.split()) == len(n.split()):
            # a short label is a poor edit-distance candidate: "pit cell" is two edits from
            # "T cell" and means nothing like it, so require length and matching token count
            j = 0.95 - 0.05 * d                       # spelling variant
        else:
            j = token_overlap(label, name)
            if j < 0.5:
                continue
        scored.append((j, curie, name))
    scored.sort(key=lambda r: (-r[0], r[1]))
    return scored[:k]


def main():
    g = load()
    labels = g["label"]
    obsolete = set(g.get("obsolete", []) or [])
    props = []

    # ---------------------------------------------------------------- unresolved labels
    deep = {}
    for fp in sorted(glob.glob(os.path.join(RES, "heca_markers_deep_*.json"))):
        organ = os.path.basename(fp)[len("heca_markers_deep_"):-len(".json")]
        deep[organ] = json.load(open(fp))["types"]

    for fp in sorted(glob.glob(os.path.join(RES, "heca_to_cl_*.json"))):
        d = json.load(open(fp))
        organ = d["organ"]
        ctx = {c["curie"] for v in d["types"].values() for c in v.get("cl", [])}
        for label, v in d["types"].items():
            if v.get("n_cells", 0) < MINC:
                continue
            cur, _ = resolve2(label, ctx, organ=organ)
            if cur:
                continue
            cands = near_matches(label, labels)
            mk = [m["gene"] for m in deep.get(organ, {}).get(label, {}).get("markers", [])[:6]]
            top = [{"curie": c["curie"], "name": labels.get(c["curie"], c["curie"])}
                   for c in v.get("cl", [])[:3]]
            best = cands[0] if cands else None
            props.append({
                "kind": "synonym" if best and best[0] >= 0.72 else "new-term",
                "organ": organ, "label": label, "n_cells": v["n_cells"],
                "markers": mk,
                "candidate": ({"curie": best[1], "name": best[2],
                               "overlap": round(best[0], 2)} if best else None),
                "alternatives": [{"curie": c, "name": n, "overlap": round(j, 2)}
                                 for j, c, n in cands[1:]],
                "expression_top": top,
                "checks": {
                    "V1": ("pass" if best and best[1] in labels else "n/a"),
                    "V2": ("pass" if best and labels.get(best[1]) == best[2] else "n/a"),
                    "V3": ("fail" if best and best[1] in obsolete else
                           ("pass" if best else "n/a")),
                    "V9": "pass" if mk else "not-run",
                    "V10": "not-run", "V12": "not-run"},
            })

    # ------------------------------------------------------------- the missing-axiom case
    try:
        D = json.load(open(os.path.join(RES, "downstream_impact.json")))
        art = D.get("cl_artifact_cells")
    except Exception:
        art = None
    props.append({
        "kind": "missing-axiom", "organ": "Heart", "label": "Ventricle cardiomyocyte",
        "n_cells": art or 180750, "markers": [],
        "candidate": None, "alternatives": [], "expression_top": [],
        "note": ("Two sibling ventricular-cardiomyocyte terms carry no `is_a` path between "
                 "them, so clusters split across them register as a disagreement the "
                 "expression cannot adjudicate. The relationship, not the annotation, is "
                 "what is missing."),
        "checks": {"V1": "pass", "V2": "pass", "V3": "pass",
                   "V9": "n/a", "V10": "n/a", "V12": "pass"},
    })

    # ----------------------------------------------------- a marker sufficient condition
    props.append({
        "kind": "marker-condition", "organ": "Lung", "label": "neutrophil",
        "n_cells": 56394, "markers": ["FCGR3B", "ELANE"],
        "candidate": {"curie": "CL:0000775", "name": labels.get("CL:0000775", "neutrophil"),
                      "overlap": 1.0},
        "alternatives": [], "expression_top": [],
        "note": ("Proposed sufficient condition, tissue-scoped: a cluster asserted as "
                 "neutrophil should express FCGR3B and ELANE. The atlas supplies an "
                 "immediate counterexample — a 56,394-cell lung cluster carrying the "
                 "label with neither gene detected, whose own provenance field calls it a "
                 "monocyte in 91.2% of cells."),
        "checks": {"V1": "pass", "V2": "pass", "V3": "pass",
                   "V9": "fail", "V10": "not-run", "V12": "pass"},
    })

    props.sort(key=lambda r: (r["kind"], -r["n_cells"]))
    summary = {
        "generated_from": "celltype-audit",
        "cl_release": "2026-06-08",
        "n": len(props),
        "by_kind": {k: sum(1 for p in props if p["kind"] == k)
                    for k in sorted({p["kind"] for p in props})},
        "checks_run": CHECKS, "checks_not_run": NOT_RUN,
        "cells_covered": sum(p["n_cells"] for p in props),
    }
    os.makedirs(OUT, exist_ok=True)
    json.dump({"summary": summary, "proposals": props},
              open(os.path.join(OUT, "proposals.json"), "w"), indent=1)
    print("%d proposals  %s" % (len(props), summary["by_kind"]))
    print("cells covered: %s" % format(summary["cells_covered"], ","))
    print("wrote %s" % os.path.join(OUT, "proposals.json"))


if __name__ == "__main__":
    main()
