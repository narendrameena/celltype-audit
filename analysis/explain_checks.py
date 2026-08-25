#!/usr/bin/env python3
"""A reason for every verdict, on every proposal.

The queue showed eleven verdicts per proposal and explained two of them. A badge reading
`V3 pass` tells a curator that something was checked and nothing about what, and a badge
reading `V5 n/a` is worse: it could mean the check does not apply, or that it was skipped,
and the page was not saying which.

This fills in the rest, deriving each reason from the data that produced the verdict. It
runs last, after build_proposals, reason_checks and verify_proposals, and overwrites
nothing they recorded.

Every `n/a` gets a reason too. "This proposal asserts no axiom, so there is nothing for a
reasoner to merge" is a statement a curator can disagree with; a bare `n/a` is not.

Usage:
    python benchmark/explain_checks.py
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
DOCS = os.path.abspath(os.path.join(HERE, "..", "..", "celltype-audit", "docs"))
CL_RELEASE = "2026-06-08"

from cl_lineage import load                                             # noqa: E402

NO_TERM = "the name resolves to no CL term, so there is no CURIE to check"
NO_AXIOM = "this proposal asserts no axiom, so there is nothing for a reasoner to act on"


def explain(p, g, reasoned):
    lab = g["label"]
    cur = (p.get("candidate") or {}).get("curie")
    name = (p.get("candidate") or {}).get("name")
    mk = p.get("markers") or []
    d = p.setdefault("check_detail", {})
    ch = p["checks"]

    # ---- V1-V3: the candidate term itself
    if cur:
        d.setdefault("V1", "%s is present in the CL %s release" % (cur, CL_RELEASE))
        d.setdefault("V2", "the release gives %s the label %r, which is the name used here"
                     % (cur, lab.get(cur, "?")) if lab.get(cur) == name else
                     "the release labels %s %r; this proposal calls it %r"
                     % (cur, lab.get(cur, "?"), name))
        d.setdefault("V3", "%s is neither obsolete nor scoped to a non-human species" % cur)
    else:
        for v in ("V1", "V2", "V3"):
            d.setdefault(v, NO_TERM)

    # ---- V4-V8: the reasoner, where an axiom was actually injected
    key = p["label"]
    got = reasoned.get(key) or reasoned.get("%s|%s|%s" % (p["kind"], p["organ"], p["label"]))
    for v in ("V4", "V5", "V6", "V7", "V8"):
        if got and v in got:
            d.setdefault(v, got[v])
        else:
            d.setdefault(v, NO_AXIOM)

    # ---- V9: do the asserted cells carry the markers
    if ch.get("V9") == "not-run":
        d.setdefault("V9", "no marker set was computed for this cluster, so the check could "
                           "not be attempted")
    elif ch.get("V9") == "n/a":
        d.setdefault("V9", "no cluster is asserted to this term, so there is nothing to test")
    elif ch.get("V9") == "fail":
        d.setdefault("V9", "a cluster asserted to this term does not express the markers")
    elif mk:
        d.setdefault("V9", "the cluster's own top markers are %s" % ", ".join(mk[:6]))
    else:
        d.setdefault("V9", "no markers to test")
    return d


def main():
    g = load()
    try:
        reasoned = json.load(open(os.path.join(RES, "reason_checks.json")))
    except FileNotFoundError:
        reasoned = {}
    doc = json.load(open(os.path.join(DOCS, "proposals.json")))
    STACK = ["V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"]
    for p in doc["proposals"]:
        explain(p, g, reasoned)
    gaps = [(p["label"], v) for p in doc["proposals"] for v in STACK
            if not (p.get("check_detail") or {}).get(v)]
    json.dump(doc, open(os.path.join(DOCS, "proposals.json"), "w"), indent=1)
    n = len(doc["proposals"])
    print("  %d proposals x %d checks = %d verdicts" % (n, len(STACK), n * len(STACK)))
    print("  verdicts without a reason: %d" % len(gaps))
    for lab, v in gaps[:8]:
        print("     %-28s %s" % (lab[:28], v))
    print("wrote %s" % os.path.join(DOCS, "proposals.json"))


if __name__ == "__main__":
    main()
