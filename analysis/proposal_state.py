#!/usr/bin/env python3
"""Stage 6: give every proposal a life, and re-check it when CL moves.

build_proposals.py emits what the audit found *today*. On its own that is a pipe: it says
nothing about whether a gap was ever closed, and re-running it silently replaces yesterday's
answer with today's. This adds the return path.

Each proposal gets a stable id from (kind, organ, label) and a status:

    drafted     the audit found the gap
    submitted   a term request was opened; the issue number is recorded
    resolved    the gap is GONE in a later CL release -- measured, by re-running the
                same check that produced it
    rejected    a curator declined it; recorded by hand with a reason
    withdrawn   the audit no longer finds the gap for a reason other than CL changing
                (e.g. the atlas was re-annotated), so the proposal is moot

One distinction is kept scrupulously. `resolved` means the gap closed and is measured.
It does NOT mean the proposal caused it -- CL gains terms constantly and a curator may have
added one independently. Attribution requires the issue outcome, which is why `accepted` is
a separate field set from the tracker rather than inferred from a release diff. Claiming
credit for every closure would be exactly the sort of unearned inference this project
exists to catch.

The re-check reads only the RELEASED ontology. Proposals never enter the reference that
judges them; that separation is what stops the loop agreeing with itself.

Usage:
    python benchmark/proposal_state.py            # sync + re-check against the current CL
    python benchmark/proposal_state.py --metrics  # just print the programme numbers
"""
import argparse
import glob
import hashlib
import json
import os
import sys
from datetime import date

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
DOCS = os.path.abspath(os.path.join(HERE, "..", "..", "celltype-audit", "docs"))
STATE = os.path.join(DOCS, "state.json")
SCHEMA = "celltype-audit/proposal-state/1"

from cl_lineage import load, ancestors                                 # noqa: E402
from cl_resolve import resolve as resolve2                             # noqa: E402
import cl_lineage as _cl                                               # noqa: E402

OPEN = ("drafted", "submitted")


def cl_release():
    """The version CL states about itself, so a diff between runs is meaningful."""
    d = json.load(open(_cl.CLJSON))
    for g in d.get("graphs", []):
        for b in (g.get("meta") or {}).get("basicPropertyValues", []) or []:
            if b.get("pred", "").endswith("versionInfo"):
                return b["val"]
    return "unknown"


def pid(p):
    """Stable across regenerations: the gap, not the run that found it."""
    key = "%s|%s|%s" % (p["kind"], p["organ"], p["label"])
    return hashlib.sha1(key.encode()).hexdigest()[:12]


def load_state():
    if os.path.exists(STATE):
        return json.load(open(STATE))
    return {"schema": SCHEMA, "first_run": str(date.today()), "runs": [], "proposals": {}}


def gap_is_closed(p):
    """Re-run the check that produced this proposal. True when CL now covers it.

    Returns None where closure cannot be decided from the ontology alone -- a marker
    sufficient condition needs the axioms parsed, not just the label graph, so it is
    reported as undecidable rather than guessed at.
    """
    if p["kind"] in ("new-term", "synonym"):
        try:
            d = json.load(open(os.path.join(RES, "heca_to_cl_%s.json" % p["organ"])))
        except FileNotFoundError:
            return None
        ctx = {c["curie"] for v in d["types"].values() for c in v.get("cl", [])}
        cur, _how = resolve2(p["label"], ctx, organ=p["organ"])
        return bool(cur)
    if p["kind"] == "missing-axiom":
        g = load()
        sibs = [c for c, n in g["label"].items()
                if "ventricular cardiomyocyte" in n.lower() or
                   "cardiac muscle cell of ventric" in n.lower()]
        for a in sibs:                       # closed once any is_a path joins the pair
            for b in sibs:
                if a != b and (b in ancestors(a) or a in ancestors(b)):
                    return True
        return False
    return None


def sync(state, proposals, rel):
    today = str(date.today())
    seen = set()
    added = 0
    for p in proposals:
        i = pid(p)
        seen.add(i)
        if i not in state["proposals"]:
            state["proposals"][i] = {
                "kind": p["kind"], "organ": p["organ"], "label": p["label"],
                "n_cells": p["n_cells"], "status": "drafted", "issue": None,
                "accepted": None, "first_seen": today, "first_seen_release": rel,
                "closed_in_release": None,
                "history": [{"date": today, "status": "drafted", "release": rel}]}
            added += 1
        else:
            state["proposals"][i]["n_cells"] = p["n_cells"]
    # a gap the audit no longer reports, for reasons other than CL moving
    gone = 0
    for i, s in state["proposals"].items():
        if i not in seen and s["status"] in OPEN:
            s["status"] = "withdrawn"
            s["history"].append({"date": today, "status": "withdrawn", "release": rel,
                                 "note": "the audit no longer reports this gap"})
            gone += 1
    return added, gone


def recheck(state, rel):
    """Stage 6 proper: has the released ontology closed any open gap?"""
    today = str(date.today())
    closed, undecidable = [], 0
    for i, s in state["proposals"].items():
        if s["status"] not in OPEN:
            continue
        v = gap_is_closed(s)
        if v is None:
            undecidable += 1
        elif v:
            s["status"] = "resolved"
            s["closed_in_release"] = rel
            s["history"].append({"date": today, "status": "resolved", "release": rel,
                                 "note": "gap closed in CL; attribution not implied"})
            closed.append(s)
    return closed, undecidable


def metrics(state, rel):
    st = {}
    for s in state["proposals"].values():
        st[s["status"]] = st.get(s["status"], 0) + 1
    sub = [s for s in state["proposals"].values() if s["issue"]]
    acc = [s for s in sub if s.get("accepted") is True]
    return {"cl_release": rel, "by_status": st,
            "open": sum(st.get(k, 0) for k in OPEN),
            "submitted": len(sub), "accepted": len(acc),
            # only meaningful once something has been submitted AND adjudicated
            "acceptance_rate": (round(100.0 * len(acc) / len(sub), 1) if sub else None)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics", action="store_true", help="print the numbers and exit")
    a = ap.parse_args()

    rel = cl_release()
    state = load_state()
    props = json.load(open(os.path.join(DOCS, "proposals.json")))["proposals"]

    if a.metrics:
        print(json.dumps(metrics(state, rel), indent=1))
        return

    added, gone = sync(state, props, rel)
    closed, undec = recheck(state, rel)
    m = metrics(state, rel)
    state["runs"].append({"date": str(date.today()), "release": rel,
                          "added": added, "closed": len(closed), "withdrawn": gone})
    state["latest"] = m
    json.dump(state, open(STATE, "w"), indent=1)

    print("CL release %s\n" % rel)
    print("  proposals tracked : %d" % len(state["proposals"]))
    print("  newly drafted     : %d" % added)
    print("  gaps closed       : %d" % len(closed))
    print("  withdrawn         : %d" % gone)
    print("  undecidable here  : %d (need the axioms parsed, not just the label graph)"
          % undec)
    for s in closed:
        print("     closed: %-13s %s" % (s["organ"], s["label"]))
    print("\n  by status: %s" % json.dumps(m["by_status"]))
    print("  submitted %d, accepted %d, acceptance rate %s"
          % (m["submitted"], m["accepted"],
             "%.0f%%" % m["acceptance_rate"] if m["acceptance_rate"] is not None
             else "not yet measurable"))
    print("\nwrote %s" % STATE)


if __name__ == "__main__":
    main()
