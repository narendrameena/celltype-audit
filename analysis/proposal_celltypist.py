#!/usr/bin/env python3
"""Attach an independent annotator's opinion to every proposal, or say why there is none.

A proposal claims CL lacks a term for a population. A CL-native automated annotator is an
independent way to test that: if CellTypist confidently names the cluster with a CL term,
the gap is doubtful; if it cannot be run at all, the reader should be told that rather than
shown a blank.

The result is mostly the second case, and that is the finding. Of the thirteen organs where
this audit proposes a gap, the HRA CellTypist crosswalk covers two. Proposals concentrate
where the annotation tooling does not reach -- eye, salivary gland, ureter, spinal cord,
pleura, nose -- which is an argument for the method rather than an omission in it.

Follows the page's existing distinction: `n/a` is a check that does not apply, `not-run` is
one that applies and was not done. Here "no model" is neither -- it is a check that cannot
be run because the comparator does not cover this organ, and it is labelled as such.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
DOCS = os.path.abspath(os.path.join(HERE, "..", "..", "celltype-audit", "docs"))

from cl_lineage import load                                 # noqa: E402


def main():
    g = load()
    L = g["label"]
    CT = json.load(open(os.path.join(RES, "baseline_celltypist.json")))
    doc = json.load(open(os.path.join(DOCS, "proposals.json")))
    try:
        WIDE = json.load(open(os.path.join(RES, "proposal_celltypist.json")))
    except FileNotFoundError:
        WIDE = {}
    scored = 0
    for p in doc["proposals"]:
        organ, label = p["organ"], p["label"]
        # Where the proposal came from, stated on the proposal rather than only in the
        # page header: a reader looking at one card should not have to infer the atlas.
        p["source"] = {"atlas": "hECA v2.0", "organ": organ,
                       "file": "RNA-%s.h5ad" % organ,
                       "n_cells": p.get("n_cells"),
                       "note": ("hECA aggregates published datasets and rewrites their "
                                "annotations into one vocabulary, so the label under test "
                                "belongs to that harmonisation step, not to the groups who "
                                "generated the cells")}
        wide = WIDE.get(organ) or {}
        if wide.get("status") == "scored" and label in (wide.get("types") or {}):
            r = wide["types"][label]
            p["celltypist"] = {
                "status": "scored", "model": wide["model"], "fit": wide.get("fit"),
                "fit_note": wide.get("fit_note"), "cluster": label,
                "call": r["celltypist"], "vote_fraction": r["vote_fraction"],
                "cl": r.get("cl"), "cl_label": r.get("cl_label"), "cl_via": r.get("cl_via"),
                "note": ("CellTypist (%s) calls this cluster %r on %s of votes%s"
                         % (wide["model"].replace(".pkl", ""), r["celltypist"],
                            "{:.0%}".format(r["vote_fraction"]),
                            "; that label resolves to %s" % r["cl_label"] if r.get("cl")
                            else "; that label is study nomenclature with no CL term"))}
            scored += 1
            continue
        if wide.get("status") == "no-model":
            p["celltypist"] = {"status": "no-model", "note": wide["note"]}
            continue
        blob = CT.get(organ)
        if not blob:
            p["celltypist"] = {
                "status": "no-model",
                "note": ("CellTypist has no model for %s that the HRA crosswalk can "
                         "translate to CL, so it cannot be asked about this proposal" % organ)}
            continue
        types = blob.get("types") or {}
        rec = types.get(label)
        if rec is None and p["kind"] == "marker-condition":
            # the proposal is about a TERM; the cluster that bears on it is the one the
            # counterexample search found, not a cluster of the same name
            viol = [k for k in types if "eutrophil" in k or "ranulocyte" in k]
            rec = types.get(viol[0]) if viol else None
            label = viol[0] if viol else label
        if rec is None:
            p["celltypist"] = {
                "status": "cluster-not-scored",
                "model": blob.get("model"),
                "note": ("CellTypist ran on %s but did not score a cluster named %r"
                         % (organ, label))}
            continue
        cl = (rec.get("cl") or [None])[0]
        p["celltypist"] = {
            "status": "scored",
            "model": blob.get("model"),
            "cluster": label,
            "call": rec.get("celltypist"),
            "vote_fraction": rec.get("vote_fraction"),
            "cl": cl,
            "cl_label": L.get(cl, cl) if cl else None,
            "note": ("CellTypist calls this cluster %r (%s of votes)%s"
                     % (rec.get("celltypist"), "{:.0%}".format(rec.get("vote_fraction", 0)),
                        ", which the crosswalk maps to %s" % L.get(cl, cl) if cl else
                        ", which the crosswalk maps to no CL term"))}
        scored += 1
    json.dump(doc, open(os.path.join(DOCS, "proposals.json"), "w"), indent=1)
    from collections import Counter
    c = Counter(p["celltypist"]["status"] for p in doc["proposals"])
    organs = {p["organ"] for p in doc["proposals"]}
    covered = {p["organ"] for p in doc["proposals"] if p["celltypist"]["status"] != "no-model"}
    print("  proposals with a CellTypist call : %d of %d" % (scored, len(doc["proposals"])))
    print("  status: %s" % dict(c))
    print("  organs where the audit proposes a gap : %d" % len(organs))
    print("  of those, ones CellTypist covers      : %d (%s)"
          % (len(covered), ", ".join(sorted(covered))))
    for p in doc["proposals"]:
        if p["celltypist"]["status"] == "scored":
            print("     %-8s %-28s -> %s" % (p["organ"], p["label"][:28], p["celltypist"]["note"]))
    print("wrote %s" % os.path.join(DOCS, "proposals.json"))


if __name__ == "__main__":
    main()
