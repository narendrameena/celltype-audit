#!/usr/bin/env python3
"""Ask CellTypist about every proposal, in the organs where it can be asked at all.

The head-to-head against the gold uses only organs the HRA crosswalk covers, because that
comparison needs CellTypist's answer as a CL term. A proposal needs less: what a curator
wants to know is whether an independent annotator, shown the same cells, names a population
CL already has. Its raw label carries that even where the crosswalk does not reach, so this
runs wider and is explicit about the two things that weakens.

  * CL mapping. For Lung and Heart the HRA CellTypist crosswalk supplies the CL term. For
    the rest there is no crosswalk level, so the label is resolved by this study's own
    lexical resolver -- the same route the paper shows is circular when used as a
    benchmark. Recorded as `cl_via` so the two are never confused.
  * Model fit. Some organs have only a fetal or developmental model. A first-trimester
    gonad model asked about adult spermatogenic cells is being misapplied, and saying so is
    more useful than a confident wrong label. Recorded as `fit`.

Writes results/proposal_celltypist.json.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
DOCS = os.path.abspath(os.path.join(HERE, "..", "..", "celltype-audit", "docs"))

import baselines as B                                       # noqa: E402
from cl_lineage import load                                 # noqa: E402
from cl_resolve import resolve as resolve2                  # noqa: E402

# organ -> (model, fit, why)
PLAN = {
    "Brain":       ("Adult_Human_MTG.pkl", "matched",
                    "adult human cortex; the proposals are cortical interneurons"),
    "Breast":      ("Cells_Adult_Breast.pkl", "matched", "adult human breast"),
    "Spinal_cord": ("Developing_Human_Brain.pkl", "matched",
                    "first-trimester brain; floor, roof and midplate ARE first-trimester "
                    "structures, so the developmental model is the right one here"),
    "Nose":        ("Cells_Lung_Airway.pkl", "approximate",
                    "airway model covering nasal locations, not a dedicated nasal model"),
    "Eye":         ("Fetal_Human_Retina.pkl", "stage-mismatched",
                    "only fetal retina models exist; the proposals are adult retina"),
    "Adrenal_gland": ("Fetal_Human_AdrenalGlands.pkl", "stage-mismatched",
                      "only a fetal adrenal model exists; the proposal is adult"),
    "Testis":      ("Developing_Human_Gonads.pkl", "stage-mismatched",
                    "first-trimester gonads, where spermatogenesis has not begun"),
}
NO_MODEL = {
    "Salivary_gland": "CellTypist ships no salivary gland model",
    "Ureter":         "CellTypist ships no ureter or urothelium model",
    "Stomach":        "CellTypist ships no stomach model (the intestinal model is not one)",
    "Pleura":         "CellTypist ships no pleura or mesothelium model",
}


def main():
    import anndata as ad
    import celltypist
    from celltypist import models

    g = load()
    L = g["label"]
    doc = json.load(open(os.path.join(DOCS, "proposals.json")))
    want = {}
    for p in doc["proposals"]:
        want.setdefault(p["organ"], []).append(p["label"])
    out = {}
    for organ, labels in sorted(want.items()):
        if organ in NO_MODEL:
            out[organ] = {"status": "no-model", "note": NO_MODEL[organ]}
            print("   %-15s no model -- %s" % (organ, NO_MODEL[organ]), flush=True)
            continue
        if organ not in PLAN:
            continue                                   # already covered by baselines.py
        mname, fit, why = PLAN[organ]
        path = os.path.abspath(B.H5 % organ)
        if not os.path.exists(path):
            out[organ] = {"status": "no-data", "note": "no hECA file"}
            continue
        print("=== %s  model=%s  (%s)" % (organ, mname, fit), flush=True)
        model = models.Model.load(model=mname)
        feats = list(model.classifier.features)
        X, lab, nfound = B.read_subset(path, feats)
        if X is None:
            out[organ] = {"status": "no-clusters", "note": "no well-powered clusters"}
            continue
        # CellTypist refuses a matrix whose log1p values exceed log1p(1e4), as a guard
        # against unnormalised input. hECA is CP10K throughout (row sums 9,216-10,000
        # across organs, checked directly), but a few individual cells carry one gene
        # above that budget and trip the guard for the whole run. Those cells are dropped
        # and counted rather than the guard being bypassed: a cell where a single gene
        # exceeds the entire CP10K total is aberrant on its own terms.
        bad = (X > 9.21) .any(axis=1)
        ndrop = int(bad.sum())
        if ndrop:
            X, lab = X[~bad], lab[~bad]
            print("   dropped %d of %d cells with a gene above 1e4 CP10K"
                  % (ndrop, ndrop + X.shape[0]), flush=True)
        a = ad.AnnData(X=X, obs={"atlas_label": lab.astype(str)},
                       var={"gene": np.array(feats, dtype=object)})
        a.var_names = [str(x) for x in feats]
        pred = celltypist.annotate(a, model=mname, majority_voting=True,
                                   over_clustering="atlas_label")
        pl = pred.predicted_labels
        col = "majority_voting" if "majority_voting" in pl.columns else "predicted_labels"
        from collections import Counter
        ctx = set()
        types = {}
        for t in sorted(set(lab)):
            m = lab == t
            call, k = Counter(pl[col].values[m]).most_common(1)[0]
            cur, how = resolve2(str(call), ctx, organ=organ)
            types[t] = {"celltypist": str(call), "vote_fraction": round(k / m.sum(), 3),
                        "cl": cur, "cl_label": L.get(cur, cur) if cur else None,
                        "cl_via": "lexical resolver (no HRA crosswalk level for this organ)",
                        "n_sampled": int(m.sum())}
        out[organ] = {"status": "scored", "model": mname, "fit": fit, "fit_note": why,
                      "cells_dropped_over_cp10k": ndrop, "types": types}
        for t in labels:
            r = types.get(t)
            print("   %-28s -> %s" % (t[:28], (r or {}).get("celltypist", "cluster not scored")),
                  flush=True)
        json.dump(out, open(os.path.join(RES, "proposal_celltypist.json"), "w"), indent=1)
    json.dump(out, open(os.path.join(RES, "proposal_celltypist.json"), "w"), indent=1)
    print("\nwrote results/proposal_celltypist.json")


if __name__ == "__main__":
    main()
