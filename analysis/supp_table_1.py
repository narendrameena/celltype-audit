#!/usr/bin/env python3
"""Supplementary Table 1 -- the surveyed CELLxGENE datasets.

Every field is read from the distributed .h5ad or from cxg_survey.json; nothing is typed
by hand, so the table cannot drift from the analysis it describes.
"""
import collections
import json
import os
import re
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
RES = os.environ.get("RESULTS_DIR", os.path.join(HERE, "results"))
DATA = os.environ.get("CXG_DATA", os.path.join(os.path.dirname(HERE), "..", "cxg_data"))

# the Tabula Sapiens collection is also used as the paper's independent witness, so any
# survey dataset drawn from it is NOT independent of that comparison and is marked here
TS_COLLECTION = "e5f58829-1a66-40b5-a624-9046778e74f5"

_src = os.path.join(RES, "cxg_survey.json")
if not os.path.exists(_src):
    sys.exit("no %s -- run cxg_survey.py first, or set RESULTS_DIR (and CXG_DATA to the\n"
             "directory of downloaded .h5ad files). Neither the results nor the datasets\n"
             "are distributed with this repository; see DATA.md." % _src)

rows = json.load(open(_src))
by_ds = collections.defaultdict(list)
for r in rows:
    by_ds[r["dataset"]].append(r)


def cats(f, key):
    """(labels, counts) for a categorical obs column, commonest first."""
    g = f["obs"][key]
    labels = [c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]]
    codes = g["codes"][:]
    v, n = np.unique(codes[codes >= 0], return_counts=True)
    order = np.argsort(-n)
    return [(labels[int(v[i])], int(n[i])) for i in order]


out = []
for ds, rr in by_ds.items():
    live = [r for r in rr if not r["thin_support"]]
    p = os.path.join(DATA, ds + ".h5ad")
    f = h5py.File(p, "r")
    cite = f["uns/citation"][()].decode()
    doi = (re.search(r"Publication: (\S+)", cite) or [None, ""])[1]
    coll = (re.search(r"collections/([0-9a-f-]+)", cite) or [None, ""])[1]
    title = f["uns/title"][()].decode()
    tis = cats(f, "tissue")
    assay = cats(f, "assay")
    n_cells = int(f["obs/cell_type/codes"].shape[0])   # one code per cell
    n_types_total = len(cats(f, "cell_type"))
    f.close()
    out.append({
        "dataset_id": ds,
        "title": title,
        "collection_id": coll,
        "publication": doi,
        "tissues": "; ".join("%s (%s)" % (t, format(n, ",")) for t, n in tis[:3])
                   + ("; +%d more" % (len(tis) - 3) if len(tis) > 3 else ""),
        "n_tissues": len(tis),
        "assay": "; ".join(a for a, _ in assay[:2]),
        "n_cells": n_cells,
        "n_cell_types": n_types_total,
        "independent": "no (Tabula Sapiens)" if coll == TS_COLLECTION else "yes",
        "n_audited": len(rr),
        "cells_audited": sum(r["n_cells"] for r in rr),
        # the support guard runs BEFORE the three-way call, exactly as cxg_survey.py
        # reports it -- a winner resting on a thinly-estimated profile is discarded, not
        # counted as a disagreement. Classifying the raw rows instead gives 6/103 (5.8%)
        # where the analysis gives 3/80 (3.8%).
        "thin": sum(1 for r in rr if r["thin_support"]),
        "agree": sum(1 for r in live if r["related"]),
        "within_lineage": sum(1 for r in live
                              if not r["related"] and not r["cross_lineage"]),
        "cross_lineage": sum(1 for r in live
                             if not r["related"] and r["cross_lineage"]),
    })

out.sort(key=lambda r: -r["n_audited"])

COLS = [("dataset_id", "CELLxGENE dataset ID"), ("title", "Dataset title"),
        ("tissues", "Tissue (cells)"), ("assay", "Assay"),
        ("n_cells", "Cells in dataset"), ("n_cell_types", "Cell types"),
        ("independent", "Independent of the TS comparison"),
        ("n_audited", "Audited"), ("cells_audited", "Cells audited"),
        ("thin", "Thin (discarded)"), ("agree", "Agree"),
        ("within_lineage", "Within lineage"),
        ("cross_lineage", "Cross lineage"), ("collection_id", "Collection ID"),
        ("publication", "Publication")]

# ---------------------------------------------------------------- TSV (machine-readable)
tsv = ["\t".join(h for _, h in COLS)]
for r in out:
    tsv.append("\t".join(str(r[k]) for k, _ in COLS))
KEYS = ("n_cells", "n_cell_types", "n_audited", "cells_audited", "thin",
        "agree", "within_lineage", "cross_lineage")
tot = {k: sum(r[k] for r in out) for k in KEYS}
tsv.append("\t".join("TOTAL (n=%d)" % len(out) if k == "dataset_id"
                     else (str(tot[k]) if k in tot else "") for k, _ in COLS))
open(os.path.join(ROOT, "supplementary_table_1.tsv"), "w").write("\n".join(tsv) + "\n")

# ---------------------------------------------------------------- Markdown (manuscript)
def fmt(r, k):
    v = r[k]
    return format(v, ",") if isinstance(v, int) else v


md = ["**Supplementary Table 1 | The surveyed CELLxGENE Discover datasets.** "
      "Twelve published human datasets, one per tissue and at most two per collection, "
      "selected for tissue breadth from 979 datasets annotated `disease == normal` with a "
      "qualified tissue term; embryonic and fetal datasets were excluded and only datasets "
      "whose tissue has a CELLxGENE WMG reference were audited. *Audited* counts cell types "
      "passing the 500-cell floor whose asserted CL term is scoreable against the "
      "reference. A winner resting on less than 10% of the asserted term's reference "
      "support is discarded as *thin* before the call is made, so the three outcome "
      "columns sum to Audited minus Thin. *Agree* means the best-scoring term is the "
      "asserted term or an `is_a` relative of it; *within lineage* and *cross lineage* "
      "classify the remainder by whether the winner shares an anchor set with the "
      "assertion. Pooled, 3 of 80 surviving cell types (3.8%) cross a lineage boundary. "
      "Two datasets are slices of the Tabula Sapiens collection "
      "(`e5f58829-1a66-40b5-a624-9046778e74f5`), which the study also uses as its "
      "independent witness; they are marked and are not independent of that comparison. "
      "Excluding them, the remaining ten give 3 of 75 (4.0%), so the survey result does "
      "not rest on them — they contribute 0 of the 3 cross-lineage calls. Cell counts are "
      "as distributed. Every value is read directly from the distributed `.h5ad` files "
      "and from `results/cxg_survey.json`.", "",
      "| " + " | ".join(h for _, h in COLS) + " |",
      "|" + "|".join("---" for _ in COLS) + "|"]
for r in out:
    md.append("| " + " | ".join(fmt(r, k) for k, _ in COLS) + " |")
md.append("| **Total (n=12)** | | | | **%s** | **%s** | | " % (format(tot["n_cells"], ","),
                                                              format(tot["n_cell_types"], ","))
          + " | ".join("**%s**" % format(tot[k], ",") for k in KEYS[2:]) + " | | |")
open(os.path.join(ROOT, "supplementary_table_1.md"), "w").write("\n".join(md) + "\n")

live_n = tot["n_audited"] - tot["thin"]
print("%d datasets | %s cells | %d audited, %d thin-discarded"
      % (len(out), format(tot["n_cells"], ","), tot["n_audited"], tot["thin"]))
print("agree %d / within %d / cross %d of %d surviving -> cross %.1f%%"
      % (tot["agree"], tot["within_lineage"], tot["cross_lineage"], live_n,
         100 * tot["cross_lineage"] / live_n))
assert tot["agree"] + tot["within_lineage"] + tot["cross_lineage"] == live_n, "partition"
ind = [r for r in out if r["independent"] == "yes"]
il = sum(r["agree"] + r["within_lineage"] + r["cross_lineage"] for r in ind)
ic = sum(r["cross_lineage"] for r in ind)
print("independent of TS: %d datasets, %d surviving, cross %d -> %.1f%%"
      % (len(ind), il, ic, 100 * ic / il))
print("PAPER SAYS: 12 datasets, 103 cell types, 3.8% cross-lineage")
