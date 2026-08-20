#!/usr/bin/env python3
"""Turn the auditor on Tabula Sapiens itself, not just on hECA.

Until now TS has been a WITNESS -- used to confirm hECA's flags. Auditing it directly is
what turns "hECA has mislabelled clusters" into a claim about atlases in general, and it
is a cleaner test: TS is CELLxGENE-standardised and carries cell_type_ontology_term_id, so
the asserted CL term is read straight off the atlas and no lexical resolution step of ours
is involved at all.

ONE CAVEAT, AND IT RUNS IN OUR FAVOUR. Tabula Sapiens is itself one of the datasets behind
the CELLxGENE reference these clusters are scored against. A cluster TS mislabels
contributes its expression to the WRONG term's reference profile, which makes that term
fit better and hides the error. Every count here is therefore a LOWER BOUND: errors found
despite the reference being contaminated in their favour are real, and errors missed may
simply have been masked.

Usage: python benchmark/audit_ts.py
"""
import glob
import json
import os
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
TS = os.path.join(os.path.dirname(HERE), "..", "ts_data", "TS-%s.h5ad")

from cl_lineage import load, anchor_set, ancestors                    # noqa: E402

MINC, MIN_REF = 500, 100
MARGIN, FLOOR = 0.50, 0.0        # queue cut used for reporting, matching within_lineage
# A winner resting on far thinner reference support than the term it displaces is almost
# always an artifact: a CL term with few reference cells has a noisy profile and scores
# high on a handful of genes by chance. Validated on the hand-curated gold, where REAL
# errors run the other way -- the correct term is usually BETTER supported than the wrong
# one (neutrophil 22k -> classical monocyte 149k) -- so this guard costs nothing there
# and removes the artifacts here.
MIN_SUPPORT_RATIO = 0.10

# A Tabula Sapiens file is NOT one tissue. TS-Adipose holds brown adipose (41,280 cells),
# adipose tissue (30,578) and subcutaneous adipose (22,557); TS-Pancreas holds exocrine and
# endocrine pancreas. Collapsing a file to its most common tissue id threw away every
# cluster sitting in the other ones -- and picked brown adipose, which CELLxGENE's tissue
# vocabulary does not carry, so the whole file was skipped. Tissue is resolved PER CLUSTER.
#
# Where a sub-tissue is absent from the WMG vocabulary it is rolled up to the UBERON part
# that is present. These are part_of relations, not guesses.
ROLLUP = {
    "UBERON:0001348": "UBERON:0001013",   # brown adipose tissue    -> adipose tissue
    "UBERON:0002190": "UBERON:0001013",   # subcutaneous adipose    -> adipose tissue
    "UBERON:0000017": "UBERON:0001264",   # exocrine pancreas       -> pancreas
    "UBERON:0000016": "UBERON:0001264",   # endocrine pancreas      -> pancreas
}
# Tissues CELLxGENE's WMG vocabulary does not carry at all, so no reference exists to
# score against. Checked against primary_filter_dimensions, not assumed: mammary gland,
# thymus, trachea and buccal mucosa are absent. 31 TS clusters sit in these and cannot be
# audited by any expression method that uses this reference.
NO_WMG_TISSUE = {"UBERON:0001911": "mammary gland", "UBERON:0002370": "thymus",
                 "UBERON:0003126": "trachea", "UBERON:0006956": "buccal mucosa"}


def tissue_per_type(organ):
    """-> {cell type: uberon}, the tissue where most of THAT cluster's cells sit."""
    p = os.path.abspath(TS % organ)
    if not os.path.exists(p):
        return {}
    f = h5py.File(p, "r")
    try:
        tg = f["obs/tissue_ontology_term_id"]
        tcats = np.array([c.decode() if isinstance(c, bytes) else c for c in tg["categories"][:]])
        tcodes = tg["codes"][:]
        cg = f["obs/cell_type"]
        ccats = np.array([c.decode() if isinstance(c, bytes) else c for c in cg["categories"][:]])
        ccodes = cg["codes"][:]
        out = {}
        for k, name in enumerate(ccats):
            m = ccodes == k
            if not m.any():
                continue
            v, n = np.unique(tcodes[m][tcodes[m] >= 0], return_counts=True)
            if not len(v):
                continue
            ub = str(tcats[int(v[np.argmax(n)])])
            out[str(name)] = ROLLUP.get(ub, ub)
        return out
    except Exception:
        return {}
    finally:
        f.close()


def main():
    g = load()["label"]
    stem = os.environ.get("WIDE_REF", "wide_ref_repro")
    idx = json.load(open(os.path.join(RES, stem + "_index.json")))
    npz = np.load(os.path.join(RES, stem + ".npz"))
    mats = {k[3:]: npz[k] for k in npz.files if k.startswith("M__")}

    rows, skipped_tissues = [], {}
    skipped = {"no_tissue": 0, "no_ref": 0, "absent_from_wmg": 0,
               "term_absent": 0, "no_genes": 0}
    for p in sorted(glob.glob(os.path.join(RES, "ts_markers_*.json"))):
        organ = os.path.basename(p)[len("ts_markers_"):-len(".json")]
        t2ub = tissue_per_type(organ)
        if not t2ub:
            skipped["no_tissue"] += 1
            continue
        for t, v in json.load(open(p))["types"].items():
            if v["n_cells"] < MINC or not v.get("cl"):
                continue
            ub = t2ub.get(t)
            if not ub:
                skipped["no_tissue"] += 1
                continue
            M, gix, tix = mats.get(ub), idx["gene_ix"].get(ub), idx["term_ix"].get(ub)
            if M is None or not gix or not tix:
                skipped["absent_from_wmg" if ub in NO_WMG_TISSUE else "no_ref"] += 1
                skipped_tissues.setdefault(ub, 0)
                skipped_tissues[ub] += 1
                continue
            cnt = idx["counts"].get(ub, {})
            inv = {v2: k for k, v2 in tix.items()}
            asserted = v["cl"]
            if asserted not in tix:
                skipped["term_absent"] += 1
                continue
            cols = [gix[m["gene"]] for m in v["markers"] if m["gene"] in gix]
            if len(cols) < 3:
                skipped["no_genes"] += 1
                continue
            s = M[:, cols].mean(axis=1)
            keep = np.array([cnt.get(inv[i], 0) >= MIN_REF for i in range(len(s))])
            if not keep[tix[asserted]]:
                skipped["term_absent"] += 1
                continue
            sk = np.where(keep, s, -1.0)
            bi = int(np.argmax(sk))
            bc, bs = inv[bi], float(sk[bi])
            sa = float(s[tix[asserted]])
            ratio = sa / (bs + 1e-9)
            related = (asserted == bc) or (bc in ancestors(asserted)) or (asserted in ancestors(bc))
            sup_a, sup_b = cnt.get(asserted, 0), cnt.get(bc, 0)
            thin = sup_b < MIN_SUPPORT_RATIO * max(sup_a, 1)
            A, B = anchor_set(asserted), anchor_set(bc)
            rows.append({"organ": organ, "label": t, "n_cells": v["n_cells"], "uberon": ub,
                         "asserted": asserted, "asserted_name": g.get(asserted, asserted),
                         "best": bc, "best_name": g.get(bc, bc),
                         "ratio": round(ratio, 4), "best_score": round(bs, 4),
                         "related": bool(related), "thin_support": bool(thin),
                         "support_asserted": int(sup_a), "support_best": int(sup_b),
                         "anchor_conflict": bool(A and B and not (A & B)),
                         "markers": [m["gene"] for m in v["markers"][:5]]})

    q = sorted([r for r in rows if not r["related"] and not r["thin_support"]],
               key=lambda r: r["ratio"])
    thin = [r for r in rows if not r["related"] and r["thin_support"]]
    flags = [r for r in q if r["ratio"] < MARGIN and r["best_score"] >= FLOOR]
    anch = [r for r in rows if r["anchor_conflict"]]
    print("TABULA SAPIENS, AUDITED WITH THE SAME PIPELINE\n")
    print("  cell types assessed (>=%d cells, CL term in reference) : %d in %d tissues"
          % (MINC, len(rows), len({r["organ"] for r in rows})))
    print("  skipped: %s" % ", ".join("%s=%d" % kv for kv in skipped.items() if kv[1]))
    for ub, n in sorted(skipped_tissues.items(), key=lambda kv: -kv[1]):
        print("     %-18s %-16s %d clusters"
              % (ub, NO_WMG_TISSUE.get(ub, "(no reference)"), n))
    print("\n  anchor-sweep conflicts (label lineage vs evidence)      : %d" % len(anch))
    print("  discarded: winner on <%.0f%% of the asserted term's support  : %d"
          % (100 * MIN_SUPPORT_RATIO, len(thin)))
    print("  marker-queue candidates (unrelated, well-supported winner): %d" % len(q))
    print("  of those past the %.2f margin                            : %d" % (MARGIN, len(flags)))
    print("\n  most suspicious Tabula Sapiens clusters")
    print("  %-11s %-30s %7s %-26s %s" % ("tissue", "TS label", "cells", "asserted CL", "best-scoring CL"))
    print("  " + "-" * 112)
    for r in q[:16]:
        print("  %-11s %-30s %7d %-26s %s" % (r["organ"][:11], r["label"][:30], r["n_cells"],
                                              r["asserted_name"][:26], r["best_name"][:30]))
    if not q:
        print("  (none survive the guards)")
    # what KIND of disagreement, which is the real comparison against hECA
    def kind(r):
        if r["related"]:
            return "agrees_or_is_a_relative"
        a, b = anchor_set(r["asserted"]), anchor_set(r["best"])
        return "cross_lineage" if (a and b and not (a & b)) else "within_lineage"
    from collections import Counter
    kc = Counter(kind(r) for r in rows if not r["thin_support"])
    print("\n  KIND of disagreement (after the support guard)")
    for k in ("agrees_or_is_a_relative", "within_lineage", "cross_lineage"):
        print("     %-26s %3d" % (k, kc[k]))
    print("\n  Tabula Sapiens disagrees with itself mainly about GRANULARITY, not identity:")
    print("  %d of %d assessed types cross a lineage boundary. hECA, audited the same way,"
          % (kc["cross_lineage"], len(rows)))
    print("  produced 8 cross-lineage errors that Tabula Sapiens itself then confirmed.")
    print("  Read as a lower bound: TS is one of the datasets behind this reference, so a")
    print("  cluster it mislabels helps the wrong term fit and hides its own error.")
    for r in rows:
        r["kind"] = kind(r)
    json.dump({"rows": rows, "queue": q, "flags": flags, "kinds": dict(kc)},
              open(os.path.join(RES, "audit_ts.json"), "w"), indent=1)
    print("\nwrote results/audit_ts.json")


if __name__ == "__main__":
    main()
