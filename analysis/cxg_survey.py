#!/usr/bin/env python3
"""Audit many published CELLxGENE datasets, not just two atlases.

Two atlases support a claim about two atlases. CELLxGENE Discover carries hundreds of
public human datasets, every one natively CL-annotated with UBERON tissue ids -- the exact
shape audit_ts.py already consumes -- so the same question can be asked of the field.

Applies the lessons the Tabula Sapiens audit forced:
  * tissue is resolved PER CLUSTER, never per file: a CELLxGENE dataset routinely spans
    several tissues and collapsing it to the commonest one silently discards the rest
  * sub-tissues absent from the WMG vocabulary roll up to the UBERON part that is present
  * a winner resting on <10%% of the asserted term's reference support is discarded as an
    artifact of a thinly-estimated profile
  * embryonic and fetal datasets are EXCLUDED: the reference is adult-dominated, so they
    would show disagreement that reflects developmental stage rather than annotation

Every count is a lower bound. These datasets are themselves part of the CELLxGENE corpus
behind the reference, so a cluster a dataset mislabels helps the wrong term fit.

Usage: python benchmark/cxg_survey.py [--markers-only]
"""
import argparse
import glob
import json
import os
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
DATA = os.path.join(os.path.dirname(HERE), "..", "cxg_data")
OUT = os.path.join(RES, "cxg_markers")

from ts_markers import markers as compute_markers                     # noqa: E402
from audit_ts import ROLLUP, MIN_SUPPORT_RATIO                        # noqa: E402
from cl_lineage import anchor_set, ancestors, load                    # noqa: E402

MINC, MIN_REF = 500, 100
FETAL = ("embryonic", "fetal", "carnegie", "gestation", "week post-fertilization")


def stage_ok(path):
    """Exclude embryonic/fetal datasets; the reference is adult-dominated."""
    try:
        f = h5py.File(path, "r")
        if "development_stage" not in f["obs"]:
            f.close()
            return True, "unknown"
        g = f["obs/development_stage"]
        cats = [c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]]
        codes = g["codes"][:]
        v, n = np.unique(codes[codes >= 0], return_counts=True)
        dom = str(cats[int(v[np.argmax(n)])]).lower()
        f.close()
        return (not any(k in dom for k in FETAL)), dom
    except Exception:
        return True, "unknown"


def tissue_per_type(path):
    f = h5py.File(path, "r")
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


def build_markers():
    os.makedirs(OUT, exist_ok=True)
    files = sorted(glob.glob(os.path.join(DATA, "*.h5ad")))
    print("datasets on disk: %d" % len(files), flush=True)
    kept = 0
    for p in files:
        did = os.path.basename(p)[:-len(".h5ad")]
        dst = os.path.join(OUT, "%s.json" % did)
        if os.path.exists(dst):
            kept += 1
            continue
        ok, stage = stage_ok(p)
        if not ok:
            print("   skip %s (%s)" % (did[:12], stage[:34]), flush=True)
            json.dump({"skipped": stage}, open(dst, "w"))
            continue
        try:
            m = compute_markers(p)
        except Exception as ex:
            print("   FAIL %s: %s" % (did[:12], str(ex)[:60]), flush=True)
            continue
        t2ub = tissue_per_type(p)
        for t in m:
            m[t]["uberon"] = t2ub.get(t)
        json.dump({"dataset": did, "stage": stage, "types": m}, open(dst, "w"))
        kept += 1
        print("   %s  %d types" % (did[:12], len(m)), flush=True)
    print("marker files: %d" % kept, flush=True)


def audit(stem="survey_ref"):
    """stem: which reference to score against. The survey uses one fetched from the API
    for exactly these tissues and genes (fetch_reference.py), not the cache-assembled
    matrix, which is sparse wherever no earlier query happened to ask for a gene."""
    g = load()["label"]
    if not os.path.exists(os.path.join(RES, stem + ".npz")):
        stem = "wide_ref"
        print("  (survey reference absent; falling back to %s)" % stem)
    idx = json.load(open(os.path.join(RES, stem + "_index.json")))
    npz = np.load(os.path.join(RES, stem + ".npz"))
    mats = {k[3:]: npz[k] for k in npz.files if k.startswith("M__")}
    rows, skipped = [], {"no_ref": 0, "term_absent": 0, "no_genes": 0, "fetal": 0}
    for p in sorted(glob.glob(os.path.join(OUT, "*.json"))):
        d = json.load(open(p))
        if d.get("skipped"):
            skipped["fetal"] += 1
            continue
        for t, v in d.get("types", {}).items():
            if v["n_cells"] < MINC or not v.get("cl") or not v.get("uberon"):
                continue
            ub = v["uberon"]
            M, gix, tix = mats.get(ub), idx["gene_ix"].get(ub), idx["term_ix"].get(ub)
            if M is None or not gix or not tix:
                skipped["no_ref"] += 1
                continue
            a = v["cl"]
            if a not in tix:
                skipped["term_absent"] += 1
                continue
            cols = [gix[m["gene"]] for m in v["markers"] if m["gene"] in gix]
            if len(cols) < 3:
                skipped["no_genes"] += 1
                continue
            cnt = idx["counts"].get(ub, {})
            inv = {vv: k for k, vv in tix.items()}
            s = M[:, cols].mean(axis=1)
            keep = np.array([cnt.get(inv[i], 0) >= MIN_REF for i in range(len(s))])
            if not keep[tix[a]]:
                skipped["term_absent"] += 1
                continue
            sk = np.where(keep, s, -1.0)
            bi = int(np.argmax(sk))
            bc, bs = inv[bi], float(sk[bi])
            sa = float(s[tix[a]])
            sup_a, sup_b = cnt.get(a, 0), cnt.get(bc, 0)
            related = (a == bc) or (bc in ancestors(a)) or (a in ancestors(bc))
            thin = sup_b < MIN_SUPPORT_RATIO * max(sup_a, 1)
            A, B = anchor_set(a), anchor_set(bc)
            rows.append({"dataset": d["dataset"], "label": t, "n_cells": v["n_cells"],
                         "uberon": ub, "asserted": a, "asserted_name": g.get(a, a),
                         "best": bc, "best_name": g.get(bc, bc),
                         "ratio": round(sa / (bs + 1e-9), 4), "related": bool(related),
                         "thin_support": bool(thin),
                         "cross_lineage": bool(A and B and not (A & B)),
                         "markers": [m["gene"] for m in v["markers"][:5]]})
    live = [r for r in rows if not r["thin_support"]]
    cross = [r for r in live if not r["related"] and r["cross_lineage"]]
    within = [r for r in live if not r["related"] and not r["cross_lineage"]]
    print("\nCELLxGENE SURVEY  (reference: %s)\n" % stem)
    print("  datasets audited              : %d" % len({r["dataset"] for r in rows}))
    print("  cell types assessed           : %d" % len(rows))
    print("  cells covered                 : %.1fM" % (sum(r["n_cells"] for r in rows) / 1e6))
    print("  skipped: %s" % ", ".join("%s=%d" % kv for kv in skipped.items() if kv[1]))
    print("\n  after the support guard (%d discarded as thin)" % (len(rows) - len(live)))
    print("     agrees with the label or an is_a relative : %d" % (len(live) - len(cross) - len(within)))
    print("     a different term, SAME lineage            : %d" % len(within))
    print("     a different term, DIFFERENT lineage       : %d (%.1f%%)"
          % (len(cross), 100 * len(cross) / max(len(live), 1)))
    print("\n  largest cross-lineage disagreements")
    for r in sorted(cross, key=lambda x: -x["n_cells"])[:14]:
        print("     %-9s %-30s %7d  %-24s -> %s"
              % (r["uberon"].split(":")[-1], r["label"][:30], r["n_cells"],
                 r["asserted_name"][:24], r["best_name"][:26]))
    json.dump(rows, open(os.path.join(RES, "cxg_survey.json"), "w"), indent=1)
    print("\nwrote results/cxg_survey.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--markers-only", action="store_true")
    ap.add_argument("--audit-only", action="store_true")
    a = ap.parse_args()
    if not a.audit_only:
        build_markers()
    if not a.markers_only:
        audit()
