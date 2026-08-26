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

from ts_markers import SingleCellType, markers as compute_markers                     # noqa: E402
from audit_ts import ROLLUP, MIN_SUPPORT_RATIO                        # noqa: E402
from cl_lineage import anchor_set, ancestors, load                    # noqa: E402

MINC, MIN_REF = 500, 100
SUBK = 20                      # markers per cluster into the space, as within_lineage
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
        except SingleCellType as ex:
            # recorded like the fetal skips rather than dropped silently: a dataset that
            # cannot yield a discriminative marker is out of scope, not a failure
            print("   skip %s (%s)" % (did[:12], str(ex)[:44]), flush=True)
            json.dump({"skipped": "single cell type"}, open(dst, "w"))
            continue
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


def audit(stem="survey_ref_wide", out="cxg_survey.json", use_subspace=True, scope="dataset"):
    """stem: which reference to score against. The survey uses one fetched from the API
    for exactly these tissues and genes (fetch_survey_ref.py), not the cache-assembled
    matrix, which is sparse wherever no earlier query happened to ask for a gene.

    survey_ref_wide replaced survey_ref because the latter was fetched for each cluster's
    own five markers, which is all the mean scorer asks for. It carried 14 tissues; the
    subspace scorer needs the union of a tissue's markers and found a median 36% of them
    present, so it fell back to the mean and reproduced the mean's number exactly. The
    replacement covers 32 tissues at a median 96% of the subspace, and agrees with the old
    one at r = 1.0000 over 840,599 overlapping values -- the same data, more of it."""
    g = load()["label"]
    if not os.path.exists(os.path.join(RES, stem + ".npz")):
        stem = "wide_ref"
        print("  (survey reference absent; falling back to %s)" % stem)
    idx = json.load(open(os.path.join(RES, stem + "_index.json")))
    npz = np.load(os.path.join(RES, stem + ".npz"))
    mats = {k[3:]: npz[k] for k in npz.files if k.startswith("M__")}
    rows, skipped = [], {"no_ref": 0, "term_absent": 0, "no_genes": 0, "fetal": 0}
    # One subspace per tissue: every cluster's markers in that tissue, pooled. Built once
    # over the whole corpus rather than per file, because a tissue's clusters are spread
    # across datasets and the space has to be the same for every candidate scored in it.
    # scope="dataset" keys the space by (dataset, tissue), which is what the rest of the
    # pipeline does: within_lineage builds one space per marker file, i.e. per experiment.
    # Pooling a whole tissue across every dataset looked more principled and is not: it
    # inflates lymph node to 508 genes while a cluster still contributes ~16 non-zero
    # rates, and cosine with a query that is zero on 97% of its dimensions ranks by how
    # little a candidate expresses OUTSIDE the query, which is noise. That put a
    # 97,130-cell CD4 T cluster's best match at "endothelial cell".
    SUBSPACE = {}
    if use_subspace:
        for _p in sorted(glob.glob(os.path.join(OUT, "*.json"))):
            try:
                _d = json.load(open(_p))
            except Exception:
                continue
            for _t, _v in (_d.get("types") or {}).items():
                _ub = _v.get("uberon")
                if _ub:
                    _k = _ub if scope == "tissue" else (_d.get("dataset", "?"), _ub)
                    SUBSPACE.setdefault(_k, set()).update(
                        m["gene"] for m in _v.get("markers", [])[:SUBK])
        SUBSPACE = {k: sorted(v) for k, v in SUBSPACE.items()}
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
            RATES = {m["gene"]: float(m.get("pc_in", 0.0)) for m in v["markers"]}
            cols = [gix[m["gene"]] for m in v["markers"] if m["gene"] in gix]
            if len(cols) < 3:
                skipped["no_genes"] += 1
                continue
            cnt = idx["counts"].get(ub, {})
            inv = {vv: k for k, vv in tix.items()}
            # Score in the tissue's shared discriminative subspace, the same space the
            # rest of the pipeline uses. Leaving this on the per-cluster mean would have
            # compared hECA under one scorer against these datasets under another, and
            # the section reporting them is called "the same oracle".
            skey = ub if scope == "tissue" else (d["dataset"], ub)
            space = SUBSPACE.get(skey) or []
            scols = [gix[g] for g in space if g in gix]
            qv = np.array([RATES.get(g, 0.0) for g in space if g in gix],
                          dtype=np.float32)
            if scols and float(np.linalg.norm(qv)) > 0:
                Rs = M[:, scols]
                s = Rs.dot(qv) / ((np.linalg.norm(Rs, axis=1) + 1e-9) * float(np.linalg.norm(qv)))
            else:
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
    json.dump(rows, open(os.path.join(RES, out), "w"), indent=1)
    print("\nwrote results/%s" % out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--markers-only", action="store_true")
    ap.add_argument("--audit-only", action="store_true")
    ap.add_argument("--ref", default="survey_ref_wide", help="reference stem to score against")
    ap.add_argument("--out", default="cxg_survey.json")
    ap.add_argument("--scope", default="dataset", choices=("dataset", "tissue"))
    ap.add_argument("--no-subspace", action="store_true",
                    help="score on each cluster's own markers, as before the subspace scorer")
    a = ap.parse_args()
    if not a.audit_only:
        build_markers()
    if not a.markers_only:
        audit(stem=a.ref, out=a.out, use_subspace=not a.no_subspace, scope=a.scope)
