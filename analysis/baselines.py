#!/usr/bin/env python3
"""Run an established annotator (CellTypist) and audit ITS output, not just ours.

The claim under test is that the contradiction sweep catches real annotation errors. That
claim is only interesting if it catches errors made by tools people actually use, so this
runs CellTypist on the same hECA data, translates its predictions to CL through the HRA
CTann crosswalk, and puts them through the same auditor.

hECA X is already exactly log1p(CP10K) -- verified, expm1 row sums are 10000.0 -- which is
what CellTypist expects, so the matrix is passed through unmodified. Only the model's own
feature genes are read, which keeps a 25 GB organ to a few hundred MB and lets the CSC
organs be read by column instead of by row.

CellTypist runs with over_clustering set to the ATLAS label, so it votes within each
existing cluster and returns one prediction per atlas cell type -- directly comparable to
the label, with no re-clustering of our own in between.

Usage: python benchmark/baselines.py [Organ ...]
"""
import json
import os
import sys
from collections import Counter, defaultdict

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
H5 = os.path.join(os.path.dirname(HERE), "..", "heca_data", "RNA-%s.h5ad")
MAXPER = int(os.environ.get("MAXPER", 300))     # cells sampled per atlas cell type
MINC = 500

# organ -> (CellTypist model, crosswalk Organ_Level key in crosswalks/celltypist.csv)
MODELS = {
    "Liver":       ("Healthy_Human_Liver.pkl", "Healthy_Human_Liver_pkl"),
    "Lung":        ("Human_Lung_Atlas.pkl", "Human_Lung_Atlas_pkl"),
    "Skin":        ("Adult_Human_Skin.pkl", "Adult_Human_Skin_pkl"),
    "Heart":       ("Healthy_Adult_Heart.pkl", "Healthy_Adult_Heart_pkl"),
    # Pancreas is deliberately absent: Adult_Human_PancreaticIslet has only 12 labels, all
    # islet types, so a whole-pancreas run forces fibroblasts and T cells into "delta".
    # That is a misapplied model, not a fair test of CellTypist, so it is not scored.
    "Blood":       ("Immune_All_Low.pkl", "blood_L1"),
    "Bone_marrow": ("Immune_All_Low.pkl", "bone marrow_L1"),
}


def read_subset(path, feats, maxper=MAXPER):
    """-> (matrix cells x feats, atlas label per cell). Reads only the model's genes."""
    f = h5py.File(path, "r")
    gs = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var/Gene_symbol"][:]])
    gi = {}
    for j, s in enumerate(gs):
        gi.setdefault(s.upper(), j)
    cols = np.array([gi.get(s.upper(), -1) for s in feats])
    have = cols >= 0
    ctg = f["obs/cell_type"]
    cats = np.array([c.decode() if isinstance(c, bytes) else c for c in ctg["categories"][:]])
    codes = ctg["codes"][:].astype(np.int64)

    rng = np.random.default_rng(0)
    keep = []
    for k in range(len(cats)):
        idx = np.where(codes == k)[0]
        if len(idx) < MINC:                     # only well-powered clusters
            continue
        keep.append(rng.choice(idx, size=min(maxper, len(idx)), replace=False))
    if not keep:
        f.close()
        return None, None, None
    rows = np.sort(np.concatenate(keep))
    lab = cats[codes[rows]]

    X = np.zeros((len(rows), len(feats)), dtype=np.float32)
    enc = dict(f["X"].attrs).get("encoding-type")
    D, I = f["X/data"], f["X/indices"]
    indptr = f["X/indptr"][:]
    n_genes = len(gs)
    if enc == "csr_matrix":
        col2k = np.full(n_genes, -1, dtype=np.int64)      # gene column -> feature slot
        for k, c in enumerate(cols):
            if c >= 0:
                col2k[c] = k
        for r, ri in enumerate(rows):
            lo, hi = int(indptr[ri]), int(indptr[ri + 1])
            if hi <= lo:
                continue
            k = col2k[I[lo:hi]]
            m = k >= 0
            if m.any():
                X[r, k[m]] = D[lo:hi][m]
    else:                                        # csc: read only the wanted columns
        row2k = np.full(int(codes.shape[0]), -1, dtype=np.int64)
        row2k[rows] = np.arange(len(rows))
        for k, c in enumerate(cols):
            if c < 0:
                continue
            lo, hi = int(indptr[c]), int(indptr[c + 1])
            if hi <= lo:
                continue
            rr = I[lo:hi]
            kk = row2k[rr]
            m = kk >= 0
            if m.any():
                X[kk[m], k] = D[lo:hi][m]
    f.close()
    return X, lab, have.sum()


def main():
    import anndata as ad
    import celltypist
    from celltypist import models
    from ctann import load as load_xw

    organs = sys.argv[1:] or list(MODELS)
    _, ct_lv = load_xw("celltypist")
    out = {}
    for organ in organs:
        mname, lvl = MODELS[organ]
        p = os.path.abspath(H5 % organ)
        if not os.path.exists(p):
            print("missing %s" % p, flush=True)
            continue
        print("=== %s  model=%s ===" % (organ, mname), flush=True)
        model = models.Model.load(model=mname)
        feats = list(model.classifier.features)
        X, lab, nfound = read_subset(p, feats)
        if X is None:
            print("   no well-powered clusters", flush=True)
            continue
        print("   %d cells x %d model genes (%d found in atlas)" % (X.shape[0], X.shape[1], nfound), flush=True)
        a = ad.AnnData(X=X, obs={"atlas_label": lab.astype(str)},
                       var={"gene": feats})
        a.var_names = feats
        pred = celltypist.annotate(a, model=mname, majority_voting=True,
                                   over_clustering="atlas_label")
        pl = pred.predicted_labels
        col = "majority_voting" if "majority_voting" in pl.columns else "predicted_labels"
        res = {}
        for t in sorted(set(lab)):
            m = lab == t
            call = Counter(pl[col].values[m]).most_common(1)[0]
            frac = call[1] / m.sum()
            cl = sorted(ct_lv.get(lvl, {}).get(call[0].lower(), set()))
            res[t] = {"celltypist": call[0], "vote_fraction": round(float(frac), 3),
                      "cl": cl, "n_sampled": int(m.sum())}
        out[organ] = {"model": mname, "crosswalk_level": lvl, "types": res}
        mapped = sum(1 for v in res.values() if v["cl"])
        print("   %d atlas types -> CellTypist; %d mapped to CL via crosswalk"
              % (len(res), mapped), flush=True)
        json.dump(out, open(os.path.join(RES, "baseline_celltypist.json"), "w"), indent=1)
    print("\nwrote results/baseline_celltypist.json")


if __name__ == "__main__":
    main()
