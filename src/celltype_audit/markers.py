"""Marker genes per cell type, computed from the data rather than looked up.

The statistic is the NS-Forest-style binary score, pc_in * (1 - max pc_out): a gene scores
well when it is detected in this cell type and in no other. It is computed from detection
RATES, so it does not assume a normalisation.

One consequence is worth stating because it defeats the obvious cross-atlas comparison: a
marker chosen this way discriminates a cell type from the OTHER types present in THAT
atlas, so the same cell type gets different top markers in different atlases. Two atlases
here share no top markers even for macrophage. To compare across atlases, score one
atlas's markers against the other's cells rather than comparing the marker lists.
"""
import re

import h5py
import numpy as np

#: Genes that dominate marker lists without identifying anything.
JUNK = re.compile(r"^(MT-|MTRNR|RPL|RPS|RP[0-9]|LINC[0-9]|AC[0-9]{6}|AL[0-9]{6}|"
                  r"AP[0-9]{6}|Z[0-9]{5}|CT[AB]-|MIR[0-9])")
SEX = {"XIST", "TSIX", "RPS4Y1", "DDX3Y", "UTY", "USP9Y", "EIF1AY", "KDM5D", "NLGN4Y"}
NOT_A_CELL_TYPE = {"", "unknown", "unclassified", "doublet", "doublets", "mixed", "nan"}


def _decode(a):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def _cat(f, key):
    """Read a categorical or plain obs/var column as (categories, codes).

    Three layouts, all of them in the wild:

      current    the column is a Group holding `categories` and `codes`
      legacy     the column is an integer Dataset of codes and the categories live in a
                 sibling `__categories/<name>` group. AnnData wrote this before 0.8 and
                 HuBMAP still publishes it. Reading such a column as values rather than
                 codes does not raise -- it silently returns integers where gene symbols
                 were expected, which is worse than a crash.
      plain      the column is an array of strings

    A code of -1 means "no category", which is how a missing HUGO symbol is stored. Left
    as a negative index it would quietly select the LAST category, so it is mapped to an
    empty string and dropped downstream.
    """
    g = f[key]
    if isinstance(g, h5py.Group) and "categories" in g:
        return _decode(g["categories"][:]), g["codes"][:]
    parent, _, name = key.rpartition("/")
    legacy = f.get("%s/__categories" % parent) if parent else None
    if legacy is not None and name in legacy and np.issubdtype(g.dtype, np.integer):
        cats = np.append(_decode(legacy[name][:]), "")     # index -1 -> ""
        return cats, g[:]
    v = _decode(g[:])
    return np.unique(v, return_inverse=True)


def marker_table(path, gene_key=None, type_key="cell_type", topk=20, min_cells=50,
                 block=20000):
    """-> {cell type: {n_cells, markers:[{gene,score,pc_in}], cl, tissue}}

    gene_key defaults to whichever of var/feature_name, var/Gene_symbol, var/_index
    exists, so CELLxGENE-standard and other layouts both work.
    """
    f = h5py.File(path, "r")
    try:
        if gene_key is None:
            for k in ("feature_name", "Gene_symbol", "gene_symbols", "hugo_symbol",
                      "_index"):
                if k in f["var"]:
                    gene_key = k
                    break
        gcats, gcodes = _cat(f, "var/%s" % gene_key)
        genes = gcats[gcodes]
        if type_key not in f["obs"]:
            # A bare KeyError here reads as a broken file. It is almost always a naming
            # difference -- annotations live in `celltype`, `Cell_type`, `annotation`,
            # `free_annotation` about as often as in `cell_type` -- so say what is there.
            cols = [k for k in f["obs"] if not k.startswith("_")]
            raise KeyError(
                "obs/%s not found in %s. Pass --type-key with the column holding the "
                "cell-type labels. Available: %s"
                % (type_key, path, ", ".join(sorted(cols)) or "(none)"))
        tcats, tcodes = _cat(f, "obs/%s" % type_key)
        codes = tcodes.astype(np.int64)
        nt, n_genes = len(tcats), len(genes)
        npt = np.bincount(codes[codes >= 0], minlength=nt).astype(float)

        # optional CELLxGENE extras
        cl_of, ub_of = {}, {}
        for key, dest in (("cell_type_ontology_term_id", cl_of),
                          ("tissue_ontology_term_id", ub_of)):
            if key in f["obs"]:
                c2, k2 = _cat(f, "obs/%s" % key)
                for i in range(nt):
                    m = codes == i
                    if m.any():
                        v, n = np.unique(k2[m][k2[m] >= 0], return_counts=True)
                        if len(v):
                            dest[tcats[i]] = str(c2[int(v[np.argmax(n)])])

        X = f["X"]
        n_cells = X.attrs["shape"][0] if "shape" in X.attrs else len(codes)
        det = np.zeros((nt, n_genes), dtype=np.float32)
        enc = dict(X.attrs).get("encoding-type")
        if enc == "csr_matrix":
            indptr = f["X/indptr"][:]
            D, I = f["X/data"], f["X/indices"]
            for s in range(0, int(n_cells), block):
                e = min(s + block, int(n_cells))
                lo, hi = int(indptr[s]), int(indptr[e])
                if hi <= lo:
                    continue
                gcol, vals = I[lo:hi], D[lo:hi]
                cnt = np.diff(indptr[s:e + 1]).astype(np.int64)
                tt = np.repeat(codes[s:e], cnt)
                m = (vals > 0) & (tt >= 0)
                np.add.at(det, (tt[m], gcol[m]), 1.0)
        elif enc == "csc_matrix":
            indptr = f["X/indptr"][:]
            D, I = f["X/data"], f["X/indices"]
            for j in range(n_genes):
                lo, hi = int(indptr[j]), int(indptr[j + 1])
                if hi <= lo:
                    continue
                rr, vv = I[lo:hi], D[lo:hi]
                m = vv > 0
                np.add.at(det, (codes[rr[m]], j), 1.0)
        else:
            for s in range(0, int(n_cells), 2000):
                e = min(s + 2000, int(n_cells))
                B = X[s:e, :]
                cb = codes[s:e]
                for t in np.unique(cb[cb >= 0]):
                    det[t] += (B[cb == t] > 0).sum(axis=0)
        det = det / np.maximum(npt, 1)[:, None]

        keep = np.array([bool(g) and not JUNK.search(g) and g not in SEX for g in genes])
        out = {}
        for k in range(nt):
            name = str(tcats[k])
            if npt[k] < min_cells or name.lower() in NOT_A_CELL_TYPE:
                continue
            inn = det[k]
            oth = np.delete(det, k, axis=0)
            score = inn * (1.0 - (oth.max(axis=0) if len(oth) else 0.0))
            score = np.where(keep, score, -1.0)
            top = np.argsort(-score)[:topk]
            out[name] = {"n_cells": int(npt[k]),
                         "cl": cl_of.get(name), "tissue": ub_of.get(name),
                         "markers": [{"gene": str(genes[j]), "score": round(float(score[j]), 4),
                                      "pc_in": round(float(inn[j]), 3)}
                                     for j in top if score[j] > 0]}
        return out
    finally:
        f.close()
