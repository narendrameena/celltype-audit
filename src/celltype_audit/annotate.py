"""Propose Cell Ontology terms for clusters that have no annotation yet.

READ THIS BEFORE USING IT. Measured against 297 hand-curated cell types in ten organs,
with no label involved at any point, picking the top-scoring CL term is right **71%** of
the time -- 51% in bone marrow, 100% in muscle. A companion experiment (calibration) found
that the confidence ranks better than chance and now gates, but only over a narrow band:
leave-one-organ-out AUC 0.72, with 80% precision over the most confident 51% of calls and
90% over the most confident 17%. A curator can work that top slice unattended; there is
still no threshold that makes the whole output safe to automate.

So this is NOT an annotator, and it is not competitive with a trained classifier such as
CellTypist, Azimuth or popV. What it is good at is the thing those tools do not do: it
returns a **ranked CL shortlist with the marker evidence attached**, and the correct term
is in the top five 81-100% of the time depending on tissue. That is a curation aid --
it turns "name this cluster" into "choose among five, with evidence" -- and it grounds the
result in the ontology from the start rather than mapping to CL afterwards.

The intended workflow for a fresh dataset:

    1. cluster it (scanpy, Seurat -- not this package's job)
    2. `celltype-audit annotate` to get a CL shortlist per cluster, OR annotate with a
       trained classifier, which will be more accurate
    3. `celltype-audit audit` to check whatever annotation you settled on

Step 3 is where this package earns its keep; step 2 is a convenience.
"""
import os

from .markers import marker_table
from .ontology import Ontology
from .reference import ROLLUP, Reference

MIN_CELLS = 500
MIN_REF = 100
TOPN = 5

#: Cluster columns commonly produced by the standard pipelines, tried in order.
CLUSTER_KEYS = ("cell_type", "leiden", "louvain", "cluster", "clusters",
                "seurat_clusters", "annotation", "cell_ontology_class")

ACCURACY_NOTE = (
    "top-1 is right ~71% of the time (297 curated cell types, 10 organs; 51-100% by tissue) "
    "and confidence gates only narrowly (AUC 0.72; 90% precision over the top 17%). Read "
    "the shortlist, not the top hit."
)


def find_cluster_key(path, prefer=None):
    """Which obs column holds the grouping? Raises if there is none."""
    import h5py
    with h5py.File(path, "r") as f:
        keys = list(f["obs"].keys())
    if prefer:
        if prefer not in keys:
            raise ValueError("obs/%s not in the file; obs has: %s"
                             % (prefer, ", ".join(sorted(keys)[:20])))
        return prefer
    for k in CLUSTER_KEYS:
        if k in keys:
            return k
    raise ValueError(
        "no cluster column found in obs (looked for %s).\n"
        "celltype-audit does not cluster -- run scanpy/Seurat first, then pass "
        "--cluster-key with the resulting column." % ", ".join(CLUSTER_KEYS))


def annotate_h5ad(path, tissue=None, cluster_key=None, ontology=None, reference=None,
                  min_cells=MIN_CELLS, topn=TOPN, min_ref=MIN_REF, cache=None,
                  verbose=True):
    """-> {"clusters": [...], "meta": {...}}; each cluster gets a ranked CL shortlist."""
    o = ontology or Ontology.load(cache=cache)
    key = find_cluster_key(path, cluster_key)
    if verbose:
        print("grouping on obs/%s" % key, flush=True)
        print("computing markers ...", flush=True)
    tbl = marker_table(path, type_key=key, min_cells=min_cells)
    tbl = {k: v for k, v in tbl.items() if v["n_cells"] >= min_cells}
    if not tbl:
        raise ValueError("no clusters with >= %d cells" % min_cells)

    # a dataset routinely spans several tissues; resolve PER CLUSTER, rolling a
    # sub-structure up to the part CELLxGENE actually has a reference for
    for v in tbl.values():
        raw = tissue or v.get("tissue")
        v["tissue_raw"] = raw
        v["tissue"] = ROLLUP.get(raw, raw) if raw else None
    tissues = {v["tissue"] for v in tbl.values() if v.get("tissue")}
    if not tissues:
        raise ValueError("no tissue found; pass tissue='UBERON:...' -- there is no "
                         "reference to score against without one")

    if reference is None:
        genes = sorted({m["gene"] for v in tbl.values() for m in v["markers"]})
        if verbose:
            print("fetching reference for %d tissue(s), %d genes ..."
                  % (len(tissues), len(genes)), flush=True)
        reference = Reference.fetch(sorted(tissues), genes, cache=cache, verbose=verbose)

    out = []
    # One subspace per tissue, built from every cluster's markers, so each candidate term
    # is asked the same question. Scoring a cluster only on its own markers asks a
    # different question of each, and cost 9.1 points of top-1 across ten organs.
    SUB = {}
    for _v in tbl.values():
        _ub = _v.get("tissue")
        if _ub:
            SUB.setdefault(_ub, []).append([m["gene"] for m in _v["markers"]])
    SUB = {u: reference.subspace(u, ms) for u, ms in SUB.items()}

    for name, v in sorted(tbl.items(), key=lambda kv: -kv[1]["n_cells"]):
        ub = v.get("tissue")
        markers = v["markers"]
        sc = reference.score(ub, markers, min_ref=min_ref,
                             subspace=SUB.get(ub)) if ub else {}
        ranked = sorted(sc, key=sc.get, reverse=True)[:topn]
        top = sc[ranked[0]] if ranked else 0.0
        # say WHY a cluster got nothing, rather than returning a silent blank
        if ranked:
            status = "ok"
        elif not ub:
            status = "no tissue for this cluster"
        elif ub not in reference.mats:
            status = "no CELLxGENE reference for %s" % ub
        elif len([g for g in markers if g in reference.gene_ix.get(ub, {})]) < 3:
            status = "fewer than 3 marker genes present in the reference"
        else:
            status = "no CL term in %s reaches %d reference cells" % (ub, min_ref)
        out.append({
            "cluster": name, "n_cells": v["n_cells"], "tissue": ub,
            "tissue_declared": v.get("tissue_raw"),
            "rolled_up": bool(v.get("tissue_raw") and v["tissue_raw"] != ub),
            "status": status,
            "markers": markers[:10],
            "proposals": [{"rank": i + 1, "curie": c, "label": o.label(c),
                           "score": round(sc[c], 4),
                           "relative": round(sc[c] / (top + 1e-9), 3),
                           "reference_cells": reference.support(ub, c),
                           "lineage": sorted(o.anchors(c))}
                          for i, c in enumerate(ranked)],
            "lineage_agreement": (len({frozenset(o.anchors(c)) for c in ranked}) == 1
                                  if ranked else None),
        })
    return {"schema": "celltype-audit-proposal/1",
            "meta": {"input": os.path.basename(path), "cluster_key": key,
                     "tissues": sorted(tissues), "min_cells": min_cells, "topn": topn,
                     "cell_ontology_release": o.version,
                     "accuracy": ACCURACY_NOTE,
                     "not_an_annotator": "These are PROPOSALS for curation, not calls. "
                                         "For automated annotation use a trained "
                                         "classifier; then audit it with celltype_audit."},
            "clusters": out}
