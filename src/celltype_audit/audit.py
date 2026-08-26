"""The audit: where a label and its own expression disagree.

Two tests, because one is not enough and neither subsumes the other.

  LINEAGE SWEEP. Compare the anchor set of the term the label asserts against the anchor
  sets of the terms the evidence points to. Disjoint on every candidate means the two
  disagree about lineage. High precision (~83% against hand-curated gold) but structurally
  blind to any error that stays inside one lineage -- which is most of them: 9 of 14
  placed errors in the reference study sat inside a single lineage.

  MARKER QUEUE. Ask which CL term best explains the cluster, guarded so that a merely
  coarser or finer name is not called an error and a winner resting on a thinly-supported
  term is discarded. As a binary flag this does not reach usable precision (~33% at a 12%
  base rate) and it is NOT shipped as one; ranked, it puts real errors near the top.

Together they are a triage system, not an oracle: read the sweep's flags, then work down
the queue. In the reference study that surfaced 9 of 14 known errors from reviewing 33 of
177 cell types.
"""
import json
import os

from .markers import marker_table
from .ontology import Ontology
from .reference import Reference
from .resolve import Resolver

MIN_CELLS = 500
MIN_REF = 100
MIN_SUPPORT_RATIO = 0.10     # a winner on <10% of the asserted term's support is noise
TOPK = 5


class Report:
    """The audit's output: one record per cell type, plus the ranked review queue."""

    def __init__(self, records, meta):
        self.records, self.meta = records, meta

    @property
    def flagged(self):
        return [r for r in self.records if r["audit"]["contradicted"]]

    @property
    def queue(self):
        """Most-suspicious first, by how much better the winner explains the cluster.

        A record with no `ratio` cannot be ranked and is not a candidate: the ratio is the
        asserted term's score over the winner's, so it is None exactly when the asserted
        term has no reference profile to score. That is common rather than exotic -- the
        reference reaches about 12% of human CL classes -- and sorting None against a float
        raised TypeError on the first genuinely external atlas this was pointed at, an
        86,530-cell liver endothelial atlas whose types are finer than the reference
        carries. Excluded here rather than coerced, because "unrankable" is not "rank 0".
        """
        q = [r for r in self.records
             if r["audit"]["best_term"] and not r["audit"]["related_to_best"]
             and not r["audit"]["thin_support"]
             and r["audit"].get("ratio") is not None]
        return sorted(q, key=lambda r: r["audit"]["ratio"])

    @property
    def unscoreable(self):
        """Records the reference cannot score, which a caller should see rather than lose."""
        return [r for r in self.records if r["audit"].get("ratio") is None]

    def summary(self):
        n = len(self.records)
        ab = sum(1 for r in self.records if r["abstained"])
        gen = sum(1 for r in self.records
                  if not r["abstained"] and not r["assignment"]["exact"])
        L = ["%d cell types (>=%d cells)" % (n, self.meta["min_cells"]),
             "  assigned from the label : %d" % (n - ab),
             "     of which superclass only : %d" % gen,
             "  abstained               : %d" % ab,
             "  lineage-sweep flags     : %d" % len(self.flagged),
             "  marker-queue candidates : %d" % len(self.queue),
             "  unscoreable against the reference : %d" % len(self.unscoreable)]
        return "\n".join(L)

    def to_json(self, path):
        json.dump({"schema": "celltype-audit-annotation/1", "meta": self.meta,
                   "annotations": self.records,
                   "queue": [{"rank": i + 1, "organ": r["organ"], "label": r["atlas_label"],
                              "ratio": r["audit"]["ratio"],
                              "best_term": r["audit"]["best_term_label"]}
                             for i, r in enumerate(self.queue)]},
                  open(path, "w"), indent=1)

    def to_tsv(self, path):
        cols = ("organ", "atlas_label", "n_cells", "cl_curie", "cl_label", "source",
                "exact", "abstained", "contradicted", "best_term", "ratio")
        with open(path, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in self.records:
                a, u = r["assignment"], r["audit"]
                fh.write("\t".join(str(x) for x in [
                    r["organ"], r["atlas_label"], r["n_cells"], a["curie"] or "",
                    a["label"] or "", a["source"], a["exact"], r["abstained"],
                    u["contradicted"], u["best_term_label"] or "", u["ratio"]]) + "\n")


def audit_h5ad(path, organ=None, tissue=None, ontology=None, reference=None,
               min_cells=MIN_CELLS, topk=TOPK, min_ref=MIN_REF,
               support_ratio=MIN_SUPPORT_RATIO, cache=None, verbose=True,
               type_key="cell_type"):
    """Audit one h5ad and return a Report.

    tissue: UBERON id. Taken from the file when it carries one (CELLxGENE layout),
    otherwise it must be given -- there is no reference to score against without it.
    type_key: obs column holding the labels under audit. Only `cell_type` was accepted
    before, which silently excluded every file not written in the CELLxGENE layout.
    """
    o = ontology or Ontology.load(cache=cache)
    res = Resolver(o)
    if verbose:
        print("computing markers ...", flush=True)
    tbl = marker_table(path, type_key=type_key, min_cells=min_cells)
    tbl = {k: v for k, v in tbl.items() if v["n_cells"] >= min_cells}
    if not tbl:
        raise ValueError("no cell types with >= %d cells in %s" % (min_cells, path))

    tissues = {v.get("tissue") for v in tbl.values() if v.get("tissue")}
    if tissue:
        tissues = {tissue}
    if not tissues:
        raise ValueError("no tissue_ontology_term_id in the file; pass tissue='UBERON:...'")

    if reference is None:
        genes = sorted({m["gene"] for v in tbl.values() for m in v["markers"]})
        if verbose:
            print("fetching reference for %d tissue(s), %d genes ..."
                  % (len(tissues), len(genes)), flush=True)
        reference = Reference.fetch(sorted(tissues), genes, cache=cache, verbose=verbose)

    ctx = {c for t in tissues for c in reference.term_ix.get(t, {})}
    records = []
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
        ub = tissue or v.get("tissue")
        native = v.get("cl")
        cur, how = (native, "native") if native else res.resolve(name, ctx, organ=organ)
        markers = v["markers"]
        sc = reference.score(ub, markers, min_ref=min_ref,
                             subspace=SUB.get(ub)) if ub else {}
        best = max(sc, key=sc.get) if sc else None
        A = o.anchors(cur) if cur else frozenset()
        cand = sorted(sc, key=sc.get, reverse=True)[:topk]
        anc = [o.anchors(c) for c in cand]
        # Unanimity is over the candidates that ASSERT a lineage. A term with no anchor
        # set says nothing about lineage -- which is why an unanchored asserted term is
        # never flagged -- and the same reasoning has to apply when it appears among the
        # candidates, or one such term silently vetoes an otherwise unanimous
        # disagreement. Measured on the ten curated organs: 6 flags -> 9, precision 100%
        # either way, and the flags it recovers include a kidney cluster labelled
        # "Neuroglial cell" whose gold term is epithelial. Found by running the shipped
        # tool on an independent liver atlas, where four hepatocyte candidates and one
        # unanchored `hepatoblast` left a 22,202-cell cluster carrying APOA1, APOC3,
        # SERPINA1 and HP unflagged.
        anchored = [B for B in anc if B]
        contradicted = bool(A) and bool(anchored) and all(not (A & B) for B in anchored)
        ratio = (sc.get(cur, 0.0) / (sc[best] + 1e-9)) if (best and cur in sc) else None
        sup_a = reference.support(ub, cur) if cur else 0
        sup_b = reference.support(ub, best) if best else 0
        records.append({
            "organ": organ or ub, "atlas_label": name, "n_cells": v["n_cells"], "tissue": ub,
            "assignment": {"curie": cur, "label": o.label(cur) if cur else None,
                           "source": ("atlas-native" if how == "native"
                                      else ("abstained" if not cur else "lexical")),
                           "resolution_status": how,
                           "exact": bool(cur) and (how == "native" or res.is_exact(how))},
            "abstained": cur is None, "abstain_reason": None if cur else how,
            "evidence": {"markers": markers[:10],
                         "candidates": [{"curie": c, "label": o.label(c),
                                         "score": round(sc[c], 4),
                                         "reference_cells": reference.support(ub, c)}
                                        for c in cand]},
            "audit": {"contradicted": contradicted,
                      "asserted_lineage": sorted(A),
                      "evidence_lineage": sorted(set().union(*anc) if anc else set()),
                      "best_term": best, "best_term_label": o.label(best) if best else None,
                      "ratio": round(ratio, 4) if ratio is not None else None,
                      "related_to_best": bool(cur and best and o.related(cur, best)),
                      "thin_support": bool(best and sup_b < support_ratio * max(sup_a, 1))},
        })
    meta = {"input": os.path.basename(path), "organ": organ,
            "tissues": sorted(tissues), "min_cells": min_cells, "topk": topk,
            "min_ref": min_ref, "support_ratio": support_ratio,
            "cell_ontology_release": o.version}
    return Report(records, meta)
