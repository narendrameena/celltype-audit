# Tutorial: auditing an atlas in five minutes

## 1. Point it at an h5ad

```bash
celltype-audit audit TS-Lung.h5ad --organ Lung -o out.json
```

Two things must be present in the file:

- **cell type labels** — `obs/cell_type` by default
- **a tissue** — `obs/tissue_ontology_term_id` (CELLxGENE layout), or pass `--tissue UBERON:0002048`

If the file is CELLxGENE-standardised it will also carry
`obs/cell_type_ontology_term_id`, and celltype-audit uses that CL term directly rather than
re-resolving the label. Gene symbols are read from whichever of `var/feature_name`,
`var/Gene_symbol` or `var/_index` exists.

## 2. Read the two outputs differently

**Lineage-sweep flags** are high-precision and few. Treat each as a claim to check.

**The queue** is ranked triage, not a set of assertions. Work down it until the hit rate
stops justifying the time. In the reference study the top 5 were 60% real and the top 30
were 30% real, against a 7.4% base rate.

## 3. Understand an entry

```json
{
  "atlas_label": "Neutrophilic granulocyte",
  "n_cells": 56394,
  "assignment": {"curie": "CL:0000775", "label": "neutrophil", "source": "lexical"},
  "evidence": {"markers": ["FCN1", "CD14", "S100A12", "MAFB"]},
  "audit": {"best_term_label": "classical monocyte", "ratio": 0.84,
            "related_to_best": false, "thin_support": false}
}
```

The label asserts *neutrophil*; the markers are a classical monocyte's; `related_to_best`
is false so this is not merely a coarser name; `thin_support` is false so the winner is not
resting on a poorly-estimated profile. `ratio` near 1 means the asserted term also fits
reasonably — neutrophils and monocytes genuinely resemble each other — which is why the
*identity of the winner* is more trustworthy than the margin.

## 4. When it abstains

`abstained: true` means the label could not be resolved to a CL term. That is a real
answer, not a failure: an unresolved label is recoverable, a wrongly resolved one is not.
`assignment.exact: false` means only a superclass was recoverable — usable for a lineage
check, **not** for deciding a subtype.

## 5. Making it reproducible

```bash
celltype-audit reference --tissues UBERON:0002048 --genes markers.txt -o lung_ref
```

Ship `lung_ref.npz` and `lung_ref_index.json` alongside a paper and the audit can be
re-run exactly.


## 6. A freshly generated dataset

If the data has never been annotated, there is nothing to audit yet — so start with
proposals:

```bash
celltype-audit annotate fresh.h5ad --cluster-key leiden --tissue UBERON:0002048
```

```
20 clusters (>=500 cells), grouped on obs/leiden

  !! PROPOSALS, NOT CALLS: top-1 is right ~71% of the time (297 curated cell types, 10 organs; ...)

  3      2869  epithelial cell of proximal tubule segment 1
                 2. kidney proximal convoluted tubule epithelial cell   0.86
                 3. epithelial cell of proximal tubule segment 2        0.84
      markers: ALDH7A1, GK, GGACT, DHRS4, FMO1, SLC16A9
```

Three things to read:

- **`relative`** — each proposal's score against the top one. Values near 1.0 mean the
  shortlist does not discriminate, and the cluster needs a human.
- **`lineage_agreement`** — whether the whole top-N sits in one lineage. `false` means the
  evidence is not even settled on the broad class.
- **`status`** — why a cluster got nothing, when it does. There is no silent blank: it
  will say `no CELLxGENE reference for UBERON:...`, or that too few marker genes were
  found, or that no term in the tissue reaches the reference-cell floor.

`tissue_declared` and `rolled_up` record when a cluster's own sub-tissue (say kidney
cortex) had no reference and was rolled up to the part that does (kidney). That is a
part_of relation, and it is reported rather than applied silently.

**The honest workflow for fresh data** is: cluster → annotate with a trained classifier →
`celltype-audit audit`. Use `annotate` to speed manual curation, not to replace it.
