# celltype-audit

[![tests](https://github.com/narendrameena/celltype-audit/actions/workflows/tests.yml/badge.svg)](https://github.com/narendrameena/celltype-audit/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/celltype-audit.svg)](https://pypi.org/project/celltype-audit/)
[![Python](https://img.shields.io/pypi/pyversions/celltype-audit.svg)](https://pypi.org/project/celltype-audit/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**Audit single-cell atlas annotations against the Cell Ontology and expression evidence.**

A cell type's *label* and its *expression* are two independent claims about the same
cells. `celltype-audit` reports where they disagree. It does not re-annotate an atlas — it
tells you where an atlas disagrees with itself, and hands a curator a ranked queue.

In the study this package comes from, it found **8 mislabelled clusters in a published
10.8M-cell atlas**, each independently confirmed against a second atlas, including a
56,394-cell cluster labelled *neutrophil* whose markers (FCN1, CD14, S100A12; no FCGR3B,
no ELANE) are a **classical monocyte**. Correcting that one cluster changes the atlas's
lung neutrophil:monocyte ratio from 1.25:1 to 0:1, against 0.07:1 in an independent atlas
— an 18-fold error in a statistic anyone might publish.

## Install

```bash
pip install celltype-audit
```

Releases are published from a version tag by the [publish workflow](.github/workflows/publish.yml),
which refuses to run if the tag and `pyproject.toml` disagree, and runs the tests before
anything is uploaded.

Only `numpy` and `h5py` are required. The Cell Ontology and the expression reference are
fetched from public endpoints on first use and cached under `~/.cache/celltype-audit`.

## Quickstart

```bash
celltype-audit audit lung.h5ad --organ Lung -o annotations.json --tsv annotations.tsv
```

```
22 cell types (>=500 cells)
  assigned from the label : 22
  abstained               : 0
  lineage-sweep flags     : 0
  marker-queue candidates : 2

review queue, most suspicious first:
   1. alveolar adventitial fibroblast   1113  -> alveolar type 1 fibroblast cell
   2. lung multiciliated epithelial     1209  -> multiciliated columnar cell of tracheobronchial tree
```

```python
from celltype_audit import audit_h5ad
report = audit_h5ad("lung.h5ad", organ="Lung")
for r in report.queue[:10]:
    print(r["atlas_label"], "->", r["audit"]["best_term_label"])
```

### A dataset with no annotation yet

```bash
celltype-audit annotate fresh.h5ad --cluster-key seurat_clusters --tissue UBERON:0002113
```

returns a ranked CL shortlist per cluster with the marker evidence attached. **This is a
curation aid, not an annotator.** Measured against 200 hand-curated cell types with no
label involved, the top-scoring term is right **55%** of the time (35% in bone marrow, 86%
in pancreas) and a right call cannot be told from a wrong one (leave-one-organ-out AUC
0.563). There is no threshold that makes it safe to automate.

What it is good for is that the correct term is in the **top five 73-100%** of the time,
so it turns "name this cluster" into "choose among five, with evidence". For automated
annotation of fresh data use a trained classifier (CellTypist, Azimuth, popV) — then audit
the result with `celltype-audit audit`, which is where this package actually earns its keep.

`celltype-audit` does not cluster. Cluster with scanpy or Seurat first; if no cluster column
is found it says so rather than guessing.

Labels can also be mapped to CL on their own:

```bash
celltype-audit resolve "Neutrophilic granulocyte" "Alpha cell" --organ Pancreas
```

## What it actually does

**1. Resolve the label to a CL term.** Exact lookup finds only about half of real atlas
labels; the misses are conventions, not exotica (`Haematopoietic stem cell` vs CL's `-e-`
spelling, `Cardiomyocyte cell` doubling a suffix, `CD8 T cell` vs
`CD8-positive, alpha-beta T cell`). A normalisation cascade plus organ-qualified lookup
raises that to **80% of labels and 99% of cells**. Terms scoped to other species (CL
carries 152 mouse-only classes) and obsolete terms are never returned.

**2. Lineage sweep.** Compare the anchor set of the asserted term against the terms the
markers point to. Disjoint on every candidate means the two disagree about *lineage*.
**~83% precision**, but structurally blind to any error inside one lineage — and most
errors are: 9 of 14 in the reference study.

**3. Marker queue.** Ask which CL term best explains the cluster, guarded so a merely
coarser or finer name is not called an error. As a binary flag this does **not** reach
usable precision and is **not** shipped as one. Ranked, it puts real errors near the top.

Together: reviewing **33 of 177** cell types surfaced **9 of 14** known errors.

## Limitations, stated

- **Recall is not high.** The two tiers together found 64% of known errors at a 19% review
  budget. This is triage, not an oracle.
- **The reference can hide errors in atlases that built it.** If your dataset is part of
  CELLxGENE, its own mislabelled clusters help the wrong term fit. Counts are a lower bound.
- **Tissue coverage is the reference's, not yours.** CELLxGENE's expression vocabulary has
  no thymus, trachea or mammary gland, so clusters there cannot be scored at all.
- **Thresholds move the numbers.** A sensitivity sweep over five constants is included;
  the atlas-discrimination result holds in every setting, the small-cluster effect holds
  in direction (1.7–2.5×) but its p-value does not.
- **Markers are atlas-relative.** A marker discriminates against whatever else is in *that*
  atlas, so marker lists are not comparable across atlases. Score one atlas's markers
  against another's cells instead.

## Reproducibility

`celltype-audit reference` rebuilds the expression reference from the public API for a stated
list of tissues and genes, with a stable cache key. Nothing depends on a local corpus.

## Tests

```bash
pip install "celltype-audit[dev]" && pytest
```

25 tests, no network required — they run against a synthetic ontology that reproduces the
structural features the code depends on (polyhierarchy, synonym collisions, species-scoped
and obsolete terms).

## What else is in this repository

- **`docs/`** — a review queue for Cell Ontology contributions drafted from the audit's
  by-products, served by GitHub Pages. Every proposal shows which checks were actually run
  against it and which were not; approval happens on the
  [cell-ontology](https://github.com/obophenotype/cell-ontology/issues) and
  [provisional Cell Ontology](https://github.com/obophenotype/provisional_cell_ontology)
  issue trackers, not here. Regenerate with `python analysis/build_proposals.py`, then
  `python analysis/proposal_state.py` to record each proposal and re-check, against the
  current CL release, whether an earlier gap has since closed. A closure is measured, not
  attributed: CL gains terms constantly, so `resolved` says the gap is gone, never that we
  closed it.
- **`curation/`** — evidence sheets for the 29 organs that have no hand-curated gold:
  388 cell types with their markers, sizes and candidate CL terms, and the decision
  columns left empty. These are not a gold standard and must not be used as one until a
  curator has filled them; the study's own result is that a gold cannot be assembled
  automatically. Regenerate with `python analysis/make_curation_sheets.py`.
- **`analysis/`** — the 43 scripts behind the study, so every number is traceable to the
  code that produced it (`cl_resolve.py`, `audit_ts.py`, `cxg_survey.py`,
  `sensitivity.py`, `downstream_impact.py`, ...).
- **`gold/`** — 7 hand-curated gold standards, 200 cell types across 7 organs, each
  assigned from marker evidence and **not** from the label. Files record which entries
  deliberately disagree with the atlas label, and why.
- **`figures/`** — figure code, source data and rendered panels for 13 themes.
- **`examples/`** — three notebooks: quickstart, unannotated data, and the worked
  case study of the 56,394-cell mislabelled cluster.

Data is not committed; see [DATA.md](DATA.md) for every input and where to get it.

## Citation

If this is useful, please cite the repository until the paper is out.

## Licence

MIT.
