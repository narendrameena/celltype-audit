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

It is not one atlas. The same pipeline runs across **every public CELLxGENE Discover
dataset it can score** — 208 met the eligibility criteria, 141 had a scoreable cell type —
and finds **21 of 670 cell types (3.1%, 95% CI 2.1–4.7) carry a label their own markers
contradict at the lineage level**. A twelve-dataset pilot put that at 3.8% [1.3–10.5]; the
estimate barely moved and the interval narrowed from nine points wide to under three.

### Whose label is being tested

Everything audited here is public data, but "public dataset" and "the authors' annotation"
are not the same thing, and the distinction changes what a flag means:

| kind of resource | example | the label under test belongs to |
|---|---|---|
| **aggregating atlas** | hECA v2.0 | the harmonisation pipeline, which rewrote the contributing datasets into one vocabulary |
| **primary atlas** | Tabula Sapiens | the consortium that generated and annotated the cells, natively in CL |
| **published dataset** | CELLxGENE Discover | the depositing authors, mapped into CL under the CELLxGENE schema |

That is not a caveat, it is the finding. The 56,394-cell *neutrophil* cluster above was
**not** an error by the people who generated the cells — hECA keeps each cell's original
name, and 91.2% of them were called a monocyte, macrophage, myeloid or dendritic cell by
their original authors. 22,075 carry a mangled name containing `Monocytenocyte`, across
nine variants: the signature of a string substitution misfiring on a label that already
matched. The harmonised label contradicts the provenance in the same file, through a
text-processing fault rather than a disputed biological judgement.

So what the audit catches depends on what you point it at — **harmonisation faults** in an
aggregating atlas, **annotation judgements** in a primary one, and a **base rate** across
a set of published datasets. It never establishes that the original authors were wrong.

The three read consistently. hECA gives 3.4% (4 of 116) and the 141-dataset survey 3.1%,
which are indistinguishable (Fisher *p* = 0.78) — the aggregating atlas is not worse than
the field it was assembled from. Tabula Sapiens, expert-curated and natively CL-annotated,
gives **0 of 81**; that ordering is what one would expect, but at these counts it is not
established (*p* = 0.15), and scaling the survey eightfold did not change it, because the
limit is Tabula Sapiens' own 81 scoreable cell types rather than the comparator.

## Install

```bash
pip install celltype-audit
```

Only `numpy` and `h5py` are required. The Cell Ontology and the expression reference are
fetched from public endpoints on first use and cached under `~/.cache/celltype-audit`.

### From git

To install `main`, ahead of the latest release:

```bash
pip install "git+https://github.com/narendrameena/celltype-audit.git"
```

A tag, branch or commit pins it:

```bash
pip install "celltype-audit @ git+https://github.com/narendrameena/celltype-audit.git@v0.1.5"
```

### From source

Clone if you want the study as well as the tool. The gold standards, curation sheets,
analysis scripts and figure code are in the repository, **not** in the wheel — `pip
install` gives you `celltype_audit`, nothing else:

```bash
git clone https://github.com/narendrameena/celltype-audit.git
cd celltype-audit
pip install -e ".[dev]"      # editable, with pytest
pytest -q                     # 39 tests; no network and no data needed
```

Use `pip install .` instead of `-e` when you want the package but not a working copy, or
`python -m build` to produce a wheel and sdist under `dist/`.

### Checking what you got

```bash
celltype-audit --version
python -c "import celltype_audit; print(celltype_audit.__version__)"
```

Releases are published from a version tag by the
[publish workflow](.github/workflows/publish.yml), which refuses to run if the tag and
`pyproject.toml` disagree, and runs the tests before anything is uploaded. Uploads use
PyPI Trusted Publishing, so no API token is stored in this repository or in CI.

## Configuration

There is nothing to configure. No account, no API key, no config file — the tool talks to
two public endpoints and caches what it fetches.

**The first run downloads two things** into `~/.cache/celltype-audit`:

| what | from | size |
|---|---|---|
| the Cell Ontology (`cl.json`) | `purl.obolibrary.org` | ~40 MB, once |
| the expression reference, per tissue | `api.cellxgene.cziscience.com` | a few MB per tissue |

Everything already cached is reused, so only the first run for a given tissue needs the
network.

**Moving the cache** — worth doing on HPC, where `$HOME` is often small or read-only from
compute nodes:

```bash
export CELLTYPE_AUDIT_CACHE=/scratch/$USER/celltype-audit    # or: --cache /scratch/...
```

If your compute nodes have no outbound network, run once on a login node with that same
path to warm the cache, then run the real job offline.

**What your `.h5ad` needs.** Two things, and both have an escape hatch:

- **a column of labels to audit.** `obs/cell_type` by default; if yours is called something
  else, `--type-key annotation`. Getting this wrong lists the columns that *are* present
  rather than failing obscurely.
- **a tissue**, so there is something to score against. Taken from the file if it carries
  `tissue_ontology_term_id` (the CELLxGENE layout); otherwise pass `--tissue UBERON:0002048`.
  `--organ Lung` is separate and optional — it only helps resolve organ-qualified labels
  like `Bronchial smooth muscle cell`.

Gene symbols are read from `var/feature_name`, `var/Gene_symbol`, `var/gene_symbols` or
`var/_index`, whichever exists. `celltype-audit` does not cluster; cluster with scanpy or
Seurat first.

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
curation aid, not an annotator.** Measured against 297 hand-curated cell types with no
label involved, the top-scoring term is right **56%** of the time (35% in bone marrow, 87%
in muscle). Its confidence ranks better than chance but not well enough to act on
unattended: leave-one-organ-out AUC 0.72, and precision reaches 80% only over the most
confident 37% of calls and never reaches 90%. There is no threshold that makes it safe to
automate.

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
raises that to **94% of labels and 99% of cells**. Terms scoped to other species (CL
carries 152 mouse-only classes) and obsolete terms are never returned.

**2. Lineage sweep.** Compare the anchor set of the asserted term against the terms the
markers point to. Disjoint on every candidate means the two disagree about *lineage*.
**71% precision** (5 of 7 flags), but structurally blind to any error inside one lineage —
and most errors are: 21 of 28 in the reference study.

**3. Marker queue.** Ask which CL term best explains the cluster, guarded so a merely
coarser or finer name is not called an error. As a binary flag this does **not** reach
usable precision and is **not** shipped as one. Ranked, it puts real errors near the top.

Together: reviewing **33 of 256** cell types surfaced **9 of 19** known errors — 27% of
those reviews landing on a real error, against a 7.4% base rate.

## Limitations, stated

- **It is triage, not an oracle.** The top 33 candidates carry 9 of the 19 known errors;
  reaching 80% of them takes 123 reviews of 256. Yield is what stays fixed as the gold
  grows — recall at a fixed budget falls as more errors are curated into the denominator,
  which says nothing about the method.
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
- **`curation/`** — evidence sheets for the 26 organs that have no hand-curated gold:
  402 cell types with their markers, sizes and candidate CL terms, and the decision
  columns left empty. An organ's sheet is deleted once it has a gold, so the directory is
  always the work that remains. These are not a gold standard and must not be used as one until a
  curator has filled them; the study's own result is that a gold cannot be assembled
  automatically. Regenerate with `python analysis/make_curation_sheets.py`.
- **`analysis/`** — the 47 scripts behind the study, so every number is traceable to the
  code that produced it (`cl_resolve.py`, `audit_ts.py`, `cxg_survey.py`,
  `sensitivity.py`, `downstream_impact.py`, ...).
- **`gold/`** — 10 hand-curated gold standards, 388 cell types across 10 organs — 353
  assigned from marker evidence and **not** from the label, and 35 explicit abstentions
  where no term is defensible.  Files record which entries
  deliberately disagree with the atlas label, and why.
- **`examples/`** — three notebooks: quickstart, unannotated data, and the worked
  case study of the 56,394-cell mislabelled cluster.

Data is not committed; see [DATA.md](DATA.md) for every input and where to get it.

## Citation

If this is useful, please cite the repository until the paper is out.

## Licence

MIT.
