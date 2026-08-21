# Figures

One numbered folder per theme. Inside each: `figure/` (svg + pdf + png @400 dpi),
`script/` (the script that made it), `sourceData/` (the exact numbers plotted, TSV).

```
figures/
├── _shared/                       figstyle.py · validate_palette.py
├── 01_atlas_and_qc/               figure/ script/ sourceData/
├── 02_mapping_performance/        figure/ script/ sourceData/
├── 03_scorer_development/         figure/ script/ sourceData/
├── 04_cell_ontology_structure/    figure/ script/ sourceData/
├── 05_granularity_gap/            figure/ script/ sourceData/
├── 06_annotation_refinement/      figure/ script/ sourceData/
├── 07_misannotation_detection/    figure/ script/ sourceData/
├── 08_external_validation/        figure/ script/ sourceData/
├── 09_gold_standard_limits/       figure/ script/ sourceData/
└── 10_error_epidemiology/         figure/ script/ sourceData/
```

Themes 02 and 07 depend on two modules in `cellscribe_tool/benchmark/`: `cl_resolve.py`
(atlas label → CL term; 80.4 % of labels, 98.7 % of cells) and `cl_lineage.py` (lineage as
the CL `is_a` anchor set, never a keyword match). `calibrate.py` records why the mapper is
an auditor and not an assigner.

| # | theme | what it shows |
|---|---|---|
| **01** | Atlas & QC | hECA v2.0 composition: cells and cell types per organ, annotated vs retained vs mapped, cluster-size ECDF (70 % of types have ≥500 cells), CELLxGENE reference depth, exclusions (0.40 M cells in catch-all categories), mapping coverage (36 organs mapped, 4 with no CELLxGENE tissue) |
| **02** | Mapping performance | **agreement** with HuBMAP HRA crosswalks — 315 cell types, 22 organs, top-1 56.8 %, top-5 85.7 %. This is *not* accuracy against truth: the crosswalk is itself a lexical label→CL map, so it trusts the label, and `cl_resolve` reproduces it on **312 of 315**. Panel b is the label-independent measure — hand-curated gold across **7 organs, pooled top-1 55.5 % (n=200)** — against the achievable ceiling. Lung, Kidney and Heart (marked ▲) were curated after the method was fixed and were never tuned on: on them top-1 is **59.8 %** and top-5 **94.6 %**, so the method holds up held-out; panel c is the top-1/top-5 gap (median 28 points) |
| **03** | Scorer development | four ideas failed and one worked; marker count is not the driver; the subspace scorer converts shortlist hits into rank-1 hits (gap narrows in 4 of 5 organs — blood is the exception). `fig3_scorer_top1only` is the same figure with top-5 suppressed, for slides |
| **04** | Cell Ontology structure | what CL actually encodes (textual definitions 95 %, but surface markers 9 % and disjointness 1 %); cell identity is not a partition (122 of 253 lineage pairs blocked by real overlap); a reasoner cannot refute a wrong genus (≤7 % catchable in principle, **0 of 28 actually refuted** after injecting all 131 assertable axioms) |
| **05** | Granularity gap | of CL's 3,537 terms, **2,491 (70%) are human, mature, in vivo and NOT reached** by a 10.8 M-cell atlas — only ~20% of CL is out of scope, so the gap is real; every organ sits above the equal-granularity line; 61% of reached terms are used in one organ and only six are ubiquitous |
| **06** | Annotation refinement | refinement under gates, taken from the term the **label** resolves to — never the expression top-1, which is right about half the time. Of **602** well-powered cell types: 85 lack an exact CL identity (G1), 11 are disputed by the auditor (G2), 301 have no evidenced better-scoring subtype (G3) and 113 are mixed clusters (G4), leaving **92 confident refinements**. Evidence, not CL, is the binding limit — **191 of 491 (39%) terms have no subtype the reference can support** (median 33 exist, 1 supported). Validation: **41 of 42** proposals falling in a hand-curated gold organ land inside the gold term's subtree |
| **07** | Mis-annotation detection | flagging atlas labels whose *own* markers imply a different lineage than the label asserts. The detector is only as good as its lineage oracle: a keyword classifier over label strings is **wrong on 3 and blind on 3** of 26 common CL terms, and every one of its 5 unique flags is refuted by the flagged cells' own markers (a cDC1 called by *CLEC9A*, *BATF3*, *WDFY4*). Reading lineage off CL's `is_a` graph instead — as the anchor set a class descends from, which keeps polyhierarchy — yields **14 flags in 611** resolved cell types. Manual review splits them three ways, because a disagreement has more than one possible culprit: **5** where the markers belong to the evidence lineage and the *label* is wrong (adipose "fibroblast" → *C1QA*, *CD163*, *S100A8*), **4** where they belong to the label's own lineage and *our mapping* is wrong (skin cDC1 → *CLEC9A*, *BATF3*, *WDFY4*), and **5** boundary cases between adjacent lineages. Panel d is the recall measurement, possible only once 7 organs were hand-curated. The anchor sweep compares lineage ANCHOR sets, so an error that stays inside one lineage is invisible to it — it reaches only 2 of 12 known errors on the scorable subset. A second test asks instead **which CL term best explains the cluster**, with an `is_a` guard so a merely coarse label is not called wrong. As a binary flag it caps near 33 % precision and is *not* shipped as one; ranked, it gives **7.6× enrichment at the top 5**. Two-tier — sweep first, then queue — a curator **reviews 33 of 177 cell types (19 %) and finds 9 of 14 errors (64 % recall)**, up from 29 %. Its flagship catch is lung "neutrophilic granulocyte", **56,394 cells**, whose best-scoring term is *classical monocyte* |
| **08** | External validation | the flags tested against evidence that never saw this pipeline. **8 of 10** flags queryable against **Tabula Sapiens** are confirmed and none is refuted — hECA's markers for each flagged cluster light up an OFF-lineage TS cell type (adipose "fibroblast" → *macrophage*, 0.53 vs 0.31; lung "tuft cell" → *mature NK T cell*, 0.31 vs 0.00). **CellTypist**, run on the same cells and mapped to CL through the HRA crosswalk, agrees with the atlas label on **179 of 187** well-powered types; 5 of its 8 disagreements are a smooth-muscle label it reads as pericyte or fibroblast. Two cases — lung tuft cell and bronchial smooth muscle cell — are confirmed by all three lines of evidence independently. Panel d turns the audit on **Tabula Sapiens itself** rather than using it only as a witness: of 115 assessed cell types in 12 tissues, **69 agree with the atlas or differ only by an `is_a` step, 12 differ within a lineage, and **none at all** crosses a lineage boundary** — against **3.4 %** for hECA under the identical pipeline, guard and reference. The method separates a CL-native, expert-curated atlas from a label-derived one. TS is a **lower bound**: it is one of the datasets behind the CELLxGENE reference it is scored against, so an error of its own helps the wrong term fit; a further 32 TS clusters sit in tissues CELLxGENE carries no reference for (thymus, mammary gland, trachea) and cannot be audited at all |
| **09** | Gold-standard limits | a label-independent gold **cannot** be built automatically from expert marker databases. Matching data-derived markers to ASCT+B and CellMarker 2.0 agrees with the 108 hand-curated types only **44–58 %** (single pool) or **57–75 %** (requiring both resources to concur), and tightening the thresholds moves coverage down without moving agreement up. The cause is **size bias**: expert sets run from 2 to 1,329 genes, the set chosen as best has a median of 29 against a typical candidate's 10, and when the call is *wrong* the winner has a median of **94**. Scaling validation needs real curation |
| **10** | Error epidemiology | where annotation errors concentrate, from the 18 curated errors across 190 cell types. The one comparison with the power to reach significance does: clusters of **500–2,000 cells are mislabelled 20.0 % of the time against 6.0 % for larger ones**, 3.3× (Fisher **p = 0.011**), so annotation confidence should scale with cluster size. The confusion matrix shows **9 of 14 placed errors sit on the diagonal** — label and evidence share a lineage — which is exactly the class an anchor-set test cannot see and the reason the sweep needs a second tier. A smooth-muscle label is the commonest single trap (**4/8** TS-confirmed errors, **5/8** CellTypist disagreements, **3/17** hand-curated), counted per source rather than pooled because the denominators are different populations |
| **11** | Downstream impact | what one confirmed error does to a statistic anyone might report. The lung cluster labelled "Neutrophilic granulocyte" (**56,394 cells**) expresses *FCN1*, *CD14*, *S100A12* and *MAFB* and lacks *FCGR3B* and *ELANE* — a classical monocyte. As published, neutrophils **outnumber** monocytes in human lung (**1.25 : 1**); correcting that single cluster takes the ratio to **0 : 1**, and an independent atlas gives **0.07 : 1** — an **18-fold** discrepancy traceable to one label. Atlas-wide, **249,264 cells (2.47 %)** carry a label the evidence contradicts, of which **26,023** are confirmed by an independent atlas; a further 180,750 sit under a CL naming artifact (two sibling ventricular-cardiomyocyte terms with no `is_a` path between them) and are excluded as a labelling question, not an error |
| **12** | Threshold sensitivity | whether the conclusions were chosen by the constants. Four tuned thresholds — support floor, minimum reference cells, minimum cluster size, shortlist depth — are swept one at a time through **the same functions the main figures call**, not a re-implementation (an earlier re-implementation silently dropped the flagship pancreatic "B cell" catch and disagreed with the figure by one error). Across **16 settings**: atlas discrimination holds everywhere (hECA 4.0–7.9 % vs TS 0.0–1.7 %); flag precision is 83 % except where the size floor leaves 2–3 flags and the denominator is the story; the fixed 33-type review budget recovers **44–64 %** of known errors against **18–23 %** at random; and the small-cluster enrichment holds in **direction** everywhere (1.7–2.5×) but **never** in significance (p = 0.14–0.34) — so the enrichment is the claim, not the p-value |
| **13** | Calibration | whether a confidence score can gate automated assignment. It cannot, and the way it fails is the point. A model over the score distribution plus the predicted term's ontology shape reaches **AUC 0.707** against the hand-curated gold and 0.684 against the crosswalks (0.653 on score shape alone) — real ranking signal — yet **90% precision is unreachable on either set**, 80% costs all but the top 22% of calls, and the crosswalk curve peaks at 73.9%. Five of the seven most confident WRONG calls carry the correct term in their own top-5, so the model fails by ranking adjacent types (megakaryocyte/platelet, mast cell/basophil) against each other, not by failing to retrieve. An earlier version of this experiment used only the four organs curated first and reported AUC 0.563; it now uses all seven |

## Reproducing

```bash
cd 0N_<theme>/script && python3 <script>.py          # writes ../figure and ../sourceData
cd 03_scorer_development/script && TOP5=0 python3 fig3_scorer.py   # top-1-only variant
```

`_shared/figstyle.py` holds the print spec (89/120/183 mm, 6.5 pt sans, 0.5 pt rules,
`pdf.fonttype=42` so text stays editable), the palette, and the `save()` helper that writes
all three formats plus source data into the theme layout.

## Palette — validated, do not hand-edit

Okabe-Ito subset `#0072B2` `#D55E00` `#009E73` `#CC79A7`, plus a deliberate achromatic neutral
`#4D4D4D` for baselines and "not applicable".

```bash
python3 _shared/validate_palette.py "#0072B2,#D55E00,#009E73,#CC79A7,#4D4D4D" --pairs all
```

All 10 pairs clear the normal-vision floor (OKLab ΔE ≥ 15) **and** the colour-vision-deficiency
floor (ΔE ≥ 8) under deuteranopia, protanopia and tritanopia. The grey intentionally fails the
chroma check — it is the neutral, not a series hue. Two attempts to "improve" Okabe-Ito by hand
were tested and **failed** (adjacent pairs collapsed to ΔE 2–7 under CVD); use it as published.

## Chart-form choices

Form follows the data's job, so nothing has to be decoded twice:

- comparison against a **baseline** → delta lollipop with a zero line (03a)
- **paired change** → dumbbell / span with an arrow (03c)
- **magnitude across categories** → sorted bars (01a, 02a, 04a)
- **distribution** → ECDF with the analysis threshold marked (01c)

Layout rules that took several passes to get right: legends sit **outside** the axes; value
labels go **above** the marker, never beside it where they collide with tick text; sample sizes
belong on the axis label rather than floating in the panel.
