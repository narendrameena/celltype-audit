# Inputs, and where to get them

No data is committed to this repository — the full working set is ~294 GB. Every input is
public, and each is fetched by a script here.

| input | size | how to obtain |
|---|---|---|
| **Cell Ontology** | ~50 MB | fetched automatically on first use, cached in `~/.cache/celltype-audit`. Source: `http://purl.obolibrary.org/obo/cl.json` |
| **CELLxGENE expression reference** | small | `celltype-audit reference --tissues UBERON:... --genes markers.txt -o ref`, or fetched automatically. Public WMG API, no key required |
| **hECA v2.0** | ~223 GB | <https://doi.org/10.5281/zenodo.7908782> — 40 organ `.h5ad` files |
| **Tabula Sapiens** | ~23 GB | CELLxGENE collection `e5f58829-1a66-40b5-a624-9046778e74f5` |
| **CELLxGENE survey datasets** | ~26 GB | selected via the Discover API; see `cellscribe_tool/benchmark/cxg_survey.py` |
| **HuBMAP ASCT+B tables** | ~1 MB | <https://humanatlas.io/asctb-tables> |
| **HRA CTann crosswalks** | ~250 KB | Azimuth / CellTypist / popV → CL, from <https://humanatlas.io> |
| **CellMarker 2.0** | ~8 MB | `Cell_marker_Human.xlsx` from the CellMarker site |

## Reproducing the reference

The expression reference is the one input everything else rests on, so it is rebuildable
from scratch with a stable cache key:

```bash
python cellscribe_tool/benchmark/fetch_reference.py --from-results --out ref
python cellscribe_tool/benchmark/compare_reference.py
```
