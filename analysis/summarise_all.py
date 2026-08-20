#!/usr/bin/env python3
"""Final report: every hECA v2.0 organ mapped to the Cell Ontology by expression.

Per organ reports the things that actually govern whether mapping can work:
  - cells / uHAF cell types / types kept (>=50 cells) / types mapped to CL
  - the CELLxGENE tissue used, and that tissue's REFERENCE DEPTH (CL terms with >=100 cells)
  - how many uHAF types are well-powered (>=500 cells) and how many are parent categories
    whose children are separately annotated (those cannot get exclusive markers)
  - scored accuracy where a hand-curated gold standard exists

Organs with no CELLxGENE tissue equivalent are listed separately: unmappable for lack of a
reference, which is a coverage limit of the reference, NOT a failure of the method.

Usage: python benchmark/summarise_all.py
"""
import glob
import json
import os
import sys
from collections import defaultdict, deque

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
SHORT = lambda i: i.rsplit("/", 1)[-1].replace("_", ":")

UH = json.load(open(os.path.join(RES, "uhaf_parsed.json")))
UH_CHILD = defaultdict(set)
for a, b in UH["edges"]:
    UH_CHILD[a].add(b)


def has_child(label, present):
    seen, q = set(), deque(UH_CHILD.get(label, ()))
    while q:
        c = q.popleft()
        if c in seen:
            continue
        seen.add(c)
        if c in present:
            return True
        q.extend(UH_CHILD.get(c, ()))
    return False


_REF = None


def _load_ref():
    """term_id_labels covers every tissue, so ONE cached response is enough.
    Scanning the whole cache per organ was O(organs x cache) and unusably slow."""
    global _REF
    if _REF is not None:
        return _REF
    _REF = {}
    for p in sorted(glob.glob(os.path.join(HERE, ".wmg_cache", "*.json")),
                    key=os.path.getsize, reverse=True):
        try:
            d = json.load(open(p))
        except Exception:
            continue
        ct = (d.get("term_id_labels") or {}).get("cell_types")
        if ct:
            _REF = ct
            break
    return _REF


def ref_depth(uberon):
    out = _load_ref().get(uberon, {})
    n100 = sum(1 for k, v in out.items()
               if k != "CL:0000000" and ((v.get("aggregated") or {}).get("total_count", 0)) >= 100)
    cells = sum(((v.get("aggregated") or {}).get("total_count", 0))
                for k, v in out.items() if k != "CL:0000000")
    return n100, cells


# organs the mapper deliberately declines to map (no CELLxGENE tissue equivalent)
sys.path.insert(0, HERE)
from lineage import crowding
try:
    from heca_to_cl import ORGAN2TISSUE
    NO_TISSUE = {k: (v is None) for k, v in ORGAN2TISSUE.items()}
except Exception:
    NO_TISSUE = {}

scores = {}
sp = os.path.join(RES, "organ_scores.json")
if os.path.exists(sp):
    for s in json.load(open(sp)):
        scores[s["organ"]] = s

rows, unmapped = [], []
for mp in sorted(glob.glob(os.path.join(RES, "heca_markers_*.json"))):
    organ = os.path.basename(mp)[len("heca_markers_"):-len(".json")]
    M = json.load(open(mp))
    meta, types = M["meta"], M["types"]
    cp = os.path.join(RES, "heca_to_cl_%s.json" % organ)
    if not os.path.exists(cp):
        # distinguish "no CELLxGENE tissue exists" from "not mapped yet"
        why = ("no CELLxGENE tissue equivalent" if NO_TISSUE.get(organ, False)
               else "not mapped yet (run heca_to_cl.py)")
        unmapped.append((organ, meta["n_cells"], meta["types_kept"], why))
        continue
    C = json.load(open(cp))
    present = set(types)
    n500 = sum(1 for v in types.values() if v["n_cells"] >= 500)
    npar = sum(1 for t in types if has_child(t, present))
    d100, dcells = ref_depth(C["uberon"])
    # a-priori difficulty: how finely the annotation subdivides a lineage (Liver showed this is
    # what governs accuracy). Computed on well-powered types only, from labels alone.
    powered = [t for t, v in types.items() if v["n_cells"] >= 500]
    crowd, biggest, _ = crowding(powered)
    rows.append({"organ": organ, "cells": meta["n_cells"], "types": meta["n_types"],
                 "kept": meta["types_kept"], "mapped": len(C["types"]), "tissue": C["tissue"],
                 "ref_terms": d100, "ref_cells": dcells, "n500": n500, "parents": npar,
                 "crowd": round(crowd, 1), "biggest": biggest,
                 "dropped": len(meta.get("dropped_not_cell_types", []))})

rows.sort(key=lambda r: -r["cells"])
tot_cells = sum(r["cells"] for r in rows)
tot_types = sum(r["kept"] for r in rows)
tot_mapped = sum(r["mapped"] for r in rows)

print("hECA v2.0 -> Cell Ontology by EXPRESSION — all organs\n")
print("%-16s %10s %5s %5s %6s  %-16s %6s %5s %6s %7s %s"
      % ("organ", "cells", "types", "kept", "mapped", "CELLxGENE tissue",
         "refCL", ">=500", "parent", "crowd", "predicted difficulty"))
print("-" * 126)
for r in rows:
    c = r["crowd"]
    pred = "EASY" if c < 6 else ("MODERATE" if c < 12 else "HARD")
    print("%-16s %10d %5d %5d %6d  %-16s %6d %5d %6d %7.1f %s (max lineage %d)"
          % (r["organ"], r["cells"], r["types"], r["kept"], r["mapped"], r["tissue"][:16],
             r["ref_terms"], r["n500"], r["parents"], c, pred, r["biggest"]))
print("-" * 126)
print("crowd = mean number of co-annotated types sharing a lineage (>=500-cell types, from labels")
print("alone). Calibration: pancreas 3.9 -> 89.5% top-1; liver 11.2 -> 54.2%; blood 23.0 -> 45.8%.")
print("%-16s %10d %5s %5d %6d" % ("TOTAL (%d organs)" % len(rows), tot_cells, "",
                                  tot_types, tot_mapped))

if unmapped:
    print("\nOrgans with NO CELLxGENE tissue equivalent (reference coverage limit, not a method failure):")
    for o, c, k, why in unmapped:
        print("   %-16s %10d cells, %d types — %s" % (o, c, k, why))

if scores:
    print("\nScored organs (hand-curated gold standards):")
    print("%-14s %5s %8s %8s %9s  %s" % ("organ", "n", "top-1", "top-5", "ceiling", "stratum"))
    for o, s in scores.items():
        k = "n>=500 & achievable"
        st = s.get(k) or s.get("n>=0 & achievable") or {}
        print("%-14s %5d %7.1f%% %7.1f%% %8.1f%%  %s"
              % (o, st.get("n", 0), st.get("top1", 0), st.get("top5", 0), s.get("ceiling", 0), k))

json.dump({"organs": rows, "unmapped": unmapped}, open(os.path.join(RES, "all_organs_summary.json"), "w"), indent=1)
print("\nwrote results/all_organs_summary.json")
