#!/usr/bin/env python3
"""A LABEL-INDEPENDENT gold standard, built from expert marker sets rather than by hand.

The hand-curated gold is 108 cell types across 4 organs -- too small to carry a claim.
This scales it by reusing the property that makes it valid: the gold term is decided from
MARKER EVIDENCE and never from the label string, so scoring a label-derived mapping
against it is not circular the way the HRA crosswalk comparison is.

Two independent expert resources supply the candidates, because ASCT+B alone covers only
15 organs and leaves most cell types with no comparable term. CellMarker 2.0 carries its
own `cellontology_id`, so bringing it in adds no lexical mapping step of ours.

Construction, deliberately conservative:
  * candidates are expert marker sets (ASCT+B and CellMarker) tagged to that organ
  * a cell type's own data-derived markers (50 per type, deep_markers.py) are intersected
    with each candidate's expert set
  * a gold term is issued only when the best candidate shares >=MIN_OVERLAP markers AND
    beats the runner-up by >=MARGIN, so ambiguous cases abstain rather than guess

This is NOT expert hand-curation and is not reported as such: it is marker-evidence gold,
semi-automated, and its rule is stated so a reader can audit every assignment. Its two
known limits are that ASCT+B covers 15 organs, and that its marker sets are themselves
literature-derived, so it inherits whatever the literature got wrong.

Usage: python benchmark/marker_gold.py
"""
import glob
import json
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")

MIN_OVERLAP = 3
MARGIN = 2.0
MIN_EXPERT = 4          # ignore ASCT+B terms with too few markers to be discriminating

ORGAN2ASCTB = {"Lung": "lung", "Kidney": "kidney", "Thymus": "thymus",
               "Lymph_node": "lymph-node", "Bone_marrow": "bone-marrow",
               "Brain": "allen-brain", "Skin": "skin", "Liver": "liver",
               "Heart": "heart", "Placenta": "placenta-full-term",
               "Pancreas": "pancreas", "Bronchi": "main-bronchus",
               "Uterine_tube": "fallopian-tube", "Prostate": "prostate"}

ORGAN2CELLMARKER = {
    "Blood": "blood", "Brain": "brain", "Bone_marrow": "bone marrow", "Liver": "liver",
    "Lung": "lung", "Kidney": "kidney", "Skin": "skin", "Eye": "eye",
    "Intestine": "intestine", "Heart": "heart", "Pancreas": "pancreas",
    "Spleen": "spleen", "Thymus": "thymus", "Breast": "breast", "Ovary": "ovary",
    "Adipose": "adipose tissue", "Oesophagus": "esophagus",
    "Salivary_gland": "salivary gland", "Uterus": "uterus", "Placenta": "placenta",
    "Testis": "testis", "Muscle": "muscle", "Stomach": "stomach",
    "Prostate": "prostate", "Lymph_node": "lymph node", "Bladder": "bladder",
    "Nose": "nose", "Bronchi": "airway"}


def asctb_by_organ():
    d = json.load(open(os.path.join(RES, "asctb_markers.json")))
    by = defaultdict(dict)
    for c, v in d.items():
        if len(v.get("genes", [])) < MIN_EXPERT:
            continue
        for o in v.get("organs", []):
            by[o][c] = (set(v["genes"]), v.get("label", c))
    return by


def _cl_label(c):
    global _LAB
    if _LAB is None:
        sys.path.insert(0, HERE)
        from cl_lineage import load
        _LAB = load()["label"]
    return _LAB.get(c, c)


_LAB = None


def cellmarker_by_organ():
    p = os.path.join(RES, "cellmarker_human.json")
    if not os.path.exists(p):
        return {}
    d = json.load(open(p))
    by = defaultdict(dict)
    for c, v in d.items():
        for t, genes in (v.get("by_tissue") or {}).items():
            if len(genes) >= MIN_EXPERT:
                by[t][c] = (set(genes), _cl_label(c))
    return by


def candidates_for(organ, ab, cm):
    """Union of the two expert resources for this organ, ASCT+B labels preferred."""
    out = dict(cm.get(ORGAN2CELLMARKER.get(organ, ""), {}))
    out.update(ab.get(ORGAN2ASCTB.get(organ, ""), {}))
    return out


def build():
    by = asctb_by_organ()
    cm = cellmarker_by_organ()
    gold, stats = {}, defaultdict(int)
    for p in sorted(glob.glob(os.path.join(RES, "heca_markers_deep_*.json"))):
        organ = os.path.basename(p)[len("heca_markers_deep_"):-len(".json")]
        cands = candidates_for(organ, by, cm)
        if not cands:
            stats["organ_uncovered"] += 1
            continue
        for t, v in json.load(open(p))["types"].items():
            if v["n_cells"] < 500:
                continue
            stats["assessed"] += 1
            mine = {m["gene"] for m in v["markers"]}
            sc = sorted(((len(mine & gs), c, lab) for c, (gs, lab) in cands.items()),
                        reverse=True)
            if not sc or sc[0][0] < MIN_OVERLAP:
                stats["no_overlap"] += 1
                continue
            best = sc[0]
            second = sc[1][0] if len(sc) > 1 else 0
            if second and best[0] < MARGIN * second:
                stats["ambiguous"] += 1
                continue
            stats["gold"] += 1
            gold["%s|%s" % (organ, t)] = {
                "organ": organ, "label": t, "n_cells": v["n_cells"],
                "gold_curie": best[1], "gold_label": best[2],
                "overlap": best[0], "runner_up": int(second),
                "shared_markers": sorted(mine & cands[best[1]][0]),
                "n_candidates": len(cands)}
    return gold, stats


if __name__ == "__main__":
    gold, stats = build()
    print("MARKER-EVIDENCE GOLD (label-independent; >=%d shared expert markers, >=%.1fx margin)\n"
          % (MIN_OVERLAP, MARGIN))
    for k in ("assessed", "no_overlap", "ambiguous", "gold"):
        print("  %-14s %5d" % (k, stats[k]))
    byorg = defaultdict(int)
    for v in gold.values():
        byorg[v["organ"]] += 1
    print("\n  organs: %s" % ", ".join("%s=%d" % (o, n) for o, n in sorted(byorg.items())))
    print("\n  sample assignments")
    print("  %-11s %-26s %-32s %s" % ("organ", "atlas label", "marker-evidence gold", "shared"))
    print("  " + "-" * 96)
    for v in sorted(gold.values(), key=lambda r: -r["overlap"])[:14]:
        print("  %-11s %-26s %-32s %d: %s"
              % (v["organ"][:11], v["label"][:26], v["gold_label"][:32], v["overlap"],
                 ",".join(v["shared_markers"][:4])))
    json.dump(gold, open(os.path.join(RES, "marker_gold.json"), "w"), indent=1)
    print("\nwrote results/marker_gold.json (%d entries)" % len(gold))
