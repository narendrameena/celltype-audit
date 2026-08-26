#!/usr/bin/env python3
"""Which gold organs can CellTypist actually be run on, and why not the rest.

The head-to-head covers seven of the ten curated organs. That is a property of the
comparator, not a choice: a fair run needs a CellTypist model whose label vocabulary the
HRA crosswalk can translate for that organ. This measures the overlap rather than
asserting it, so "no runnable model" is a number.

A level is servable by a model when the model can emit its labels. Blood, bone marrow and
spleen are served by the cross-tissue Immune_All_Low; heart, liver, lung and skin have
dedicated models. Kidney, skeletal muscle and pancreas have crosswalk levels whose
vocabularies are organ-specific -- nephron segments, myofibre types, islet and acinar
types -- and CellTypist ships no adult human model that emits them.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")

import gold_organs                                          # noqa: E402
from ctann import load as load_xw                           # noqa: E402

LEVEL = {"Blood": "blood_L1", "Bone_marrow": "bone marrow_L1", "Spleen": "spleen_L1",
         "Heart": "Healthy_Adult_Heart_pkl", "Liver": "Healthy_Human_Liver_pkl",
         "Lung": "Human_Lung_Atlas_pkl", "Skin": "Adult_Human_Skin_pkl",
         "Kidney": "kidney_L1", "Muscle": "skeletal muscle_L1", "Pancreas": "pancreas_L1"}
MODEL = {"Blood": "Immune_All_Low.pkl", "Bone_marrow": "Immune_All_Low.pkl",
         "Spleen": "Immune_All_Low.pkl", "Heart": "Healthy_Adult_Heart.pkl",
         "Liver": "Healthy_Human_Liver.pkl", "Lung": "Human_Lung_Atlas.pkl",
         "Skin": "Adult_Human_Skin.pkl",
         # the three with no adult human model that speaks their vocabulary
         "Kidney": "Immune_All_Low.pkl", "Muscle": "Immune_All_Low.pkl",
         "Pancreas": "Adult_Human_PancreaticIslet.pkl"}
RUNNABLE = 0.50     # a model serving under half a level's labels is misapplied, not applied


def main():
    from celltypist import models
    _, lv = load_xw("celltypist")
    rows = []
    for organ in gold_organs.curated():
        level, mname = LEVEL.get(organ), MODEL.get(organ)
        if not level or level not in lv:
            rows.append({"organ": organ, "level": level, "model": mname,
                         "labels": 0, "emittable": 0, "frac": 0.0, "runnable": False})
            continue
        try:
            labs = {str(x).lower() for x in models.Model.load(mname).cell_types}
        except Exception as exc:
            rows.append({"organ": organ, "level": level, "model": mname,
                         "error": str(exc)[:60], "runnable": False})
            continue
        keys = set(lv[level])
        hit = len(keys & labs)
        rows.append({"organ": organ, "level": level, "model": mname,
                     "labels": len(keys), "emittable": hit,
                     "frac": round(hit / max(len(keys), 1), 3),
                     "runnable": hit / max(len(keys), 1) >= RUNNABLE})
    json.dump(rows, open(os.path.join(RES, "celltypist_coverage.json"), "w"), indent=1)
    print("can CellTypist be run on each curated organ?\n")
    print("   %-13s %-24s %-30s %s" % ("organ", "crosswalk level", "model", "labels it can emit"))
    for r in sorted(rows, key=lambda x: (not x["runnable"], x["organ"])):
        mark = "yes" if r["runnable"] else "NO "
        print("   %-13s %-24s %-30s %3d/%-3d = %3.0f%%  %s"
              % (r["organ"], r["level"], r["model"], r.get("emittable", 0),
                 r.get("labels", 0), 100 * r.get("frac", 0), mark))
    ok = [r["organ"] for r in rows if r["runnable"]]
    print("\n   runnable on %d of %d curated organs: %s" % (len(ok), len(rows), ", ".join(ok)))
    print("   wrote results/celltypist_coverage.json")


if __name__ == "__main__":
    main()
