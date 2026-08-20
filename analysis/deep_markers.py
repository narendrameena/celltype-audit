#!/usr/bin/env python3
"""Recompute markers at greater depth for the organs ASCT+B covers.

The live marker files hold 5 genes per cell type, but ASCT+B expert marker sets carry a
median of 8, so a 5-gene query overlaps them too rarely to build a gold from. This writes
DEEPER marker sets to a SEPARATE prefix -- heca_markers_deep_<Organ>.json -- because
heca_to_cl.py and figures 2/6/7 read the 5-gene files and must not see different input.

Usage: python benchmark/deep_markers.py [--topk 50]
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
from heca_markers import marker_table                               # noqa: E402

# ASCT+B organ tag -> hECA organ file. palatine-tonsil has no hECA counterpart.
# ASCT+B covers 15 organs; CellMarker covers ~27 of hECA's, so the deep pass is run over
# the union -- the gold can only reach an organ that at least one expert resource covers.
ORGANS = ["Lung", "Kidney", "Thymus", "Lymph_node", "Bone_marrow", "Brain", "Skin",
          "Liver", "Heart", "Placenta", "Pancreas", "Bronchi", "Uterine_tube", "Prostate",
          "Blood", "Eye", "Intestine", "Spleen", "Breast", "Ovary", "Adipose", "Oesophagus",
          "Salivary_gland", "Uterus", "Testis", "Muscle", "Stomach", "Bladder", "Nose"]
H5 = os.path.join(os.path.dirname(HERE), "..", "heca_data", "RNA-%s.h5ad")

if __name__ == "__main__":
    topk = 50
    if "--topk" in sys.argv:
        topk = int(sys.argv[sys.argv.index("--topk") + 1])
    for o in ORGANS:
        out = os.path.join(RES, "heca_markers_deep_%s.json" % o)
        if os.path.exists(out):
            print("skip %s (done)" % o, flush=True)
            continue
        p = os.path.abspath(H5 % o)
        if not os.path.exists(p):
            print("MISSING %s" % p, flush=True)
            continue
        print("computing %s (topk=%d) ..." % (o, topk), flush=True)
        tbl, meta = marker_table(p, topk=topk, min_cells=50, verbose=False)
        meta["organ"] = o
        meta["topk"] = topk
        json.dump({"meta": meta, "types": tbl}, open(out, "w"), indent=1)
        print("   %s: %d types written" % (o, len(tbl)), flush=True)
    print("DEEP MARKERS DONE")
