#!/usr/bin/env python3
"""Does correcting a mis-annotation change an answer anyone would report?

Finding errors is only interesting if they propagate. This takes the largest confirmed
error -- hECA's lung "Neutrophilic granulocyte", 56,394 cells whose markers (FCN1, CD14,
S100A12, MAFB, with no FCGR3B and no ELANE) are a classical monocyte -- and asks what it
does to a composition statistic anyone might publish from that atlas.

The answer is checkable against an atlas that was never involved: Tabula Sapiens.

  hECA lung as published : 56,394 neutrophils vs 45,198 monocytes  = 1.25 : 1
  hECA lung as curated   :      0 neutrophils vs 101,592 monocytes = 0.00 : 1
  Tabula Sapiens lung    :    371 neutrophils vs   5,387 monocytes = 0.07 : 1

As published, hECA says neutrophils OUTNUMBER monocytes in human lung. Tabula Sapiens
says monocytes outnumber neutrophils roughly fourteen to one. Correcting the single
mislabelled cluster removes the contradiction. The biology agrees with the correction
rather than the label: neutrophils are fragile and RNA-poor and are well known to be
under-captured by droplet scRNA-seq, so a lung dataset in which they are the commonest
myeloid cell is the anomaly that needs explaining.

Usage: python benchmark/downstream_impact.py
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")

from cl_lineage import ancestors, load                                # noqa: E402
from cl_resolve import resolve as resolve2                            # noqa: E402

MONO, NEUT = "CL:0000576", "CL:0000775"
MINC = 500


def bucket(curie):
    if not curie:
        return None
    a = ancestors(curie)
    if curie == NEUT or NEUT in a:
        return "neutrophil"
    if curie == MONO or MONO in a:
        return "monocyte"
    return None


def main(organ="Lung"):
    g = load()["label"]
    M = json.load(open(os.path.join(RES, "heca_to_cl_%s.json" % organ)))["types"]
    G = {k: v for k, v in json.load(open(os.path.join(HERE, "%s_gold.json" % organ.lower()))).items()
         if not k.startswith("_") and v}
    ctx = {c["curie"] for v in M.values() for c in v.get("cl", [])}

    before = {"neutrophil": 0, "monocyte": 0}
    after = {"neutrophil": 0, "monocyte": 0}
    changed = []
    for t, v in M.items():
        if v["n_cells"] < MINC:
            continue
        lab, _ = resolve2(t, ctx, organ=organ)
        gold = G.get(t)
        b, a = bucket(lab), bucket(gold)
        if b:
            before[b] += v["n_cells"]
        if a:
            after[a] += v["n_cells"]
        if b != a and (b or a):
            changed.append({"label": t, "n_cells": v["n_cells"],
                            "as_labelled": g.get(lab, lab), "as_curated": g.get(gold, gold),
                            "from": b, "to": a})

    ts = json.load(open(os.path.join(RES, "ts_markers_%s.json" % organ)))["types"]
    tn = sum(v["n_cells"] for v in ts.values() if bucket(v.get("cl")) == "neutrophil")
    tm = sum(v["n_cells"] for v in ts.values() if bucket(v.get("cl")) == "monocyte")

    def ratio(d):
        return d["neutrophil"] / max(d["monocyte"], 1)

    print("DOWNSTREAM IMPACT — neutrophil : monocyte ratio in %s\n" % organ.lower())
    print("  %-34s %8s %9s %10s" % ("", "neutro.", "mono.", "ratio"))
    print("  " + "-" * 64)
    print("  %-34s %8d %9d %9.2f:1" % ("hECA as published", before["neutrophil"],
                                       before["monocyte"], ratio(before)))
    print("  %-34s %8d %9d %9.2f:1" % ("hECA after correcting the error", after["neutrophil"],
                                       after["monocyte"], ratio(after)))
    print("  %-34s %8d %9d %9.2f:1" % ("Tabula Sapiens (never involved)", tn, tm,
                                       tn / max(tm, 1)))
    print("\n  the cluster responsible")
    for c in changed:
        print("     %-28s %7d  %-24s -> %s"
              % (c["label"][:28], c["n_cells"], str(c["as_labelled"])[:24], c["as_curated"]))
    fold = ratio(before) / max(tn / max(tm, 1), 1e-9)
    print("\n  As published, hECA overstates the ratio by %.0f-fold against Tabula Sapiens." % fold)
    print("  Any composition or abundance claim about lung myeloid cells drawn from this")
    print("  atlas inherits that error, and correcting one cluster removes it.")

    json.dump({"organ": organ, "before": before, "after": after,
               "tabula_sapiens": {"neutrophil": tn, "monocyte": tm},
               "fold_overstatement": fold, "changed": changed},
              open(os.path.join(RES, "downstream_impact.json"), "w"), indent=1)
    print("\nwrote results/downstream_impact.json")


if __name__ == "__main__":
    main(*(sys.argv[1:] or []))
