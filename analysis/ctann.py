#!/usr/bin/env python3
"""HuBMAP HRA CTann crosswalks: annotation-tool label -> Cell Ontology.

Azimuth, CellTypist and popV each emit their own label vocabulary. HRA publishes an
expert crosswalk from each vocabulary to CL, which is what lets a tool's output be
compared with ours in one shared vocabulary.

Note what this means for evaluation: these crosswalks are LEXICAL, expert-curated maps
from a label string to a CL term. They are the right instrument for translating a tool's
output, and the wrong instrument for deciding whether the underlying annotation is
correct -- see calibrate.py, where a purely lexical resolver reproduces the HRA
crosswalk on 312 of 315 cell types.

Usage: python benchmark/ctann.py
"""
import csv
import os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
XW = os.path.join(os.path.dirname(HERE), "..", "crosswalks")
TOOLS = ("azimuth", "celltypist", "popv")


def load(tool):
    """-> {label_lower: {CL curies}}, {level: {label_lower: {CL curies}}}"""
    p = os.path.abspath(os.path.join(XW, "%s.csv" % tool))
    rows = list(csv.reader(open(p)))
    h = next(i for i, r in enumerate(rows) if r and r[0] == "Organ_Level")
    flat, bylevel = defaultdict(set), defaultdict(lambda: defaultdict(set))
    for r in rows[h + 1:]:
        if len(r) < 6 or not r[5].startswith("CL:"):
            continue
        lvl, lab, cl = r[0].strip(), r[2].strip().lower(), r[5].strip()
        if not lab:
            continue
        flat[lab].add(cl)
        bylevel[lvl][lab].add(cl)
    return dict(flat), {k: dict(v) for k, v in bylevel.items()}


def load_all():
    return {t: load(t) for t in TOOLS}


if __name__ == "__main__":
    A = load_all()
    print("%-11s %8s %8s %10s" % ("tool", "labels", "levels", "ambiguous"))
    print("-" * 42)
    for t in TOOLS:
        flat, bylevel = A[t]
        amb = sum(1 for v in flat.values() if len(v) > 1)
        print("%-11s %8d %8d %10d" % (t, len(flat), len(bylevel), amb))
    ct = A["celltypist"][1]
    print("\ncelltypist levels: %s" % ", ".join(sorted(ct)[:12]))
    print("\nsample celltypist L1 blood entries:")
    for k, v in list(sorted(ct.get("blood_L1", {}).items()))[:6]:
        print("   %-34s %s" % (k[:34], ", ".join(sorted(v))))
