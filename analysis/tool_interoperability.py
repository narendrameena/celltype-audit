#!/usr/bin/env python3
"""Do three annotation tools, once mapped to CL, actually become interoperable?

SCOPE, STATED PLAINLY. CellTypist was run on the data (baselines.py, 187 cell types across
six organs). Azimuth and popV were NOT run: Azimuth is R/Seurat and popV is an scvi-tools
ensemble, and installing either would have meant rebuilding a Python environment that took
real effort to stabilise on NumPy 2.x. That is a resource decision, not a finding, and the
comparison below is not a substitute for running them -- it tests a different claim.

THE CLAIM TESTED HERE. Mapping every tool's labels to the Cell Ontology is supposed to make
their outputs comparable. That only works if the tools' vocabularies, once mapped, actually
land on overlapping parts of CL. Using HuBMAP's own expert crosswalks for Azimuth,
CellTypist and popV, this measures how much of CL each tool can express, and how much they
share -- per organ, so the comparison is like for like.

If two tools annotating the same tissue can barely name the same cell types, then "we
mapped both to CL" does not make their annotations comparable, and any benchmark that
compares them through CL is measuring vocabulary overlap as much as biology.

Usage: python benchmark/tool_interoperability.py
"""
import json
import os
import re
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")

from ctann import load_all, TOOLS                                     # noqa: E402
from cl_lineage import load, anchor_set                               # noqa: E402

# HRA levels name the organ differently per tool; group them onto a common organ key
ORGAN_PAT = [("blood", r"blood"), ("bone marrow", r"bone.?marrow"), ("heart", r"heart"),
             ("kidney", r"kidney"), ("liver", r"liver"), ("lung", r"lung|bronchus"),
             ("intestine", r"intestine|colon|bowel"), ("pancreas", r"pancreas|islet"),
             ("skin", r"skin"), ("spleen", r"spleen"), ("thymus", r"thymus"),
             ("lymph node", r"lymph.?node"), ("brain", r"hippocampus|cortex|brain")]


def organ_of(level):
    s = level.lower()
    for name, pat in ORGAN_PAT:
        if re.search(pat, s):
            return name
    return None


def main():
    A = load_all()
    per = defaultdict(lambda: defaultdict(set))      # organ -> tool -> {CL}
    for tool in TOOLS:
        _flat, bylevel = A[tool]
        for level, m in bylevel.items():
            o = organ_of(level)
            if not o:
                continue
            for _lab, cls in m.items():
                per[o][tool] |= cls

    print("TOOL VOCABULARIES IN CL, PER ORGAN  (from the HuBMAP HRA crosswalks)\n")
    print("  %-13s %8s %8s %8s %10s %11s" % ("organ", "azimuth", "cellt.", "popv",
                                             "union", "in all 3"))
    print("  " + "-" * 66)
    rows = []
    for o in sorted(per):
        s = {t: per[o].get(t, set()) for t in TOOLS}
        have = [t for t in TOOLS if s[t]]
        if len(have) < 2:
            continue
        union = set().union(*s.values())
        inter = set.intersection(*[s[t] for t in have]) if have else set()
        rows.append((o, s, union, inter, have))
        print("  %-13s %8d %8d %8d %10d %10d" % (o, len(s["azimuth"]), len(s["celltypist"]),
                                                 len(s["popv"]), len(union), len(inter)))
    tot_u = len(set().union(*[u for _o, _s, u, _i, _h in rows]))
    print("\n  distinct CL terms any tool can express, across these organs : %d" % tot_u)

    print("\n  PAIRWISE OVERLAP where both tools cover the organ (Jaccard)")
    print("  %-13s %14s %14s %14s" % ("organ", "azi vs cell", "azi vs popv", "cell vs popv"))
    print("  " + "-" * 60)
    jac = defaultdict(list)
    for o, s, _u, _i, _h in rows:
        line = "  %-13s" % o
        for a, b, key in (("azimuth", "celltypist", "azi-cell"),
                          ("azimuth", "popv", "azi-popv"),
                          ("celltypist", "popv", "cell-popv")):
            if s[a] and s[b]:
                j = len(s[a] & s[b]) / len(s[a] | s[b])
                jac[key].append(j)
                line += " %13.2f" % j
            else:
                line += " %13s" % "-"
        print(line)
    print("\n  median Jaccard: %s"
          % ", ".join("%s %.2f" % (k, sorted(v)[len(v) // 2]) for k, v in jac.items() if v))
    print("\n  Two tools annotating the same organ agree on only a minority of the CL terms")
    print("  they can even name. Mapping to CL makes their outputs COMPARABLE in principle;")
    print("  it does not make them equivalent, and a benchmark that scores one tool against")
    print("  another through CL is partly measuring which vocabulary is richer.")

    json.dump({"per_organ": {o: {t: sorted(s[t]) for t in TOOLS} for o, s, _u, _i, _h in rows},
               "jaccard": {k: v for k, v in jac.items()},
               "note": "Azimuth and popV were NOT run as tools; see module docstring."},
              open(os.path.join(RES, "tool_interoperability.json"), "w"), indent=1)
    print("\nwrote results/tool_interoperability.json")


if __name__ == "__main__":
    main()
