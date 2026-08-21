#!/usr/bin/env python3
"""Assemble the evidence a curator needs to extend the gold standard.

This does NOT assign gold terms. The study's own result is that a gold cannot be built
automatically -- matching data-derived markers to expert databases agrees with hand
curation only 44-58% -- and assigning terms from the same expression evidence the method
is scored against would make the benchmark circular by construction, which is the failure
the manuscript exists to document.

What it does is remove the part that is mechanical: pulling each cell type's marker set,
its size, and the CL terms in play, into one sheet per organ with the decision columns
left empty. Judgement stays with the curator; fetching does not.

The protocol matches the existing seven organs so the new entries are comparable: the
label is visible, and the rule is that where the markers contradict it, the markers decide.

Usage: python benchmark/make_curation_sheets.py [Organ ...]
Out:   results/curation/<Organ>.tsv
"""
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
OUT = os.path.join(RES, "curation")

from cl_lineage import load                                            # noqa: E402
from cl_resolve import resolve as resolve2                             # noqa: E402

CURATED = {"Pancreas", "Liver", "Blood", "Bone_marrow", "Lung", "Kidney", "Heart"}
MINC = 500
NMARK = 50

HEADER = ["cell_type", "n_cells", "markers_top50", "label_resolves_to", "resolved_name",
          "expression_candidates", "GOLD_CURIE", "GOLD_NAME", "NOTE"]


def main():
    want = set(sys.argv[1:])
    os.makedirs(OUT, exist_ok=True)
    g = load()
    lab = g["label"]
    made = []
    for fp in sorted(glob.glob(os.path.join(RES, "heca_to_cl_*.json"))):
        d = json.load(open(fp))
        organ = d["organ"]
        if organ in CURATED or (want and organ not in want):
            continue
        deep = {}
        dp = os.path.join(RES, "heca_markers_deep_%s.json" % organ)
        if os.path.exists(dp):
            deep = json.load(open(dp))["types"]
        ctx = {c["curie"] for v in d["types"].values() for c in v.get("cl", [])}
        rows = []
        for t, v in sorted(d["types"].items(), key=lambda x: -x[1].get("n_cells", 0)):
            n = v.get("n_cells", 0)
            if n < MINC:
                continue
            cur, _how = resolve2(t, ctx, organ=organ)
            mk = [m["gene"] for m in deep.get(t, {}).get("markers", [])[:NMARK]]
            cands = ["%s %s" % (c["curie"], lab.get(c["curie"], "")) for c in v.get("cl", [])[:5]]
            rows.append([t, str(n), ";".join(mk), cur or "", lab.get(cur, "") if cur else "",
                         " | ".join(cands), "", "", ""])
        if not rows:
            continue
        p = os.path.join(OUT, "%s.tsv" % organ)
        with open(p, "w") as fh:
            fh.write("# %s -- %d cell types >= %d cells. Fill GOLD_CURIE from the markers.\n"
                     "# Where the markers contradict the label, the markers decide; say so in NOTE.\n"
                     "# Leave GOLD_CURIE empty to abstain -- abstention is explicit and is not\n"
                     "# scored as a miss. Do not copy label_resolves_to: that is what is being tested.\n"
                     % (organ, len(rows), MINC))
            fh.write("\t".join(HEADER) + "\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")
        made.append((organ, len(rows)))
    made.sort(key=lambda r: -r[1])
    print("%d sheets, %d cell types awaiting a decision\n" % (len(made), sum(n for _o, n in made)))
    for o, n in made:
        print("   %-16s %3d" % (o, n))
    print("\nwrote %s/" % OUT)


if __name__ == "__main__":
    main()
