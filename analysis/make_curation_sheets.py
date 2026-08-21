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

from gold_organs import curated
CURATED = set(curated())
MINC = 500
NMARK = 50

# The pipeline's own candidate ranking is deliberately ABSENT. It is one of the two things
# the gold exists to score, so putting it in front of the curator makes the measurement
# partly circular -- the same failure the study documents for the label. The curator gets
# the markers, the size and the label, and looks up CL directly.
HEADER = ["cell_type", "n_cells", "markers_top50", "label_resolves_to", "resolved_name",
          "GOLD_CURIE", "GOLD_NAME", "NOTE"]


def main():
    want = set(sys.argv[1:])
    os.makedirs(OUT, exist_ok=True)
    g = load()
    lab = g["label"]
    made = []
    for fp in sorted(glob.glob(os.path.join(RES, "heca_to_cl_*.json"))):
        d = json.load(open(fp))
        organ = d["organ"]
        # naming an organ explicitly regenerates it even if a gold already exists, which is
        # what re-curating a contaminated sheet requires
        if (organ in CURATED and not want) or (want and organ not in want):
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
            rows.append([t, str(n), ";".join(mk), cur or "", lab.get(cur, "") if cur else "",
                         "", "", ""])
        if not rows:
            continue
        p = os.path.join(OUT, "%s.tsv" % organ)
        with open(p, "w") as fh:
            fh.write("# %s -- %d cell types >= %d cells. Fill GOLD_CURIE from the markers.\n"
                     "# Where the markers contradict the label, the markers decide; say so in NOTE.\n"
                     "# Leave GOLD_CURIE empty to abstain -- abstention is explicit and is not\n"
                     "# scored as a miss. Do not copy label_resolves_to: that is what is being\n"
                     "# tested. The pipeline's own ranking is not shown here, for the same reason.\n"
                     % (organ, len(rows), MINC))
            fh.write("\t".join(HEADER) + "\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")
        made.append((organ, len(rows)))

    # Drop sheets for organs that have since been curated. Skipping them on write is not
    # enough: the file from before the gold existed stays on disk, and the directory goes
    # on advertising "organs with no hand-curated gold" while three of them have one --
    # inviting a curator to redo finished work. Only ever removes a sheet whose organ now
    # has a gold; a full regeneration is not a licence to delete anything else.
    pruned = []
    if not want:
        for fp in sorted(glob.glob(os.path.join(OUT, "*.tsv"))):
            organ = os.path.basename(fp)[:-4]
            if organ in CURATED:
                os.remove(fp)
                pruned.append(organ)

    made.sort(key=lambda r: -r[1])
    print("%d sheets, %d cell types awaiting a decision\n" % (len(made), sum(n for _o, n in made)))
    for o, n in made:
        print("   %-16s %3d" % (o, n))
    if pruned:
        print("\n   pruned (now curated): %s" % ", ".join(pruned))
    print("\nwrote %s/" % OUT)


if __name__ == "__main__":
    main()
