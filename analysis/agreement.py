#!/usr/bin/env python3
"""Is the gold standard one person's opinion?

The obvious referee objection to a single-annotator gold is that it measures the
annotator, not the ontology. The usual answer is an inter-annotator kappa, which needs a
second curator this study does not have. It does, however, have something better than a
second curator recruited for the purpose: two expert consortia have already done the same
task, independently and for their own reasons.

The HRA / ASCT+B crosswalks assign a Cell Ontology term to a cell-type label, organ by
organ, curated by organ-specific expert panels. That is exactly what the gold does. Where
a label appears in both, the two assignments are an independent double annotation -- made
years apart, by different people, with no knowledge of this study.

Two things have to be got right for the comparison to mean anything.

FIRST, granularity is not disagreement. `effector memory CD8-positive alpha-beta T cell`
and `effector CD8-positive alpha-beta T cell` are one is_a step apart; two curators
choosing different depths on the same path agree about what the cell IS. Those are counted
separately from assignments on unrelated branches, which are real disagreements.

SECOND, and this is the part that inverts the result: HRA maps the NAME. It takes the
label a curator wrote and finds the CL term that label denotes. The gold maps the CELLS --
it reads the markers and assigns the term the expression supports. On a correctly labelled
cluster the two tasks coincide. On a MIS-labelled one they cannot, by construction: HRA
will faithfully return the term for the wrong name. So the audit's own error calls have to
be held out before the residual is read as annotator noise, or the study's central finding
is laundered into evidence against its own gold.

Cohen's kappa is reported because it will be asked for, but it is close to uninformative
here and the number shows why: the category set is the whole of CL, so chance agreement is
about 0.02 and kappa lands within 0.05 of raw agreement. Percent agreement, split by
whether the difference is a granularity step, carries the information kappa compresses away.

Usage:
    python benchmark/agreement.py
"""
import collections
import glob
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
GOLD = os.path.abspath(os.path.join(HERE, "..", "..", "celltype-audit", "gold"))

from cl_lineage import load, ancestors                                  # noqa: E402


def norm(s):
    """Match labels across two vocabularies without inventing equivalences."""
    return re.sub(r"[^a-z0-9]+", " ", s.lower()).strip()


def hra_index():
    """normalised label -> {CL curie}, from the HRA / ASCT+B crosswalks."""
    idx = collections.defaultdict(set)
    tools = set()
    for r in json.load(open(os.path.join(RES, "hra_crosswalks.json"))):
        tools.add(r.get("tool"))
        if r.get("cl"):
            idx[norm(r["label"])].add(r["cl"])
    return idx, sorted(t for t in tools if t)


def gold_index():
    """(organ, normalised label) -> (organ, label, curie), over every curated organ.

    Keyed by the PAIR. The same label is curated separately in each organ and the calls
    can differ -- `Conventional dendritic cell` is an error in both skin and spleen but a
    different one -- so collapsing on the label alone silently discards annotations and
    makes the surviving one depend on directory order.
    """
    idx = {}
    for p in sorted(glob.glob(os.path.join(GOLD, "*_gold.json"))):
        organ = os.path.basename(p)[:-len("_gold.json")]
        for k, v in json.load(open(p)).items():
            if k.startswith("_"):                      # _abstentions, _disagreements, ...
                continue
            curie = v.get("curie") if isinstance(v, dict) else v
            if curie:
                idx[(organ, norm(k))] = (organ, k, curie)
    return idx


def classify(mine, theirs):
    if mine in theirs:
        return "exact"
    if any(t in ancestors(mine) or mine in ancestors(t) for t in theirs):
        return "granularity"
    return "unrelated"


def kappa(counts, n):
    """Cohen's kappa with the category set taken to be the terms actually used.

    Every label draws a near-unique term, so the marginals are flat and chance agreement
    is ~1/n. Reported to show that it is negligible, not because it adds information.
    """
    po = counts["exact"] / n
    pe = 1.0 / n
    return po, pe, (po - pe) / (1 - pe)


def n_for_ci(po, halfwidth=0.10, pe=0.02):
    """Cell types a real double-curation study would need for a given CI half-width."""
    for n in range(20, 2000):
        if 1.96 * math.sqrt(po * (1 - po) / (n * (1 - pe) ** 2)) <= halfwidth:
            return n
    return None


def main():
    g = load()
    H, tools = hra_index()
    G = gold_index()
    wl = json.load(open(os.path.join(RES, "within_lineage.json")))["queue"]
    flagged = {(r["organ"].lower(), r["label"].lower()) for r in wl if r["error"]}

    both = sorted(k for k in G if k[1] in H)
    rows = []
    for key in both:
        organ, label, mine = G[key]
        n = key[1]
        rows.append({"organ": organ, "label": label, "mine": mine,
                     "mine_name": g["label"].get(mine, mine),
                     "theirs": sorted(H[n]),
                     "theirs_name": [g["label"].get(t, t) for t in sorted(H[n])],
                     "class": classify(mine, H[n]),
                     "audit_flags_label": (organ.lower(), label.lower()) in flagged})

    def block(rs, title):
        n = len(rs)
        c = collections.Counter(r["class"] for r in rs)
        po, pe, k = kappa(c, n)
        print("  %-46s n=%3d" % (title, n))
        print("     exact %d (%.1f%%)   granularity %d   unrelated %d"
              % (c["exact"], 100 * po, c["granularity"], c["unrelated"]))
        print("     identity agreement %.1f%%    kappa %.3f (chance %.3f)\n"
              % (100 * (n - c["unrelated"]) / n, k, pe))
        return {"n": n, "exact": c["exact"], "granularity": c["granularity"],
                "unrelated": c["unrelated"], "pct_exact": round(100 * po, 1),
                "pct_identity": round(100 * (n - c["unrelated"]) / n, 1),
                "kappa": round(k, 3), "chance": round(pe, 3)}

    print("Independent double annotation: gold vs HRA/ASCT+B (%s)\n" % ", ".join(tools))
    print("  %d gold (organ, label) annotations, %d HRA labels, %d annotated by both\n"
          % (len(G), len(H), len(both)))
    all_ = block(rows, "all doubly-annotated cell types")
    kept = [r for r in rows if not r["audit_flags_label"]]
    held = block(kept, "excluding types the audit calls mis-labelled")

    print("  the %d unrelated disagreements:" % all_["unrelated"])
    for r in rows:
        if r["class"] == "unrelated":
            print("     %-8s %-30s gold=%-30s HRA=%s%s"
                  % (r["organ"], r["label"][:30], r["mine_name"][:30],
                     ", ".join(r["theirs_name"])[:34],
                     "   <- audit flags this label" if r["audit_flags_label"] else ""))
    n_flagged = sum(1 for r in rows if r["class"] == "unrelated" and r["audit_flags_label"])
    print("\n  %d of %d fall on labels the audit independently calls errors, where HRA"
          % (n_flagged, all_["unrelated"]))
    print("  mapped the name and the gold mapped the cells.")
    need = n_for_ci(all_["pct_exact"] / 100.0)
    print("\n  a purpose-built double curation would need n=%d types for a +/-0.10 CI" % need)

    out = {"tools": tools, "n_gold": len(G), "n_hra": len(H), "n_both": len(both),
           "all": all_, "excluding_flagged": held,
           "unrelated_on_flagged_labels": n_flagged,
           "n_for_ci_0.10": need, "pairs": rows}
    p = os.path.join(RES, "agreement.json")
    json.dump(out, open(p, "w"), indent=1)
    print("\nwrote %s" % p)


if __name__ == "__main__":
    main()
