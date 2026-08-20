#!/usr/bin/env python3
"""Parse HuBMAP ASCT+B tables into a CL-term -> expert-curated marker-gene resource.

ASCT+B (Borner et al., Nat Cell Biol 2021; HRA v2.0, Borner et al., Nat Methods 2024) is an
EXPERT-CURATED mapping of Anatomical Structures -> Cell Types -> Biomarkers, with CL and Uberon
CURIEs. It is independent of both CL's own logical definitions and of CELLxGENE expression, so
it is the natural third-party gold standard for marker-based cell-type identity.

Row layout: AS/n[,LABEL,ID] ... CT/n[,LABEL,ID] ... BGene/n[,LABEL,ID] ...
The DEEPEST populated CT/n/ID in a row is the cell type the row's biomarkers characterise.

Usage: python benchmark/parse_asctb.py [<asctb_dir>]
Output: results/asctb_markers.json  {CL curie: {label, genes[], organs[], n_rows}}
"""
import csv
import glob
import json
import os
import re
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
DEFAULT_DIR = os.path.join(os.path.dirname(os.path.dirname(HERE)), "asctb")
os.makedirs(RES, exist_ok=True)

CT_ID = re.compile(r"^CT/(\d+)/ID$")
CT_LAB = re.compile(r"^CT/(\d+)/LABEL$")
# NOTE: ASCT+B tables are INCONSISTENT about which column holds the gene SYMBOL.
# Pancreas: BGene/n = "Pancreatic alpha-amylase" (name), BGene/n/LABEL = "AMY2A" (symbol).
# Heart/Brain: the two are SWAPPED. So scan BOTH and keep whichever looks like a symbol.
B_ANY = re.compile(r"^B[A-Za-z]*/(\d+)(/LABEL)?$")
SYMBOL = re.compile(r"[A-Z][A-Z0-9\-]{1,14}")


def parse(path):
    rows = list(csv.reader(open(path, encoding="utf-8", errors="replace")))
    hi = None
    for i, r in enumerate(rows[:40]):
        if r and r[0].strip() == "AS/1":
            hi = i
            break
    if hi is None:
        return []
    hdr = [c.strip() for c in rows[hi]]
    ct_id = {int(m.group(1)): i for i, h in enumerate(hdr) if (m := CT_ID.match(h))}
    ct_lab = {int(m.group(1)): i for i, h in enumerate(hdr) if (m := CT_LAB.match(h))}
    b_lab = [i for i, h in enumerate(hdr) if B_ANY.match(h)]
    out = []
    for r in rows[hi + 1:]:
        if not any(c.strip() for c in r):
            continue
        deepest = None
        for n in sorted(ct_id, reverse=True):                 # deepest populated CT level
            i = ct_id[n]
            if i < len(r) and r[i].strip().startswith("CL:"):
                deepest = (r[i].strip(),
                           r[ct_lab[n]].strip() if ct_lab.get(n, 10**9) < len(r) else "")
                break
        if not deepest:
            continue
        genes = set()
        for i in b_lab:
            if i < len(r):
                g = r[i].strip()
                # gene SYMBOLS only: reject prose, keep HGNC-style tokens
                if g and SYMBOL.fullmatch(g):
                    genes.add(g)
        if genes:
            out.append((deepest[0], deepest[1], genes))
    return out


if __name__ == "__main__":
    d = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_DIR
    files = sorted(glob.glob(os.path.join(d, "*.csv")))
    if not files:
        raise SystemExit("no CSVs in %s" % d)
    agg = defaultdict(lambda: {"label": "", "genes": set(), "organs": set(), "n_rows": 0})
    for f in files:
        organ = os.path.basename(f).replace("asct-b-", "").replace("vh-", "").replace(".csv", "")
        recs = parse(f)
        for curie, label, genes in recs:
            a = agg[curie]
            a["label"] = a["label"] or label
            a["genes"] |= genes
            a["organs"].add(organ)
            a["n_rows"] += 1
        print("  %-34s %4d rows with CL+markers, %3d distinct CL terms"
              % (organ, len(recs), len({x[0] for x in recs})))
    out = {k: {"label": v["label"], "genes": sorted(v["genes"]),
               "organs": sorted(v["organs"]), "n_rows": v["n_rows"]}
           for k, v in sorted(agg.items())}
    json.dump(out, open(os.path.join(RES, "asctb_markers.json"), "w"), indent=1)
    ng = sum(len(v["genes"]) for v in out.values())
    print("\nASCT+B: %d CL terms with expert-curated markers, %d (term,gene) pairs, %d distinct genes"
          % (len(out), ng, len({g for v in out.values() for g in v["genes"]})))
    print("wrote results/asctb_markers.json")
