#!/usr/bin/env python3
"""Does scoring in the DISCRIMINATIVE SUBSPACE beat the 5-marker production scorer?

Production: score(CL term) = mean over the type's 5 markers of (pc x me) in the reference.
SUBSPACE  : build the union of ALL types' top-k markers for the organ; represent both the query
            (its own marker detection rates) and every candidate CL term in that subspace; rank
            by cosine.

Rationale: the full 1,958-gene profile is dominated by abundance/housekeeping and scores WORSE
than 5 markers (pancreas 68.2% vs 81.8% top-1). Five markers is too few to be robust. The union
of everyone's discriminative markers is the space where cell identity actually lives.

Reads the WMG cache ONCE for all tissues (it was 4 rescans of 3.3GB before).

Usage: python benchmark/test_subspace.py [topk]
"""
import glob
import json
import os
import sys
import urllib.request
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from scoring_variants import ok                                    # noqa: E402
from heca_markers import marker_table                              # noqa: E402

RES = os.path.join(HERE, "results")
TOPK = int(sys.argv[1]) if len(sys.argv) > 1 else 20
# Every curated organ, with its tissue id read from the atlas mapping. This was the four
# organs that existed when the experiment was written, and the comparison it reports --
# the subspace scorer beating production by 8.4 points -- had never been tried on the six
# curated since. A scorer chosen on four organs is a scorer chosen on four organs.
def _org():
    from gold_organs import curated
    out = {}
    for o in curated():
        m = os.path.join(RES, "heca_to_cl_%s.json" % o)
        h5 = os.path.join(os.path.dirname(HERE), "..", "heca_data", "RNA-%s.h5ad" % o)
        if os.path.exists(m) and os.path.exists(h5):
            ub = json.load(open(m)).get("uberon")
            if ub:
                out[o] = (ub, h5)
    return out


ORG = _org()
WANT = {v[0] for v in ORG.values()}

print("scanning WMG cache once for %d tissues ..." % len(WANT), flush=True)
ref = {u: defaultdict(dict) for u in WANT}
for p in glob.glob(os.path.join(HERE, ".wmg_cache", "*.json")):
    try:
        d = json.load(open(p))
    except Exception:
        continue
    es = d.get("expression_summary")
    if not es:
        continue
    for gene, by_t in es.items():
        for u in WANT:
            ct = by_t.get(u)
            if not ct:
                continue
            for c, v in ct.items():
                a = v.get("aggregated") or {}
                if a.get("pc") is not None:
                    ref[u][c][gene] = a["pc"] * a.get("me", 0.0)
ref = {u: {c: v for c, v in r.items() if len(v) >= 100 and c != "CL:0000000"}
       for u, r in ref.items()}
for u, r in ref.items():
    print("   %-18s %d CL terms" % (u, len(r)), flush=True)

dims = json.loads(urllib.request.urlopen(
    "https://api.cellxgene.cziscience.com/wmg/v2/primary_filter_dimensions", timeout=180).read())
E2S = {}
for e in dims["gene_terms"]["NCBITaxon:9606"]:
    for k, v in e.items():
        E2S[k] = v.upper()

print("\n%-13s %5s | %-17s | %-17s" % ("organ", "n", "PRODUCTION (M5)", "SUBSPACE(top%d)" % TOPK))
print("%-13s %5s | %8s %8s | %8s %8s" % ("", "", "top-1", "top-5", "top-1", "top-5"))
print("-" * 64)
tot = [0, 0, 0, 0, 0]
rows = {}
for organ, (ub, h5) in ORG.items():
    r = ref[ub]
    if not r:
        print("  %s: no reference" % organ)
        continue
    GEN = sorted({g for v in r.values() for g in v})
    RT = sorted(r)
    R = np.array([[r[c].get(g, 0.0) for g in GEN] for c in RT])
    sym = [E2S.get(x, "") for x in GEN]
    s2j = {s: j for j, s in enumerate(sym) if s}
    GOLD = {k: v for k, v in json.load(open(os.path.join(HERE, "%s_gold.json" % organ.lower()))).items()
            if not k.startswith("_") and v}
    PROD = json.load(open(os.path.join(RES, "heca_to_cl_%s.json" % organ)))["types"]
    tbl, _ = marker_table(h5, topk=TOPK, min_cells=50, verbose=False)
    uni = [g for g in sorted({m["gene"] for v in tbl.values() for m in v["markers"]}) if g in s2j]
    js = [s2j[g] for g in uni]
    Rs = R[:, js]
    rn = np.linalg.norm(Rs, axis=1) + 1e-9
    a, n = [0, 0, 0, 0], 0
    for t, gd in GOLD.items():
        if t not in PROD or PROD[t]["n_cells"] < 500 or t not in tbl:
            continue
        n += 1
        p = [c["curie"] for c in PROD[t]["cl"]]
        a[0] += ok(p[0], gd)
        a[1] += any(ok(x, gd) for x in p[:5])
        mine = {m["gene"]: m["pc_in"] for m in tbl[t]["markers"]}
        q = np.array([mine.get(g, 0.0) for g in uni])
        s = Rs.dot(q) / (rn * (np.linalg.norm(q) + 1e-9))
        top = [RT[i] for i in np.argsort(-s)[:5]]
        a[2] += ok(top[0], gd)
        a[3] += any(ok(x, gd) for x in top)
    if not n:
        continue
    print("%-13s %5d | %7.1f%% %7.1f%% | %7.1f%% %7.1f%%"
          % (organ, n, 100 * a[0] / n, 100 * a[1] / n, 100 * a[2] / n, 100 * a[3] / n), flush=True)
    rows[organ] = {"n": n, "prod_top1": round(100 * a[0] / n, 1), "prod_top5": round(100 * a[1] / n, 1),
                   "sub_top1": round(100 * a[2] / n, 1), "sub_top5": round(100 * a[3] / n, 1)}
    for i in range(4):
        tot[i] += a[i]
    tot[4] += n
n = tot[4]
if n:
    print("-" * 64)
    print("%-13s %5d | %7.1f%% %7.1f%% | %7.1f%% %7.1f%%"
          % ("POOLED", n, 100 * tot[0] / n, 100 * tot[1] / n, 100 * tot[2] / n, 100 * tot[3] / n))
    rows["POOLED"] = {"n": n, "prod_top1": round(100 * tot[0] / n, 1), "prod_top5": round(100 * tot[1] / n, 1),
                      "sub_top1": round(100 * tot[2] / n, 1), "sub_top5": round(100 * tot[3] / n, 1)}
json.dump({"topk": TOPK, "rows": rows}, open(os.path.join(RES, "subspace_test.json"), "w"), indent=1)
print("\nwrote results/subspace_test.json")
