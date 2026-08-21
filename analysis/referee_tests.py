#!/usr/bin/env python3
"""Two tests a referee will ask for, run before a referee asks.

TEST A -- is the small-cluster enrichment an artefact of reference influence?
  The audited atlases are part of the corpus behind the expression reference, so a large
  cluster helps shape the profile of the term it is then scored against and should be
  harder to flag. If so, the size effect in the epidemiology is manufactured by the
  design rather than observed in the data. Tested by modelling error against cluster size
  and against the reference support for the asserted term, separately and together.

TEST B -- is the flagship lung cluster genuine neutrophils that lost their markers?
  Neutrophils are RNA-poor and fragile, so a degraded neutrophil cluster could lack
  FCGR3B and ELANE for technical rather than biological reasons. Tested by counting genes
  detected per cell (from the CSR indptr, so no expression values are read), and by
  reading what the SOURCE datasets called these cells before harmonisation.

Usage: python benchmark/referee_tests.py
"""
import json
import os
import re
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
LUNG = os.path.join(os.path.dirname(HERE), "..", "heca_data", "RNA-Lung.h5ad")

from cl_resolve import resolve as resolve2                            # noqa: E402
import within_lineage as wl                                           # noqa: E402

MONO = re.compile(r"monocyte|macrophage|\bm[φΦ]|myeloid|\bdc\d?\b|cdc\d?|dendritic"
                  r"|as-dc|promonocyte|gmp\b|\bmono\.", re.I)
NEU = re.compile(r"neutrophil|granulocyte|myelocyte", re.I)


def test_a():
    """Error ~ cluster size vs error ~ reference support, on the epidemiology population."""
    import statsmodels.api as sm
    from scipy.stats import chi2, fisher_exact, spearmanr
    idx, _mats = wl.load_ref()
    rows = []
    for r in json.load(open(os.path.join(RES, "auditor_recall.json"))):
        o = idx["organs"].get(r["organ"])
        if not o:
            continue
        d = json.load(open(os.path.join(RES, "heca_to_cl_%s.json" % r["organ"])))
        ctx = {c["curie"] for v in d["types"].values() for c in v.get("cl", [])}
        cur, _ = resolve2(r["label"], ctx, organ=r["organ"])
        rows.append({"n_cells": r["n_cells"], "error": bool(r["error"]),
                     "ref_support": int(idx["counts"].get(o["uberon"], {}).get(cur, 0)) if cur else 0})
    json.dump(rows, open(os.path.join(RES, "size_vs_support_epi.json"), "w"), indent=1)

    def model(rr, tag):
        n = np.log10([x["n_cells"] for x in rr])
        s = np.log10([max(x["ref_support"], 1) for x in rr])
        y = np.array([1.0 if x["error"] else 0.0 for x in rr])
        f = lambda X: sm.Logit(y, sm.add_constant(np.column_stack(X))).fit(disp=0)
        m1, m2, m3 = f([n]), f([s]), f([n, s])
        rho, prho = spearmanr(n, s)
        print("\n  %s  (n=%d, %d errors)" % (tag, len(rr), int(y.sum())))
        print("    Spearman(log size, log support)   rho = %+.3f   p = %.2g" % (rho, prho))
        print("    error ~ size                      beta = %+.3f   p = %.4f" % (m1.params[1], m1.pvalues[1]))
        print("    error ~ support                   beta = %+.3f   p = %.4f" % (m2.params[1], m2.pvalues[1]))
        print("    error ~ size + support            size p = %.4f | support p = %.4f"
              % (m3.pvalues[1], m3.pvalues[2]))
        print("    LR, adding size to support-only   p = %.4f" % chi2.sf(2 * (m3.llf - m2.llf), 1))
        print("    LR, adding support to size-only   p = %.4f" % chi2.sf(2 * (m3.llf - m1.llf), 1))
        return {"n": len(rr), "errors": int(y.sum()), "rho": rho, "p_rho": prho,
                "p_size": m1.pvalues[1], "p_support": m2.pvalues[1],
                "p_size_adj": m3.pvalues[1], "p_support_adj": m3.pvalues[2],
                "p_lr_size": chi2.sf(2 * (m3.llf - m2.llf), 1)}

    print("TEST A  does reference support explain the small-cluster enrichment?")
    out = {"all": model(rows, "all resolvable curated cell types"),
           "support_gt0": model([r for r in rows if r["ref_support"] > 0],
                                "excluding rows with no reference support")}
    nc = np.array([r["n_cells"] for r in rows]); y = np.array([r["error"] for r in rows])
    sup = np.array([r["ref_support"] for r in rows])
    a = int(y[nc < 2000].sum()); b = int((nc < 2000).sum() - a)
    c = int(y[nc >= 2000].sum()); dd = int((nc >= 2000).sum() - c)
    print("\n  dichotomised at 2,000 cells: small %d/%d vs large %d/%d, Fisher p = %.4f"
          % (a, a + b, c, c + dd, fisher_exact([[a, b], [c, dd]])[1]))
    print("  stratified by reference support (median split):")
    strat = []
    for lab, mask in (("low support", sup <= np.median(sup)), ("high support", sup > np.median(sup))):
        s1 = (nc < 2000) & mask; l1 = (nc >= 2000) & mask
        a2 = int(y[s1].sum()); b2 = int(s1.sum() - a2)
        c2 = int(y[l1].sum()); d2 = int(l1.sum() - c2)
        p = fisher_exact([[a2, b2], [c2, d2]])[1]
        print("    %-13s small %d/%-3d (%5.1f%%)   large %d/%-3d (%5.1f%%)   p = %.3f"
              % (lab, a2, a2 + b2, 100 * a2 / max(a2 + b2, 1), c2, c2 + d2,
                 100 * c2 / max(c2 + d2, 1), p))
        strat.append({"stratum": lab, "small_k": a2, "small_n": a2 + b2,
                      "large_k": c2, "large_n": c2 + d2, "p": p})
    out["dichotomised"] = {"small_k": a, "small_n": a + b, "large_k": c, "large_n": c + dd,
                           "p": fisher_exact([[a, b], [c, dd]])[1]}
    out["stratified"] = strat
    return out


def test_b():
    """Genes detected per cell, and what the source datasets called these cells."""
    import h5py
    print("\n\nTEST B  is the lung 'Neutrophilic granulocyte' cluster RNA-poor?")
    f = h5py.File(LUNG, "r")
    ct = f["obs/cell_type"]
    cats = np.array([c.decode() if isinstance(c, bytes) else c for c in ct["categories"][:]])
    codes = ct["codes"][:].astype(np.int64)
    ngene = np.diff(f["X/indptr"][:]).astype(np.int64)       # nonzeros per cell = genes detected
    per = []
    for i, c in enumerate(cats):
        m = codes == i
        if m.sum() < 500:
            continue
        per.append((c, int(m.sum()), float(np.median(ngene[m]))))
    per.sort(key=lambda r: r[2])
    tgt = next(r for r in per if r[0] == "Neutrophilic granulocyte")
    mono = next((r for r in per if r[0] == "Monocyte"), None)
    print("    genes detected per cell, median")
    print("      the flagged cluster        %5.0f   (rank %d of %d cell types, 1 = lowest)"
          % (tgt[2], per.index(tgt) + 1, len(per)))
    if mono:
        print("      lung 'Monocyte' cluster    %5.0f" % mono[2])
    print("      median across cell types   %5.0f" % np.median([r[2] for r in per]))
    print("      whole lung, per cell       %5.0f" % np.median(ngene))
    print("      lowest three cell types    %s"
          % ", ".join("%s %.0f" % (r[0][:24], r[2]) for r in per[:3]))

    g = f["obs/original_name"]
    ocats = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]])
    ocodes = g["codes"][:].astype(np.int64)
    m = codes == int(np.where(cats == "Neutrophilic granulocyte")[0][0])
    vals, cnts = np.unique(ocodes[m], return_counts=True)
    n = int(m.sum())
    names = {ocats[v]: int(c) for v, c in zip(vals, cnts)}
    mono_n = sum(k for nm, k in names.items() if MONO.search(nm) or not NEU.search(nm) and False)
    mono_n = sum(k for nm, k in names.items() if MONO.search(nm) and not (NEU.search(nm) and not MONO.search(nm)))
    neu_n = sum(k for nm, k in names.items() if NEU.search(nm) and not MONO.search(nm))
    mangled = {nm: k for nm, k in names.items() if "Monocytenocyte" in nm}
    print("\n    what the SOURCE datasets called these %s cells (obs/original_name, %d distinct)"
          % (format(n, ","), len(names)))
    print("      mononuclear phagocyte      %6s  %5.1f%%" % (format(mono_n, ","), 100 * mono_n / n))
    print("      neutrophil / granulocyte   %6s  %5.1f%%" % (format(neu_n, ","), 100 * neu_n / n))
    print("      other or unspecified       %6s  %5.1f%%"
          % (format(n - mono_n - neu_n, ","), 100 * (n - mono_n - neu_n) / n))
    print("      top five: %s" % "; ".join("%s %d" % (k, v) for k, v in
                                           sorted(names.items(), key=lambda x: -x[1])[:5]))
    print("\n      cells whose original name is MANGLED ('Monocytenocyte'): %s in %d strings"
          % (format(sum(mangled.values()), ","), len(mangled)))
    f.close()
    return {"median_genes_flagged": tgt[2], "median_genes_monocyte": mono[2] if mono else None,
            "rank": per.index(tgt) + 1, "n_types": len(per),
            "median_across_types": float(np.median([r[2] for r in per])),
            "n_cells": n, "mononuclear_phagocyte": mono_n, "neutrophil": neu_n,
            "mangled_cells": sum(mangled.values()), "mangled_strings": len(mangled),
            "original_names": names}


if __name__ == "__main__":
    out = {"test_a": test_a(), "test_b": test_b()}
    json.dump(out, open(os.path.join(RES, "referee_tests.json"), "w"), indent=1)
    print("\nwrote results/referee_tests.json")
