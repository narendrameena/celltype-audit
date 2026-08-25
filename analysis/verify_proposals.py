#!/usr/bin/env python3
"""V10 and V12: the two checks that make a proposal falsifiable rather than merely stated.

A proposal that carries only V1-V3 says its CURIE is real and not obsolete. That is
bookkeeping. What makes a proposed axiom worth a curator's attention is evidence that
someone tried to break it and reports whether they succeeded:

  V10  EXCLUSIVITY. Do the proposed markers pick this population out, or does something
       else in the same tissue express them too? Computed against the CELLxGENE reference
       for the tissue, over every CL term it carries.

       Ancestors and descendants of the target are excluded from the competitor set, and
       that is not a convenience. Asked naively for neutrophil on FCGR3B+ELANE in lung,
       the top three terms are neutrophil, granulocyte and blood cell -- its own parents,
       which share the profile by construction. A test that counted those as competitors
       would report every well-specified marker set as non-exclusive.

  V12  COUNTEREXAMPLES. A marker-based sufficient condition claims every cluster of this
       type expresses these genes. That is falsifiable, so falsify it: search every
       audited cluster asserted to the term and report the ones that do not.

Both report `n/a` where they do not apply rather than `not-run`, because "this check does
not bear on this kind of proposal" and "we did not do it" are different statements and a
curator is entitled to tell them apart.

The reference is queried for the specific marker genes, not read from the pre-built
matrices: those carry a 2,000-gene marker subset per tissue and do not contain FCGR3B or
ELANE, so the flagship proposal could not have been checked against them at all.

Usage:
    python benchmark/verify_proposals.py          # re-verify docs/proposals.json in place
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "celltype-audit", "src")))
RES = os.path.join(HERE, "results")
DOCS = os.path.abspath(os.path.join(HERE, "..", "..", "celltype-audit", "docs"))

from cl_lineage import load, ancestors                                  # noqa: E402

# A competitor must beat this share of the target's own signal to count as one; below it
# the markers are doing their job. Reported alongside the verdict so it can be argued with.
EXCL_RATIO = 0.50
MINC = 500
# A competitor estimated from a handful of cells can win on noise. Nose "Serous gland
# cell" lost to endothelium at 78% on MUC7, ODAM and LPO -- salivary enzymes an
# endothelial cell does not make, from a profile too thin to mean anything.
MIN_RIVAL_CELLS = 500


def _ref():
    from celltype_audit.reference import Reference
    return Reference


def _uberon(organ):
    idx = json.load(open(os.path.join(RES, "wide_ref_index.json")))
    o = idx["organs"].get(organ)
    return o["uberon"] if isinstance(o, dict) else o


def v10(term, markers, organ, g):
    """-> (verdict, detail). Do these markers pick the term out within its tissue?"""
    ub = _uberon(organ)
    if not ub or not markers:
        return "n/a", "no tissue reference or no markers"
    try:
        ref = _ref().fetch([ub], markers, verbose=False)
    except Exception as exc:
        return "not-run", "reference query failed: %s" % str(exc)[:60]
    M, gi, ti = ref.mats.get(ub), ref.gene_ix.get(ub, {}), ref.term_ix.get(ub, {})
    support = (ref.counts or {}).get(ub, {})
    cols = [gi[m] for m in markers if m in gi]
    if M is None or not cols or term not in ti:
        return "not-run", "term or markers absent from the tissue reference"
    got = [m for m in markers if m in gi]
    strength = M[:, cols].min(axis=1)          # a term must express ALL of them
    inv = {v: k for k, v in ti.items()}
    mine = float(strength[ti[term]])
    kin = set(ancestors(term)) | {term}
    obsolete = {c for c, n in g["label"].items() if str(n).lower().startswith("obsolete")}
    rivals = [(float(strength[i]), inv[i]) for i in range(len(strength))
              # a competitor must be a real, current, adequately-supported alternative:
              # not kin (they share the profile by construction), not obsolete (CL has
              # withdrawn it, so it cannot be what the population is), and not estimated
              # from too few cells to distinguish signal from noise
              if inv[i] not in kin and term not in ancestors(inv[i])
              and inv[i] not in obsolete
              and support.get(inv[i], 0) >= MIN_RIVAL_CELLS]
    if not rivals or mine <= 0:
        return "not-run", "no comparable terms in this tissue"
    best, who = max(rivals)
    ratio = best / mine
    detail = ("markers %s; strongest unrelated, current, well-supported term in %s is %s "
              "(%s cells) at %.0f%% of the target's signal"
              % ("+".join(got), organ, g["label"].get(who, who),
                 format(support.get(who, 0), ","), 100 * ratio))
    return ("pass" if ratio < EXCL_RATIO else "fail"), detail


def v12(term, markers, g):
    """-> (verdict, detail). Clusters asserted to be this type that lack the markers."""
    if not markers:
        return "n/a", "no markers to violate"
    kin = {term}
    viol, checked = [], 0
    for p in sorted(__import__("glob").glob(os.path.join(RES, "heca_to_cl_*.json"))):
        d = json.load(open(p))
        organ = d["organ"]
        try:
            deep = json.load(open(os.path.join(RES, "heca_markers_deep_%s.json" % organ)))["types"]
        except FileNotFoundError:
            continue
        for label, v in d["types"].items():
            if v.get("n_cells", 0) < MINC:
                continue
            if not any(c["curie"] in kin for c in v.get("cl", [])[:1]):
                continue
            checked += 1
            have = {m["gene"] for m in deep.get(label, {}).get("markers", [])[:50]}
            if not (set(markers) & have):
                viol.append((organ, label, v["n_cells"]))
    if not checked:
        return "not-run", "no audited cluster is asserted to this term"
    if not viol:
        return "pass", "%d clusters asserted to this term; none lack the markers" % checked
    viol.sort(key=lambda r: -r[2])
    return "fail", ("%d of %d clusters asserted to this term lack every proposed marker; "
                    "largest is %s %s at %s cells"
                    % (len(viol), checked, viol[0][0], viol[0][1], format(viol[0][2], ",")))


def main():
    g = load()
    doc = json.load(open(os.path.join(DOCS, "proposals.json")))
    ran = {"V10": {}, "V12": {}}
    for p in doc["proposals"]:
        mk, ch = p.get("markers") or [], p["checks"]
        term = (p.get("candidate") or {}).get("curie")
        if p["kind"] == "marker-condition" and term:
            ch["V10"], d10 = v10(term, mk, p["organ"], g)
            ch["V12"], d12 = v12(term, mk, g)
        elif p["kind"] == "new-term":
            # There is no term yet, so exclusivity is asked of the population instead:
            # does anything CL already carries in this tissue express all these markers?
            near = (p.get("expression_top") or [{}])[0].get("curie")
            if near and mk:
                ch["V10"], d10 = v10(near, mk, p["organ"], g)
            else:
                ch["V10"], d10 = "n/a", "no scoreable term to compare against"
            ch["V12"], d12 = "n/a", "a missing term is not falsified by a counterexample"
        else:
            ch["V10"], d10 = "n/a", "not a marker-based proposal"
            ch["V12"] = ch.get("V12", "n/a"); d12 = "structural, not marker-based"
        p.setdefault("check_detail", {})["V10"] = d10
        p["check_detail"]["V12"] = d12
        ran["V10"][ch["V10"]] = ran["V10"].get(ch["V10"], 0) + 1
        ran["V12"][ch["V12"]] = ran["V12"].get(ch["V12"], 0) + 1
        print("  %-16s %-9s %-28s V10=%-7s V12=%s"
              % (p["kind"], p["organ"], p["label"][:28], ch["V10"], ch["V12"]))
    json.dump(doc, open(os.path.join(DOCS, "proposals.json"), "w"), indent=1)
    print("\n  V10 %s\n  V12 %s" % (ran["V10"], ran["V12"]))
    print("wrote %s" % os.path.join(DOCS, "proposals.json"))


if __name__ == "__main__":
    main()
