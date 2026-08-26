#!/usr/bin/env python3
"""V10 and V11: the two checks that make a proposal falsifiable rather than merely stated.

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

  V11  COUNTEREXAMPLES. A marker-based sufficient condition claims every cluster of this
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
from cl_resolve import resolve as resolve2                              # noqa: E402

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
    """Organ -> UBERON, from the atlas mapping rather than the wide reference.

    The wide reference covers 28 organs; the atlas maps 36. Reading the tissue id from it
    made V10 unrunnable for Adrenal_gland, Pleura, Spinal_cord and Ureter -- reported as
    "no tissue reference" when the reference was never the problem, since the check queries
    WMG directly and only needs the id. heca_to_cl_<organ>.json carries it for every organ
    the atlas maps.
    """
    p = os.path.join(RES, "heca_to_cl_%s.json" % organ)
    if os.path.exists(p):
        ub = json.load(open(p)).get("uberon")
        if ub:
            return ub
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


def v11(term, markers, g):
    """-> (verdict, detail). Clusters ASSERTED to be this type that lack the markers.

    "Asserted" means the atlas LABEL resolves to the term. An earlier version matched on
    heca_to_cl's `cl` field, which is the expression-ranked candidate list -- the opposite
    claim. It reported Skin "Monocyte" as a neutrophil counterexample, when what that row
    records is the audit finding neutrophil to be the best expression match for a cluster
    labelled monocyte. Counting those as violations of a neutrophil condition inverts the
    check: it collects clusters the markers point AT and calls them clusters that fail the
    markers.
    """
    if not markers:
        return "n/a", "no markers to violate"
    viol, checked = [], 0
    for p in sorted(__import__("glob").glob(os.path.join(RES, "heca_to_cl_*.json"))):
        d = json.load(open(p))
        organ = d["organ"]
        try:
            deep = json.load(open(os.path.join(RES, "heca_markers_deep_%s.json" % organ)))["types"]
        except FileNotFoundError:
            continue
        ctx = {c["curie"] for v in d["types"].values() for c in v.get("cl", [])}
        for label, v in d["types"].items():
            if v.get("n_cells", 0) < MINC:
                continue
            asserted, _how = resolve2(label, ctx, organ=organ)
            if asserted != term:
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


def _name_contained(label, rival):
    """The rival's name already contains the label, or the reverse.

    "Luminal cell" against "luminal cell of prostate epithelium" is not two cell types that
    happen to share markers; it is one cell type the resolver failed to reach. Distinguish
    that from "Pit cell" against "goblet cell", which are genuinely different types that
    both make mucus -- shared markers there weaken a proposal without settling it.
    """
    import re as _re
    a = _re.sub(r"[^a-z0-9]+", " ", label.lower()).strip()
    b = _re.sub(r"[^a-z0-9]+", " ", rival.lower()).strip()
    return a in b or b in a or (set(a.split()) and set(a.split()) <= set(b.split()))


def _toks(x):
    """Label as an unordered bag of stems, `cell` dropped -- the same normalisation the
    already-known gate uses, so the two agree about when CL covers a population."""
    import re as _re
    w = _re.sub(r"[^a-z0-9]+", " ", x.lower()).split()
    w = [t[:-1] if len(t) > 3 and t.endswith("s") and not t.endswith("ss") else t for t in w]
    return frozenset(t for t in w if t not in ("cell", "cells"))


def v12(p, g):
    """-> (verdict, detail). Under the improved scorer, is this population already named?

    A proposal says CL has no term for what the atlas calls X. The scorer changed
    underneath it -- the shared-subspace ranking is 14.5 points more accurate against the
    hand-curated gold -- so the claim is re-asked with the better evidence: does the term
    the improved scorer now puts first actually name this population?

    Deterministic on purpose. The obvious form of this check is a threshold on the
    scorer's confidence, and that is measured rather than assumed (`v12_signal.py`, which
    writes results/v12_signal.json over the same 297 hand-curated cell types): the top-1
    cosine separates correct from incorrect calls so poorly that 82% of wrong calls exceed
    the 10th percentile of right ones; the margin over the best unrelated runner-up reaches
    AUC 0.683, weak signal rather than none; and the term's is_a depth reaches 0.509, which
    is chance. Thresholding the margin would buy a little and cost a confident verdict on a
    proposal, which is the trade this stack exists to refuse. Name matching carries no
    threshold.
    """
    top = ((p.get("expression_top") or [{}])[0] or {}).get("curie")
    if p["kind"] != "new-term" or not top:
        return "n/a", "only a new-term proposal claims CL lacks a term for a population"
    name = g["label"].get(top, top)
    if _toks(p["label"]) and _toks(p["label"]) <= _toks(name):
        return "fail", ("the improved scorer puts %s first, which already names this "
                        "population; CL may not be missing a term" % name)
    return "pass", ("the improved scorer puts %s first, which does not name this "
                    "population, so the gap stands under the better evidence" % name)


def main():
    g = load()
    doc = json.load(open(os.path.join(DOCS, "proposals.json")))
    ran = {"V10": {}, "V11": {}, "V12": {}}
    for p in doc["proposals"]:
        mk, ch = p.get("markers") or [], p["checks"]
        term = (p.get("candidate") or {}).get("curie")
        if p["kind"] == "marker-condition" and term:
            ch["V10"], d10 = v10(term, mk, p["organ"], g)
            ch["V11"], d12 = v11(term, mk, g)
        elif p["kind"] == "new-term":
            # There is no term yet, so exclusivity is asked of the population instead:
            # does anything CL already carries in this tissue express all these markers?
            near = (p.get("expression_top") or [{}])[0].get("curie")
            if near and mk:
                ch["V10"], d10 = v10(near, mk, p["organ"], g)
            else:
                ch["V10"], d10 = "n/a", "no scoreable term to compare against"
            ch["V11"], d12 = "n/a", "a missing term is not falsified by a counterexample"
        else:
            ch["V10"], d10 = "n/a", "not a marker-based proposal"
            ch["V11"], d12 = "n/a", ("this proposal asserts no marker condition, so there "
                                     "is nothing a cluster could violate")
        ch["V12"], d12b = v12(p, g)
        p.setdefault("check_detail", {})["V12"] = d12b
        p["check_detail"]["V10"] = d10
        p["check_detail"]["V11"] = d12

        # A proposal that fails exclusivity is not ready to put in front of a curator, and
        # saying so is the point of running the check. Two outcomes, deliberately not one:
        #   * the winning term's name already contains the label -> the population is
        #     covered and the resolver simply missed it; drop it as satisfied
        #   * a different term merely shares the markers -> the proposal is WEAKENED, not
        #     refuted, and is kept where a curator can see why it is doubted
        p["readiness"] = "ready"
        if ch["V10"] == "fail":
            import re as _re
            m = _re.search(r"term in \S+ is (.+?) \([\d,]+ cells\) at ", d10)
            rival = m.group(1) if m else ""
            if rival and _name_contained(p["label"], rival):
                p["readiness"] = "covered"
                p["covered_by"] = rival
            else:
                p["readiness"] = "weakened"
                p["weakened_by"] = rival
        ran["V10"][ch["V10"]] = ran["V10"].get(ch["V10"], 0) + 1
        ran["V11"][ch["V11"]] = ran["V11"].get(ch["V11"], 0) + 1
        ran["V12"][ch["V12"]] = ran["V12"].get(ch["V12"], 0) + 1
        print("  %-16s %-9s %-28s V10=%-7s V11=%s"
              % (p["kind"], p["organ"], p["label"][:28], ch["V10"], ch["V11"]))
    covered = [p for p in doc["proposals"] if p.get("readiness") == "covered"]
    doc["proposals"] = [p for p in doc["proposals"] if p.get("readiness") != "covered"]
    doc.setdefault("covered_after_exclusivity", []).extend(
        {"organ": p["organ"], "label": p["label"], "covered_by": p.get("covered_by")}
        for p in covered)
    r = doc["summary"]
    r["n"] = len(doc["proposals"])
    r["ready"] = sum(1 for p in doc["proposals"] if p.get("readiness") == "ready")
    r["weakened"] = sum(1 for p in doc["proposals"] if p.get("readiness") == "weakened")
    r["covered_after_exclusivity"] = len(doc["covered_after_exclusivity"])
    print("\n  ready %d | weakened by exclusivity %d | dropped as already covered %d"
          % (r["ready"], r["weakened"], r["covered_after_exclusivity"]))
    for c in covered:
        print("     dropped: %-9s %-24s -> %s" % (c["organ"], c["label"][:24], c.get("covered_by")))
    json.dump(doc, open(os.path.join(DOCS, "proposals.json"), "w"), indent=1)
    print("\n  V10 %s\n  V11 %s\n  V12 %s" % (ran["V10"], ran["V11"], ran["V12"]))
    print("wrote %s" % os.path.join(DOCS, "proposals.json"))


if __name__ == "__main__":
    main()
