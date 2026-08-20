#!/usr/bin/env python3
"""Atlas label -> Cell Ontology term. Identity resolution, stage 1 of the pipeline.

Everything downstream (assignment, contradiction, refinement) needs one unambiguous
CL term per atlas label. Raw exact-string lookup gets only about half of them, and the
misses are not exotic -- they are spelling and word-order conventions:

  'Haematopoietic stem cell'  CL indexes only the -e- spelling
  'Cardiomyocyte cell'        'cardiomyocyte' is already a synonym; the suffix is doubled
  'Myofibroblast'             CL's label is 'myofibroblast cell'; the suffix is missing
  'CD8 T cell'                CL spells it 'CD8-positive, alpha-beta T cell'
  'Proximal convoluted tubule' an anatomical structure, not a cell type at all

so the resolver is a normalisation cascade, then a disambiguation rule for labels that
match several terms, then a quarantine for labels that are not cell types.

DISAMBIGUATION. A label matching several CL terms is usually a general term colliding
with its own specialisations ('Monocyte' -> monocyte + CD14-positive monocyte). When one
hit is an is_a ancestor of every other hit, that hit is the intended reading and is
returned. Organ context breaks the remaining ties ('Beta cell' is the pancreatic B cell
in pancreas, the pars distalis basophil in pituitary). Anything still ambiguous abstains
rather than guessing -- an unresolved label is recoverable, a wrong one is not.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cl_lineage import load, ancestors                                # noqa: E402

# Labels that name an anatomical structure, a state, or a bucket -- not a cell type.
QUARANTINE = {
    "unknown", "unclassified", "unassigned", "other", "others", "doublet", "doublets",
    "mixed", "mixture", "ambiguous", "low quality", "lowquality", "debris", "nan", "na",
    "none", "contamination", "artefact", "artifact",
}

# Cell-state and qualifier words. A cluster called 'Activated fibroblast' is a fibroblast
# in a state; CL types the cell, not the state, so these are dropped to reach the type.
STATE = ("activated", "committed", "proliferating", "proliferative", "cycling",
         "differentiating", "quiescent", "resting", "early", "late", "primed")

# Abbreviations and lab shorthand CL does not carry as synonyms. Only entries that a
# curator would sign off; anything requiring a judgement call is left to abstain.
ALIAS = {
    "cd4 t cell": "CD4-positive, alpha-beta T cell",
    "cd8 t cell": "CD8-positive, alpha-beta T cell",
    "naive cd4 t cell": "naive thymus-derived CD4-positive, alpha-beta T cell",
    "naive cd8 t cell": "naive thymus-derived CD8-positive, alpha-beta T cell",
    "memory cd4 t cell": "central memory CD4-positive, alpha-beta T cell",
    "memory cd8 t cell": "central memory CD8-positive, alpha-beta T cell",
    "effector cd8 t cell": "effector CD8-positive, alpha-beta T cell",
    "cd4 treg": "CD4-positive, CD25-positive, alpha-beta regulatory T cell",
    "treg": "regulatory T cell",
    "nkt cell": "mature NK T cell",
    "nk cell": "natural killer cell",
    "cytotoxic t cell": "CD8-positive, alpha-beta cytotoxic T cell",
    "neutrophilic granulocyte": "neutrophil",
    "eosinophilic granulocyte": "eosinophil",
    "basophilic granulocyte": "basophil",
    "haematopoietic stem and progenitor cell": "hematopoietic precursor cell",
    "hematopoietic stem and progenitor cell": "hematopoietic precursor cell",
    "cardiomyocyte cell": "cardiac muscle cell",
    "atrial cardiomyocyte cell": "regular atrial cardiac myocyte",
    "ventricle cardiomyocyte cell": "regular ventricular cardiac myocyte",
    "ventricular cardiomyocyte cell": "regular ventricular cardiac myocyte",
    "excitatory neuron": "glutamatergic neuron",
    "inhibitory neuron": "GABAergic neuron",
    "adipocyte stem cell": "preadipocyte",
    "epithelial progenitor cell": "epithelial cell",
    "artery endothelial cell": "endothelial cell of artery",
    "vein endothelial cell": "vein endothelial cell",
    "lymphoid cell": "lymphocyte",
    "myeloid cell": "myeloid cell",
    "spermatogonia": "spermatogonium",
    "cd8 trm": "CD8-positive, alpha-beta memory T cell",
    "megakaryocyte-erythrocyte progenitor (mep)": "megakaryocyte-erythroid progenitor cell",
    "megakaryocyte-erythrocyte progenitor": "megakaryocyte-erythroid progenitor cell",
    "monocyte-derived macrophage": "macrophage",
    "enterocyte progenitor": "enterocyte",
    "nk t cell": "mature NK T cell",
    "ilc": "innate lymphoid cell",
    "type i alveolar cell": "type I pneumocyte",
    "type ii alveolar cell": "type II pneumocyte",
    "type 1 alveolar cell": "type I pneumocyte",
    "type 2 alveolar cell": "type II pneumocyte",
    "ductal cell": "duct epithelial cell",
    "myoid cell": "peritubular myoid cell",
    "sinusoidal endothelial cell": "endothelial cell of sinusoid",
    "steroidogenic cortical cell": "adrenal cortex cell",
    "granulocyte-monocyte progenitor (gmp)": "granulocyte monocyte progenitor cell",
    "granulocyte-monocyte progenitor": "granulocyte monocyte progenitor cell",
    # hECA names some clusters after the structure; the cluster is that structure's cells
    "proximal convoluted tubule": "epithelial cell of proximal tubule",
    "distal convoluted tubule": "kidney distal convoluted tubule epithelial cell",
    "loop of henle": "kidney loop of Henle epithelial cell",
    "collecting duct": "kidney collecting duct epithelial cell",
}

_AE = re.compile(r"ae")

# CL carries terms scoped to other species (152 mouse-only, plus 'sensu ...' clades) and
# obsolete terms. They are legitimate CL classes but must never be returned for human
# atlas data: 'Alpha cell' has exactly ONE synonym hit in CL and it is
# 'alpha retinal ganglion cell (Mmus)', which the resolver used to accept as exact.
_NOT_HUMAN = re.compile(r"\((?:Mmus|Rnor|Dmel|Cele|Xlae|Drer|sensu [^)]+)\)", re.I)


def _admissible(c):
    lab = load()["label"].get(c, "")
    return bool(lab) and not _NOT_HUMAN.search(lab) and not lab.lower().startswith("obsolete")


# Atlas labels routinely drop the organ that CL keeps: hECA writes "Alpha cell" where CL
# has only "pancreatic alpha cell". The organ is known to the caller, so the qualified form
# is tried before giving up. Adjectives are the ones CL itself uses in term names.
ORGAN_ADJ = {
    "Pancreas": ("pancreatic",), "Kidney": ("kidney", "renal"),
    "Liver": ("liver", "hepatic"), "Lung": ("lung", "pulmonary"),
    "Heart": ("cardiac", "heart"), "Brain": ("brain",), "Skin": ("skin",),
    "Bone_marrow": ("bone marrow",), "Blood": ("blood",), "Thymus": ("thymic",),
    "Intestine": ("intestinal",), "Stomach": ("stomach", "gastric"),
    "Testis": ("testicular",), "Ovary": ("ovarian",), "Prostate": ("prostate",),
    "Bladder": ("bladder",), "Spleen": ("splenic",), "Placenta": ("placental",),
    "Breast": ("breast", "mammary"), "Uterus": ("uterine",),
}


def _norm(s):
    return re.sub(r"\s+", " ", (s or "").strip().lower())


def _variants(label):
    """Normalisation cascade, most to least conservative."""
    s = _norm(label)
    out = [s]
    if s in ALIAS:
        out.append(_norm(ALIAS[s]))
    b = _AE.sub("e", s)                      # haematopoietic -> hematopoietic
    if b != s:
        out.append(b)
    # drop cell-STATE qualifiers to reach the underlying type
    for v in list(out):
        w = v.split()
        while w and w[0] in STATE:
            w = w[1:]
        if w and " ".join(w) != v:
            out.append(" ".join(w))
    # trailing subset digits: 'conventional dendritic cell 2' -> 'conventional dendritic cell'
    for v in list(out):
        m = re.match(r"^(.*?)[ -]*(?:\d+|i{1,3}|iv)$", v)
        if m and m.group(1).strip() and len(m.group(1).strip()) > 4:
            out.append(m.group(1).strip())
    for v in list(out):
        if v.endswith(" cell"):
            out.append(v[:-5])               # 'cardiomyocyte cell' -> 'cardiomyocyte'
        else:
            out.append(v + " cell")          # 'myofibroblast' -> 'myofibroblast cell'
        if v.endswith("s"):
            out.append(v[:-1])               # crude plural
    seen, uniq = set(), []
    for v in out:
        if v and v not in seen:
            seen.add(v)
            uniq.append(v)
    return uniq


def _primary(hits, v):
    """A term whose PRIMARY label is the query beats one that merely lists it as a
    synonym -- 'macrophage' is the label of CL:0000235 but a Drosophila synonym of
    CL:0000394 plasmatocyte."""
    g = load()
    p = [c for c in hits if _norm(g["label"].get(c, "")) == v]
    return p[0] if len(p) == 1 else None


def _most_general(hits):
    """If one hit is an is_a ancestor of all others, it is the intended reading."""
    hits = list(hits)
    for c in hits:
        if all(c in ancestors(o) for o in hits if o != c):
            return c
    return None


def _lcs(hits):
    """Least common subsumer: the most specific CL term subsuming every hit. Turns an
    abstention into a correct-but-general call -- 'Alveolar fibroblast' matches two
    sibling subtypes, and their subsumer is the fibroblast class both belong to."""
    common = None
    for c in hits:
        a = ancestors(c)
        common = a if common is None else (common & a)
    common = (common or set()) - set(hits)
    if not common:
        return None
    # most specific = the one no other member subsumes. CL is a polyhierarchy, so two
    # subsumers can tie without either subsuming the other; both are then valid
    # superclasses of every hit, and the more general of the tie is the safer call.
    best = [c for c in common if not any(c in ancestors(o) for o in common if o != c)]
    if not best:
        return None
    return min(sorted(best), key=lambda c: len(ancestors(c)))


def resolve(label, organ_terms=None, organ=None):
    """-> (curie|None, status). status in exact/normalised/alias/general/organ/
    organ-qualified/ambiguous/absent/quarantined."""
    s = _norm(label)
    if s in QUARANTINE:
        return None, "quarantined"
    g = load()
    # If the organ-qualified form exists AND the organ is observed to contain it, that is
    # stronger evidence than any bare-string match, so it is tried first. "Delta cell"
    # otherwise normalises onto "retinal ganglion cell C1" before anything else is checked.
    if organ and organ_terms:
        for adj in ORGAN_ADJ.get(organ, ()):
            for cand in ("%s %s" % (adj, s),
                         ("%s %s" % (adj, s[:-5])) if s.endswith(" cell") else None):
                if not cand:
                    continue
                hits = {c for c in g["syn"].get(_norm(cand), ())
                        if _admissible(c)} & set(organ_terms)
                if len(hits) == 1:
                    return next(iter(hits)), "organ-qualified"
                if hits:
                    pick = _primary(hits, _norm(cand)) or _most_general(hits)
                    if pick:
                        return pick, "organ-qualified"
    for i, v in enumerate(_variants(label)):
        hits = {c for c in g["syn"].get(v, ()) if _admissible(c)}
        if not hits:
            continue
        # Organ context is a VALIDITY CHECK, not only a tie-break. A single hit that the
        # organ has never been observed to contain is worse evidence than no hit at all.
        if organ_terms:
            inorg = hits & set(organ_terms)
            if inorg:
                hits = inorg
        how = "exact" if i == 0 else ("alias" if _norm(ALIAS.get(s, "")) == v else "normalised")
        if len(hits) == 1:
            return next(iter(hits)), how
        pri = _primary(hits, v)
        if pri:
            return pri, how
        gen = _most_general(hits)
        if gen:
            return gen, "general"
        if organ_terms:                       # break the tie with organ context
            inorg = [c for c in hits if c in organ_terms]
            if len(inorg) == 1:
                return inorg[0], "organ"
        sub = _lcs(hits)
        if sub:
            return sub, "generalised"
        return None, "ambiguous"
    # the label may simply have dropped the organ CL keeps in the term name
    if organ:
        for adj in ORGAN_ADJ.get(organ, ()):
            for cand in ("%s %s" % (adj, s), "%s %s" % (adj, s[:-5]) if s.endswith(" cell") else None):
                if not cand:
                    continue
                hits = {c for c in g["syn"].get(_norm(cand), ()) if _admissible(c)}
                if organ_terms:
                    inorg = hits & set(organ_terms)
                    if inorg:
                        hits = inorg
                if len(hits) == 1:
                    return next(iter(hits)), "organ-qualified"
                if hits:
                    pick = _primary(hits, _norm(cand)) or _most_general(hits)
                    if pick:
                        return pick, "organ-qualified"
    # last resort: drop leading modifiers to reach the head noun. This yields a genuine
    # CL superclass rather than the exact type, so it is reported as "generalised" and
    # callers needing exact identity (refinement) must not accept it.
    w = _norm(label).split()
    for i in range(1, max(1, len(w) - 1)):
        for tail in (" ".join(w[i:]), " ".join(w[i:]) + " cell"):
            hits = set(g["syn"].get(tail, ()))
            if not hits:
                continue
            if len(hits) == 1:
                return next(iter(hits)), "generalised"
            pick = _primary(hits, tail) or _most_general(hits)
            if pick:
                return pick, "generalised"
    return None, "absent"


def selftest():
    """Every ALIAS target must itself resolve, or the alias is dead weight."""
    g, dead = load(), []
    for k, v in ALIAS.items():
        if not g["syn"].get(_norm(v)):
            dead.append((k, v))
    return dead


if __name__ == "__main__":
    d = selftest()
    print("dead aliases: %d %s\n" % (len(d), d if d else ""))
    for t in ["Macrophage", "Monocyte", "Schwann cell", "Cardiomyocyte cell",
              "Myofibroblast", "Haematopoietic stem cell", "CD8 T cell", "NKT cell",
              "Neutrophilic granulocyte", "Beta cell", "Proximal convoluted tubule",
              "Naive CD4 T cell", "Excitatory neuron", "Fibroblast",
              "Activated fibroblast", "Conventional dendritic cell 2", "ILC",
              "Basal decidual cell", "Type I alveolar cell", "Vascular epithelial cell",
              "Committed preadipocyte", "Unknown"]:
        c, how = resolve(t)
        g = load()
        print("  %-30s %-12s %-11s %s" % (t, how, c or "-", g["label"].get(c, "") if c else ""))
