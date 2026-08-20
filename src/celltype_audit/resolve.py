"""Atlas label -> Cell Ontology term.

Everything downstream needs one unambiguous CL term per label, and plain exact-string
lookup finds only about half of them. The misses are not exotic; they are conventions:

    'Haematopoietic stem cell'   CL indexes only the -e- spelling
    'Cardiomyocyte cell'         'cardiomyocyte' is already a synonym; suffix doubled
    'Myofibroblast'              CL's label is 'myofibroblast cell'; suffix missing
    'CD8 T cell'                 CL spells it 'CD8-positive, alpha-beta T cell'
    'Alpha cell'                 matches ONE CL term, and it is a MOUSE retinal neuron

so this is a normalisation cascade, then a disambiguation rule, then a quarantine for
labels that are not cell types at all.

The resolver translates the label FAITHFULLY. It deliberately does not second-guess a
label that looks wrong for its organ: tested on seven hand-curated organs, a rule that
abstained whenever the term was not observed in that tissue cost 46 extra abstentions for
3 points of precision, and removed exactly the mis-annotations the audit exists to find.
"""
import re

#: Labels naming a structure, a state, or a bucket -- not a cell type.
QUARANTINE = {
    "unknown", "unclassified", "unassigned", "other", "others", "doublet", "doublets",
    "mixed", "mixture", "ambiguous", "low quality", "lowquality", "debris", "nan", "na",
    "none", "contamination", "artefact", "artifact",
}

#: Cell-STATE qualifiers. CL types the cell, not the state, so these are stripped.
STATE = ("activated", "committed", "proliferating", "proliferative", "cycling",
         "differentiating", "quiescent", "resting", "early", "late", "primed")

#: Abbreviations and lab shorthand CL does not carry as synonyms.
ALIAS = {
    "cd4 t cell": "CD4-positive, alpha-beta T cell",
    "cd8 t cell": "CD8-positive, alpha-beta T cell",
    "naive cd4 t cell": "naive thymus-derived CD4-positive, alpha-beta T cell",
    "naive cd8 t cell": "naive thymus-derived CD8-positive, alpha-beta T cell",
    "memory cd4 t cell": "central memory CD4-positive, alpha-beta T cell",
    "memory cd8 t cell": "central memory CD8-positive, alpha-beta T cell",
    "cd4 treg": "CD4-positive, CD25-positive, alpha-beta regulatory T cell",
    "treg": "regulatory T cell", "nkt cell": "mature NK T cell",
    "nk t cell": "mature NK T cell", "nk cell": "natural killer cell",
    "ilc": "innate lymphoid cell", "cytotoxic t cell": "CD8-positive, alpha-beta cytotoxic T cell",
    "neutrophilic granulocyte": "neutrophil", "eosinophilic granulocyte": "eosinophil",
    "basophilic granulocyte": "basophil", "cardiomyocyte cell": "cardiac muscle cell",
    "atrial cardiomyocyte cell": "regular atrial cardiac myocyte",
    "ventricle cardiomyocyte cell": "regular ventricular cardiac myocyte",
    "excitatory neuron": "glutamatergic neuron", "inhibitory neuron": "GABAergic neuron",
    "adipocyte stem cell": "preadipocyte", "type i alveolar cell": "type I pneumocyte",
    "type ii alveolar cell": "type II pneumocyte", "ductal cell": "duct epithelial cell",
    "spermatogonia": "spermatogonium", "monocyte-derived macrophage": "macrophage",
    "haematopoietic stem and progenitor cell": "hematopoietic precursor cell",
    "hematopoietic stem and progenitor cell": "hematopoietic precursor cell",
    "proximal convoluted tubule": "epithelial cell of proximal tubule",
    "loop of henle": "kidney loop of Henle epithelial cell",
    "collecting duct": "kidney collecting duct epithelial cell",
}

#: Atlas labels drop the organ CL keeps: 'Alpha cell' vs CL's 'pancreatic alpha cell'.
ORGAN_ADJ = {
    "pancreas": ("pancreatic",), "kidney": ("kidney", "renal"),
    "liver": ("liver", "hepatic"), "lung": ("lung", "pulmonary"),
    "heart": ("cardiac", "heart"), "brain": ("brain",), "skin": ("skin",),
    "bone_marrow": ("bone marrow",), "blood": ("blood",), "thymus": ("thymic",),
    "intestine": ("intestinal",), "stomach": ("stomach", "gastric"),
    "testis": ("testicular",), "ovary": ("ovarian",), "prostate": ("prostate",),
    "bladder": ("bladder",), "spleen": ("splenic",), "placenta": ("placental",),
    "breast": ("breast", "mammary"), "uterus": ("uterine",),
}

_AE = re.compile(r"ae")
#: CL carries 152 mouse-only terms plus 'sensu <clade>' classes and obsoletes. They are
#: valid CL classes that must never be returned for human data.
_NOT_HUMAN = re.compile(r"\((?:Mmus|Rnor|Dmel|Cele|Xlae|Drer|sensu [^)]+)\)", re.I)


def _norm(s):
    return re.sub(r"\s+", " ", (s or "").strip().lower())


class Resolver:
    """Map atlas cell-type labels onto CL terms.

    >>> r = Resolver(ontology)                                    # doctest: +SKIP
    >>> r.resolve("Myofibroblast")                                # doctest: +SKIP
    ('CL:0000186', 'normalised')
    """

    #: statuses that give an EXACT identity; 'generalised' is deliberately absent
    EXACT = ("exact", "normalised", "alias", "general", "organ", "organ-qualified")

    def __init__(self, ontology):
        self.o = ontology

    # ------------------------------------------------------------------ helpers
    def _admissible(self, c):
        lab = self.o.labels.get(c, "")
        return bool(lab) and not _NOT_HUMAN.search(lab) and not lab.lower().startswith("obsolete")

    def _variants(self, label):
        s = _norm(label)
        out = [s]
        if s in ALIAS:
            out.append(_norm(ALIAS[s]))
        b = _AE.sub("e", s)
        if b != s:
            out.append(b)
        for v in list(out):
            w = v.split()
            while w and w[0] in STATE:
                w = w[1:]
            if w and " ".join(w) != v:
                out.append(" ".join(w))
        for v in list(out):
            m = re.match(r"^(.*?)[ -]*(?:\d+|i{1,3}|iv)$", v)
            if m and len(m.group(1).strip()) > 4:
                out.append(m.group(1).strip())
        for v in list(out):
            out.append(v[:-5] if v.endswith(" cell") else v + " cell")
            if v.endswith("s"):
                out.append(v[:-1])
        seen, uniq = set(), []
        for v in out:
            if v and v not in seen:
                seen.add(v)
                uniq.append(v)
        return uniq

    def _primary(self, hits, v):
        """A term whose PRIMARY label is the query beats one merely listing it as a
        synonym: 'macrophage' is the label of CL:0000235 but a Drosophila synonym of
        CL:0000394 plasmatocyte."""
        p = [c for c in hits if _norm(self.o.labels.get(c, "")) == v]
        return p[0] if len(p) == 1 else None

    def _most_general(self, hits):
        for c in hits:
            if all(c in self.o.ancestors(o) for o in hits if o != c):
                return c
        return None

    def _lcs(self, hits):
        """Least common subsumer: turns an abstention into a correct-but-general call."""
        common = None
        for c in hits:
            a = self.o.ancestors(c)
            common = a if common is None else (common & a)
        common = (common or set()) - set(hits)
        if not common:
            return None
        best = [c for c in common if not any(c in self.o.ancestors(o) for o in common if o != c)]
        return min(sorted(best), key=lambda c: len(self.o.ancestors(c))) if best else None

    # ------------------------------------------------------------------ resolve
    def resolve(self, label, organ_terms=None, organ=None):
        """-> (curie|None, status).

        status: exact | normalised | alias | general | organ | organ-qualified |
                generalised | ambiguous | absent | quarantined
        """
        s = _norm(label)
        if s in QUARANTINE:
            return None, "quarantined"
        syn = self.o.synonyms
        okey = _norm(organ).replace(" ", "_") if organ else None

        # An organ-qualified form beats a bare match, because the atlas label has simply
        # dropped an organ that CL keeps in the term name. 'Delta cell' otherwise
        # normalises onto 'retinal ganglion cell C1'; with organ='Pancreas' the qualified
        # name 'pancreatic delta cell' exists and is unambiguous. Where a reference is
        # available its observed terms narrow the hits further, but the qualified name is
        # sufficient evidence on its own.
        if okey:
            for adj in ORGAN_ADJ.get(okey, ()):
                for cand in ("%s %s" % (adj, s),
                             ("%s %s" % (adj, s[:-5])) if s.endswith(" cell") else None):
                    if not cand:
                        continue
                    hits = {c for c in syn.get(_norm(cand), ()) if self._admissible(c)}
                    if organ_terms:
                        inorg = hits & set(organ_terms)
                        if inorg:
                            hits = inorg
                    pick = (next(iter(hits)) if len(hits) == 1
                            else (self._primary(hits, _norm(cand)) or self._most_general(hits)
                                  if hits else None))
                    if pick:
                        return pick, "organ-qualified"

        for i, v in enumerate(self._variants(label)):
            hits = {c for c in syn.get(v, ()) if self._admissible(c)}
            if not hits:
                continue
            if organ_terms:
                inorg = hits & set(organ_terms)
                if inorg:
                    hits = inorg
            how = "exact" if i == 0 else ("alias" if _norm(ALIAS.get(s, "")) == v else "normalised")
            if len(hits) == 1:
                return next(iter(hits)), how
            pick = self._primary(hits, v)
            if pick:
                return pick, how
            gen = self._most_general(hits)
            if gen:
                return gen, "general"
            if organ_terms:
                inorg = [c for c in hits if c in organ_terms]
                if len(inorg) == 1:
                    return inorg[0], "organ"
            sub = self._lcs(hits)
            if sub:
                return sub, "generalised"
            return None, "ambiguous"

        # last resort: drop leading modifiers to the head noun. Yields a genuine CL
        # SUPERCLASS, so it is reported as 'generalised' and callers needing exact
        # identity (refinement) must not accept it.
        w = _norm(label).split()
        for i in range(1, max(1, len(w) - 1)):
            for tail in (" ".join(w[i:]), " ".join(w[i:]) + " cell"):
                hits = {c for c in syn.get(tail, ()) if self._admissible(c)}
                if not hits:
                    continue
                if len(hits) == 1:
                    return next(iter(hits)), "generalised"
                pick = self._primary(hits, tail) or self._most_general(hits)
                if pick:
                    return pick, "generalised"
        return None, "absent"

    def is_exact(self, status):
        return status in self.EXACT

    def dead_aliases(self):
        """Every ALIAS target must itself resolve, or the alias is dead weight."""
        return [(k, v) for k, v in ALIAS.items() if not self.o.synonyms.get(_norm(v))]
