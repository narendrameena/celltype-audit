"""The Cell Ontology, as this package uses it: an is_a graph and lineage anchor sets.

Lineage is read off CL's own is_a hierarchy, never guessed from a label string. That
distinction is the whole reason this module exists: a keyword classifier over CL term
names is wrong on terms whose names mislead ('pericyte' is connective tissue in CL, not
muscle) and blind to terms whose names carry no lineage word at all ('type B pancreatic
cell'). Measured over 26 common terms it was wrong on 3 and blind on 3, and because the
audit flags DISAGREEMENT, every one of those becomes a false flag.

Anchor sets are compared as SETS, never as single labels, because CL is a polyhierarchy:
a microglial cell descends from both hematopoietic cell and neural cell, and must conflict
with neither.
"""
import json
import os
import urllib.request
from collections import defaultdict, deque

CL_URL = "http://purl.obolibrary.org/obo/cl.json"

#: CL classes used as lineage anchors. Chosen to be the coarse divisions a marker set can
#: actually distinguish; finer ones are left to the ontology graph.
ANCHORS = {
    "CL:0000988": "haematopoietic",
    "CL:0000066": "epithelial",
    "CL:0000115": "endothelial",
    "CL:0002320": "connective",
    "CL:0000187": "muscle",
    "CL:0002319": "neural",
    "CL:0000586": "germ",
    "CL:0000039": "germline",
}


def default_cache():
    return os.environ.get("CELLTYPE_AUDIT_CACHE",
                          os.path.join(os.path.expanduser("~"), ".cache", "celltype-audit"))


def _curie(u):
    return u.rsplit("/", 1)[-1].replace("_", ":") if "/" in u else u


class Ontology:
    """CL is_a graph, label/synonym index, and lineage anchor sets.

    >>> o = Ontology.load()                      # doctest: +SKIP
    >>> sorted(o.anchors("CL:0000057"))          # doctest: +SKIP
    ['connective']
    """

    def __init__(self, parents, labels, synonyms, version=None):
        self.parents = parents
        self.labels = labels
        self.synonyms = synonyms
        self.version = version
        self._anc = {}

    # ------------------------------------------------------------------ loading
    @classmethod
    def load(cls, path=None, cache=None, download=True):
        """Load CL from a local obographs JSON, or fetch and cache it."""
        cache = cache or default_cache()
        path = path or os.path.join(cache, "cl.json")
        if not os.path.exists(path):
            if not download:
                raise FileNotFoundError(path)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            urllib.request.urlretrieve(CL_URL, path)
        return cls.from_json(json.load(open(path)))

    @classmethod
    def from_json(cls, doc):
        par, lab, syn = defaultdict(set), {}, defaultdict(set)
        version = None
        for g in doc.get("graphs", []):
            meta = g.get("meta") or {}
            for bp in meta.get("basicPropertyValues", []) or []:
                if "version" in str(bp.get("pred", "")).lower():
                    version = bp.get("val")
            for n in g.get("nodes", []):
                c = _curie(n.get("id", ""))
                if not c.startswith("CL:"):
                    continue
                if n.get("lbl"):
                    lab[c] = n["lbl"]
                    syn[n["lbl"].lower()].add(c)
                for s in ((n.get("meta") or {}).get("synonyms") or []):
                    if s.get("val"):
                        syn[s["val"].lower()].add(c)
            for e in g.get("edges", []):
                if e.get("pred") in ("is_a", "rdfs:subClassOf"):
                    s, o = _curie(e.get("sub", "")), _curie(e.get("obj", ""))
                    if s.startswith("CL:") and o.startswith("CL:"):
                        par[s].add(o)
        return cls(dict(par), lab, dict(syn), version)

    # ------------------------------------------------------------------ queries
    def ancestors(self, curie):
        """All is_a ancestors, inclusive of the term itself."""
        if curie in self._anc:
            return self._anc[curie]
        seen, q = {curie}, deque([curie])
        while q:
            x = q.popleft()
            for p in self.parents.get(x, ()):
                if p not in seen:
                    seen.add(p)
                    q.append(p)
        self._anc[curie] = seen
        return seen

    def anchors(self, curie):
        """Lineage anchor set. May contain several members -- CL is a polyhierarchy."""
        return frozenset(ANCHORS[a] for a in self.ancestors(curie) & set(ANCHORS))

    def conflict(self, a, b):
        """True when two terms' anchor sets are DISJOINT and both are non-empty.

        An empty anchor set asserts nothing, which is the correct reading for terms such
        as 'Schwann cell' or 'melanocyte' that sit outside every anchor.
        """
        A, B = self.anchors(a), self.anchors(b)
        return bool(A) and bool(B) and not (A & B)

    def related(self, a, b):
        """is_a-related in either direction: a coarser or finer name for the same thing."""
        return a == b or b in self.ancestors(a) or a in self.ancestors(b)

    def label(self, curie):
        return self.labels.get(curie, curie)

    def __repr__(self):
        return "<Ontology %s terms, release %s>" % (len(self.labels), self.version)
