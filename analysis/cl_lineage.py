#!/usr/bin/env python3
"""Ontology-grounded lineage oracle.

Replaces the regex keyword classifier in lineage.py. A cell type's lineage is not
guessed from its label string but read off CL's own is_a hierarchy: the set of
LINEAGE ANCHOR terms it descends from. Polyhierarchy is handled natively -- in CL
an endothelial cell IS an epithelial cell, so anchor sets overlap and are compared
as SETS, never as single labels. Two terms conflict only when their anchor sets are
disjoint, which is the closest thing CL supports to an assertion of distinctness.
"""
import json
import os
from collections import defaultdict, deque

HERE = os.path.dirname(os.path.abspath(__file__))
CLJSON = os.path.join(HERE, "..", "cl-full.json")

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

_G = {}


def _curie(u):
    return u.rsplit("/", 1)[-1].replace("_", ":") if "/" in u else u


def load():
    """Build is_a parent map + label/synonym index once."""
    if _G:
        return _G
    d = json.load(open(CLJSON))
    par = defaultdict(set)
    lab, syn = {}, defaultdict(set)
    for g in d.get("graphs", []):
        for n in g.get("nodes", []):
            c = _curie(n.get("id", ""))
            if not c.startswith("CL:"):
                continue
            m = n.get("meta") or {}
            if n.get("lbl"):
                lab[c] = n["lbl"]
                syn[n["lbl"].lower()].add(c)
            for s in m.get("synonyms", []) or []:
                if s.get("val"):
                    syn[s["val"].lower()].add(c)
        for e in g.get("edges", []):
            if e.get("pred") in ("is_a", "rdfs:subClassOf"):
                s, o = _curie(e.get("sub", "")), _curie(e.get("obj", ""))
                if s.startswith("CL:") and o.startswith("CL:"):
                    par[s].add(o)
    _G.update(parents=par, label=lab, syn=syn)
    return _G


def ancestors(c):
    """All is_a ancestors of c, inclusive."""
    g = load()
    seen, q = {c}, deque([c])
    while q:
        x = q.popleft()
        for p in g["parents"].get(x, ()):
            if p not in seen:
                seen.add(p)
                q.append(p)
    return seen


def anchor_set(c):
    """Which lineage anchors this CL term descends from (may be several)."""
    return frozenset(ANCHORS[a] for a in ancestors(c) & set(ANCHORS))


def resolve(label):
    """Exact label/synonym -> CL term(s). Returns a set (may be empty/ambiguous)."""
    g = load()
    s = (label or "").strip().lower()
    return set(g["syn"].get(s, ()))


def conflict(a, b):
    """True when two CL terms' anchor sets are disjoint AND both are non-empty."""
    A, B = anchor_set(a), anchor_set(b)
    return bool(A) and bool(B) and not (A & B)


if __name__ == "__main__":
    g = load()
    print("CL terms with is_a parents : %d" % len(g["parents"]))
    print("labels+synonyms indexed    : %d\n" % len(g["syn"]))
    for lbl in ["fibroblast", "type B pancreatic cell", "astrocyte", "macrophage",
                "T cell", "endothelial cell", "neutrophil", "Leydig cell"]:
        cs = resolve(lbl)
        for c in sorted(cs)[:1]:
            print("%-24s %-12s %s" % (lbl, c, sorted(anchor_set(c)) or "-"))
    print()
    for a, b in [("CL:0000057", "CL:0000775"), ("CL:0000169", "CL:0000171"),
                 ("CL:0000115", "CL:0000066"), ("CL:0000540", "CL:0000235")]:
        print("conflict(%s,%s) = %s" % (g["label"].get(a, a)[:22], g["label"].get(b, b)[:22], conflict(a, b)))
