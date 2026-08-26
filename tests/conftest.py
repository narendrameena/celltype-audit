"""A tiny synthetic Cell Ontology, so the tests need no network and no 50 MB download.

It reproduces the structural features the code actually depends on: an is_a chain, a
polyhierarchy (microglia under two anchors), a species-qualified term, an obsolete term,
a synonym that collides with another term's primary label, and an organ-qualified name
whose bare form is absent.
"""
import os
import sys

# Test the repository, not whatever happens to be installed. A stale celltype-audit in
# site-packages shadowed src/ here and every local run silently exercised 0.1.0 -- including
# test_stated_numbers, the guard written to catch exactly this kind of drift. Editable
# installs are not available on this machine's pip, so the path is set explicitly.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import pytest

from celltype_audit.ontology import Ontology


def _node(cid, label, syns=()):
    return {"id": "http://purl.obolibrary.org/obo/" + cid.replace(":", "_"), "lbl": label,
            "meta": {"synonyms": [{"val": s} for s in syns]}}


def _edge(sub, obj):
    return {"sub": "http://purl.obolibrary.org/obo/" + sub.replace(":", "_"),
            "pred": "is_a",
            "obj": "http://purl.obolibrary.org/obo/" + obj.replace(":", "_")}


@pytest.fixture(scope="session")
def ontology():
    nodes = [
        _node("CL:0000988", "hematopoietic cell", ["haematopoietic cell"]),
        _node("CL:0000066", "epithelial cell"),
        _node("CL:0002320", "connective tissue cell"),
        _node("CL:0000187", "muscle cell"),
        _node("CL:0002319", "neural cell"),
        _node("CL:0000235", "macrophage"),
        _node("CL:0000394", "plasmatocyte", ["macrophage"]),      # synonym collision
        _node("CL:0000057", "fibroblast"),
        _node("CL:0000186", "myofibroblast cell"),
        _node("CL:0000129", "microglial cell"),                    # polyhierarchy
        _node("CL:0004117", "alpha retinal ganglion cell (Mmus)", ["alpha cell"]),
        _node("CL:0000171", "pancreatic A cell", ["pancreatic alpha cell"]),
        _node("CL:9999999", "obsolete thing", ["ghost cell"]),
        _node("CL:0000775", "neutrophil"),
        _node("CL:0000860", "classical monocyte"),
        _node("CL:0000576", "monocyte"),
        _node("CL:0001054", "CD14-positive monocyte"),
    ]
    edges = [
        _edge("CL:0000235", "CL:0000988"), _edge("CL:0000775", "CL:0000988"),
        _edge("CL:0000576", "CL:0000988"), _edge("CL:0000860", "CL:0000576"),
        _edge("CL:0001054", "CL:0000576"),
        _edge("CL:0000057", "CL:0002320"), _edge("CL:0000186", "CL:0002320"),
        _edge("CL:0000129", "CL:0000988"), _edge("CL:0000129", "CL:0002319"),
        _edge("CL:0000171", "CL:0000066"),
    ]
    return Ontology.from_json({"graphs": [{"nodes": nodes, "edges": edges,
                                           "meta": {"basicPropertyValues":
                                                    [{"pred": "owl:versionInfo",
                                                      "val": "test"}]}}]})
