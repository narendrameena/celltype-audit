def test_ancestors_are_inclusive_and_transitive(ontology):
    a = ontology.ancestors("CL:0000860")           # classical monocyte
    assert "CL:0000860" in a and "CL:0000576" in a and "CL:0000988" in a


def test_anchor_sets_are_read_from_the_graph(ontology):
    assert ontology.anchors("CL:0000057") == frozenset({"connective"})
    assert ontology.anchors("CL:0000235") == frozenset({"haematopoietic"})


def test_polyhierarchy_is_preserved(ontology):
    """A microglial cell is both haematopoietic and neural and must conflict with neither."""
    assert ontology.anchors("CL:0000129") == frozenset({"haematopoietic", "neural"})
    assert not ontology.conflict("CL:0000129", "CL:0000235")
    assert not ontology.conflict("CL:0000129", "CL:0002319")


def test_conflict_needs_both_sides_non_empty(ontology):
    assert ontology.conflict("CL:0000057", "CL:0000775")       # connective vs haemato
    assert not ontology.conflict("CL:9999999", "CL:0000775")   # no anchors -> asserts nothing


def test_related_is_symmetric_over_is_a(ontology):
    assert ontology.related("CL:0000576", "CL:0000860")
    assert ontology.related("CL:0000860", "CL:0000576")
    assert not ontology.related("CL:0000775", "CL:0000057")
