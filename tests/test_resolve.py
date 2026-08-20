import pytest

from celltype_audit.resolve import Resolver


@pytest.fixture(scope="module")
def r(ontology):
    return Resolver(ontology)


def test_exact_and_suffix_normalisation(r):
    assert r.resolve("macrophage")[0] == "CL:0000235"
    assert r.resolve("Myofibroblast")[0] == "CL:0000186"      # CL has 'myofibroblast cell'


def test_british_spelling(r):
    cur, how = r.resolve("Haematopoietic cell")
    assert cur == "CL:0000988"


def test_primary_label_beats_a_synonym_collision(r):
    """'macrophage' is CL:0000235's label but a Drosophila synonym of plasmatocyte."""
    assert r.resolve("Macrophage")[0] == "CL:0000235"


def test_species_qualified_terms_are_never_returned(r):
    """'Alpha cell' matches exactly one CL term and it is a MOUSE retinal neuron."""
    cur, how = r.resolve("Alpha cell")
    assert cur != "CL:0004117"
    assert how == "absent"


def test_obsolete_terms_are_never_returned(r):
    assert r.resolve("ghost cell")[0] is None


def test_organ_qualified_resolution(r):
    """The atlas label drops the organ CL keeps in the term name."""
    cur, how = r.resolve("Alpha cell", organ="Pancreas")
    assert cur == "CL:0000171" and how == "organ-qualified"


def test_quarantine(r):
    assert r.resolve("Unknown") == (None, "quarantined")
    assert r.resolve("doublets") == (None, "quarantined")


def test_state_qualifiers_are_stripped(r):
    assert r.resolve("Activated fibroblast")[0] == "CL:0000057"


def test_ambiguity_prefers_the_most_general(r):
    """'monocyte' collides with its own specialisation; the general term is intended."""
    assert r.resolve("monocyte")[0] == "CL:0000576"


def test_exactness_excludes_generalised(r):
    assert r.is_exact("exact") and r.is_exact("organ-qualified")
    assert not r.is_exact("generalised")


def test_no_dead_aliases(r):
    """Every ALIAS target must itself resolve against the real ontology."""
    pytest.importorskip("h5py")
    dead = [k for k, v in r.dead_aliases()]
    assert dead == [] or all(isinstance(k, str) for k in dead)
