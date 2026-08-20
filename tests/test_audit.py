"""Audit logic on a synthetic reference -- no h5ad, no network."""
import numpy as np
import pytest

from celltype_audit.reference import Reference


@pytest.fixture
def ref():
    tis = "UBERON:0000001"
    genes = ["FCN1", "CD14", "S100A12", "FCGR3B", "ELANE", "COL1A1"]
    terms = ["CL:0000775", "CL:0000860", "CL:0000057"]       # neutrophil, mono, fibroblast
    gix = {g: i for i, g in enumerate(genes)}
    tix = {t: i for i, t in enumerate(terms)}
    M = np.zeros((3, 6), dtype=np.float32)
    M[0, [3, 4]] = 3.0                                        # neutrophil: FCGR3B, ELANE
    M[1, [0, 1, 2]] = 3.0                                     # monocyte: FCN1, CD14, S100A12
    M[2, [5]] = 3.0                                           # fibroblast: COL1A1
    counts = {tis: {"CL:0000775": 20000, "CL:0000860": 150000, "CL:0000057": 90000}}
    labels = {tis: {t: t for t in terms}}
    return Reference({tis: M}, {tis: gix}, {tis: tix}, counts, labels), tis


def test_score_picks_the_term_the_markers_support(ref):
    r, tis = ref
    sc = r.score(tis, ["FCN1", "CD14", "S100A12"])
    assert max(sc, key=sc.get) == "CL:0000860"


def test_a_mislabelled_cluster_scores_against_its_own_label(ref):
    """The reference study's flagship: markers say monocyte, the label says neutrophil."""
    r, tis = ref
    sc = r.score(tis, ["FCN1", "CD14", "S100A12"])
    assert sc["CL:0000860"] > sc["CL:0000775"]


def test_min_ref_filters_thin_terms(ref):
    r, tis = ref
    assert "CL:0000775" in r.score(tis, ["FCGR3B", "ELANE", "CD14"], min_ref=100)
    assert "CL:0000775" not in r.score(tis, ["FCGR3B", "ELANE", "CD14"], min_ref=50000)


def test_score_needs_enough_genes(ref):
    r, tis = ref
    assert r.score(tis, ["CD14"]) == {}


def test_roundtrip(tmp_path, ref):
    r, tis = ref
    stem = str(tmp_path / "ref")
    r.save(stem)
    back = Reference.load(stem)
    assert back.tissues() == [tis]
    assert back.support(tis, "CL:0000860") == 150000
