"""The subspace scorer must actually run, and must not silently fall back to the mean.

Three separate times while adopting this, the cosine path looked adopted and was not:
the genes it needs were missing from the reference, `score` quietly took the mean branch,
and the run reproduced the previous number to the digit. Identical output is what a silent
fallback looks like, so the only safe guard is data where the two paths disagree -- if the
subspace branch stops executing, the assertion moves, rather than the result staying
plausible.

The fixture reproduces the phenomenon the scorer exists for. T2 is abundant everywhere and
T1 is specific to the cluster's markers. The mean over the cluster's own two markers ranks
the abundant term first, because abundance raises every gene it touches. The cosine over
the shared subspace ranks the specific term first, because the genes the cluster does NOT
express carry weight zero in the query and full weight in T2's norm.
"""
import numpy as np
import pytest

from celltype_audit.reference import Reference

TIS = "UBERON:0000000"
GENES = ["AAA", "BBB", "CCC", "DDD"]


@pytest.fixture
def ref():
    #                AAA   BBB   CCC   DDD
    M = np.array([[10.0, 10.0, 0.0, 0.0],      # T1: specific
                  [12.0, 12.0, 12.0, 12.0]],   # T2: abundant everywhere
                 dtype=np.float32)
    return Reference({TIS: M},
                     {TIS: {g: i for i, g in enumerate(GENES)}},
                     {TIS: {"CL:0000001": 0, "CL:0000002": 1}},
                     {TIS: {"CL:0000001": 5000, "CL:0000002": 5000}},
                     {TIS: {"CL:0000001": "specific", "CL:0000002": "abundant"}})


def markers(*names):
    """Marker records, which is what the subspace path needs; rates weight the query."""
    return [{"gene": g, "pc_in": 0.9} for g in names]


def test_subspace_keeps_only_genes_the_reference_has(ref):
    sub = ref.subspace(TIS, [["AAA", "BBB"], ["CCC", "ZZZ"]])
    assert sub == ["AAA", "BBB", "CCC"], "ZZZ is absent from the reference and must drop"


def test_mean_path_prefers_the_abundant_term(ref):
    """Not a wish -- this is the failure the subspace scorer was built to fix."""
    s = ref.score(TIS, markers("AAA", "BBB", "CCC"))
    assert max(s, key=s.get) == "CL:0000002"


def test_subspace_path_prefers_the_specific_term(ref):
    sub = ref.subspace(TIS, [["AAA", "BBB"], ["CCC", "DDD"]])
    s = ref.score(TIS, markers("AAA", "BBB"), subspace=sub)
    assert max(s, key=s.get) == "CL:0000001", "cosine must beat abundance here"


def test_the_two_paths_disagree(ref):
    """The guard against a silent fallback: if these ever agree, the branch stopped running."""
    sub = ref.subspace(TIS, [["AAA", "BBB"], ["CCC", "DDD"]])
    mean = ref.score(TIS, markers("AAA", "BBB", "CCC"))
    cos = ref.score(TIS, markers("AAA", "BBB"), subspace=sub)
    assert max(mean, key=mean.get) != max(cos, key=cos.get)


def test_bare_gene_names_fall_back_and_say_so(ref):
    """No rates means no query vector, so the subspace path cannot run; the mean must."""
    sub = ref.subspace(TIS, [["AAA", "BBB"], ["CCC", "DDD"]])
    s = ref.score(TIS, ["AAA", "BBB", "CCC"], subspace=sub)
    assert max(s, key=s.get) == "CL:0000002"


def test_too_few_reference_genes_falls_back(ref):
    """Fewer than three usable genes is not a subspace; the mean is the honest answer."""
    s = ref.score(TIS, markers("AAA", "BBB", "CCC"), subspace=["AAA", "BBB"])
    assert max(s, key=s.get) == "CL:0000002"


def test_min_ref_gates_thin_terms(ref):
    """A term below the support floor must not be rankable on either path."""
    ref.counts[TIS]["CL:0000001"] = 10
    sub = ref.subspace(TIS, [["AAA", "BBB"], ["CCC", "DDD"]])
    s = ref.score(TIS, markers("AAA", "BBB"), subspace=sub, min_ref=100)
    assert "CL:0000001" not in s
