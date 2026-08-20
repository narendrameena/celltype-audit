"""De-novo proposal path: cluster discovery, tissue rollup, and honest failure."""
import pytest

from celltype_audit.annotate import find_cluster_key, CLUSTER_KEYS
from celltype_audit.reference import ROLLUP


def test_rollup_maps_substructures_to_a_referenced_part():
    """CELLxGENE's expression vocabulary is coarser than UBERON."""
    assert ROLLUP["UBERON:0001225"] == "UBERON:0002113"   # kidney cortex -> kidney
    assert ROLLUP["UBERON:0000017"] == "UBERON:0001264"   # exocrine pancreas -> pancreas
    assert all(v not in ROLLUP for v in ROLLUP.values()), "rollup targets must be terminal"


def test_cluster_key_candidates_cover_the_usual_pipelines():
    for k in ("leiden", "louvain", "seurat_clusters", "cell_type"):
        assert k in CLUSTER_KEYS


def test_missing_cluster_column_fails_loudly(tmp_path):
    """celltype-audit does not cluster; it must say so rather than guess."""
    h5py = pytest.importorskip("h5py")
    p = tmp_path / "x.h5ad"
    with h5py.File(p, "w") as f:
        f.create_group("obs").create_dataset("barcode", data=[b"a", b"b"])
    with pytest.raises(ValueError) as e:
        find_cluster_key(str(p))
    assert "does not cluster" in str(e.value)


def test_explicit_cluster_key_must_exist(tmp_path):
    h5py = pytest.importorskip("h5py")
    p = tmp_path / "y.h5ad"
    with h5py.File(p, "w") as f:
        f.create_group("obs").create_dataset("leiden", data=[b"0", b"1"])
    assert find_cluster_key(str(p), "leiden") == "leiden"
    with pytest.raises(ValueError):
        find_cluster_key(str(p), "not_there")
