"""The label column is not always called `cell_type`.

`audit` hardcoded `obs/cell_type` and offered no way to change it, so any file not written
in the CELLxGENE layout could not be audited at all -- and the failure was a bare KeyError
from deep inside h5py, which reads as a corrupt file rather than a naming difference.
"""
import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from celltype_audit.markers import marker_table          # noqa: E402


def _h5ad(path, type_col):
    with h5py.File(path, "w") as f:
        obs, var = f.create_group("obs"), f.create_group("var")
        g = obs.create_group(type_col)
        g.create_dataset("categories", data=np.array([b"T cell", b"B cell"]))
        g.create_dataset("codes", data=np.array([0, 0, 1, 1], dtype=np.int8))
        obs.create_dataset("donor", data=np.zeros(4))
        gv = var.create_group("feature_name")
        gv.create_dataset("categories", data=np.array([b"CD3E", b"MS4A1"]))
        gv.create_dataset("codes", data=np.array([0, 1], dtype=np.int8))
        # a dense 4x2 counts matrix: CD3E in the T cells, MS4A1 in the B cells
        X = f.create_dataset("X", data=np.array([[3., 0.], [2., 0.],
                                                 [0., 4.], [0., 1.]], dtype=np.float32))
        X.attrs["shape"] = np.array([4, 2])
    return str(path)


def test_a_missing_column_names_what_is_there(tmp_path):
    p = _h5ad(tmp_path / "a.h5ad", "annotation")
    with pytest.raises(KeyError) as e:
        marker_table(p)
    msg = str(e.value)
    assert "obs/cell_type not found" in msg
    assert "--type-key" in msg, "the error must say how to fix it"
    assert "annotation" in msg and "donor" in msg, "it must list the columns present"


def test_type_key_selects_the_column(tmp_path):
    """The fix itself: naming the column makes the file readable."""
    p = _h5ad(tmp_path / "b.h5ad", "annotation")
    tbl = marker_table(p, type_key="annotation", min_cells=1)
    assert set(tbl) == {"T cell", "B cell"}


def test_audit_h5ad_accepts_type_key():
    """Threaded all the way through, not just available on the low-level function."""
    import inspect

    from celltype_audit import audit_h5ad
    assert inspect.signature(audit_h5ad).parameters["type_key"].default == "cell_type"
