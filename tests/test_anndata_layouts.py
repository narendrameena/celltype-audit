"""An .h5ad may store a categorical three ways, and two of them were read wrongly.

AnnData before 0.8 wrote a categorical column as an INTEGER array of codes with the
categories in a sibling `__categories/<name>` group. HuBMAP still publishes that layout.
The reader assumed the current one -- a group holding `categories` and `codes` -- and on
the legacy form fell through to treating the codes as values, so `var/hugo_symbol` would
have yielded the integers 29035, 28677, ... where gene symbols were expected. It never got
that far, because a dead guard called `.get()` on an h5py Dataset and raised first, but a
crash was the lucky outcome: the silent path returns numbers that look like data.

A code of -1 means "no category" and is how a missing HUGO symbol is stored. Left alone it
is a negative index and selects the LAST category, quietly mislabelling those genes.

These build each layout in a temporary file and check the decoded values, because the bug
was in decoding and only data can show it.
"""
import os
import tempfile

import h5py
import numpy as np
import pytest

from celltype_audit.markers import _cat

GENES = ["CD3E", "MS4A1", "LYZ", "PTPRC"]


@pytest.fixture
def h5(tmp_path):
    p = str(tmp_path / "layouts.h5ad")
    with h5py.File(p, "w") as f:
        # current: a group with categories + codes
        g = f.create_group("var/modern")
        g.create_dataset("categories", data=np.array(GENES, dtype="S8"))
        g.create_dataset("codes", data=np.array([0, 1, 2, 3], dtype="i1"))
        # legacy: an integer dataset of codes, categories in a sibling __categories group
        f.create_dataset("var/legacy", data=np.array([3, 2, 1, 0], dtype="i2"))
        f.create_dataset("var/legacy_missing", data=np.array([0, -1, 2, -1], dtype="i2"))
        lc = f.create_group("var/__categories")
        lc.create_dataset("legacy", data=np.array(GENES, dtype="S8"))
        lc.create_dataset("legacy_missing", data=np.array(GENES, dtype="S8"))
        # plain: an array of strings
        f.create_dataset("var/plain", data=np.array(GENES, dtype="S8"))
    return p


def decode(path, key):
    with h5py.File(path, "r") as f:
        cats, codes = _cat(f, key)
        return list(np.asarray(cats)[codes])


def test_current_layout(h5):
    assert decode(h5, "var/modern") == GENES


def test_legacy_layout_is_decoded_not_read_as_values(h5):
    """The failure this test exists for: integers where gene symbols belong."""
    got = decode(h5, "var/legacy")
    assert got == ["PTPRC", "LYZ", "MS4A1", "CD3E"], got
    assert not any(str(x).lstrip("-").isdigit() for x in got), (
        "codes were read as values; the categories were never applied")


def test_a_missing_category_does_not_wrap_to_the_last_one(h5):
    """-1 means no category. As an index it silently selects PTPRC."""
    got = decode(h5, "var/legacy_missing")
    assert got[1] == "" and got[3] == "", got
    assert "PTPRC" not in got, "a -1 code wrapped round to the last category"


def test_plain_string_column(h5):
    assert sorted(decode(h5, "var/plain")) == sorted(GENES)
