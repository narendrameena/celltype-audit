"""Properties of the independent-annotator comparison that are easy to break silently.

agreement.py answers the single-annotator objection by comparing the gold against the
HRA / ASCT+B crosswalks. Three things have to stay true or the number it reports stops
meaning what the paper says it means.
"""
import importlib.util
import json
import os

import pytest

ANALYSIS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "analysis", "agreement.py")


@pytest.fixture
def mod():
    if not os.path.exists(ANALYSIS):
        pytest.skip("analysis/agreement.py not present")
    spec = importlib.util.spec_from_file_location("agreement", ANALYSIS)
    m = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(m)
    except Exception as exc:                       # needs the CL graph, absent in CI
        pytest.skip("agreement needs the ontology on disk: %s" % exc)
    return m


def test_the_same_label_in_two_organs_is_two_annotations(mod, tmp_path, monkeypatch):
    """The bug this guards: keying on the label alone drops one of them, and which one
    survives depends on directory order. `Conventional dendritic cell` is curated in both
    skin and spleen and is an error in both."""
    for organ in ("skin", "spleen"):
        (tmp_path / ("%s_gold.json" % organ)).write_text(json.dumps(
            {"Conventional dendritic cell": {"curie": "CL:0000784"}}))
    monkeypatch.setattr(mod, "GOLD", str(tmp_path))
    idx = mod.gold_index()
    assert len(idx) == 2, "one annotation per (organ, label), not per label"
    assert {k[0] for k in idx} == {"skin", "spleen"}


def test_underscore_keys_are_not_annotations(mod, tmp_path, monkeypatch):
    """_abstentions / _disagreements are curator notes, not assignments."""
    (tmp_path / "lung_gold.json").write_text(json.dumps(
        {"_abstentions": {"T cell": "cell-cycle programme only"},
         "Tuft cell": {"curie": "CL:0002204"}}))
    monkeypatch.setattr(mod, "GOLD", str(tmp_path))
    assert len(mod.gold_index()) == 1


def test_a_granularity_step_is_not_a_disagreement(mod):
    """Two curators choosing different depths on one is_a path agree about identity."""
    from cl_lineage import ancestors
    narrow = "CL:0000895"                          # naive thymus-derived CD4+ ab T cell
    try:                                           # the CL graph is a data file, not in the repo
        anc = [a for a in ancestors(narrow) if a.startswith("CL:")]
    except FileNotFoundError:
        pytest.skip("lineage graph unavailable")
    if not anc:
        pytest.skip("lineage graph unavailable")
    assert mod.classify(narrow, {anc[0]}) == "granularity"
    assert mod.classify(narrow, {narrow}) == "exact"


def test_kappa_chance_is_negligible_over_an_open_vocabulary(mod):
    """The reason the paper reports percent agreement and treats kappa as decoration."""
    po, pe, k = mod.kappa({"exact": 141}, 169)
    assert pe < 0.01
    assert abs(k - po) < 0.05
