"""The return path: a gap must close when CL covers it, and must not claim credit.

analysis/proposal_state.py is stage 6 -- it re-runs the check that produced each proposal
against a later Cell Ontology release and closes what the ontology now covers. Two
properties matter enough to pin down, because both are easy to break silently:

  * a closure is DETECTED, so the loop actually turns; and
  * a closure is not ATTRIBUTED, because CL gains terms constantly and a curator may have
    added one independently of anything we submitted.

The module lives in analysis/ rather than the package, so it is loaded by path.
"""
import importlib.util
import json
import os
from datetime import date

import pytest

ANALYSIS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "analysis", "proposal_state.py")


@pytest.fixture
def mod(monkeypatch, tmp_path):
    """Load stage 6 with its state file redirected into a temp directory."""
    if not os.path.exists(ANALYSIS):
        pytest.skip("analysis/proposal_state.py not present")
    spec = importlib.util.spec_from_file_location("proposal_state", ANALYSIS)
    m = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(m)
    except Exception as exc:                      # needs the CL graph, absent in CI
        pytest.skip("stage 6 needs the ontology on disk: %s" % exc)
    monkeypatch.setattr(m, "STATE", str(tmp_path / "state.json"))
    return m


def _seed(m, kind="new-term", organ="Kidney", label="Principle cell"):
    state = {"schema": m.SCHEMA, "first_run": "2026-08-21", "runs": [], "proposals": {}}
    m.sync(state, [{"kind": kind, "organ": organ, "label": label, "n_cells": 11753}],
           "2026-06-08")
    return state


def test_a_new_gap_is_recorded_once(mod):
    state = _seed(mod)
    assert len(state["proposals"]) == 1
    (s,) = state["proposals"].values()
    assert s["status"] == "drafted"
    assert s["first_seen_release"] == "2026-06-08"
    # re-running the audit must not duplicate the same gap
    added, _gone = mod.sync(state, [{"kind": "new-term", "organ": "Kidney",
                                     "label": "Principle cell", "n_cells": 11753}],
                            "2026-06-08")
    assert added == 0 and len(state["proposals"]) == 1


def test_closure_is_detected_when_the_ontology_covers_the_gap(mod, monkeypatch):
    state = _seed(mod)
    monkeypatch.setattr(mod, "gap_is_closed", lambda p: True)
    closed, _undec = mod.recheck(state, "2026-09-01")
    (s,) = state["proposals"].values()
    assert len(closed) == 1
    assert s["status"] == "resolved"
    assert s["closed_in_release"] == "2026-09-01"
    assert s["history"][-1]["status"] == "resolved"


def test_closure_does_not_claim_credit(mod, monkeypatch):
    """The property the whole design rests on: measured, not attributed."""
    state = _seed(mod)
    monkeypatch.setattr(mod, "gap_is_closed", lambda p: True)
    mod.recheck(state, "2026-09-01")
    (s,) = state["proposals"].values()
    assert s["accepted"] is None, "a release diff must never set accepted"
    assert "attribution not implied" in s["history"][-1]["note"]
    m = mod.metrics(state, "2026-09-01")
    assert m["accepted"] == 0
    assert m["acceptance_rate"] is None, "no acceptance rate without an adjudicated issue"


def test_an_open_gap_stays_open(mod, monkeypatch):
    state = _seed(mod)
    monkeypatch.setattr(mod, "gap_is_closed", lambda p: False)
    closed, _ = mod.recheck(state, "2026-09-01")
    (s,) = state["proposals"].values()
    assert closed == [] and s["status"] == "drafted"


def test_undecidable_gaps_are_counted_not_guessed(mod, monkeypatch):
    state = _seed(mod, kind="marker-condition", organ="Lung", label="neutrophil")
    monkeypatch.setattr(mod, "gap_is_closed", lambda p: None)
    closed, undec = mod.recheck(state, "2026-09-01")
    (s,) = state["proposals"].values()
    assert closed == [] and undec == 1
    assert s["status"] == "drafted", "undecidable must not be silently resolved"


def test_a_gap_the_audit_drops_is_withdrawn_not_resolved(mod):
    state = _seed(mod)
    _added, gone = mod.sync(state, [], "2026-06-08")
    (s,) = state["proposals"].values()
    assert gone == 1
    assert s["status"] == "withdrawn", "vanishing from the audit is not the same as closure"


def test_acceptance_rate_needs_an_adjudicated_issue(mod):
    state = _seed(mod)
    (s,) = state["proposals"].values()
    assert mod.metrics(state, "x")["acceptance_rate"] is None
    s["issue"], s["status"] = 4321, "submitted"
    assert mod.metrics(state, "x")["acceptance_rate"] == 0.0
    s["accepted"] = True
    assert mod.metrics(state, "x")["acceptance_rate"] == 100.0
