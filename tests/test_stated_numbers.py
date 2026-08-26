"""The numbers the package states about itself must match the results that produced them.

Twice now a figure or a docstring has drifted from its own source and nothing caught it,
because a printed string is not exercised by any test that checks a value. ACCURACY_NOTE is
the worst case of the two: it goes to every user at run time, and 0.1.2 shipped it saying
"~55% ... 200 curated cell types, 7 organs ... AUC 0.563" long after the gold reached ten
organs and the AUC reached 0.72 -- which was not merely a stale number but a different
claim, since 0.563 says the confidence is worthless and 0.72 says it ranks but cannot gate.

This parses the numbers back out of the shipped strings and compares them with the results
files. It skips where those files are absent, which is the normal case in CI: the results
are study outputs, not package data. That means the test only bites where it can, on a
machine that has the pipeline -- which is precisely where the strings get edited.
"""
import json
import os
import re

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.abspath(os.path.join(HERE, "..", "..", "cellscribe_tool",
                                       "benchmark", "results"))


def _load(name):
    p = os.path.join(RESULTS, name)
    if not os.path.exists(p):
        pytest.skip("study results not present (%s); nothing to check against" % name)
    return json.load(open(p))


@pytest.fixture(scope="module")
def gold():
    """Pooled top-1 over every curated organ, the way fig2 panel b pools it."""
    d = _load("organ_scores.json")
    n = sum(s["n>=500"]["n"] for s in d)
    top1 = sum(s["n>=500"]["n"] * s["n>=500"]["top1"] / 100.0 for s in d) / n * 100
    per = [s["n>=500"]["top1"] for s in d]
    return {"n": n, "top1": top1, "organs": len(d), "lo": min(per), "hi": max(per)}


@pytest.fixture(scope="module")
def note():
    from celltype_audit.annotate import ACCURACY_NOTE
    return ACCURACY_NOTE


def test_top1_percentage(note, gold):
    m = re.search(r"top-1 is right ~(\d+)%", note)
    assert m, "ACCURACY_NOTE no longer states a top-1 percentage"
    # exact round-match, not a tolerance: the string is written by rounding the result,
    # so anything else is drift. A +/-1 window let "35-86%" survive a true 35-86.7%.
    assert int(m.group(1)) == round(gold["top1"]), (
        "note says ~%s%%, organ_scores.json pools to %.1f%%" % (m.group(1), gold["top1"]))


def test_gold_size_and_organ_count(note, gold):
    m = re.search(r"\((\d+) curated cell types, (\d+) organs", note)
    assert m, "ACCURACY_NOTE no longer states the gold's size"
    assert int(m.group(1)) == gold["n"], (
        "note says %s cell types, organ_scores.json has %d" % (m.group(1), gold["n"]))
    assert int(m.group(2)) == gold["organs"], (
        "note says %s organs, organ_scores.json has %d" % (m.group(2), gold["organs"]))


def test_per_tissue_range(note, gold):
    m = re.search(r"(\d+)-(\d+)% by tissue", note)
    assert m, "ACCURACY_NOTE no longer states a per-tissue range"
    lo, hi = int(m.group(1)), int(m.group(2))
    assert (lo, hi) == (round(gold["lo"]), round(gold["hi"])), (
        "note says %d-%d%%, organs run %.1f-%.1f%%" % (lo, hi, gold["lo"], gold["hi"]))


def test_auc_matches_the_calibration_experiment(note):
    """The one that changed meaning, not just value."""
    auc = _load("calibration.json")["auc_gold_organs"]
    m = re.search(r"AUC (\d\.\d+)", note)
    assert m, "ACCURACY_NOTE no longer states an AUC"
    assert abs(float(m.group(1)) - auc) < 0.01, (
        "note says AUC %s, calibration.json has %.3f" % (m.group(1), auc))


def test_the_module_docstring_agrees_with_the_note(gold):
    """The docstring is the long form of the same claim and drifted independently."""
    import celltype_audit.annotate as a
    doc = a.__doc__
    m = re.search(r"against (\d+) hand-curated cell types in (\w+) organs", doc)
    assert m, "annotate.__doc__ no longer states the gold's size"
    assert int(m.group(1)) == gold["n"]
    words = {"four": 4, "five": 5, "six": 6, "seven": 7, "eight": 8, "nine": 9, "ten": 10}
    assert words.get(m.group(2)) == gold["organs"], (
        "docstring says %r organs, organ_scores.json has %d" % (m.group(2), gold["organs"]))


# --------------------------------------------------------------------------- survey
# The README's survey figures went stale exactly once and it took a reference refetch to
# notice: survey_ref had been fetched for each cluster's own five markers, the subspace
# scorer silently fell back to the mean, and the published 3.1% survived a rerun that
# reproduced it to the digit. These read the numbers back out of the README and compare
# them with cxg_survey.json, so the prose cannot outlive the result again.
@pytest.fixture(scope="module")
def survey():
    rows = _load("cxg_survey.json")
    live = [r for r in rows if not r["thin_support"]]
    cross = [r for r in live if r["cross_lineage"]]
    return {"datasets": len({r["dataset"] for r in rows}), "surviving": len(live),
            "cross": len(cross), "pct": 100.0 * len(cross) / max(len(live), 1)}


@pytest.fixture(scope="module")
def readme():
    p = os.path.join(HERE, "..", "README.md")
    if not os.path.exists(p):
        pytest.skip("README not present")
    return open(p).read()


def test_readme_survey_counts(readme, survey):
    m = re.search(r"\*\*(\d+) of ([\d,]+) cell types \(([\d.]+)%", readme)
    assert m, "README no longer states the survey's cross-lineage count"
    assert int(m.group(1)) == survey["cross"], (
        "README says %s cross-lineage, cxg_survey.json has %d" % (m.group(1), survey["cross"]))
    assert int(m.group(2).replace(",", "")) == survey["surviving"], (
        "README says %s surviving, cxg_survey.json has %d" % (m.group(2), survey["surviving"]))
    assert abs(float(m.group(3)) - survey["pct"]) < 0.05, (
        "README says %s%%, cxg_survey.json gives %.1f%%" % (m.group(3), survey["pct"]))


def test_readme_dataset_count(readme, survey):
    m = re.search(r"(\d+) had a scoreable cell type", readme)
    assert m, "README no longer states how many datasets were scoreable"
    assert int(m.group(1)) == survey["datasets"], (
        "README says %s datasets, cxg_survey.json has %d" % (m.group(1), survey["datasets"]))
