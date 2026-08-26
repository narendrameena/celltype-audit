"""V12 must not claim CL lacks a term for a population CL demonstrably names.

V12 asks whether the improved scorer now returns a term that already names the proposal's
population. Its first form compared bags of word stems against the single top candidate,
and that was wrong in two ways at once: it read only rank 1, so `VIP GABAergic interneuron`
sitting at rank 2 went unseen; and stems alone cannot equate `PV` with `pvalb`,
`inhibitory` with `GABAergic`, or `Mesothelium` with `mesothelial cell`. Seven of nineteen
proposals asserted a gap in CL over a term CL already had.

The direction of the error is what makes it worth a test. A false "gap stands" sends a
curator to write a term that exists; a false "already named" only withdraws a proposal the
audit could have made. These pin the pairs that were actually wrong, and the negatives that
must keep failing so the equivalence table cannot be widened into uselessness.
"""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "analysis"))

verify = pytest.importorskip("verify_proposals",
                             reason="analysis scripts need the study's CL release")
toks = verify._toks


def names_same(label, candidate):
    return bool(toks(label)) and toks(label) <= toks(candidate)


@pytest.mark.parametrize("label,candidate", [
    ("PV inhibitory neuron", "pvalb GABAergic interneuron"),
    ("VIP inhibitory neuron", "VIP GABAergic interneuron"),
    ("Mesothelium", "mesothelial cell"),
    ("Serous gland cell", "serous secreting cell"),
    ("OFF-BC", "OFF-bipolar cell"),
    ("ON-BC", "ON-bipolar cell"),
    ("Bipolar cell", "rod bipolar cell"),
    ("Umbrella cell", "umbrella cell of urothelium"),
])
def test_the_same_population_written_two_ways(label, candidate):
    assert names_same(label, candidate), (
        "%r and %r name one population; V12 would assert a gap CL does not have"
        % (label, candidate))


@pytest.mark.parametrize("label,candidate", [
    ("Capsular cell", "mesenchymal stem cell"),
    ("Roof plate (RP) cell", "astrocyte"),
    ("Floor plate (FP) cell", "ependymal cell"),
    ("Midplate (MP) cell", "astrocyte"),
    ("Spermatogenic cell", "primordial germ cell"),
    ("Serous gland cell", "endothelial cell of vascular tree"),
])
def test_different_populations_stay_different(label, candidate):
    assert not names_same(label, candidate), (
        "%r and %r are different populations; matching them would withdraw a real proposal"
        % (label, candidate))


def test_the_equivalence_table_stays_small_and_justified():
    """Every entry is a curated claim about naming, so the table must not drift into a
    general synonym list -- that is what the ontology is for."""
    assert len(verify.EQUIV) <= 20, "EQUIV is a short curated list, not a thesaurus"
    for k, v in verify.EQUIV.items():
        assert k != v or k in ("sst",), "an entry that maps a token to itself does nothing"
