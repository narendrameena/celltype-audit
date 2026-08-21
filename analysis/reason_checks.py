#!/usr/bin/env python3
"""V4-V8: the checks that need a reasoner, run atomically.

Each check is an independent function of (axiom, ontology) -> (verdict, evidence). None
depends on another having run, none shares mutable state, and each records the evidence it
decided on, so a verdict can be re-derived rather than trusted. Given the same CL release
and the same axiom the verdicts are identical, which is what lets them sit beside the
deterministic checks rather than under a caveat.

  V4  parse          the axiom is well-formed OWL and merges into CL
  V5  consistency    no class becomes unsatisfiable
  V6  entailment     the inferred class hierarchy gains nothing outside the axiom's own
                     signature -- an axiom that silently re-parents something else is the
                     failure a drafter cannot see
  V7  conservativity nothing new is entailed BETWEEN classes that already existed
  V8  non-triviality something new is entailed at all; an axiom that changes no entailment
                     is noise, and LLM drafts produce these constantly

Only proposals that assert an axiom are in scope. A request for a new term or a synonym
asserts nothing for a reasoner to check, so those record `n/a` -- which is not a pass, and
is not a failure either.

Usage: python benchmark/reason_checks.py
Env:   CL_OWL, ROBOT_JAR, JAVA_XMX (default 8g)
"""
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
DOCS = os.path.abspath(os.path.join(HERE, "..", "..", "celltype-audit", "docs"))
CL_OWL = os.environ.get("CL_OWL", os.path.join(ROOT, ".tools", "cl-base.owl"))
ROBOT = os.environ.get("ROBOT_JAR", os.path.join(ROOT, ".tools", "robot.jar"))
XMX = os.environ.get("JAVA_XMX", "8g")

IRI = lambda c: "http://purl.obolibrary.org/obo/%s" % c.replace(":", "_")
# ROBOT writes functional syntax with PREFIXED names (obo:CL_0002131), not full IRIs.
# Matching only <...> parsed nothing, so the entailment diff compared two empty sets and
# V6/V7 passed vacuously -- a green badge with no evidence behind it, which is worse than
# a red one. Accept both forms and normalise.
_TERM = r"(?:<([^>]+)>|([A-Za-z][\w.-]*:[\w.-]+))"
SUBCLASS = re.compile(r"SubClassOf\(\s*%s\s+%s\s*\)" % (_TERM, _TERM))


def _norm_term(iri, curie):
    """obo:CL_0002131 and <http://purl.obolibrary.org/obo/CL_0002131> are the same class."""
    if iri:
        return iri.rsplit("/", 1)[-1].replace("_", ":")
    if curie.startswith("obo:"):
        return curie[4:].replace("_", ":")
    return curie


def robot(args, timeout=1800):
    p = subprocess.run(["java", "-Xmx" + XMX, "-jar", ROBOT] + args,
                       capture_output=True, text=True, timeout=timeout)
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def axiom_for(p):
    """The OWL a proposal asserts, or None when it asserts nothing checkable.

    A marker condition is written as a NECESSARY condition (the type implies the markers),
    not a sufficient one. CL has no marker vocabulary to point at, so the genes are carried
    as an annotation and the logical content is the subclass assertion; claiming a
    sufficient condition we cannot express would be the sort of overreach this stack exists
    to stop.
    """
    if p["kind"] == "marker-condition" and p.get("candidate"):
        c = p["candidate"]["curie"]
        return {"kind": "annotation-only", "subject": c, "ofn": None,
                "note": ("CL carries no marker property, so this proposal cannot be stated "
                         "as an axiom against the current release; it is a request for one "
                         "rather than an axiom to check.")}
    if p["kind"] == "missing-axiom":
        # the pair is named in the proposal, not guessed from a label pattern: matching on
        # "ventricular cardiomyocyte" found neither term, since CL calls them "regular
        # ventricular cardiac myocyte" and "ventricular cardiac muscle cell"
        a = (p.get("candidate") or {}).get("curie")
        b = ((p.get("alternatives") or [{}])[0]).get("curie")
        if not a or not b:
            return None
        return {"kind": "subclass", "subject": a, "object": b,
                "ofn": "SubClassOf(<%s> <%s>)" % (IRI(a), IRI(b)),
                "note": "asserts an is_a path between the two sibling terms"}
    return None


_LAB = {}


def _cl_labels():
    if _LAB:
        return _LAB
    from cl_lineage import load
    _LAB.update(load()["label"])
    return _LAB


def inferred_pairs(extra_ofn=None, cache={}):
    """SubClassOf pairs ELK infers, with the axiom merged in if given. Deterministic."""
    key = extra_ofn or "__baseline__"
    if key in cache:
        return cache[key]
    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, "reasoned.ofn")
        args = ["merge", "--input", CL_OWL]
        if extra_ofn:
            args += ["--input", extra_ofn]
        args += ["reason", "--reasoner", "ELK", "--axiom-generators", "SubClass",
                 "--output", out]
        code, log = robot(args)
        if code != 0 or not os.path.exists(out):
            cache[key] = (None, log[-800:])
            return cache[key]
        pairs = set()
        with open(out) as fh:
            for line in fh:
                m = SUBCLASS.search(line)
                if m:
                    a = _norm_term(m.group(1), m.group(2))
                    b = _norm_term(m.group(3), m.group(4))
                    pairs.add((a, b))
        cache[key] = (pairs, "")
    return cache[key]


# ------------------------------------------------------------------ atomic checks
def v4_parse(ax, td):
    ofn = os.path.join(td, "ax.ofn")
    open(ofn, "w").write("Prefix(:=<http://x/>)\nOntology(<http://x/ax>\n%s\n)\n" % ax["ofn"])
    code, log = robot(["merge", "--input", CL_OWL, "--input", ofn,
                       "--output", os.path.join(td, "merged.owl")])
    return ("pass" if code == 0 else "fail",
            "merged cleanly" if code == 0 else log[-300:], ofn)


def v5_consistency(ofn, td):
    unsat = os.path.join(td, "unsat.owl")
    code, log = robot(["merge", "--input", CL_OWL, "--input", ofn, "reason",
                       "--reasoner", "ELK", "--dump-unsatisfiable", unsat,
                       "--output", os.path.join(td, "r.owl")])
    if code != 0 and "unsatisfiable" in log.lower():
        return "fail", "a class became unsatisfiable"
    if os.path.exists(unsat) and os.path.getsize(unsat) > 0:
        n = sum(1 for line in open(unsat) if "Class(" in line)
        if n:
            return "fail", "%d unsatisfiable class(es)" % n
    return ("pass", "no class becomes unsatisfiable") if code == 0 else ("fail", log[-300:])


def v6_v7_v8(ax, ofn):
    """One reasoner pair, three independent verdicts read off it."""
    base, err1 = inferred_pairs()
    withx, err2 = inferred_pairs(ofn)
    if base is None or withx is None:
        e = (err1 or err2)[-200:]
        return {"V6": ("not-run", e), "V7": ("not-run", e), "V8": ("not-run", e)}
    if not base:
        # nothing parsed: any diff would be vacuous, so refuse to report a verdict
        e = "the reasoned ontology yielded no SubClassOf axioms to compare"
        return {"V6": ("not-run", e), "V7": ("not-run", e), "V8": ("not-run", e)}
    new = withx - base
    sig = {ax["subject"]} | ({ax["object"]} if ax.get("object") else set())
    outside = {(a, b) for a, b in new if a not in sig and b not in sig}
    existing = {(a, b) for a, b in new if a not in sig or b not in sig}
    return {
        "V6": ("pass" if not outside else "fail",
               "%d new entailment(s), none outside the axiom's signature" % len(new)
               if not outside else "%d entailment(s) outside the signature" % len(outside)),
        "V7": ("pass" if not existing else "fail",
               "nothing new entailed between pre-existing classes" if not existing
               else "%d new entailment(s) touch pre-existing classes" % len(existing)),
        "V8": ("pass" if new else "fail",
               "%d new entailment(s)" % len(new) if new else "entails nothing new"),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="recompute and fail on disagreement instead of writing; "
                         "this is what CI runs, so a CL release that changes a verdict "
                         "is caught rather than silently absorbed")
    args = ap.parse_args()
    for f, what in ((CL_OWL, "CL OWL"), (ROBOT, "robot.jar")):
        if not os.path.exists(f):
            sys.exit("missing %s: %s" % (what, f))
    path = os.path.join(DOCS, "proposals.json")
    doc = json.load(open(path))
    props = doc["proposals"]
    # keyed by the gap, not the label: two organs share "Multipotent hematopoietic
    # progenitor", and keying on the label alone had one proposal overwrite the other's
    # snapshot, producing a phantom drift report
    key = lambda q: "%s|%s|%s" % (q["kind"], q["organ"], q["label"])
    before = {key(p): dict(p["checks"]) for p in props}

    scoped = [p for p in props if axiom_for(p)]
    print("V4-V8 apply to %d of %d proposals; the rest assert nothing for a reasoner\n"
          % (len([p for p in scoped if axiom_for(p)["ofn"]]), len(props)))

    results = {}
    for p in props:
        ax = axiom_for(p)
        if not ax or not ax["ofn"]:
            for v in ("V4", "V5", "V6", "V7", "V8"):
                p["checks"][v] = "n/a"
            if ax:
                p["reason_note"] = ax["note"]
            continue
        print("  %s / %s" % (p["kind"], p["label"]))
        with tempfile.TemporaryDirectory() as td:
            v4, why4, ofn = v4_parse(ax, td)
            p["checks"]["V4"] = v4
            ev = {"V4": why4}
            if v4 == "pass":
                v5, why5 = v5_consistency(ofn, td)
                p["checks"]["V5"] = v5
                ev["V5"] = why5
                rest = v6_v7_v8(ax, ofn)
                for k, (verdict, why) in rest.items():
                    p["checks"][k] = verdict
                    ev[k] = why
            else:
                for k in ("V5", "V6", "V7", "V8"):
                    p["checks"][k] = "not-run"
        p["reason_note"] = ax["note"]
        p["reason_evidence"] = ev
        results[p["label"]] = ev
        for k in ("V4", "V5", "V6", "V7", "V8"):
            print("     %-3s %-8s %s" % (k, p["checks"][k], ev.get(k, "")[:70]))

    # a check that could not run is never reported as a pass; if the reasoner produced
    # nothing to compare, say so rather than banking a vacuous green
    vacuous = [p["label"] for p in props
               if axiom_for(p) and axiom_for(p)["ofn"]
               and any(p["checks"].get(v) == "not-run" for v in ("V6", "V7", "V8"))]
    if vacuous:
        print("\n  WARNING: entailment checks could not run for: %s" % ", ".join(vacuous))

    if args.check:
        drift = {}
        for p in props:
            was = before[key(p)]
            now = {k: p["checks"][k] for k in was}
            if was != now:
                drift[key(p)] = (was, now)
        if drift:
            print("\nVERDICTS CHANGED against the committed file:")
            for k_, (was, now) in drift.items():
                for c in was:
                    if was[c] != now[c]:
                        print("   %-38s %-4s %s -> %s" % (k_[:38], c, was[c], now[c]))
            print("\nThe ontology moved under a recorded verdict. Re-run without --check,"
                  "\nreview the change, and commit it.")
            sys.exit(1)
        print("\nno verdict changed against the committed file")
        return sys.exit(1) if vacuous else None

    doc["summary"]["checks_run"] = sorted(set(doc["summary"]["checks_run"])
                                          | {"V4", "V5", "V6", "V7", "V8"})
    doc["summary"]["checks_not_run"] = []
    doc["summary"]["reasoner"] = {"engine": "ELK via ROBOT", "ontology": os.path.basename(CL_OWL)}
    json.dump(doc, open(path, "w"), indent=1)
    json.dump(results, open(os.path.join(RES, "reason_checks.json"), "w"), indent=1)
    print("\nwrote %s" % path)


if __name__ == "__main__":
    main()
