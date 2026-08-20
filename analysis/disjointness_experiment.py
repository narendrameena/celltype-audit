#!/usr/bin/env python3
"""Does injecting disjointness axioms give a reasoner power to refute wrong genus proposals?

Answer: no. Four independent measurements, in order of increasing cost:

  D1  Which major CL lineages can be asserted pairwise disjoint at all?
      (criterion: zero overlap in the is_a descendant closure -- asserting disjointness
       for an overlapping pair would contradict CL's own content.)
  D2  UPPER BOUND on refutation power: for each observed genus error, could ANY assertable
      disjointness axiom separate the proposed genus from the gold genus?
  D3  EMPIRICAL: inject every assertable disjointness axiom into CL, add the candidate
      definitions, run ELK, count what actually becomes unsatisfiable.
      (Doubles as a sanity check: no *existing* CL class may become unsatisfiable.)
  D4  Could *sufficient* conditions (marker -> type) rescue it? Measures how many errors
      even have a differentium to hang such an axiom on, and how many are sibling-level.

Prereqs: benchmark/results/genus_experiment.csv (run genus_experiment.py first),
         cl-full.json, .tools/cl-base.owl, .tools/robot.jar, Java.

Usage:   CL_JSON=cl-full.json python benchmark/disjointness_experiment.py
"""
import csv
import itertools
import json
import os
import subprocess
import sys
import tempfile
from collections import defaultdict, deque

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CL_JSON = os.environ.get("CL_JSON", os.path.join(ROOT, "cl-full.json"))
CL_OWL = os.environ.get("CL_OWL", os.path.join(ROOT, ".tools", "cl-base.owl"))
ROBOT = os.environ.get("ROBOT_JAR", os.path.join(ROOT, ".tools", "robot.jar"))
RES = os.path.join(HERE, "results")
GENUS_CSV = os.path.join(RES, "genus_experiment.csv")
os.makedirs(RES, exist_ok=True)

SHORT = lambda i: i.rsplit("/", 1)[-1].replace("_", ":")
IRI = lambda c: "<http://purl.obolibrary.org/obo/%s>" % c.replace(":", "_")

print("loading %s ..." % os.path.basename(CL_JSON), flush=True)
G = json.load(open(CL_JSON))["graphs"][0]
LAB = {SHORT(n["id"]): n["lbl"] for n in G["nodes"] if n.get("lbl")}
CHILD, PARENT = defaultdict(set), defaultdict(set)
for e in G["edges"]:
    if e.get("pred") == "is_a":
        CHILD[SHORT(e["obj"])].add(SHORT(e["sub"]))
        PARENT[SHORT(e["sub"])].add(SHORT(e["obj"]))
LD = defaultdict(lambda: {"genus": set(), "diff": set()})
for a in G.get("logicalDefinitionAxioms", []):
    cur = SHORT(a["definedClassId"])
    for gid in a.get("genusIds", []):
        LD[cur]["genus"].add(SHORT(gid))
    for r in a.get("restrictions") or []:
        if r and r.get("fillerId"):
            LD[cur]["diff"].add((SHORT(r["propertyId"]), SHORT(r["fillerId"])))


def closure(c, adj):
    s, q = set(), deque([c])
    while q:
        for k in adj.get(q.popleft(), ()):
            if k not in s:
                s.add(k)
                q.append(k)
    return s


# ---------------------------------------------------------------- D1
print("\n=== D1  which lineages can be asserted disjoint at all? ===")
kids = sorted(CHILD["CL:0000000"])
POOL = sorted({k for k in kids if len(closure(k, CHILD)) >= 20} |
              {c for c in ["CL:0000540", "CL:0000066", "CL:0000738", "CL:0000057",
                           "CL:0000115", "CL:0000187", "CL:0000125", "CL:0000232",
                           "CL:0000039", "CL:0000034", "CL:0000988", "CL:0002320"] if c in LAB})
DESC = {c: closure(c, CHILD) | {c} for c in POOL}
SAFE, UNSAFE = [], []
for a, b in itertools.combinations(POOL, 2):
    (SAFE if not (DESC[a] & DESC[b]) else UNSAFE).append((a, b))
print("pool=%d lineages | assertable-disjoint pairs=%d | blocked by overlap=%d"
      % (len(POOL), len(SAFE), len(UNSAFE)))

MEMBER = defaultdict(set)
for c in POOL:
    for x in DESC[c]:
        MEMBER[x].add(c)
SAFESET = {tuple(sorted(p)) for p in SAFE}

# ---------------------------------------------------------------- D2
print("\n=== D2  upper bound on refutation power ===")
rows = list(csv.DictReader(open(GENUS_CSV)))
d2 = {}
for arm in ("A0", "A1", "A2"):
    errs = [r for r in rows if r.get(arm) == "miss" and r.get(arm + "_curie")]
    catch = 0
    for r in errs:
        lp = MEMBER.get(r[arm + "_curie"], set())
        lg = set().union(*[MEMBER.get(x, set()) for x in r["gold"].split(";")]) or set()
        if any(tuple(sorted((a, b))) in SAFESET for a in lp for b in lg if a != b):
            catch += 1
    d2[arm] = {"errors": len(errs), "catchable": catch,
               "rate": round(catch / float(len(errs)), 4) if errs else None}
    print("  %s: %3d errors -> at most %d catchable (%.0f%%)"
          % (arm, len(errs), catch, 100 * catch / max(1, len(errs))))

need = defaultdict(int)
for r in rows:
    if r.get("A1") == "miss" and r.get("A1_curie"):
        lg = set().union(*[MEMBER.get(x, set()) for x in r["gold"].split(";")]) or set()
        for a in MEMBER.get(r["A1_curie"], set()):
            for b in lg:
                if a != b:
                    need[tuple(sorted((a, b)))] += 1
print("  most-needed discriminations:")
for (a, b), n in sorted(need.items(), key=lambda x: -x[1])[:8]:
    print("    %-26s vs %-26s n=%3d  %s" % (LAB.get(a, a)[:26], LAB.get(b, b)[:26], n,
                                            "ASSERTABLE" if tuple(sorted((a, b))) in SAFESET else "BLOCKED (overlaps)"))

# ---------------------------------------------------------------- D3
print("\n=== D3  empirical: inject every assertable disjointness axiom, run ELK ===")
catch_rows, ctrl_rows = [], []
for r in rows:
    p = r.get("A1_curie")
    if not p or not LD.get(r["curie"], {}).get("diff"):
        continue
    lp = MEMBER.get(p, set())
    lg = set().union(*[MEMBER.get(x, set()) for x in r["gold"].split(";")]) or set()
    hit = any(tuple(sorted((a, b))) in SAFESET for a in lp for b in lg if a != b)
    if r["A1"] == "miss" and hit:
        catch_rows.append(r)
    elif r["A1"] in ("exact", "ancestor") and len(ctrl_rows) < 25:
        ctrl_rows.append(r)
cand = catch_rows + ctrl_rows
print("  candidates: %d catchable-in-principle errors + %d correct controls"
      % (len(catch_rows), len(ctrl_rows)))

lines = ["Prefix(owl:=<http://www.w3.org/2002/07/owl#>)", "Ontology(<http://example.org/disjtest>"]
for i, r in enumerate(cand):
    nid = "<http://example.org/NEW_%d>" % i
    parts = [IRI(r["A1_curie"])] + ["ObjectSomeValuesFrom(%s %s)" % (IRI(p), IRI(f))
                                    for p, f in sorted(LD[r["curie"]]["diff"])]
    lines.append("Declaration(Class(%s))" % nid)
    lines.append("EquivalentClasses(%s ObjectIntersectionOf(%s))" % (nid, " ".join(parts)))
for a, b in sorted(SAFE):
    lines.append("DisjointClasses(%s %s)" % (IRI(a), IRI(b)))
lines.append(")")

with tempfile.TemporaryDirectory() as td:
    ofn, unsat, out = (os.path.join(td, f) for f in ("cand.ofn", "unsat.owl", "reasoned.owl"))
    open(ofn, "w").write("\n".join(lines))
    cmd = ["java", "-Xmx8g", "-jar", ROBOT, "merge", "--input", CL_OWL, "--input", ofn,
           "reason", "--reasoner", "ELK", "--dump-unsatisfiable", unsat, "--output", out]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    blob = p.communicate()[0].decode("utf8", "replace")
    n_unsat = 0
    if os.path.exists(unsat):
        body = open(unsat).read()
        n_unsat = body.count("NEW_")
        print("  ELK: %d unsatisfiable (existing CL classes affected: %s)"
              % (n_unsat, "YES -- axioms invalid!" if "CL_" in body else "none"))
    else:
        print("  ELK: NOTHING unsatisfiable -> injected disjointness fires on 0/%d candidates"
              % len(cand))
    print("  (robot exit=%s)" % p.returncode)

# ---------------------------------------------------------------- D4
print("\n=== D4  could sufficient conditions (marker -> type) rescue it? ===")
errs = [r for r in rows if r.get("A1") == "miss"]
hasd = lambda r, pre: any(f.startswith(pre) for _, f in LD.get(r["curie"], {}).get("diff", set()))
sib = 0
for r in errs:
    p = r.get("A1_curie")
    if p and (PARENT.get(p, set()) & set(r["gold"].split(";")) or
              PARENT.get(p, set()) & PARENT.get(r["curie"], set())):
        sib += 1
d4 = {"errors": len(errs),
      "any_differentium": sum(bool(LD.get(r["curie"], {}).get("diff")) for r in errs),
      "pro_marker": sum(hasd(r, "PR:") for r in errs),
      "go_function": sum(hasd(r, "GO:") for r in errs),
      "uberon_location": sum(hasd(r, "UBERON:") for r in errs),
      "sibling_level": sib}
for k in ("any_differentium", "pro_marker", "go_function", "uberon_location", "sibling_level"):
    print("  %-20s %3d / %d (%.0f%%)" % (k, d4[k], len(errs), 100.0 * d4[k] / len(errs)))

json.dump({"D1": {"pool": POOL, "safe_pairs": len(SAFE), "blocked_pairs": len(UNSAFE)},
           "D2": d2, "D3": {"candidates": len(cand), "errors": len(catch_rows),
                            "controls": len(ctrl_rows), "unsatisfiable": n_unsat},
           "D4": d4},
          open(os.path.join(RES, "disjointness_experiment.json"), "w"), indent=2)
print("\nwrote results/disjointness_experiment.json")
