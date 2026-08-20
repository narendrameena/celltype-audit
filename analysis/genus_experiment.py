#!/usr/bin/env python3
"""Go/no-go experiment: can we beat the 50.7% genus-derivation baseline?

Arms (all scored with the SAME rubric as benchmark B2: exact / ancestor / defaulted / miss,
on the SAME seeded sample, so numbers are directly comparable to the published 50.7%):

  A0  keyword         14 hand-written rules  (`agent._derive_parent`)      -- the baseline
  A1  suffix-ladder   longest proper suffix of the label that grounds to a CL term
  A2  defn-phrase     leading noun phrase of the curated textual definition ("A <genus> that ...")
  A3  differentia-nbr genera of terms sharing >=1 differentium             -- known-weak, for the record
  L1  lexical-exact    longest label suffix that EXACTLY matches a CL label/synonym (offline)
  L2  lexical-vote     L1 + synonym surfaces + head-noun preference + surface voting (offline)

Rubric note: a pick of the root `cell` (CL:0000000) is scored "defaulted", never "ancestor" --
a root pick is technically a superclass of everything and would trivially inflate every arm.

Leakage guards: a candidate is rejected if it grounds to the target term itself or to any
descendant of it; A1 never uses the full label; A3 excludes the target and its subtree.

Usage:  CL_JSON=cl-full.json python benchmark/genus_experiment.py
Env:    N (sample size, default 400), WORKERS (default 12), ARMS (comma list, default all)
"""
import csv
import json
import os
import random
import re
import sys
from collections import defaultdict, deque
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
os.environ.setdefault("CELLSCRIBE_CACHE", os.path.join(HERE, ".bench_cache"))

from cellscribe.agent import _derive_parent
from cellscribe.tools.ontology import OLSSearchTool

CL_JSON = os.environ.get("CL_JSON", os.path.join(ROOT, "cl-full.json"))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)
WORKERS = int(os.environ.get("WORKERS", 12))
N = int(os.environ.get("N", 400))
ARMS = os.environ.get("ARMS", "A0,A1,A2,A3,L1,L2").split(",")
SHORT = lambda iri: iri.rsplit("/", 1)[-1].replace("_", ":")


# ---------------------------------------------------------------- gold standard
def load_cl(path):
    g = json.load(open(path))["graphs"][0]
    label, defn, syn = {}, {}, defaultdict(list)
    for n in g["nodes"]:
        cur = SHORT(n.get("id", ""))
        if n.get("lbl"):
            label[cur] = n["lbl"]
        meta = n.get("meta", {}) or {}
        for sv in meta.get("synonyms", []) or []:
            if sv.get("val"):
                syn[cur].append(sv["val"])
        if (meta.get("definition") or {}).get("val"):
            defn[cur] = meta["definition"]["val"]
    parents, children = defaultdict(set), defaultdict(set)
    for e in g["edges"]:
        if e.get("pred") == "is_a":
            parents[SHORT(e["sub"])].add(SHORT(e["obj"]))
            children[SHORT(e["obj"])].add(SHORT(e["sub"]))
    ld = defaultdict(lambda: {"genus": set(), "diff": set()})
    for a in g.get("logicalDefinitionAxioms", []):
        cur = SHORT(a["definedClassId"])
        if not cur.startswith("CL:"):
            continue
        for gid in a.get("genusIds", []):
            ld[cur]["genus"].add(SHORT(gid))
        for r in a.get("restrictions") or []:
            if r and r.get("fillerId"):
                ld[cur]["diff"].add((SHORT(r.get("propertyId", "")), SHORT(r["fillerId"])))
    return label, defn, dict(syn), dict(parents), dict(children), dict(ld)


CL_ROOT = "CL:0000000"


def _head(s):
    t = s.lower().split()
    return t[-1] if t else ""


def _ladder(s):
    """Progressively shorter PROPER suffixes, plus the pre-'of' head. Longest first."""
    s = s.lower().strip()
    toks = s.split()
    out = []
    if " of " in s:
        h = s.split(" of ")[0].strip()
        if h and h != s:
            out.append(h)
    for i in range(1, len(toks)):
        x = _STOP.sub("", " ".join(toks[i:]))
        if x and x not in out:
            out.append(x)
    return out


def _closure(seed, adj):
    seen, q = set(), deque([seed])
    while q:
        for p in adj.get(q.popleft(), ()):
            if p not in seen:
                seen.add(p)
                q.append(p)
    return seen


print("loading gold standard ...", flush=True)
LABEL, DEFN, SYN, PARENTS, CHILDREN, LD = load_cl(CL_JSON)

# exact lexical index: every CL label and synonym -> the terms bearing it (offline, no network)
LEX = defaultdict(set)
for _c in LABEL:
    if _c.startswith("CL:"):
        LEX[LABEL[_c].lower()].add(_c)
        for _s in SYN.get(_c, ()):
            LEX[_s.lower().strip()].add(_c)
cl_terms = [c for c in LABEL if c.startswith("CL:")]
print("CL labelled: %d | with logical def: %d | lexical index: %d strings"
      % (len(cl_terms), len(LD), len(LEX)), flush=True)

# index: differentium -> terms carrying it (for A3)
BY_DIFF = defaultdict(set)
for c, v in LD.items():
    for d in v["diff"]:
        BY_DIFF[d].add(c)

ols = OLSSearchTool()


def cl_hits(q):
    try:
        return [h for h in ols(q, ontology="cl", rows=5, offline=False) if h.curie.startswith("CL:")]
    except Exception:
        return []


# ---------------------------------------------------------------- candidate generators
_STOP = re.compile(r"^[\s,;:]+|[\s,;:]+$")


def suffixes(lbl):
    """Progressively shorter PROPER suffixes of the label (longest first).

    'striatal parvalbumin-positive GABAergic interneuron'
      -> 'parvalbumin-positive GABAergic interneuron', 'GABAergic interneuron', 'interneuron'
    Also yields the pre-'of' head for 'X of Y' labels ('epithelial cell of lung' -> 'epithelial cell').
    """
    toks = lbl.split()
    out = []
    if " of " in lbl:
        head = _STOP.sub("", lbl.split(" of ")[0])
        if head and head != lbl:
            out.append(head)
    for i in range(1, len(toks)):
        s = _STOP.sub("", " ".join(toks[i:]))
        if s and s not in out:
            out.append(s)
    return out


def defn_phrase(text):
    """Leading genus phrase of a CL house-style definition: 'A(n) <genus> that/which/, ...'."""
    if not text:
        return None
    t = text.strip()
    t = re.sub(r"^(An?|The)\s+", "", t, flags=re.I)
    m = re.split(r"\b(that|which|whose|with|in which|located|expressing|derived|capable)\b|[,.;(]", t, 1)
    head = _STOP.sub("", m[0]) if m else None
    if not head or len(head.split()) > 8:
        return None
    return head


# ---------------------------------------------------------------- scoring
def score_pick(cur, picked_curie, picked_query, gold, anc, defaulted):
    if defaulted or picked_curie == CL_ROOT:
        return "defaulted"
    if picked_curie and picked_curie in gold:
        return "exact"
    if picked_curie and picked_curie in anc:
        return "ancestor"
    return "miss"


def run_one(cur):
    gold = {g for g in (LD.get(cur, {}).get("genus", set()) | PARENTS.get(cur, set()))
            if g.startswith("CL:")}
    if not gold:
        return None
    anc = _closure(cur, PARENTS)
    banned = {cur} | _closure(cur, CHILDREN)          # leakage guard
    lbl = LABEL[cur]
    row = {"curie": cur, "label": lbl, "gold": ";".join(sorted(gold))}

    def ground_first(queries):
        """First query that grounds to a CL term outside the banned set."""
        for q in queries:
            for h in cl_hits(q):
                if h.curie not in banned:
                    return h.curie, q
        return None, None

    # --- A0 baseline: 14 keyword rules
    if "A0" in ARMS:
        pq = _derive_parent(lbl + " " + DEFN.get(cur, ""))
        hits = cl_hits(pq)
        P = hits[0].curie if hits else None
        row["A0_query"], row["A0_curie"] = pq, P
        row["A0"] = score_pick(cur, P, pq, gold, anc, defaulted=(pq == "cell"))

    # --- A1 suffix ladder (label only; no curated definition, no LLM, no reasoner)
    if "A1" in ARMS:
        P, q = ground_first(suffixes(lbl))
        row["A1_query"], row["A1_curie"] = q, P
        row["A1"] = score_pick(cur, P, q, gold, anc, defaulted=(P is None))

    # --- A2 leading phrase of the curated textual definition
    if "A2" in ARMS:
        ph = defn_phrase(DEFN.get(cur))
        qs = ([ph] + suffixes(ph)) if ph else []
        P, q = ground_first(qs)
        row["A2_query"], row["A2_curie"] = q, P
        row["A2"] = score_pick(cur, P, q, gold, anc, defaulted=(P is None))

    # --- A3 differentia-matched neighbours (structural; known weak)
    if "A3" in ARMS:
        nbr = set()
        for d in LD.get(cur, {}).get("diff", set()):
            nbr |= BY_DIFF.get(d, set())
        nbr -= banned
        cand = defaultdict(float)
        mine = LD.get(cur, {}).get("diff", set())
        for x in nbr:
            xd = LD[x]["diff"]
            w = len(mine & xd) / float(len(mine | xd)) if (mine | xd) else 0.0
            for g in LD[x]["genus"]:
                cand[g] += w
        P = max(cand, key=cand.get) if cand else None
        row["A3_query"], row["A3_curie"] = ("%d nbrs" % len(nbr)), P
        row["A3"] = score_pick(cur, P, None, gold, anc, defaulted=(P is None))

    # --- L1/L2 exact-lexical ladder: OFFLINE, no OLS. A suffix only counts if it EXACTLY
    #     matches a CL label or synonym -- fuzzy search is what manufactured the sibling errors.
    if "L1" in ARMS or "L2" in ARMS:
        surfaces = [lbl] + list(SYN.get(cur, ()))
        cands = []                                    # (surface_idx, ladder_idx, curie)
        for si, surf in enumerate(surfaces):
            for li, q in enumerate(_ladder(surf)):
                for m in LEX.get(q, ()):
                    if m not in banned and m != CL_ROOT:
                        cands.append((si, li, m))
        if "L1" in ARMS:                              # longest suffix of the primary label
            prim = [x for x in cands if x[0] == 0]
            P = min(prim, key=lambda x: x[1])[2] if prim else None
            row["L1_query"], row["L1_curie"] = "", P
            row["L1"] = score_pick(cur, P, None, gold, anc, defaulted=(P is None))
        if "L2" in ARMS:                              # + synonyms, head-noun preference, voting
            h = _head(lbl)
            pref = [x for x in cands if _head(LABEL.get(x[2], "")) == h]
            votes = {}
            for x in cands:
                votes[x[2]] = votes.get(x[2], 0) + 1
            pool = pref or cands
            P = max(pool, key=lambda x: (votes[x[2]], -x[1]))[2] if pool else None
            row["L2_query"], row["L2_curie"] = ("%d cands" % len(cands)), P
            row["L2"] = score_pick(cur, P, None, gold, anc, defaulted=(P is None))

    return row


random.seed(0)
sample = random.sample(cl_terms, min(N, len(cl_terms)))
print("running arms %s on n=%d (scored on those with a gold genus) ..." % (ARMS, len(sample)), flush=True)
with ThreadPoolExecutor(max_workers=WORKERS) as ex:
    rows = [r for r in ex.map(run_one, sample) if r is not None]

# ---------------------------------------------------------------- report
summary = {"n": len(rows), "cl_json": os.path.basename(CL_JSON), "arms": {}}
print("\nn = %d terms with a curated genus\n" % len(rows))
print("%-16s %8s %8s %10s %7s %10s" % ("arm", "exact", "ancestor", "defaulted", "miss", "HIER-VALID"))
print("-" * 64)
NAMES = {"A0": "A0 keyword", "A1": "A1 fuzzy-suffix", "A2": "A2 defn-phrase",
         "A3": "A3 differentia-nbr", "L1": "L1 lexical-exact", "L2": "L2 lexical-vote"}
for a in ARMS:
    if a not in rows[0]:
        continue
    cnt = defaultdict(int)
    for r in rows:
        cnt[r[a]] += 1
    n = len(rows)
    hv = (cnt["exact"] + cnt["ancestor"]) / float(n)
    summary["arms"][a] = {"exact": cnt["exact"], "ancestor": cnt["ancestor"],
                          "defaulted": cnt["defaulted"], "miss": cnt["miss"],
                          "exact_rate": round(cnt["exact"] / float(n), 4),
                          "hierarchically_valid_rate": round(hv, 4)}
    print("%-16s %8d %8d %10d %7d %9.1f%%" %
          (NAMES.get(a, a), cnt["exact"], cnt["ancestor"], cnt["defaulted"], cnt["miss"], 100 * hv))

# ---------------------------------------------------------------- cascades & oracle
def _fires(r, a):
    return r.get(a) not in (None, "defaulted") and bool(r.get(a + "_curie"))


def _hit(r, a):
    return r.get(a) in ("exact", "ancestor")


def cascade(rs, order):
    hit = cov = 0
    for r in rs:
        for a in order:
            if _fires(r, a):
                cov += 1
                hit += _hit(r, a)
                break
    return 100.0 * hit / len(rs), 100.0 * cov / len(rs)


CASCADES = [("A0",), ("A0", "A1"), ("L1", "A0", "A1"), ("L2", "A0", "A1"),
            ("L2", "A0", "A2", "A1")]
CASCADES = [c for c in CASCADES if all(a in ARMS for a in c)]
if CASCADES:
    print("\n%-40s %9s %10s %11s" % ("cascade (first arm that fires wins)", "coverage", "precision", "hier-valid"))
    print("-" * 74)
    for c in CASCADES:
        hv, cov = cascade(rows, c)
        print("%-40s %8.1f%% %9.1f%% %10.1f%%" % (" -> ".join(c), cov, hv * 100 / max(cov, 1e-9), hv))
        summary.setdefault("cascades", {})[" -> ".join(c)] = {
            "hierarchically_valid_rate": round(hv / 100, 4), "coverage": round(cov / 100, 4)}

    best = max(CASCADES, key=lambda c: cascade(rows, c)[0])
    random.seed(0)
    bs = sorted(cascade([rows[random.randrange(len(rows))] for _ in range(len(rows))], best)[0]
                for _ in range(2000))
    lo, hi = bs[50], bs[1949]
    print("\nBEST: %-30s %.1f%%  95%% CI [%.1f, %.1f]   (published baseline 50.7%%)"
          % (" -> ".join(best), cascade(rows, best)[0], lo, hi))
    summary["best_cascade"] = {"arms": list(best), "rate": round(cascade(rows, best)[0] / 100, 4),
                               "ci95": [round(lo / 100, 4), round(hi / 100, 4)]}
    oracle = 100.0 * sum(any(_hit(r, a) for a in ARMS) for r in rows) / len(rows)
    print("ORACLE over all arms (perfect selector): %.1f%%  -> %.1f pts remain in SELECTION"
          % (oracle, oracle - cascade(rows, best)[0]))
    summary["oracle"] = round(oracle / 100, 4)

with open(os.path.join(OUT, "genus_experiment.csv"), "w") as f:
    w = csv.DictWriter(f, fieldnames=sorted(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
json.dump(summary, open(os.path.join(OUT, "genus_experiment.json"), "w"), indent=2)
print("\nwrote results/genus_experiment.{csv,json}")
