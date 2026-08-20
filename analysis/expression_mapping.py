#!/usr/bin/env python3
"""Map atlas cell-type labels to the Cell Ontology by EXPRESSION, not by name.

Motivation: lexical mapping asks "do the names agree?"; biology asks "do the cells look the
same?". On a real atlas vocabulary (uHAF/hECA v2.0) the lexical method mapped `Alveolar cell`
(lung) to a *mammary gland* term, because CL synonyms are not organ-disambiguated. Expression
has no such failure mode.

Method (name-free):
  uHAF label -> its curated marker genes + its organ
    -> CELLxGENE WMG API (whole corpus, CL-annotated by mandate)
    -> for each CL term in that tissue, score = mean over markers of (pc * me)
       pc = fraction of that cell type's cells expressing the gene
       me = mean expression among expressing cells
    -> top-ranked CL term = the biological match

Evaluated WITHOUT a gold standard, exactly as the lexical arm was: map both ends of each uHAF
parent-child edge and ask whether CL's is_a closure preserves the relation.

Usage: python benchmark/expression_mapping.py            (needs network; ~1 API call per gene set)
Env:   N (labels to try, default 120), MINCELLS (default 100)
"""
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict, deque

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(HERE, "results")
WMG = "https://api.cellxgene.cziscience.com/wmg/v2/query"
DIMS = "https://api.cellxgene.cziscience.com/wmg/v2/primary_filter_dimensions"
HUMAN = "NCBITaxon:9606"
N = int(os.environ.get("N", 120))
MINCELLS = int(os.environ.get("MINCELLS", 100))
CACHE = os.path.join(HERE, ".wmg_cache")
os.makedirs(CACHE, exist_ok=True)
os.makedirs(RES, exist_ok=True)
SHORT = lambda i: i.rsplit("/", 1)[-1].replace("_", ":")

# uHAF organ sheet -> CELLxGENE tissue label
ORGAN = {"Bronchi": "lung", "Lung": "lung", "Blood": "blood", "Bone marrow": "bone marrow",
         "Brain": "brain", "Breast": "breast", "Heart": "heart", "Kidney": "kidney",
         "Liver": "liver", "Pancreas": "pancreas", "Skin": "skin", "Stomach": "stomach",
         "Testis": "testis", "Thymus": "thymus", "Spleen": "spleen", "Uterus": "uterus",
         "Ovary": "ovary", "Prostate": "prostate", "Eye": "eye", "Placenta": "placenta",
         "Small intestine": "small intestine", "Large intestine": "large intestine",
         "Muscle": "muscle organ", "Thyroid": "thyroid gland", "Bladder": "urinary bladder",
         "Lymph node": "lymph node", "Oesophagus": "esophagus", "Vessel": "blood vessel"}


def _get(url, body=None, tries=4):
    key = os.path.join(CACHE, "%x.json" % abs(hash((url, json.dumps(body, sort_keys=True)))))
    if os.path.exists(key):
        return json.load(open(key))
    for a in range(tries):
        try:
            if body is None:
                r = json.loads(urllib.request.urlopen(url, timeout=300).read())
            else:
                req = urllib.request.Request(url, data=json.dumps(body).encode(),
                                             headers={"Content-Type": "application/json"})
                r = json.loads(urllib.request.urlopen(req, timeout=300).read())
            json.dump(r, open(key, "w"))
            return r
        except urllib.error.HTTPError:
            if a == tries - 1:
                return None
            time.sleep(4)
        except Exception:
            if a == tries - 1:
                return None
            time.sleep(4)


print("loading CELLxGENE dimensions ...", flush=True)
dims = _get(DIMS)
SYM2ENS, TISSUE = {}, {}
for e in dims["gene_terms"][HUMAN]:
    for k, v in e.items():
        SYM2ENS.setdefault(v.upper(), k)
for e in dims["tissue_terms"][HUMAN]:
    for k, v in e.items():
        TISSUE[v] = k
print("  %d human genes, %d tissues" % (len(SYM2ENS), len(TISSUE)), flush=True)

# ---- CL structure (for the hierarchy-preservation test)
CLJ = os.environ.get("CL_JSON", os.path.join(ROOT, "cl-full.json"))
g = json.load(open(CLJ))["graphs"][0]
LAB = {SHORT(n["id"]): n["lbl"] for n in g["nodes"] if n.get("lbl")}
PAR = defaultdict(set)
for e in g["edges"]:
    if e.get("pred") == "is_a":
        PAR[SHORT(e["sub"])].add(SHORT(e["obj"]))


def anc(c):
    s, q = set(), deque([c])
    while q:
        for p in PAR.get(q.popleft(), ()):
            if p not in s:
                s.add(p)
                q.append(p)
    return s


# ---- uHAF
UP = os.path.join(RES, "uhaf_parsed.json")
U = json.load(open(UP))
LABELS, EDGES, MARK = U["labels"], U["edges"], U["markers"]
GENE = re.compile(r"\b[A-Z][A-Z0-9\-]{1,14}\b")
STOP = {"AND", "OR", "NOT", "THE"}


def genes_of(L):
    return sorted({x for x in GENE.findall(" ".join(MARK.get(L, []))) if x not in STOP and x in SYM2ENS})


def organ_of(L):
    for o in LABELS.get(L, []):
        if o in ORGAN and ORGAN[o] in TISSUE:
            return ORGAN[o]
    return None


def expression_match(L, topn=5):
    """Return ranked CL candidates for a uHAF label, from expression alone."""
    gs, org = genes_of(L), organ_of(L)
    if not gs or not org:
        return None
    ids = sorted(SYM2ENS[x] for x in gs)
    r = _get(WMG, {"filter": {"gene_ontology_term_ids": ids,
                              "organism_ontology_term_id": HUMAN}, "is_rollup": True})
    if not r:
        return None
    u = TISSUE[org]
    es = r["expression_summary"]
    lbl = r["term_id_labels"]["cell_types"].get(u, {})
    sc = defaultdict(list)
    for gid in ids:
        for ct, v in es.get(gid, {}).get(u, {}).items():
            a = v.get("aggregated") or {}
            sc[ct].append(a.get("pc", 0.0) * a.get("me", 0.0))
    out = []
    for ct, vals in sc.items():
        if ct == "CL:0000000":
            continue
        info = (lbl.get(ct, {}).get("aggregated") or {})
        if info.get("total_count", 0) < MINCELLS:
            continue
        out.append({"score": round(sum(vals) / len(ids), 4), "curie": ct,
                    "label": info.get("name", ct), "n": info.get("total_count", 0)})
    out.sort(key=lambda x: -x["score"])
    return out[:topn] or None


cands = [L for L in LABELS if genes_of(L) and organ_of(L)][:N]
print("uHAF labels with markers AND a mappable organ: %d (querying %d)" % (
    len([L for L in LABELS if genes_of(L) and organ_of(L)]), len(cands)), flush=True)

EXPR = {}
for i, L in enumerate(cands):
    m = expression_match(L)
    if m:
        EXPR[L] = m
    if (i + 1) % 20 == 0:
        print("  %d/%d ..." % (i + 1, len(cands)), flush=True)
print("expression matches obtained: %d" % len(EXPR))

LEXP = os.path.join(RES, "uhaf_to_cl_DRAFT_unvalidated.json")
LEX = json.load(open(LEXP)) if os.path.exists(LEXP) else {}


def hier(pick):
    both = ok = 0
    for a, b in EDGES:
        A, B = pick(a), pick(b)
        if not A or not B:
            continue
        both += 1
        ok += (A in anc(B) or A == B)
    return ok, both


lex_pick = lambda L: (LEX.get(L, {}).get("cl") or [None])[0]
exp_pick = lambda L: EXPR[L][0]["curie"] if L in EXPR else None
sub = set(EXPR)
lex_sub = lambda L: lex_pick(L) if L in sub else None

print("\n%-46s %s" % ("hierarchy preservation (uHAF edge -> CL is_a)", "result"))
print("-" * 66)
for name, fn in [("name-based, all labels", lex_pick),
                 ("name-based, restricted to expression subset", lex_sub),
                 ("EXPRESSION-based", exp_pick)]:
    ok, both = hier(fn)
    print("%-46s %3d/%3d = %.1f%%" % (name, ok, both, 100.0 * ok / max(1, both)))

agree = diff = 0
examples = []
for L in EXPR:
    lx, ex = lex_pick(L), exp_pick(L)
    if not lx:
        continue
    if lx == ex:
        agree += 1
    else:
        diff += 1
        if len(examples) < 12:
            examples.append((L, lx, LAB.get(lx, "?"), ex, LAB.get(ex, EXPR[L][0]["label"])))
print("\nname vs expression on the same labels: agree %d, differ %d" % (agree, diff))
for L, lx, ll, ex, el in examples:
    print("   %-26s name=%-11s %-30s  expr=%-11s %s" % (L[:26], lx, ll[:30], ex, el[:30]))

json.dump({"expression": EXPR}, open(os.path.join(RES, "uhaf_expression_mapping.json"), "w"), indent=1)
print("\nwrote results/uhaf_expression_mapping.json")
