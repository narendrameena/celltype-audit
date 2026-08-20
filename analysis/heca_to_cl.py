#!/usr/bin/env python3
"""Map hECA/uHAF cell types to the Cell Ontology BY EXPRESSION (no vocabulary in the loop).

Input : results/heca_markers_<Organ>.json  (markers computed from the data by heca_markers.py)
Method: markers + organ -> CELLxGENE WMG API -> rank CL terms by mean over markers of (pc * me)
Output: results/heca_to_cl_<Organ>.json

Usage: python benchmark/heca_to_cl.py <Organ> [<Organ> ...]     e.g. Pancreas
       python benchmark/heca_to_cl.py --all
"""
import glob
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
CACHE = os.path.join(HERE, ".wmg_cache")
os.makedirs(CACHE, exist_ok=True)
WMG = "https://api.cellxgene.cziscience.com/wmg/v2/query"
DIMS = "https://api.cellxgene.cziscience.com/wmg/v2/primary_filter_dimensions"
HUMAN = "NCBITaxon:9606"
MINCELLS = int(os.environ.get("MINCELLS", 100))
TOPN = int(os.environ.get("TOPN", 5))

# hECA organ file name -> CELLxGENE WMG tissue label (only where they differ / need a choice)
ORGAN2TISSUE = {
    # explicit hECA organ -> CELLxGENE WMG tissue label. Verified against the 64 human tissues
    # the API exposes; a wrong guess here silently maps a whole organ to the wrong reference.
    "Adipose": "adipose tissue", "Adrenal_gland": "adrenal gland", "Bladder": "urinary bladder",
    "Blood": "blood", "Bone_marrow": "bone marrow", "Brain": "brain", "Breast": "breast",
    "Bronchi": "lung",                       # bronchi are part of the lung
    "Eye": "eye", "Femur": "bone marrow",    # hECA femur samples are marrow
    "Gallbladder": "gallbladder", "Heart": "heart", "Intestine": "intestine",
    "Kidney": "kidney", "Liver": "liver", "Lung": "lung", "Lymph_node": "lymph node",
    "Muscle": "musculature", "Nose": "nose", "Oesophagus": "esophagus", "Ovary": "ovary",
    "Pancreas": "pancreas", "Placenta": "placenta", "Pleura": "pleura",
    "Prostate": "prostate gland", "Rib": "bone marrow",
    "Salivary_gland": "exocrine gland", "Skin": "skin of body", "Spinal_cord": "spinal cord",
    "Spleen": "spleen", "Stomach": "stomach", "Testis": "testis", "Ureter": "ureter",
    "Uterine_tube": "fallopian tube", "Uterus": "uterus", "Vessel": "vasculature",
    # No CELLxGENE tissue equivalent -> deliberately NOT mapped rather than mapped to a wrong
    # proxy. Reported as a reference-coverage limit, not a method failure.
    "Thymus": None, "Thyroid": None, "Trachea": None, "Oral_cavity": None,
}


def _fetch(url, body=None, tries=4):
    # NB: this used to key the cache on hash(), which Python salts per process, so every
    # run computed fresh keys, missed every entry and re-downloaded the lot. sha1 is
    # stable across processes; entries written before this fix are unreachable by key.
    tag = hashlib.sha1(("%s|%s" % (url, json.dumps(body, sort_keys=True))).encode()).hexdigest()[:16]
    key = os.path.join(CACHE, tag + ".json")
    if os.path.exists(key):
        try:
            return json.load(open(key))
        except Exception:
            pass
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
        except Exception:
            if a == tries - 1:
                return None
            time.sleep(4)


dims = _fetch(DIMS)
SYM2ENS, TISSUE = {}, {}
for e in dims["gene_terms"][HUMAN]:
    for k, v in e.items():
        SYM2ENS.setdefault(v.upper(), k)
for e in dims["tissue_terms"][HUMAN]:
    for k, v in e.items():
        TISSUE[v] = k


def tissue_for(organ):
    if organ in ORGAN2TISSUE:
        cand = ORGAN2TISSUE[organ]
        if cand is None:                              # explicitly no equivalent
            return None, None
    else:
        cand = organ.replace("_", " ").lower()
    if cand in TISSUE:
        return cand, TISSUE[cand]
    return None, None                                 # never guess by substring


def rank_cl(markers, uberon, topn=TOPN):
    ids = sorted({SYM2ENS[m] for m in markers if m in SYM2ENS})
    if not ids:
        return None
    r = _fetch(WMG, {"filter": {"gene_ontology_term_ids": list(ids),
                                "organism_ontology_term_id": HUMAN}, "is_rollup": True})
    if not r:
        return None
    es = r["expression_summary"]
    lbl = r["term_id_labels"]["cell_types"].get(uberon, {})
    sc = defaultdict(list)
    for gid in ids:
        for ct, v in es.get(gid, {}).get(uberon, {}).items():
            a = v.get("aggregated") or {}
            sc[ct].append(a.get("pc", 0.0) * a.get("me", 0.0))
    out = []
    for ct, vals in sc.items():
        if ct == "CL:0000000":
            continue
        info = (lbl.get(ct, {}).get("aggregated") or {})
        if info.get("total_count", 0) < MINCELLS:
            continue
        out.append({"curie": ct, "label": info.get("name", ct),
                    "score": round(sum(vals) / len(ids), 4), "n": info.get("total_count", 0)})
    out.sort(key=lambda x: -x["score"])
    return out[:topn] or None


def run(organ):
    p = os.path.join(RES, "heca_markers_%s.json" % organ)
    if not os.path.exists(p):
        print("  no marker file for %s" % organ)
        return None
    M = json.load(open(p))
    tname, uberon = tissue_for(organ)
    if not uberon:
        print("  %-16s NO CELLxGENE tissue match -> skipped" % organ)
        return None
    res, miss = {}, 0
    for t, v in M["types"].items():
        genes = [m["gene"] for m in v["markers"]]
        hits = rank_cl(genes, uberon)
        if hits:
            res[t] = {"n_cells": v["n_cells"], "markers": genes, "cl": hits}
        else:
            miss += 1
    json.dump({"organ": organ, "tissue": tname, "uberon": uberon, "types": res},
              open(os.path.join(RES, "heca_to_cl_%s.json" % organ), "w"), indent=1)
    print("  %-16s tissue=%-18s mapped %d/%d types" % (organ, tname, len(res), len(M["types"])))
    return res


if __name__ == "__main__":
    args = sys.argv[1:]
    if not args:
        raise SystemExit(__doc__)
    if args[0] == "--all":
        organs = sorted(os.path.basename(p)[len("heca_markers_"):-len(".json")]
                        for p in glob.glob(os.path.join(RES, "heca_markers_*.json")))
    else:
        organs = args
    print("mapping %d organ(s) to CL by expression ...\n" % len(organs))
    allres = {}
    for o in organs:
        r = run(o)
        if r:
            allres[o] = r
    if len(organs) == 1 and allres:
        o = organs[0]
        print("\n%s — top CL term per uHAF cell type\n" % o)
        print("%-40s %-8s %-13s %s" % ("uHAF label", "n_cells", "CL", "CL label (expression-derived)"))
        print("-" * 100)
        for t, v in sorted(allres[o].items(), key=lambda x: -x[1]["n_cells"]):
            top = v["cl"][0]
            print("%-40s %-8d %-13s %s" % (t[:40], v["n_cells"], top["curie"], top["label"][:44]))
