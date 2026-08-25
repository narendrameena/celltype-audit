#!/usr/bin/env python3
"""Select and download the CELLxGENE Discover datasets the survey audits.

The survey answers "does this generalise beyond the two atlases", and the strength of that
answer is the size and breadth of the population it runs on. The first pass used twelve
datasets picked for tissue breadth -- one per tissue, at most two per collection -- which
was the right shape for a demonstration and the wrong shape for an estimate. This selects
every dataset that CAN be audited and downloads it, so the population is defined by the
eligibility criteria rather than by a sampling decision.

Eligibility, in the order applied:

  organism      exactly Homo sapiens; the reference is human
  disease       exactly {normal}; a diseased cluster may legitimately fail to match a
                healthy reference profile, which is not an annotation error
  stage         not embryonic or fetal; the reference is adult-dominated, so a fetal
                cluster would disagree for developmental reasons
  tissue        at least one UBERON term the WMG reference actually carries, else there
                is nothing to score against
  asset         an H5AD to download
  tombstone     excluded

`is_primary_data` is deliberately NOT filtered on. It marks whether cells are original to
the dataset rather than re-analysed from another, which is a provenance question and not an
auditability one: a re-analysed dataset still asserts CL terms, and those assertions are
exactly what the audit tests. Filtering on it would have dropped seven of the twelve
datasets already published, changing the population rather than growing it. Duplication is
handled where it belongs -- reported per collection, as Supplementary Table 1 already does
for the two Tabula Sapiens slices.

Downloads are resumable and skip anything already present at the expected size, so this can
be re-run after an interruption. Smallest first, so the survey has something to chew on
early rather than blocking on a 45 GB file.

Usage:
    python benchmark/fetch_cxg_survey.py --list          # selection only, no download
    python benchmark/fetch_cxg_survey.py --max-gb 200    # cap the total
    python benchmark/fetch_cxg_survey.py
"""
import argparse
import json
import os
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
DATA = os.path.abspath(os.path.join(HERE, "..", "..", "cxg_data"))
INDEX = "https://api.cellxgene.cziscience.com/curation/v1/datasets"
WMG = "https://api.cellxgene.cziscience.com/wmg/v2/primary_filter_dimensions"
FETAL = ("embryonic", "fetal", "carnegie", "gestation", "post-fertilization")


def _get(url):
    with urllib.request.urlopen(url, timeout=300) as r:
        return json.load(r)


def wmg_tissues():
    """UBERON ids the expression reference actually carries, for human."""
    d = _get(WMG)
    return {list(t)[0] for t in d["tissue_terms"]["NCBITaxon:9606"]}


def eligible(index=None, tissues=None):
    index = index if index is not None else _get(INDEX)
    tissues = tissues if tissues is not None else wmg_tissues()
    out = []
    for r in index:
        if r.get("tombstone"):
            continue
        if [o["label"] for o in r.get("organism") or []] != ["Homo sapiens"]:
            continue
        if {x["label"] for x in r.get("disease") or []} != {"normal"}:
            continue
        stage = " ".join(x["label"].lower() for x in r.get("development_stage") or [])
        if any(k in stage for k in FETAL):
            continue
        assets = [a for a in r.get("assets") or [] if a.get("filetype") == "H5AD"]
        tis = {t["ontology_term_id"] for t in r.get("tissue") or []} & tissues
        if not assets or not tis:
            continue
        out.append({"dataset_id": r["dataset_id"], "url": assets[0]["url"],
                    "filesize": assets[0]["filesize"], "tissues": sorted(tis),
                    "cells": r.get("cell_count") or 0,
                    "collection_id": r["collection_id"],
                    "title": (r.get("title") or "")[:80]})
    return sorted(out, key=lambda c: c["filesize"])


def fetch(c, dest):
    """Resumable; a complete file at the expected size is left alone."""
    if os.path.exists(dest) and os.path.getsize(dest) == c["filesize"]:
        return "have"
    tmp = dest + ".part"
    have = os.path.getsize(tmp) if os.path.exists(tmp) else 0
    req = urllib.request.Request(c["url"])
    if have:
        req.add_header("Range", "bytes=%d-" % have)
    with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "ab" if have else "wb") as fh:
        while True:
            b = r.read(1 << 22)
            if not b:
                break
            fh.write(b)
    if os.path.getsize(tmp) != c["filesize"]:
        raise IOError("size mismatch for %s" % c["dataset_id"])
    os.rename(tmp, dest)
    return "fetched"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true", help="print the selection and exit")
    ap.add_argument("--max-gb", type=float, help="stop once this much has been downloaded")
    a = ap.parse_args()

    cand = eligible()
    total = sum(c["filesize"] for c in cand) / 1e9
    print("%d eligible datasets | %.0f GB | %d tissues | %d collections\n"
          % (len(cand), total, len({t for c in cand for t in c["tissues"]}),
             len({c["collection_id"] for c in cand})))
    os.makedirs(RES, exist_ok=True)
    json.dump(cand, open(os.path.join(RES, "cxg_selection.json"), "w"), indent=1)
    if a.list:
        for c in cand:
            print("  %-38s %7.2f GB  %-9s %s"
                  % (c["dataset_id"], c["filesize"] / 1e9, c["tissues"][0], c["title"][:44]))
        return

    os.makedirs(DATA, exist_ok=True)
    got = 0.0
    for i, c in enumerate(cand, 1):
        dest = os.path.join(DATA, "%s.h5ad" % c["dataset_id"])
        if a.max_gb and got >= a.max_gb:
            print("stopping at --max-gb %.0f" % a.max_gb)
            break
        try:
            how = fetch(c, dest)
        except Exception as exc:                       # one bad asset must not end the run
            print("  [%3d/%d] FAIL %s: %s" % (i, len(cand), c["dataset_id"][:8], exc),
                  flush=True)
            continue
        got += c["filesize"] / 1e9
        print("  [%3d/%d] %-7s %s  %6.2f GB  (%.0f GB total)"
              % (i, len(cand), how, c["dataset_id"][:8], c["filesize"] / 1e9, got), flush=True)


if __name__ == "__main__":
    main()
