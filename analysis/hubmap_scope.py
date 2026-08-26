#!/usr/bin/env python3
"""Does HuBMAP publish cell-type annotations the audit could test? Measured, not assumed.

The external check in this paper rests on one atlas, and a referee is entitled to ask why
not more. The answer is that an audit needs an ASSERTED cell-type label to test against
expression, and the largest public single-cell repository outside CELLxGENE does not
publish one. That is a claim about a live resource, so it is derived here rather than
stated.

The method: query HuBMAP's portal search for published, centrally processed datasets whose
data are public (no human genetic sequences), and enumerate the `obs` columns their
distributed AnnData objects carry -- readable from the zarr paths in the file manifest
without downloading anything. A column that asserts a cell type is one a curator or a
classifier wrote: cell_type, predicted_label, azimuth, annotation. Leiden clusters are not
such a column: an unsupervised cluster number asserts nothing that could be wrong.

Writes results/hubmap_scope.json.
"""
import json
import os
import re
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
SEARCH = "https://search.api.hubmapconsortium.org/v3/portal/search"
# the assets CDN refuses a bare request but serves public data to a browser agent
UA = "Mozilla/5.0 (X11; Linux x86_64)"
ASSERTS = re.compile(r"cell.?type|predicted|azimuth|annotation|^label$", re.I)


def _get(url):
    return urllib.request.urlopen(
        urllib.request.Request(url, headers={"User-Agent": UA}), timeout=300).read()


def post(body):
    """A large result is handed back as a presigned URL, either in the body or via a 303."""
    req = urllib.request.Request(SEARCH, data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json", "User-Agent": UA})
    try:
        raw = urllib.request.urlopen(req, timeout=300).read()
    except urllib.error.HTTPError as exc:
        if exc.code not in (301, 302, 303, 307, 308):
            raise
        # the presigned URL arrives in the BODY of the 303, not in a Location header
        loc = exc.headers.get("Location") or exc.read().decode().strip()
        return json.loads(_get(loc))
    if raw[:5] == b"https":
        raw = _get(raw.decode().strip())
    return json.loads(raw)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    d = post({"size": n,
              "query": {"bool": {"must": [
                  {"term": {"entity_type.keyword": "Dataset"}},
                  {"term": {"status.keyword": "Published"}},
                  {"term": {"creation_action.keyword": "Central Process"}},
                  {"term": {"contains_human_genetic_sequences": False}}]}},
              "_source": ["hubmap_id", "uuid", "dataset_type",
                          "origin_samples.mapped_organ", "files.rel_path"]})
    hits = d["hits"]["hits"]
    cols, organs, per_col = {}, set(), {}
    for h in hits:
        s = h["_source"]
        organs.add((s.get("origin_samples") or [{}])[0].get("mapped_organ"))
        for f in (s.get("files") or []):
            rp = str(f.get("rel_path", ""))
            if "/obs/" not in rp or ".zarr" not in rp:
                continue
            seg = rp.split("/obs/", 1)[1].split("/")[0]
            if not seg or seg.startswith("."):
                continue
            cols[seg] = cols.get(seg, 0) + 1
            per_col.setdefault(seg, set()).add(s["hubmap_id"])
    asserting = sorted(c for c in cols if ASSERTS.search(c))
    out = {"datasets_examined": len(hits),
           "total_published": d["hits"]["total"]["value"],
           "organs": sorted(str(o) for o in organs if o),
           "obs_columns": dict(sorted(cols.items(), key=lambda kv: -kv[1])),
           "columns_asserting_a_cell_type": asserting}
    json.dump(out, open(os.path.join(RES, "hubmap_scope.json"), "w"), indent=1)
    print("public, centrally processed HuBMAP datasets examined : %d" % len(hits))
    print("organs represented                                   : %d" % len(out["organs"]))
    print("\nobs columns their distributed objects carry:")
    for c, k in list(out["obs_columns"].items())[:20]:
        print("   %-24s in %3d datasets" % (c[:24], k))
    print("\ncolumns asserting a cell type: %s"
          % (", ".join(asserting) if asserting else "NONE"))
    print("wrote results/hubmap_scope.json")


if __name__ == "__main__":
    main()
