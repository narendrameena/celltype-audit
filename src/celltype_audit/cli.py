"""Command line for celltype_audit.

    celltype-audit audit lung.h5ad --organ Lung -o annotations.json
    celltype-audit annotate new.h5ad --tissue UBERON:0002048 -o proposals.json
    celltype-audit resolve "Neutrophilic granulocyte" --organ Lung
    celltype-audit reference --tissues UBERON:0002048 --genes markers.txt -o ref
"""
import argparse
import json
import sys

from . import __version__


def _audit(a):
    from .audit import audit_h5ad
    rep = audit_h5ad(a.h5ad, organ=a.organ, tissue=a.tissue, min_cells=a.min_cells,
                     type_key=a.type_key,
                     topk=a.topk, cache=a.cache, verbose=not a.quiet)
    print()
    print(rep.summary())
    if rep.flagged:
        print("\nlineage-sweep flags (the label and the evidence disagree on lineage):")
        for r in rep.flagged:
            print("   %-30s %7d  %-26s -> %s"
                  % (r["atlas_label"][:30], r["n_cells"],
                     (r["assignment"]["label"] or "-")[:26],
                     r["audit"]["best_term_label"] or "-"))
    q = rep.queue[:a.queue]
    if q:
        print("\nreview queue, most suspicious first (top %d of %d):" % (len(q), len(rep.queue)))
        for i, r in enumerate(q, 1):
            print("   %2d. %-30s %7d  %-26s -> %s"
                  % (i, r["atlas_label"][:30], r["n_cells"],
                     (r["assignment"]["label"] or "-")[:26],
                     r["audit"]["best_term_label"] or "-"))
    if a.out:
        rep.to_json(a.out)
        print("\nwrote %s" % a.out)
    if a.tsv:
        rep.to_tsv(a.tsv)
        print("wrote %s" % a.tsv)
    return 0


def _annotate(a):
    from .annotate import annotate_h5ad, ACCURACY_NOTE
    doc = annotate_h5ad(a.h5ad, tissue=a.tissue, cluster_key=a.cluster_key,
                        min_cells=a.min_cells, topn=a.topn, cache=a.cache,
                        verbose=not a.quiet)
    cl = doc["clusters"]
    print("\n%d clusters (>=%d cells), grouped on obs/%s\n"
          % (len(cl), a.min_cells, doc["meta"]["cluster_key"]))
    print("  !! PROPOSALS, NOT CALLS: %s\n" % ACCURACY_NOTE)
    for c in cl[:a.show]:
        head = c["proposals"][0]["label"] if c["proposals"] else "(%s)" % c["status"]
        agree = "" if c["lineage_agreement"] is not False else "   [top-%d span >1 lineage]" % a.topn
        print("  %-24s %7d  %s%s" % (str(c["cluster"])[:24], c["n_cells"], head[:40], agree))
        for p in c["proposals"][1:]:
            print("  %-24s %7s    %d. %-38s %.2f" % ("", "", p["rank"], p["label"][:38],
                                                     p["relative"]))
        print("      markers: %s" % ", ".join(c["markers"][:6]))
    if a.out:
        json.dump(doc, open(a.out, "w"), indent=1)
        print("\nwrote %s" % a.out)
    return 0


def _resolve(a):
    from .ontology import Ontology
    from .resolve import Resolver
    o = Ontology.load(cache=a.cache)
    r = Resolver(o)
    for label in a.labels:
        cur, how = r.resolve(label, organ=a.organ)
        print("%-34s %-16s %s" % (label, how, o.label(cur) if cur else "-"))
    return 0


def _reference(a):
    from .reference import Reference
    genes = sorted(set(open(a.genes).read().split())) if a.genes else []
    ref = Reference.fetch(a.tissues.split(","), genes, cache=a.cache)
    ref.save(a.out)
    print("wrote %s.npz and %s_index.json (%d tissues)" % (a.out, a.out, len(ref.mats)))
    return 0


def main(argv=None):
    # prog matches the console script, not the module: the help text is what a user reads
    # after `pip install`, and it should show the command they actually have.
    p = argparse.ArgumentParser(prog="celltype-audit", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version",
                   version="celltype-audit %s" % __version__)
    p.add_argument("--cache", help="cache directory (default ~/.cache/celltype-audit)")
    sub = p.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("audit", help="audit an h5ad's cell-type labels")
    a.add_argument("h5ad")
    a.add_argument("--organ", help="organ name, used to resolve organ-qualified labels")
    a.add_argument("--type-key", default="cell_type",
                   help="obs column holding the cell-type labels (default: cell_type)")
    a.add_argument("--tissue", help="UBERON id, if the file does not carry one")
    a.add_argument("--min-cells", type=int, default=500)
    a.add_argument("--topk", type=int, default=5)
    a.add_argument("--queue", type=int, default=15, help="how many queue rows to print")
    a.add_argument("-o", "--out", help="write the full report as JSON")
    a.add_argument("--tsv", help="write a flat TSV")
    a.add_argument("-q", "--quiet", action="store_true")
    a.set_defaults(fn=_audit)

    n = sub.add_parser("annotate",
                       help="propose CL shortlists for UNANNOTATED clusters (curation aid)")
    n.add_argument("h5ad")
    n.add_argument("--tissue", help="UBERON id, if the file does not carry one")
    n.add_argument("--cluster-key", help="obs column holding the clustering")
    n.add_argument("--min-cells", type=int, default=500)
    n.add_argument("--topn", type=int, default=5)
    n.add_argument("--show", type=int, default=12, help="how many clusters to print")
    n.add_argument("-o", "--out")
    n.add_argument("-q", "--quiet", action="store_true")
    n.set_defaults(fn=_annotate)

    r = sub.add_parser("resolve", help="map labels onto CL terms")
    r.add_argument("labels", nargs="+")
    r.add_argument("--organ")
    r.set_defaults(fn=_resolve)

    f = sub.add_parser("reference", help="fetch and cache an expression reference")
    f.add_argument("--tissues", required=True, help="comma-separated UBERON ids")
    f.add_argument("--genes", required=True, help="file of gene symbols, one per line")
    f.add_argument("-o", "--out", default="reference")
    f.set_defaults(fn=_reference)

    args = p.parse_args(argv)
    try:
        return args.fn(args)
    except KeyboardInterrupt:
        return 130
    except Exception as ex:
        print("error: %s" % ex, file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
