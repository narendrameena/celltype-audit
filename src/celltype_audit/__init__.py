"""celltype-audit -- audit single-cell atlas annotations against the Cell Ontology.

The premise is that a cell type's LABEL and its EXPRESSION are two independent claims
about the same cells, and that disagreement between them is worth surfacing. The package
does not try to re-annotate an atlas; it reports where an atlas disagrees with itself.

    from celltype_audit import Ontology, Resolver, audit_h5ad
    report = audit_h5ad("lung.h5ad", organ="Lung")
    print(report.summary())

For data that has no annotation yet, `annotate_h5ad` proposes a ranked CL shortlist per
cluster -- a curation aid, not an annotator: top-1 is right about 56%% of the time and its
confidence is not calibratable, so read the shortlist rather than the top hit.

or from the command line:

    celltype-audit audit lung.h5ad --organ Lung -o annotations.json
    celltype-audit annotate new.h5ad --tissue UBERON:0002048
"""
from .ontology import Ontology
from .resolve import Resolver
from .annotate import annotate_h5ad
from .audit import audit_h5ad, Report

__version__ = "0.1.1"
__all__ = ["Ontology", "Resolver", "audit_h5ad", "annotate_h5ad", "Report", "__version__"]
