"""Coarse lineage assignment from a cell-type LABEL, and the crowding metric it supports.

Deliberately label-based and a priori: it needs no mapping, no reference and no gold standard,
so it can be reported for every organ. Coarse lineage is the one thing a label reliably tells
you ("CD8 T cell" is obviously haematopoietic even if which T cell is unclear) — which is
precisely the regime where string matching is trustworthy.

CROWDING = for each annotated cell type, how many OTHER types in the same organ share its
lineage. Liver showed that accuracy is governed by how finely a lineage is subdivided, so
crowding is the a-priori predictor of difficulty.
"""
import re

LINEAGE = [
    ("haematopoietic", r"\bt[ -]cell|\bb[ -]cell|\bnk\b|natural killer|lymph|myeloid|"
                       r"macrophage|monocyte|dendritic|\bdc\b|mast|granulocyte|neutrophil|"
                       r"eosinophil|basophil|plasma|kupffer|microglia|megakaryocyt|erythro|"
                       r"platelet|haematopoiet|hematopoiet|leukocyte|treg|mait|\bilc|thymocyte|"
                       r"langerhans|osteoclast|mononuclear|gdt|prebnk|\bcd4|\bcd8"),
    ("neural",        r"neuron|neural|glia|astrocyte|oligodendro|schwann|purkinje|photoreceptor|"
                      r"bipolar cell|amacrine|horizontal cell|satellite|nerval|ganglion"),
    ("endothelial",   r"endotheli|aerocyte"),
    ("muscle",        r"muscle|myocyte|cardiomyocyt|pericyte|myofibro"),
    ("stromal",       r"fibroblast|stromal|stellate|mesenchym|adipocyte|chondrocyte|osteo|"
                      r"interstitial|capsular|mesothel"),
    ("epithelial",    r"epitheli|keratinocyte|hepatocyte|cholangiocyte|acinar|ductal|alveolar|"
                      r"goblet|enterocyte|secretory|basal cell|club cell|ciliated|urothel|"
                      r"podocyte|tubule|principal cell|intercalated|paneth|tuft|ionocyte|"
                      r"endocrine|exocrine|alpha cell|beta cell|delta cell|gamma cell|islet|"
                      r"chromaffin|steroidogenic|pneumocyte|hepatoblast|enteroendocrine"),
    ("germ",          r"germ|sperm|oocyte|spermato|sertoli|leydig"),
]


def lineage_of(label):
    L = label.lower()
    for name, pat in LINEAGE:
        if re.search(pat, L):
            return name
    return "other"


def crowding(labels):
    """(mean co-annotated siblings per type, largest lineage size, per-lineage counts)."""
    counts = {}
    for L in labels:
        k = lineage_of(L)
        counts[k] = counts.get(k, 0) + 1
    if not labels:
        return 0.0, 0, counts
    mean = sum((counts[lineage_of(L)] - 1) for L in labels) / float(len(labels))
    return mean, (max(counts.values()) if counts else 0), counts
