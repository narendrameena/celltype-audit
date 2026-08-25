#!/usr/bin/env python3
"""A graph per proposal, drawn from the real Cell Ontology.

Every card linked to the same general explainer, which tells a curator what an anchor set
is but nothing about the proposal in front of them. What they need is the specific
structure the proposal concerns: where the name lands in CL, where the expression lands,
which anchors each reaches, and -- where exclusivity failed -- the term that beat it.

So each proposal gets its own diagram, built from the pinned release rather than drawn:
the is_a path from each term up to its anchor, both columns side by side, and the point
where they meet or fail to. Nothing here is illustrative; every node and edge is read out
of cl-full.json.

Emitted as inline SVG into proposals.json so the page needs no diagramming library and no
second request. Colours are CSS custom properties, so the diagrams follow the reader's
theme with the rest of the page.

Usage:
    python benchmark/proposal_graphs.py
"""
import html
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
DOCS = os.path.abspath(os.path.join(HERE, "..", "..", "celltype-audit", "docs"))

from cl_lineage import load, ancestors, anchor_set, ANCHORS                # noqa: E402

MAXUP = 3          # is_a steps drawn before eliding; deep chains are the norm in CL
W, BW, BH, VG = 520, 210, 40, 56


def _chain(cur, g, stop):
    """The is_a path upward, first parent at each step, halted at an anchor."""
    out, seen = [cur], {cur}
    while len(out) <= MAXUP:
        ps = [p for p in g["parents"].get(out[-1], ()) if p.startswith("CL:")]
        ps = [p for p in ps if p not in seen]
        if not ps:
            break
        out.append(ps[0]); seen.add(ps[0])
        if ps[0] in stop:
            break
    return out


def _box(x, y, lines, kind="n"):
    fill = {"n": "var(--d-fill)", "ok": "var(--d-ok)", "bad": "var(--d-bad)",
            "anc": "var(--d-anchor)"}[kind]
    stroke = {"n": "var(--d-stroke)", "ok": "var(--d-okline)", "bad": "var(--d-badline)",
              "anc": "var(--d-okline)"}[kind]
    o = ['<rect x="%d" y="%d" width="%d" height="%d" rx="4" fill="%s" stroke="%s" '
         'stroke-width="1.1"/>' % (x, y, BW, BH, fill, stroke)]
    ty = y + BH / 2 - (len(lines) - 1) * 7 + 4
    for i, (t, small) in enumerate(lines):
        o.append('<text x="%d" y="%.1f" text-anchor="middle" font-size="%d" '
                 'fill="var(--d-ink)" font-weight="%s">%s</text>'
                 % (x + BW / 2, ty + i * 14, 9 if small else 11,
                    "400" if small else "600", html.escape(t[:34])))
    return "".join(o)


def _arrow(x, y1, y2, label=None):
    o = ['<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="var(--d-line)" stroke-width="1.1" '
         'marker-end="url(#pa)"/>' % (x, y1, x, y2)]
    if label:
        o.append('<text x="%d" y="%.1f" text-anchor="middle" font-size="8" '
                 'fill="var(--d-mute)">%s</text>' % (x + 16, (y1 + y2) / 2 + 3, label))
    return "".join(o)


def graph(p, g):
    lab = g["label"]
    anchors = set(ANCHORS)
    left = (p.get("candidate") or {}).get("curie")
    right = ((p.get("expression_top") or [{}])[0] or {}).get("curie")
    rival = p.get("weakened_by") or p.get("covered_by")
    rival_c = next((c for c, n in lab.items() if n == rival), None) if rival else None
    if not left and not right:
        return None

    cols = []
    if left:
        cols.append(("the name resolves to", left, "bad" if p.get("readiness") != "ready" else "n"))
    if right and right != left:
        cols.append(("the markers point to", right, "ok"))
    if rival_c and rival_c not in (left, right):
        cols.append(("exclusivity lost to", rival_c, "bad"))
    if not cols:
        return None

    n = len(cols)
    width = max(W, n * (BW + 30))
    chains = [(_chain(c, g, anchors), head, kind) for head, c, kind in cols]
    depth = max(len(ch) for ch, _, _ in chains)
    height = 26 + depth * VG + 10
    body = []
    for i, (ch, head, kind) in enumerate(chains):
        x = 15 + i * (width - 30 - BW) // max(n - 1, 1) if n > 1 else (width - BW) // 2
        body.append('<text x="%d" y="14" text-anchor="middle" font-size="8.5" '
                    'fill="var(--d-mute)" letter-spacing=".06em">%s</text>'
                    % (x + BW / 2, html.escape(head.upper())))
        for j, c in enumerate(ch):
            y = 26 + j * VG
            k = "anc" if c in anchors else (kind if j == 0 else "n")
            body.append(_box(x, y, [(lab.get(c, c), False), (c, True)], k))
            if j:
                body.append(_arrow(x + BW / 2, y + BH + VG - BH - 2, y - 2, "is_a"))
    return ('<svg viewBox="0 0 %d %d" width="100%%" style="max-width:%dpx" '
            'xmlns="http://www.w3.org/2000/svg" role="img" aria-label="%s in the Cell Ontology" '
            'font-family="ui-monospace,SFMono-Regular,Menlo,monospace">'
            '<defs><marker id="pa" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5" '
            'markerHeight="5" orient="auto-start-reverse">'
            '<path d="M0,0 L10,5 L0,10 z" fill="var(--d-line)"/></marker></defs>%s</svg>'
            % (width, height, width, html.escape(p["label"]), "".join(body)))


def main():
    g = load()
    doc = json.load(open(os.path.join(DOCS, "proposals.json")))
    made = 0
    for p in doc["proposals"]:
        s = graph(p, g)
        if s:
            p["graph"] = s
            made += 1
        print("  %-9s %-28s %s" % (p["organ"], p["label"][:28], "graph" if s else "-"))
    json.dump(doc, open(os.path.join(DOCS, "proposals.json"), "w"), indent=1)
    print("\n  %d of %d proposals carry an ontology graph" % (made, len(doc["proposals"])))


if __name__ == "__main__":
    main()
