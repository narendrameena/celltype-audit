#!/usr/bin/env python3
"""Colour-blind safety validator for a categorical palette (Python port; no node available).

Checks, per the dataviz method:
  1 lightness band      - all swatches in a usable L* range for the surface
  2 chroma floor        - no swatch so grey it reads as "disabled"
  3 CVD separation      - adjacent pairs stay apart under deuter/protan/tritan (OKLab dE >= 8)
  4 normal-vision floor - adjacent pairs apart for full-colour readers (dE >= 15)
  5 contrast            - each swatch vs the chart surface (>= 3:1 for marks)

Usage: python validate_palette.py "#hex,#hex,..." [--surface "#ffffff"] [--pairs all]
"""
import itertools
import sys


def hex2rgb(h):
    h = h.strip().lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


def _lin(c):
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def rgb2oklab(rgb):
    r, g, b = (_lin(c) for c in rgb)
    l = 0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b
    m = 0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b
    s = 0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b
    l_, m_, s_ = l ** (1 / 3), m ** (1 / 3), s ** (1 / 3)
    return (0.2104542553 * l_ + 0.7936177850 * m_ - 0.0040720468 * s_,
            1.9779984951 * l_ - 2.4285922050 * m_ + 0.4505937099 * s_,
            0.0259040371 * l_ + 0.7827717662 * m_ - 0.8086757660 * s_)


def dE(a, b):
    """OKLab distance x100 (the scale the dataviz thresholds are quoted in)."""
    A, B = rgb2oklab(a), rgb2oklab(b)
    return 100.0 * sum((x - y) ** 2 for x, y in zip(A, B)) ** 0.5


# Viénot/Brettel-style LMS simulation matrices
_SIM = {
    "deuteranopia": ((1.0, 0.0, 0.0), (0.9513092, 0.0, 0.04866992), (0.0, 0.0, 1.0)),
    "protanopia":   ((0.0, 1.05118294, -0.05116099), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
    "tritanopia":   ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (-0.86744736, 1.86727089, 0.0)),
}
_RGB2LMS = ((0.31399022, 0.63951294, 0.04649755),
            (0.15537241, 0.75789446, 0.08670142),
            (0.01775239, 0.10944209, 0.87256922))
_LMS2RGB = ((5.47221206, -4.6419601, 0.16963708),
            (-1.1252419, 2.29317094, -0.1678952),
            (0.02980165, -0.19318073, 1.16364789))


def _mv(M, v):
    return tuple(sum(M[i][j] * v[j] for j in range(3)) for i in range(3))


def simulate(rgb, kind):
    lms = _mv(_RGB2LMS, [_lin(c) for c in rgb])
    out = _mv(_LMS2RGB, _mv(_SIM[kind], lms))
    def g(c):
        c = max(0.0, min(1.0, c))
        return 12.92 * c if c <= 0.0031308 else 1.055 * c ** (1 / 2.4) - 0.055
    return tuple(g(c) for c in out)


def relL(rgb):
    r, g, b = (_lin(c) for c in rgb)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def contrast(a, b):
    la, lb = relL(a), relL(b)
    hi, lo = max(la, lb), min(la, lb)
    return (hi + 0.05) / (lo + 0.05)


def validate(hexes, surface="#ffffff", all_pairs=False):
    rgbs = [hex2rgb(h) for h in hexes]
    surf = hex2rgb(surface)
    rows, fails, warns = [], 0, 0
    print("palette: %s   surface %s\n" % (", ".join(hexes), surface))
    print("%-4s %-9s %6s %6s %7s  %s" % ("#", "hex", "L*", "chroma", "contr", "status"))
    for i, (h, c) in enumerate(zip(hexes, rgbs)):
        L, a_, b_ = rgb2oklab(c)
        chroma = 100 * (a_ ** 2 + b_ ** 2) ** 0.5
        cr = contrast(c, surf)
        st = []
        if not (0.35 <= L <= 0.75):
            st.append("L-band FAIL")
            fails += 1
        if chroma < 5:
            st.append("chroma FAIL")
            fails += 1
        if cr < 3.0:
            st.append("contrast WARN")
            warns += 1
        print("%-4d %-9s %6.3f %6.1f %7.2f  %s" % (i + 1, h, L, chroma, cr, ", ".join(st) or "ok"))
    pairs = list(itertools.combinations(range(len(rgbs)), 2)) if all_pairs \
        else [(i, i + 1) for i in range(len(rgbs) - 1)]
    print("\n%-9s %8s %10s %10s %10s  %s" % ("pair", "normal", "deuter", "protan", "tritan", "status"))
    for i, j in pairs:
        n = dE(rgbs[i], rgbs[j])
        d = {k: dE(simulate(rgbs[i], k), simulate(rgbs[j], k)) for k in _SIM}
        st = []
        if n < 15:
            st.append("NORMAL FAIL")
            fails += 1
        worst = min(d.values())
        if worst < 6:
            st.append("CVD FAIL")
            fails += 1
        elif worst < 8:
            st.append("CVD floor (needs 2nd encoding)")
            warns += 1
        print("%-9s %8.1f %10.1f %10.1f %10.1f  %s"
              % ("%d-%d" % (i + 1, j + 1), n, d["deuteranopia"], d["protanopia"], d["tritanopia"],
                 ", ".join(st) or "ok"))
    print("\n%d FAIL, %d WARN" % (fails, warns))
    return fails


if __name__ == "__main__":
    hexes = [x for x in sys.argv[1].split(",") if x.strip()]
    surface = sys.argv[sys.argv.index("--surface") + 1] if "--surface" in sys.argv else "#ffffff"
    sys.exit(1 if validate(hexes, surface, "--pairs" in sys.argv) else 0)
