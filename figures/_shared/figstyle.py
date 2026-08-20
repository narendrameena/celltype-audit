"""Shared figure style — Nature-journal print specification.

Widths: 89 mm (single column), 120 mm (1.5), 183 mm (double).
Type: sans-serif, 5-7 pt. Rules 0.5 pt. No gridlines heavier than the data.

Palette: Okabe-Ito subset + a deliberate achromatic neutral. Validated with
`validate_palette.py --pairs all`: all 10 pairs clear the normal-vision floor (dE >= 15)
and the CVD floor (dE >= 8) under deuteranopia, protanopia and tritanopia. The grey
intentionally fails the chroma check — it is the neutral/reference, not a series hue.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MM = 1 / 25.4
W1, W15, W2 = 89 * MM, 120 * MM, 183 * MM          # column widths in inches

BLUE   = "#0072B2"   # series 1 / primary result
VERM   = "#D55E00"   # series 2 / comparison
GREEN  = "#009E73"   # series 3
PINK   = "#CC79A7"   # series 4
GREY   = "#4D4D4D"   # neutral: baseline, reference, "not applicable"
SERIES = [BLUE, VERM, GREEN, PINK]

INK    = "#1A1A1A"   # primary text
INK2   = "#595959"   # secondary text
RULE   = "#BFBFBF"   # axes / rules
FAINT  = "#E8E8E8"   # very light fill


def apply():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 6.5,
        "axes.labelsize": 6.5,
        "axes.titlesize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "axes.edgecolor": RULE,
        "axes.linewidth": 0.5,
        "axes.labelcolor": INK,
        "axes.titlecolor": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": INK2, "ytick.color": INK2,
        "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2, "ytick.major.size": 2,
        "text.color": INK,
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 400,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.01,
        "pdf.fonttype": 42,          # editable text in Illustrator
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def panel(ax, letter, dx=-0.16, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=8,
            fontweight="bold", va="top", ha="left", color=INK)


def save(fig, outdir, name, source_df=None):
    """Write svg + pdf + png and the source-data TSV.

    In the themed layout (`<NN_theme>/script/*.py`) figures go to `<NN_theme>/figure/`
    and data to `<NN_theme>/sourceData/`. Anywhere else, everything lands in `outdir`.
    """
    import os
    if os.path.basename(os.path.abspath(outdir)) == "script":
        theme = os.path.dirname(os.path.abspath(outdir))
        figdir, datadir = os.path.join(theme, "figure"), os.path.join(theme, "sourceData")
    else:
        figdir = datadir = outdir
    os.makedirs(figdir, exist_ok=True)
    os.makedirs(datadir, exist_ok=True)
    paths = []
    for ext in ("svg", "pdf", "png"):
        p = os.path.join(figdir, "%s.%s" % (name, ext))
        fig.savefig(p, format=ext)
        paths.append(p)
    if source_df is not None:
        p = os.path.join(datadir, "%s_source_data.tsv" % name)
        with open(p, "w") as fh:
            fh.write(source_df)
        paths.append(p)
    plt.close(fig)
    return paths
