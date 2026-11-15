"""House style of the errbudget figures (CPC, elsarticle 5p two-column).

Widths are the measured LaTeX widths: \\columnwidth = 252 pt, \\textwidth = 522 pt in the
article, \\textwidth = 472 pt in the supplement (11 pt, A4, 2.2 cm margins). Text is serif
(STIX, a Times-metric face), 8 pt, at most 9 pt. Colours: Okabe--Ito, with fixed roles.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLUMN = 252 / 72.27          # 3.49 in
TEXT = 522 / 72.27            # 7.22 in
SUPP = 472.3 / 72.27          # 6.54 in

# Okabe--Ito, fixed roles across all figures
BLUE = "#0072B2"              # errbudget (estimate, returned setting)
VERM = "#D55E00"              # tutorial settings and procedures of current practice
BLACK = "#000000"             # true or reference error
GREY = "#777777"              # target line, annotations
SKY, GREEN, ORANGE, PURPLE, YELLOW = "#56B4E9", "#009E73", "#E69F00", "#CC79A7", "#F0E442"
SHADE = "#EEEEEE"


def use_style():
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["STIXGeneral"], "mathtext.fontset": "stix",
        "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
        "ytick.major.width": 0.6, "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
        "xtick.major.size": 3, "ytick.major.size": 3, "lines.linewidth": 1.2,
        "lines.markersize": 4, "legend.frameon": False, "axes.spines.top": False,
        "axes.spines.right": False, "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.dpi": 300,
    })


def panel_label(ax, s, x=-0.02, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight="bold",
            ha="right", va="bottom")


def target_line(ax, y, **kw):
    ax.axhline(y, color=GREY, lw=0.8, ls=(0, (4, 2.5)), zorder=1, **kw)
