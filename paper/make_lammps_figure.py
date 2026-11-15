"""Figure: the method on LAMMPS, against LAMMPS's own documented practice.

Regenerated for the revision without version history: the marker of the setting that an
earlier version returned ("this method") is removed.  The released program (errbudget.V121)
is run from the documented setting only to report its outcome, which the caption states.
The figure is written to ERRBUDGET_FIG_DIR (default: figures/)."""
import itertools
import os
import pickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from errbudget import run_protocol, V121
from lammps_madelung import LammpsMadelung

OUT = os.environ.get("ERRBUDGET_FIG_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures") + "/")
C = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7"]
INK, MUTED = "#1a1a1a", "#6b6b6b"
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["STIXGeneral"], "mathtext.fontset": "stix", "font.size": 9,
    "axes.labelsize": 9, "axes.titlesize": 9.5, "legend.fontsize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.spines.top": False, "axes.spines.right": False,
    "grid.color": "#dddddd", "grid.linewidth": 0.4, "lines.linewidth": 1.6,
    "figure.dpi": 150, "savefig.bbox": "tight", "savefig.pad_inches": 0.02})

LAD = [[6.0, 8.0, 10.0, 12.0, 14.0], [1e-2, 1e-3, 1e-4, 1e-5, 1e-6], [4, 5, 6, 7]]
KINDS = ["discretisation", "discretisation", "cost-only"]
HS = [[1.0 / c for c in LAD[0]], list(LAD[1]), [1.0 / o for o in LAD[2]]]
ALL = list(itertools.product(*[range(len(l)) for l in LAD]))

p = LammpsMadelung(n_cells=20, n_eval=3)
QS = p.exact_observables()["energy_per_ion"]
from paths import DATA
tensor = pickle.load(open(os.path.join(DATA,
                                       "lammps_tensor_n20.pkl"), "rb"))
rel = lambda i: abs(tensor[i][0]["energy_per_ion"] - QS) / abs(QS)
cost = lambda i: tensor[i][1]
look = lambda i: (tensor[tuple(i)][0]["energy_per_ion"], tensor[tuple(i)][1])

fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.9))

# ---- left: the whole parameter space, and four ways of choosing within it ----
target = 1e-5
eps = target * abs(QS)
feas = [i for i in ALL if rel(i) <= target]
best = min(feas, key=cost)
doc = (LAD[0].index(10.0), LAD[1].index(1e-5), LAD[2].index(5))
cautious = (LAD[0].index(14.0), LAD[1].index(1e-6), LAD[2].index(7))
r = run_protocol(look, LAD, doc, eps, hs=HS, kinds=KINDS, **V121)
print(f"released program from the documented setting: accepted={r['certified']}, "
      f"declined={r.get('declined')}, proposal={r['proposal']}, error={rel(r['proposal']):.2e}")

ax[0].scatter([cost(i) for i in ALL], [rel(i) for i in ALL], s=7,
              color="#c9c9c9", edgecolor="none", zorder=1,
              label="all 100 settings")
ax[0].axhline(target, color=INK, lw=1.0, ls="--", zorder=2)
ax[0].annotate("accuracy requested", xy=(0.02, target), xycoords=("axes fraction", "data"),
               xytext=(0, 3), textcoords="offset points", fontsize=7.5, color=INK)
for idx, col, mark, size, lab in [
        (doc,      C[1], "s", 62, "LAMMPS documented practice"),
        (cautious, C[4], "D", 46, "cautious user"),
        (best,     C[0], "o", 62, "cheapest that attains it")]:
    ax[0].scatter([cost(idx)], [rel(idx)], s=size, marker=mark, color=col,
                  edgecolor="white", linewidth=0.8, zorder=5, label=lab)
ax[0].set_xscale("log"); ax[0].set_yscale("log")
ax[0].set(xlabel="cost per evaluation (s)", ylabel="relative error in the energy")
ax[0].set_title("(a) All 100 settings", loc="left")
ax[0].grid(alpha=.5, lw=.4)
fig.legend(*ax[0].get_legend_handles_labels(), frameon=False, fontsize=7.5, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.16), handletextpad=0.3, columnspacing=1.2)

# ---- right: the accuracy-neutral parameter LAMMPS does not give advice on ----
orders = LAD[2]
for ci, (cut, acc) in enumerate([(10.0, 1e-5), (8.0, 1e-6)]):
    ii = [(LAD[0].index(cut), LAD[1].index(acc), j) for j in range(len(orders))]
    cs = [cost(i) for i in ii]
    ax[1].plot(orders, cs, "o-", color=C[ci], ms=4.5, markeredgecolor="white",
               markeredgewidth=0.5,
               label=f"cutoff {cut:.0f} \u00c5, accuracy $10^{{{int(np.log10(acc))}}}$")
    j5 = orders.index(5)
    ax[1].scatter([5], [cs[j5]], s=90, facecolor="none", edgecolor=C[ci],
                  linewidth=1.4, zorder=4)
ax[1].annotate("LAMMPS default", xy=(5, cs[orders.index(5)]), xytext=(10, 12),
               textcoords="offset points", ha="left", fontsize=7.5, color=MUTED,
               arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7))
ax[1].set_yscale("log")
ax[1].set(xlabel="PPPM interpolation order", ylabel="cost per evaluation (s)")
ax[1].set_xticks(orders)
ax[1].set_title("(b) Cost against interpolation order", loc="left")
ax[1].grid(alpha=.5, lw=.4)
ax[1].legend(frameon=False, fontsize=7.5, loc="upper right")

fig.subplots_adjust(wspace=0.34)
fig.savefig(OUT + "fig_lammps.pdf")
print(f"lammps figure done: documented={cost(doc):.4f}s err={rel(doc):.2e}, "
      f"cautious={cost(cautious):.4f}s err={rel(cautious):.2e}, best={cost(best):.4f}s err={rel(best):.2e}")
