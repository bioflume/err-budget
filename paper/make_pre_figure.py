"""Figure of REPORT.md: PRE's 95% band is too narrow on 2D-1, and its selector misses or declines.

Copy of studies/pre_comparison/make_figure.py for the article's supplement (Fig. S-4.1), with the
comparisons renamed: "Arm 1" is the same-runs comparison (PRE fitted to the runs the program
used) and "Arm 2" the matched-budget comparison (PRE chooses its own runs at the program's
search cost).  It reads studies/pre_comparison/results/ and writes figures/fig_sup_pre.pdf; nothing in
studies/pre_comparison/ is changed.

(a) Share of PRE fits whose band m +- 2s contains the exact value, per test case: same-runs
    comparison (608 fits on errbudget's runs) and matched-budget comparison (50 designs at
    budget B).  Nominal: 95%.
(b) Outcome of the 608 searches of the same-runs comparison per test case at z = 2.

    PYTHONDONTWRITEBYTECODE=1 python make_pre_figure.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "studies", "pre_comparison"))

import numpy as np
import pandas as pd

from orx_figstyle import PALETTE, figure_grid, panel_labels, save, use_style
import common as K

MM = 1 / 25.4
LAB = {"cd": "$c_D$", "dp": r"$\Delta p$", "cl": "$c_L$"}
TG = {1e-2: r"$10^{-2}$", 3e-3: r"$3{\times}10^{-3}$", 1e-3: r"$10^{-3}$", 3e-2: r"$3{\times}10^{-2}$",
      1e-1: r"$10^{-1}$"}


def coverage(g):
    eps = np.array([K.eps_abs(o, t) for o, t in zip(g.ob, g.target)])
    Q = np.array([K.CASE.exact[o] for o in g.ob])
    return 100 * np.mean(np.abs(g.pre_m - Q) <= 2 * g.pre_s)


def main():
    a1 = pd.read_csv(os.path.join(K.RESULTS, "arm1_rows.csv"))
    a2 = pd.read_csv(os.path.join(K.RESULTS, "arm2_rows.csv"))
    a2 = a2[a2.budget_factor == 1]
    tests = K.TESTS[::-1]                                   # first test case at the top
    y = np.arange(len(tests))
    labels = [f"{LAB[o]} {TG[t]}" for o, t in tests]

    use_style()
    fig, (ax1, ax2) = figure_grid(1, 2, width=89 * MM, ratio=0.82, sharey=True)

    c1 = [coverage(a1[(a1.ob == o) & (a1.target == t)]) for o, t in tests]
    c2 = [coverage(a2[(a2.ob == o) & (a2.target == t)]) for o, t in tests]
    ax1.axvline(95, color="0.35", lw=0.7, ls="--", zorder=1)
    ax1.text(93, len(tests) - 0.35, "nominal", ha="right", va="center", fontsize=6, color="0.35")
    ax1.scatter(c1, y + 0.12, s=14, color=PALETTE["orange"], zorder=3, label="same runs")
    ax1.scatter(c2, y - 0.12, s=14, facecolor="white", edgecolor=PALETTE["orange"], lw=0.9, marker="s",
                zorder=3, label="matched budget")
    ax1.set_xlim(0, 104)
    ax1.set_xticks([0, 25, 50, 75, 100])
    ax1.set_xlabel(r"PRE fits with $Q^*$ in $m\pm2s$ (%)")
    ax1.set_yticks(y, labels)
    ax1.set_ylabel("Observable, target")
    ax1.set_ylim(-0.6, len(tests) - 0.4 + 0.3)
    ax1.grid(axis="x", color="#D9D9D9", lw=0.4)
    ax1.grid(axis="y", visible=False)

    h = 0.34
    for k, (o, t) in enumerate(tests):
        g = a1[(a1.ob == o) & (a1.target == t)]
        n = len(g)
        met_eb = 100 * np.mean(g.eb_err <= 1)
        ax2.barh(k + h / 2 + 0.02, met_eb, height=h, color=PALETTE["blue"], lw=0)
        miss = 100 * np.mean(g.pre_err_z2 > 1)
        dec = 100 * np.mean(g.pre_err_z2.isna())
        met = 100 - miss - dec
        ax2.barh(k - h / 2 - 0.02, met, height=h, color=PALETTE["orange"], lw=0)
        ax2.barh(k - h / 2 - 0.02, miss, left=met, height=h, color=PALETTE["red"], lw=0)
        ax2.barh(k - h / 2 - 0.02, dec, left=met + miss, height=h, color="#BBBBBB", lw=0)
        if miss > 0:
            ax2.text(101, k - h / 2 - 0.02, f"{int((g.pre_err_z2 > 1).sum())}", va="center", ha="left",
                     fontsize=5.5, color=PALETTE["red"])
    ax2.set_xlim(0, 100)
    ax2.set_xlabel("Same-runs searches (%)")
    ax2.grid(axis="y", visible=False)
    ax2.grid(axis="x", color="#D9D9D9", lw=0.4)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [Line2D([], [], ls="none", marker="o", ms=4, color=PALETTE["orange"], label="(a) same runs"),
               Line2D([], [], ls="none", marker="s", ms=4, mfc="white", mec=PALETTE["orange"], label="(a) matched budget"),
               Patch(color=PALETTE["blue"], label="errbudget met"),
               Patch(color=PALETTE["orange"], label="PRE met"),
               Patch(color=PALETTE["red"], label="PRE missed"), Patch(color="#BBBBBB", label="PRE declined")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, fontsize=6, handlelength=1.0,
               columnspacing=0.8, handletextpad=0.3)
    panel_labels([ax1, ax2])
    save(fig, os.path.join(HERE, "figures", "fig_sup_pre"), formats=("pdf",))


if __name__ == "__main__":
    main()
