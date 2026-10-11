"""Figures of revision 3:
  fig_dfg_mesh.pdf    OpenFOAM 2D-1: error against mesh level at the tightest and at the
                      tutorial solver settings (drag, pressure difference, lift), with
                      the version-1.1 post-processing for the drag
  fig_search.pdf      one search: the error budget at each step, against the true error
  fig_l5_stall.pdf    L5: error of the y-coordinate of the minimum against N
"""
import os
import pickle
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cases as C
import grids
from knobs import LocalBudget, run_protocol, V12

OUT = os.environ.get("KNOBS_FIG_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures"))
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True,
                     "grid.color": GRID, "grid.linewidth": 0.5, "lines.linewidth": 1.5,
                     "legend.frameon": False, "font.family": "serif"})


def fig_dfg():
    c = C._dfg1(cost="deterministic")
    T, LAD = c.tensor, c.LAD
    na = LAD[0]
    tight = lambda a: (a, len(LAD[1]) - 1, len(LAD[2]) - 1, 0)
    tut = lambda a: (a, 0, LAD[2].index(1e-6), 0)
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.6))
    nf = []
    for ax, ob, lab in zip(axs, ["cd", "dp", "cl"], ["drag $c_D$", "pressure difference $\\Delta p$",
                                                     "lift $c_L$"]):
        ex = c.exact[ob]
        e1 = [abs(T[tight(a)][0][ob] - ex) / abs(ex) for a in range(len(na))]
        e2 = [abs(T[tut(a)][0][ob] - ex) / abs(ex) for a in range(len(na))]
        ax.loglog(na, e1, "-o", color=BLUE, ms=3.5, label="tightest solver settings")
        ax.loglog(na, e2, "-s", color=ORANGE, ms=3.5, label="tutorial solver settings")
        if ob == "cd" and "cd_vm" in T[tight(0)][0]:
            e3 = [abs(T[tight(a)][0]["cd_vm"] - ex) / abs(ex) for a in range(len(na))]
            ax.loglog(na, e3, ":", color=INK2, lw=1.0, label="tightest, v1.1 post-processing")
        fine = os.path.join(C.D, "dfg_fine.json")
        if os.path.exists(fine):
            import json
            fj = json.load(open(fine))
            nf = [int(n) for n in sorted(fj, key=int) if int(n) > max(na)]
            if nf:
                ef = [abs(fj[str(n)]["obs"][ob] - ex) / abs(ex) for n in nf]
                ax.loglog([na[-1]] + nf, [e1[-1]] + ef, "--o", color=BLUE, ms=3.5, mfc="white",
                          lw=1.0, label="tightest, beyond the study grid")
        ax.set_title(lab, fontsize=8, color=INK)
        ax.set_xlabel("cells per quarter cylinder, $n_\\mathrm{arc}$")
        allx = sorted(set(na) | set(nf if os.path.exists(fine) and nf else []))
        ax.set_xticks(allx); ax.set_xticklabels([str(n) for n in allx], fontsize=6)
        ax.minorticks_off()
    axs[0].set_ylabel("relative error")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, fontsize=6.5, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(os.path.join(OUT, "fig_dfg_mesh.pdf"))
    plt.close(fig)


def fig_search():
    """The search from the tutorial settings at n_arc = 32 for the lift at 3e-2."""
    c = C._dfg1(cost="deterministic")
    T, LAD, ob = c.tensor, c.LAD, "cl"
    eps = 3e-2 * abs(c.exact[ob])
    look = lambda i: ((T[tuple(i)][0][ob], T[tuple(i)][1]) if c.converged(tuple(i))
                      else (np.nan, T[tuple(i)][1]))
    start = (LAD[0].index(32), 0, LAD[2].index(1e-6), 0)
    r = run_protocol(look, LAD, start, eps, hs=c.HS, kinds=c.KINDS, **V12)
    steps = []
    for a in r["trace"]:
        lb = LocalBudget(look, LAD, a, hs=c.HS, kinds=c.KINDS, local_order=True,
                         resid_margin=5.0, resolution=0.1 * eps, cancel_rule="none")
        lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
        if np.isfinite(lb.E0) and abs(lb.E0) <= eps:      # as the program does
            lb.joint_correct(look, eps, cancel_rule="none")
        steps.append((a, [abs(x) / eps if np.isfinite(x) else np.nan for x in lb.rk],
                      abs(lb.E0) / eps, abs(T[a][0][ob] - c.exact[ob]) / eps))
    names = ["mesh $n_\\mathrm{arc}$", "residualControl", "GAMG tolerance"]
    cols = [BLUE, ORANGE, AQUA]
    marks = ["o", "s", "^"]
    fig, ax = plt.subplots(figsize=(3.4, 2.7))
    x = np.arange(len(steps))
    TOP, BOT = 3e3, 1e-3
    clip = lambda v: np.clip(v, BOT, TOP)
    for k in range(3):
        v = np.array([s[1][k] for s in steps], float)
        fin = np.isfinite(v) & (v > 0)
        ax.plot(x[fin], clip(v[fin]), marks[k], color=cols[k], ms=4, mfc="white", mew=1.2,
                label=f"budget: {names[k]}")
        inf = ~np.isfinite(v)
        if inf.any():
            ax.plot(x[inf], np.full(inf.sum(), TOP * 0.6), marker=r"$\uparrow$", ls="none",
                    color=cols[k], ms=8)
    est = np.array([s[2] for s in steps], float)
    ax.plot(x[np.isfinite(est)], clip(est[np.isfinite(est)]), "D", color=INK, ms=4.5,
            label="estimate $\\widehat{E}$ (with margin)")
    if (~np.isfinite(est)).any():
        ax.plot(x[~np.isfinite(est)], np.full((~np.isfinite(est)).sum(), TOP * 0.6), marker=r"$\uparrow$",
                ls="none", color=INK, ms=9)
    ax.plot(x, clip(np.array([s[3] for s in steps])), "x", color=INK, ms=6, mew=1.5,
            label="true error")
    ax.axhline(1.0, color=INK2, lw=0.8, ls="--")
    ax.text(len(steps) - 0.6, 1.15, "target", fontsize=6, color=INK2, ha="right", va="bottom")
    ax.set_yscale("log"); ax.set_ylim(BOT, TOP)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{LAD[0][a[0]]}\n{LAD[1][a[1]]:.0e}\n{LAD[2][a[2]]:.0e}" for a, *_ in steps],
                       fontsize=5)
    ax.set_xlabel("setting at each step ($n_\\mathrm{arc}$, residualControl, GAMG tolerance)", fontsize=6.5)
    ax.set_ylabel("error / target")
    ax.legend(fontsize=5.5, loc="upper center", bbox_to_anchor=(0.5, 1.32), ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_search.pdf"))
    plt.close(fig)
    return steps, r


def fig_l5():
    P = grids.L5
    t = grids.load(P)
    ex = grids.exact(P)["y_min"]
    N = P["lad"][0]
    fin = (len(P["lad"][1]) - 1, len(P["lad"][2]) - 1, len(P["lad"][3]) - 1)
    e = [abs(t[(a,) + fin][0]["y_min"] - ex) / abs(ex) for a in range(len(N))]
    fig, ax = plt.subplots(figsize=(3.2, 2.2))
    ax.loglog(N, e, "-o", color=BLUE, ms=3.5)
    for j in range(len(N) - 1):
        p = np.log(e[j] / e[j + 1]) / np.log(N[j + 1] / N[j])
        ax.annotate(f"{p:.1f}", (np.sqrt(N[j] * N[j + 1]), np.sqrt(e[j] * e[j + 1])),
                    fontsize=6, color=INK2, textcoords="offset points", xytext=(0, 4), ha="center")
    ax.axhline(3e-4, color=INK2, lw=0.8, ls="--")
    ax.text(N[0], 3.3e-4, "tightest target, $3\\times10^{-4}$", fontsize=6, color=INK2)
    ax.set_xticks(N); ax.set_xticklabels([str(n) for n in N]); ax.minorticks_off()
    ax.set_xlabel("grid points per side, $N$")
    ax.set_ylabel("relative error of $y_\\mathrm{min}$")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_l5_stall.pdf"))
    plt.close(fig)
    return e


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    print("L5 errors:", [f"{x:.2e}" for x in fig_l5()])
    fig_dfg()
    steps, r = fig_search()
    print("search:", r["certified"], r["proposal"], r["runs"], "runs")
    for a, rk, E, tr in steps:
        print("  ", a, [f"{x:.2f}" for x in rk], f"estimate {E:.2f}", f"true {tr:.2f}")
