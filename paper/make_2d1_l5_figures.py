"""Figs. 1 (fig_teaser.pdf) and 3 (fig_l5_stall.pdf) and Figs. S-5.1 (fig_sup_dfg.pdf) and
S-5.2 (fig_sup_search.pdf), from figures/data_2d1_l5.json (extract_2d1_l5.py); prints the
plotted values.

    python make_2d1_l5_figures.py       # writes figures/*.pdf and a .png of each
"""
import json
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from figstyle import (use_style, panel_label, target_line, COLUMN, TEXT, SUPP,
                      BLUE, VERM, BLACK, GREY, GREEN, ORANGE, PURPLE, SHADE)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
D = json.load(open(os.path.join(OUT, "data_2d1_l5.json")))


def sci(x):
    e = int(np.floor(np.log10(x)))
    return f"10^{{{e}}}" if abs(x / 10 ** e - 1) < 1e-9 else f"{x / 10 ** e:g}\\times10^{{{e}}}"


def narc_axis(ax, ticks):
    ax.set_xscale("log")
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(n) for n in ticks])
    ax.minorticks_off()
    ax.set_xlabel(r"cells per quarter cylinder, $n_\mathrm{arc}$")


def save(fig, name):
    fig.savefig(os.path.join(OUT, name))
    fig.savefig(os.path.join(OUT, name.replace(".pdf", ".png")), dpi=200)
    plt.close(fig)
    print("wrote", name)


def teaser():
    f, s = D["dfg"], D["search"]
    na = f["n_arc"]
    tgt = s["target_rel"]
    steps = s["steps"]
    fig = plt.figure(figsize=(TEXT, 2.45))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.55], wspace=0.22,
                          left=0.075, right=0.995, bottom=0.17, top=0.93)
    # (a) lift error against the mesh
    ax = fig.add_subplot(gs[0])
    ax.plot(na, f["tut_cl"], "-s", color=VERM, ms=3.5, zorder=3)
    ax.plot(na, f["tight_cl"], "-o", color=BLACK, ms=3.2, zorder=3)
    target_line(ax, tgt)
    fin = steps[-1]
    e_fin = fin["true"] * tgt
    ax.plot([fin["setting"][0]], [e_fin], "*", color=BLUE, ms=9, mec="white", mew=0.5, zorder=5)
    ax.annotate("errbudget's setting\n" r"($r_c=10^{-4}$, $t=10^{-7}$)", (fin["setting"][0], e_fin),
                xytext=(14, -15), textcoords="offset points", fontsize=7, color=BLUE,
                ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=BLUE, lw=0.6, shrinkA=0, shrinkB=4))
    ax.text(30, 0.42, "tutorial tolerances\n" r"($r_c=10^{-2}$, $t=10^{-6}$)", color=VERM,
            fontsize=7, ha="right", va="center")
    ax.text(9.5, 5.2e-3, "tightest tolerances", color=BLACK, fontsize=7, ha="left", va="top")
    ax.text(64, tgt * 1.12, r"target $3\times10^{-2}$", color=GREY, fontsize=7, ha="right",
            va="bottom")
    narc_axis(ax, na)
    ax.set_yscale("log")
    ax.set_ylim(2e-3, 1.5)
    ax.set_ylabel("relative error of the lift")
    panel_label(ax, "(a)", x=-0.17)
    # (b) the search: estimate and true error per step
    ax = fig.add_subplot(gs[1])
    x = np.arange(len(steps))
    p1 = s["phase1_end"]
    ax.axvspan(p1 + 0.5, len(steps) - 0.5, color=SHADE, lw=0, zorder=0)
    TOP = 2e3
    ax.text(p1 / 2, TOP * 0.95, "Phase 1: refine", fontsize=7, color=GREY, ha="center", va="top")
    ax.text((p1 + len(steps)) / 2, TOP * 0.95, "Phase 2: save cost", fontsize=7, color=GREY,
            ha="center", va="top")
    est = np.array([np.nan if v is None else v for v in (st["estimate"] for st in steps)])
    tru = np.array([st["true"] for st in steps])
    ok = np.isfinite(est)
    ax.plot(x, tru, "-x", color=BLACK, ms=5, mew=1.2, lw=0.7, zorder=3, label="true error")
    ax.plot(x, est, "-D", color=BLUE, ms=3.6, lw=0.7, zorder=4, label=r"estimate $\widehat{E}$")
    for i in np.where(~ok)[0]:      # the estimate cannot be bounded
        ax.annotate("", (i, TOP * 0.55), xytext=(i, 120), zorder=4,
                    arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1.0, mutation_scale=7))
    target_line(ax, 1.0)
    ax.text(len(steps) - 0.6, 1.12, "target", fontsize=7, color=GREY, ha="right", va="bottom")
    # short labels on the moves that change the error; GAMG and coarsest-level moves unlabeled
    lab = {2: (r"$n_\mathrm{arc}$ 32$\to$64", 1.6, 2.6, "center"),
           6: (r"$r_c$ $10^{-2}\to10^{-3}$", 5.55, 2.2, "right"),
           8: (r"$r_c$ $10^{-3}\to10^{-4}$", 7.85, 0.06, "right"),
           10: (r"$n_\mathrm{arc}$ 64$\to$32", 9.6, 0.56, "center")}
    for i, (t, xt, yt, ha) in lab.items():
        ax.text(xt, yt, t, fontsize=7, color=BLACK, ha=ha, va="center")
    ax.set_yscale("log")
    ax.set_ylim(1e-2, TOP)
    ax.set_xlim(-0.6, len(steps) - 0.4)
    ax.set_xticks(x)
    ax.set_xlabel("step of the search")
    ax.set_ylabel("error / target")
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.86), handlelength=1.6)
    panel_label(ax, "(b)", x=-0.09)
    save(fig, "fig_teaser.pdf")
    print("teaser (a) n_arc, tutorial, tightest:", list(zip(na, f["tut_cl"], f["tight_cl"])))
    print("teaser (a) errbudget setting:", fin["setting"], "relative lift error", e_fin)
    print("teaser (b) step, setting, estimate/target, true/target:")
    for i, st in enumerate(steps):
        print("  ", i, st["setting"], st["estimate"], st["true"])


def sup_dfg():
    f = D["dfg"]
    na, nf = f["n_arc"], f["n_arc_fine"]
    fig, axs = plt.subplots(1, 3, figsize=(SUPP, 2.5))
    fig.subplots_adjust(left=0.085, right=0.995, bottom=0.33, top=0.89, wspace=0.34)
    names = {"cd": r"drag $c_D$", "dp": r"pressure difference $\Delta p$", "cl": r"lift $c_L$"}
    for n, (ax, ob) in enumerate(zip(axs, ["cd", "dp", "cl"])):
        ax.plot(na, f["tight_" + ob], "-o", color=BLACK, ms=3.2, zorder=3)
        ax.plot(na, f["tut_" + ob], "-s", color=VERM, ms=3.2, zorder=3)
        if ob == "cd":
            ax.plot(na, f["tight_cd_vm"], ":", color=GREY, lw=1.1, zorder=2)
        ax.plot([na[-1]] + nf, [f["tight_" + ob][-1]] + f["fine_" + ob], "--o", color=BLACK,
                ms=3.2, mfc="white", lw=0.8, markevery=list(range(1, len(nf) + 1)), zorder=3)
        narc_axis(ax, na + nf)
        ax.set_xlabel("")
        ax.tick_params(axis="x", labelsize=6.5)
        ax.set_yscale("log")
        ax.set_ylabel(f"relative error of the {names[ob].split(' $')[0]}")
        panel_label(ax, f"({'abc'[n]})", x=-0.05, y=1.07)
    fig.text(0.54, 0.195, r"cells per quarter cylinder, $n_\mathrm{arc}$", ha="center", va="center")
    h = [Line2D([], [], color=BLACK, marker="o", ms=3.2),
         Line2D([], [], color=VERM, marker="s", ms=3.2),
         Line2D([], [], color=GREY, ls=":", lw=1.1),
         Line2D([], [], color=BLACK, ls="--", marker="o", mfc="white", ms=3.2, lw=0.8)]
    fig.legend(h, ["tightest solver settings", "tutorial solver settings",
                   "tightest, vertex-mean cell centers (drag)", "tightest, beyond the study"],
               loc="lower center", ncol=2, bbox_to_anchor=(0.54, -0.01), columnspacing=2.0)
    save(fig, "fig_sup_dfg.pdf")
    for ob in ["cd", "dp", "cl"]:
        print("sup_dfg", ob, "tightest", f["tight_" + ob], "tutorial", f["tut_" + ob],
              "fine", f["fine_" + ob])
    print("sup_dfg cd_vm", f["tight_cd_vm"])


def sup_search():
    s = D["search"]
    steps = s["steps"]
    p1 = s["phase1_end"]
    fig, ax = plt.subplots(figsize=(SUPP, 2.9))
    fig.subplots_adjust(left=0.115, right=0.995, bottom=0.30, top=0.83)
    x = np.arange(len(steps))
    TOP, BOT = 3e3, 1e-3
    ax.axvspan(p1 + 0.5, len(steps) - 0.5, color=SHADE, lw=0, zorder=0)
    ax.text(p1 / 2, TOP * 0.7, "Phase 1", fontsize=7, color=GREY, ha="center", va="top")
    ax.text((p1 + len(steps)) / 2, TOP * 0.7, "Phase 2", fontsize=7, color=GREY, ha="center",
            va="top")
    names = [r"budget: mesh $n_\mathrm{arc}$", "budget: residualControl", "budget: GAMG tolerance",
             "budget: coarsest level (cost-only)"]
    cols, marks = [GREEN, PURPLE, ORANGE, GREY], ["o", "s", "^", "v"]
    for k in range(4):
        v = np.array([np.inf if st["budget"][k] is None else st["budget"][k] for st in steps])
        fin = np.isfinite(v) & (v > 0)
        ax.plot(x[fin], np.clip(v[fin], BOT, TOP), marks[k], color=cols[k], ms=4, mfc="white",
                mew=1.0, label=names[k], zorder=3)
        for i in np.where(~np.isfinite(v))[0]:
            ax.annotate("", (i, TOP * 0.9), xytext=(i, TOP * 0.12),
                        arrowprops=dict(arrowstyle="-|>", color=cols[k], lw=1.0, mutation_scale=7))
    est = np.array([np.nan if st["estimate"] is None else st["estimate"] for st in steps])
    ok = np.isfinite(est)
    ax.plot(x[ok], est[ok], "D", color=BLUE, ms=4.2, label=r"estimate $\widehat{E}$", zorder=4)
    for i in np.where(~ok)[0]:
        ax.annotate("", (i, TOP * 0.9), xytext=(i, TOP * 0.12),
                    arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=1.0, mutation_scale=7))
    ax.plot(x, [st["true"] for st in steps], "x", color=BLACK, ms=5.5, mew=1.3, label="true error",
            zorder=5)
    target_line(ax, 1.0)
    ax.text(len(steps) - 0.55, 1.15, "target", fontsize=7, color=GREY, ha="right", va="bottom")
    ax.set_yscale("log"); ax.set_ylim(BOT, TOP); ax.set_xlim(-0.6, len(steps) - 0.4)
    ax.set_xticks(x)

    def e(v):
        return f"{v:.0e}".replace("e-0", "e-")
    ax.set_xticklabels([f"{i}\n\n{st['setting'][0]}\n{e(st['setting'][1])}\n{e(st['setting'][2])}"
                        f"\n{st['setting'][3]}" for i, st in enumerate(steps)], fontsize=6)
    ax.text(-0.02, -0.098, "step", transform=ax.transAxes, fontsize=6, ha="right", va="center")
    for j, t in enumerate([r"$n_\mathrm{arc}$", r"$r_c$", "GAMG", "coarsest"]):
        ax.text(-0.02, -0.205 - j * 0.0705, t, transform=ax.transAxes, fontsize=6, ha="right",
                va="center")
    ax.set_ylabel("error / target")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=7,
              columnspacing=1.2, handletextpad=0.3)
    save(fig, "fig_sup_search.pdf")


def l5():
    d = D["l5"]
    N, e = d["N"], d["err"]
    fig, ax = plt.subplots(figsize=(COLUMN, 2.3))
    fig.subplots_adjust(left=0.19, right=0.97, bottom=0.19, top=0.96)
    ax.plot(N, e, "-o", color=BLACK, ms=3.2, zorder=3)
    for j in range(len(N) - 1):
        p = np.log(e[j] / e[j + 1]) / np.log(N[j + 1] / N[j])
        ax.annotate(f"{p:.1f}", (np.sqrt(N[j] * N[j + 1]), np.sqrt(e[j] * e[j + 1])),
                    fontsize=7, color=GREY, textcoords="offset points", xytext=(3, 3),
                    ha="left", va="bottom")
    target_line(ax, d["target"])
    ax.text(N[0], d["target"] * 0.85, r"tightest target, $3\times10^{-4}$", fontsize=7,
            color=GREY, va="top")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks(N); ax.set_xticklabels([str(n) for n in N]); ax.minorticks_off()
    ax.set_xlabel(r"grid points per side, $N$")
    ax.set_ylabel(r"relative error of $y_\mathrm{min}$")
    save(fig, "fig_l5_stall.pdf")
    print("l5 N, error, observed order:")
    for j in range(len(N)):
        p = np.log(e[j - 1] / e[j]) / np.log(N[j] / N[j - 1]) if j else float("nan")
        print("  ", N[j], e[j], round(p, 2))


if __name__ == "__main__":
    use_style()
    teaser()
    sup_dfg()
    sup_search()
    l5()
