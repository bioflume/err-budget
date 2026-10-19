"""Publication figures. Palette: Okabe-Ito, ordered so the worst adjacent CVD
pair clears dE 11 (validated); ordered magnitudes use a single-hue ramp."""
import sys, os, itertools, pickle, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
from knobs import run_protocol, LocalBudget
from knobs.problems import AdvectionDiffusion
from knobs.problems import Burgers2D

OUT = os.environ.get("KNOBS_FIG_DIR", "figures/") 
C = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7"]
INK, MUTED = "#1a1a1a", "#6b6b6b"
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"], "font.size": 9,
    "axes.labelsize": 9, "axes.titlesize": 9.5, "legend.fontsize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.spines.top": False, "axes.spines.right": False,
    "grid.color": "#dddddd", "grid.linewidth": 0.4, "lines.linewidth": 1.6,
    "figure.dpi": 150, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})

# ------------------------------------------------------------------ Fig 1
def fig_masking():
    ad = AdvectionDiffusion(nu=0.02, warm_start=False)
    QS = ad.exact_Q()
    tols = [10.0 ** -k for k in range(1, 11)]
    Ns = [64, 128, 256, 512, 1024]
    err = np.zeros((len(Ns), len(tols))); cost = np.zeros_like(err)
    for i, N in enumerate(Ns):
        for j, t in enumerate(tols):
            Q, Cc, _ = ad.run(dict(N=N, steps=400, tol=t))
            err[i, j], cost[i, j] = abs(Q - QS), Cc
    ramp = plt.cm.Blues(np.linspace(0.42, 0.95, len(Ns)))

    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.7))
    for i, N in enumerate(Ns):
        floor = err[i, -1]
        j = int(np.argmax(err[i] < 2 * floor))
        ax[0].loglog(tols, err[i], "-", color=ramp[i], marker="o", ms=2.6,
                     markerfacecolor=ramp[i], markeredgecolor="white", markeredgewidth=0.4)
        ax[0].axhline(floor, color=ramp[i], ls=(0, (1, 2)), lw=0.7)
        ax[0].plot(tols[j], err[i, j], "*", color=ramp[i], ms=11,
                   markeredgecolor="white", markeredgewidth=0.5, zorder=5)
        ax[0].plot([], [], "-", color=ramp[i], label=f"$N={N}$")
        ax[1].loglog(tols, cost[i], "-", color=ramp[i], marker="o", ms=2.6,
                     markerfacecolor=ramp[i], markeredgecolor="white", markeredgewidth=0.4)
        ax[1].plot(tols[j], cost[i, j], "*", color=ramp[i], ms=11,
                   markeredgecolor="white", markeredgewidth=0.5, zorder=5)
    ax[0].set(xlabel="linear-solver tolerance", ylabel=r"error in the observable")
    ax[0].set_title("accuracy stops improving", loc="left")
    ax[1].set(xlabel="linear-solver tolerance", ylabel="cost")
    ax[1].set_title("cost does not stop rising", loc="left")
    for a in ax:
        a.invert_xaxis(); a.grid(alpha=.5, which="major", lw=.4)
    ax[0].legend(frameon=False, loc="lower left", fontsize=7.5, ncol=2,
                 handlelength=1.2, columnspacing=1.0)
    ax[1].text(0.97, 0.06, "$\\star$ = last useful setting", transform=ax[1].transAxes,
               ha="right", fontsize=7.5, color=MUTED)
    fig.subplots_adjust(wspace=0.32)
    fig.savefig(OUT + "fig1_masking.pdf"); plt.close(fig)
    print("fig1 done")


# ------------------------------------------------------------------ Fig 2
def fig_budget(tensor, LAD, HS, NAMES):
    fig, ax = plt.subplots(figsize=(6.6, 2.2))
    anchor = (2, 2, 1, 1, 1)
    for oi, ob in enumerate(["energy", "gradient"]):
        look = lambda i, ob=ob: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
        lb = LocalBudget(look, LAD, anchor, hs=HS)
        lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
        vals = [abs(r) if np.isfinite(r) else 0.0 for r in lb.rk]
        tot = sum(vals) or 1.0
        y = np.arange(len(NAMES)) + (0.16 if oi == 0 else -0.16)
        ax.scatter([max(v / tot, 1e-26) for v in vals], y, s=34, color=C[oi],
                   edgecolor="white", linewidth=0.6, zorder=3,
                   label={"energy": "energy", "gradient": "gradient"}[ob])
    ax.set_yticks(np.arange(len(NAMES)))
    ax.set_yticklabels([n.replace("_", r"\_") if False else n for n in NAMES])
    ax.set_xscale("log"); ax.invert_yaxis()
    ax.axvline(0.05, color=MUTED, ls="--", lw=0.8)
    ax.annotate("5% of the budget", xy=(0.05, 1.0), xycoords=("data", "axes fraction"),
                xytext=(-4, -10), textcoords="offset points", ha="right",
                fontsize=7.5, color=MUTED)
    ax.set_xlabel("share of the observable's error budget")
    ax.set_title("one parameter carries everything; a tolerance carries nothing",
                 loc="left")
    ax.set_xlim(1e-26, 3)
    ax.xaxis.set_major_locator(LogLocator(numticks=8))
    ax.grid(axis="x", alpha=.5, lw=.4); ax.legend(frameon=False, loc="lower left")
    fig.savefig(OUT + "fig2_budget.pdf"); plt.close(fig)
    print("fig2 done")


# ------------------------------------------------------------------ Fig 3
def fig_trajectory(tensor, LAD, HS):
    ob = "gradient"; QS = Burgers2D(nu=0.05, T=1.0).exact_observables()[ob]
    ALL = list(itertools.product(*[range(len(l)) for l in LAD]))
    look = lambda i: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
    tr = lambda i: abs(tensor[i][0][ob] - QS)
    eps = 1e-2 * abs(QS)
    start = (0, 3, 0, 2, 2)
    r = run_protocol(look, LAD, start, eps, hs=HS)
    trace = [t for i, t in enumerate(r["trace"]) if i == 0 or t != r["trace"][i - 1]]
    feas = [i for i in ALL if tr(i) <= eps]
    opt = min(feas, key=lambda i: tensor[i][1])

    fig, ax = plt.subplots(figsize=(6.6, 3.1))
    ax.scatter([tensor[i][1] for i in ALL], [tr(i) for i in ALL], s=5,
               color="#c9c9c9", edgecolor="none", zorder=1, label="all 432 settings")
    ax.axhline(eps, color=INK, lw=1.0, ls="--", zorder=2)
    xs = [tensor[t][1] for t in trace]; ys = [tr(t) for t in trace]
    ax.plot(xs, ys, "-", color=C[0], lw=1.3, alpha=.85, zorder=3)
    ax.scatter(xs[1:-1], ys[1:-1], s=16, color=C[0], edgecolor="white",
               linewidth=0.6, zorder=4, label="steps taken by the method")
    ax.scatter([xs[0]], [ys[0]], s=70, marker="s", color=C[1], edgecolor="white",
               linewidth=0.8, zorder=5, label="where the user started")
    ax.scatter([xs[-1]], [ys[-1]], s=95, marker="*", color=C[2], edgecolor="white",
               linewidth=0.8, zorder=6, label="what the method recommends")
    ax.scatter([tensor[opt][1]], [tr(opt)], s=70, marker="o", facecolor="none",
               edgecolor=C[3], linewidth=1.5, zorder=5, label="true best (exhaustive search)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.annotate("accuracy the user asked for", xy=(0.015, eps),
                xycoords=("axes fraction", "data"), xytext=(0, 3),
                textcoords="offset points", fontsize=7.5, color=INK, va="bottom")
    ax.set_xlabel("cost of one run"); ax.set_ylabel("error in the observable")
    ax.set_title("the method walks from the user's guess to a setting it can certify",
                 loc="left")
    ax.grid(alpha=.5, lw=.4); ax.legend(frameon=False, loc="lower left", fontsize=7.5)
    fig.savefig(OUT + "fig3_trajectory.pdf"); plt.close(fig)
    print(f"fig3 done  start={start} -> {r['proposal']}  opt={opt}  steps={len(trace)}")


# ------------------------------------------------------------------ Fig 4
def fig_outcomes(tensor, LAD, HS):
    ALL = list(itertools.product(*[range(len(l)) for l in LAD]))
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.5))
    for oi, ob in enumerate(["energy", "gradient"]):
        QS = Burgers2D(nu=0.05, T=1.0).exact_observables()[ob]
        look = lambda i, ob=ob: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
        tr = lambda i: abs(tensor[i][0][ob] - QS)
        eps = 1e-2 * abs(QS)
        c_opt = min(tensor[i][1] for i in ALL if tr(i) <= eps)
        saved, pen = [], []
        for a in ALL:
            r = run_protocol(look, LAD, a, eps, hs=HS)
            p = r["proposal"]
            pen.append(tensor[p][1] / c_opt)
            if tr(a) <= eps:
                saved.append(tensor[a][1] / tensor[p][1])
        lab = {"energy": "energy", "gradient": "gradient"}[ob]
        for k, (dat, a_) in enumerate([(saved, ax[0]), (pen, ax[1])]):
            d = np.sort(np.array(dat))
            a_.step(d, np.arange(1, len(d) + 1) / len(d), where="post",
                    color=C[oi], label=lab)
        print(f"  {ob}: median saved {np.median(saved):.2f}x, median penalty {np.median(pen):.2f}x")
    ax[0].set_xscale("log")
    ax[0].set(xlabel="cost saved vs the user's own setting",
              ylabel="fraction of users")
    ax[0].set_title("users already accurate: what they save", loc="left")
    ax[0].axvline(1, color=MUTED, lw=0.8, ls=":")
    ax[1].set(xlabel=r"cost relative to the true best setting", ylabel="")
    ax[1].set_title("how close to optimal the advice is", loc="left")
    ax[1].axvline(1, color=MUTED, lw=0.8, ls=":")
    ax[1].set_xlim(0.9, 8)
    for a in ax:
        a.grid(alpha=.5, lw=.4); a.set_ylim(0, 1.02)
        a.legend(frameon=False, loc="lower right", fontsize=7.5)
    fig.savefig(OUT + "fig4_outcomes.pdf"); plt.close(fig)
    print("fig4 done")


# ------------------------------------------------------------------ Fig 5
def fig_transfer():
    labels = ["grid\nspacing", "time\nstep", "Newton\ntolerance",
              "linear-solver\ntolerance"]
    kind = ["discretisation", "discretisation", "nested solver", "nested solver"]
    fine_to_coarse = [0, 0, 0, 88 / 125]         # gmres_tol threshold too coarse
    coarse_to_fine = [0, 0, 9 / 24, 0]           # newton_tol missed by cheap screen
    x = np.arange(4); w = 0.36
    fig, ax = plt.subplots(figsize=(6.6, 2.4))
    b1 = ax.bar(x - w / 2, fine_to_coarse, w, color=C[0], edgecolor="white", linewidth=1.2,
                label="calibrated on a fine run, used on a coarse one")
    b2 = ax.bar(x + w / 2, coarse_to_fine, w, color=C[1], edgecolor="white", linewidth=1.2,
                label="calibrated on a coarse run, used on a fine one")
    for bars in (b1, b2):
        for b in bars:
            if b.get_height() > 0:
                ax.annotate(f"{b.get_height()*100:.0f}%",
                            (b.get_x() + b.get_width() / 2, b.get_height()),
                            xytext=(0, 2), textcoords="offset points",
                            ha="center", fontsize=8, color=INK)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    for i, k in enumerate(kind):
        ax.annotate(k, (i, -0.30), xycoords=("data", "axes fraction"),
                    ha="center", fontsize=7, color=MUTED, style="italic",
                    annotation_clip=False)
    ax.set_ylabel("settings where the\nverdict did not carry over")
    ax.set_title("what a parameter does in one regime need not hold in another",
                 loc="left")
    ax.set_ylim(0, 0.85); ax.grid(axis="y", alpha=.5, lw=.4)
    ax.legend(frameon=False, loc="upper left", fontsize=7.5)
    fig.subplots_adjust(bottom=0.32)
    fig.savefig(OUT + "fig5_transfer.pdf"); plt.close(fig)
    print("fig5 done")


if __name__ == "__main__":
    tensor = pickle.load(open(os.path.join(os.path.dirname(__file__), "data", "burgers2d_tensor.pkl"), "rb"))
    LAD = [[16, 32, 64, 128], [10, 20, 40, 80], [1e-2, 1e-5, 1e-8],
           [1e-2, 1e-5, 1e-8], [1e-1, 1e-2, 1e-3]]
    NAMES = ["grid spacing", "time step", "Newton tol.", "linear-solver tol.",
             "preconditioner"]
    HS = [[2 * np.pi / n for n in LAD[0]], [1.0 / s for s in LAD[1]]] + \
         [list(l) for l in LAD[2:]]
    fig_masking()
    fig_budget(tensor, LAD, HS, NAMES)
    fig_trajectory(tensor, LAD, HS)
    fig_outcomes(tensor, LAD, HS)
    fig_transfer()
