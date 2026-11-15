"""fig_baselines.pdf: the program against the three procedures of current practice.

Data source: paper/tab_baselines_full.tex (Table S-4.2), parsed as written; no other input.
Columns of each row: group, procedure, Met, Missed, Declined, Cost ratio, Runs, K.
  (a) fraction of searches from converged starting settings whose returned setting misses
      = Missed / (Met + Missed + Declined); Declined is "---" (zero) for the procedures.
  (b) median cost ratio of the program and of SJ at its largest threshold with no miss,
      with K beside each point; the legend title says so.
  python3 figures/make_baselines_figure.py
"""
import os
import re

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap

from figstyle import use_style, panel_label, TEXT, BLUE, VERM, BLACK, GREY

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
SRC = os.path.join(OUT, "tab_baselines_full.tex")


def num(s):
    s = s.replace("{,}", "").strip()
    return 0.0 if s == "---" else float(s)


def parse():
    rows, group = {}, None
    for line in open(SRC):
        line = line.strip()
        if not line or line.startswith("\\addlinespace"):
            continue
        cells = [c.strip() for c in line.rstrip("\\").split("&")]
        if cells[0]:
            group = re.sub(r"\s*\(\d+ cases\)", "", cells[0])
            rows[group] = {}
        met, missed, declined, ratio, runs, K = map(num, cells[2:8])
        rows[group][cells[1]] = dict(met=met, missed=missed, declined=declined, ratio=ratio, K=K)
    return rows


COLS = [("S", r"$\varepsilon$"), ("S", r"$\varepsilon/3$"), ("S", r"$\varepsilon/10$"),
        ("SJ", r"$\varepsilon$"), ("SJ", r"$\varepsilon/3$"), ("SJ", r"$\varepsilon/10$"),
        ("Ref.", None), ("This program", None)]
ETAS = {r"$\varepsilon$": 1, r"$\varepsilon/3$": 3, r"$\varepsilon/10$": 10}


def key(proc, eta):
    return proc if eta is None else f"{proc}, $\\eta={eta.strip('$')}$"


def main():
    use_style()
    R = parse()
    groups = list(R)
    ng = len(groups)
    frac = np.zeros((ng, len(COLS)))
    for i, g in enumerate(groups):
        for j, (p, eta) in enumerate(COLS):
            r = R[g][key(p, eta)]
            frac[i, j] = r["missed"] / (r["met"] + r["missed"] + r["declined"])

    fig = plt.figure(figsize=(TEXT, 3.0))
    ax = fig.add_axes([0.155, 0.20, 0.43, 0.62])
    cax = fig.add_axes([0.24, 0.115, 0.26, 0.028])
    bx = fig.add_axes([0.63, 0.20, 0.36, 0.62], sharey=ax)

    # (a) dot matrix
    # four bins of the fraction, light to dark; the first bin is (0, 1%]
    cmap = ListedColormap(["#F7CBA9", "#EE9154", "#C24F00", "#5A2400"])
    bounds = [0.0, 0.01, 0.1, 0.3, 1.0]
    norm = BoundaryNorm(bounds, cmap.N)
    y = np.arange(ng)
    for j in range(len(COLS)):
        for i in range(ng):
            f = frac[i, j]
            if f == 0:
                ax.plot(j, i, "o", ms=6.5, mfc="white", mec=GREY, mew=0.8)
            else:
                ax.plot(j, i, "o", ms=6.5, mfc=cmap(norm(f)), mec="none")
    ax.set_xlim(-0.6, len(COLS) - 0.4)
    ax.set_ylim(ng - 0.5, -0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([g.replace("--", "\u2013") for g in groups])
    ax.tick_params(axis="y", length=0)
    ax.set_xticks(range(len(COLS)))
    ax.set_xticklabels([e if e else "" for _, e in COLS])
    ax.tick_params(axis="x", length=0)
    ax.xaxis.tick_top()
    for sp in ax.spines.values():
        sp.set_visible(False)
    for xc, lab in [(1, "S, $\\eta=$"), (4, "SJ, $\\eta=$"), (6, "Ref."), (7, "errbudget")]:
        ax.text(xc, -1.55 if lab.endswith("=$") else -1.1, lab, ha="center", va="bottom",
                fontsize=8, color=BLUE if lab == "errbudget" else BLACK)
    for xc in (0.5 + 2, 0.5 + 5, 0.5 + 6):
        ax.axvline(xc, color="#DDDDDD", lw=0.6, zorder=0)
    panel_label(ax, "(a)", x=-0.30, y=1.12)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal",
                      ticks=bounds, spacing="uniform")
    cb.ax.set_xticklabels(["0", "1%", "10%", "30%", "100%"])
    cb.set_label("fraction of searches that miss (open: none)", fontsize=7.5, labelpad=2)
    cb.ax.tick_params(labelsize=7)
    cb.outline.set_linewidth(0.4)

    # (b) cost ratio of the program and of SJ at its largest safe threshold
    out = []
    for i, g in enumerate(groups):
        pr = R[g]["This program"]
        safe = [e for e in (r"$\varepsilon$", r"$\varepsilon/3$", r"$\varepsilon/10$")
                if R[g][key("SJ", e)]["missed"] == 0]
        bx.plot(pr["ratio"], i, "o", color=BLUE, ms=5, zorder=5)
        if safe:
            e = safe[0]
            sj = R[g][key("SJ", e)]
            left = pr["ratio"] <= sj["ratio"]          # program's K away from SJ's point
            bx.plot([pr["ratio"], sj["ratio"]], [i, i], color="#BBBBBB", lw=1.0, zorder=2)
            bx.plot(sj["ratio"], i, "s", color=VERM, ms=4.5, zorder=4)
            bx.annotate(f"{pr['K']:,g}", (pr["ratio"], i), xytext=(-5 if left else 5, 0),
                        textcoords="offset points", fontsize=6.5, color=BLUE,
                        ha="right" if left else "left", va="center")
            sjlab = f"{sj['K']:g}" + ("" if e == r"$\varepsilon/3$" else f" (at $\\eta=${e})")
            bx.annotate(sjlab, (sj["ratio"], i), xytext=(5 if left else -5, 0),
                        textcoords="offset points", fontsize=6.5, color=VERM,
                        ha="left" if left else "right", va="center")
            out.append((g, pr["ratio"], pr["K"], e, sj["ratio"], sj["K"]))
        else:
            bx.annotate(f"{pr['K']:,g}", (pr["ratio"], i), xytext=(-5, 0),
                        textcoords="offset points", fontsize=6.5, color=BLUE, ha="right",
                        va="center")
            bx.annotate("SJ: no safe threshold", (pr["ratio"], i), xytext=(6, 0),
                        textcoords="offset points", fontsize=7, color=VERM, ha="left",
                        va="center", style="italic")
            out.append((g, pr["ratio"], pr["K"], None, None, None))
    bx.set_xscale("log")
    bx.set_xlim(0.55, 250)
    bx.set_xticks([1, 3, 10, 30, 100])
    bx.set_xticklabels(["1", "3", "10", "30", "100"])
    bx.minorticks_off()
    bx.xaxis.tick_top()
    bx.xaxis.set_label_position("top")
    bx.set_xlabel("median cost ratio (cost / optimal cost)", labelpad=4)
    bx.tick_params(axis="y", left=False, labelleft=False)
    bx.spines["left"].set_visible(False)
    bx.spines["bottom"].set_visible(False)
    bx.spines["top"].set_visible(True)
    for i in range(ng):
        bx.axhline(i, color="#F0F0F0", lw=0.6, zorder=0)
    panel_label(bx, "(b)", x=-0.02, y=1.12)
    h = [plt.Line2D([], [], color=BLUE, marker="o", ls="none", ms=5),
         plt.Line2D([], [], color=VERM, marker="s", ls="none", ms=4.5)]
    bx.legend(h, ["errbudget", r"SJ, largest $\eta$ with no miss"], loc="upper center",
              bbox_to_anchor=(0.5, -0.03), ncol=2, fontsize=7.5, handletextpad=0.3,
              title="numbers: search cost $K$ (production runs)", title_fontsize=7.5)

    fig.savefig(os.path.join(OUT, "fig_baselines.pdf"))
    fig.savefig(os.path.join(OUT, "fig_baselines.png"), dpi=200)
    plt.close(fig)
    print("wrote fig_baselines.pdf")
    print("(a) fraction missed; columns:", [f"{p} {e or ''}" for p, e in COLS])
    for g, row in zip(groups, frac):
        print(f"  {g:18s}", " ".join(f"{v:.4f}" for v in row))
    print("(b) group, program ratio, program K, SJ threshold, SJ ratio, SJ K, SJ/program:")
    for g, a, ka, e, b, kb in out:
        print(f"  {g:18s} {a:6.3g} {ka:6g}  {e or '-':18s} {b if b else '-'} {kb if kb else '-'}"
              f" {b / a if b else '-'}")


if __name__ == "__main__":
    main()
