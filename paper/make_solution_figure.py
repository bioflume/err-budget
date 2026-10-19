"""Figure: the problems actually being solved, and what under-resolution does."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from knobs.problems import Burgers1D, Burgers2D, ReactionDiffusion2D, DrivenCavity

OUT = os.environ.get("KNOBS_FIG_DIR", "figures/")
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

fig, ax = plt.subplots(1, 4, figsize=(9.4, 2.35))

# --- L2: 1D Burgers, exact vs a coarse run -----------------------------------
b1 = Burgers1D(nu=0.05, T=1.0)
xf, uf, _ = b1.exact(b1.T)
ax[0].plot(xf, uf, "-", color=INK, lw=1.4, label="exact")
for N, col, ms in [(16, C[1], 4.2), (64, C[0], 3.0)]:
    x = 2 * np.pi * np.arange(N) / N
    Dx, Lx = b1._ops(N)
    u = np.sin(x)
    dt = b1.T / 100
    for _ in range(100):
        un = u.copy(); Fn = b1._F(un, Dx, Lx); v = un.copy()
        for _ in range(20):
            G = v - un - 0.5 * dt * (b1._F(v, Dx, Lx) + Fn)
            if np.linalg.norm(G) <= 1e-10 * (np.linalg.norm(un) + 1e-30):
                break
            import scipy.sparse.linalg as spla
            dv, _ = spla.gmres(b1._J(v, Dx, Lx, dt, N), -G, rtol=1e-10, atol=0.0,
                               restart=min(50, N), maxiter=400)
            v = v + dv
        u = v
    ax[0].plot(np.append(x, 2*np.pi), np.append(u, u[0]), "o--", ms=ms, lw=0.8,
               color=col, markeredgecolor="white", markeredgewidth=0.4,
               label=f"$N={N}$")
ax[0].set(xlabel="$x$", ylabel="$u$", xlim=(0, 2 * np.pi))
ax[0].set_title("(a) L2: Burgers, 1D", loc="left")
ax[0].set_xticks([0, np.pi, 2 * np.pi]); ax[0].set_xticklabels(["0", r"$\pi$", r"$2\pi$"])
ax[0].legend(frameon=False, fontsize=7.5, loc="lower left")

# --- L3: 2D Burgers field ----------------------------------------------------
b2 = Burgers2D(nu=0.05, T=1.0)
N = 128
x = 2 * np.pi * np.arange(N) / N
X, Y = np.meshgrid(x, x, indexing="ij")
U = b2.u_ex(X, Y, b2.T)
im = ax[1].pcolormesh(X, Y, U, cmap="RdBu_r", shading="auto",
                      vmin=-abs(U).max(), vmax=abs(U).max(), rasterized=True)
ax[1].set(xlabel="$x$", ylabel="$y$", aspect="equal")
ax[1].set_title("(b) L3: Burgers, 2D", loc="left")
plt.colorbar(im, ax=ax[1], fraction=0.046, pad=0.03).ax.tick_params(labelsize=7)

# --- L4: reaction-diffusion inhibitor field ----------------------------------
rd = ReactionDiffusion2D()
V = rd.v_ex(X, Y, rd.T)
im2 = ax[2].pcolormesh(X, Y, V, cmap="viridis", shading="auto", rasterized=True)
ax[2].set(xlabel="$x$", ylabel="$y$", aspect="equal")
ax[2].set_title("(c) L4: Gray–Scott, $v$", loc="left")
plt.colorbar(im2, ax=ax[2], fraction=0.046, pad=0.03).ax.tick_params(labelsize=7)

for a in (ax[1], ax[2]):
    a.set_xticks([0, np.pi, 2 * np.pi]); a.set_xticklabels(["0", r"$\pi$", r"$2\pi$"])
    a.set_yticks([0, np.pi, 2 * np.pi]); a.set_yticklabels(["0", r"$\pi$", r"$2\pi$"])
# --- L5: driven cavity, streamfunction with the boundary layer --------------
cv = DrivenCavity()
xs = np.linspace(0, 1, 201)
XC, YC = np.meshgrid(xs, xs, indexing="ij")
PS = cv.psi_ex(XC, YC)
im3 = ax[3].contourf(XC, YC, PS, levels=14, cmap="Blues_r")
ex = cv.exact_observables()
ax[3].plot([0.5], [ex["y_min"]], "x", color=C[1], ms=6, mew=1.4)
ax[3].set(xlabel="$x$", ylabel="$y$", aspect="equal")
ax[3].set_title(r"(d) L5: cavity, $\psi$", loc="left")
plt.colorbar(im3, ax=ax[3], fraction=0.046, pad=0.03).ax.tick_params(labelsize=7)
fig.subplots_adjust(wspace=0.6)
fig.savefig(OUT + "fig0_problems.pdf")
print("solution figure done")
