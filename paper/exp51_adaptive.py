"""Item: a solver with an adaptive, error-controlled time stepper.

L2's Burgers equation (nu = 0.05, T = 1, exact by Cole-Hopf), discretised in
space by second-order central differences on N points and integrated in time by
SciPy's variable-order, variable-step BDF method (solve_ivp), whose step size
is chosen by a local-error controller.  The time step is no longer a parameter;
its role is taken by the controller's tolerances:

    N          grid points                         discretisation
    rtol       controller relative tolerance       ?
    atol       controller absolute tolerance       ?
    max_step   cap on the step size                ?

The question is which class the controller tolerances behave like.  Two
observables: the energy at T (an integral) and max |u_x| at T (an extremum,
refined by a quadratic through the three nodes around the grid maximum).
Cost is counted deterministically: right-hand-side evaluations times N plus
LU factorisations times 10 N (the Jacobian is banded/periodic).
"""
import itertools
import os
import pickle
import sys
import time

import numpy as np
import scipy.sparse as sp
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(__file__))
from errbudget.problems import Burgers1D
from replay import grid, method_cell, baseline_cell, fmt_method
from exp47_transfer import failure_rate, PHIS

NU, T = 0.05, 1.0
NAMES = ["N", "rtol", "atol", "max_step"]
LAD = [[32, 64, 128, 256, 512], [1e-2, 1e-4, 1e-6, 1e-8],
       [1e-4, 1e-6, 1e-8, 1e-10], [1.0, 0.1, 0.03, 0.01]]
HS = [[2 * np.pi / n for n in LAD[0]], list(LAD[1]), list(LAD[2]), list(LAD[3])]
CACHE = os.path.join(os.path.dirname(__file__), "data", "adaptive_tensor.pkl")


def peak(v):
    i = int(np.argmax(v))
    fm, f0, fp = v[i - 1], v[i], v[(i + 1) % len(v)]
    den = fm - 2 * f0 + fp
    return float(f0 - 0.125 * (fm - fp) ** 2 / den) if den < 0 else float(f0)


def run(theta):
    N = int(theta["N"])
    h = 2 * np.pi / N
    x = h * np.arange(N)
    def f(t, u):
        up, um = np.roll(u, -1), np.roll(u, 1)
        return -(up**2 - um**2) / (4 * h) + NU * (up - 2 * u + um) / h**2
    o = np.ones(N)
    S = sp.diags([o, o, o], [-1, 0, 1], shape=(N, N), format="lil")
    S[0, -1] = S[-1, 0] = 1
    sol = solve_ivp(f, (0.0, T), np.sin(x), method="BDF", rtol=theta["rtol"],
                    atol=theta["atol"], max_step=theta["max_step"],
                    jac_sparsity=S.tocsc())
    u = sol.y[:, -1]
    ux = (np.roll(u, -1) - np.roll(u, 1)) / (2 * h)
    obs = {"energy": 0.5 * float(np.sum(u**2)) * h, "peak_gradient": peak(np.abs(ux))}
    cost = sol.nfev * N + sol.nlu * 10 * N
    return obs, float(cost), dict(nfev=sol.nfev, nlu=sol.nlu, steps=len(sol.t) - 1,
                                  ok=bool(sol.success))


b = Burgers1D(nu=NU, T=T)
xe, ue, uxe = b.exact(T)
EX = {"energy": b.exact_observables()["energy"], "peak_gradient": float(np.max(np.abs(uxe)))}

if os.path.exists(CACHE):
    tensor = pickle.load(open(CACHE, "rb"))
else:
    tensor, t0 = {}, time.time()
    for idx in grid(LAD):
        obs, c, info = run({NAMES[k]: LAD[k][j] for k, j in enumerate(idx)})
        tensor[idx] = (obs, c if info["ok"] else float("inf"), info)
    pickle.dump(tensor, open(CACHE, "wb"))
    print(f"tensor: {len(tensor)} solve_ivp runs in {time.time()-t0:.0f}s")

print("=" * 100)
print("Burgers (nu = 0.05) with SciPy's adaptive BDF integrator;  exact values", EX)
print("=" * 100)

print("\nerror at the finest controller settings, by N:")
for ob in EX:
    e = [abs(tensor[(j, 3, 3, 3)][0][ob] - EX[ob]) / abs(EX[ob]) for j in range(5)]
    print(f"  {ob:<14} " + "  ".join(f"{v:.1e}" for v in e))
print("\nerror against rtol at N = 512, atol = 1e-10, max_step = 0.01:")
for ob in EX:
    e = [abs(tensor[(4, j, 3, 3)][0][ob] - EX[ob]) / abs(EX[ob]) for j in range(4)]
    print(f"  {ob:<14} " + "  ".join(f"{v:.1e}" for v in e))

def get(idx, ob):
    v = tensor[idx]
    return v[0][ob] if np.isfinite(v[1]) else np.nan
get.exact = EX
print("\ntransfer of thresholds (as Table 3): failure rate per parameter")
print(f"  {'parameter':<10} " + " ".join(f"phi={p:<6g}" for p in PHIS))
for k in range(4):
    row = []
    for phi in PHIS:
        fl, n = failure_rate(get, LAD, k, list(EX), phi)
        row.append(f"{100*fl/max(n,1):>9.0f}%")
    print(f"  {NAMES[k]:<10} " + " ".join(row))

for decl in (["discretisation", "nested-solver", "nested-solver", "discretisation"],
             ["discretisation", "discretisation", "discretisation", "discretisation"]):
    print(f"\ndeclared classes: {dict(zip(NAMES, decl))}")
    for ob in EX:
        for tg in [1e-2, 1e-3]:
            eps = tg * abs(EX[ob])
            o = method_cell(tensor, LAD, HS, decl, ob, EX[ob], eps)
            print(f"  {ob:<14} {tg:.0e}  " + fmt_method(o))
            if decl[1] == "nested-solver":
                for lab, r in baseline_cell(tensor, LAD, ob, EX[ob], eps).items():
                    print(f"      {lab:<28} attained {r['attained']}/{r['n']}, "
                          f"penalty median {r['pen']:.2f}, runs median {r['runs']:.0f}")
