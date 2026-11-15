"""Two measurements on L2 (S-4.1): the change of the energy and the dissipation when the
Newton tolerance is tightened from 1e-8 to 1e-10 at the finest mesh and time step, and the
dissipation error of the optimal setting for the energy at 1e-3.

    python exp40_section5.py
"""
import numpy as np
from errbudget import LocalBudget
from grids import L2, L3, hs, load, exact, grid

t, ex = load(L2), exact(L2)
fin = (4, 4, 4, 4)
print("=" * 90)
print("5.1  masking on L2 at the finest mesh and time step (N=512, 400 steps)")
print("=" * 90)
for ob in L2["obs"]:
    ref = t[fin][0][ob]
    disc = abs(ref - ex[ob]) / abs(ex[ob])
    print(f"  {ob}: remaining discretisation error {disc:.2e} relative")
    for k, nm in [(2, "Newton"), (3, "GMRES")]:
        for j, tol in enumerate(L2["lad"][k][:-1]):
            idx = list(fin); idx[k] = j
            d = abs(t[tuple(idx)][0][ob] - ref) / abs(ref)
            print(f"    {nm:<7} tolerance {tol:.0e} vs 1e-10: change {d:.2e} relative"
                  f"  ({100*d/disc:.3f}% of the discretisation error)")

print("\n" + "=" * 90)
print("5.1  masking on L3: budget at the anchor N=64, 40 steps, tolerances 1e-5, drop 1e-2")
print("=" * 90)
t3, H3 = load(L3), hs(L3)
for ob in L3["obs"]:
    look = lambda i, ob=ob: (t3[tuple(i)][0][ob], t3[tuple(i)][1])
    lb = LocalBudget(look, L3["lad"], (2, 2, 1, 1, 1), hs=H3)
    lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
    r = np.abs(np.nan_to_num(np.array(lb.rk, float)))
    print(f"  {ob}: " + ", ".join(f"{n} {100*x/r.sum():.3f}%" for n, x in zip(L3["names"], r)))

print("\n" + "=" * 90)
print("5.3  the same runs, two observables (L2, relative targets)")
print("=" * 90)
ALL = grid(L2)
for rel in [1e-3, 1e-4]:
    best = {}
    for ob in L2["obs"]:
        eps = rel * abs(ex[ob])
        feas = [i for i in ALL if abs(t[i][0][ob] - ex[ob]) <= eps]
        if not feas:
            print(f"  target {rel:.0e}, {ob}: not attainable on the grid"); continue
        b = min(feas, key=lambda i: t[i][1]); best[ob] = b
        tight = [i for i in feas if i[2] == 4 and i[3] == 4]
        print(f"  target {rel:.0e}, {ob}: optimum {b}; cheapest with both tolerances at"
              f" 1e-10 costs {min(t[i][1] for i in tight)/t[b][1]:.1f} times the optimum")
    if "energy" in best:
        b = best["energy"]; e = abs(t[b][0]["dissipation"] - ex["dissipation"])
        print(f"    energy-optimal setting: dissipation error {e:.3e} against "
              f"{rel*ex['dissipation']:.3e} permitted -> "
              f"{'FAILS' if e > rel*ex['dissipation'] else 'meets'}")
