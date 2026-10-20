"""Is the explanation offered for limited transferability actually true?

The paper explains transfer failure by saying that coarsening the time step makes
each inner solve harder, so a tolerance's own error grows just as the
discretisation error that would bury it grows.  That is testable and it has not
been tested.  Two quantities decide it, both measured per SOLVE rather than per
run, since a run with more steps trivially performs more work:

    kappa_proxy(dt) = GMRES matrix-vector products per Newton iteration
    e_tol(dt; tol)  = |Q(tol) - Q(tol_ref)| at fixed mesh and step

If the mechanism is real, both rise as dt rises, and e_tol/tol rises with the
iteration count.  If they do not, the transfer finding stands as a measured
frequency with no mechanism, and the paper must say so.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from errbudget.problems import Burgers1D

p = Burgers1D(nu=0.05, T=1.0)
TOLREF = 1e-12
STEPS = [25, 50, 100, 200, 400]
TOLS = [1e-3, 1e-5, 1e-7]

print("=" * 96)
print("conditioning of the inner solve against the outer time step")
print("=" * 96)
print("  1D Burgers, N = 256, dissipation.  J = I - (dt/2) F'(u): a larger dt")
print("  moves J further from the identity, so the Krylov solve should be harder.")
print()
print(f"  {'steps':>6} {'dt':>8} {'Newton its':>11} {'matvecs':>9} "
      f"{'matvecs/Newton it':>18}")
rows = []
for s in STEPS:
    _, _, _, info = p.run(dict(N=256, steps=s, newton_tol=1e-10,
                               gmres_tol=1e-8))
    per = info["matvecs"] / info["newton_its"]
    rows.append((s, 1.0 / s, info["newton_its"], info["matvecs"], per))
    print(f"  {s:>6} {1.0/s:>8.4f} {info['newton_its']:>11} "
          f"{info['matvecs']:>9} {per:>18.2f}")
dts = np.array([r[1] for r in rows])
per = np.array([r[4] for r in rows])
sl = np.polyfit(np.log(dts), np.log(per), 1)[0]
print(f"\n  matvecs per Newton iteration ~ dt^{sl:.2f}"
      f"   ({'rises' if sl > 0.05 else 'does not rise'} with dt)")

print()
print("=" * 96)
print("the tolerance's own error contribution against the time step")
print("=" * 96)
print("  e_tol = |Q(tol) - Q(1e-12)| at fixed N and steps; the discretisation")
print("  error at the same point is shown for comparison.")
print()
print(f"  {'steps':>6} {'dt':>8} " +
      "  ".join(f"e_tol @{t:.0e}" for t in TOLS) + "   e_discr")
for s in STEPS:
    obs_ref, _, _, _ = p.run(dict(N=256, steps=s, newton_tol=1e-10,
                                  gmres_tol=TOLREF))
    base = obs_ref["dissipation"]
    cells = []
    for t in TOLS:
        o, _, _, _ = p.run(dict(N=256, steps=s, newton_tol=1e-10, gmres_tol=t))
        cells.append(abs(o["dissipation"] - base))
    ed = abs(base - p.exact_observables()["dissipation"])
    print(f"  {s:>6} {1.0/s:>8.4f} " +
          "  ".join(f"{c:>12.3e}" for c in cells) + f"   {ed:>9.3e}")

print()
print("=" * 96)
print("the quantity that decides transfer: is the tolerance MASKED here?")
print("=" * 96)
print("  A threshold transfers only if the ratio e_tol/e_discr keeps the same")
print("  side of 1 at both discretisations.  That ratio is what we tabulate.")
print()
print(f"  {'steps':>6} {'dt':>8} " + "  ".join(f"ratio @{t:.0e}" for t in TOLS))
for s in STEPS:
    obs_ref, _, _, _ = p.run(dict(N=256, steps=s, newton_tol=1e-10,
                                  gmres_tol=TOLREF))
    base = obs_ref["dissipation"]
    ed = abs(base - p.exact_observables()["dissipation"])
    cells = []
    for t in TOLS:
        o, _, _, _ = p.run(dict(N=256, steps=s, newton_tol=1e-10, gmres_tol=t))
        cells.append(abs(o["dissipation"] - base) / max(ed, 1e-300))
    print(f"  {s:>6} {1.0/s:>8.4f} " +
          "  ".join(f"{c:>12.3e}" for c in cells))
print()
print("  A ratio that crosses 1 as dt changes is a threshold that does not")
print("  transfer.  A ratio that stays on one side is one that does.")
