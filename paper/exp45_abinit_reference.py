"""Section 8.3: how far the converged ABINIT reference can be trusted.

The reference (ecut = 48 Ha, 8x8x8 k-points, toldfe = 1e-14) lies beyond the
ladder in every parameter.  It is compared with ecut = 36 and 42 Ha and with a
10x10x10 k-point mesh; the largest change from the last increment estimates its
residual uncertainty, and targets are admitted only if they lie at least an
order of magnitude above it.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from abinit_silicon import AbinitSilicon

p = AbinitSilicon()
ref = p.reference()
vals = {}
for lab, upd in [("ecut 36", dict(ecut=36.0)), ("ecut 42", dict(ecut=42.0)),
                 ("reference", {}), ("ngkpt 10", dict(ngkpt=10))]:
    th = dict(ref); th.update(upd)
    obs, cost, wall, info = p.run(th)
    vals[lab] = obs
    print(f"{lab:<10} {th}: E = {obs['energy']:.10f} Ha, |F| = {obs['force']:.8e} Ha/Bohr ({wall:.0f} s)")
R = vals["reference"]
u = {}
for ob in ["energy", "force"]:
    de = abs(vals["ecut 42"][ob] - R[ob]) / abs(R[ob])
    dk = abs(vals["ngkpt 10"][ob] - R[ob]) / abs(R[ob])
    u[ob] = max(de, dk)
    print(f"{ob}: change 42 -> 48 Ha {de:.2e}, 8 -> 10 k-points {dk:.2e}; uncertainty {u[ob]:.2e}")
for ob, tgs in [("energy", [1e-2, 3e-3, 1e-3]), ("force", [1e-1, 3e-2, 1e-2])]:
    for tg in tgs:
        print(f"  {ob:<6} target {tg:.0e}: {tg/u[ob]:8.1f} times the uncertainty -> "
              f"{'admitted' if tg >= 10*u[ob] else 'not admitted'}")
