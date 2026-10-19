"""Does the certificate survive a learned component?"""
import numpy as np
from knobs import Knob, tune
from knobs.problems import MLHybridBurgers1D

p = MLHybridBurgers1D()
EXACT = p.exact_observables()["dissipation"]

def solver(N, steps, gmres_tol, model_level):
    obs, cost, _, _ = p.run(dict(N=N, steps=steps, newton_tol=1e-10,
                                 gmres_tol=gmres_tol, model_level=model_level))
    return obs, cost

KN = [Knob("N", [32, 64, 128, 256, 512]),
      Knob("steps", [25, 50, 100, 200, 400]),
      Knob("gmres_tol", [1e-2, 1e-4, 1e-6, 1e-8, 1e-10], kind="nested-solver"),
      Knob("model_level", [0, 1, 2, 3, 4, 5], kind="discretisation",
           h=lambda L: p.surrogate_error(L))]

print("the surrogate's own error does not go to zero:")
for L in [0, 1, 2, 3, 4, 5]:
    print(f"  level {L}: relative error in the learned viscosity "
          f"{p.surrogate_error(L):.3e}")
print(f"  ... floor = {p.delta_floor:.1e}, reached at level ~5\n")

print("error attainable with the surrogate fully trained (level 5), "
      "classical parameters at their finest:")
obs, _, _, _ = p.run(dict(N=512, steps=400, newton_tol=1e-10, gmres_tol=1e-10,
                          model_level=5))
floor_err = abs(obs["dissipation"] - EXACT)
print(f"  |Q - Q*| = {floor_err:.3e}   <- no configuration can do better\n")

start = dict(N=64, steps=50, gmres_tol=1e-6, model_level=1)
for target in [3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5]:
    r = tune(solver, KN, target=target, observable="dissipation", start=start)
    obs, _, _, _ = p.run(dict(**r.settings, newton_tol=1e-10))
    err = abs(obs["dissipation"] - EXACT)
    verdict = "certified" if r.certified else "REFUSED (target unreachable)"
    reach = "attainable" if target > floor_err else "BELOW THE SURROGATE FLOOR"
    ok = "met" if err <= target else "missed"
    print(f"target {target:.0e} ({reach:<25}) -> {verdict:<30} "
          f"true error {err:.2e} ({ok})")
    if r.certified and r.limiting:
        print(f"   limited by: {r.limiting}")
