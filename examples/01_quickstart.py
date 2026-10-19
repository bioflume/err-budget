"""Quickstart: hand `tune` a solver and ask for an accuracy.

The solver here is 1D viscous Burgers, which has an exact solution via the
Cole-Hopf transformation, so we can check afterwards whether the answer really
met the target.  `tune` itself never sees the exact solution.
"""
from knobs import Knob, tune
from knobs.problems import Burgers1D

problem = Burgers1D(nu=0.05, T=1.0)


def solver(N, steps, newton_tol, gmres_tol):
    """Anything with this shape will do: keywords in, observables out."""
    obs, cost, wall, _ = problem.run(
        dict(N=N, steps=steps, newton_tol=newton_tol, gmres_tol=gmres_tol))
    return obs, cost


knobs = [
    Knob("N",          [32, 64, 128, 256, 512],            kind="discretisation"),
    Knob("steps",      [25, 50, 100, 200, 400],            kind="discretisation"),
    Knob("newton_tol", [1e-2, 1e-4, 1e-6, 1e-8, 1e-10],    kind="nested-solver"),
    Knob("gmres_tol",  [1e-2, 1e-4, 1e-6, 1e-8, 1e-10],    kind="nested-solver"),
]

# A configuration a cautious user might well have picked: a middling mesh, and
# both tolerances wound right down because tighter feels safer.
start = dict(N=128, steps=100, newton_tol=1e-10, gmres_tol=1e-10)

result = tune(solver, knobs, target=1e-3, observable="dissipation", start=start)
print(result)

# --- check it against the exact solution (the tool did not get to see this) ---
exact = problem.exact_observables()["dissipation"]
obs, cost, _, _ = problem.run(result.settings)
before, cost_before, _, _ = problem.run(start)
print(f"\ntrue error at the recommended setting : {abs(obs['dissipation']-exact):.3e}"
      f"   (target {result.target:.3e})")
print(f"true error at the user's own setting   : "
      f"{abs(before['dissipation']-exact):.3e}")
print(f"cost of one run, before -> after       : {cost_before:.3g} -> {cost:.3g}"
      f"   ({cost_before/cost:.1f}x cheaper)")
