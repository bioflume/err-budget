"""What happens when you declare a parameter's kind wrongly.

The only declaration that is unsafe to get wrong is calling an influential
parameter 'cost-only'.  Here we do exactly that to the mesh, which in this
problem carries the entire error budget, and watch the tool contradict us.
"""
from errbudget import Parameter, tune
from errbudget.problems import Burgers1D

problem = Burgers1D(nu=0.05, T=1.0)
EXACT = problem.exact_observables()["dissipation"]

def solver(N, steps, newton_tol, gmres_tol):
    obs, cost, _, _ = problem.run(dict(N=N, steps=steps,
                                       newton_tol=newton_tol, gmres_tol=gmres_tol))
    return obs, cost

start = dict(N=64, steps=100, newton_tol=1e-6, gmres_tol=1e-6)
for kind in ["discretisation", "cost-only"]:
    params = [Parameter("N", [32, 64, 128, 256, 512], kind=kind),
              Parameter("steps", [25, 50, 100, 200, 400]),
              Parameter("newton_tol", [1e-2, 1e-4, 1e-6, 1e-8, 1e-10], kind="nested-solver"),
              Parameter("gmres_tol", [1e-2, 1e-4, 1e-6, 1e-8, 1e-10], kind="nested-solver")]
    r = tune(solver, params, target=3e-3, observable="dissipation", start=start,
             screen=True)       # the preliminary sweep is where declarations are tested
    obs, _, _, _ = problem.run(r.settings)
    err = abs(obs["dissipation"] - EXACT)
    print(f"N declared {kind!r:>18}: certified={r.certified}, "
          f"true error {err:.2e} vs target {r.target:.2e} "
          f"-> {'MET' if err <= r.target else 'MISSED'}"
          f"{'   [declaration contradicted]' if r.violations else ''}")
