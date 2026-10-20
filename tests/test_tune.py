"""Behavioural guarantees. These are the claims the paper makes."""
import numpy as np, pytest
from errbudget import Parameter, tune
from errbudget.problems import Burgers1D

problem = Burgers1D(nu=0.05, T=1.0)
EXACT = problem.exact_observables()

def solver(N, steps, newton_tol, gmres_tol):
    obs, cost, _, _ = problem.run(dict(N=N, steps=steps,
                                       newton_tol=newton_tol, gmres_tol=gmres_tol))
    return obs, cost

PARAMS = [Parameter("N", [32, 64, 128, 256, 512]),
          Parameter("steps", [25, 50, 100, 200, 400]),
          Parameter("newton_tol", [1e-2, 1e-4, 1e-6, 1e-8, 1e-10], kind="nested-solver"),
          Parameter("gmres_tol", [1e-2, 1e-4, 1e-6, 1e-8, 1e-10], kind="nested-solver")]

@pytest.mark.parametrize("start", [
    dict(N=32,  steps=25,  newton_tol=1e-2,  gmres_tol=1e-2),
    dict(N=128, steps=100, newton_tol=1e-10, gmres_tol=1e-10),
    dict(N=512, steps=400, newton_tol=1e-6,  gmres_tol=1e-6),
])
def test_recommendation_meets_target(start):
    r = tune(solver, PARAMS, target=3e-3, observable="dissipation", start=start)
    assert r.certified
    obs, _, _, _ = problem.run(r.settings)
    assert abs(obs["dissipation"] - EXACT["dissipation"]) <= r.target

def test_unreachable_target_is_not_certified():
    r = tune(solver, PARAMS, target=1e-9, observable="energy",
             start=dict(N=64, steps=50, newton_tol=1e-4, gmres_tol=1e-4))
    assert not r.certified

def test_diagnostic_agrees_with_recommendation():
    """Nothing reported as wasted may be something the tool declined to loosen
    for cost reasons -- the advice must not contradict itself."""
    r = tune(solver, PARAMS, target=3e-3, observable="dissipation",
             start=dict(N=256, steps=200, newton_tol=1e-10, gmres_tol=1e-10))
    for line in r.wasted:
        assert "would still meet the target" in line

def test_misdeclared_cost_only_parameter_is_reported():
    bad = list(PARAMS)
    bad[0] = Parameter("N", [32, 64, 128, 256, 512], kind="cost-only")
    r = tune(solver, bad, target=3e-3, observable="dissipation",
             start=dict(N=64, steps=100, newton_tol=1e-6, gmres_tol=1e-6))
    assert "N" in r.violations or "N" not in r.frozen


def test_reaction_diffusion_recommendation_meets_target():
    """A different domain and a different solver architecture: operator
    splitting, no Newton iteration, three competing discretisation errors."""
    from errbudget.problems import ReactionDiffusion2D
    rd = ReactionDiffusion2D()
    ex = rd.exact_observables()["inhibitor_mass"]

    def rd_solver(N, steps, n_react, gmres_tol, ilu_drop):
        obs, cost, _, _ = rd.run(dict(N=N, steps=steps, n_react=n_react,
                                      gmres_tol=gmres_tol, ilu_drop=ilu_drop))
        return obs, cost

    ks = [Parameter("N", [16, 32, 64]),
          Parameter("steps", [5, 10, 20]),
          Parameter("n_react", [1, 2, 4]),
          Parameter("gmres_tol", [1e-2, 1e-5, 1e-8], kind="nested-solver"),
          Parameter("ilu_drop", [1e-1, 1e-2, 1e-3], kind="cost-only")]
    r = tune(rd_solver, ks, target=1e-2, relative=True, observable="inhibitor_mass",
             start=dict(N=32, steps=10, n_react=4, gmres_tol=1e-8, ilu_drop=1e-3))
    assert r.certified
    obs, _, _, _ = rd.run(r.settings)
    assert abs(obs["inhibitor_mass"] - ex) <= r.target
