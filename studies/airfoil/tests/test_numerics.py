"""The numerics variants: the registered pair unchanged, the sweep's variants as described,
and the readers the sweep's diagnostics use."""
import json
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.dont_write_bytecode = True
from af import foamcase, foamio, study  # noqa: E402

# fvSchemes and fvSolution at the middle setting (residualControl 1e-5, GAMG 1e-7): the
# variants as registered (da6e9e0; the primary is now "primary_v1", the fallback, withdrawn
# by A3, keeps its name) and the primary as amended by A3 (PROTOCOL.md 14)
SNAPSHOT = json.load(open(os.path.join(HERE, "registered_numerics.json")))


@pytest.mark.parametrize("name,now", [("primary", "primary_v1"), ("fallback", "fallback")])
def test_registered_variants_reproduced(name, now):
    run = study.run_at((0.03, 0.12), study.MID)
    v1 = SNAPSHOT["registered_da6e9e0"][name]
    assert foamcase.fv_schemes(now) == v1["schemes"]
    assert foamcase.fv_solution(run, now) == v1["solution"]


def test_amended_primary():
    run = study.run_at((0.03, 0.12), study.MID)
    a3 = SNAPSHOT["amended_A3"]["primary"]
    assert foamcase.fv_schemes("primary") == a3["schemes"]
    assert foamcase.fv_solution(run, "primary") == a3["solution"]


def test_a3_primary_differs_from_v1_only_in_k_omega_convection():
    old = foamcase.fv_schemes("primary_v1").splitlines()
    new = foamcase.fv_schemes("primary").splitlines()
    changed = [(a, b) for a, b in zip(old, new) if a != b]
    assert len(old) == len(new) and [b for _, b in changed] == [
        "    div(phi,k)      bounded Gauss upwind;", "    div(phi,omega)  bounded Gauss upwind;"]
    run = study.run_at((0.03, 0.12), study.MID)
    assert foamcase.fv_solution(run, "primary") == foamcase.fv_solution(run, "primary_v1")


def test_free_stream_of_a3():
    assert study.K_INF == 1e-6 and study.OMEGA_INF == 5.0          # TMR SST-sust, U = c = 1
    assert abs(study.K_INF / study.OMEGA_INF / study.NU - 0.6) < 1e-12
    text = foamcase.turbulence_properties()
    assert "decayControl    yes;" in text and "kInf            1e-06;" in text and "omegaInf        5.0;" in text


@pytest.mark.parametrize("name", sorted(foamcase.VARIANTS))
def test_every_variant_is_complete(name):
    run = study.run_at((0.03, 0.12), study.MID)
    v = foamcase.VARIANTS[name]
    sch, sol = foamcase.fv_schemes(name), foamcase.fv_solution(run, name)
    turb = ["k", "omega"] if foamcase.model_of(name) == "kOmegaSST" else ["nuTilda"]
    for f in ["U"] + turb:
        assert f"div(phi,{f})" in sch
    assert ("grad(U)         cellLimited Gauss linear 1;" in sch) == (v["grad_U"] == "limited")
    assert ("linearUpwind grad(U)" in sch) == (v["div_U"] == "linearUpwind")
    assert f"nNonOrthogonalCorrectors {v['nonorth']};" in sol
    assert ("consistent      yes;" in sol) == v["consistent"]
    assert ("p               0.3;" in sol) == (v["p"] == 0.3)
    assert sch.count("{") == sch.count("}") and sol.count("{") == sol.count("}")


def test_convergence_map():
    runs = study.sweep_convmap()
    assert len(runs) == 4 * 45 + 2 and len({r.tag for r in runs}) == len(runs)
    assert all(r.res_ctrl == 1e-6 and r.gamg_tol == 1e-8 and r.num == "primary" for r in runs[:-2])
    assert runs[-1].setting() == study.reference_setting() and (runs[-2].m, runs[-2].t) == (0.06, 0.16)


def test_tags_of_variants():
    tags = [r.tag for r in study.sweep("numerics")]
    assert len(tags) == len(set(tags)) == 33
    assert sum(t.endswith("_nvsa") for t in tags) == 6
    assert "m6.0_t16_n64_y2e-05_R50_rc1e-06_gt1e-08" in tags            # the primary: no suffix


def test_internal_field_reader(tmp_path):
    p = tmp_path / "U"
    p.write_text("FoamFile { class volVectorField; }\ninternalField   nonuniform List<vector> \n3\n(\n"
                 "(1 0 0)\n(0.5 -1e-3 0)\n(2e-1 3 0)\n)\n;\nboundaryField { }\n")
    u = foamio.internal_field(str(p))
    assert u.shape == (3, 3) and u[1, 1] == -1e-3 and u[2, 0] == 0.2
    p.write_text("internalField   uniform (1 0 0);\nboundaryField { }\n")
    assert foamio.internal_field(str(p), 4).shape == (4, 3)
    p.write_text("internalField   uniform 4e-07;\n")
    assert np.allclose(foamio.internal_field(str(p), 2), 4e-7)
    p.write_text("internalField   nonuniform List<scalar> \n2\n(\n1\n2\n3\n)\n;\n")
    with pytest.raises(ValueError):
        foamio.internal_field(str(p))


def test_residual_fields_by_equation():
    st = foamio.LogState(10**9)
    st.feed("Time = 1\nsmoothSolver:  Solving for Ux, Initial residual = 0.5, Final residual = 0.01, No Iterations 2\n"
            "GAMG:  Solving for p, Initial residual = 1, Final residual = 0.1, No Iterations 3\n"
            "Time = 2\nsmoothSolver:  Solving for Ux, Initial residual = 0.25, Final residual = 0.01, No Iterations 2\n"
            "GAMG:  Solving for p, Initial residual = 0.125, Final residual = 0.1, No Iterations 3\n")
    st.finish()
    rf = foamio.residual_fields(st.field_history)
    assert rf == {"it": [1, 2], "Ux": [0.5, 0.25], "p": [1.0, 0.125]}


def test_settling_iteration_of_a_damped_oscillation():
    sys.path.insert(0, os.path.join(os.path.dirname(HERE), "analysis"))
    import settle
    it = np.arange(1, 40001)
    ld = 78.0 * (1 + 0.02 * np.exp(-it / 6000.0) * np.cos(np.pi * it / 3500.0))
    s = settle.settled_from(it, ld, 2000, 1e-4)
    spread = settle.rolling_spread(ld, 2000) / ld
    assert s is not None and np.all(spread[s - 1:] <= 1e-4) and spread[s - 2] > 1e-4
    assert settle.settled_from(it, ld, 2000, 1e-9) is None
    assert settle.settled_from(it, ld, 2000, 1e-3) < s
