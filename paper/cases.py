"""Registry of the development studies that the replays search: for each case the stored
runs (setting -> observables, cost, convergence information), the candidate values,
refinement scalars, declared classes and error models, the reference value and its
relative uncertainty for each observable, and the targets.

A target is graded (admitted) only if it is at least ten times the uncertainty of its
reference (Section 5.3).  converged(i) tells whether run i met its own convergence
criterion; the released configuration treats a run that did not as failed.
"""
from __future__ import annotations

import itertools
import os
import pickle
from dataclasses import dataclass, field

import numpy as np

import grids

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


@dataclass
class Case:
    key: str                     # short name used in tables
    label: str
    tensor: dict                 # idx -> (obs dict, cost, info dict)
    LAD: list
    HS: list
    KINDS: list
    exact: dict                  # observable -> exact or reference value
    targets: dict                # observable -> list of relative targets
    models: list | None = None
    unc: dict = field(default_factory=dict)   # observable -> relative uncertainty
    floor: float | None = None   # round-off floor measured by the probe (absolute)
    cost_note: str = ""

    def grid(self):
        return list(itertools.product(*[range(len(l)) for l in self.LAD]))

    def usable(self, i):
        return np.isfinite(self.tensor[i][1])

    def converged(self, i):
        return bool(self.tensor[i][2].get("converged", True))

    def admitted(self, ob, tg):
        return tg >= 10 * self.unc.get(ob, 0.0)


def _norm(t):
    """(obs, cost[, info]) with info a dict."""
    out = {}
    for i, v in t.items():
        info = v[2] if len(v) > 2 and isinstance(v[2], dict) else {}
        out[i] = (v[0], v[1], dict(info))
    return out


def _with_flags(name, t):
    path = os.path.join(D, f"flags_{name}.pkl")
    if not os.path.exists(path):
        return t
    fl = pickle.load(open(path, "rb"))
    for i, v in t.items():
        v[2].update({k: fl[i][k] for k in ("converged", "newton_fail", "gmres_fail",
                                           "newton_converged") if k in fl[i]})
    return t


def _manufactured():
    out = []
    for P in (grids.L2, grids.L3, grids.L4, grids.L5):
        t = _with_flags(P["name"], _norm(grids.load(P)))
        ex = grids.exact(P)
        tg = P["targets"] if isinstance(P["targets"], dict) else {o: P["targets"] for o in P["obs"]}
        if not P["rel"]:
            # L2 targets are absolute in the article; express them relative to
            # the exact value so that every case is graded the same way
            tg = {o: [x / abs(ex[o]) for x in v] for o, v in tg.items()}
        out.append(Case(P["name"], P["name"], t, P["lad"], grids.hs(P), P["kinds"], ex, tg,
                        cost_note="work units"))
    return out


def _dfg1(cost="deterministic"):
    from exp46_dfg_build import LAD
    from openfoam_cylinder import REF
    # forces from the exact wall-cell centroids and the second-order wall pressure;
    # dfg_tensor.pkl holds the earlier post-processing (vertex-mean cell centres)
    raw = pickle.load(open(os.path.join(D, os.environ.get("ERRBUDGET_DFG_TENSOR", "dfg_tensor_v12.pkl")), "rb"))
    t = {}
    for i, v in raw.items():
        info = dict(v[2]) if len(v) > 2 else {}
        c = v[1]
        if cost == "deterministic" and np.isfinite(c) and info.get("cells") and info.get("iters"):
            # cells x outer iterations: a deterministic work measure (CPU time of the
            # cheapest runs is a few hundredths of a second and mostly noise)
            c = float(info["cells"] * info["iters"])
        t[i] = (v[0], c, info)
    return Case("2D-1", "OpenFOAM DFG 2D-1", t, LAD,
                [[1.0 / n for n in LAD[0]], list(LAD[1]), list(LAD[2]), [1.0 / c for c in LAD[3]]],
                ["discretisation", "nested-solver", "nested-solver", "cost-only"], dict(REF),
                {"cd": [1e-2, 3e-3, 1e-3], "dp": [3e-2, 1e-2], "cl": [1e-1, 3e-2]},
                cost_note="cells x iterations" if cost == "deterministic" else "CPU s")


def _dfg2():
    z = pickle.load(open(os.path.join(D, "dfg2d2_tensor.pkl"), "rb"))
    LAD, ref = z["LAD"], z["ref"]
    t = {}
    for i, v in z["tensor"].items():
        info = dict(v[2])
        # PIMPLE stops at 50 outer iterations; a step that reaches it did not meet
        # residualControl
        info["converged"] = not (info.get("max_outer", 0) >= 50)
        t[i] = (v[0], v[1], info)
    R = ref[0][1]
    unc = {o: max(abs(ref[k][1][o] - R[o]) / abs(R[o]) for k in (1, 2)) for o in R}
    return Case("2D-2", "OpenFOAM DFG 2D-2", t, LAD,
                [[1.0 / n for n in LAD[0]], list(LAD[1]), list(LAD[2]), list(LAD[3])],
                ["discretisation", "discretisation", "nested-solver", "nested-solver"], dict(R),
                {"cd_max": [1e-1, 3e-2, 1e-2, 3e-3], "cl_max": [1e-1, 3e-2, 1e-2],
                 "st": [1e-1, 3e-2, 1e-2]}, unc=unc, cost_note="CPU s")


def _lammps():
    from lammps_madelung import LammpsMadelung
    out = []
    for n, ions in [(20, "8,000"), (40, "64,000"), (64, "262,144")]:
        t = _norm(pickle.load(open(os.path.join(D, f"lammps_tensor_n{n}.pkl"), "rb")))
        LAD = ([[6.0, 8.0, 10.0, 12.0, 14.0], [1e-2, 1e-3, 1e-4, 1e-5, 1e-6], [4, 5, 6, 7]]
               if n <= 20 else [[8.0, 10.0, 12.0, 14.0], [1e-3, 1e-4, 1e-5, 1e-6], [4, 5, 6, 7]])
        QS = LammpsMadelung(n_cells=n).exact_observables()["energy_per_ion"]
        out.append(Case(f"LMP-{n}", f"LAMMPS {ions} ions", t, LAD,
                        [[1.0 / c for c in LAD[0]], list(LAD[1]), [1.0 / o for o in LAD[2]]],
                        ["discretisation", "discretisation", "cost-only"],
                        {"energy_per_ion": QS}, {"energy_per_ion": [1e-4, 1e-5, 1e-6]},
                        cost_note="wall s"))
    return out


def _gromacs():
    LAD = [[0.9, 1.0, 1.1, 1.2, 1.4], [0.20, 0.16, 0.12, 0.10, 0.08], [4, 5, 6], [1e-4, 1e-5, 1e-6]]
    HS = [[1 / c for c in LAD[0]], list(LAD[1]), [1 / o for o in LAD[2]], list(LAD[3])]
    K = ["discretisation", "discretisation", "cost-only", "discretisation"]
    td, QD = pickle.load(open(os.path.join(D, "gromacs_tensor_gmx_d.pkl"), "rb"))
    tm, _ = pickle.load(open(os.path.join(D, "gromacs_tensor_gmx.pkl"), "rb"))
    # round-off floors measured by the probe of three converged settings
    # (exp44_gromacs_floor.py; values in figures/sec82_floor.txt)
    return [Case("GMX-d", "GROMACS double", _norm(td), LAD, HS, K, {"potential": QD},
                 {"potential": [1e-5, 1e-6]}, floor=2.85e-10 * abs(QD), cost_note="wall s"),
            Case("GMX-s", "GROMACS mixed", _norm(tm), LAD, HS, K, {"potential": QD},
                 {"potential": [1e-4, 3e-5, 1e-5]}, floor=0.4375, cost_note="wall s")]


def abinit_uncertainty(REF):
    """Relative uncertainty of the 48-Ha ABINIT reference (S-5.5).

    exp75_abinit_reference.py measures it with the recipe used for Quantum ESPRESSO and
    writes data/abinit_reference_v6.json; that value is used if the reference it measured
    is the one in the tensor.  Otherwise, and with ERRBUDGET_ABINIT_UNC=rev2, the value of
    the earlier 30-Ha reference (figures/sec83_reference.txt) is used."""
    unc = {"energy": 1.22e-4, "force": 4.45e-3}
    f = os.path.join(D, "abinit_reference_v6.json")
    if os.path.exists(f) and os.environ.get("ERRBUDGET_ABINIT_UNC", "measured") == "measured":
        import json
        m = json.load(open(f))
        same = all(abs(m["ref"][o] - REF[o]) <= 1e-8 * abs(REF[o]) for o in ("energy", "force"))
        if same:
            unc = {o: float(m["unc"][o]) for o in ("energy", "force")}
        else:
            print(f"cases.py: {f} measured a different reference ({m['ref']} against {dict(REF)}); "
                  "keeping the revision-2 uncertainty", flush=True)
    return unc


def _abinit():
    LAD = [[8.0, 12.0, 16.0, 20.0, 24.0, 28.0, 32.0], [2, 3, 4, 6],
           [1e-6, 1e-8, 1e-10, 1e-12], [2, 4, 6]]
    t, REF = pickle.load(open(os.path.join(D, "abinit_tensor_v11.pkl"), "rb"))
    return Case("ABINIT", "ABINIT Si", _norm(t), LAD,
                [list(LAD[0]), [1.0 / k for k in LAD[1]], list(LAD[2]), [1.0 / n for n in LAD[3]]],
                ["discretisation", "discretisation", "nested-solver", "cost-only"], dict(REF),
                {"energy": [1e-2, 3e-3, 1e-3], "force": [1e-1, 3e-2, 1e-2]},
                models=["exponential", "power", "power", "power"],
                unc=abinit_uncertainty(REF), cost_note="CPU s")


def _dfg2_refined():
    """The refined OpenFOAM 2D-2 study (Hamilton job J2, case 2D-2h): n_arc 24-96, time
    step, PIMPLE residualControl, GAMG tolerance, discarded transient R and averaging
    window W; forces from the forceCoeffs function object.  The reference is the run at
    n_arc = 128, dt = 1.25e-4 (the run at 192 diverged); hamilton/J2_dfg2d2/analyse_j2.py
    measures its uncertainty and writes the case specification."""
    p = os.path.normpath(os.path.join(os.path.dirname(D), "..", "hamilton", "results", "j2_tensor.pkl"))
    if not os.path.exists(p):
        raise FileNotFoundError(p)
    z = pickle.load(open(p, "rb"))
    if "case" not in z:
        raise FileNotFoundError(f"{p} was written by the revision-3 analyse_j2.py; re-run the "
                                "revision-4 analyse_j2.py to add the case specification")
    s = z["case"]
    t = {tuple(i): (v[0], v[1], dict(v[2])) for i, v in z["tensor"].items()}
    return Case("2D-2h", "OpenFOAM DFG 2D-2 (refined, current release)", t, s["LAD"], s["HS"],
                s["KINDS"], dict(s["exact"]), {k: list(v) for k, v in s["targets"].items()},
                models=list(s["models"]), unc=dict(s["unc"]), cost_note="core-seconds")


def _qe():
    LAD = [[12.0, 16.0, 20.0, 24.0, 28.0, 32.0, 40.0, 48.0], [2, 3, 4, 6],
           [1e-4, 1e-6, 1e-8, 1e-10, 1e-12], [0.3, 0.5, 0.7]]
    t, REF, UNC = pickle.load(open(os.path.join(D, "qe_tensor.pkl"), "rb"))
    return Case("QE", "Quantum ESPRESSO Si", _norm(t), LAD,
                [list(LAD[0]), [1.0 / k for k in LAD[1]], list(LAD[2]), [1.0 / b for b in LAD[3]]],
                ["discretisation", "discretisation", "nested-solver", "cost-only"], dict(REF),
                {"energy": [1e-3, 1e-4, 1e-5], "force": [1e-1, 3e-2, 1e-2]},
                models=["exponential", "power", "power", "power"], unc=dict(UNC),
                cost_note="CPU s")


def all_cases(dfg1_cost="deterministic", refined_2d2=True):
    """Every development case.  With refined_2d2 the refined 2D-2 study (case 2D-2h) is
    appended when hamilton/results/j2_tensor.pkl is present; the coarser 2D-2 study (case
    2D-2) stays in the list, and the tables show the refined one."""
    out = (_manufactured() + [_dfg1(dfg1_cost), _dfg2()] + _lammps() + _gromacs()
           + [_abinit(), _qe()])
    if refined_2d2:
        try:
            out.append(_dfg2_refined())
        except FileNotFoundError:
            pass
    return out
