"""Every study of the article as one registry, for the pooled results table.

Each case carries what the replay needs: the exhaustive tensor (configuration
-> observables, cost, per-run information), the candidate values, refinement
scalars, declared classes and error models, the exact or reference value of
each observable with its uncertainty, the graded targets, and a flag saying
whether each run met its own convergence criterion.

The declarations are those a user following the article would make; they are
the same for every configuration of the program that is compared.  Targets are
graded only if they are at least ten times the relative uncertainty of the
reference (`admitted`).
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
    for P, tg_abs in [(grids.L2, None), (grids.L3, None), (grids.L4, None), (grids.L5, None)]:
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
    # version-1.2 post-processing (exact cell centroids at the wall, second-order
    # wall pressure); the tensor of revision 2 is dfg_tensor.pkl
    raw = pickle.load(open(os.path.join(D, os.environ.get("KNOBS_DFG_TENSOR", "dfg_tensor_v12.pkl")), "rb"))
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
    """Relative uncertainty of the ABINIT reference used for grading.

    Until revision 4 this was the value measured on the 30-Ha reference of revision 2
    (figures/sec83_reference.txt: 1.22e-4 in the energy, 4.45e-3 in the force), although
    the tensor is graded against the 48-Ha reference (abinit_silicon.py, reference()).
    paper/exp75_abinit_reference.py measures the uncertainty of the 48-Ha reference with
    the recipe used for Quantum ESPRESSO and writes data/abinit_reference_v6.json; that
    value is used if the reference it measured is the one in the tensor."""
    unc = {"energy": 1.22e-4, "force": 4.45e-3}
    f = os.path.join(D, "abinit_reference_v6.json")
    if os.path.exists(f) and os.environ.get("KNOBS_ABINIT_UNC", "measured") == "measured":
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
    """Revision 4: DFG 2D-2 refined on Hamilton (job J2).  Parameters n_arc 24-96, time
    step, PIMPLE residualControl, GAMG tolerance, relaxation time and averaging window;
    forces from OpenFOAM's forceCoeffs function object; reference n_arc 192, dt 1.25e-4,
    whose uncertainty is measured against two other reference runs and the window
    variants (hamilton/J2_dfg2d2/analyse_j2.py, which writes the case specification)."""
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
    """Every development case.  With refined_2d2 (revision 4) the refined 2D-2 study of
    the Hamilton job J2 is appended as case "2D-2h" when its results are present; the
    revision-3 2D-2 case stays in the list, and the tables choose which one to show."""
    out = (_manufactured() + [_dfg1(dfg1_cost), _dfg2()] + _lammps() + _gromacs()
           + [_abinit(), _qe()])
    if refined_2d2:
        try:
            out.append(_dfg2_refined())
        except FileNotFoundError:
            pass
    return out
