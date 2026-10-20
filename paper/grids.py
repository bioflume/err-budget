"""Ladders, classes and tensors of the three test problems, shared by the scripts."""
import itertools
import os
import pickle

import numpy as np

D = os.path.join(os.path.dirname(__file__), "data")

L2 = dict(name="L2", file="burgers_tensor.pkl", rel=False,
          names=["N", "steps", "newton_tol", "gmres_tol"],
          lad=[[32, 64, 128, 256, 512], [25, 50, 100, 200, 400],
               [1e-2, 1e-4, 1e-6, 1e-8, 1e-10], [1e-2, 1e-4, 1e-6, 1e-8, 1e-10]],
          kinds=["discretisation", "discretisation", "nested-solver", "nested-solver"],
          obs=["energy", "dissipation"], targets=[1e-2, 3e-3, 1e-3])
L3 = dict(name="L3", file="burgers2d_tensor.pkl", rel=True,
          names=["N", "steps", "newton_tol", "gmres_tol", "ilu_drop"],
          lad=[[16, 32, 64, 128], [10, 20, 40, 80], [1e-2, 1e-5, 1e-8],
               [1e-2, 1e-5, 1e-8], [1e-1, 1e-2, 1e-3]],
          kinds=["discretisation", "discretisation", "nested-solver",
                 "nested-solver", "cost-only"],
          obs=["energy", "gradient"], targets=[1e-2, 3e-3])
L4 = dict(name="L4", file="reactdiff_tensor.pkl", rel=True,
          names=["N", "steps", "n_react", "gmres_tol", "ilu_drop"],
          lad=[[16, 32, 64, 128], [5, 10, 20, 40], [1, 2, 4, 8],
               [1e-2, 1e-5, 1e-8], [1e-1, 1e-2, 1e-3]],
          kinds=["discretisation", "discretisation", "discretisation",
                 "nested-solver", "cost-only"],
          obs=["inhibitor_mass", "pattern_sharpness"], targets=[1e-2, 3e-3])


L5 = dict(name="L5", file="cavity_tensor.pkl", rel=True,
          names=["N", "newton_tol", "gmres_tol", "ilu_drop"],
          lad=[[16, 24, 32, 48, 64, 96, 128], [1e-2, 1e-4, 1e-6, 1e-8, 1e-10],
               [1e-1, 1e-2, 1e-4, 1e-6, 1e-8], [1e-3, 1e-4, 1e-5, 1e-6]],
          kinds=["discretisation", "nested-solver", "nested-solver", "cost-only"],
          obs=["psi_min", "y_min", "omega_lid"],
          targets={"psi_min": [1e-2, 3e-3], "y_min": [1e-3, 3e-4],
                   "omega_lid": [1e-2, 3e-3]})


def hs(P):
    lad = P["lad"]
    if P["name"] == "L5":
        return [[1.0 / n for n in lad[0]]] + [list(l) for l in lad[1:]]
    out = [[2 * np.pi / n for n in lad[0]], [1.0 / s for s in lad[1]]]
    for k in range(2, len(lad)):
        out.append([1.0 / r for r in lad[k]] if P["names"][k] == "n_react"
                   else list(lad[k]))
    return out


def load(P):
    t = pickle.load(open(os.path.join(D, P["file"]), "rb"))
    return t[0] if isinstance(t, tuple) else t


def exact(P):
    from errbudget.problems import Burgers1D, Burgers2D, ReactionDiffusion2D
    from errbudget.problems import DrivenCavity
    p = {"L2": lambda: Burgers1D(nu=0.05, T=1.0), "L3": lambda: Burgers2D(nu=0.05, T=1.0),
         "L4": ReactionDiffusion2D, "L5": DrivenCavity}[P["name"]]()
    return p.exact_observables()


def grid(P):
    return list(itertools.product(*[range(len(l)) for l in P["lad"]]))
