"""Test problems with exact ground truth, for validating the method.

`run(theta)` takes the numerical parameters as a dict and returns
`(observables, cost, wall, info)`: a dict of observables, a deterministic work
count, the wall-clock time in seconds, and a dict of diagnostics whose
`"converged"` is False if a Newton or GMRES solve missed its tolerance;
`exact_observables()` gives the exact values.  Two differ: `AdvectionDiffusion`
has `exact_Q()`, and its `run` returns `(Q, cost, wall)`; `make_synthetic(...)`
returns `(solver, Qstar)`, where `solver(theta)` returns `(Q, cost)`.
`errbudget.tune` calls a solver with keywords, so wrap `run` in a function that
returns `(observables, cost)`.
"""
from .burgers1d import Burgers1D
from .burgers2d import Burgers2D
from .advdiff import AdvectionDiffusion
from .synthetic import make_synthetic
from .reactdiff import ReactionDiffusion2D
from .mlhybrid import MLHybridBurgers1D
from .cavity import DrivenCavity

__all__ = ["Burgers1D", "Burgers2D", "AdvectionDiffusion", "make_synthetic", "ReactionDiffusion2D", "MLHybridBurgers1D", "DrivenCavity"]
