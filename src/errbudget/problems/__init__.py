"""Test problems with exact ground truth, for validating the method.

`run(theta)` takes the numerical parameters as a dict and returns
`(observables, cost, wall, info)`: a dict of observables, a deterministic work
count, the wall-clock time in seconds, and a dict of diagnostics whose
`"converged"` is False if a Newton or GMRES solve missed its tolerance;
`exact_observables()` gives the exact values.  `errbudget.tune` calls a solver
with keywords, so wrap `run` in a function that returns `(observables, cost)`.
"""
from .burgers1d import Burgers1D
from .burgers2d import Burgers2D
from .reactdiff import ReactionDiffusion2D
from .mlhybrid import MLHybridBurgers1D
from .cavity import DrivenCavity

__all__ = ["Burgers1D", "Burgers2D", "ReactionDiffusion2D", "MLHybridBurgers1D", "DrivenCavity"]
