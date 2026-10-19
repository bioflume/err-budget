"""Test problems with exact ground truth, for validating the method.

Each exposes `exact_observables()` and a `run(theta)` returning observables and
a cost, so any of them can be driven by `knobs.tune`.
"""
from .burgers1d import Burgers1D
from .burgers2d import Burgers2D
from .advdiff import AdvectionDiffusion
from .synthetic import make_synthetic
from .reactdiff import ReactionDiffusion2D
from .mlhybrid import MLHybridBurgers1D
from .cavity import DrivenCavity

__all__ = ["Burgers1D", "Burgers2D", "AdvectionDiffusion", "make_synthetic", "ReactionDiffusion2D", "MLHybridBurgers1D", "DrivenCavity"]
