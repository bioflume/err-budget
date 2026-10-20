"""A solver with a learned component, to test what happens to the certificate.

Machine-learning-augmented solvers -- neural closures, learned subgrid models,
surrogate operators, learned preconditioners used as coarse propagators --
introduce a parameter that does not behave like any classical one.  Training a
larger model, or training it longer, reduces its error only down to an
irreducible floor set by model bias and by the data.  There is no ladder along
which that error tends to zero.

Here the viscous term of the Burgers problem is supplied by a surrogate whose
relative error follows

    delta(level) = delta_floor + A * 2^(-level),

so refining the surrogate improves the answer until `delta_floor` and then stops
permanently.  In the taxonomy of the method this is a discretisation-like
parameter with a hard plateau, and the question is whether the method's
certificate survives: it should refuse to certify any target below the floor
rather than refining the classical parameters indefinitely against an error they
cannot reach.
"""
from __future__ import annotations

import numpy as np

from .burgers1d import Burgers1D


class MLHybridBurgers1D(Burgers1D):
    """Burgers with the diffusion coefficient supplied by a surrogate model."""

    def __init__(self, nu=0.05, T=1.0, max_newton=25,
                 delta_floor=2.0e-3, amplitude=6.0e-2):
        super().__init__(nu=nu, T=T, max_newton=max_newton)
        self.delta_floor, self.amplitude = delta_floor, amplitude

    def surrogate_error(self, level: int) -> float:
        """Relative error of the learned viscosity at this training level."""
        return self.delta_floor + self.amplitude * 2.0 ** (-float(level))

    def run(self, theta: dict):
        level = int(theta.get("model_level", 8))
        delta = self.surrogate_error(level)
        nu_true = self.nu
        try:
            # the surrogate reports a biased viscosity; everything else is exact
            self.nu = nu_true * (1.0 + delta)
            obs, cost, wall, info = super().run(
                {k: v for k, v in theta.items() if k != "model_level"})
        finally:
            self.nu = nu_true
        # a larger surrogate costs more to evaluate, as a network would
        cost = cost * (1.0 + 0.15 * level)
        info = dict(info)
        info["surrogate_error"] = delta
        return obs, cost, wall, info

    def exact_observables(self, Nf: int = 1 << 15) -> dict:
        """Truth is the true physics, not the surrogate's version of it."""
        return super().exact_observables(Nf)
