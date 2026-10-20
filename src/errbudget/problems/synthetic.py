"""L0: a synthetic solver whose error and cost are imposed analytically."""
from __future__ import annotations
import numpy as np


def make_synthetic(c=(1.0, -0.6, 0.25), p=(2.0, 1.0, 1.0), q=(1.0, 1.0, 0.0),
                   interaction=0.0, noise=0.0, seed=0):
    """A 'solver' whose error structure we impose analytically.

    Q(h) = Qstar + sum_k c_k h_k^p_k + interaction * h_1^p_1 h_2^p_2 + noise
    cost = prod_k h_k^-q_k   (q_k = 0 means a log-cost parameter)
    """
    rng = np.random.default_rng(seed)
    Qstar = 1.0

    def runner(theta):
        h = np.array([theta["h1"], theta["h2"], theta["h3"]])
        e = sum(c[k] * h[k] ** p[k] for k in range(3))
        e += interaction * h[0] ** p[0] * h[1] ** p[1]
        if noise:
            e += noise * rng.standard_normal()
        cost = 1.0
        for k in range(3):
            cost *= h[k] ** (-q[k]) if q[k] > 0 else (1.0 + 0.15 * np.log(1.0 / h[k]))
        return Qstar + e, cost

    return runner, Qstar


