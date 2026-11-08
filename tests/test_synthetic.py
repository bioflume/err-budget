"""On a solver whose error we impose, the measurement must recover it."""
import numpy as np
from errbudget import LocalBudget


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


def test_recovers_imposed_orders_and_residual():
    P, Q = (2.0, 1.0, 1.0), (1.0, 1.0, 0.0)
    lad = [list(2.0 ** -np.arange(1, 9))] * 3
    runner, Qstar = make_synthetic(p=P, q=Q)
    look = lambda idx: runner({f"h{k+1}": lad[k][j] for k, j in enumerate(idx)})
    lb = LocalBudget(look, lad, (7, 7, 7), hs=lad)
    lb.sweep(fine_reach=1); lb.per_knob_error()
    for k in range(3):
        assert abs(lb.pk[k] - P[k]) < 0.05
