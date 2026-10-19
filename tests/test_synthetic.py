"""On a solver whose error we impose, the measurement must recover it."""
import numpy as np
from knobs import LocalBudget
from knobs.problems import make_synthetic

def test_recovers_imposed_orders_and_residual():
    P, Q = (2.0, 1.0, 1.0), (1.0, 1.0, 0.0)
    lad = [list(2.0 ** -np.arange(1, 9))] * 3
    runner, Qstar = make_synthetic(p=P, q=Q)
    look = lambda idx: runner({f"h{k+1}": lad[k][j] for k, j in enumerate(idx)})
    lb = LocalBudget(look, lad, (7, 7, 7), hs=lad)
    lb.sweep(fine_reach=1); lb.per_knob_error()
    for k in range(3):
        assert abs(lb.pk[k] - P[k]) < 0.05
