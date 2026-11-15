"""Trace of the stepwise selection on the one V2 data set where the two implementations disagreed
(c_D 3e-3, start 10): the LOOCV criterion of every candidate in both implementations."""
import numpy as np, jax.numpy as jnp
import common as K, spre_white as SW
from sparse_pre.SPRE import SPRE
r, used = K.errbudget_search("cd", 3e-3, K.starts()[10])
X, Y, _ = K.pre_dataset(used, "cd")
nX = SW.EP + (X.max(0) - X.min(0)); nY = SW.EP + (Y.max() - Y.min())
Xn, Yn = X / nX, Y / nY
sp = SPRE("white", 3); sp.set_normalised_data(jnp.array(X), jnp.array(Y)); sp.prepare_jit_for_extrapolation_optimization()
def both(A):
    A = np.array(A, dtype=int)
    mine, a = SW.loocv(Xn, Yn, A)
    fit = sp.perform_extrapolation_optimization(jnp.array(A), do_jit=False)
    return mine, float(fit["cv"]), a, float(np.log1p(np.exp(fit["x"][0])) + 1e-16)
A = [[0, 0, 0]]
for order in range(1, 5):
    base = both(A)
    print(f"order {order} base A={A}: fast L={base[0]:.6f} authors L={base[1]:.6f} amp fast {base[2]:.3e} authors {base[3]:.3e}")
    cands = SW.next_candidates(np.array(A), order)
    keep = []
    for c in cands:
        An = A + [list(c)]
        if not SW.unisolvent(Xn, np.array(An)):
            print("   ", list(c), "not unisolvent"); continue
        m, au, am, aa = both(An)
        print(f"    + {list(c)}: fast L={m:.6f} authors L={au:.6f}  amp fast {am:.3e} authors {aa:.3e}  {'keep' if m < base[0] else ''}/{'keep' if au < base[1] else ''}")
        if m < base[0]: keep.append(list(c))
    if not keep: break
    A = A + keep
