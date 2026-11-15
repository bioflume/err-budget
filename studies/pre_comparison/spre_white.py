"""Sparse probabilistic Richardson extrapolation (SPRE) with the white-noise kernel, in NumPy.

This is the model and the algorithm of Oates, Howey and Karvonen, "Sparse probabilistic
Richardson extrapolation" (arXiv 2604.02072), Sections 3.1 and 3.3, as implemented in the
authors' code (github.com/NewcastleRSE/sparse-probabilistic-richardson-extrapolation,
commit 3918107, src/sparse_pre/SPRE.py): universal kriging with mean sum_{a in A} beta_a x^a
(flat prior on beta), covariance sigma^2 delta(x, x'), the amplitude sigma^2 chosen by
leave-one-out cross-validation (LOOCV), and the index set A chosen by the authors' stepwise
procedure.  Inputs and outputs are normalised by their ranges, as in the authors' code.

With the white-noise kernel every quantity has a closed form, which is what makes 4,256
fits affordable (the authors' JAX code recompiles for every data size):

  * the posterior mean at 0 is the intercept of the least-squares fit on the basis A;
  * the leave-one-out mean of point i is y_i - e_i / (1 - h_ii) and its variance is
    sigma^2 / (1 - h_ii) (e: least-squares residuals, h_ii: leverages);
  * the LOOCV criterion is minimised over sigma^2 at sigma^2 = mean(e_i^2 / (1 - h_ii)), where it equals
    n/2 log(2 pi sigma^2) - 1/2 sum log(1 - h_ii) + n/2;
  * the posterior variance at 0 is sigma^2 (1 + e_1^T (V^T V)^{-1} e_1), the first term being
    the kernel's own value k(0, 0) (the authors' code includes it).

validate_spre.py checks this file against the authors' code.
"""
from __future__ import annotations

import numpy as np

EP = 1e-16          # smallest normalising constant and amplitude, as in sparse_pre.SPRE


def vandermonde(Xn: np.ndarray, A: np.ndarray) -> np.ndarray:
    """V[i, j] = prod_k Xn[i, k] ** A[j, k] (0 ** 0 = 1)."""
    return np.prod(Xn[:, None, :] ** A[None, :, :], axis=2)


def next_candidates(A: np.ndarray, order: int) -> np.ndarray:
    """sparse_pre.helper_functions.stepwise: raise each row of degree order - 1 by one in
    every coordinate; unique rows in lexicographic order."""
    d = A.shape[1]
    rows = A[A.sum(axis=1) == order - 1]
    if len(rows) == 0:
        return np.empty((0, d), dtype=int)
    out = (rows[:, None, :] + np.eye(d, dtype=int)[None, :, :]).reshape(-1, d)
    return np.unique(out, axis=0)


def loocv(Xn: np.ndarray, Yn: np.ndarray, A: np.ndarray):
    """Minimal LOOCV criterion L(A) (negative log predictive density, normalised data) and the
    amplitude that attains it.  Returns (inf, nan) where a leave-one-out fit is not defined."""
    n = len(Yn)
    V = vandermonde(Xn, A)
    if n < A.shape[0] + 1 or np.linalg.matrix_rank(V) < A.shape[0]:
        return np.inf, np.nan
    Q, R = np.linalg.qr(V)
    hat = np.sum(Q ** 2, axis=1)
    if np.any(hat > 1.0 - 1e-10):
        return np.inf, np.nan                       # leaving this point out loses rank
    beta = np.linalg.solve(R, Q.T @ Yn)
    e = Yn - V @ beta
    a = float(np.mean(e ** 2 / (1.0 - hat)))
    a = max(a, EP)
    L = 0.5 * n * np.log(2 * np.pi * a) - 0.5 * np.sum(np.log(1.0 - hat)) \
        + 0.5 * float(np.sum(e ** 2 / (1.0 - hat))) / a
    return float(L), a


def unisolvent(Xn: np.ndarray, A: np.ndarray) -> bool:
    return np.linalg.matrix_rank(vandermonde(Xn, A)) == A.shape[0]


def fit(X: np.ndarray, Y: np.ndarray, max_order: int = 0, A_fixed: np.ndarray | None = None) -> dict:
    """SPRE with the white-noise kernel and stepwise selection of A (SPRE.stepwise_selection).

    X: (n, d) distinct inputs (x -> 0 is the limit); Y: (n,) outputs.  Returns the posterior
    mean and variance of f(0), the index set A, the amplitude (normalised units) and L(A)."""
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    n, d = X.shape
    nX = EP + (X.max(axis=0) - X.min(axis=0))
    nY = EP + (Y.max() - Y.min())
    Xn, Yn = X / nX, Y / nY
    A = np.zeros((1, d), dtype=int)
    if A_fixed is not None:
        A = np.asarray(A_fixed, dtype=int)
        cv, a = loocv(Xn, Yn, A)
    else:
        cv, a = loocv(Xn, Yn, A)
        order, carry_on = 0, True
        while carry_on and (max_order == 0 or order < max_order):
            m = A.shape[0]
            order += 1
            extra = next_candidates(A, order)
            keep = np.zeros(len(extra), dtype=bool)
            for i, row in enumerate(extra):
                A_new = np.vstack([A, row])
                if unisolvent(Xn, A_new):
                    cv_new, _ = loocv(Xn, Yn, A_new)
                    if cv_new < cv:
                        keep[i] = True
            if keep.any() and (m + keep.sum()) < (n - 1):
                A_up = np.vstack([A, extra[keep]])
                cv_up, a_up = loocv(Xn, Yn, A_up)
                if cv_up >= cv:
                    carry_on = False
                else:
                    A, cv, a = A_up, cv_up, a_up
            else:
                carry_on = False
    if not np.isfinite(cv):
        return dict(mu=np.nan, var=np.nan, sd=np.nan, A=A, amp=np.nan, cv=cv, n=n)
    V = vandermonde(Xn, A)
    beta, *_ = np.linalg.lstsq(V, Yn, rcond=None)
    G = np.linalg.inv(V.T @ V)
    var_n = a * (1.0 + G[0, 0])
    return dict(mu=float(nY * beta[0]), var=float(nY ** 2 * var_n), sd=float(nY * np.sqrt(var_n)),
                A=A, amp=a, cv=cv, n=n)


def design_variance(Xn_all: np.ndarray, A: np.ndarray) -> float:
    """Design criterion of SPRE (eq. 11 of the paper) up to the fixed amplitude: the posterior
    variance of f(0) for the design Xn_all, which does not depend on the outputs.  For the
    white-noise kernel it is a (1 + e_1^T (V^T V)^{-1} e_1), so only the bracket is returned.
    Infinite if the design is not P_A-unisolvent."""
    V = vandermonde(Xn_all, A)
    if V.shape[0] < V.shape[1] or np.linalg.matrix_rank(V) < V.shape[1]:
        return np.inf
    G = np.linalg.inv(V.T @ V)
    return float(1.0 + G[0, 0])
