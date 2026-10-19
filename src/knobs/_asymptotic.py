"""
knobs: observable-targeted error budgets for black-box scientific solvers.

Structural model
----------------
A solver depends on numerical knobs theta = (theta_1, ..., theta_p).  Each knob
has an ordered ladder of settings (coarse -> fine) and a refinement scalar h_k
that decreases monotonically toward the fully-refined limit.

We assume, and then *test*:

    error is additive across knobs      E(theta)  ~ sum_k e_k(theta_k)
    cost is multiplicative across knobs log C(theta) ~ sum_k g_k(theta_k)

This duality is what makes the problem tractable: both structures are
identified by a one-at-a-time (OAT) design costing O(sum_k m_k) runs rather
than the O(prod_k m_k) runs of a tensor sweep.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from typing import Any, Callable, Sequence

import numpy as np


@dataclass
class Knob:
    """One numerical parameter, with an ordered ladder from coarse to fine."""

    name: str
    settings: Sequence[Any]
    h: Sequence[float]          # refinement scalar; strictly decreasing
    cost_law: str = "power"     # "power" (cost ~ h^-q) or "log" (cost ~ log(1/h))
    label: str = ""

    def __post_init__(self):
        self.h = np.asarray(self.h, dtype=float)
        if not np.all(np.diff(self.h) < 0):
            raise ValueError(f"knob {self.name!r}: h must be strictly decreasing")
        if len(self.h) != len(self.settings):
            raise ValueError(f"knob {self.name!r}: h and settings length mismatch")
        if not self.label:
            self.label = self.name

    @property
    def n(self) -> int:
        return len(self.settings)


class BudgetStudy:
    """Measure, model and optimize an observable's numerical error budget."""

    def __init__(self, runner: Callable[[dict], tuple], knobs: Sequence[Knob]):
        self.runner = runner
        self.knobs = list(knobs)
        self.cache: dict[tuple, tuple] = {}
        self.e: list[np.ndarray] | None = None
        self.g: list[np.ndarray] | None = None

    # ------------------------------------------------------------------ runs
    def theta(self, idx) -> dict:
        return {k.name: k.settings[i] for k, i in zip(self.knobs, idx)}

    def run(self, idx) -> tuple[float, float]:
        idx = tuple(int(i) for i in idx)
        if idx not in self.cache:
            self.cache[idx] = self.runner(self.theta(idx))
        return self.cache[idx]

    @property
    def ref_idx(self) -> tuple:
        """Reference configuration: every knob at its finest available setting."""
        return tuple(k.n - 1 for k in self.knobs)

    @property
    def n_runs(self) -> int:
        return len(self.cache)

    # -------------------------------------------------------------- measuring
    def measure_axes(self):
        """OAT relaxation away from the reference. O(sum_k m_k) runs."""
        ref = self.ref_idx
        self.Qref, self.Cref = self.run(ref)
        self.e, self.g = [], []
        for kk, knob in enumerate(self.knobs):
            e_k = np.zeros(knob.n)
            g_k = np.zeros(knob.n)
            for j in range(knob.n):
                idx = list(ref)
                idx[kk] = j
                Q, C = self.run(idx)
                e_k[j] = Q - self.Qref          # signed: sign carries cancellation
                g_k[j] = np.log(C) - np.log(self.Cref)
            self.e.append(e_k)
            self.g.append(g_k)
        self.n_runs_axes = self.n_runs
        return self

    # ---------------------------------------------------------------- models
    def predict_error(self, idx, target: str = "reference") -> float:
        """Signed error, relative to the reference config or to the true value."""
        v = float(sum(self.e[k][i] for k, i in enumerate(idx)))
        if target == "truth":
            v += getattr(self, "ref_residual", 0.0)
        return v

    def predict_abs_error(self, idx, target: str = "reference") -> float:
        """Triangle-inequality (cancellation-free) error bound."""
        v = float(sum(abs(self.e[k][i]) for k, i in enumerate(idx)))
        if target == "truth":
            v += abs(getattr(self, "ref_residual", 0.0))
        return v

    def predict_logcost(self, idx) -> float:
        return float(np.log(self.Cref) + sum(self.g[k][i] for k, i in enumerate(idx)))

    def predict_cost(self, idx) -> float:
        return float(np.exp(self.predict_logcost(idx)))

    # ------------------------------------------------------------- power laws
    def fit_power_laws(self, p_grid=None) -> list[dict]:
        """Fit the OAT curves and recover the reference's own residual error.

        Critical subtlety: the measured OAT curve is *not* the error.  It is the
        error minus the reference's error along that axis,

            e_k(h) = c_k h^p_k - c_k h_ref^p_k ,

        so a naive log-log fit of |e_k| against h is biased whenever h_ref is not
        negligible -- and it is never negligible near the reference, which is
        exactly where the curve is being sampled.  Fitting the offset form
        instead both de-biases p_k and yields r_k = c_k h_ref^p_k, the error the
        reference configuration still carries.  Summing r_k over knobs converts
        errors-relative-to-reference into errors-relative-to-truth without ever
        knowing the truth.
        """
        if p_grid is None:
            p_grid = np.linspace(0.25, 6.0, 231)
        fits = []
        for kk, knob in enumerate(self.knobs):
            h, e, href = knob.h, self.e[kk], knob.h[-1]
            best = (np.inf, np.nan, np.nan)
            for p in p_grid:
                basis = h**p - href**p
                denom = basis @ basis
                if denom <= 0:
                    continue
                c = (basis @ e) / denom
                res = np.sum((e - c * basis) ** 2)
                if res < best[0]:
                    best = (res, p, c)
            _, p, c = best
            rel_res = np.sqrt(best[0] / max(np.sum(e**2), 1e-300))
            gk = self.g[kk]
            gd = np.abs(gk) > 1e-12
            q = -np.polyfit(np.log(h[gd]), gk[gd], 1)[0] if gd.sum() >= 2 else np.nan
            fits.append(dict(name=knob.name, p=p, c=c, q=q,
                             ref_residual=c * href**p, fit_rel_residual=rel_res))
        self.fits = fits
        self.ref_residual = float(sum(f["ref_residual"] for f in fits))
        return fits

    # ------------------------------------------------------------- regimes
    def classify(self, kappa: float = 1e-2) -> list[dict]:
        """Classify each knob's OAT curve as converging or masked.

        The OAT curve e_k(h) is measured with every *other* knob at its finest
        setting, and it vanishes at the reference by construction -- so it can
        never show saturation as a nonzero floor.  What it shows instead is a
        *plateau at zero*: a run of settings over which refining knob k changes
        the observable not at all.  The coarsest setting at which that plateau
        begins is the knob's masking threshold.

        This threshold is conservative.  Masking can only strengthen as the
        other knobs coarsen (their errors grow, so knob k's contribution is
        buried deeper), hence a threshold measured with everything else at its
        finest is a setting the user can always safely adopt.  `verify_
        conservative` tests that claim against the tensor sweep.
        """
        out = []
        for kk, knob in enumerate(self.knobs):
            e = np.abs(self.e[kk])
            emax = e.max()
            tau = kappa * emax
            above = np.where(e > tau)[0]
            j_sat = int(above[-1]) if len(above) else 0
            thr = min(j_sat + 1, knob.n - 1)
            masked = thr <= knob.n - 2          # plateau covers >1 setting
            out.append(dict(name=knob.name, threshold_idx=thr,
                            threshold=knob.settings[thr], masked=masked,
                            n_wasted_steps=knob.n - 1 - thr, emax=float(emax)))
        self.regimes = out
        return out

    def refresh_residual(self):
        """Reference residual, counting only knobs that are still converging."""
        reg = {r["name"]: r for r in getattr(self, "regimes", [])}
        tot = 0.0
        for f in self.fits:
            if not reg or not reg[f["name"]]["masked"]:
                tot += f["ref_residual"]
        self.ref_residual = float(tot)
        return tot

    def advice(self, user_idx, eps: float, mode: str = "conservative",
               target: str = "reference") -> dict:
        """What to tell a user sitting at user_idx who wants accuracy eps."""
        rec, c_rec = self.optimize(eps, mode=mode, target=target)
        c_user = self.predict_cost(user_idx)
        shares = self.budget_shares(user_idx)
        dom = max(shares, key=shares.get) if shares else None
        lines = []
        for kk, knob in enumerate(self.knobs):
            r = self.regimes[kk]
            over = user_idx[kk] - r["threshold_idx"]
            if over > 0:
                lines.append(f"{knob.name}: {over} ladder step(s) past its masking "
                             f"threshold ({knob.settings[r['threshold_idx']]}) "
                             f"- refinement here is not reaching the observable")
        return dict(recommended=rec, cost_user=c_user, cost_rec=c_rec,
                    overpayment=c_user / c_rec if rec is not None else np.inf,
                    dominant=dom, shares=shares, notes=lines)

    # ----------------------------------------------------------- separability
    def check_separability(self, n_samples: int = 40, rng=None) -> dict:
        """Compare the additive prediction against off-axis runs."""
        rng = np.random.default_rng(0 if rng is None else rng)
        pred, actual, costs_p, costs_a = [], [], [], []
        seen = set()
        for _ in range(n_samples):
            idx = tuple(int(rng.integers(0, k.n)) for k in self.knobs)
            if idx in seen or sum(i != r for i, r in zip(idx, self.ref_idx)) < 2:
                continue
            seen.add(idx)
            Q, C = self.run(idx)
            actual.append(Q - self.Qref)
            pred.append(self.predict_error(idx))
            costs_a.append(np.log(C))
            costs_p.append(self.predict_logcost(idx))
        pred, actual = np.asarray(pred), np.asarray(actual)
        scale = np.max(np.abs(actual)) if len(actual) else 1.0
        return dict(
            n=len(actual),
            pred=pred,
            actual=actual,
            interaction_rel=float(np.max(np.abs(pred - actual)) / (scale + 1e-300)),
            interaction_rms=float(np.sqrt(np.mean((pred - actual) ** 2)) / (scale + 1e-300)),
            logcost_max_err=float(np.max(np.abs(np.asarray(costs_p) - np.asarray(costs_a)))),
        )

    # -------------------------------------------------------------- analytics
    def budget_shares(self, idx) -> dict:
        contrib = {k.name: abs(self.e[i][j]) for i, (k, j) in enumerate(zip(self.knobs, idx))}
        tot = sum(contrib.values())
        if tot == 0:
            return {k: 0.0 for k in contrib}
        return {k: v / tot for k, v in contrib.items()}

    def all_idx(self):
        return itertools.product(*[range(k.n) for k in self.knobs])

    def _errfun(self, mode: str, target: str):
        if mode == "conservative":
            return lambda i: self.predict_abs_error(i, target)
        return lambda i: abs(self.predict_error(i, target))

    def optimize(self, eps: float, mode: str = "conservative", target: str = "reference"):
        """Cheapest configuration whose modelled error is within eps.

        mode="conservative": sum_k |e_k| <= eps   (no reliance on cancellation)
        mode="optimistic":  |sum_k e_k| <= eps    (exploits cancellation; fragile)
        """
        err = self._errfun(mode, target)
        best, best_cost = None, np.inf
        for idx in self.all_idx():
            if err(idx) <= eps:
                c = self.predict_logcost(idx)
                if c < best_cost:
                    best, best_cost = idx, c
        return best, (np.exp(best_cost) if best is not None else np.inf)

    def equimarginal(self, eps: float, mode: str = "conservative",
                     target: str = "reference"):
        """Greedy refinement by error-reduction per unit log-cost.

        This is the algorithmic form of the equimarginal principle: refine
        whichever knob currently buys the most accuracy per unit of extra cost.
        """
        err = self._errfun(mode, target)
        idx = [0] * len(self.knobs)
        while err(idx) > eps:
            best_gain, best_k = -np.inf, None
            for kk, knob in enumerate(self.knobs):
                if idx[kk] == knob.n - 1:
                    continue
                trial = list(idx)
                trial[kk] += 1
                d_err = err(idx) - err(trial)
                d_cost = self.g[kk][trial[kk]] - self.g[kk][idx[kk]]
                gain = d_err / max(d_cost, 1e-12)
                if gain > best_gain:
                    best_gain, best_k = gain, kk
            if best_k is None:
                break
            idx[best_k] += 1
        return tuple(idx), self.predict_cost(idx)

    def overpayment(self, idx, mode: str = "conservative",
                    target: str = "reference") -> float:
        """How many times more expensive than necessary, at the same accuracy."""
        eps = self._errfun(mode, target)(idx)
        _, c_opt = self.optimize(eps, mode=mode, target=target)
        return self.predict_cost(idx) / c_opt


def balance_law(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Optimal ratio of error contributions across knobs: e_k proportional to q_k/p_k.

    Derived from min log C s.t. sum_k e_k = eps with e_k = c_k h_k^p_k and
    log C = const - sum_k q_k log h_k.  The folklore rule "balance all error
    contributions equally" is the special case q_k/p_k = const.
    """
    r = np.asarray(q, float) / np.asarray(p, float)
    return r / r.sum()
