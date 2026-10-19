"""Local error budget: characterise the solver where the user actually runs it.

A fine-reference design fails for nested solver tolerances, because a knob's own
error contribution is not invariant to the other knobs: coarsening the time step
makes every Newton/GMRES solve harder, so the tolerance-induced error grows just
as the discretisation error that would have buried it grows.  A threshold
measured with everything else at its finest is therefore not transferable.  So we
anchor locally, and assume nothing about the solver anywhere we did not sample.

One anchor is enough.  An earlier design added a second -- a jointly refined
Richardson path -- to supply the absolute error.  Once the per-knob fit was
sharpened (`per_knob_error`), the path became redundant: measured over 10
(problem, observable, target) combinations in 1D and 2D, dropping it left safety
at 100% and penalties unchanged or better, for 1.8-3.3x less calibration cost.
It had also been the source of the worst over-refinement we saw (5.68x), because
taking the larger of two estimators inherits the weaker one's pessimism.  It is
kept behind `path_levels` for the ablation, and off by default.

Everything now comes from one asymmetric sweep around the anchor: sensitivities,
marginal costs, masking thresholds, the per-knob error budget, and the absolute
error estimate.
"""
from __future__ import annotations

import itertools

import numpy as np


class LocalBudget:
    def __init__(self, lookup, ladders, anchor, hs=None, fs_cap=None,
                 resid_margin=0.0, models=None, unfit="order1", kinds=None,
                 local_order=False, resolution=0.0, resid_mode="relative",
                 cancel_rule="credit", joint_res=3.0):
        """lookup(idx) -> (Q, cost);  ladders: settings;  hs: refinement scalars.

        models[k] is "power" (default: e_k(h) = c (h^p - h_a^p)) or
        "exponential" (e_k(x) = c (exp(-b x) - exp(-b x_a)), with x the setting
        itself, e.g. a plane-wave cutoff, for which algebraic convergence in a
        refinement scalar is the wrong model).
        """
        self.models = list(models) if models is not None else None
        self.unmeasured_rule = "fit"      # "max": also charge the last step
        # An axis with fewer than three successful runs cannot be fitted (a
        # start at the coarse end of a ladder sees only one finer value).
        # "floor": it enters only through the consistency floor (v1.0);
        # "add": its measured one-step change is added to the budget;
        # "order1" (default): that change is extrapolated assuming the axis
        # converges at least to first order in its refinement scalar, i.e.
        # multiplied by h_finer / (h - h_finer) + 1.  Without it, first-order
        # observables on closely spaced ladders were certified up to 2.7 eps
        # (OpenFOAM, exp46_dfg).
        self.unfit = unfit
        # Version 1.2.  kinds: declared classes (the local-order check and the
        # joint check act on discretisation parameters only).  local_order: never
        # claim less remaining error along a discretisation axis than the
        # three-point extrapolation through the anchor's neighbours with their
        # observed order (Section 3.2).  resolution: differences smaller than
        # this (e.g. the measured round-off floor) are treated as zero.
        self.kinds = list(kinds) if kinds is not None else None
        self.local_order = bool(local_order)
        self.resolution = float(resolution or 0.0)
        self.resid_mode = resid_mode
        self.cancel_rule = cancel_rule
        # a parameter takes part in the joint check if its one-step change exceeds
        # joint_res times the resolution (3 in the frozen version 1.2, 1 from 1.2.1)
        self.joint_res = float(joint_res)
        self.cancelled = False
        self.local = {}
        self.joint = None
        self.lookup = lookup
        self.ladders = [list(l) for l in ladders]
        self.m = [len(l) for l in self.ladders]
        self.hs = [np.asarray(h, float) for h in hs] if hs is not None else None
        self.anchor = tuple(anchor)
        self.fs_cap = fs_cap
        self.resid_margin = float(resid_margin)
        self.runs_used = 0
        self.protocol_cost = 0.0

    def _get(self, idx):
        idx = tuple(idx)
        cb = getattr(self, "on_charge", None)
        if cb is not None:
            cb(idx, getattr(self, "_tag", "?"))
        Q, C = self.lookup(idx)
        self.runs_used += 1
        self.protocol_cost += C
        return Q, C

    # ------------------------------------------------------------- anchor 1
    def sweep(self, fine_reach: int = 1, active=None):
        """OAT sweep, asymmetric by design.

        Coarsening is nearly free, so we sample the whole coarse side.  Refining
        is what costs, so we step only `fine_reach` settings finer than the
        anchor and let the refinement path carry the rest.  Sampling the full
        ladder instead makes calibration cost hundreds of times the production
        run -- which is why an unconstrained sweep is not a usable protocol.
        Settings beyond the sampled window are extrapolated log-linearly.
        """
        self._tag = "anchor"
        self.Q0, self.C0 = self._get(self.anchor)
        self.e, self.g = [], []
        self.sampled = []
        self.fine_failed = []
        act = set(range(len(self.m))) if active is None else set(active)
        for k in range(len(self.m)):
            e_k, g_k = np.zeros(self.m[k]), np.zeros(self.m[k])
            if k not in act:
                # Frozen by screening: never probed, never moved.  Zero
                # sensitivity and zero cost delta keep it out of both phases.
                self.e.append(e_k); self.g.append(g_k)
                self.sampled.append([self.anchor[k]])
                self.fine_failed.append(False)
                continue
            reach = fine_reach
            if (self.local_order and self.anchor[k] == 0 and self._disc(k)):
                # no coarse side: a second finer value gives the three points
                # the local-order check needs
                reach = max(reach, 2)
            hi = min(self.anchor[k] + reach, self.m[k] - 1)
            js = [j for j in range(0, hi + 1)]
            failed_fine = False
            for j in js:
                if j == self.anchor[k]:
                    continue
                idx = list(self.anchor)
                idx[k] = j
                self._tag = "sweep_fine" if j > self.anchor[k] else "sweep_coarse"
                Q, C = self._get(idx)
                if not (np.isfinite(Q) and np.isfinite(C)):
                    # A failed run (out of memory, crash) carries no information.
                    e_k[j], g_k[j] = np.nan, np.inf
                    if j > self.anchor[k]:
                        failed_fine = True
                    continue
                e_k[j] = Q - self.Q0
                g_k[j] = np.log(C) - np.log(self.C0)
            # extrapolate the unsampled fine tail
            if hi < self.m[k] - 1 and hi >= 1:
                de = e_k[hi] - e_k[hi - 1]
                dg = g_k[hi] - g_k[hi - 1]
                for j in range(hi + 1, self.m[k]):
                    e_k[j] = e_k[hi] + de * (j - hi) * 0.5   # damped: errors flatten
                    g_k[j] = g_k[hi] + dg * (j - hi)         # cost keeps growing
            self.e.append(e_k)
            self.g.append(g_k)
            self.sampled.append(js)
            self.fine_failed.append(failed_fine)
        return self

    def _disc(self, k):
        return self.kinds is None or self.kinds[k] == "discretisation"

    def _hk(self, k):
        """Refinement scalar of axis k in the form h^p (exponential axes are
        mapped as in the fit, so one formula serves both models)."""
        if self.models is not None and self.models[k] == "exponential":
            x = np.asarray(self.hs[k], float)
            span = max(x.max() - x.min(), 1e-300)
            return np.exp(-(x - x.min()) / span)
        return np.asarray(self.hs[k], float)

    def _local_check(self, k):
        """Three-point extrapolation along axis k through the anchor.

        Uses the anchor's coarser and finer neighbours, or, at the coarse end of
        the ladder, the two finer values.  With d1, d2 the two successive changes,
        the observed order p solves (h0^p - h1^p)/(h1^p - h2^p) = d1/d2 and the
        error left at the anchor is c h_a^p with c fixed by d2.  If the changes do
        not decrease (no positive order), the remaining error is unbounded; if
        they alternate in sign with decreasing size, it is at most 2|d2| beyond
        the last point (alternating-series bound).
        """
        a = self.anchor[k]
        ok = lambda j: 0 <= j < self.m[k] and j in self.sampled[k] and np.isfinite(self.e[k][j])
        if ok(a - 1) and ok(a + 1):
            j0, j1, j2 = a - 1, a, a + 1
        elif a == 0 and ok(1) and ok(2):
            j0, j1, j2 = 0, 1, 2
        elif self.fine_failed[k] and ok(a - 2) and ok(a - 1) and np.isfinite(self.e[k][a]):
            # the finer run failed, so nothing is known on the refining side: the
            # order is taken from the three values ending at the anchor
            j0, j1, j2 = a - 2, a - 1, a
        else:
            return None
        e = self.e[k]
        d1, d2 = e[j1] - e[j0], e[j2] - e[j1]
        res = 3.0 * self.resolution
        if abs(d1) <= res and abs(d2) <= res:
            return dict(p=np.inf, E=0.0, unbounded=False, pts=(j0, j1, j2))
        h = self._hk(k)
        tail_from = e[j2] - e[a]              # change from the anchor to the last point
        if d1 * d2 < 0 or abs(d1) <= res:
            # alternating (or a first change below resolution): bounded only if
            # the changes shrink
            if abs(d2) < abs(d1):
                E = abs(tail_from) + 2.0 * abs(d2)
                return dict(p=np.nan, E=E, unbounded=False, pts=(j0, j1, j2))
            return dict(p=np.nan, E=np.inf, unbounded=True, pts=(j0, j1, j2))
        if abs(d2) <= res:
            # the second change is below resolution: nothing is left beyond the
            # last point, so the anchor carries what lies between it and that point
            return dict(p=np.inf, E=abs(tail_from), unbounded=False, pts=(j0, j1, j2))
        ratio = d1 / d2
        f = lambda p: (h[j0]**p - h[j1]**p) / (h[j1]**p - h[j2]**p) - ratio
        lo, hi = 1e-6, 30.0
        if f(lo) >= 0:                         # changes do not decrease
            return dict(p=0.0, E=np.inf, unbounded=True, pts=(j0, j1, j2))
        if f(hi) < 0:
            p = hi
        else:
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
            p = 0.5 * (lo + hi)
        c = d2 / (h[j2]**p - h[j1]**p)
        E = abs(c * h[a]**p)
        return dict(p=p, E=E, unbounded=False, pts=(j0, j1, j2))

    # ------------------------------------------- anchor 2a: per-knob residual
    def per_knob_error(self, p_grid=None):
        """Error remaining at the anchor, knob by knob, from the coarse side.

        The sweep already samples each knob's whole coarse side, and those runs
        are cheap.  Fitting

            e_k(h) = c_k ( h^{p_k} - h_anchor^{p_k} )

        to them recovers that knob's own convergence order -- one rate per knob,
        rather than the single blended ratio a jointly-refined path yields -- and
        r_k = c_k h_anchor^{p_k} is the error the anchor still carries along that
        axis.  Summing r_k gives the absolute error, and the r_k *are* the error
        budget, so the estimate and the decomposition come from the same fit.
        A masked knob fits c_k ~ 0 and contributes nothing, which is correct.
        """
        if p_grid is None:
            p_grid = np.linspace(0.25, 6.0, 121)
        self.rk, self.pk, self.fit_q = [], [], []
        self.fit_rms = {}
        for k in range(len(self.m)):
            # only runs that succeeded enter the fit
            js = [j for j in self.sampled[k] if np.isfinite(self.e[k][j])]
            if self.hs is None or len(js) < 3:
                r = np.nan
                a = self.anchor[k]
                if self.unfit in ("add", "order1") and len(self.sampled[k]) > 1:
                    if a + 1 in js:
                        r = -self.e[k][a + 1]       # remaining error >= one finer step
                        if self.unfit == "order1" and self.hs is not None:
                            # if the axis converges at least to first order, the
                            # error left after the step is at most step/(ratio-1)
                            ha, hb = self.hs[k][a], self.hs[k][a + 1]
                            if self.models is not None and self.models[k] == "exponential":
                                pass
                            elif ha > hb > 0:
                                r = r * (1.0 + hb / (ha - hb))
                    elif a - 1 in js:
                        r = self.e[k][a - 1]        # bounded by the last coarser step
                self.rk.append(r); self.pk.append(np.nan); self.fit_q.append(np.inf)
                continue
            model = self.models[k] if self.models is not None else "power"
            if model == "exponential":
                # x = setting; h = exp(-(x - x_min)/span) turns exp(-b x) into
                # h^p with p = b * span, so the same grid search applies.
                x = np.asarray(self.hs[k], float)
                span = max(x.max() - x.min(), 1e-300)
                hk = np.exp(-(x - x.min()) / span)
                grid = np.linspace(0.25, 40.0, 160)
            elif model == "power":
                hk, grid = self.hs[k], p_grid
            else:
                raise ValueError(f"unknown error model {model!r}")
            h = hk[js]
            e = self.e[k][js]
            ha = hk[self.anchor[k]]
            best = (np.inf, np.nan, np.nan)
            for p in grid:
                basis = h**p - ha**p
                d = basis @ basis
                if d <= 0:
                    continue
                c = (basis @ e) / d
                res = np.sum((e - c * basis) ** 2)
                if res < best[0]:
                    best = (res, p, c)
            res, p, c = best
            self.rk.append(c * ha**p)
            self.pk.append(p)
            self.fit_q.append(np.sqrt(res / max(np.sum(e**2), 1e-300)))
            self.fit_rms[k] = float(np.sqrt(res / len(js)))
        # An axis whose finer run failed has no measurement on the refining
        # side.  Treating it as converged is what certified two settings that
        # missed on LAMMPS at 262,144 ions.  Such an axis is charged at least
        # its last measured one-step change; if it has none, the setting is
        # declined (E_knob = inf).
        self.unmeasured = []
        for k in range(len(self.m)):
            if not self.fine_failed[k]:
                continue
            a = self.anchor[k]
            step = abs(self.e[k][a - 1]) if a >= 1 and np.isfinite(self.e[k][a - 1]) \
                else np.nan
            if not np.isfinite(step):
                self.unmeasured.append(k)
                self.rk[k] = np.inf
                continue
            r = self.rk[k]
            if not np.isfinite(r):
                # too few successful runs on this axis for a fit
                self.rk[k] = step
            elif self.unmeasured_rule == "max":
                if abs(r) < step:
                    self.rk[k] = np.sign(r or 1.0) * step
        if self.local_order:
            for k in range(len(self.m)):
                if not self._disc(k) or self.hs is None:
                    continue
                L = self._local_check(k)
                if L is None:
                    continue
                self.local[k] = L
                r = self.rk[k]
                if L["unbounded"]:
                    self.rk[k] = np.inf
                elif not np.isfinite(r) or L["E"] > abs(r):
                    sgn = np.sign(r) if np.isfinite(r) and r != 0 else \
                        -np.sign(self.e[k][L["pts"][2]] - self.e[k][self.anchor[k]] or 1.0)
                    self.rk[k] = sgn * L["E"]
        self.E_knob = float(np.nansum(self.rk))
        return self.E_knob

    # ------------------------------------------------------- joint check
    def joint_correct(self, lookup, eps=None, cancel_rule="credit"):
        """Measure the interaction of the parameters at this anchor.

        Refine every discretisation and nested-solver parameter whose one-step
        change is resolved by one position together (one extra run), and compare
        the observed change with the sum of the one-at-a-time changes.  The
        difference I is the interaction term of Appendix A, measured where the
        estimate is used.  The estimate becomes

            Ehat = Fs * max( |E_fit| + 2|I|,  max_k |delta_k|,  |delta_joint| ),

        and if |I| exceeds half the largest one-at-a-time change, the errors of
        different parameters are not allowed to cancel (|E_fit| is replaced by
        sum_k |r_k|).  If the joint run fails, cancellation is not credited.
        """
        a = self.anchor
        R = []
        for k in range(len(self.m)):
            if self.kinds is not None and self.kinds[k] == "cost-only":
                continue
            j = a[k] + 1
            if j >= self.m[k] or j not in self.sampled[k] or not np.isfinite(self.e[k][j]):
                continue
            if abs(self.e[k][j]) <= self.joint_res * self.resolution or self.e[k][j] == 0.0:
                continue
            R.append(k)
        if len(R) < 2 or not np.isfinite(self.E0):
            self.joint = dict(axes=R, skipped=True)
            return self.E0
        idx = tuple(a[k] + 1 if k in R else a[k] for k in range(len(self.m)))
        self._tag = "joint"
        Q, C = self._get(idx)
        rk = [r for r in self.rk if np.isfinite(r)]
        E_fit, E_abs = float(np.sum(rk)), float(np.sum(np.abs(rk)))
        steps = [abs(self.e[k][a[k] + 1]) for k in R]
        floor = max(steps)
        if not (np.isfinite(Q) and np.isfinite(C)):
            I, dJ, strong = np.nan, 0.0, True
            core = E_abs
        else:
            dJ = Q - self.Q0
            I = dJ - sum(self.e[k][a[k] + 1] for k in R)
            strong = abs(I) > 0.5 * floor
            core = (E_abs if strong else abs(E_fit)) + 2.0 * abs(I)
        E_path = np.nan
        # cancellation: the signed sum is smaller than the largest single term
        big = max((abs(r) for r in rk), default=0.0)
        cancels = big > 0 and abs(E_fit) < 0.8 * big
        if cancels and cancel_rule == "none":
            core = E_abs + 2.0 * (abs(I) if np.isfinite(I) else 0.0)
        elif cancels and cancel_rule == "path" and np.isfinite(dJ):
            # verify the cancellation along the jointly refined path: one more run
            # two positions finer in every parameter of the joint step
            idx2 = tuple(min(a[k] + 2, self.m[k] - 1) if k in R else a[k]
                         for k in range(len(self.m)))
            if idx2 != idx:
                self._tag = "joint2"
                Q2, C2 = self._get(idx2)
                if np.isfinite(Q2) and np.isfinite(C2):
                    d1, d2 = dJ, Q2 - Q
                    if d1 * d2 > 0 and abs(d2) < abs(d1):
                        q = d2 / d1
                        E_path = abs(d1 + d2 / (1.0 - q))
                    elif d1 * d2 < 0 and abs(d2) < abs(d1):
                        E_path = abs(d1) + abs(d2)
                    else:
                        E_path = np.inf
                else:
                    E_path = np.inf
                core = max(core, E_path)
        E = max(core, floor, abs(dJ))
        sgn = np.sign(self.E0) if self.E0 != 0 else 1.0
        E0_new = sgn * self.Fs * (E + getattr(self, "add_margin", 0.0))
        self.joint = dict(axes=R, idx=idx, I=I, dJ=dJ, strong=strong, E_path=E_path,
                          E0_before=self.E0, E0_after=E0_new, skipped=False)
        if abs(E0_new) > abs(self.E0):
            self.E0 = E0_new
            self.E0_reliable = False
        return self.E0

    # ------------------------------------------------------------- anchor 2b
    def refine_path(self, n_levels: int = 3):
        """Absolute error at the anchor, from a jointly-refined path.

        n_levels=0 skips the path entirely and relies on the per-knob fit plus
        the consistency floor, both of which come free with the sweep.
        """
        if n_levels < 2:
            E, rel = float(getattr(self, "E_knob", 0.0)), False
            self.path_ratio = np.nan
            if self.cancel_rule == "none" and np.isfinite(E):
                # Version 1.2: no credit for cancellation between parameters.  If
                # the signed sum is smaller than the largest single term, the
                # terms are summed in magnitude (the sign of the largest is kept).
                rk = [r for r in self.rk if np.isfinite(r)]
                big = max(rk, key=abs) if rk else 0.0
                if rk and abs(E) < 0.8 * abs(big):
                    E = float(np.sign(big) * np.sum(np.abs(rk)))
                    self.cancelled = True
        else:
            Qs = []
            for l in range(n_levels):
                idx = tuple(min(self.anchor[k] + l, self.m[k] - 1)
                            for k in range(len(self.m)))
                self._tag = f"path_L{l}"
                Qs.append(self._get(idx)[0])
            d1 = Qs[1] - Qs[0]
            d2 = (Qs[2] - Qs[1]) if n_levels > 2 else None
            if abs(d1) < 1e-300:
                E, rel = 0.0, True
            elif d2 is None:
                E, rel = 2.0 * d1, False        # no ratio available: crude margin
            else:
                r = d2 / d1
                if 0 < r < 1:                   # clean geometric convergence
                    E, rel = d1 / (r - 1), True
                else:                           # non-monotone or stalled
                    E, rel = d1 / (0.5 - 1), False
            self.path_ratio = (d2 / d1) if (d2 is not None and abs(d1) > 1e-300) else np.nan

        # Consistency floor.  A single Richardson ratio is meaningless when the
        # knobs converge at very different rates -- refining N by 4x and a
        # tolerance by 100x in the same step does not produce one geometric
        # sequence.  But the sweep already bounds the error from below: if
        # *refining* knob k alone moved the observable by delta, the anchor's
        # error is at least about delta.  Only the refining direction counts:
        # that coarsening a knob changes the answer says the coarse point is
        # bad, not that the anchor is.  (Using the whole curve here makes the
        # floor enormous and the protocol recommends the finest setting always.)
        floor = 0.0
        for k in range(len(self.m)):
            fine = [j for j in self.sampled[k]
                    if j > self.anchor[k] and np.isfinite(self.e[k][j])]
            if fine:
                floor = max(floor, abs(self.e[k][max(fine)]))
        if abs(E) < floor:
            E, rel = np.sign(E if E != 0 else 1.0) * floor, False

        # Two independent estimators: this jointly-refined path, and the
        # per-knob fit to the coarse side.  Take the larger magnitude, and let
        # the *disagreement between them* set the margin.  A single estimator
        # has no way to know it is in a pre-asymptotic regime; two estimators
        # that rely on different structure do -- when they agree the anchor is
        # asymptotic and a small margin suffices, and when they diverge the
        # ratio itself measures how little we should trust either.
        E_path = E
        E_knob = getattr(self, "E_knob", np.nan)
        if np.isfinite(E_knob):
            if abs(E_knob) > abs(E):
                E = E_knob
            a, bmag = abs(E_path), abs(E_knob)
            lo, hi = min(a, bmag), max(a, bmag)
            disagree = hi / lo if lo > 1e-300 else np.inf
        else:
            disagree = np.inf
        self.disagreement = disagree
        if rel and disagree < 2.0:
            self.Fs = 1.25
        elif disagree < 5.0:
            self.Fs = 2.0
        else:
            self.Fs = 3.0
        if self.fs_cap is not None:
            self.Fs = min(self.Fs, self.fs_cap)
        # Residual-gated margin.  A fixed margin is calibrated against ladders
        # of a given coarseness and silently loses its protection on finer
        # ones: as the sampled range of h shrinks, the offset fit is less
        # determined and the estimate, when it is optimistic, is optimistic by
        # far more (measured: the 10th percentile of Ehat/Etrue falls from
        # 0.93 to 0.39 as the spacing closes from 2x to 1.2x).  The fit's own
        # normalised residual tracks that degradation -- it rises 47-fold over
        # the same range -- and the method already computes it, so the margin
        # is scaled by it instead of being a constant chosen once.
        add = 0.0
        if self.resid_margin and self.resid_mode == "additive":
            # Version 1.2: the misfit of each axis's error model, in the units of
            # the observable, is an uncertainty of its extrapolation; it is added
            # (root-sum-square over the axes, cost-only axes excluded) instead of
            # scaling the whole estimate, so an axis whose changes are far below
            # the target cannot inflate it.
            rms = [v for k, v in getattr(self, "fit_rms", {}).items()
                   if np.isfinite(v) and (self.kinds is None or self.kinds[k] != "cost-only")]
            add = self.resid_margin * float(np.sqrt(np.sum(np.square(rms)))) if rms else 0.0
            self.fit_residual = add
        elif self.resid_margin:
            good = [q for q in getattr(self, "fit_q", []) if np.isfinite(q)]
            r = float(np.median(good)) if good else 0.0
            self.fit_residual = r
            self.Fs *= 1.0 + self.resid_margin * r
        else:
            self.fit_residual = np.nan
        self.add_margin = add
        self.E0 = (abs(E) + add) * self.Fs * (np.sign(E) if E != 0 else 1.0)
        self.E0_reliable = bool(rel and disagree < 2.0)
        return self.E0

    # -------------------------------------------------------------- models
    def model_error(self, idx):
        return float(self.E0 + sum(self.e[k][i] for k, i in enumerate(idx)))

    def model_abs_error(self, idx):
        return float(abs(self.E0) + sum(abs(self.e[k][i]) for k, i in enumerate(idx)))

    def model_logcost(self, idx):
        return float(np.log(self.C0) + sum(self.g[k][i] for k, i in enumerate(idx)))

    # ---------------------------------------------------------- diagnostics
    def thresholds(self, kappa: float = 1e-2):
        """Local masking threshold: coarsest setting past which the knob is inert."""
        out = []
        for k in range(len(self.m)):
            dev = np.abs(self.e[k] - self.e[k][self.m[k] - 1])
            emax = dev.max()
            if emax <= 0:
                thr = 0
            else:
                above = np.where(dev > kappa * emax)[0]
                thr = min(int(above[-1]) + 1, self.m[k] - 1) if len(above) else 0
            out.append(dict(k=k, threshold_idx=thr, setting=self.ladders[k][thr],
                            masked=thr <= self.m[k] - 2,
                            past=self.anchor[k] - thr))
        return out

    def shares(self, idx=None):
        idx = self.anchor if idx is None else idx
        c = np.array([abs(self.e[k][idx[k]] - self.e[k][self.m[k] - 1])
                      for k in range(len(self.m))])
        return c / c.sum() if c.sum() > 0 else c

    # ------------------------------------------------------------- proposal
    def propose(self, eps: float, mode: str = "optimistic"):
        err = (lambda i: abs(self.model_error(i))) if mode == "optimistic" \
            else self.model_abs_error
        # Trust region: only configurations we actually sampled.  Anywhere else
        # we would be extrapolating the very couplings that made the global
        # design fail, and a one-shot jump to the finest setting -- the obvious
        # fallback -- costs two orders of magnitude more than the target needs.
        best, bc = None, np.inf
        for idx in itertools.product(*self.sampled):
            if err(idx) <= eps:
                c = self.model_logcost(idx)
                if c < bc:
                    best, bc = idx, c
        if best is not None:
            self.infeasible = False
            return best, np.exp(bc)

        # Target not reachable inside the window: take the single most
        # cost-effective refinement step and re-anchor there.
        self.infeasible = True
        e_now = err(self.anchor)
        cand, gain_best = None, -np.inf
        for k in range(len(self.m)):
            j = self.anchor[k] + 1
            if j >= self.m[k] or j not in self.sampled[k]:
                continue
            trial = tuple(j if i == k else self.anchor[i] for i in range(len(self.m)))
            d_cost = self.g[k][j] - self.g[k][self.anchor[k]]
            gain = (e_now - err(trial)) / max(d_cost, 1e-12)
            if gain > gain_best:
                cand, gain_best = trial, gain
        if cand is None:
            return self.anchor, self.C0
        return cand, np.exp(self.model_logcost(cand))


def screen(cached, ladders, anchor, kinds, tau=0.05, drop=2, at="tight"):
    """Decide which knobs enter the search, cheaply.

    Influence is measured at a configuration whose DISCRETISATION knobs have been
    coarsened by `drop` rungs, which makes every probe roughly an order of
    magnitude cheaper.  Knobs holding less than `tau` of the influence there are
    frozen.

    Nested-solver knobs are exempt and always enter the search.  Measured over
    24 parameter combinations on the 8-knob problem, every screening miss fell
    on that class (newton_tol, 9/24) and none on discretisation or cost-only
    knobs: coarsening the grid changes the conditioning of the problem those
    tolerances are solving, so their irrelevance at low cost says nothing about
    production resolution.

    A declaration can also be falsified here: a knob declared cost-only that
    visibly moves the observable is reported and forced back into the search.
    Inference in the other direction is not possible -- a masked nested-solver
    knob is observationally identical to a cost-only one -- so classes are
    declared, then falsified, never inferred.
    """
    nk = len(ladders)
    # Screen where nothing masks what we are trying to see: coarsen the
    # DISCRETISATION knobs (expensive, and what we want to make cheap) while
    # setting the NESTED-SOLVER knobs to their finest (cheap, and otherwise the
    # dominant error source at a coarse configuration).  Screening at the
    # anchor's own tolerances instead lets a loose newton_tol swamp everything
    # and the grid itself measures as irrelevant -- which is how screening
    # first lost 2.1% of its safety.
    if at == "user":
        # Ablation only: screen at the user's own configuration.
        sc = tuple(anchor)
    else:
        sc = tuple(max(0, anchor[k] - drop) if kinds[k] == "discretisation"
                   else (len(ladders[k]) - 1 if kinds[k] == "nested-solver" else anchor[k])
                   for k in range(nk))
    base = cached(sc)[0]
    infl = {}
    for k in range(nk):
        if kinds[k] == "nested-solver":
            continue
        j = min(sc[k] + 1, len(ladders[k]) - 1)
        if j == sc[k]:
            infl[k] = 0.0
            continue
        cc = tuple(j if i == k else sc[i] for i in range(nk))
        infl[k] = abs(cached(cc)[0] - base)
    tot = sum(infl.values()) or 1.0
    frozen, violations = [], []
    for k, v in infl.items():
        share = v / tot
        if share < tau:
            frozen.append(k)
        elif kinds[k] == "cost-only":
            violations.append(k)            # declared inert, demonstrably is not
    active = [k for k in range(nk) if k not in frozen]
    return active, frozen, violations, sc


# The default configuration of version 1.2 (Section 3 of the article).  The
# keyword defaults of run_protocol itself reproduce version 1.1, so that the
# ablations and the comparisons with version 1.1 can be re-run exactly.
V12 = dict(resid_margin=5.0, joint_check=True, local_order_check=True, res_frac=0.1,
           cancel_rule="none")
V11 = dict(resid_margin=0.0, joint_check=False, local_order_check=False, res_frac=0.0,
           cancel_rule="credit")
# Version 1.2.1: the joint check includes every parameter whose one-step change exceeds
# the resolution, not three times it.  Found by the held-out problem H2, where the check
# was skipped at a setting that missed (PREREGISTRATION.md, amendment).
V121 = dict(V12, joint_res=1.0)


def run_protocol(lookup, ladders, anchor, eps, mode="optimistic", max_iter=40,
                 fine_reach=1, hs=None, path_levels=0, fs_cap=1.25,
                 kinds=None, use_screen=False, tau=0.05, drop=2, audit=True,
                 resid_margin=0.0, walk_mode="monotone", phase2="both",
                 screen_at="tight", audit_factor=2.0, models=None,
                 noise_floor=None, floor_factor=10.0, unfit="order1",
                 joint_check=False, local_order_check=False, cancel_rule="credit",
                 resid_mode="relative", res_frac=0.0, joint_res=3.0):
    """Screen, walk, audit.

    Every configuration is evaluated at most once for the whole protocol.

    The last four keyword arguments exist only to reproduce the design
    ablations of the article (Table 4); the defaults are the method:
      walk_mode="propose"   jump to the model's predicted optimum instead of the
                            monotone refine/improve walk;
      phase2="coarsen"      Phase 2 may only coarsen;
      screen_at="user"      screen at the user's own setting;
      audit_factor=1.0      the audit uses the one-step difference itself.
    Successive anchors overlap heavily, so without the cache each step re-pays
    for most of the previous sweep.
    """
    # Round-off floor.  If the solver cannot resolve the observable to better
    # than `noise_floor` (measured beforehand, e.g. by `measure_floor`), no
    # parameter setting can certify a target near it: decline at once.
    if noise_floor is not None and eps < floor_factor * noise_floor:
        return dict(proposal=tuple(anchor), runs=0, protocol_cost=0.0, iters=0,
                    hit_cap=False, cycled=False, reliable=False, certified=False,
                    infeasible=True, budget=None, breakdown={}, frozen=[],
                    violations=[], repairs=0, trace=[], screen_cfg=None,
                    declined="below_floor")
    cache: dict = {}
    charged = [0.0]
    distinct = [0]
    breakdown: dict = {}
    pending = ["?"]

    def _note_tag(idx, tag):
        pending[0] = tag

    def cached(idx):
        idx = tuple(idx)
        if idx not in cache:
            val = lookup(idx)                 # exactly one evaluation per config
            cache[idx] = val
            charged[0] += val[1]
            distinct[0] += 1
            breakdown[pending[0]] = breakdown.get(pending[0], 0.0) + val[1]
        return cache[idx]

    def analyse(idx, active):
        lb = LocalBudget(cached, ladders, idx, hs=hs, fs_cap=fs_cap,
                         resid_margin=resid_margin, models=models, unfit=unfit,
                         kinds=kinds, local_order=local_order_check,
                         resolution=max(noise_floor or 0.0, res_frac * eps),
                         resid_mode=resid_mode, cancel_rule=cancel_rule, joint_res=joint_res)
        lb.on_charge = _note_tag
        lb.sweep(fine_reach=fine_reach, active=active)
        lb.per_knob_error()
        lb.refine_path(n_levels=path_levels)
        if not np.isfinite(lb.Q0):
            lb.E0 = np.inf                     # the anchor itself failed
        elif joint_check and np.isfinite(lb.E0) and abs(lb.E0) <= eps:
            # only where the estimate would be accepted: one extra run
            lb.joint_correct(cached, eps, cancel_rule=cancel_rule)
        return lb

    def finite_err(lb, trial):
        v = lb.model_error(trial)
        return v if np.isfinite(v) else None

    nk = len(ladders)
    frozen, violations, screen_cfg = [], [], None
    active = list(range(nk))
    if use_screen and kinds is not None:
        pending[0] = "screen"
        active, frozen, violations, screen_cfg = screen(cached, ladders, anchor,
                                                        kinds, tau=tau, drop=drop,
                                                        at=screen_at)

    trace = []

    def walk_propose(start, active, budget):
        """Ablation: repeatedly jump to the model's predicted optimum."""
        cur, steps = tuple(start), 0
        trace.append(cur)
        lb = analyse(cur, active)
        seen = {cur}
        while steps < budget:
            nxt, _ = lb.propose(eps, mode="optimistic")
            nxt = tuple(nxt)
            if nxt == cur:
                break
            cur, lb, steps = nxt, analyse(nxt, active), steps + 1
            trace.append(cur)
            if cur in seen:
                walk_propose.cycled = True
            seen.add(cur)
        return cur, lb, abs(lb.E0) <= eps, steps

    def walk(start, active, budget):
        if walk_mode == "propose":
            return walk_propose(start, active, budget)
        cur, steps = tuple(start), 0
        trace.append(tuple(start))
        lb = analyse(cur, active)
        if not np.isfinite(lb.Q0):
            walk.start_failed = True
            return cur, lb, False, steps
        # Phase 1 -- refine while the anchor's own estimate misses the target.
        # Every step increases an index, so this terminates.
        while abs(lb.E0) > eps and steps < budget:
            best, gain_best = None, -np.inf
            if not np.isfinite(lb.E0):
                # an axis whose error is unbounded at this anchor (its changes do
                # not decrease): refine it, the cheapest first, to learn more
                cand = [k for k in active
                        if not np.isfinite(lb.rk[k]) and cur[k] + 1 < lb.m[k]
                        and (cur[k] + 1) in lb.sampled[k]
                        and np.isfinite(lb.e[k][cur[k] + 1])]
                if cand:
                    k = min(cand, key=lambda k: lb.g[k][cur[k] + 1] - lb.g[k][cur[k]])
                    best = tuple(cur[k] + 1 if i == k else cur[i] for i in range(nk))
            else:
                for k in active:
                    j = cur[k] + 1
                    if j >= lb.m[k] or j not in lb.sampled[k]:
                        continue
                    trial = tuple(j if i == k else cur[i] for i in range(nk))
                    me = finite_err(lb, trial)
                    if me is None or not np.isfinite(lb.g[k][j]):
                        continue
                    d_err = abs(lb.E0) - abs(me)
                    d_cost = lb.g[k][j] - lb.g[k][cur[k]]
                    gain = d_err / max(d_cost, 1e-12)
                    if gain > gain_best:
                        best, gain_best = trial, gain
            if best is None:
                break
            cur, lb, steps = best, analyse(best, active), steps + 1
            trace.append(cur)
        certified = abs(lb.E0) <= eps
        # Phase 2 -- any single-knob move, either direction, that strictly lowers
        # modelled cost and keeps the estimate in target.  Cost is a strict
        # potential, so this terminates.  Coarsen-only cannot reach the interior
        # cost minimum of an accuracy-neutral knob.
        blocked = set()
        while steps < budget:
            c_now = lb.model_logcost(cur)
            cands = []
            for k in active:
                if k in blocked:
                    continue
                for d in ((-1,) if phase2 == "coarsen" else (-1, +1)):
                    j = cur[k] + d
                    if j < 0 or j >= lb.m[k] or j not in lb.sampled[k]:
                        continue
                    trial = tuple(j if i == k else cur[i] for i in range(nk))
                    c_try = lb.model_logcost(trial)
                    me = finite_err(lb, trial)
                    if me is None or not np.isfinite(c_try) or abs(me) > eps:
                        continue
                    cheaper = c_try < c_now - 1e-12
                    # A move that coarsens at no measurable cost saving is still
                    # worth taking: it is the difference between telling a user
                    # a tolerance is six orders past useful and actually
                    # loosening it.  Termination is preserved because the
                    # potential is lexicographic -- cost first, then the sum of
                    # ladder indices, which such a move strictly decreases.
                    tidy = (not cheaper) and j < cur[k] and c_try <= c_now + 1e-12
                    if cheaper or tidy:
                        cands.append((c_try - (1e-9 if tidy else 0.0), k, trial))
            if not cands:
                break
            cands.sort()
            _, k, trial = cands[0]
            lb2 = analyse(trial, active)
            steps += 1
            if abs(lb2.E0) <= eps:
                cur, lb = trial, lb2
                trace.append(cur)
            else:
                blocked.add(k)
        return cur, lb, certified, steps

    pending[0] = "walk"
    walk.start_failed = False
    cur, lb, certified, steps = walk(anchor, active, max_iter)
    if walk.start_failed:
        return dict(proposal=tuple(anchor), runs=distinct[0], protocol_cost=charged[0],
                    iters=0, hit_cap=False, cycled=False, reliable=False, certified=False,
                    infeasible=True, budget=lb, breakdown=breakdown, frozen=frozen,
                    violations=violations, repairs=0, trace=trace, screen_cfg=screen_cfg,
                    declined="start_failed")

    # Audit -- probe each frozen knob once at the answer and ADD its measured
    # deviation to the error budget.  Testing each frozen knob individually
    # against the full target is wrong: several knobs each contributing a
    # fraction of eps pass one at a time and bust the budget together, which is
    # exactly how screening lost 2.1% of its safety.  Frozen knobs are part of
    # the budget, so they are summed into it.
    repairs = 0
    frozen_contrib = 0.0
    if audit and frozen:
        for _ in range(3):
            pending[0] = "audit"
            q0 = cached(cur)[0]
            dev = {}
            for k in list(frozen):
                j = min(cur[k] + 1, len(ladders[k]) - 1)
                if j == cur[k]:
                    dev[k] = 0.0
                    continue
                cc = tuple(j if i == k else cur[i] for i in range(nk))
                # A one-step difference is the FIRST TERM of what the knob still
                # owes, not the whole of it.  For a converging knob the tail sums
                # to dev/(1-r); with r bounded by 1/2 (first order under 2x
                # refinement) a factor of 2 is the conservative reading.  Using
                # the raw difference certified a run 25% over budget.
                dev[k] = audit_factor * abs(cached(cc)[0] - q0)
            frozen_contrib = sum(dev.values())
            if abs(lb.E0) + frozen_contrib <= eps:
                break
            # unfreeze the largest contributors until the budget could close
            order = sorted(dev, key=lambda k: -dev[k])
            moved = False
            for k in order:
                if dev[k] <= 0:
                    continue
                frozen.remove(k); active.append(k); moved = True
                if abs(lb.E0) + sum(dev[j] for j in frozen) <= eps:
                    break
            if not moved:
                break
            repairs += 1
            pending[0] = "repair"
            cur, lb, certified, s2 = walk(cur, sorted(active), max_iter)
            steps += s2
        certified = certified and (abs(lb.E0) + frozen_contrib <= eps)

    # Cost-only pass -- knobs declared accuracy-neutral are swept once at the
    # answer and set to their cheapest verified value.  Freezing them at the
    # user's setting is safe but leaves the interior cost minimum unclaimed,
    # which is most of the recommendation penalty screening otherwise incurs.
    if kinds is not None:
        pending[0] = "cost_pass"
        for k in range(nk):
            if kinds[k] != "cost-only":
                continue
            q0, c0 = cached(cur)
            best_j, best_c = cur[k], c0
            for j in range(len(ladders[k])):
                if j == cur[k]:
                    continue
                cc = tuple(j if i == k else cur[i] for i in range(nk))
                q, c = cached(cc)
                if c < best_c and abs(q - q0[0] if isinstance(q0, tuple) else q - q0) \
                        <= 0.05 * eps:
                    best_j, best_c = j, c
            if best_j != cur[k]:
                cur = tuple(best_j if i == k else cur[i] for i in range(nk))

    return dict(proposal=cur, runs=distinct[0], protocol_cost=charged[0],
                iters=steps, hit_cap=steps >= max_iter,
                cycled=bool(getattr(walk_propose, "cycled", False)), reliable=lb.E0_reliable, certified=certified,
                infeasible=not certified, budget=lb, breakdown=breakdown,
                frozen=frozen, violations=violations, repairs=repairs, trace=trace,
                screen_cfg=screen_cfg,
                declined=None if certified else (
                    "unmeasured_axis" if getattr(lb, "unmeasured", []) else
                    ("unbounded_axis" if not np.isfinite(lb.E0) else "not_reached")))


def measure_floor(values):
    """Round-off floor from repeated evaluations that should agree.

    `values` are the observable at two or more settings refined beyond every
    candidate value (or differing only in thread count or another parameter
    with no discretisation effect).  Their spread is a lower bound on the
    error any setting can be certified to.
    """
    v = np.asarray([x for x in values if np.isfinite(x)], float)
    if len(v) < 2:
        raise ValueError("need at least two finite values to measure a floor")
    return float(v.max() - v.min())
