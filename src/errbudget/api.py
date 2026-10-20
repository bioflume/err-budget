"""User-facing API.

The contract is deliberately small.  You supply a function that runs your solver,
a description of each numerical parameter, the quantity you care about and how
accurate you need it.  You get back a setting, an error estimate, and a statement
of whether the target was met -- plus, importantly, which parameters are being
refined for nothing.

    from errbudget import Parameter, tune

    def run(N, dt, rtol):
        ...                      # your solver
        return {"drag": cd}      # or just: return cd

    result = tune(
        run,
        [Parameter("N",    [32, 64, 128, 256],       kind="discretisation", h=lambda v: 1/v),
         Parameter("dt",   [1e-2, 5e-3, 2.5e-3],     kind="discretisation"),
         Parameter("rtol", [1e-4, 1e-6, 1e-8, 1e-10], kind="nested-solver")],
        observable="drag", target=1e-3, start={"N": 64, "dt": 5e-3, "rtol": 1e-10})

    print(result)
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any, Callable, Sequence

import numpy as np

from ._core import run_protocol, V12, V121

KINDS = ("discretisation", "nested-solver", "cost-only")


class RunFailed(Exception):
    """Raise this from the solver function to report a run that failed (crashed,
    diverged, or whose nested solver did not converge).  The program records it
    as unusable and never proposes or fits to it."""


class Parameter:
    """One numerical parameter of the solver.

    Parameters
    ----------
    name : the keyword the solver function expects.
    values : candidate settings, ordered COARSE (cheap, inaccurate) to FINE.
    kind : one of

        "discretisation"  refining improves the answer and costs more
                          (mesh spacing, time step, basis size, cutoff)
        "nested-solver"   a tolerance or iteration limit of a solver nested
                          inside the outer problem; what it does depends on the
                          conditioning that the discretisation sets, so it is
                          never screened out cheaply
        "cost-only"       changes cost but not the answer, often with an
                          interior cost optimum (preconditioner fill or drop
                          tolerance, how often a Jacobian is refreshed)

        The class cannot be inferred reliably from data -- a masked
        nested-solver parameter is indistinguishable from an inert one -- so it
        is declared here and falsified by the tool, which reports any
        declaration its own measurements contradict.
    h : optional map from a setting to a refinement scalar that decreases as the
        setting gets finer.  Defaults to the value itself for tolerances, and to
        1/value for counts (inferred when the values increase).
    model : "power" (default) fits e(h) = c (h^p - h_a^p) in the refinement
        scalar h.  "exponential" fits e(x) = c (exp(-b x) - exp(-b x_a)) in the
        setting x itself, which must increase as it refines; use it for
        spectral or plane-wave cutoffs, whose error decays exponentially.
    """

    def __init__(self, name: str, values: Sequence[Any], kind: str = "discretisation",
                 h: Callable[[Any], float] | None = None, model: str = "power"):
        if model not in ("power", "exponential"):
            raise ValueError(f"model must be 'power' or 'exponential', got {model!r}")
        self.model = model
        if model == "exponential":
            x = [float(v) for v in values]
            if not all(x[i + 1] > x[i] for i in range(len(x) - 1)):
                raise ValueError(f"parameter {name!r}: an exponential-model parameter must "
                                 f"increase from coarse to fine, got {values!r}")
            if len(values) < 2:
                raise ValueError(f"parameter {name!r}: need at least two candidate values")
            if kind not in KINDS:
                raise ValueError(f"kind must be one of {KINDS}, got {kind!r}")
            self.name, self.values, self.kind, self.h = name, list(values), kind, x
            return
        if kind not in KINDS:
            raise ValueError(f"kind must be one of {KINDS}, got {kind!r}")
        if len(values) < 2:
            raise ValueError(f"parameter {name!r}: need at least two candidate values")
        self.name, self.values, self.kind = name, list(values), kind
        if h is not None:
            self.h = [float(h(v)) for v in self.values]
        else:
            v = np.asarray(self.values, dtype=float)
            increasing = v[-1] > v[0]
            if increasing and np.any(np.abs(v) < 1.0):
                # An increasing ladder of sub-unit values is ambiguous: it looks
                # like a count ladder but is far more often a tolerance ladder
                # written fine-to-coarse.  Guessing here silently accepts a
                # reversed ladder, so refuse and ask.
                raise ValueError(
                    f"parameter {name!r}: cannot infer the refinement direction from "
                    f"{self.values!r}. Values increase but are below 1, which "
                    f"looks like a tolerance ladder written fine-to-coarse. "
                    f"Reverse it, or pass h=... explicitly.")
            # counts increase as they refine; tolerances decrease
            self.h = list(1.0 / v) if increasing else list(v)
        if not all(self.h[i + 1] < self.h[i] for i in range(len(self.h) - 1)):
            raise ValueError(
                f"parameter {name!r}: values must run coarse -> fine, so the refinement "
                f"scalar decreases. Got h = {self.h}. Pass h=... if the default "
                f"(value, or 1/value for increasing values) is wrong.")

    def __len__(self):
        return len(self.values)

    def __repr__(self):
        return f"Parameter({self.name!r}, {self.values!r}, kind={self.kind!r})"


@dataclass
class TuneResult:
    """What the tool returns."""
    settings: dict
    certified: bool
    error_estimate: float
    target: float
    observable: str
    limiting: str | None
    budget_shares: dict
    wasted: list
    frozen: list
    violations: list
    runs_used: int
    calibration_cost: float
    recommended_cost: float
    start: dict
    trace: list = field(default_factory=list)

    def __str__(self):
        w = 74
        L = ["=" * w,
             f"errbudget: recommendation for '{self.observable}'",
             "=" * w]
        if self.certified:
            L.append(f"  Target {self.target:.3g} MET. Estimated error "
                     f"{abs(self.error_estimate):.3g}.")
        else:
            L.append(f"  Target {self.target:.3g} NOT REACHABLE with the values given.")
            L.append(f"  Best estimate at the returned setting: "
                     f"{abs(self.error_estimate):.3g}.")
        L.append("")
        L.append("  setting:")
        for k, v in self.settings.items():
            was = self.start.get(k)
            mark = "" if was == v else f"   (was {was})"
            L.append(f"    {k:<20} {v!r}{mark}")
        if self.limiting:
            L.append("")
            L.append(f"  what limits the accuracy: {self.limiting}")
        if self.wasted:
            L.append("")
            L.append("  refined past the point of having any effect:")
            for s in self.wasted:
                L.append(f"    {s}")
        if self.violations:
            L.append("")
            L.append("  declarations contradicted by measurement:")
            for s in self.violations:
                L.append(f"    {s} was declared 'cost-only' but moves the observable")
        L.append("")
        L.append(f"  {self.runs_used} solver runs used; calibration cost "
                 f"{self.calibration_cost / max(self.recommended_cost, 1e-300):.0f}x "
                 f"one run at the recommended setting")
        L.append("=" * w)
        return "\n".join(L)


def tune(solver: Callable, params: Sequence[Parameter], target: float,
         observable: str | None = None, start: dict | None = None,
         relative: bool = False, screen: bool = False, verbose: bool = False,
         floor_probe: Sequence[dict] | None = None,
         accept_unconverged: bool = False,
         **protocol_kwargs) -> TuneResult:
    """Choose parameter settings that meet `target` at close to minimal cost.

    Parameters
    ----------
    solver : callable taking the parameter names as keywords.  It may return
        a scalar, a dict of named observables, or a ``(value, cost)`` pair.  If
        no cost is returned, wall-clock time is measured.
    params : the parameters to tune.
    target : allowed error in the observable.  Absolute unless ``relative``.
    observable : which key of the returned dict to target.  Required if the
        solver returns a dict with more than one entry.
    start : the setting you would otherwise have used.  Defaults to the middle
        of each ladder.  The method calibrates locally, so this matters: it is
        where measurement begins, not merely a seed.
    relative : interpret ``target`` as a fraction of the observable's magnitude.
    screen : optional preliminary sweep that freezes parameters which
        demonstrably do not matter, using cheap runs.  Lowers calibration cost
        (2-3x in our tests) at the price of a more expensive recommendation,
        because frozen parameters stay at the user's value.  Off by default.
    accept_unconverged : by default a run is treated as failed if the solver
        reports that it did not converge, by returning a dict with
        ``"converged": False`` or by raising ``RunFailed``.  Set True to use such
        runs as answers (the behaviour of version 1.1).
    The remaining keyword arguments override the default configuration of
    version 1.2.1 (``errbudget.V121``: residual-gated safety factor, joint-refinement
    check, local-order check, no credit for cancellation between parameters).
    floor_probe : two or more settings (dicts of keyword arguments for
        ``solver``) refined beyond every candidate value, or differing only in
        something with no discretisation effect such as the thread count.  The
        spread of the observable over them is the round-off floor; targets
        below ``floor_factor`` (default 10) times it are declined without
        searching.  Recommended for single-precision or otherwise noisy codes.
    """
    names = [k.name for k in params]
    ladders = [k.values for k in params]
    hs = [k.h for k in params]
    kinds = [k.kind for k in params]
    models = [getattr(k, "model", "power") for k in params]

    if start is None:
        start_idx = tuple(len(k) // 2 for k in params)
    else:
        missing = set(start) - set(names)
        if missing:
            raise ValueError(f"start mentions unknown parameters: {sorted(missing)}")
        start_idx = []
        for k in params:
            v = start.get(k.name, k.values[len(k) // 2])
            if v not in k.values:
                raise ValueError(f"start value {v!r} for {k.name!r} is not among "
                                 f"its candidate values {k.values!r}")
            start_idx.append(k.values.index(v))
        start_idx = tuple(start_idx)

    probe = {"n": 0}

    def lookup(idx):
        kw = {names[i]: ladders[i][j] for i, j in enumerate(idx)}
        t0 = time.perf_counter()
        try:
            out = solver(**kw)
        except RunFailed:
            probe["n"] += 1
            if verbose:
                print(f"  run {probe['n']:3d}  {kw}  -> failed")
            return float("nan"), max(time.perf_counter() - t0, 1e-9)
        wall = time.perf_counter() - t0
        probe["n"] += 1
        cost = None
        if isinstance(out, tuple) and len(out) == 2:
            out, cost = out
        if isinstance(out, dict) and "converged" in out:
            out = dict(out)
            ok = bool(out.pop("converged"))
            if not ok and not accept_unconverged:
                if verbose:
                    print(f"  run {probe['n']:3d}  {kw}  -> not converged (treated as failed)")
                return float("nan"), float(cost if cost is not None else max(wall, 1e-9))
        if isinstance(out, dict):
            if observable is None:
                if len(out) != 1:
                    raise ValueError("solver returned several observables; say "
                                     "which one with observable=...")
                val = float(next(iter(out.values())))
            else:
                if observable not in out:
                    raise KeyError(f"observable {observable!r} not in solver output "
                                   f"{sorted(out)}")
                val = float(out[observable])
        else:
            val = float(out)
        if verbose:
            print(f"  run {probe['n']:3d}  {kw}  -> {val:.8g}")
        return val, float(cost if cost is not None else max(wall, 1e-9))

    if relative:
        q0, _ = lookup(start_idx)
        eps = abs(target * q0)
    else:
        eps = float(target)

    floor = None
    if floor_probe:
        from ._core import measure_floor
        vals = []
        for th in floor_probe:
            out = solver(**th)
            if isinstance(out, tuple) and len(out) == 2:
                out = out[0]
            if isinstance(out, dict):
                out = out[observable] if observable is not None else next(iter(out.values()))
            vals.append(float(out))
        floor = measure_floor(vals)
        if verbose:
            print(f"  round-off floor: {floor:.3g} (from {len(vals)} runs)")

    cfg = dict(V121)
    cfg.update(protocol_kwargs)
    r = run_protocol(lookup, ladders, start_idx, eps, hs=hs, kinds=kinds,
                     use_screen=screen, models=models, noise_floor=floor, **cfg)
    if r["budget"] is None:          # declined below the round-off floor
        return TuneResult(settings={names[i]: ladders[i][j] for i, j in enumerate(start_idx)},
                          certified=False, error_estimate=float("nan"), target=eps,
                          observable=observable or "", limiting=None,
                          budget_shares={}, wasted=[], frozen=[], violations=[],
                          runs_used=probe["n"], calibration_cost=0.0,
                          recommended_cost=1.0,
                          start={names[i]: ladders[i][j] for i, j in enumerate(start_idx)})

    idx = r["proposal"]
    lb = r["budget"]
    settings = {names[i]: ladders[i][j] for i, j in enumerate(idx)}

    # which parameter holds the error, and which are being refined for nothing
    rk = getattr(lb, "rk", [np.nan] * len(params))
    finite = [abs(v) for v in rk if np.isfinite(v)]
    tot = sum(finite) or 1.0
    shares = {names[i]: (abs(rk[i]) / tot if np.isfinite(rk[i]) else 0.0)
              for i in range(len(params))}
    limiting = max(shares, key=shares.get) if shares else None
    if limiting is not None and shares[limiting] < 0.5:
        limiting = None

    # "Refined for nothing" must mean the same thing the recommendation means:
    # how far this parameter could be coarsened, with the others held at the
    # recommended setting, before the error estimate leaves the target.  A
    # criterion relative to the parameter's own range instead (a natural choice)
    # reports parameters as wasted that the method then declines to loosen,
    # which is advice contradicting itself.
    wasted = []
    for k in range(len(params)):
        coarsest = idx[k]
        for j in range(idx[k] - 1, -1, -1):
            if j not in lb.sampled[k]:
                break
            trial = tuple(j if i == k else idx[i] for i in range(len(params)))
            if abs(lb.model_error(trial)) <= eps:
                coarsest = j
            else:
                break
        if coarsest < idx[k]:
            wasted.append(f"{names[k]:<20} {idx[k] - coarsest} step(s) finer than "
                          f"needed (={ladders[k][coarsest]!r} would still meet the "
                          f"target)")

    return TuneResult(
        settings=settings, certified=bool(r["certified"]),
        error_estimate=float(lb.E0), target=eps,
        observable=observable or "observable", limiting=limiting,
        budget_shares=shares, wasted=wasted,
        frozen=[names[k] for k in r["frozen"]],
        violations=[names[k] for k in r["violations"]],
        runs_used=r["runs"], calibration_cost=r["protocol_cost"],
        recommended_cost=float(np.exp(lb.model_logcost(idx))),
        start={names[i]: ladders[i][j] for i, j in enumerate(start_idx)},
        trace=[{names[i]: ladders[i][j] for i, j in enumerate(t)}
               for t in r.get("trace", [])],
    )
