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

import inspect
import json
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Sequence

import numpy as np

from . import report
from . import _core
from ._core import measure_floor, run_protocol, V11, V12, V121, V122

KINDS = ("discretisation", "nested-solver", "cost-only")
CONFIGURATIONS = {"v1.2.2": V122, "v1.2.1": V121, "v1.2": V12, "v1.1": V11}
DEFAULT_CONFIGURATION = "v1.2.2"


def full_configuration(name: str) -> dict:
    """The complete engine options of a named configuration.

    Each entry of CONFIGURATIONS lists only what differs from the engine's own
    defaults; the keys it leaves out take run_protocol's defaults, so that "v1.2"
    runs with v1.2's joint-check resolution (3) rather than inheriting another
    configuration's."""
    if name not in CONFIGURATIONS:
        raise ValueError(f"configuration must be one of {', '.join(CONFIGURATIONS)}, got {name!r}")
    sig = inspect.signature(_core.run_protocol).parameters
    keys = set().union(*CONFIGURATIONS.values())
    full = {k: sig[k].default for k in keys}
    full.update(CONFIGURATIONS[name])
    return full


DECLINE_REASONS = {
    "below_floor": "the target is below floor_factor (10) times the round-off floor that "
                   "the floor probe measured",
    "start_failed": "the run at the start setting failed or did not converge",
    "unmeasured_axis": "the error along {axes} cannot be measured: the runs that would "
                       "measure it failed",
    "unbounded_axis": "the changes along {axes} do not decrease, so the error left there "
                      "is unbounded; extend their values or start finer",
    "not_reached": "no setting within the values given meets the target with the "
                   "estimated error; extend the finest values{of}",
}


class RunFailed(Exception):
    """Raise this from the solver function to report a run that failed (crashed,
    diverged, or whose nested solver did not converge).  The program records it
    as unusable and never proposes or fits to it.  `cost`, if given, is what the
    failed run cost; otherwise the wall-clock time of the call is charged."""

    def __init__(self, *args, cost: float | None = None):
        super().__init__(*args)
        self.cost = cost


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
        if kind not in KINDS:
            raise ValueError(f"kind must be one of {KINDS}, got {kind!r}")
        if len(values) < 2:
            raise ValueError(f"parameter {name!r}: need at least two candidate values")
        self.name, self.values, self.kind, self.model = name, list(values), kind, model
        if model == "exponential":
            self.h = [float(v) for v in self.values]
            if not all(self.h[i + 1] > self.h[i] for i in range(len(self.h) - 1)):
                raise ValueError(f"parameter {name!r}: an exponential-model parameter must "
                                 f"increase from coarse to fine, got {values!r}")
            return
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
    declined: str | None = None           # why the target was not met (DECLINE_REASONS)

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
            if self.declined:
                L.append(f"  Why: {self.declined}.")
            if np.isfinite(self.error_estimate):
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
        if self.frozen:
            L.append("")
            L.append("  frozen by the preliminary sweep (kept at your value): "
                     + ", ".join(self.frozen))
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

    # ------------------------------------------------------------ persistence
    def to_dict(self) -> dict:
        """A plain dict (JSON-serialisable) with every field of the result."""
        def clean(v):
            if isinstance(v, dict):
                return {str(k): clean(x) for k, x in v.items()}
            if isinstance(v, (list, tuple)):
                return [clean(x) for x in v]
            if isinstance(v, (np.floating, float)):
                v = float(v)
                return v if np.isfinite(v) else (None if np.isnan(v) else ("inf" if v > 0 else "-inf"))
            if isinstance(v, (np.integer,)):
                return int(v)
            if isinstance(v, np.bool_):
                return bool(v)
            if isinstance(v, np.ndarray):
                return clean(v.tolist())
            return v
        d = {k: clean(getattr(self, k)) for k in self.__dataclass_fields__}
        d["errbudget_version"] = __import__("errbudget").__version__
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "TuneResult":
        """Rebuild a result written by `to_dict` (unknown keys are ignored)."""
        def unclean(v):
            if v is None:
                return float("nan")
            if v == "inf":
                return float("inf")
            if v == "-inf":
                return float("-inf")
            return v
        fields = cls.__dataclass_fields__
        kw = {k: d[k] for k in fields if k in d}
        for k in ("error_estimate", "target", "calibration_cost", "recommended_cost"):
            if k in kw:
                kw[k] = unclean(kw[k])
        return cls(**kw)

    def to_json(self, path: str, indent: int = 1) -> None:
        """Write the result as JSON; read it back with `TuneResult.load`."""
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=indent)

    @classmethod
    def load(cls, path: str) -> "TuneResult":
        with open(path) as f:
            return cls.from_dict(json.load(f))

    def markdown(self) -> str:
        """The same report as `str(result)`, as Markdown with a settings table."""
        return report.markdown(self)


def tune(solver: Callable, params: Sequence[Parameter], target: float,
         observable: str | None = None, start: dict | None = None,
         relative: bool = False, screen: bool = False, verbose: bool = False,
         floor_probe: Sequence[dict] | None = None,
         accept_unconverged: bool = False,
         configuration: str = DEFAULT_CONFIGURATION,
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
        The target is first scaled by the value at the start setting; if the
        search ends at a setting whose value, less its estimated error, is smaller
        in magnitude, the search is repeated with the smaller absolute target
        (repeated runs come from the cache), so the target holds relative to the
        converged value as far as the estimate can tell.
    screen : optional preliminary sweep that freezes parameters which
        demonstrably do not matter, using cheap runs.  Lowers calibration cost
        (2-3x in our tests) at the price of a more expensive recommendation,
        because frozen parameters stay at the user's value.  Off by default.
    accept_unconverged : by default a run is treated as failed if the solver
        reports that it did not converge, by returning a dict with
        ``"converged": False`` or by raising ``RunFailed``.  Set True to use such
        runs as answers (the behaviour of version 1.1).
    configuration : the engine's rules, by name: "v1.2.2" (default; ``errbudget.V122``),
        "v1.2.1" (the article's program, ``V121``), "v1.2" or "v1.1".
    The remaining keyword arguments override single options of that configuration
    (they are added to it, not used instead of it: to run version 1.2, pass
    ``configuration="v1.2"``, not ``**V12``).
    floor_probe : two or more settings (dicts of keyword arguments for
        ``solver``) refined beyond every candidate value, or differing only in
        something with no discretisation effect such as the thread count.  The
        spread of the observable over them is the round-off floor; targets
        below ``floor_factor`` (default 10) times it are declined without
        searching.  Probe runs are read like every other run, so an unconverged
        one is left out.  Recommended for single-precision or otherwise noisy codes.
    """
    names = [k.name for k in params]
    ladders = [k.values for k in params]
    hs = [k.h for k in params]
    kinds = [k.kind for k in params]
    models = [k.model for k in params]

    def setting(idx):
        return {names[i]: ladders[i][j] for i, j in enumerate(idx)}

    def unflag(out):
        """A run's output without its "converged" flag, and whether the run is used:
        an unconverged run counts as failed unless accept_unconverged."""
        if isinstance(out, dict) and "converged" in out:
            out = dict(out)
            return out, bool(out.pop("converged")) or accept_unconverged
        return out, True

    def value(out):
        """The observable in a run's output: a scalar, or a dict of observables."""
        if not isinstance(out, dict):
            return float(out)
        if observable is None:
            if len(out) != 1:
                raise ValueError("solver returned several observables; say "
                                 "which one with observable=...")
            return float(next(iter(out.values())))
        if observable not in out:
            raise KeyError(f"observable {observable!r} not in solver output {sorted(out)}")
        return float(out[observable])

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
    memo: dict = {}

    def lookup(idx):
        # every setting is run once, whichever part of tune asks for it first
        idx = tuple(idx)
        if idx not in memo:
            memo[idx] = run(idx)
        return memo[idx]

    def run(idx):
        kw = setting(idx)
        t0 = time.perf_counter()
        try:
            out = solver(**kw)
        except RunFailed as exc:
            probe["n"] += 1
            if verbose:
                print(f"  run {probe['n']:3d}  {kw}  -> failed")
            cost = getattr(exc, "cost", None)
            return float("nan"), float(cost if cost is not None
                                       else max(time.perf_counter() - t0, 1e-9))
        wall = time.perf_counter() - t0
        probe["n"] += 1
        cost = None
        if isinstance(out, tuple) and len(out) == 2:
            out, cost = out
        out, used = unflag(out)
        if not used:
            if verbose:
                print(f"  run {probe['n']:3d}  {kw}  -> not converged (treated as failed)")
            return float("nan"), float(cost if cost is not None else max(wall, 1e-9))
        val = value(out)
        if verbose:
            print(f"  run {probe['n']:3d}  {kw}  -> {val:.8g}")
        return val, float(cost if cost is not None else max(wall, 1e-9))

    if relative:
        q0, _ = lookup(start_idx)
        if not np.isfinite(q0) or q0 == 0:
            raise ValueError(
                f"relative target: the run at the start setting {setting(start_idx)} "
                f"{'failed or did not converge' if not np.isfinite(q0) else 'gave 0'}, "
                f"so the target cannot be scaled by the observable; give an absolute "
                f"target or another start")
        eps = abs(target * q0)
    else:
        eps = float(target)
    if not (np.isfinite(eps) and eps > 0):
        raise ValueError(f"target must be a positive finite number, got {target!r}")

    floor = None
    if floor_probe:
        vals = []
        for th in floor_probe:
            unknown, missing = set(th) - set(names), set(names) - set(th)
            if unknown or missing:
                raise ValueError(f"floor_probe setting {th}: give every parameter "
                                 f"({', '.join(names)}) and no others")
            try:
                out = solver(**th)
            except RunFailed:
                vals.append(float("nan"))
                continue
            if isinstance(out, tuple) and len(out) == 2:
                out = out[0]
            out, used = unflag(out)
            vals.append(value(out) if used else float("nan"))
        if sum(np.isfinite(vals)) < 2:
            raise ValueError(f"floor_probe: {sum(np.isfinite(vals))} of {len(vals)} probe "
                             f"runs gave a value (the others failed or did not converge); "
                             f"the round-off floor needs two")
        floor = measure_floor(vals)
        if verbose:
            print(f"  round-off floor: {floor:.3g} (from {len(vals)} runs)")

    cfg = full_configuration(configuration)
    cfg.update(protocol_kwargs)

    def search(eps):
        return run_protocol(lookup, ladders, start_idx, eps, hs=hs, kinds=kinds,
                            use_screen=screen, models=models, noise_floor=floor, **cfg)

    r = search(eps)
    if relative:
        # the target is relative to the observable, which the start value only
        # approximates: tighten it until it holds at the returned setting's value
        # less its estimated error (a few repeats; every run comes from the cache)
        for _ in range(5):
            if not r["certified"] or r["budget"] is None:
                break
            q_fin = memo[tuple(r["proposal"])][0]
            lower = abs(q_fin) - abs(r["budget"].E0)
            eps_needed = abs(target) * lower
            if eps_needed >= eps * (1 - 1e-12):
                break
            if not eps_needed > 0:
                r = dict(r, certified=False, declined="not_reached")
                break
            eps = eps_needed
            r = search(eps)

    if r["budget"] is None:          # declined below the round-off floor
        return TuneResult(settings=setting(start_idx),
                          certified=False, error_estimate=float("nan"), target=eps,
                          observable=observable or "", limiting=None,
                          budget_shares={}, wasted=[], frozen=[], violations=[],
                          runs_used=probe["n"], calibration_cost=0.0,
                          recommended_cost=1.0, start=setting(start_idx),
                          declined=DECLINE_REASONS["below_floor"])

    idx = r["proposal"]
    lb = r["budget"]
    settings = setting(idx)

    # which parameter holds the error, and which are being refined for nothing
    rk = lb.rk
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

    declined = None
    if not r["certified"]:
        why = r.get("declined") or "not_reached"
        axes = [names[k] for k in range(len(params))
                if why == "unbounded_axis" and not np.isfinite(lb.rk[k])] if lb is not None else []
        if why == "unmeasured_axis":
            axes = [names[k] for k in getattr(lb, "unmeasured", [])]
        of = f" of {limiting}" if limiting else ""
        declined = DECLINE_REASONS.get(why, why).format(axes=", ".join(axes) or "an axis", of=of)

    return TuneResult(
        settings=settings, certified=bool(r["certified"]),
        error_estimate=float(lb.E0), target=eps,
        observable=observable or "observable", limiting=limiting,
        budget_shares=shares, wasted=wasted,
        frozen=[names[k] for k in r["frozen"]],
        violations=[names[k] for k in r["violations"]],
        runs_used=r["runs"], calibration_cost=r["protocol_cost"],
        recommended_cost=float(np.exp(lb.model_logcost(idx))),
        start=setting(start_idx),
        trace=[setting(t) for t in r["trace"]],
        declined=declined,
    )
