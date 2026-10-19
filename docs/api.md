# Python API

```python
from knobs import Knob, tune, TuneResult, RunFailed
```

## `Knob(name, values, kind="discretisation", h=None, model="power")`

One numerical parameter.

| argument | |
|---|---|
| `name` | the keyword the solver function expects |
| `values` | candidate settings, ordered coarse (cheap, inaccurate) to fine; at least two |
| `kind` | `"discretisation"`, `"nested-solver"` or `"cost-only"` (`knobs.KINDS`) |
| `h` | optional map from a value to a refinement scalar that decreases as the setting gets finer. Default: the value itself for tolerances, `1/value` for counts (inferred when the values increase). An increasing list of values below 1 is refused as ambiguous: reverse it or give `h` |
| `model` | `"power"` fits `e(h) = c (h^p - h_a^p)`; `"exponential"` fits `e(x) = c (exp(-b x) - exp(-b x_a))` in the setting `x` itself, for spectral or plane-wave cutoffs (values must increase) |

## `tune(solver, knobs, target, observable=None, start=None, relative=False, screen=False, verbose=False, floor_probe=None, accept_unconverged=False, **engine)`

| argument | |
|---|---|
| `solver` | a callable taking the knob names as keywords. It may return a number, a dict of named observables, or `(observables, cost)`. If no cost is returned, wall-clock time is used. A dict may carry `"converged": False`; the function may raise `RunFailed` |
| `knobs` | the `Knob`s |
| `target` | the allowed error in the observable, absolute unless `relative` |
| `observable` | which key of the returned dict to target; required when the dict has several entries |
| `start` | the setting you would otherwise have used, as a dict; defaults to the middle of each list. The measurement is local, so this is where it begins |
| `relative` | interpret `target` as a fraction of the observable's magnitude at `start` |
| `screen` | run the preliminary sweep that freezes parameters which demonstrably do not matter (fewer runs, more expensive recommendations) |
| `verbose` | print every run |
| `floor_probe` | two or more settings (dicts) refined beyond every candidate value; the spread of the observable over them is the round-off floor, and targets below `floor_factor` (10) times it are declined |
| `accept_unconverged` | use runs reported as not converged as answers (version 1.1 behaviour) |
| `**engine` | overrides of the default configuration `knobs.V121` (below) |

Returns a `TuneResult`.

## `TuneResult`

| field | |
|---|---|
| `settings` | the recommended setting, `{name: value}` |
| `certified` | whether the target was met |
| `error_estimate` | the estimated error at `settings` (with the safety margin) |
| `target` | the absolute target used |
| `observable` | its name |
| `limiting` | the parameter that carries most of the error budget, or `None` if no single one does |
| `budget_shares` | `{name: share}` of the error budget |
| `wasted` | parameters that, by the fitted models, could be coarsened and still meet the target |
| `frozen` | parameters frozen by the preliminary sweep |
| `violations` | parameters declared `cost-only` that move the observable |
| `runs_used`, `calibration_cost`, `recommended_cost` | the count and total cost of the runs, and the cost of one run at `settings` |
| `start` | the starting setting |
| `trace` | the settings the search went through |

Methods: `str(result)` (the text report), `result.markdown()`, `result.to_dict()`,
`result.to_json(path)`, `TuneResult.load(path)`, `TuneResult.from_dict(d)`.

## `RunFailed`

Raise it from the solver to report a run that crashed, diverged or whose nested solver did
not converge. The run is recorded with its cost and never proposed or fitted.

## Configurations

`knobs.V121` (default), `knobs.V12`, `knobs.V11` are dicts of engine options;
`knobs.CONFIGURATIONS` maps `"v1.2.1"`, `"v1.2"`, `"v1.1"` to them. The options that
distinguish them:

| option | v1.1 | v1.2 | v1.2.1 | meaning |
|---|---|---|---|---|
| `resid_margin` | 0 | 5 | 5 | the safety factor grows with the misfit of the error models |
| `joint_check` | off | on | on | one extra run with every parameter refined together: interaction term added to the estimate |
| `joint_res` | — | 3 | 1 | which parameters enter the joint check: one-step change above `joint_res` × the resolution |
| `local_order_check` | off | on | on | three-point order check along every discretisation parameter |
| `res_frac` | 0 | 0.1 | 0.1 | differences below `res_frac` × target are treated as zero |
| `cancel_rule` | credit | none | none | whether cancellation between parameters is credited |

Other engine options (`max_iter`, `fs_cap`, `floor_factor`, `unfit`, `walk_mode`,
`phase2`, `screen_at`, `audit_factor`, `tau`, `drop`) are documented in
`knobs/_core.py`; the last four exist to reproduce the article's ablations.

## Job files from Python

```python
from knobs import load_job, run_job
job = load_job("job.yaml")           # a JobSpec; JobSpec.from_dict(d, base_dir) for a dict
result = run_job(job)                # the same as `knobs run`, returns the TuneResult
```

`knobs.runner.CommandSolver(job)` is the callable solver the job file describes, usable
with `tune` directly.

## The engine

`knobs.run_protocol(lookup, ladders, anchor, eps, hs=..., kinds=..., **config)` is the
protocol `tune` wraps: `lookup(index_tuple) -> (value, cost)`, ladders as lists of values,
the start as a tuple of indices. `knobs.LocalBudget` is the local error budget at one
anchor, `knobs.measure_floor(values)` the round-off floor, `knobs.screen(...)` the
preliminary sweep. These are what the article's replay scripts call; `tune` is the
interface for everyone else.

## Test problems

`knobs.problems` holds solvers with exact solutions, used by the tests and the examples:
`Burgers1D`, `Burgers2D`, `AdvectionDiffusion`, `ReactionDiffusion2D`, `DrivenCavity`,
`MLHybridBurgers1D`, `make_synthetic`. Each has `run(theta) -> (observables, cost, wall, info)`
and `exact_observables()`.
