# Python API

```python
from errbudget import Parameter, tune, TuneResult, RunFailed
```

## `Parameter(name, values, kind="discretisation", h=None, model="power")`

One numerical parameter.

| argument | |
|---|---|
| `name` | the keyword the solver function expects |
| `values` | candidate settings, ordered coarse (cheap, inaccurate) to fine; at least two |
| `kind` | `"discretisation"`, `"nested-solver"` or `"cost-only"` (`errbudget.KINDS`) |
| `h` | optional map from a value to a refinement scalar that decreases as the setting gets finer. Default: the value itself for tolerances, `1/value` for counts (inferred when the values increase). An increasing list of values below 1 is refused as ambiguous: reverse it or give `h` |
| `model` | `"power"` fits `e(h) = c (h^p - h_a^p)`; `"exponential"` fits `e(x) = c (exp(-b x) - exp(-b x_a))` in the setting `x` itself, for spectral or plane-wave cutoffs (values must increase) |

## `tune(solver, params, target, observable=None, start=None, relative=False, screen=False, verbose=False, floor_probe=None, accept_unconverged=False, configuration="v1.2.2", **engine)`

| argument | |
|---|---|
| `solver` | a callable taking the parameter names as keywords. It may return a number, a dict of named observables, or `(observables, cost)`. If no cost is returned, wall-clock time is used. A dict may carry `"converged": False`; the function may raise `RunFailed` |
| `params` | the `Parameter`s to tune |
| `target` | the allowed error in the observable, absolute unless `relative` |
| `observable` | which key of the returned dict to target; required when the dict has several entries |
| `start` | the setting you would otherwise have used, as a dict; defaults to the middle of each list. The measurement is local, so this is where it begins |
| `relative` | interpret `target` as a fraction of the observable's magnitude. It is first scaled by the value at `start`; if the search ends where the value, less its estimated error, is smaller in magnitude, the search is repeated with the smaller target (from the cache). A failed or zero value at `start` is an error |
| `screen` | run the preliminary sweep that freezes parameters which demonstrably do not matter (fewer runs, more expensive recommendations) |
| `verbose` | print every run |
| `floor_probe` | two or more settings (dicts giving every parameter) refined beyond every candidate value; the spread of the observable over them is the round-off floor, and targets below `floor_factor` (10) times it are declined. A probe run that fails or does not converge is left out; fewer than two usable ones is an error |
| `accept_unconverged` | use runs reported as not converged as answers (version 1.1 behaviour) |
| `configuration` | the engine's rules by name: `"v1.2.2"` (default), `"v1.2.1"` (the article's program), `"v1.2"`, `"v1.1"` (below) |
| `**engine` | overrides of single options of that configuration; they are added to it. To run version 1.2, pass `configuration="v1.2"`, not `**V12` |

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
| `declined` | when the target was not met, why: a failed start, the round-off floor, a parameter whose error cannot be measured or is unbounded, or no setting within the values given; `None` when it was met |

Methods: `str(result)` (the text report), `result.markdown()`, `result.to_dict()`,
`result.to_json(path)`, `TuneResult.load(path)`, `TuneResult.from_dict(d)`.

## `RunFailed`

Raise it from the solver to report a run that crashed, diverged or whose nested solver did
not converge, optionally with what it cost: `RunFailed("diverged", cost=12.5)` (otherwise the
wall-clock time of the call is charged). The run is never proposed or fitted.

## Configurations

`errbudget.V122` (default), `errbudget.V121` (the article's program), `errbudget.V12` and
`errbudget.V11` are dicts of engine options; `errbudget.CONFIGURATIONS` maps `"v1.2.2"`,
`"v1.2.1"`, `"v1.2"`, `"v1.1"` to them, and `errbudget.full_configuration(name)` gives every
option of one, the engine's defaults filled in. The options that distinguish them:

| option | v1.1 | v1.2 | v1.2.1 | v1.2.2 | meaning |
|---|---|---|---|---|---|
| `resid_margin` | 0 | 5 | 5 | 5 | the safety factor grows with the misfit of the error models |
| `joint_check` | off | on | on | on | one extra run with every parameter refined together: interaction term added to the estimate |
| `joint_res` | 3 | 3 | 1 | 1 | which parameters enter the joint check: one-step change above `joint_res` × the resolution |
| `local_order_check` | off | on | on | on | three-point order check along every discretisation parameter |
| `res_frac` | 0 | 0.1 | 0.1 | 0.1 | differences below `res_frac` × target are treated as zero |
| `cancel_rule` | credit | none | none | none | whether cancellation between parameters is credited |
| `strict_unfit` | off | off | off | on | a parameter with too few successful runs to fit is charged its nearest coarser change (first order), or declined; off, it counted as error-free when it had no neighbouring run |
| `local_at_finest` | off | off | off | on | the order check also runs at the finest value of a list |
| `cancel_strict` | off | off | off | on | with `cancel_rule` none, any cancellation is summed in magnitude; off, only when the signed sum is below 0.8 times the largest term |

Other engine options (`max_iter`, `fs_cap`, `floor_factor`, `unfit`, `tau`, `drop`,
`walk_mode`, `phase2`, `screen_at`, `audit_factor`) are documented in
`errbudget/_core.py`; the last four exist only to reproduce design ablations.

## Job files from Python

```python
from errbudget import load_job, run_job
job = load_job("job.yaml")           # a JobSpec; JobSpec.from_dict(d, base_dir) for a dict
result = run_job(job)                # the same as `errbudget run`, returns the TuneResult
```

`errbudget.runner.CommandSolver(job)` is the callable solver the job file describes, usable
with `tune` directly.

## The engine

`errbudget.run_protocol(lookup, ladders, anchor, eps, hs=..., kinds=..., **config)` is the
protocol `tune` wraps: `lookup(index_tuple) -> (value, cost)`, ladders as lists of values,
the start as a tuple of indices. `errbudget.LocalBudget` is the local error budget at one
anchor, `errbudget.measure_floor(values)` the round-off floor, `errbudget.screen(...)` the
preliminary sweep. These are what the article's replay scripts call; `tune` is the
interface for everyone else.

## Test problems

`errbudget.problems` holds the article's test problems, solvers with exact solutions:
`Burgers1D`, `Burgers2D`, `ReactionDiffusion2D`, `DrivenCavity` and `MLHybridBurgers1D`.
Each has `run(theta) -> (observables, cost, wall, info)`
and `exact_observables()`.

`Burgers2D`, `ReactionDiffusion2D` and `DrivenCavity` need SymPy, which only the `dev` extra installs.
