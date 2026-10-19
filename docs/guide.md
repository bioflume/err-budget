# User guide

## The problem it solves

A solver has numerical parameters: a mesh, a time step, a cutoff, the tolerance of a
Newton or linear solver inside it, how often a preconditioner is rebuilt. You want one
observable (a drag coefficient, an energy, a flux) to an accuracy you can name, and you do
not want to pay for refinement that does not improve it. Choosing those parameters by hand
means either a convergence study for every parameter, at every setting of the others, or
guessing on the safe side and overpaying, often by a factor of ten.

`knobs` measures, at the setting you give it, how each parameter moves the observable and
the cost, fits a convergence model to each, and walks to the cheapest setting whose
estimated error is within your target. It then verifies that estimate at the proposed
setting and refines or declines if the verification fails. The number of runs grows with
the number of parameters, not with the product of their candidate values.

## What you have to provide

1. **Something that runs your solver for one setting and returns the observable.** From
   Python, a function `solver(**params) -> value` (or `{"name": value}`, or
   `(observables, cost)`). From the command line, a command with `{{placeholders}}` and a
   rule for reading the number from its output (`docs/job-files.md`).
2. **The parameters, each with candidate values from coarse to fine and a kind**:
   `discretisation`, `nested-solver` or `cost-only` (the README explains the three and why
   you declare them).
3. **The observable and the target accuracy**, absolute or relative.
4. **A starting setting**, the one you would otherwise have used. The measurement is local,
   so this is where it begins, not just a seed.

## Choosing candidate values

- List them from **coarse to fine**. Counts (mesh points, time steps, k-points) increase;
  tolerances decrease. The tool checks the direction and refuses an ambiguous list.
- Use **geometric spacing** (factors of 2 for meshes, 10 or 100 for tolerances). The error
  models are fitted in the refinement scalar, and a ladder with a constant ratio gives the
  cleanest fits.
- Give **enough room on both sides** of your starting setting: at least two coarser values
  (so that the tool can find the cheap setting) and two finer ones (so that it can verify
  the order of convergence through the current value).
- The finest value need not be affordable: it is where the extrapolation points, not
  where the tool goes.
- For a plane-wave or spectral cutoff, whose error decays exponentially rather than
  algebraically, say so (`model: exponential`); its values must increase.
- A parameter that takes a few discrete, unordered values (a scheme name) is not a knob
  in this sense; tune the others for each choice.

## Choosing a target

Ask for what you need, not for what feels safe. Every factor of ten in the target costs
runs and makes the recommendation more expensive; the point of the tool is that you no
longer need the margin you were adding by hand. If you need several observables, tune for
the most demanding one, or run a job per observable: the cost-optimal setting for one
observable can fail the same target on another computed from the same runs.

A relative target (`relative: true`) is a fraction of the observable's value at the
starting setting; for an observable near zero, use an absolute target.

## Reading the report

```
  Target 0.001 MET. Estimated error 2.3e-04.

  setting:
    n_cells              40     (was 80)
    p_tol                1e-06
    n_sweeps             1      (was 2)

  what limits the accuracy: n_cells

  refined past the point of having any effect:
    p_tol                2 step(s) finer than needed (=0.0001 would still meet the target)

  146 solver runs used; calibration cost 31x one run at the recommended setting
```

- **MET** means the estimated error at the recommended setting, including the safety
  margin, is within the target, and that this was verified by runs at and around the
  setting. **NOT REACHABLE** means no setting of the values offered can be certified: the
  report says which parameter could not be bounded, and the usual fix is to extend that
  parameter's list at the fine end.
- **what limits the accuracy** is the parameter that carries most of the error budget at
  the recommended setting. If you need a tighter target later, this is the one to refine.
- **refined past the point of having any effect** lists parameters that, by the error
  models fitted at the recommended setting, could be coarsened (with the others held
  there) and still meet the target. The search did not take those steps, because its
  step-by-step verification, which carries a safety margin, did not confirm them; read the
  list as where to look if you need a cheaper setting.
- **declarations contradicted by measurement** (reported when the preliminary sweep,
  `screen`, is on) means a parameter you called `cost-only` moves the observable. The
  recommendation is still safe; redeclare it and run again for a better one. Without the
  sweep a `cost-only` parameter is moved to a cheaper value only when the observable stays
  the same to within 5 % of the target.
- **calibration cost** is the total cost of the tuning in units of one run at the
  recommended setting. It pays back when you will run that setting many times.

The Markdown report adds the share of the error budget per parameter and the search path
(every setting the search went through).

## Failed runs, unconverged runs, round-off

- A run that crashes, diverges, exceeds the time limit or exits with an error is a failed
  run: its cost is recorded, it is never proposed, never fitted. Nothing stops.
- A run whose inner solver did not reach its tolerance is not an answer either. Tell the
  tool (`converged: False` from Python; `unconverged_if_regex` or `converged_if_regex` in a
  job file). Accepting such runs (`accept_unconverged`) makes the search cheaper and the
  certificate weaker.
- A single-precision or otherwise noisy code cannot resolve differences below its round-off
  floor, and no setting can then be certified to a target near it. Give two settings refined
  beyond every candidate value (`floor_probe`); the tool measures the floor and declines
  targets below ten times it instead of searching for them.

## When it declines

The tool declines rather than guesses when the observable converges irregularly along a
parameter (the error does not decrease between consecutive values), when a target could be
met only by a lucky cancellation between parameters, when the target is below the measured
round-off floor, or when the finest value offered is not fine enough. The report says which.
Extending the candidate list of the parameter it names, usually at the fine end, is the
remedy; so is checking that the parameter does what you declared.

## Cost

Calibration costs runs: about a hundred to a few hundred in the article's tests, many of
them at coarser and cheaper settings than your starting one, some at finer ones for the
verification of the convergence order. Wall-clock time is the default cost; a cost read from the output (CPU time, core-seconds, a
work count) is better when runs share a machine. For parallel codes, use core-seconds.

## Parallel and cluster use

The runner executes one run at a time, in the order the search needs (each setting
depends on the previous ones). Put the parallelism inside the run: `mpirun -n 64 ...` in
the command, and core-seconds as the cost. On a cluster, run `knobs run` inside one
allocation that is long enough for the whole tuning, or in a job that is resubmitted: the
run log resumes the tuning where it stopped.

## Reproducibility

Every run directory keeps `params.json`, `stdout.txt` and `stderr.txt`; `runs.jsonl` is
the complete log; `knobs_result.json` holds the result with the search path and the
package version. The engine is deterministic: the same solver outputs give the same
recommendation.
