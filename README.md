# knobs

Choose a solver's numerical parameters to hit a target accuracy **in the quantity you
actually care about**, at close to the cheapest cost that achieves it.

You give it something that runs your solver (a Python function or a shell command), the
numerical parameters with their candidate values, the observable you want, and how
accurate you need it. It gives you a setting, an error estimate, and a statement of which
parameter is limiting the accuracy and which ones you are refining for nothing.

It needs no adjoint, no access to the solver's source, and no analysis of your particular
problem. It costs runs proportional to the *number* of parameters rather than to their
product. It was tested on OpenFOAM, Quantum ESPRESSO, ABINIT, GROMACS and LAMMPS, driven
only through their input files, and on test problems with exact solutions; the article
reports every miss and its cause (see *Citing*).

```bash
pip install "knobs[yaml] @ git+https://github.com/bioflume/knobs"
```

## Two ways to use it

**From the command line, for any code that runs as a command** (OpenFOAM, LAMMPS, your
own binary, `python my_solver.py`):

```bash
knobs init job.yaml      # writes an annotated job file: edit the command, the parameters, the target
knobs check job.yaml     # validates it and runs your starting setting once
knobs run job.yaml       # tunes, prints the report, writes knobs_result.json and .md
```

A job file is a page of YAML (or JSON, or TOML). Here is the one from
`examples/04_command_line/`:

```yaml
command: python {{job_dir}}/solver.py --steps {{steps}} --tol {{tol}} --refresh {{refresh}}

parameters:                      # values run COARSE (cheap) to FINE (accurate)
  - {name: steps,   kind: discretisation, values: [25, 50, 100, 200, 400, 800], start: 100}
  - {name: tol,     kind: nested-solver,  values: [1e-2, 1e-4, 1e-6, 1e-8, 1e-10], start: 1e-10}
  - {name: refresh, kind: cost-only,      values: [16, 8, 4, 2, 1],               start: 1}

observable: {name: u_final, regex: "u_final\\s*=\\s*([-+0-9.eE]+)"}   # read from stdout
target: 1e-4
cost: {regex: "work\\s*=\\s*([0-9]+)"}   # optional; wall-clock time by default
```

Every run gets its own directory (template files are copied in with `{{placeholders}}`
filled), and every run is logged, so a tuning that is interrupted resumes where it was.
`docs/job-files.md` is the reference and `docs/recipes.md` has job files for OpenFOAM,
LAMMPS, GROMACS and plane-wave DFT codes.

**From Python**, for a solver you can call as a function:

```python
from knobs import Knob, tune

def solver(steps, inner_tol, inner_max, jac_every):
    ...                                  # your code, untouched
    return {"u_final": u}, work          # observables, and optionally a cost

result = tune(
    solver,
    [Knob("steps",     [25, 50, 100, 200, 400, 800], kind="discretisation"),
     Knob("inner_tol", [1e-2, 1e-4, 1e-6, 1e-8, 1e-10], kind="nested-solver"),
     Knob("inner_max", [2, 4, 8, 16],               kind="nested-solver"),
     Knob("jac_every", [16, 8, 4, 2, 1],            kind="cost-only")],
    target=1e-4, observable="u_final",
    start=dict(steps=100, inner_tol=1e-10, inner_max=16, jac_every=1))

print(result)
```

```
  Target 0.0001 MET. Estimated error 2.08e-06.

  setting:
    steps                25     (was 100)
    inner_tol            0.01   (was 1e-10)
    inner_max            2      (was 16)
    jac_every            16     (was 1)
```

That starting configuration is a reasonable, cautious guess. It costs 1000 units of work
for an error of 1.4e-7 when 1e-4 was asked for. The recommendation costs 64, **15× less**,
and still meets the target with room to spare. Run it yourself:
`python examples/02_bring_your_own_solver.py`.

## Declaring your parameters

Each parameter is declared as one of three kinds. This matters, and it is the one thing
the tool asks of you.

| kind | meaning | examples |
|---|---|---|
| `discretisation` | refining improves the answer and costs more | mesh spacing, time step, basis size, cutoff radius, k-point mesh |
| `nested-solver` | tolerance or iteration limit of a solver nested inside the outer problem | Newton tolerance, linear-solver relative tolerance, iteration caps, Krylov restart, SCF or constraint tolerance |
| `cost-only` | changes cost but not the answer, often with an interior cost optimum | preconditioner fill or drop tolerance, how often a Jacobian or preconditioner is refreshed, neighbour-list rebuild frequency |

**Why you declare it rather than the tool inferring it.** We tested inference and it scores
3–4 out of 8. The obstacle is structural: a nested-solver parameter that is *masked*,
currently buried under another parameter's error, is indistinguishable from one that does
nothing at all. In one of our tests even the mesh measured zero influence, because a loose
Newton tolerance dominated so completely that refining the mesh changed nothing.

**The tool falsifies what it cannot infer.** The only declaration that is unsafe to get
wrong is calling an influential parameter `cost-only`, and that is exactly what a sweep
detects: with `screen: true` the preliminary sweep tests every declaration and the report
lists the ones its measurements contradict. Without the sweep, a parameter declared
`cost-only` is only moved to a cheaper value when the observable stays the same to within
5 % of the target, so a wrong declaration costs efficiency, not correctness.

Values are listed from coarse to fine. For a parameter whose error decays exponentially
(a plane-wave or spectral cutoff) add `model: exponential` (`Knob(..., model="exponential")`);
its values must increase.

## What it will tell you

- the recommended setting, and what it changed from yours;
- whether the target was met, and, if it cannot be met with the values you offered, it
  says so instead of returning a confident wrong answer;
- which parameter is limiting the accuracy;
- which parameters are refined finer than they need to be, and by how much;
- any declaration its own measurements contradict (with the preliminary sweep).

## Three things worth knowing before you use it

**Calibration is not free.** In the article's tests it cost about a hundred to a few
hundred runs at the returned setting (more where convergence is irregular). You cannot
certify an answer without computing a better one, so part of that is irreducible. It pays
back if you are running parameter studies or ensembles; it does not pay back for a single
run.

**It declines rather than guesses.** Where an observable converges irregularly, or a target
can only be met through a lucky cancellation, it declines and says which parameter it could
not bound. Extend that parameter's list of values and run it again.

**The certificate is measured, not proven.** You get a computed estimate with a calibrated
margin and a stopping rule that is re-verified at every step. On four test problems with
exact solutions and on five production codes the article reports every miss and its cause.
That is stronger than a heuristic and weaker than a proof, and we would rather say so.

## Failed and unconverged runs

A run whose nested solver stopped at its iteration limit is not an answer. From Python,
return `{"drag": cd, "converged": False}` or raise `knobs.RunFailed`; from a job file, give
a pattern that marks it (`failure: {unconverged_if_regex: "did not converge"}`). A non-zero
exit code is always a failure. Such runs are recorded with their cost but never proposed
and never fitted. `accept_unconverged: true` restores the behaviour of version 1.1.

## What it checks before it accepts a setting

- the observed order along every discretisation parameter (three runs through the current
  value): a list that starts before the asymptotic range, or an error that stalls between
  two grids, is refined or declined instead of extrapolated;
- the interaction of the parameters, by refining them all together by one step (one extra
  run): errors of different parameters that cancel at one setting but not at the next are
  not credited;
- a safety factor that grows with the misfit of the error models.

The default configuration is the article's version 1.2.1 (`knobs.V121`); `v1.2` and `v1.1`
reproduce the earlier versions. For single-precision or otherwise noisy codes,
`floor_probe` measures the round-off floor from two settings refined beyond every
candidate value and declines targets below ten times it; `screen: true` runs a
preliminary sweep (fewer runs, more expensive recommendations).

## Documentation

| | |
|---|---|
| `docs/guide.md` | how it works, how to choose values and a target, how to read the report, when it declines |
| `docs/job-files.md` | the job-file reference (every field) |
| `docs/recipes.md` | job files for OpenFOAM, LAMMPS, GROMACS, Quantum ESPRESSO / ABINIT, MPI codes, Python solvers |
| `docs/api.md` | the Python API (`Knob`, `tune`, `TuneResult`, `RunFailed`, the engine) |
| `examples/` | four runnable examples, including the command-line one |
| `CHANGELOG.md` | what changed between versions |

## Install

```bash
pip install "knobs[yaml] @ git+https://github.com/bioflume/knobs"     # with YAML job files
pip install "knobs @ git+https://github.com/bioflume/knobs"           # Python API, JSON/TOML job files
```

Python 3.9 or later, NumPy and SciPy. From a clone: `pip install -e ".[dev]"` and `pytest`
(about 80 s; the tests encode the article's claims, not just code health).

## Citing

The method and its evaluation are described in

> G. Kabacaoğlu, *knobs: choosing the numerical parameters of a black-box solver for a
> target accuracy in one observable*, Computer Physics Communications, 2026.

The engine in this package is version 1.2.1, the version the article describes
(`knobs.ENGINE_VERSION`); the article's own repository, with every script, stored solver
output and the cluster jobs, is <https://github.com/bioflume/knobs-paper>. See
`CITATION.cff`.

## Licence

MIT.
