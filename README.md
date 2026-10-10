# knobs

Choose a solver's numerical parameters to hit a target accuracy **in the quantity
you actually care about**, at close to the cheapest cost that achieves it.

You give it a function that runs your solver, a list of its numerical parameters
with candidate values, the observable you want, and how accurate you need it. It
gives you a setting, an error estimate, and a statement of which parameter is
limiting the accuracy and which ones you are refining for nothing.

It needs no adjoint, no access to the solver's source, and no analysis of your
particular problem. It costs runs proportional to the *number* of parameters
rather than to their product.

```bash
pip install -e .
```

## In one example

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

That starting configuration is a reasonable, cautious guess. It costs 1000 units
of work for an error of $1.4\times10^{-7}$ when $10^{-4}$ was asked for. The
recommendation costs 64 — **15× less** — and still meets the target with room to
spare. Run it yourself: `python examples/02_bring_your_own_solver.py`.

## Declaring your parameters

Each parameter is declared as one of three kinds. This matters, and it is the one
thing the tool asks of you.

| kind | meaning | examples |
|---|---|---|
| `discretisation` | refining improves the answer and costs more | mesh spacing, time step, basis size, cutoff radius, k-point mesh |
| `nested-solver` | tolerance or iteration limit of a solver nested inside the outer problem | Newton tolerance, linear-solver relative tolerance, iteration caps, Krylov restart, constraint tolerance |
| `cost-only` | changes cost but not the answer, often with an interior cost optimum | preconditioner fill or drop tolerance, how often a Jacobian or preconditioner is refreshed, neighbour-list rebuild frequency |

**Why you declare it rather than the tool inferring it.** We tested inference and
it scores 3–4 out of 8. The obstacle is structural: a nested-solver parameter that
is *masked* — currently buried under another parameter's error — is
indistinguishable from one that does nothing at all. In one of our tests even the
mesh measured zero influence, because a loose Newton tolerance dominated so
completely that refining the mesh changed nothing.

**The tool falsifies what it cannot infer.** The only declaration that is unsafe
to get wrong is calling an influential parameter `cost-only`, and that is exactly
what a sweep detects — in our data every inert parameter moved the observable
by less than any active one. Contradicted declarations are reported. Get a
declaration wrong and you lose efficiency, not correctness.

## What it will tell you

- the recommended setting, and what it changed from yours
- whether the target was met — and, if it cannot be met with the values you
  offered, it says so instead of returning a confident wrong answer
- which parameter is limiting the accuracy
- which parameters are refined finer than they need to be, and by how much
- any declaration its own measurements contradict

## Two things worth knowing before you use it

**Calibration is not free.** It costs tens to a couple of hundred times one
production run. You cannot certify an answer without computing a better one, so
part of that is irreducible. It pays back if you are running parameter studies or
ensembles; it does not pay back for a single run.

**The certificate is measured, not proven.** You get a computed estimate with a
calibrated margin and a stopping rule that is re-verified at every step. On four
test problems with exact solutions (L2-L5, 22 cases) and on OpenFOAM, Quantum
ESPRESSO, ABINIT, GROMACS and LAMMPS, the article reports every miss and its
cause; `paper/README.md` maps each number to the script that produces it. That
is stronger than a heuristic and weaker than a proof, and we would rather say so.

**Options added in version 1.1.** `Knob(..., model="exponential")` for
spectral or plane-wave cutoffs; `tune(..., floor_probe=[...])` to measure the
round-off floor of a single-precision or noisy code and decline targets below
it; `resid_margin=5` to scale the safety factor with the quality of the fits
(recommended with the exponential model and for observables such as extremum
locations whose convergence is irregular).

## Layout

```
src/knobs/
  api.py             Knob, tune, TuneResult -- the public interface
  _core.py           measurement, the screen/walk/audit procedure, floor check
  _asymptotic.py     the balance law e_k ~ q_k/p_k (secondary result)
  problems/          test problems with exact ground truth
examples/            quickstart, bring-your-own-solver, declaration falsified
paper/               everything needed to reproduce the paper: reproduce.sh
tests/               behavioural guarantees the paper claims
DESIGN_NOTES.md      why the method is shaped this way, including what failed
```

`paper/reproduce.sh` regenerates every number and figure in the paper.

## Tests

```bash
pip install -e ".[dev]" && pytest
```

The tests encode the paper's claims, not just code health: that recommendations
meet their target from a range of starting points, that unreachable targets are
refused, that the diagnostics never contradict the recommendation, and that on a
solver whose error is imposed analytically the measurement recovers it.

## Licence

MIT.
