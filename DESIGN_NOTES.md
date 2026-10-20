# errbudget (formerly knobs) — working notes

These notes were written while the program was called knobs, and keep the names of that time.

**Question.** An external user of someone else's solver cannot choose its numerical
parameters. They lack the developer's intuition for how resolution and tolerances trade
accuracy against cost *for the observable they care about*. Can we measure that trade-off
automatically and hand it over?

## The three things we found

**Masking.** Every knob has a threshold past which refining it changes the observable by
exactly nothing, because another knob's error buries it. Cost keeps rising across that
threshold; accuracy does not. In 2D, `gmres_tol`'s residual error at a mid configuration is
**5.8e-23** — it is not a small contributor, it is absent — while N carries 100% of the
budget.

**Non-transferability.** A masking threshold measured with every other knob at its finest is
*too coarse* for `gmres_tol` in 73/125 and 88/125 of the other-knob combinations. Coarsening
the time step makes each Newton/GMRES solve harder, so the tolerance's own error grows
faster than the discretisation error that would have buried it. Discretisation knobs never
violated. **Calibration must be local.**

**Parameters belong to the (solver, observable, target) triple, not the solver.** The
configuration cost-optimal for energy at a 1e-3 relative target *fails* that target on
dissipation. "Set everything to 1e-10" overpays 10.4× on dissipation and 1.3× on energy —
same runs, same solver.

There is also a **third species of knob**. An ILU drop tolerance leaves the answer untouched
(observables agree to every printed digit across 4.5 orders of magnitude) while its cost
curve has an interior minimum: loose ⇒ iteration count explodes, tight ⇒ factorisation does.
It trades cost against cost. Users tune it as though it controlled accuracy.

## Method (`knobs/local.py`) — one anchor

At each anchor, a single **asymmetric OAT sweep**: the whole coarse side (cheap), only
`fine_reach=1` settings finer (expensive). From it:

- **sensitivities and marginal costs** per knob;
- **masking thresholds** — where a knob's curve flattens to zero;
- **per-knob residual**: fit `e_k(h) = c_k (h^{p_k} − h_a^{p_k})`; then `r_k = c_k h_a^{p_k}`
  is the error the anchor still carries on that axis. One rate per knob. The `r_k` *are* the
  budget, so the absolute estimate and the decomposition come from the same fit;
- **consistency floor**: the error is at least the largest single-knob *refinement* deviation;
- **margin** Fs = 1.25.

Then a **monotone two-phase walk**, re-verifying at every new anchor:

- **Phase 1, refine** while the anchor's own estimate misses the target, stepping the knob
  with the best error-reduction per unit log-cost. Each step increases an index → terminates.
- **Phase 2, improve** by accepting any single-knob move, in *either* direction, that
  strictly lowers modelled cost and keeps the estimate in target. Cost is a strict potential
  → terminates.

Every configuration is evaluated at most once across the whole walk.

### Seven corrections, each a real defect

1. **Biased order estimates.** The OAT curve is `c_k h^p − c_k h_ref^p`, not the error; a
   log-log fit of it is biased exactly where it is sampled. Fit the offset form.
2. **The fine reference is not transferable** (killed design v1) — see above.
3. **The consistency floor must use the refining direction only.** That coarsening a knob
   changes the answer says the coarse point is bad, not that the anchor is. Using the whole
   sweep makes the floor enormous and the protocol recommends the finest setting always.
4. **Infeasibility must not return the anchor** (silently certifies an under-resolved run) —
   nor jump to the finest, which costs 300× what the target needs.
5. **The walk must be monotone.** Jumping to the predicted optimum each iteration never
   reaches a fixed point; it stops at the iteration cap, so the answer depended on the cap.
6. **Phase 2 cannot be coarsen-only**, or it can never reach a cost minimum lying in the
   refining direction — which is where accuracy-neutral knobs keep theirs.
7. **The second anchor was scaffolding.** A jointly-refined Richardson path supplied the
   absolute error until the per-knob fit was sharpened. Dropping it leaves safety at 100%
   and penalties unchanged or better, for 1.8–3.3× less calibration — and it had caused our
   worst over-refinement (5.68×), because taking the larger of two estimators inherits the
   weaker one's pessimism.

A caution on our own instincts: three times a number contradicted theory, and three times
it was a bug in the harness, not a finding. The cancellation-free mode once measured *less*
safe than the optimistic one (23% vs 75%) — because "cannot reach this accuracy" was being
reported as "you are fine where you are."

## Results

Every starting configuration in the tensor is replayed through the protocol, so these are
distributions over plausible users, not anecdotes. Ground truth is exact in both cases
(Cole–Hopf; manufactured solution) — no fine reference stands in for truth anywhere.

**L2 — 1D viscous Burgers, 4 knobs, 625 starting points**

| observable | target | safe | penalty vs optimum (med/p90) | already accurate: cost saved | missed target: fixed |
|---|---|---|---|---|---|
| energy | 1e-2 | 100% | 1.72× / 2.20× | 13.4× (419 users) | 206/206 |
| energy | 3e-3 | 100% | 1.03× / 1.82× | 10.2× (311) | 314/314 |
| energy | 1e-3 | 100% | 1.00× / 1.64× | 6.7× (199) | 426/426 |
| dissipation | 1e-2 | 100% | 1.00× / 2.08× | 15.2× (420) | 205/205 |
| dissipation | 3e-3 | 100% | 1.03× / 1.82× | 10.3× (315) | 310/310 |
| dissipation | 1e-3 | 100% | 1.08× / 3.99× | 5.9× (211) | 414/414 |

Calibration 20–83× the user's own run.

**L3 — 2D Burgers, 5 knobs, 432 starting points** (linear solve ≈73% of cost)

| observable | target | safe | penalty (med/p90) | already accurate: saved | missed: fixed | calibration |
|---|---|---|---|---|---|---|
| energy | 1e-2 rel | 100% | 1.04× / 1.09× | 8.9× (144) | 288/288 | 148× |
| energy | 3e-3 rel | 100% | 1.00× / 1.00× | 4.3× (72) | 360/360 | 437× |
| gradient | 1e-2 rel | 100% | 1.09× / 1.09× | 8.9× (144) | 288/288 | 148× |
| gradient | 3e-3 rel | 100% | 1.00× / 1.00× | 4.3× (72) | 360/360 | 448× |

### The safety margin is a third of the literature's

The GCI margin is calibrated for *reporting* an error bar; ours drives a *decision* that is
re-verified at the next anchor, so the walk absorbs most of the risk.

| Fs cap | safety over 10 (problem, observable, target) cases | worst penalty |
|---|---|---|
| 1.0 | collapses to 22.9% in one case | — |
| **1.25** | **100% everywhere** | 5.68× (now 1.09× with the path off) |
| 3.0 (GCI default) | 100% everywhere | 8.93× |

Its price scales with dimension: one extra grid refinement costs ~2–4× in 1D, ~8× in 2D. A
conservative estimator that is nearly free in 1D is expensive in 2D and worse in 3D.

### Secondary: the balance law

Minimising `log C = const − Σ q_k log h_k` subject to `Σ c_k h_k^{p_k} = ε` gives error
contributions in ratio `e_k ∝ q_k/p_k`. The folklore "equalise all error contributions" is
the special case `q_k/p_k = const`. Verified on the synthetic ladder, but it governs only
knobs in their asymptotic regime, and the interesting knobs are masked instead.

## Scaling to real codes: screening (exp14)

The protocol is O(Σ knobs) per anchor, never O(∏), and used **48–94 distinct runs** in these
studies — the tensor sweeps exist only to grade the method, they are not part of it. The
real obstacle is not sample count but knob count: production codes have 20–50 parameters.

Screening answers that *if* a knob's irrelevance can be established cheaply. Tested by
measuring each knob's influence share at a coarse configuration and asking whether the
verdict still holds at a production-resolution one:

| knob | missed by the cheap screen (2D) | (1D) |
|---|---|---|
| N | 0/27 | 0/25 |
| steps | 0/27 | 1/25 |
| **newton_tol** | **9/27** | **6/25** |
| gmres_tol | 0/27 | 1/25 |
| ilu_drop | 0/27 | — |

Screening runs are 34–40× cheaper per run. So it works — except for the nested nonlinear
solver tolerance, which a cheap screen declares irrelevant in a third of cases when it
actually matters at production resolution.

**This is the same law as the non-transferability result, now from the opposite direction.**
Earlier, `gmres_tol`'s threshold measured at fine resolution failed when transferred to
coarse. Here `newton_tol`'s irrelevance measured at coarse resolution fails when transferred
to fine. Nested iterative-solver tolerances are the knobs whose influence does not transfer
across discretisation settings, in *either* direction, because coarsening or refining the
grid changes the conditioning of the problem they are solving. Discretisation knobs transfer
fine.

Protocol implication: **screen everything cheaply, but exempt the nested solver tolerances
and probe those at production resolution.** There are usually only two or three of them, so
the expensive part stays small even at 30 parameters.

## The exemption rule at eight knobs (exp16)

Added `gmres_restart`, `precond_freq`, `newton_max` to the 2D solver and tested the
prediction: screening misses should fall *only* on nested iterative-solver knobs.
24 parameter combinations, 328 runs, cheap (N=16, steps=10) → production (N=64, steps=40).

| knob | class | missed by the cheap screen | kept needlessly |
|---|---|---|---|
| N | discretisation | 0/24 | 9/24 |
| steps | discretisation | 0/24 | 9/24 |
| **newton_tol** | **nested solver** | **9/24** | 0/24 |
| gmres_tol | nested solver | 0/24 | 0/24 |
| ilu_drop | cost-only | 0/24 | 0/24 |
| gmres_restart | nested solver | 0/24 | 0/24 |
| precond_freq | cost-only | 0/24 | 0/24 |
| newton_max | nested solver | 0/24 | 0/24 |

Identical for both observables. **Every miss is a nested-solver knob; zero misses among
discretisation and cost-only knobs.** The class-level rule is therefore sufficient but
conservative — it exempts four knobs when one needed it here, and across all our experiments
only `newton_tol` (coarse→fine) and `gmres_tol` (fine→coarse) have actually failed. That
error asymmetry is the right one for a safety rule: it never misses a knob that matters.

Screening runs are 35× cheaper. The saving grows with knob count, because the exempt class
does not: at 8 knobs with 4 exempt, screening costs ~50% of probing everything at production
resolution; at 30 knobs with 3 exempt, ~13%, i.e. ~8× cheaper.

## Screening folded into the protocol (exp18, exp19)

The protocol is now **screen → walk → audit**, with the knob class declared per knob.

*Screen.* Measure each knob's influence at a configuration where the discretisation knobs
are coarsened by two rungs and the **nested-solver knobs are set to their finest**. Knobs
below a 5% influence share are frozen; nested-solver knobs are exempt and always enter the
search. A knob declared cost-only that visibly moves the observable is reported as a
violated declaration and forced back in.

*Audit.* Probe each frozen knob once at the answer and **add its contribution to the error
budget**, certifying only if `|E0| + Σ frozen ≤ ε`. Contributions are extrapolated, not
differenced. Failing knobs are unfrozen and the walk resumes.

*Cost-only pass.* Sweep each accuracy-neutral knob once at the answer and take its cheapest
verified setting.

### Results

| study | mode | safe | calibration / user run | distinct runs | frozen |
|---|---|---|---|---|---|
| 2D, 5 knobs, 432 starts, 1e-2 | baseline | 100% | 149× | 49 | — |
| | **screened** | **100%** | **77×** | **25** | 2 |
| 2D, 5 knobs, 3e-3 | baseline | 100% | 437× | 60 | — |
| | **screened** | **100%** | **182×** | **26** | 2 |
| 2D, 8 knobs, live solver, 20 starts | baseline | 20/20 | 335× | 100 | — |
| | **screened** | **20/20** | **108×** | **43** | 3 |

Calibration falls 2.4× at five knobs and **3.1× at eight** — the benefit grows with knob
count because the exempt class does not. The price is a more expensive recommendation:
2.4–2.7× median at five knobs, but only 1.3× at eight. So screening is a user-facing option,
not a default: it wins when few production runs follow, the full protocol when many.

### Two more defects, both the same mistake in new clothes

- **Screening was itself masked.** Nine unsafe cases, all one proposal: `N` was frozen
  because the screening configuration inherited the anchor's `newton_tol = 1e-2`, whose
  error swamped the grid — exactly what exp17 showed. Fixed by tightening the cheap knobs
  while coarsening the expensive ones, which is the masking law applied to our own procedure.
- **The audit differenced instead of extrapolating.** A one-step deviation is the first term
  of what a knob still owes, not the whole tail. Using it raw certified a run 25% over
  budget. Bounding the ratio by 1/2 gives a factor of 2, which restored 100%.

Also fixed: the cost instrumentation called `lookup` twice on a cache miss, which under a
live runner ran the solver twice per new configuration.

### Honest refusal

At targets below what the ladder can reach (best 6.66e-3 relative; targets 5e-3 and 3e-3),
**0/20 runs were certified** in both modes. The protocol distinguishes "here is your answer"
from "this target is not reachable with the parameters you gave me" rather than returning a
confident wrong answer.

### Falsification in practice

Deliberately declaring the influential `newton_tol` as cost-only: the declaration was
rejected or the knob kept active in 6/10 starts, and **0/10 answers were unsafe**. The other
four were anchors where the knob was genuinely masked, so freezing it was harmless there —
and the audit is the backstop if it stops being harmless.

## Can the knob class be inferred instead of declared? (exp17) — No.

If a tool could infer each knob's structural class from its own cheap sweep, it would need
nothing from the user but a parameter list. Tested on the 8-knob 2D problem with rules fixed
in advance, at two anchors, both observables: **3–4 correct out of 8.**

The failure is not a bad rule, it is an identifiability problem, and it is the masking law
again:

- `gmres_tol`, `gmres_restart` and `newton_max` all measure `err_range = 0.00e+00` at these
  anchors — they are masked there. **A masked nested-solver knob is observationally
  identical to a cost-only knob.** One local sweep cannot distinguish "does nothing here"
  from "does nothing ever".
- At the off-centre anchor, `N` itself measures `err_range = 0` — the *grid* looked
  irrelevant, because `newton_tol = 1e-2` dominated so completely (`err_range = 1.0e-1`)
  that refining the mesh changed nothing.

So the generic, zero-input version is not reachable this way, and the class must be declared
— by the solver's developer, who knows it, which is the right place for it anyway.

**But declarations can be falsified even when they cannot be inferred**, and the asymmetry
falls the right way. The only misdeclaration that is unsafe is calling an influential knob
cost-only, and that is exactly the one the sweep detects directly: an active knob separates
from a genuinely inert one by seven orders of magnitude here (1e-1 vs 1e-8). A knob wrongly
declared cost-only while masked does no harm *at that anchor*, and the walk re-verifies at
every new anchor — so a wrong declaration is caught precisely when it begins to matter.
Declare, then falsify; do not infer.

## What learning could buy (exp15)

Calibration cost scales almost linearly with the ladder distance between the user's starting
point and the answer, and **safety is 100% at every distance**:

| distance from optimum | 2D: median runs | median calibration | 1D: median runs | calibration |
|---|---|---|---|---|
| 0 | 9 | 1.0× | 20 | 1.0× |
| 3 | 32 | 6.0× | 34 | 2.3× |
| 6 | 59 | 17.7× | 53 | 4.3× |
| 9 | 79 | 23.2× | 83 | 8.6× |

So a perfect prior is worth **up to ~23× in 2D and ~10× in 1D, in calibration cost only, and
nothing in accuracy** — because the verification walk reaches a certified answer from every
starting point tested. That is the ceiling on any learned component, RL or otherwise.

The favourable part is that a wrong prediction costs time, never correctness: the guarantee
lives in the walk, not in the prediction. Which argues for learning the *prior* and keeping
the *certificate* — a regression from problem features to a starting configuration, trained
on calibrations users generate anyway, or at most a contextual bandit for the "which knob
next" step. Model-free RL would have to rediscover structure we already exploit in closed
form (monotone convergence, additive budgets, power laws) and would need orders of magnitude
more simulations than the 9–94 the model-based walk uses, to win a bounded prize it cannot
certify.

## Open

- **2D calibration at the tight target is still ~440× a production run.** Breakdown: sweep
  fine side 32–37%, coarse side 17–36%, anchor 1–2% (the path, formerly 26–49%, is gone).
  The remaining lever is reuse across anchors — each re-anchoring currently re-sweeps every
  knob from scratch, and a knob already shown to be masked need not be re-probed.
- Separability residual ≈0.27 rms on Burgers 1D. The additive model is genuinely violated;
  the trust region is what makes the method work anyway.
- The protocol picks the cost-minimising drop tolerance in 100% of cases on one observable
  but 65% on the other; the failures over-tighten it — the user error we set out to stop.
- Nothing yet on a code we did not write ourselves.

## Layout

```
knobs/core.py       v1 global design + balance law — kept for the asymptotic results
knobs/local.py      the method
knobs/models.py     L0 synthetic ladder, L1 advection-diffusion
knobs/burgers.py    L2 1D viscous Burgers, Cole-Hopf exact solution
knobs/burgers2d.py  L3 2D Burgers, manufactured solution, ILU drop tolerance
experiments/        exp01 synthetic · exp02 adv-diff · exp03 masking figure · exp04 budget
                    exp05 protocol from 625 starts · exp06/12 path ablation · exp07 2D
                    exp08/09 margin ablation · exp10 cost breakdown · exp11 cheap calibration
                    exp13 consolidated results
```
