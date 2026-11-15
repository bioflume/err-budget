# Protocol: errbudget against probabilistic Richardson extrapolation on OpenFOAM 2D-1

Fixed on 2026-11-09, before any PRE result on the 2D-1 study was graded against the exact values. This file
is committed before `arm1.py` and `arm2.py` are run. Any later change is listed as a deviation in REPORT.md.

## 0. Material

- **Study.** The OpenFOAM DFG 2D-1 study of the article, read from the article's repository at commit
  4d32148 (read-only) through its own loader, `paper/cases.py::_dfg1`. It has 840 settings: n_arc (7 values),
  residualControl (6), GAMG tolerance (5) and nCellsInCoarsestLevel (4, cost-only). 780 runs completed and
  608 converged. The cost of a run is mesh cells × outer iterations. The exact values are those of
  `openfoam_cylinder.REF`.
- **Test cases.** The seven graded cases: c_D at 1e-2, 3e-3 and 1e-3; Δp at 3e-2 and 1e-2; c_L at 1e-1 and
  3e-2. Every target is relative to the exact value: ε = target × |Q*|.
- **errbudget.** Version 1.2.1 (`pooled.CONFIGS["v1.2.1"]`), replayed with `run_protocol` exactly as
  `paper/exp72_baselines_pooled.py::cell` does.
- **PRE.** The authors' code, github.com/NewcastleRSE/sparse-probabilistic-richardson-extrapolation at commit
  3918107 (MIT licence; cloned into `ext/`, not committed). It implements sparse PRE (SPRE), the
  Gauss–Richardson extrapolation of the JRSS-B paper (GRE) and multivariate Richardson extrapolation.
  - **Method used: SPRE with the white-noise kernel.** This is the authors' default
    (`extrapolation(..., {"name": "SPRE", "k_name": "white"})`). They found it, with Matérn-1/2, "more
    accurate and robust" in their case studies. They also argue that GRE needs a number of runs that grows
    super-exponentially with the number of parameters, here 3, and their own comparison at d = 3 shows GRE
    failing to accelerate convergence. SPRE is the same authors' current method in the same probabilistic
    framework: a Gaussian process for x ↦ Q(x) with a posterior for Q* = Q(0).
  - **Implementation for the bulk runs.** `spre_white.py` is a closed-form NumPy copy of the authors' SPRE
    with the white-noise kernel: the same model, the same LOOCV criterion, the same stepwise selection of the
    index set A and the same normalisation. The authors' code takes about 75 s per fit and the closed form
    about 2 ms, and the study needs 4,256 fits in Arm 1 alone.
  - **Validation, done before this protocol (`validate_spre.py`, `results/validation.json`).**
    V0: the authors' 15 unit tests against their MATLAB implementation pass.
    V1: a fresh run of the authors' `run_cubature_experiments.py` reproduces the stored results from which
    Figs. 1 and 2 of the sparse PRE paper are drawn. It also reproduces the convergence rates the paper
    states for SPRE with the white-noise kernel: O(h²) for s = 0 and O(h⁴) for s = 1.
    V2: `spre_white.fit` against the authors' `extrapolation()` on the Fig. 1 data and on 9 run sets of
    errbudget's 2D-1 searches, at least one per test case. 21 were planned; the run was cut because the
    authors' code took up to 10 minutes per run set.
    - In 9 of 10 data sets the two give the same A, and the posterior mean and sd agree to < 2e-8 of the sd.
    - In the tenth (c_D 3e-3, n = 55) the LOOCV criterion agrees to 6 digits for every candidate through
      order 4. The selections then part at orders 5–8, where the basis (powers up to (1/n_arc)⁷) is
      ill-conditioned: the authors' code inverts these matrices directly, while `spre_white.py` uses QR.
      The posterior means differ by 0.18 sd.
    - We accept the closed form as the same method. Because high orders are where the two can differ, Arm 1
      also reports SPRE with the stepwise selection capped at order 4 (`max_order = 4`, an option of the
      authors' code) as a sensitivity.

## 1. What PRE sees (a)

- **Inputs.** x = (1/n_arc, residualControl, GAMG tolerance). These are errbudget's own refinement scalars
  (`CASE.HS`), and x → 0 is the limit Q*. SPRE normalises each input by its range, as in the authors' code.
- **The cost-only parameter.** It is not an input of PRE, since it has no limit to extrapolate to.
  - In **Arm 2**, PRE's candidate settings have it fixed at its default, nCellsInCoarsestLevel = 10. That
    gives 210 settings, all completed and 164 converged.
  - In **Arm 1**, PRE sees the runs errbudget used, which include other values of the cost-only parameter.
    Runs that differ only in it are averaged into one data point; their observables agree to about 1e-7
    relative. The selector may return any of these runs.
- **How errbudget is treated.**
  - The main comparison uses errbudget as in the article: it searches all four parameters, including the
    cost-only one.
  - Sensitivity "errbudget-3p": the same program on the 3-parameter problem with nCellsInCoarsestLevel
    fixed at 10, from the 164 converged starts at that value. Its cost ratios are reported against both the
    study's optimum and the optimum at nCellsInCoarsestLevel = 10.

## 2. From PRE's posterior to a setting (b)

PRE returns a posterior N(m, s²) for Q*. **Selector:** among the converged runs PRE has seen, return the
cheapest θ with |Q(θ) − m| + z·s ≤ ε.
- z = 2 is primary (the 95% band); z = 1 and z = 3 are sensitivities.
- The rule follows from |Q(θ) − Q*| ≤ |Q(θ) − m| + |m − Q*|, with |m − Q*| ≤ z·s at the band's credibility.
- **Decline.** PRE declines if no seen run qualifies. It also declines if the fit is undefined: fewer than two
  distinct converged runs, or a non-finite posterior.
- No other selector is used.

## 3. Two arms with matched information (c)

**Arm 1: same data, different inference.**
- For each test case and each of the 608 converged starts (the starts of Table 7 and `exp72`), errbudget is
  replayed and every run it asks for is recorded, including failed and unconverged runs.
- PRE is fitted to the converged runs among them and the selector picks from them.
- Both methods are charged the same search cost, the summed cost of all runs errbudget used.
- K = that cost ÷ the cost of the method's own returned setting.
- Sensitivity: SPRE with the stepwise selection capped at order 4.

**Arm 2: PRE's own design under errbudget's budget.**
- **Budget.** For each test case, the budget B is the median over the 608 searches of the summed cost of
  the runs errbudget used.
- **Design.** PRE spends B with the sequential cost-constrained design of the sparse PRE paper (Section 3.3,
  eq. 11), as in the authors' `run_sequential_design.py`, adapted to the discrete candidate set of §1:
  1. Initial design: an 8-point Latin hypercube over the level indices of the three parameters (with
     distinct settings), of which the 4 cheapest are run. The authors use the cheaper half of an 8-point
     reference design.
  2. The rest of the budget is split equally over 3 iterations (the authors' setting). Each iteration:
     - fit SPRE to the converged runs so far (stepwise A, LOOCV amplitude);
     - draw 1,000 random batches. A batch adds unused candidates, sampled uniformly, while the batch stays
       within the iteration's budget, and stops at the first sample that would exceed it once the batch is
       non-empty (the authors' `random_batch`);
     - run the batch that minimises PRE's posterior variance of Q* (it depends only on the design).
       If no batch has a finite variance, run the one with the most settings.
  3. Budget left unspent in an iteration carries over to the next. This is a deviation from the authors'
     fixed per-iteration budget, which can leave a discrete design unable to buy any run.
  4. The final SPRE fit uses every converged run of the design, and the selector of §2 picks from them.
- **Cost model.** The cost of a candidate is taken as known before it is run: its stored cost. This favours
  PRE, because a user would have to predict costs. The authors did the same in the cardiac study of the
  JRSS-B paper, where they recorded all candidate times first.
- **Repetitions.** 50 seeds per test case (seeds 0–49), each giving the Latin hypercube and the stochastic
  search. Sensitivity: budget 3B.
- **errbudget.** Its Arm 2 counterpart is its own result from the 608 starts, and errbudget-3p as in §1.

## 4. Failed and unconverged runs (d)

As in the article, a run that crashed or did not meet its residualControl within 2,000 outer iterations is
a failed run for both methods:
- it is paid for (its stored cost, for crashed runs where finite);
- it carries no value, so it is not a data point for PRE;
- it can never be returned.

Unconverged starting settings are not starts (as in Table 7; Table S-3.1's 172 extra starts are declined by
errbudget at once).

## 5. Outcome measures (e)

Per test case, arm and method (PRE at z = 1, 2, 3):
- **Misses.** Returned settings whose true relative error exceeds the target.
- **Declines.** Searches that returned no setting.
- **Cost ratio.** Cost of the returned setting ÷ cost of the cheapest converged setting of the study that
  meets the target. Median and worst over returned settings that meet the target (as `exp72`); the median
  over all returned settings is also given.
- **Search cost K.** As in §3; median over returned settings that meet the target (as `exp72`).
- **PRE's estimate of Q*.**
  - |m − Q*|/ε: median, 90th percentile and maximum;
  - coverage, the fraction of fits with |m − Q*| ≤ z·s, for z = 1, 2, 3;
  - the fraction with 2s ≤ ε, i.e. PRE's own claim that it resolves Q* to the target;
  - the fraction with |m − Q*| ≤ ε, i.e. m itself meets the target.

**Check of the replay.** errbudget's results are recomputed by the same code that grades PRE and compared with
Table S-3.1 (`paper/tab_pooled.tex` at 79edc9a) and with `data/pooled_v1.2.1.json` and
`data/baselines_v6.json`. They must match exactly:
- 780 starts, 608 accepted, 0 missed;
- the distinct settings;
- the median and worst cost ratios to three significant digits;
- the median runs;
- the K of each case.

## 6. What counts as a headline result (f), stated in advance

Judged at z = 2, in Arm 1 first (same data), with Arm 2 as the check that the finding survives when PRE
chooses its own runs.
1. **PRE is less reliable than errbudget.** In at least one test case PRE misses where errbudget does not, in
   at least 1% of the searches (≥ 7 of 608). A second form of the same finding: PRE's coverage of Q* at
   m ± 2s is below 90% in at least two test cases.
2. **PRE reaches the targets at lower cost for the same reliability.** PRE has 0 misses in all seven cases,
   declines at most 5% of searches in each, and has a median cost ratio at least 1.5 times lower than
   errbudget's in at least 5 of 7 cases. In Arm 2 the same must hold with K no larger than errbudget's.
3. **PRE is as reliable but less usable.** PRE has no misses but declines more than 10% of searches in at
   least 3 cases. This is a limitation of PRE as a selector, not a reliability result.
4. **No clear difference.** None of 1–3 holds.

A result is a headline for the paper only under 1 or 2 in Arm 1, with Arm 2 not contradicting it. Every
finding is limited to one study of one benchmark (steady 2D-1, three parameters).

## 7. Exploratory check (not part of the headline)

- **GRE.** The authors' GRE (white-noise base kernel, stepwise rate) on 1 random Arm 1 start per test case
  (7 data sets, seed 1), compared with SPRE on the same data: m, s, coverage and the selector's outcome.
- **Matérn-1/2 dropped.** A planned Matérn-1/2 SPRE check was dropped before any grading: one fit with the
  authors' code ran for over 10 minutes without finishing, and the task has a 3-hour CPU limit.
