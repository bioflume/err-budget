# Cost estimate (before Stage 0; amended with A4 and A5, 13 November 2026)

Agreed budget about 600 core-hours, ceiling 1,000 (PROTOCOL.md §12). Already used on Hamilton
before Stage 0: 215.1 core-hours (check jobs, diagnostics, numerics sweep, convergence map,
confirmation map, settling measurement). With the measured speed of the mesh runs (4.6 µs per
cell and iteration) applied to every run, Stages 0–2 come to about **464 core-hours central and
928 with a factor-2 margin**. The reference-mesh runs measured 2.1–2.9 µs on Hamilton; at that
speed for them, Stage 2's setting (d) falls from about 242 to about 150 core-hours and the total
to about 370. Stage 0 replaces these assumptions by measurements; `analysis/stage0.py` prints
the new prediction for Stages 1–2 and stops the study if it passes 900.

## Assumptions

| Quantity | Value | Basis |
|---|---|---|
| Cells per run | exact | from the mesh generator (`af/costmodel.py`) |
| Core-seconds per cell and iteration | 4.6 µs | **measured** (A4): median of the 180 mesh runs of the convergence map on Hamilton (IQR 3.5–6.6); first assumed 3 µs |
| Iterations to meet residualControl, n_surf 128 | 300, 700, 1,200, 1,800, 4,800 for 10⁻³ … 10⁻⁷ | 10⁻⁶: the convergence map's converged runs gave about 2,080 at n_surf 128; 10⁻⁷: calibrated on the one reference run (10,316 iterations at n_surf 362); the rest as first assumed |
| Growth of iterations per mesh level | ×1.29 | **measured** in the convergence map (first assumed ×1.25) |
| Reference runs (A4, A5) | settled at R = 5000 (about 11,000 iterations), + the first run at n_surf 91 (in the table) and, per design, the far-field run at R = 500 (in the table) | the settling measurement: L/D settled to 10⁻⁴ from iteration 8,626 |
| Cores | 1 per run | serial is the cheapest per core-hour; Stage 0 measures a 4-core run for the record |

Reproduce: `python3 analysis/cost_estimate.py` (options `--t`, `--hardest`).

## Per stage (hardest design taken as m 4 %, t 12 %; the result barely depends on it)

| item | runs | cells | core-h (central) | core-h (x2) | longest run (h) |
|---|---|---|---|---|---|
| Stage 0 (check job: 2 tiny runs) | 2 | 13,664-13,664 | 0.0 | 0.0 | 0.00 |
| Stage 0 | 26 | 18,000-251,430 | 4.4 | 8.8 | 0.96 |
| Stage 1, n_surf 64 (108 runs) | 108 | 12,768-23,290 | 1.5 | 3.0 | 0.03 |
| Stage 1, n_surf 91 (108 runs) | 108 | 24,964-45,396 | 3.8 | 7.5 | 0.08 |
| Stage 1, n_surf 128 (108 runs) | 108 | 47,952-87,906 | 9.3 | 18.6 | 0.20 |
| Stage 1, n_surf 181 (108 runs) | 108 | 92,022-171,264 | 23.2 | 46.4 | 0.51 |
| Stage 1, n_surf 256 (108 runs) | 108 | 177,684-336,192 | 58.4 | 116.8 | 1.29 |
| Stage 1, reference runs (6) | 6 | 34,848-1,336,720 | 54.5 | 109.0 | 24.12 |
| Stage 2, (a) errbudget, if = default start | 25 | 66,690-66,690 | 2.6 | 5.1 | 0.10 |
| Stage 2, (a) errbudget, if = cautious | 25 | 336,192-336,192 | 32.2 | 64.3 | 1.29 |
| Stage 2, (b) tutorial | 25 | 66,690-66,690 | 2.6 | 5.1 | 0.10 |
| Stage 2, (c) cautious | 25 | 336,192-336,192 | 32.2 | 64.3 | 1.29 |
| Stage 2, (d) reference | 50 | 34,848-682,630 | 242.2 | 484.5 | 9.55 |

Stage 0 4.4, Stage 1 150.6, Stage 2 309.1 (with (a) = cautious, the upper bound); total 464 core-hours central, 928 with the x2 margin

Stage 2 (d), the settled reference of all 25 designs with its far-field runs, is about half of it.

## What could change it

- **Time per cell and iteration.** At 6 µs every number doubles: 428 core-hours central,
  856 with the margin, still under the ceiling.
- **Slow or stalled convergence.** A run that stalls stops 3,000 iterations after its last
  new minimum of the residuals; one that crawls stops at the cap of 20,000 iterations
  (40,000 for the reference). If every n_surf 256 run at residualControl 10⁻⁶ (27 runs) ran
  to the cap, Stage 1 would gain about 120 core-hours. Stage 0 run A′ measures the slowest
  setting before Stage 1 is submitted.
- **The far-field amendment** of PROTOCOL.md §3.2 is withdrawn (A5): the radius list never changes.
- **The fallback numerics** were withdrawn by amendment A3 (PROTOCOL.md §14): if the
  numerics check of §3.3 triggers, the study stops for an amendment.
- **Amendment A3** (free stream, first-order convection of k and ω). Local runs on the
  coarsest mesh needed 570–1,900 iterations to reach residualControl 10⁻⁶, in line with
  the model (1,800 at n_surf 128, ×1.25 per level); on some meshes of the thick cambered
  sections the iteration levels off and stops after 4,000–6,000 iterations. The convergence
  map before Stage 0 (182 runs) costs about 45 core-hours by the model, at most about 110
  if its runs stall or the reference runs to its cap.

## Guards (PROTOCOL.md §12)

Stage 0: 8 core-hours. Stage 1 with its reference runs: 400. Stage 2: 500. `slurm/status.sh`
prints the core-hours used per stage from `sacct` and says GUARD PASSED when a stage passes
its guard; stop submitting that stage and report.

## Wall time and queue

All runs are serial tasks on `shared`. With the central estimate the longest Stage-1 run
takes under an hour, and the reference runs 4–9 hours (n_surf 362 and 512). Time limits are
set per mesh size at two to three times the predicted run, between 3 hours (n_surf ≤ 128)
and 3 days (n_surf 512); a run that reaches its limit writes its fields and resumes at the
next submission. Stage 0 should finish within a few hours of starting, Stage 1 within a day
of starting (the reference runs set the pace), Stage 2 within a day.
