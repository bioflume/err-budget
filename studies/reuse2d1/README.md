# Reuse and payback on a Reynolds-number study of 2D-1

The article's OpenFOAM case 2D-1 (steady laminar flow past a cylinder in a channel) at 13
Reynolds numbers, Re = 10 to 40: the whole 840-setting grid at each, and references at
n_arc 96, 128 and 192. errbudget searches at the hardest member; its setting is reused at
every member and graded against the references; the search cost gives the payback.

- `PROTOCOL.md`: the pre-registration (questions, runs, rules, outcomes).
- `RUNBOOK.md`: the commands for Hamilton.
- `cyl/case.py`: the article's case and post-processing, ported unchanged apart from the
  inflow velocity; `cyl/study.py`: members, lists, runs, tasks.
- `driver.py`, `slurm/`: one Slurm array element per mesh and member.
- `analysis/reuse.py`: the analysis on the Mac (needs the errbudget engine of the article).

## Checked before submission (14 November 2026, local OpenFOAM v2512 in Docker)

- At Re = 20 the port reproduces the article's stored runs (OpenFOAM v1912): c_D 5.49171559
  against 5.49170935 (n_arc 8, residualControl 10⁻⁴), 5.53936848 against 5.53936848 (n_arc 16,
  10⁻⁶), 5.51808664 against 5.51808678 (n_arc 12), with the same iteration counts (95, 123,
  81).
- Re = 10 and 40 converge on n_arc 8 and 32 (c_D 8.414 and 4.030 at n_arc 32, residualControl
  10⁻⁷); the reference setting at Re = 40 on n_arc 96 (112,704 cells) converged to
  residualControl 10⁻⁸ in 1,522 iterations, 140 CPU seconds.

## Cost

The article's 840-run grid at Re = 20 took 2.25 CPU hours of `simpleFoam` on the article's
laptop. Thirteen members: about 30 core-hours, plus the driver's meshing and
post-processing; the 39 references: about 5–10 (n_arc 192 has 442,368 cells, about four
times n_arc 96, and needs more iterations). Predicted 40–60 core-hours, guard 200.
