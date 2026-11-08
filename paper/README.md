# Reproducing the article

This folder holds every script, every stored solver output and every table and figure of
the article. There are three ways to reproduce it, from cheap to expensive:

| Tier | What it does | Needs | Time |
|---|---|---|---|
| **A. Replay** | re-runs the program over the stored solver output in `data/` and regenerates every number, table and figure | Python with NumPy, SciPy, Matplotlib, SymPy, pytest | 2–4 hours on 2 cores |
| **B. From scratch** | runs every solver again and rebuilds `data/`, then replays | tier A plus LAMMPS, GROMACS, ABINIT, OpenFOAM and Quantum ESPRESSO (versions below) | about 16 CPU-hours, plus the 2D-2 build (`dfg2d2_job/`, about 40 CPU-hours) |
| **C. Cluster studies** | the five SLURM jobs J1–J5 of revision 4, including the pre-registered held-out test H3 | Durham's Hamilton (or any SLURM cluster, see `../hamilton/config.sh`) | 1–3 days of queue time, about 400–700 core-hours |

Tier A is what a reader checking the article needs. Tier B is for someone who wants to
regenerate the solver output itself. Tier C lives entirely in `../hamilton/`, whose README
is the runbook.

## Setup

```bash
git clone <this repository> errbudget-paper && cd errbudget-paper
python -m venv .venv && source .venv/bin/activate
pip install -e ".[paper]"            # numpy, scipy, matplotlib, sympy, pytest and the errbudget package
cd paper && bash run_v6_analyses.sh check
```

`check` lists what it finds: the Python stack, the `errbudget` package (it must report
`V121`), the solvers, and the git tag `v1.2.0-frozen` (or `frozen_src/`), which the
held-out replays need. Missing solvers only matter for tier B and for the three steps of
tier A that call them (A03, A04 and the second half of A12).

The article's own environment is in the table at the end of this file. Replays with
NumPy 2.x and SciPy 1.18 give identical counts (checked on 2026-10-19).

## Tier A: replay the stored output

### Revision 4 (current): `run_v6_analyses.sh`

```bash
bash run_v6_analyses.sh before                          # 2-4 hours on 2 cores
bash run_v6_analyses.sh after hamilton_results.tar.gz   # once the cluster results are back
bash run_v6_analyses.sh collect                         # RESULTS_SUMMARY.md and errbudget_v6_results.tar.gz
```

Every step logs to `../logs/v6/<step>.txt` and records its exit status in
`../logs/v6/STATUS.txt`. A failed step is reported and the script goes on, except the
environment check and the tests. `ERRBUDGET_JOBS` (default 2) sets how many pooled
configurations run at once.

**`before`** needs nothing from the cluster:

| Step | Script | Produces | What it answers; what to check |
|---|---|---|---|
| A01 | `pytest ../tests -q` | — | all tests pass (52) |
| A02 | `../tests/test_spec_v121.py` | — | 14 lines `ok`: the program's behaviour as Section 4 describes it |
| A03 | `exp75_abinit_reference.py` | `data/abinit_reference_v6.json` | the uncertainty of the 48-Ha ABINIT reference with the Quantum ESPRESSO recipe; the log must say `REFERENCE REPRODUCED` (otherwise the revision-2 value is kept). It decides which ABINIT targets are graded, hence it runs before A05. **Needs ABINIT.** |
| A04 | `exp74_pppm_mechanism.py` | `data/pppm_mechanism_v6.json` | why PPPM order 7 is cheapest (coarser FFT grid): 24 rows, 3 sizes × 2 settings × orders 4–7, and the "cheapest order" lines. **Needs the LAMMPS Python module.** |
| A05 | `pooled.py <config>` × 20 | `data/pooled_<config>.json` | Table 3 and Table 4: v1.2.1, v1.2, v1.1, v1.1+gated, the eight ablations from v1.2.1 and the same eight from v1.2; every log ends with `TOTAL`. The 2D-1 rows of v1.2.1 show 780 starts and 608 accepted at all 7 targets, `missed` 0 everywhere. |
| A06 | `summ.py v1.2.1 v1.2 v1.1` | — | pooled totals, low-order acceptances, step-cap hits |
| A07 | `exp72_baselines_pooled.py` | `data/baselines_v6.json` | Table 6: sequential and coordinate searches against the program, with a "sweep order" line per case; LAMMPS, GROMACS double and mixed, ABINIT and QE included |
| A08 | the same, `LMP-20 LMP-40 LMP-64 GMX-d` | `data/baselines_lmp_gmx.json` | Supplement S5 |
| A09 | `exp76_accounting.py` | `data/accounting_v6.json` | Table 2 against Table 3: candidate lists m1×…×mp, settings without a value, not converged, converged, starts. 2D-1: 840 / 60 / 172 / 608 / 780. |
| A10 | `exp77_reltol_2d1.py` | `data/reltol_2d1.json` | is the GAMG tolerance inactive when 10× below residualControl? Regime line `100 t/rc < 0.1 ... 78 settings, 0 not converged; largest change ... cd 0.0e+00` |
| A11 | `exp78_gromacs_floor_v121.py` | `data/gromacs_floor_v121.json` | GROMACS mixed precision: 3e-5 and 1e-5 reachable only by chance (table below) |
| A12 | `exp58_floor_sensitivity.py v1.2.1` | — | floor sensitivity; its second half times probe runs and **needs GROMACS** |
| A13 | `exp79_unconverged.py` | `data/unconverged_v6.json` | accepting unconverged runs: 2D-1 row 5460 starts, 1204 whose own run did not converge, 4256 accepted by v1.2.1, 5460 with unconverged runs accepted; 1,514 extra acceptances in total |
| A14 | `exp73_decline_reasons.py` | — | Section 6.1; a `TOTAL` line |
| A15 | `exp62_heldout.py` | `data/heldout_replay.json` | Table 5: H1 and H2 with the frozen program; `frozen program verified (sha256 34dd3f44...)`; every row carries `distinct` |
| A16 | `exp62_heldout.py --current` | `data/heldout_replay_current.json` | the same with v1.2.1 |
| A17 | `exp71_dfg2_detail.py` | — | Section 7.2, every 2D-2 target and configuration, including v1.2.1 with cancellation credited (drag at 3e-3: `accepted 172 missed 116`) and the budget at the one setting it returns outside the target (r_k/eps 1.456 and -1.199 for the mesh and the time step, true error +1.031) |
| A18–A21 | `exp65_misdeclare_v12.py`, `exp66_estimate_v12.py`, `exp67_margin_v12.py`, `exp64_transfer_floor.py` | — | Supplement S3, Appendix A, Appendix B, Supplement S2 (rerun because the admitted ABINIT cells may change after A03) |
| A21b | `exp80_v6_numbers.py` | — | numbers that no table holds: the error budgets of ABINIT and Quantum ESPRESSO (v1.2.1, its fits alone, and v1.1), the Phase 1 gains of Fig. 3, the L5 joint check with cancellation credited, the L5 streamfunction minimum, and (part 5) the two L5 settings that miss with a constant safety factor: estimate/eps -0.997 with it and -1.112 with the gated factor, true error -1.080 |
| A21c | `exp81_bound_check.py` | — | Section 4.8 (about 3 minutes): Lemma 1 against the code, `largest relative difference 2.3e-14` at three alignments; Corollary 1 `violated at 0 starts` at all five targets; same-sign higher-order terms `missed 400` at 3e-3, worst 1.083; opposite-sign terms `\|E\| > \|Ehat\| at 0`; the bilinear interaction `at 0 with the joint check and at 16 without it` for b = -60 and -40; at a parameter's finest value `missed 240`, worst 1.281 |
| A22 | `exp63_figures.py` | `figures/fig_search.pdf`, `figures/fig_search_steps.txt` | Fig. 3: 14 steps, Phase 1 ends at step 8, final setting (32, 1e-4, 1e-7, 10), error 0.225 ε, 146 runs; Fig. 4, Fig. 5 |
| A23 | `exp70_tables.py` | `figures/tab_*.tex`, including `tab_ablation_v12.tex` and, for Section S15 of the supplement, `tab_versions_summary.tex`, `tab_versions.tex` and `tab_estimate_v11.tex` | every table of the article and of the supplement; the tables of the article carry no columns of earlier versions |

The deterministic replay of the GROMACS tensor (A11) must give exactly:

| Target | Settings within the target by chance (of 225) | v1.2.1 | v1.1-noprobe |
|---|---|---|---|
| 1e-4 | 206 | 69 accepted (7 distinct), 0 missed | — |
| 3e-5 | 160 | 0 accepted | 152 accepted, 0 missed |
| 1e-5 | 51 | 0 accepted | 9 accepted, 6 missed (worst 1.58 ε) |
| 3e-6 | 13 | 0 accepted | 6 accepted, 6 missed (worst 3.20 ε) |
| 1e-6 | 4 | 0 accepted | — |

**`after <archive>`** needs `hamilton_results.tar.gz` from `../hamilton/`:

| Step | Script | What it does |
|---|---|---|
| B01 | `exp60_hamilton.py` | unpacks into `../hamilton/results`; prints J1–J5; replays H1–H3 with the frozen program and with v1.2.1 (`heldout_replay*.json` then contain case `H3`); regenerates the tables |
| B02 | `pooled.py --only 2D-2h` × 20 | merges the refined 2D-2 study into every `pooled_<config>.json` (skipped, and reported, if `j2_tensor.pkl` lacks the case specification) |
| B03 | `exp72_baselines_pooled.py 2D-2h` | 2D-2h baselines merged into `baselines_v6.json` |
| B04–B06 | `exp76`, `exp79`, `exp73` | accounting, unconverged counts and decline reasons including 2D-2h |
| B07 | `exp70_tables.py` | final tables; Table 3 uses the refined 2D-2 study (`ERRBUDGET_TABLE3_2D2=refined`, the default) |

**`collect`** runs `exp82_collect_v6.py`: `RESULTS_SUMMARY.md` in the repository root
(headline numbers for every configuration with both 2D-2 studies, the held-out tables,
accounting, baselines, every revision-4 analysis, the tails of `j1.txt`–`j5.txt`) and
`errbudget_v6_results.tar.gz` (the summary, every result file, tables, figures, logs and
`hamilton/results` without the raw archive).

### Earlier revisions: `reproduce.sh --cached`

```bash
./reproduce.sh --cached   # re-analyse the stored solver output of revisions 2 and 3 (about 2 hours)
```

Each step writes its output to `figures/<step>.txt`; the tables further down map every
result of revisions 3 and 2 to its script. Those results are now in the supplement or in
`../CHANGELOG.md`.

## Tier B: rebuild the solver output from scratch

```bash
./reproduce.sh            # deletes data/*.pkl and data/*.json, runs every solver again (about 10 hours on 2 cores)
```

Then the DFG 2D-2 grid, which is a batch job of 243 OpenFOAM runs (about 40 CPU-hours;
SLURM or one machine): see `dfg2d2_job/README.md`. Finally `run_v6_analyses.sh before`
as in tier A.

Nothing is read from saved output unless `--cached` is given. The stored solver output is
in `data/` (pickled tensors, one per case, plus JSON result files), and `cases.py` is the
registry of every case: its tensor, declarations, reference, admission rule and targets.

- **Machine-dependent numbers:** the costs of every production code (wall-clock time for
  LAMMPS and GROMACS, CPU time for OpenFOAM, ABINIT and Quantum ESPRESSO), the penalties
  derived from them, the CPU-time replay of `exp52`, the 8-parameter study `exp54`
  (which uses work units but whose recommendations are graded by them), and running times.
- **Machine-independent numbers:** everything else, which reproduces exactly.
- **One exception, the L2 cost model.** Its work-unit cost counts GMRES iterations. At loose
  tolerances this count can differ by one between SciPy builds, which changes L2 costs by up
  to 0.25%. In our checks this changed no number in the article.

## Tier C: the cluster studies

Everything that runs on the cluster is in `../hamilton/`: the five jobs, their drivers,
the submit, status and collect scripts, and the runbook (`../hamilton/README.md`) that
takes the results back to `run_v6_analyses.sh after`.

## What is in this folder

```
run_v6_analyses.sh   revision 4: check | before | after <archive> | collect
reproduce.sh         revisions 2 and 3, from scratch or --cached
cases.py             the case registry (every tensor, declaration, reference, admission rule)
pooled.py, summ.py   the pooled replay over every case and configuration, and its totals
replay.py, grids.py, baselines.py   shared replay, grid and baseline code
exp*.py              one script per analysis; the tables below say which result each one produces
*_silicon.py, gromacs_water.py, lammps_*.py, openfoam_*.py, heldout_*.py
                     drivers of the solvers (ABINIT, QE, GROMACS, LAMMPS, OpenFOAM, held-out problems)
snfix.c              the snprintf shim for Quantum ESPRESSO 6.7 (see ../THIRD_PARTY_RESULTS.md)
data/                stored solver output and result files
figures/             tables (tab_*.tex), figures (fig_*.pdf) and the text output of every step
frozen_src/          copy of version 1.2.0, the frozen program of the held-out tests
dfg2d2_job/          batch build of the revision-2 2D-2 grid (SLURM or one machine)
```

## Environment used for the article

| Component | Version |
|---|---|
| OS | Ubuntu 24.04, 2 vCPU, 7 GB RAM |
| Python | 3.11 |
| numpy | 1.26.4 |
| scipy | 1.13.1 (1.12 or later is required) |
| Other Python packages | matplotlib, sympy, pytest |
| LAMMPS | 20240207, with its Python module (`apt install lammps python3-lammps`) |
| GROMACS | 2023.3, `gmx` (mixed precision) and `gmx_d` (double precision) (`apt install gromacs`) |
| ABINIT | 9.10.4, pseudopotentials in `/usr/share/abinit/psp` (`apt install abinit`) |
| OpenFOAM | v1912 (`apt install openfoam`); `simpleFoam`, `pimpleFoam`, `blockMesh`, `checkMesh` |
| Quantum ESPRESSO | 6.7 (`apt install quantum-espresso`), run with the preloaded shim `snfix.c` (needs `gcc`) |


## Where each result of the current article comes from (revision 5, version 1.2.1)

Section, table and figure numbers are those of revision 5 of the article and its
supplement. Every step is in `run_v6_analyses.sh before` (the step is given in brackets).
The tables of the article report the released program only; the comparisons with earlier
versions are in Section S15 of the supplement and in `../CHANGELOG.md`.

| Result | Script [step] | Output |
|---|---|---|
| Table 4 (pooled), totals in Section 6.1 and in the abstract | `pooled.py v1.2.1` [A05], then `exp70_tables.py` [A23] | `data/pooled_v1.2.1.json`, `figures/tab_pooled.tex` |
| Table 5 (ablation, including the frozen program's joint-check threshold); Sections 6.1 and 8 | `pooled.py v1.2.1-* v1.2` [A05], then `exp70_tables.py` | `figures/tab_ablation.tex` |
| Section 4.8: Lemma 1, Proposition 1, Corollary 1 and the synthetic tests of the Remark | `exp81_bound_check.py` [A21c] | stdout |
| Section 4.8, Remark: the L5 settings that miss with a constant safety factor | `exp80_v6_numbers.py`, part 5 [A21b] | stdout |
| Section 4: the central rules of the specification | `../tests/test_spec_v121.py` [A02] | — |
| Section 4.6 and Fig. 3: Phase 1 gains of the cost-only moves | `exp80_v6_numbers.py`, part 2 [A21b] | stdout |
| Section 6.1: low-order acceptances, step cap | `summ.py v1.2.1 v1.2 v1.1` [A06] | stdout |
| Section 6.1: decline reasons | `exp73_decline_reasons.py` [A14] | stdout |
| Section 6.1: accepting unconverged runs | `exp79_unconverged.py` [A13] | `data/unconverged_v6.json` |
| Section 6.1: the L5 joint check with cancellation credited; the L5 streamfunction minimum | `exp80_v6_numbers.py`, parts 3 and 4 [A21b] | stdout |
| Table 3 and Section 5.2: accounting of every grid | `exp76_accounting.py` [A09] | `data/accounting_v6.json` |
| Table 6 and Section 6.2: held-out H1, H2 (H3 after the cluster run) | `exp62_heldout.py` [A15, frozen], `exp62_heldout.py --current` [A16] | `data/heldout_replay*.json`, `figures/tab_heldout.tex` |
| Table 7, Table S6, Sections 6.3 and 6.4: baselines and calibration cost | `exp72_baselines_pooled.py` [A07], then `exp70_tables.py` | `data/baselines_v6.json`, `figures/tab_baselines_compact.tex`, `figures/tab_baselines_full.tex` |
| Table S4: LAMMPS and GROMACS double precision against the baselines | `exp72_baselines_pooled.py LMP-20 LMP-40 LMP-64 GMX-d` [A08] | `figures/tab_lmpbase.tex` |
| Fig. 3, Fig. 4, Fig. 5 | `exp63_figures.py` [A22] | `figures/fig_search.pdf`, `figures/fig_l5_stall.pdf`, `figures/fig_dfg_mesh.pdf` |
| Section 7.1: when the GAMG tolerance acts | `exp77_reltol_2d1.py` [A10] | `data/reltol_2d1.json` |
| Section 7.2 and the Remark of Section 4.8: 2D-2 below the admission threshold, with and without credit for cancellation | `exp71_dfg2_detail.py` [A17] | stdout |
| Section 7.3: PPPM grid and cost split (needs the LAMMPS Python module) | `exp74_pppm_mechanism.py` [A04] | `data/pppm_mechanism_v6.json` |
| Section 7.3: GROMACS round-off floor | `exp78_gromacs_floor_v121.py` [A11] | `data/gromacs_floor_v121.json` |
| Section 7.3, Section S4, Table S3: floor sensitivity and probe cost | `exp58_floor_sensitivity.py v1.2.1` [A12] | `figures/tab_floor.tex` |
| Section 7.4: ABINIT reference uncertainty (needs ABINIT) | `exp75_abinit_reference.py` [A03] | `data/abinit_reference_v6.json` |
| Section 7.4: error budgets of ABINIT and Quantum ESPRESSO (program and fits alone); ABINIT errors along the cutoff | `exp80_v6_numbers.py`, part 1 [A21b] | stdout |
| Section 7.4: ABINIT misses without the local-order check | `pooled.py v1.2.1-nolocal` [A05] | `data/pooled_v1.2.1-nolocal.json` (`miss_list`) |
| Appendix A, Table A.1 | `exp66_estimate_v12.py` [A19], then `exp70_tables.py` | `figures/tab_estimate.tex` |
| Appendix B | `exp67_margin_v12.py` [A20] | stdout |
| Section 3.2, Section S2, Table S1 | `exp64_transfer_floor.py` [A21] | `figures/tab_transfer.tex` |
| Section S3, Table S2: misdeclared classes | `exp65_misdeclare_v12.py` [A18] | `figures/tab_misdeclare.tex` |
| Section S15, Tables S9 and S10: earlier versions over every cell | `pooled.py v1.2 v1.1+gated v1.1` [A05], then `exp70_tables.py` | `figures/tab_versions_summary.tex`, `figures/tab_versions.tex` |
| Section S15, Table S11: accuracy of the estimate of version 1.1 | `exp66_estimate_v12.py` [A19], then `exp70_tables.py` | `figures/tab_estimate_v11.tex` |
| Section S15: held-out problems, 2D-2, GROMACS and ABINIT with earlier versions | `exp62_heldout.py` [A15], `exp71_dfg2_detail.py` [A17], `exp78_gromacs_floor_v121.py` [A11], `exp80_v6_numbers.py` [A21b] | stdout |
| Cluster studies J1–J5 | `../hamilton/` (see its README), then `run_v6_analyses.sh after` | `../hamilton/results/` |

The studies of Sections S8–S14 of the supplement were not repeated with version 1.2.1;
their scripts are listed under revision 2 below.

## Where each result of revision 3 came from (version 1.2.1)

| Result | Script | Output |
|---|---|---|
| Table 3 (pooled), Table 4 (ablation), totals in the text | `pooled.py v1.2.1 v1.2 v1.1 v1.2-nojoint …` then `exp70_tables.py` | `data/pooled_*.json`, `figures/tab_pooled.tex`, `figures/tab_ablation.tex` |
| Pooled totals, low-order acceptances, step-cap hits | `summ.py v1.2.1 v1.2 v1.1` | stdout |
| Case registry (every tensor, declaration, reference, admission rule) | `cases.py` | — |
| Section 4.1: L2–L5 convergence flags | `exp57_flags.py` | `data/flags_L*.pkl` |
| Table 5 and Section 6.2: held-out H1, H2 (H3 after the cluster run) | `exp61_heldout_build.py H1/H2`, `exp62_heldout.py` (frozen, from the tag or `frozen_src/`), `exp62_heldout.py --current` | `figures/rev3_heldout*.txt`, `figures/tab_heldout.tex` |
| Section 6.2: post-hoc check of the amendment | `posthoc_v121.py` | — |
| Table 6 and Section 6.3: baselines | `exp72_baselines_pooled.py` | `figures/rev3_baselines_pooled.txt`, `figures/tab_baselines.tex` |
| Section 6.1: decline reasons | `exp73_decline_reasons.py` | `figures/rev3_decline_reasons.txt` |
| Fig. 3: one search | `exp63_figures.py` | `figures/fig_search.pdf` |
| Fig. 4: L5 stall | `exp63_figures.py` | `figures/fig_l5_stall.pdf` |
| Fig. 5, Section 7.1: 2D-1 rebuilt with exact centroids; n_arc 96 and 128 | `exp46_dfg_build.py`, `exp69_dfg_fine.py`, `exp63_figures.py` | `data/dfg_tensor_v12.pkl`, `data/dfg_fine.json`, `figures/fig_dfg_mesh.pdf` |
| Section 7.2: 2D-2, every target and configuration | `exp71_dfg2_detail.py` | `figures/rev3_dfg2_detail.txt` |
| Section 7.3, Supplement S4: GROMACS floor sensitivity, probe cost | `exp58_floor_sensitivity.py v1.2` | `figures/rev3_floor_sensitivity.txt`, `figures/tab_floor.tex` |
| Section 7.4: ABINIT SCF convergence check | `exp68_abinit_scf.py` | `figures/rev3_abinit_scf.txt` |
| Appendix A: accuracy of the estimate | `exp66_estimate_v12.py` | `figures/rev3_estimate.txt`, `figures/tab_estimate.tex` |
| Appendix B: base safety factor | `exp67_margin_v12.py` | `figures/rev3_margin.txt` |
| Section 3.2, Supplement S2: transfer with a resolution floor | `exp64_transfer_floor.py` | `figures/rev3_transfer_floor.txt`, `figures/tab_transfer.tex` |
| Supplement S3: misdeclared classes | `exp65_misdeclare_v12.py` | `figures/rev3_misdeclare.txt`, `figures/tab_misdeclare.tex` |
| Supplement S5: LAMMPS and GROMACS baselines | `exp72_baselines_pooled.py` (cells `LMP-*`, `GMX-d`) | `data/baselines_lmp_gmx.json`, `figures/tab_lmpbase.tex` |
| Cluster studies J1–J5 | `../hamilton/` (see its README), then `exp60_hamilton.py hamilton_results.tar.gz` | `../hamilton/results/` |
| Behavioural claims (38 tests) | `pytest ../tests` | — |

The sections below map the results of earlier versions of the article. Those studies are
now in the supplement or in `../CHANGELOG.md`.

## Where each result of revision 2 came from

| Article | Script | Output (`figures/`) |
|---|---|---|
| Section 3.1: class inference (3 or 4 of 8) | `exp17_classify.py` | `sec31_classes.txt` |
| Table 2: L4 error budget | `exp20_reactdiff.py` | `build_L4.txt` |
| Section 5.1: masking on L2 and L3; Section 5.3: two observables | `exp40_section5.py` | `sec5_masking.txt` |
| Section 5.1: round-off floor of the observables | `exp25_referee.py` | `sec5_floor.txt` |
| Table 3 and Section 5.2: transfer of thresholds | `exp32_conditions.py` | `tab03_transfer.txt` |
| Section 5.2: conditioning mechanism | `exp31_mechanism.py` | `sec52_mechanism.txt` |
| Table 4: design ablations | `exp41_ablations.py` | `tab04_ablations.txt` |
| Table 5 and Section 7.1: all 11,048 replays | `exp26_outcomes.py` | `tab05_outcomes.txt` |
| Table 6: separability and estimate quality | `exp28_separability.py` | `tab06_rho.txt` |
| Table 7: ladder spacing | `exp34_spacing.py`, `exp35_whyfine.py`, `exp36_gated.py` | `tab07_*.txt` |
| Table 8: baselines | `exp38_baselines.py` | `tab08_baselines.txt` |
| Section 7.6: margin scan | `exp29_margin.py` | `sec76_margin.txt` |
| Table 9: transfer across ν | `exp30_transfer.py` | `tab09_transfer.txt` |
| Section 7.7: calibrating at the hardest member | `exp37_calhardest.py` | `sec77_hardest.txt` |
| Table 10: calibration cost and break-even; calibration runs vs distance | `exp42_calibration.py` | `tab10_calibration.txt` |
| LAMMPS table, Fig. 3: LAMMPS at 8,000 / 64,000 / 262,144 ions | `exp22_lammps.py` (`ERRBUDGET_LMP_N=20,40,64`), `make_lammps_figure.py` | `sec81_lammps_n*.txt` |
| Table 12: LAMMPS against baselines | `exp39_lmp_baselines.py` | `sec81_baselines.txt` |
| Section 8.1: the certified misses and their cause | `exp43_lammps_misses.py` | `sec81_misses.txt` |
| Table 13: GROMACS | `exp24_gromacs.py` | `sec82_gromacs_d.txt` |
| Section 8.2: precision floor of `gmx` and `gmx_d` | `exp44_gromacs_floor.py` | `sec82_floor.txt` |
| Table 14, Section 8.3: ABINIT budget and outcomes | `exp27_abinit.py` | `sec83_abinit.txt` |
| Section 8.3: uncertainty of the ABINIT reference | `exp45_abinit_reference.py` | `sec83_reference.txt` |
| Section 8.4: learned component | `exp23_mlhybrid.py` | `sec84_learned.txt` |
| Fig. 2: test problems | `make_solution_figure.py` | `fig0_problems.pdf` |
| Section 6, L5 and Table 7.x: driven cavity, non-integral observables | `exp48_cavity.py` | `sec5x_cavity.txt` |
| Table 3 on L2–L5 and OpenFOAM | `exp47_transfer.py` | `tab03_all.txt` |
| Section 3.1: misdeclared classes (L3, OpenFOAM) | `exp49_misdeclare.py` | `tab_misdeclare.txt` |
| Appendix B: margin calibrated on each problem in turn | `exp50_margin_heldout.py` | `app_margin_held.txt` |
| Section 9: adaptive BDF integrator | `exp51_adaptive.py` | `sec9x_adaptive.txt` |
| Section 10: Table 5 with CPU time as cost | `exp52_wallclock.py` | `tab05_cputime.txt` |
| Section 7: eight parameters, 200 starts | `exp54_eight.py` | `sec7x_eight.txt` |
| Section 8.1, Table: OpenFOAM DFG 2D-1 | `exp46_dfg_build.py`, `exp46_dfg.py` | `build_dfg.txt`, `sec8x_dfg.txt` |
| Section 8.1: OpenFOAM DFG 2D-2 (Re = 100) | `exp56_dfg2d2_build.py` (run through `dfg2d2_job/`, see its README), `exp56_dfg2d2.py` | `sec8x_dfg2d2.txt` |
| Section 8.4: ABINIT and Quantum ESPRESSO | `exp27_abinit.py`, `exp45_abinit_reference.py`, `exp53_qe.py` | `sec83_*.txt`, `sec8x_qe.txt` |
| Section 8.6: parallel execution | `exp55_parallel.py` (`ERRBUDGET_RANKS`, `ERRBUDGET_THREADS`) | `sec8x_parallel.txt` |
| Behavioural claims (31 tests) | `pytest ../tests` | — |

Exhaustive grids are built before the analyses that use them:
`exp00_build_L2.py`, `exp07_2d.py` (L3), `exp20_reactdiff.py` (L4),
`exp33_fineladder.py`, `exp34_spacing.py` and `exp37_calhardest.py`. These grids are used
only to grade the method. The method itself never sees more than the configurations it
chooses to evaluate.


## Version 1.1 (revision 2)

The method changed in three places (see `src/errbudget/_core.py`):
`unfit="order1"` (a parameter that cannot be fitted is charged its one-step change
extrapolated at first order), the `noise_floor` check, and `models=[..."exponential"...]`.
Every table was regenerated with version 1.1. The stored solver output in `data/` and
the text output of every script in `figures/` are included so that results can be
inspected without rerunning; `./reproduce.sh` regenerates them from scratch.
