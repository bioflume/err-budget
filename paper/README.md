# Reproducing the article

This folder regenerates every number, table and figure of the article (v12) and its
supplement. Each result below names the script that produces it, the step of
`replay_all.sh` that runs the script, and the stored data that the script reads.

There are three levels, from cheap to expensive:

| Level | What it does | Needs | Time |
|---|---|---|---|
| **Replay** (`replay_all.sh`) | searches the stored runs in `data/` and `../hamilton/results/` with the released program and regenerates every result that needs no solver | Python 3.11, `requirements-article.txt` | about 2–4 hours on two cores, plus about 6.5 CPU-hours for the eight-parameter study (`--no-eight` leaves it out) |
| **Rebuild** (`reproduce.sh`) | runs the solvers again, rebuilds `data/`, runs the analyses that need a solver, then replays | the replay, plus LAMMPS, GROMACS, ABINIT, OpenFOAM v1912 and Quantum ESPRESSO | about 16 CPU-hours, plus about 40 for the coarser 2D-2 study (`dfg2d2_job/`) |
| **Cluster studies** | the Hamilton jobs J1–J5 (`../hamilton/`) and the studies of Sections S-6 and S-7 (`../studies/`) | a SLURM cluster with OpenFOAM | days of queue time; see each runbook |

## Setup

```bash
git clone <this repository> && cd <it> && git checkout paper
python3.11 -m venv .venv && source .venv/bin/activate
pip install -e . -r paper/requirements-article.txt
bash paper/replay_all.sh                  # or: ERRBUDGET_JOBS=4 bash paper/replay_all.sh
git status paper/data paper/figures       # what the replay changed
```

The replay writes its results where the stored ones are, so `git diff` compares the two. The
printed output of each step goes to `logs/replay/<step>.txt`. Results that the replay prints
but does not store are compared with the stored printed output in `figures/<name>.txt`.
`REPRODUCTION_CHECK.md` records what was reproduced on which machine.

The searches are deterministic: on Linux x86-64 and aarch64 every count of the stored
replays is reproduced exactly. Scripts that rerun a solver (the L2 viscosity study, the
exact solutions of L3–L5) can differ from the stored values in the last digits on another
platform; no rounded number of the article changes.

## Article

| Result | Script (step of `replay_all.sh`) | Reads | Writes |
|---|---|---|---|
| Fig. 1, Section 1 and 7.1 (lift error 1.3×10⁻² to 0.79, the program's setting at 0.23 ε) | `extract_2d1_l5.py` (fig_2d1_l5_data), `make_2d1_l5_figures.py` (fig_2d1_l5) | `data/dfg_tensor_v12.pkl`, `data/dfg_fine.json` | `figures/data_2d1_l5.json`, `figures/fig_teaser.pdf`; printed values as in `figures/values_2d1_l5.txt` |
| Table 1, Table 2, Fig. 2 | none (literature, notation, schematic) | | |
| Fig. 3; 21 searches miss without the observed-order check | `make_2d1_l5_figures.py`; `pooled.py v1.2.1-nolocal` | `data/cavity_tensor.pkl` | `figures/fig_l5_stall.pdf`; `data/pooled_v1.2.1-nolocal.json` (L5, y of the minimum, 3×10⁻⁴) |
| Section 4.3, the base safety factor | `exp67_margin_v12.py` (margin) | stored runs of L2–L5, 2D-1 | printed |
| Table 3: settings and converged runs | `exp76_accounting.py` (accounting) | every case of `cases.py` | `data/accounting_v6.json` |
| Table 3: reference uncertainties | `exp75_abinit_reference.py` (ABINIT, rebuild only), `../hamilton/J2_dfg2d2/analyse_j2.py` (2D-2), `exp44_gromacs_floor.py` (GROMACS, rebuild only), `data/qe_tensor.pkl` (Quantum ESPRESSO) | | `data/abinit_reference_v6.json`, `../hamilton/results/j2.txt`, `figures/sec82_floor.txt` |
| Table 4 (`tab_pooled_compact.tex`) and Section 6.1 | `pooled.py v1.2.1` (pooled), `exp72_baselines_pooled.py` (baselines), `exp70_tables.py` (tables) | `data/`, `../hamilton/results/j2_tensor.pkl` | `data/pooled_v1.2.1.json`, `data/baselines_v6.json`, `figures/tab_pooled_compact.tex` |
| Section 6.1: 436 and 105 misses of the ablation; up to 36% underestimates on the 2D-1 lift | `pooled.py v1.2.1-constant`, `v1.2.1-nolocal`; `exp66_estimate_v12.py` (estimate) | | Table S-3.2, Table S-3.3 |
| Table 5 (`tab_preregistered.tex`) and Section 6.2 | `exp62_heldout.py` (heldout_frozen), `exp62_heldout.py --current` (heldout_released), `exp70_tables.py` | `data/heldout_H1.pkl`, `data/heldout_H2.pkl`, `../hamilton/results/j5_tensor.pkl`, the frozen program (tag `v1.2.0-frozen` or `frozen_src/`) | `data/heldout_replay.json`, `data/heldout_replay_current.json`, `figures/tab_preregistered.tex` |
| Section 6.2, the airfoil study | `../studies/airfoil/` (Section S-7) | | |
| Fig. 4 and Section 6.3 (S and SJ misses, the safe thresholds, 2.7 to 40 times cheaper, search cost 14 to 32 and six to seven times SJ's) | `exp72_baselines_pooled.py`, `exp70_tables.py`, `make_baselines_figure.py` (fig_baselines) | `data/` | `figures/tab_baselines_full.tex`, `figures/fig_baselines.pdf`; printed values as in `figures/values_baselines.txt` |
| Section 6.3, PRE (7–30% coverage, 8.6% of searches, 1.84 ε) | `../studies/pre_comparison/` (Section S-4.6) | | `results/summary.json`, `results/tables.md`, `REPORT.md` |
| Section 6.4: search cost K (40–460, 1,550 on L5, median 112; S 4–8) | Table S-4.2 | | |
| Section 6.4: reuse across Reynolds numbers | `../studies/reuse2d1/` (Section S-6) | | |
| Section 7.1: the GAMG tolerance acts through 100 t/r_c | `exp77_reltol_2d1.py` (reltol_2d1) | `data/dfg_tensor_v12.pkl` | `data/reltol_2d1.json` |
| Section 7.2: 180 runs, 40 diverged | `../hamilton/J2_dfg2d2/analyse_j2.py` | J2 raw runs | `../hamilton/results/j2.txt` |
| Section 7.2: 136 misses at two ungraded targets, up to 2.39 times | `exp83_dfg2_traces.py` (dfg2_traces) | `../hamilton/results/j2_tensor.pkl` | `data/dfg2_traces.json`, `figures/dfg2_traces.txt` |
| Section 7.3: LAMMPS within 1.0–2.6 times the optimal cost; 400 of 660 declined | Table S-4.3; Table S-3.1 | | |
| Section 7.3: GROMACS round-off 3.0×10⁻⁵ and 2.3×10⁻⁹; 450 declined | `exp44_gromacs_floor.py` (rebuild only); Table S-3.1 | | `figures/sec82_floor.txt` |
| Section 7.4: 82% and 91% of the errors; force error 27% | `exp80_v6_numbers.py` (numbers), part 1 | `data/qe_tensor.pkl`, `data/abinit_tensor_v11.pkl` | printed |
| Section 8: 19% declines; the stricter rule | Table S-3.1; `exp84_dfg2_strict.py` (dfg2_strict) | | `data/dfg2_strict.json`, `figures/dfg2_strict.txt` |
| Section 9: the interface | the package on the `main` branch (`errbudget init/check/run` is its command-line tool); `tests/test_spec_v121.py` (14 tests) | | |

## Supplement

| Result | Script (step) | Reads | Writes |
|---|---|---|---|
| S-1, Tables S-1.1–S-1.3: the specification | `../src/errbudget/_core.py`, `api.py` (the engine as it produced the results) | | |
| S-1.9: tests of the conditions of Proposition 1 | `exp81_bound_check.py` (bound_check) | synthetic problems | printed |
| S-1.10: scan of the base safety factor | `exp67_margin_v12.py` (margin) | | printed |
| S-1.11: the 14 tests | `../tests/test_spec_v121.py` (tests) | | |
| S-2.1, S-2.2, Tables S-2.1, S-2.2: candidate values and work units | `grids.py:10–38`, `../src/errbudget/problems/`, `qe_silicon.py:38–54`, `abinit_silicon.py:33–53`, `cases.py`, `exp46_dfg_build.py` | | |
| Fig. S-2.1 | `make_solution_figure.py` (fig_problems) | | `figures/fig0_problems.pdf` |
| S-2.3: the coarser 2D-2 study (240 settings, its reference); 116 misses with cancellation | `exp56_dfg2d2.py` (dfg2d2_coarse); `exp71_dfg2_detail.py` (dfg2_detail) | `data/dfg2d2_tensor.pkl` | printed, as in `figures/sec8x_dfg2d2.txt`, `figures/rev3_dfg2_detail.txt` |
| S-2.4: the pre-registered problems | `../PREREGISTRATION.md`; `heldout_fem.py`, `heldout_cavity.py`, `../hamilton/J5_abc3d/`; H3 errors at the finest setting: `../hamilton/results/j5.txt` | | |
| S-2.5: references; the remark on limits beyond the 2D-2 reference | `../hamilton/results/j2.txt`; `exp83_dfg2_traces.py`, part 1 | | `data/dfg2_traces.json` |
| S-2.6: accounting | `exp76_accounting.py` (accounting) | | `data/accounting_v6.json` |
| S-3.2, S-5.1, S-5.2: the numbers that no script printed before (LAMMPS errors along the cutoffs, the ABINIT estimate at the finest setting, momentum solves stopped by their tolerance, the Courant number and the diverged 2D-2 runs) | `exp85_unscripted.py` (unscripted) | `../hamilton/results/dfg_tensor_current.pkl`, `j2_tensor.pkl`, `data/lammps_tensor_n20.pkl`, the ABINIT tensor | printed, as in `figures/unscripted.txt` |
| Table S-3.1 (`tab_pooled.tex`) | `pooled.py v1.2.1`, `exp70_tables.py` | | `figures/tab_pooled.tex` |
| S-3.2: declines by reason; L5 orders | `exp73_decline_reasons.py` (decline_reasons); `exp80_v6_numbers.py`, part 4 | | printed |
| Table S-3.2 (`tab_sup_ablation.tex`) and S-3.3 | `pooled.py v1.2.1-<variant>` (8 variants), `exp70_tables.py`; 0.997, 1.08 and 1.11 ε: `exp80_v6_numbers.py`, part 5; 116 of 172: `exp71_dfg2_detail.py`; accepting unconverged runs: `exp79_unconverged.py` (unconverged); 2D-1 with the function object: `../hamilton/results/j1.txt` | | `figures/tab_sup_ablation.tex`, `data/unconverged_v6.json` |
| S-3.4: the additive model (57,000 designs); the coarser 2D-2 study | `exp28_separability.py` (separability); `exp56_dfg2d2.py` | | printed, as in `figures/tab06_rho.txt` |
| Table S-3.3 (`tab_estimate.tex`) | `exp66_estimate_v12.py` (estimate), `exp70_tables.py` | | `data/estimate_v12.json` |
| S-3.5, Table S-3.4: closely spaced values | `exp34_spacing.py` (spacing_released), `exp35_whyfine.py` (spacing_estimate_released), `exp34_spacing_margin.py` (spacing_margin), with `l2flags.py` | `data/burgers_sp_*.pkl` | `data/spacing_released.json`, `data/spacing_estimate_released.json`, `data/spacing_margin_check.json` |
| S-3.5: eight parameters | `exp54_eight.py` (eight_released) | `data/eight_cache.pkl` | `data/eight_released.json` |
| S-3.6, Table S-3.5: wrongly declared classes | `exp17_classify.py` (classify); `exp65_misdeclare_v12.py` (misdeclare) | `data/classify.pkl` | `figures/tab_misdeclare.tex` |
| S-3.7: traces of the pre-registered searches | `exp62_heldout.py`; `../hamilton/J5_abc3d/analyse_j5.py` | | printed; `../hamilton/results/j5.txt` |
| Table S-4.1 (`tab_transfer.tex`); two measurements on L2 | `exp64_transfer_floor.py` (transfer), `exp70_tables.py`; `exp40_section5.py` (section5) | | `data/transfer_floor.json` |
| S-4.2: the three procedures | `baselines.py` | | |
| Table S-4.2 (`tab_baselines_full.tex`) | `exp72_baselines_pooled.py` (baselines), `exp70_tables.py` | | `data/baselines_v6.json` |
| Table S-4.3 (`tab_lmpbase.tex`) | `exp72_baselines_pooled.py LMP-20 LMP-40 LMP-64 GMX-d` (baselines_md), `exp70_tables.py` | | `data/baselines_lmp_gmx.json` |
| S-4.5, Table S-4.4: reuse on L2 | `exp30_transfer.py` (reuse_released), `exp37_calhardest.py` (reuse_hardest_released) | `data/burgers_tensor.pkl`, `data/burgers_nu002_tensor.pkl` | `data/reuse_released.json`, `data/reuse_hardest_released.json` |
| S-4.6, Tables S-4.5, S-4.6: PRE | `../studies/pre_comparison/` (`replay_all.sh --pre`) | the 2D-1 runs in `data/` | `results/` |
| Fig. S-4.1 | `make_pre_figure.py` (fig_pre) | `../studies/pre_comparison/results/` | `figures/fig_sup_pre.pdf` |
| Figs. S-5.1, S-5.2 | `extract_2d1_l5.py`, `make_2d1_l5_figures.py` | | `figures/fig_sup_dfg.pdf`, `figures/fig_sup_search.pdf` |
| S-5.1: the regimes of the GAMG tolerance; `relTol` 0; why the linear solves stopped | `exp77_reltol_2d1.py`; the v2406 repeat: `../hamilton/J1_dfg2d1/analyse_j1.py` | | `data/reltol_2d1.json`, `../hamilton/results/j1.txt`, `j1_reltol.json`, `j1_stop.json` |
| S-5.1: vertex-mean bias, runs beyond the study | `exp59_dfg_plateau.py`, `exp69_dfg_fine.py` (rebuild only) | | `data/dfg_plateau.json`, `data/dfg_fine.json` |
| S-5.1: the function object of OpenFOAM v2406 | `../hamilton/J1_dfg2d1/analyse_j1.py` | | `../hamilton/results/j1.txt`, `j1_fo.json`, `j1_replay.json` |
| S-5.2: divergence, convergence along each parameter, the conservative estimate | `../hamilton/J2_dfg2d2/analyse_j2.py` | | `../hamilton/results/j2.txt` |
| S-5.2: misses at the ungraded targets, the traced search, the mesh errors | `exp83_dfg2_traces.py` (dfg2_traces) | `../hamilton/results/j2_tensor.pkl` | `data/dfg2_traces.json` |
| S-5.2: the stricter rule | `exp84_dfg2_strict.py` (dfg2_strict) | `../hamilton/results/j2_tensor.pkl`, `data/pooled_v1.2.1.json` | `data/dfg2_strict.json` |
| S-5.3: PPPM order | `exp74_pppm_mechanism.py` (rebuild only; needs LAMMPS) | | `data/pppm_mechanism_v6.json` |
| Fig. S-5.3 | `make_lammps_figure.py` (fig_lammps) | `data/lammps_tensor_n20.pkl` | `figures/fig_lammps.pdf` |
| S-5.4: round-off measurement | `exp44_gromacs_floor.py` (rebuild only); `exp78_gromacs_floor_v121.py` (gromacs_floor); `exp58_floor_sensitivity.py v1.2.1` (floor_sensitivity; its timings need GROMACS) | | `figures/sec82_floor.txt`, `data/gromacs_floor_v121.json`, `figures/tab_floor.tex` |
| S-5.5: plane-wave DFT | `exp75_abinit_reference.py`, `exp45_abinit_reference.py` (earlier 30-Ha reference), `exp68_abinit_scf.py` (rebuild only); `exp80_v6_numbers.py`, part 1 | | `data/abinit_reference_v6.json`, `figures/sec83_reference.txt`, `figures/rev3_abinit_scf.txt` |
| S-5.6: OpenFOAM v2406, Quantum ESPRESSO 7.6, MPI ranks | `../hamilton/J1_dfg2d1/`, `J3_qe/`, `J4_parallel/` analyses | | `../hamilton/results/j1.txt`, `j3.txt`, `j4.txt` |
| S-6, Tables S-6.1, S-6.2: reuse across Reynolds numbers | `../studies/reuse2d1/analysis/reuse.py` (reuse2d1, when the raw pack is unpacked) | the raw pack (`../DATA.md`) | `../studies/reuse2d1/summary/reuse2d1-20261115-1111/` |
| S-7: the airfoil study | `../studies/airfoil/analysis/stage0.py`, `stage1.py` (airfoil_stage1, when the raw pack is unpacked) | the raw packs (`../DATA.md`) | `../studies/airfoil/summary/` |

"Rebuild only": the script runs a solver and is part of `reproduce.sh`, not of the replay;
its stored output is listed.

## Cluster studies

- **Hamilton jobs J1–J5** (`../hamilton/`, runbook `../hamilton/README.md`): 2D-1 with
  OpenFOAM v2406, the refined 2D-2 study, Quantum ESPRESSO 7.6, MPI ranks and the
  pre-registered problem H3. Their summarised results are in `../hamilton/results/`; after
  a new run, `python exp60_hamilton.py hamilton_results.tar.gz` unpacks them, replays H1–H3
  and regenerates the tables.
- **The Reynolds-number reuse study** (`../studies/reuse2d1/`, runbook `RUNBOOK.md`) and
  **the airfoil study** (`../studies/airfoil/`, runbook `RUNBOOK.md`, cost model `COST.md`):
  both pre-registered. Their protocols were committed before any run, and each amendment of
  the airfoil protocol alone, in the project repository: airfoil da6e9e0 (10 November
  2026), amendments A1 5774f09, A2 9b6f670, A3 31f17e0, A4 5b8cdc5, A5 4a1d510, A6 88c1910;
  reuse2d1 8cb358c (14 November 2026, no amendment). That project repository is not
  available, and none of these commits is in this repository, so their dates cannot be
  checked here; `../TRUST.md` says what evidence of the order of events remains. The code that ran on Hamilton is that
  of the batch tarballs listed in `../DATA.md` (reuse2d1: 9b7482e). The raw packs are not in
  git; their analyses are committed under `summary/`.
- **The PRE comparison** (`../studies/pre_comparison/`): protocol committed before any PRE
  result was graded (e1d1bce, 9 November 2026); the commit is not available either). `fetch_pre.sh` clones the PRE authors' code
  from https://github.com/NewcastleRSE/sparse-probabilistic-richardson-extrapolation at
  commit 3918107 (MIT licence, © 2025–2026 Oates, Howey and Karvonen) and installs its
  environment; nothing of it is kept here. `spre_white.py` restates their model in closed
  form, and `validate_spre.py` checks it against their code. Its `README.md` gives the
  validation steps.

## Stored data that v12 does not use

The scripts removed with v12 read these files, which are kept as the record of earlier
results: `adaptive_tensor.pkl`, `burgers_fine_tensor.pkl`, `cputime_L2.pkl`,
`cputime_L3.pkl`, `screen8.pkl`, `margin_scan.json`, `outcomes.json`,
`baselines_pooled_main.json`, `pooled_v1.2.1_run1.json`, the pooled files of versions 1.2
and 1.1 (`pooled_v1.2*.json` without `.1`, `pooled_v1.1*.json`), and the local parallel
runs `gromacs_gmx_d_threads*.pkl`, `lammps_mpi_n64_r*.pkl`.

## The article's environment

Ubuntu 24.04 on a 2-vCPU, 7 GB virtual machine; Python 3.11, NumPy 1.26.4, SciPy 1.13.1
(`requirements-article.txt`); OpenFOAM v1912, LAMMPS 20240207, GROMACS 2023.3, ABINIT 9.10.4,
Quantum ESPRESSO 6.7 (Ubuntu packages); on Hamilton, OpenFOAM v2406 and v2606, Quantum
ESPRESSO 7.6. Replays with Python 3.13, NumPy 2.5 and SciPy 1.18 give identical counts.
