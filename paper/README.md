# Reproducing the article

`reproduce.sh` regenerates every number, table and figure in the article **from
scratch**. It deletes `data/` and runs every solver again, including LAMMPS, GROMACS and
ABINIT. Nothing is read from saved output unless you ask for it.

```bash
./reproduce.sh            # from scratch (about 10 hours on 2 cores)
./reproduce.sh --cached   # re-analyse the stored solver output (about 2 hours)
```

Each step writes its output to `figures/<step>.txt`. The table below maps every result in
the article to the script that produces it.

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

- **Machine-dependent numbers:** the costs of every production code (wall-clock time for
  LAMMPS and GROMACS, CPU time for OpenFOAM, ABINIT and Quantum ESPRESSO), the penalties
  derived from them, the CPU-time replay of `exp52`, the 8-parameter study `exp54`
  (which uses work units but whose recommendations are graded by them), and running times.
- **Machine-independent numbers:** everything else, which reproduces exactly.
- **One exception, the L2 cost model.** Its work-unit cost counts GMRES iterations. At loose
  tolerances this count can differ by one between SciPy builds, which changes L2 costs by up
  to 0.25%. In our checks this changed no number in the article.

## Where each result of the current article comes from (revision 3, version 1.2.1)

| Result | Script | Output |
|---|---|---|
| Table 3 (pooled), Table 4 (ablation), totals in the text | `pooled.py v1.2.1 v1.2 v1.1 v1.2-nojoint …` then `exp70_tables.py` | `data/pooled_*.json`, `figures/tab_pooled.tex`, `figures/tab_ablation.tex` |
| Pooled totals, low-order acceptances, step-cap hits | `summ.py v1.2.1 v1.2 v1.1` | stdout |
| Case registry (every tensor, declaration, reference, admission rule) | `cases.py` | — |
| Section 4.1: L2–L5 convergence flags | `exp57_flags.py` | `data/flags_L*.pkl` |
| Table 5 and Section 6.2: held-out H1, H2 (H3 after the cluster run) | `exp61_heldout_build.py H1/H2`, `exp62_heldout.py` (frozen, from the tag or `frozen_src/`), `exp62_heldout.py --current` | `figures/rev3_heldout*.txt`, `figures/tab_heldout.tex` |
| Section 6.2: post-hoc check of the amendment | `posthoc_v121.py` | — |
| Table 6 and Section 6.3: baselines | `exp72_baselines_pooled.py` | `figures/rev3_baselines_pooled.txt`, `figures/tab_baselines.tex` |
| Section 6.1: decline reasons | (inline in the article's text; reproduced by `run_protocol(...)["declined"]`) | — |
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
| LAMMPS table, Fig. 3: LAMMPS at 8,000 / 64,000 / 262,144 ions | `exp22_lammps.py` (`KNOBS_LMP_N=20,40,64`), `make_lammps_figure.py` | `sec81_lammps_n*.txt` |
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
| Section 8.6: parallel execution | `exp55_parallel.py` (`KNOBS_RANKS`, `KNOBS_THREADS`) | `sec8x_parallel.txt` |
| Behavioural claims (31 tests) | `pytest ../tests` | — |

Exhaustive grids are built before the analyses that use them:
`exp00_build_L2.py`, `exp07_2d.py` (L3), `exp20_reactdiff.py` (L4),
`exp33_fineladder.py`, `exp34_spacing.py` and `exp37_calhardest.py`. These grids are used
only to grade the method. The method itself never sees more than the configurations it
chooses to evaluate.

Scripts in `legacy/` belong to earlier versions of the article and are not needed.

## Version 1.1 (revision 2)

The method changed in three places (see `src/knobs/_core.py`):
`unfit="order1"` (a parameter that cannot be fitted is charged its one-step change
extrapolated at first order), the `noise_floor` check, and `models=[..."exponential"...]`.
Every table was regenerated with version 1.1. The stored solver output in `data/` and
the text output of every script in `figures/` are included so that results can be
inspected without rerunning; `./reproduce.sh` regenerates them from scratch.
