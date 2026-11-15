# Reproduction check of v12

Machine: Ubuntu 22.04 (aarch64) in a virtual machine on the author's Mac, 4 cores, 3 GB;
a fresh environment with Python 3.11.17 and `requirements-article.txt` (NumPy 1.26.4, SciPy
1.13.1, Matplotlib 3.11.2, SymPy 1.14.0, pandas 3.0.6); `pip install -e .`. Date: 15 November
2026, branch `paper`. The virtual machine stops any command after three minutes, so the
replay ran here step by step, not as one run of `replay_all.sh`.

Status:
- **exact**: the regenerated file is byte-identical to the stored one or to the file of v12.
- **digits**: every number printed in the article is reproduced, but stored values differ in
  later digits.
- **stored**: not rerun here; the table or figure is regenerated exactly from the stored
  results of the script.
- **not reproduced**: no script produces the number.

The engine's searches are bit-identical: all 56 graded test cases of `pooled.py v1.2.1`
were searched again from every starting setting (and the 20 test cases of L5 and ABINIT of
`v1.2.1-nolocal` and `v1.2.1-constant`, which hold the ablation's misses). Every count,
cost ratio and run count equals the stored one, and the worst error ratio agrees to 10⁻⁹.

## Article

| Result | Status | How checked |
|---|---|---|
| Tables 1 and 2, Fig. 2 | — | no data |
| Table 3 | exact (counts); stored (reference uncertainties) | `exp76_accounting.py` rewrote `data/accounting_v6.json` byte for byte; the ABINIT and GROMACS uncertainties need the solvers (`reproduce.sh`) |
| Table 4 | exact | `exp70_tables.py` writes `tab_pooled_compact.tex` byte-identical to v12; the replay it reads is reproduced (above) |
| Table 5 | exact | `exp62_heldout.py` and `--current` rewrote both replay files byte for byte (frozen program verified by its hash); `tab_preregistered.tex` byte-identical |
| Fig. 1 | exact (data) | `extract_2d1_l5.py` rewrote `figures/data_2d1_l5.json` byte for byte; the printed values equal `figures/values_2d1_l5.txt`; the drawing differs from v12 only in the rendering of the fonts |
| Fig. 3 | exact (data) | as Fig. 1 |
| Fig. 4 | exact (data) | drawn from Table S-4.2 (exact); the printed values equal `figures/values_baselines.txt` |

## Supplement

| Result | Status | How checked |
|---|---|---|
| Tables S-1.1 to S-1.3, S-2.1, S-2.2 | — | from the code; the files cited by line (`grids.py`, `qe_silicon.py`, `abinit_silicon.py`, `src/errbudget/problems/`) are unchanged |
| Fig. S-2.1 | exact | pixel-identical to v12 |
| Table S-3.1 | exact | byte-identical; replay reproduced |
| Table S-3.2 | exact | byte-identical from the stored ablation; 20 of its test cases searched again, identical |
| Table S-3.3 | stored | `exp66_estimate_v12.py` not rerun; the table is byte-identical from `data/estimate_v12.json` |
| Table S-3.4 | stored | every cell equals `data/spacing_released.json` and `spacing_estimate_released.json` (the v12 run); `exp34_spacing.py` (about 9 minutes and its flag reruns) not rerun here |
| Table S-3.5 | stored | the stored output, relabelled by the new `exp65_misdeclare_v12.py`, equals v12 byte for byte; not rerun (long) |
| Table S-4.1 | stored | byte-identical from `data/transfer_floor.json` |
| Tables S-4.2, S-4.3 | stored | byte-identical from `data/baselines_v6.json`, `baselines_lmp_gmx.json`; `exp72_baselines_pooled.py` not rerun |
| Table S-4.4 | digits | `exp37_calhardest.py` and `exp30_transfer.py` rerun: same settings, acceptances, members met and worst ratios to the printed digits; the errors of the reruns of the solver differ from the stored ones from the 4th significant digit (at most 4×10⁻³ relative, at errors near 10⁻⁴) |
| Tables S-4.5, S-4.6 | stored | from `studies/pre_comparison/results/`; the comparison needs network access (`replay_all.sh --pre`) and was not rerun |
| Fig. S-4.1 | exact (data) | redrawn from the stored CSVs; same content as v12, fonts differ |
| Figs. S-5.1, S-5.2 | exact (data) | as Fig. 1 |
| Fig. S-5.3 | exact (data) | redrawn from `data/lammps_tensor_n20.pkl`; fonts differ |
| S-2.5 remark, S-5.2 (new `exp83`, `exp84`) | digits, with two exceptions | every number reproduced except the Strouhal limit at the observed order (6.8×10⁻⁴ here, 6.7×10⁻⁴ in the supplement) and the change of the lift from n_arc 24 to 32 (−0.09 to +0.11 ε over the four missed settings, −0.09 to +0.09 ε in the supplement, which holds for the two settings at GAMG 10⁻⁸ that 54 of the 64 searches return). "The mesh error there, 1.3–1.44 ε" is read as the change of the drag from n_arc 32 to 96 (1.25–1.33 ε) and its error at 32 against the reference (1.44–1.45 ε) |
| Tables S-6.1, S-6.2 | digits | `studies/reuse2d1/analysis/reuse.py` rerun from the raw pack: the printed results are identical; 19 floats of `reuse_summary.json` differ in the last one or two digits (order of summation) |
| S-7 | digits | `studies/airfoil/analysis/stage1.py` rerun from the raw pack: `stage1_summary.json` identical except two sums in the last digit |
| The other analyses rerun here | exact | `exp77_reltol_2d1.py`, `exp78_gromacs_floor_v121.py`: stored JSON rewritten byte for byte; `exp56_dfg2d2.py`: printed output equals `figures/sec8x_dfg2d2.txt`; tests: 52 passed |

## Not reproduced: numbers without a script

No script in the repository prints these numbers of the supplement; they appear to come
from analyses of the writing sessions:

- S-5.1: the shares of momentum solves stopped by their tolerance alone over the 608
  converged runs (26% at 100 t/r_c = 0.1, 94% at 1). `../hamilton/results/j1.txt` prints them
  over all 840 runs (26% and 93%); the 48 runs with residual histories are in
  `../hamilton/results/j1_reltol.json`, but no script tallies them. (The change with `relTol`
  0, at most 3.9×10⁻³ and 9.6×10⁻² for the lift, is in `j1.txt`.)
- S-5.2: "22 of 54 runs diverged" and Co_max = 27.8 n_arc Δt. The Courant numbers are in
  `hamilton/results/j2_tensor.pkl` ("courant"), but no script tallies them.
- S-3.2: the LAMMPS errors along the five cutoffs, (3.1, 3.6, 0.77, 2.4, 0.95)×10⁻²; the ABINIT
  estimate at the finest setting (1.15 and 3.5 times the target) and the "about 50 times".

## Not run here

- The analyses that run a solver (`reproduce.sh`): `exp44`, `exp45`, `exp59`, `exp68`,
  `exp69`, `exp74`, `exp75`. Their stored outputs are what the article quotes.
- The long replay steps: the ablation beyond the 20 test cases above, `exp72`, `exp64`,
  `exp65`, `exp66`, `exp67`, `exp73`, `exp79`, `exp28`, `exp81`, `exp17`, `exp34`, `exp35`
  and `exp54`.
- The run time of `replay_all.sh` as a whole. It is to be measured on the author's Mac
  with `time bash paper/replay_all.sh`; then `git status paper/data paper/figures` lists
  every stored result the replay changed.
