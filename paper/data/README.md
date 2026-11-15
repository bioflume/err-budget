# Stored solver output

These grids are used to **grade** the method, not by it: the method never sees more
than the settings it asks for. They are included so that every table can be checked
without rerunning the solvers; `../reproduce.sh` (without `--cached`) deletes them and
regenerates everything from scratch.

| file | what it is | script |
|---|---|---|
| `burgers_tensor.pkl` | L2, 625 settings | `exp00_build_L2.py` |
| `burgers2d_tensor.pkl` | L3, 432 settings | `exp07_2d.py` |
| `reactdiff_tensor.pkl` | L4, 576 settings | `exp20_reactdiff.py` |
| `cavity_tensor.pkl` | L5 driven cavity, 700 settings | `exp48_cavity.py` |
| `burgers_sp_*.pkl` | L2 at four spacings of the candidate values | `exp34_spacing.py` |
| `burgers_fine_tensor.pkl` | L2 on a fine ladder (not used by v12) | `exp33_fineladder.py`, removed |
| `burgers_nu002_tensor.pkl` | L2 at nu = 0.02 | `exp37_calhardest.py` |
| `classify.pkl`, `eight_cache.pkl` | 8-parameter L3 runs (classes, 200-start study) | `exp17`, `exp54` |
| `screen8.pkl` | 8-parameter L3 runs of the screening study (not used by v12) | `exp16`, removed |
| `cputime_L2.pkl`, `cputime_L3.pkl` | CPU time of every L2 and L3 setting (not used by v12) | `exp52_wallclock.py`, removed |
| `adaptive_tensor.pkl` | L2 with SciPy's adaptive BDF, 320 settings (not used by v12) | `exp51_adaptive.py`, removed |
| `dfg_tensor.pkl` | OpenFOAM DFG 2D-1, 840 settings | `exp46_dfg_build.py` |
| `lammps_tensor_n{20,40,64}.pkl` | LAMMPS at 8,000 / 64,000 / 262,144 ions | `exp22_lammps.py` |
| `lammps_mpi_n64_r{1,2}.pkl` | LAMMPS 262,144 ions on 1 and 2 MPI ranks (not used by v12) | `exp55_parallel.py` |
| `gromacs_tensor_gmx{,_d}.pkl` | GROMACS mixed and double precision, 225 settings | `exp24_gromacs.py` |
| `gromacs_gmx_d_threads{1,2}.pkl` | GROMACS gmx_d on 1 and 2 threads (not used by v12) | `exp55_parallel.py` |
| `abinit_tensor_v11.pkl` | ABINIT, 336 settings, reference at 48 Ha | `exp27_abinit.py` |
| `qe_tensor.pkl` | Quantum ESPRESSO, 480 settings, reference at 90 Ry | `exp53_qe.py` |
| `outcomes.json`, `baselines.json`, `margin_scan.json` | scored outcomes of earlier revisions (not used by v12) | `exp26`, `exp38`, `exp29`, removed |

The replay results (`pooled_*.json`, `baselines_*.json`, `*_released.json` and the other
`.json` files written by the scripts) are listed with their scripts in `../README.md`.
The printed output of the scripts is in `../figures/`.
