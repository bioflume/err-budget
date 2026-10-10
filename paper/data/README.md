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
| `burgers_sp_*.pkl`, `burgers_fine_tensor.pkl` | L2 at four ladder spacings, fine ladder | `exp34_spacing.py`, `exp33_fineladder.py` |
| `burgers_nu002_tensor.pkl` | L2 at nu = 0.02 | `exp37_calhardest.py` |
| `classify.pkl`, `screen8.pkl`, `eight_cache.pkl` | 8-parameter L3 runs (classes, screening, 200-start study) | `exp17`, `exp16`, `exp54` |
| `cputime_L2.pkl`, `cputime_L3.pkl` | CPU time of every L2 and L3 setting | `exp52_wallclock.py` |
| `adaptive_tensor.pkl` | L2 with SciPy's adaptive BDF, 320 settings | `exp51_adaptive.py` |
| `dfg_tensor.pkl` | OpenFOAM DFG 2D-1, 840 settings | `exp46_dfg_build.py` |
| `lammps_tensor_n{20,40,64}.pkl` | LAMMPS at 8,000 / 64,000 / 262,144 ions | `exp22_lammps.py` |
| `lammps_mpi_n64_r{1,2}.pkl` | LAMMPS 262,144 ions on 1 and 2 MPI ranks | `exp55_parallel.py` |
| `gromacs_tensor_gmx{,_d}.pkl` | GROMACS mixed and double precision, 225 settings | `exp24_gromacs.py` |
| `gromacs_gmx_d_threads{1,2}.pkl` | GROMACS gmx_d on 1 and 2 threads | `exp55_parallel.py` |
| `abinit_tensor_v11.pkl` | ABINIT, 336 settings, reference at 48 Ha | `exp27_abinit.py` |
| `qe_tensor.pkl` | Quantum ESPRESSO, 480 settings, reference at 90 Ry | `exp53_qe.py` |
| `outcomes.json`, `baselines.json`, `margin_scan.json` | scored outcomes behind Tables 5, 8 and Appendix B | `exp26`, `exp38`, `exp29` |

The text output of every script is in `../figures/`.
