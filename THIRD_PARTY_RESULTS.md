# Running the method on solvers we did not write

Five production codes, all driven only through their input files.

| code | version | observable | ground truth | script |
|---|---|---|---|---|
| LAMMPS | 20240207 (Ubuntu 24.04) | electrostatic energy per ion, rocksalt lattice | exact (Madelung constant) | `exp22_lammps.py`, `exp39`, `exp43`, `exp55` |
| GROMACS | 2023.3 (`gmx`, `gmx_d`) | potential energy, 1,728 SPC waters | converged `gmx_d` run, uncertainty 2.3e-9 | `exp24_gromacs.py`, `exp44`, `exp55` |
| ABINIT | 9.10.4 | total energy and force, Si with a displaced atom | converged run (48 Ha, 8x8x8), uncertainty measured in `exp45` | `exp27_abinit.py`, `exp45` |
| Quantum ESPRESSO | 6.7 (Ubuntu 24.04) | total energy and force, same Si system | converged run (90 Ry, 8x8x8), uncertainty measured | `exp53_qe.py` |
| OpenFOAM | v1912 (Ubuntu 24.04), `simpleFoam` | drag, lift, pressure difference, DFG benchmark 2D-1 | published reference values | `exp46_dfg_build.py`, `exp46_dfg.py` |

## Packaging problems met, and how they were handled

- **Quantum ESPRESSO 6.7** (Ubuntu package) aborts with "buffer overflow
  detected" when it hashes the pseudopotential file: glibc's fortified
  `snprintf` rejects an off-by-one terminator write that an unfortified build
  performs silently. `paper/snfix.c`, compiled and preloaded by
  `qe_silicon.py`, restores the unfortified behaviour. No physics is affected.
- **OpenFOAM v1912** (Ubuntu package) aborts with an IOstream "sha1" error
  whenever a function object (`forces`, `forceCoeffs`, `probes`,
  `writeCellCentres`) is constructed. `openfoam_cylinder.py` therefore computes
  the drag and lift from the written `p` and `U` fields the way the `forces`
  function object does (pressure times face area vector, plus the wall shear
  nu U_t / d from the wall-adjacent cell), and the pressure difference from the
  wall pressure interpolated in angle.
- **Memory.** At 262,144 ions some LAMMPS settings need a PPPM grid larger than
  the 5 GB allowed per run; each LAMMPS setting runs in its own process so that
  an out-of-memory kill is recorded as a failed run.

## Reference uncertainties (admission rule: target >= 10 x uncertainty)

| code | energy | force |
|---|---|---|
| ABINIT, 48 Ha / 8x8x8 | 1.3e-6 (42 -> 48 Ha) | 4.9e-3 (8 -> 10 k-points) |
| Quantum ESPRESSO, 90 Ry / 8x8x8 | 1.2e-6 | 4.2e-3 |
| GROMACS `gmx_d` | 2.3e-9 | -- |

Force targets of 3e-2 and 1e-2 are therefore not admitted for either DFT code:
the force converges slowly in the k-point mesh, and a k-converged reference is
beyond the budget of this study.

All numbers are in `paper/figures/` after `paper/reproduce.sh`.
