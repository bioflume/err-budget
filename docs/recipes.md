# Recipes: job files for common codes

These are starting points, not finished cases: the parameter lists, the output patterns
and the file layouts are the ones the article used for the same codes (its repository,
<https://github.com/bioflume/errbudget-paper>, has the complete drivers). Adapt the command,
the template files and the observable to your case, then `errbudget check job.yaml`.

The pattern is always the same: a `template_dir` holding the case with `{{placeholders}}`
in its input files, a `command` that runs the code inside the run directory, an
`observable` read from a log or output file, and a `cost` read from the code's own timing
when runs share a machine.

## OpenFOAM (steady flow, `simpleFoam`)

The article's 2D-1 cylinder study tuned the mesh, the SIMPLE residual control, the GAMG
pressure tolerance and the GAMG agglomeration (`nCellsInCoarsestLevel`), for the drag
coefficient. Put the case in `case/` with placeholders in `system/blockMeshDict`
(the cell counts) and `system/fvSolution`:

```
solvers
{
    p { solver GAMG; smoother GaussSeidel; tolerance {{lin_tol}}; relTol 0.1;
        nCellsInCoarsestLevel {{n_coarsest}}; }
    U { solver smoothSolver; smoother symGaussSeidel; tolerance {{lin_tol_U}}; relTol 0.1; }
}
SIMPLE { residualControl { p {{res_ctrl}}; U {{res_ctrl_U}}; } }
```

Placeholders can only take parameter values, so a value derived from another one (the U
tolerance ten times the p tolerance) is best written into the template by a one-line
script in the command, or declared as its own `cost-only` parameter if it has no effect.

```yaml
name: cylinder drag, OpenFOAM
command: >
  blockMesh > log.blockMesh &&
  simpleFoam > log.simpleFoam &&
  postProcess -func forceCoeffs -latestTime > log.post
template_dir: case
parameters:
  - {name: n_arc,      kind: discretisation, values: [8, 12, 16, 24, 32, 48, 64], start: 24}
  - {name: res_ctrl,   kind: nested-solver,  values: [1e-2, 1e-3, 1e-4, 1e-5, 1e-6], start: 1e-4}
  - {name: lin_tol,    kind: nested-solver,  values: [1e-4, 1e-5, 1e-6, 1e-7, 1e-8], start: 1e-6}
  - {name: n_coarsest, kind: cost-only,      values: [10, 50, 250], start: 10}
observable:
  name: Cd
  source: file
  path: postProcessing/forceCoeffs/0/coefficient.dat    # columns: Time Cd Cd(f) Cd(r) Cl ...
  column: 2
target: 3e-3
relative: true
cost: {source: file, path: log.simpleFoam, regex: "ExecutionTime = ([0-9.]+) s"}
failure:
  failed_if_regex: "FOAM FATAL"
  converged_if_regex: "SIMPLE solution converged"       # residualControl met before endTime
env: {OMP_NUM_THREADS: "1"}
timeout: 14400
keep_runs: false
```

The `converged_if_regex` line is what makes a run that stopped at `endTime` without meeting
`residualControl` a non-answer: without it, such runs would be fitted as if converged. The
patterns are searched in stdout and stderr and in the files the observable and the cost
come from, here `log.simpleFoam`; for a log that is neither, add `path: log.simpleFoam`
under `failure`.
For a parallel run, `decomposePar > log.decompose && mpirun -n 8 simpleFoam -parallel > log.simpleFoam`
and a cost of `ExecutionTime` times 8 (put the factor in a small script, or accept that
costs are compared only within one rank count). For a transient case (`pimpleFoam`), the
time step is a discretisation parameter and the observable is usually a time average
computed by a function object or by a short script in the command.

## LAMMPS (PPPM electrostatics)

The article tuned the real-space cutoff, the reciprocal-space accuracy and the PPPM
interpolation order for the Coulomb energy per ion, with the neighbour-list rebuild
frequency as a cost-only parameter. With variables on the command line nothing in the
input deck needs placeholders:

```yaml
name: rocksalt electrostatics, LAMMPS
command: >
  lmp -in {{job_dir}}/in.madelung
      -var cutoff {{cutoff}} -var accuracy {{accuracy}} -var order {{order}} -var every {{every}}
      -log log.lammps > stdout.lammps
parameters:
  - {name: cutoff,   kind: discretisation, values: [6.0, 8.0, 10.0, 12.0, 14.0], start: 10.0}
  - {name: accuracy, kind: discretisation, values: [1e-2, 1e-3, 1e-4, 1e-5, 1e-6], start: 1e-4}
  - {name: order,    kind: discretisation, values: [4, 5, 6, 7], start: 5}
  - {name: every,    kind: cost-only,      values: [1, 2, 5, 10], start: 1}
observable:
  name: ecoul_per_atom
  source: file
  path: log.lammps
  regex: "^\\s*ecoul/atom\\s*=\\s*([-0-9.eE+]+)"     # printed by the deck: print "ecoul/atom = ${e}"
target: 1e-5
relative: true
cost: {source: file, path: log.lammps, regex: "Loop time of ([0-9.]+)"}
failure:
  failed_if_regex: ["ERROR", "Killed"]
  path: stdout.lammps                       # the patterns are also searched here
timeout: 3600
```

In the deck, `pair_style lj/cut/coul/long ${cutoff}`, `kspace_style pppm ${accuracy}`,
`kspace_modify order ${order}` and `neigh_modify every ${every} delay 0`. The order is
declared as a discretisation parameter: LAMMPS adjusts the FFT grid to the requested
accuracy, so the order changes cost and, through the grid, the answer. Out-of-memory
kills (large grids at tight accuracies) are failed runs; give each run its own process, as
the command does.

## GROMACS (PME, mixed and double precision)

Tuned in the article: `rcoulomb`, `fourierspacing`, `pme_order` and `ewald_rtol` for the
potential energy of a water box. The `.mdp` file is the template:

```
rcoulomb        = {{rcoulomb}}
rvdw            = {{rcoulomb}}
fourierspacing  = {{fourierspacing}}
pme_order       = {{pme_order}}
ewald_rtol      = {{ewald_rtol}}
```

```yaml
name: water box potential energy, GROMACS
command: >
  gmx_d grompp -f md.mdp -c {{job_dir}}/conf.gro -p {{job_dir}}/topol.top -o run.tpr > log.grompp 2>&1 &&
  gmx_d mdrun -deffnm run -nt 4 > log.mdrun 2>&1 &&
  echo Potential | gmx_d energy -f run.edr -o energy.xvg > log.energy 2>&1
template_dir: case                      # holds md.mdp
parameters:
  - {name: rcoulomb,       kind: discretisation, values: [0.9, 1.0, 1.1, 1.2, 1.4], start: 1.0, h: inverse}   # nm; larger is finer
  - {name: fourierspacing, kind: discretisation, values: [0.20, 0.16, 0.12, 0.10, 0.08], start: 0.12}
  - {name: pme_order,      kind: cost-only,      values: [4, 5, 6], start: 4}
  - {name: ewald_rtol,     kind: discretisation, values: [1e-4, 1e-5, 1e-6], start: 1e-5}
observable:
  name: potential
  source: file
  path: energy.xvg
  column: 2                              # xvg comment lines start with # and @
  comment: "#@"
target: 1e-5
relative: true
cost: {source: file, path: log.mdrun, regex: "Time:\\s+[0-9.]+\\s+([0-9.]+)"}   # wall seconds
options:
  floor_probe:                           # needed for the mixed-precision gmx: measures its round-off floor
    - {rcoulomb: 1.4, fourierspacing: 0.06, pme_order: 8, ewald_rtol: 1e-8}
    - {rcoulomb: 1.4, fourierspacing: 0.05, pme_order: 8, ewald_rtol: 1e-9}
failure:
  failed_if_regex: ["Fatal error", "Segmentation fault"]
  path: [log.grompp, log.energy]            # log.mdrun is searched already (the cost's file)
timeout: 3600
```

`comment: "#@"` is read as a set of comment markers: lines starting with either are
skipped. With the mixed-precision `gmx`, the floor probe is what tells the tool that
targets of 1e-5 and below are at the round-off floor: the article found they can be met
only by chance, and the probe makes the tool decline them instead.

## Quantum ESPRESSO and ABINIT (plane-wave DFT)

The cutoff converges exponentially: declare `model: exponential`. The k-point mesh and
the SCF threshold are the other two accuracy parameters; the mixing parameter (QE) or the
number of line minimisations (ABINIT) only changes cost.

```yaml
name: silicon total energy, Quantum ESPRESSO
command: mpirun -n 4 pw.x -in pw.in > pw.out
template_dir: case                       # pw.in with the placeholders below; pseudopotentials beside it
parameters:
  - {name: ecutwfc,  kind: discretisation, model: exponential, values: [12, 16, 20, 24, 28, 32, 40, 48], start: 24}
  - {name: nk,       kind: discretisation, values: [2, 3, 4, 6], start: 4}
  - {name: conv_thr, kind: nested-solver,  values: [1e-4, 1e-6, 1e-8, 1e-10, 1e-12], start: 1e-8}
  - {name: mixing,   kind: cost-only,      values: [0.3, 0.5, 0.7], start: 0.7, h: inverse}
observable:
  name: energy
  source: file
  path: pw.out
  regex: "^!\\s+total energy\\s+=\\s+([-0-9.]+) Ry"
target: 1e-4
relative: true
cost: {source: file, path: pw.out, regex: "PWSCF\\s+:\\s+([0-9.]+)s CPU"}
failure:
  failed_if_regex: "Error in routine"
  converged_if_regex: "convergence has been achieved"
```

In `pw.in`: `ecutwfc = {{ecutwfc}}`, `conv_thr = {{conv_thr}}`, `mixing_beta = {{mixing}}`
and `K_POINTS automatic` followed by `{{nk}} {{nk}} {{nk}} 0 0 0`. For ABINIT the
parameters are `ecut`, `ngkpt`, `toldfe` and `nline`, the energy is read with
`regex: "etotal\\s+([-0-9.E+]+)"` and the cost from `Proc.\s+0 individual time \(sec\): cpu=\s+([0-9.]+)`;
an SCF cycle that stops at `nstep` without meeting `toldfe` is marked by
`unconverged_if_regex: "nstep=.*was not enough SCF cycles"` (adapt the wording to your
version's warning).

The force needs a converged k-point mesh and a reference that most studies cannot afford;
tune for the energy, or admit only force targets well above the reference uncertainty.

## An MPI code of your own

```yaml
command: mpirun -n {{ranks}} {{job_dir}}/bin/solve --mesh {{mesh}} --tol {{tol}} --out result.json
parameters:
  - {name: mesh,  kind: discretisation, values: [32, 64, 128, 256], start: 64}
  - {name: tol,   kind: nested-solver,  values: [1e-3, 1e-5, 1e-7, 1e-9], start: 1e-7}
  - {name: ranks, kind: cost-only,      values: [1, 2, 4, 8], start: 4}
observable: {name: flux, source: file, path: result.json, key: observables.flux}
target: 1e-3
cost: {source: file, path: result.json, key: core_seconds}
```

The rank count as a cost-only parameter is legitimate when it does not change the answer;
the tool checks this and reports a contradiction if domain decomposition moves the
observable beyond the resolution. Report core-seconds as the cost so that rank counts are
comparable.

## Checking a recipe

1. `errbudget check job.yaml` runs the starting setting once and prints the number it read and
   the cost. If it prints `FAILED`, look at `stdout.txt` and `stderr.txt` in the run
   directory it names: the command, the placeholders and the pattern are almost always the
   problem, and the fix takes a minute.
2. Make the coarsest and the finest value of each parameter run once by hand (edit `start`
   and run `errbudget check`), so that you know the candidate lists are feasible.
3. `errbudget run job.yaml`. The first lines show the runs as they happen; an interrupted
   tuning resumes with the same command.
