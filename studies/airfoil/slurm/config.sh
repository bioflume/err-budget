# Settings of the airfoil study on Hamilton.  The only file to edit (RUNBOOK.md, step 3).
# Every script sources it through slurm/env.sh.

# where everything lives: code (this tarball), the run store, cases, logs, packs
export AIRFOIL_WORK=${AIRFOIL_WORK:-/nobackup/$USER/airfoil}
export AIRFOIL_CODE=${AIRFOIL_CODE:-$AIRFOIL_WORK/code}

# OpenFOAM: the exact module name that `module avail openfoam` prints for ESI v2606
# (RUNBOOK.md, step 2).  Some modules also need their bashrc sourced: env.sh does it
# when the module sets FOAM_BASH or leaves simpleFoam off the PATH.
export AIRFOIL_OF_MODULE=${AIRFOIL_OF_MODULE:-openfoam/v2606}
# OpenFOAM's etc/bashrc, if the module asks for it to be sourced and env.sh does not
# find it under $OPENFOAM_HOME or $WM_PROJECT_DIR (absolute path; set in config.local.sh)
export AIRFOIL_OF_BASHRC=${AIRFOIL_OF_BASHRC:-}
# modules to load before it, if the OpenFOAM module needs them (e.g. "gcc/11.2 openmpi")
export AIRFOIL_PRE_MODULES=${AIRFOIL_PRE_MODULES:-}

# Python for the driver (NumPy) and the live errbudget run (NumPy, SciPy, PyYAML).
# Default: the cluster's python module and a virtual environment made by setup.sh.
export AIRFOIL_PY_MODULE=${AIRFOIL_PY_MODULE:-python}
export AIRFOIL_VENV=${AIRFOIL_VENV:-$HOME/venvs/airfoil}

# Slurm
export AIRFOIL_PARTITION=${AIRFOIL_PARTITION:-shared}
export AIRFOIL_TEST_PARTITION=${AIRFOIL_TEST_PARTITION:-test}
export AIRFOIL_SB_EXTRA=${AIRFOIL_SB_EXTRA:-}          # e.g. "-A myproject" if an account is needed
# parallel runs only (Stage 0, run B): ranks are not bound to cores (see the article's hamilton/README.md)
export AIRFOIL_MPIRUN=${AIRFOIL_MPIRUN:-mpirun}
# --oversubscribe: the job is one task with --cpus-per-task CPUs, which OpenMPI counts as one
# slot (Stage 0, run B, 13 November: "not enough slots available")
export AIRFOIL_MPIRUN_OPTS=${AIRFOIL_MPIRUN_OPTS:---bind-to none --oversubscribe}

# cost guards in core-hours (PROTOCOL.md 12)
export AIRFOIL_GUARD_stage0=8 AIRFOIL_GUARD_stage1=400 AIRFOIL_GUARD_stage2=500
