# Settings of the reuse study on Hamilton.  Every script sources it through slurm/env.sh;
# config.local.sh (not in the tarball, not overwritten by later uploads) can override any.

# where everything lives: code (this tarball), the run store, cases, logs, packs
export REUSE_WORK=${REUSE_WORK:-/nobackup/$USER/reuse2d1}
export REUSE_CODE=${REUSE_CODE:-$REUSE_WORK/code}

# OpenFOAM ESI v2606, the module of the airfoil study (its RUNBOOK.md, step 2).  Hamilton's
# module only sets OPENFOAM_HOME; env.sh sources $OPENFOAM_HOME/OpenFOAM/etc/bashrc itself.
export REUSE_OF_MODULE=${REUSE_OF_MODULE:-openfoam/v2606}
export REUSE_OF_BASHRC=${REUSE_OF_BASHRC:-}          # absolute path, if env.sh does not find it
export REUSE_PRE_MODULES=${REUSE_PRE_MODULES:-}      # modules to load first, if any

# Python with NumPy: the airfoil study's virtual environment (made by its slurm/setup.sh)
export REUSE_PY_MODULE=${REUSE_PY_MODULE:-python}
export REUSE_VENV=${REUSE_VENV:-$HOME/venvs/airfoil}

# Slurm
export REUSE_PARTITION=${REUSE_PARTITION:-shared}
export REUSE_TEST_PARTITION=${REUSE_TEST_PARTITION:-test}
export REUSE_SB_EXTRA=${REUSE_SB_EXTRA:-}             # e.g. "-A myproject" if an account is needed

# cost guard in core-hours (PROTOCOL.md 11)
export REUSE_GUARD=200
