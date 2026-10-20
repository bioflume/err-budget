# Settings for running the heavy studies of the errbudget article on Durham's Hamilton.
# Edit this file once; every job script sources it.  Nothing else needs editing.
#
# Hamilton (Hamilton8): 128 cores and 246 GB per standard node; partitions
#   shared (default, jobs share nodes, 3 days), multi (whole nodes, 3 days),
#   long (7 days), test (15 minutes).

# ---------------------------------------------------------------- paths
# the errbudget-paper repository (this file lives in <repository>/hamilton)
export ERRBUDGET_ROOT=${ERRBUDGET_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
# scratch for the cases (large, temporary); results are small and are copied
# into $ERRBUDGET_ROOT/hamilton/results by collect.sh
export WORK=${WORK:-/nobackup/$USER/errbudget}

# ---------------------------------------------------------------- SLURM
# Hamilton does not normally need an account; add e.g. "-A myproject" if yours does.
export SB_EXTRA=${SB_EXTRA:-""}
export SB_SHARED=${SB_SHARED:-shared}      # serial and few-core jobs
export SB_WHOLE=${SB_WHOLE:-multi}         # whole-node jobs (timing studies)

# ---------------------------------------------------------------- Python
# setup_env.sh creates a micromamba environment with Python, NumPy, SciPy,
# Quantum ESPRESSO, LAMMPS and the errbudget package.
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$HOME/micromamba}
export ERRBUDGET_ENV=${ERRBUDGET_ENV:-errbudget}
activate_python() {
    eval "$("$HOME/bin/micromamba" shell hook -s bash)"
    micromamba activate "$ERRBUDGET_ENV"
    export PYTHONPATH="$ERRBUDGET_ROOT/paper:$ERRBUDGET_ROOT/src:$ERRBUDGET_ROOT/hamilton/common:${PYTHONPATH:-}"
    export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
}

# ---------------------------------------------------------------- OpenFOAM
# Either the OpenCFD container (default; any current ESI release works) or a
# Hamilton module.  The container is pulled by setup_env.sh.
export OF_MODE=${OF_MODE:-container}                 # container | module
export OF_VERSION=${OF_VERSION:-2406}
export OF_SIF=${OF_SIF:-$WORK/openfoam-v$OF_VERSION.sif}
export OF_BASHRC=${OF_BASHRC:-/usr/lib/openfoam/openfoam$OF_VERSION/etc/bashrc}
export OF_MODULE=${OF_MODULE:-openfoam}              # see: module avail openfoam
load_openfoam() {
    if [ "$OF_MODE" = "container" ]; then
        # every OpenFOAM executable is run as: singularity exec ... bash -c 'source bashrc; exec "$@"'
        export ERRBUDGET_OF_WRAP="singularity exec --bind /nobackup,$ERRBUDGET_ROOT $OF_SIF bash -c 'source $OF_BASHRC >/dev/null 2>&1; exec \"\$@\"' bash"
    else
        module load "$OF_MODULE"
        [ -n "${FOAM_BASH:-}" ] && source "$FOAM_BASH"
        export ERRBUDGET_OF_WRAP=""
    fi
}

# ---------------------------------------------------------------- other codes
# Quantum ESPRESSO and LAMMPS come from the conda environment; GROMACS in double
# precision is built from source by setup_env.sh (conda-forge builds are single
# precision, whose round-off floor is too high for the targets of the study).
export ERRBUDGET_PW=${ERRBUDGET_PW:-pw.x}
export ERRBUDGET_QE_SHIM=0
export ERRBUDGET_QE_PSEUDO_DIR=${ERRBUDGET_QE_PSEUDO_DIR:-$ERRBUDGET_ROOT/hamilton/J3_qe/pseudo}
export ERRBUDGET_LMP=${ERRBUDGET_LMP:-lmp}
export GMX_DOUBLE_PREFIX=${GMX_DOUBLE_PREFIX:-$WORK/gromacs-2024-double}
