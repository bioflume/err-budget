# Settings for running the heavy studies of the knobs article on Durham's Hamilton.
# Edit this file once; every job script sources it.  Nothing else needs editing.
#
# Hamilton (Hamilton8): 128 cores and 246 GB per standard node; partitions
#   shared (default, jobs share nodes, 3 days), multi (whole nodes, 3 days),
#   long (7 days), test (15 minutes).

# ---------------------------------------------------------------- paths
# the unpacked knobs package (this file lives in <package>/hamilton)
export KNOBS_ROOT=${KNOBS_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
# scratch for the cases (large, temporary); results are small and are copied
# into $KNOBS_ROOT/hamilton/results by collect.sh
export WORK=${WORK:-/nobackup/$USER/knobs}

# ---------------------------------------------------------------- SLURM
# Hamilton does not normally need an account; add e.g. "-A myproject" if yours does.
export SB_EXTRA=${SB_EXTRA:-""}
export SB_SHARED=${SB_SHARED:-shared}      # serial and few-core jobs
export SB_WHOLE=${SB_WHOLE:-multi}         # whole-node jobs (timing studies)

# ---------------------------------------------------------------- Python
# setup_env.sh creates a micromamba environment with Python, NumPy, SciPy,
# Quantum ESPRESSO, LAMMPS and the knobs package.
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$HOME/micromamba}
export KNOBS_ENV=${KNOBS_ENV:-knobs}
activate_python() {
    eval "$("$HOME/bin/micromamba" shell hook -s bash)"
    micromamba activate "$KNOBS_ENV"
    export PYTHONPATH="$KNOBS_ROOT/paper:$KNOBS_ROOT/src:$KNOBS_ROOT/hamilton/common:${PYTHONPATH:-}"
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
        export KNOBS_OF_WRAP="singularity exec --bind /nobackup,$KNOBS_ROOT $OF_SIF bash -c 'source $OF_BASHRC >/dev/null 2>&1; exec \"\$@\"' bash"
    else
        module load "$OF_MODULE"
        [ -n "${FOAM_BASH:-}" ] && source "$FOAM_BASH"
        export KNOBS_OF_WRAP=""
    fi
}

# ---------------------------------------------------------------- other codes
# Quantum ESPRESSO and LAMMPS come from the conda environment; GROMACS in double
# precision is built from source by setup_env.sh (conda-forge builds are single
# precision, whose round-off floor is too high for the targets of the study).
export KNOBS_PW=${KNOBS_PW:-pw.x}
export KNOBS_QE_SHIM=0
export KNOBS_QE_PSEUDO_DIR=${KNOBS_QE_PSEUDO_DIR:-$KNOBS_ROOT/hamilton/J3_qe/pseudo}
export KNOBS_LMP=${KNOBS_LMP:-lmp}
export GMX_DOUBLE_PREFIX=${GMX_DOUBLE_PREFIX:-$WORK/gromacs-2024-double}
