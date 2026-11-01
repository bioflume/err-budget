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

# ---------------------------------------------------------------- local overrides
# config.local.sh (not in the repository) holds what differs on one installation, e.g.
# the absolute paths of codes installed in their own environments:
#   export ERRBUDGET_LMP=$HOME/micromamba/envs/lmp/bin/lmp
# It is read before the defaults below, which apply only to what it leaves unset.  A
# config.local.sh written for the copy of the repository from before the renaming
# (KNOBS_PW, KNOBS_LMP, KNOBS_MPIRUN) works too.
if [ -f "$(dirname "${BASH_SOURCE[0]}")/config.local.sh" ]; then
    source "$(dirname "${BASH_SOURCE[0]}")/config.local.sh"
fi
for _v in PW LMP MPIRUN; do
    _o="KNOBS_$_v"; _n="ERRBUDGET_$_v"
    if [ -z "${!_n:-}" ] && [ -n "${!_o:-}" ]; then export "$_n=${!_o}"; fi
done
unset _v _o _n

# ---------------------------------------------------------------- Python
# setup_env.sh creates a micromamba environment with Python, NumPy, SciPy,
# Quantum ESPRESSO, LAMMPS and the errbudget package.
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$HOME/micromamba}
export ERRBUDGET_ENV=${ERRBUDGET_ENV:-errbudget}
# Modules loaded by every job: the GROMACS build of setup_env.sh needs the C++ runtime
# of the compiler that built it (gcc/11.2 on Hamilton; without it "gmx_d genconf" exits
# with status 1, as every GROMACS run of J4 did in October 2026).  Empty: load none.
export ERRBUDGET_MODULES=${ERRBUDGET_MODULES-gcc/11.2}
load_modules() {
    [ -n "$ERRBUDGET_MODULES" ] || return 0
    type module > /dev/null 2>&1 || source /etc/profile.d/lmod.sh 2> /dev/null \
        || source /etc/profile.d/modules.sh 2> /dev/null || return 0
    module load $ERRBUDGET_MODULES > /dev/null 2>&1 \
        || echo "config.sh: could not load the modules $ERRBUDGET_MODULES" >&2
}
activate_python() {
    load_modules
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
        # Open MPI writes its session directory under $TMPDIR, which Slurm on Hamilton sets
        # to node-local /local/...; /local does not exist in the container, so every mpirun
        # failed with "mkdir /local: Read-only file system" (October 2026).  Bind /local (on
        # the compute nodes); where it is absent, use /tmp inside the container.
        local B="/nobackup,$ERRBUDGET_ROOT"
        if [ -d /local ]; then B="$B,/local"; else export SINGULARITYENV_TMPDIR=/tmp APPTAINERENV_TMPDIR=/tmp; fi
        export ERRBUDGET_OF_WRAP="singularity exec --bind $B $OF_SIF bash -c 'source $OF_BASHRC >/dev/null 2>&1; exec \"\$@\"' bash"
    else
        module load "$OF_MODULE"
        [ -n "${FOAM_BASH:-}" ] && source "$FOAM_BASH"
        export ERRBUDGET_OF_WRAP=""
    fi
}
# Parallel OpenFOAM runs (hamilton/common/ofmodern.py).  In October 2026 the J2 spin-up of
# n_arc 192 on 32 ranks ran at about 4 s per linear-solver iteration instead of a few ms,
# and the reference at n_arc 128 and dt 2.5e-4 on 16 ranks was as slow; both used up their
# time limits.  The ranks are therefore not bound to cores.  J2 also stops a run whose
# linear solves and their iterations take longer than ERRBUDGET_OF_MAX_S_PER_ITER seconds
# each (J2's default 0.1; normal: 2-6 ms) over ERRBUDGET_OF_WATCH_WINDOW seconds (900), and
# retries it from its last checkpoint at the next submission; set them here to change them
# (0: never stop).
export ERRBUDGET_OF_MPIRUN_OPTS=${ERRBUDGET_OF_MPIRUN_OPTS-"--bind-to none"}

# ---------------------------------------------------------------- other codes
# Quantum ESPRESSO and LAMMPS come from the conda environment; GROMACS in double
# precision is built from source by setup_env.sh (conda-forge builds are single
# precision, whose round-off floor is too high for the targets of the study).
export ERRBUDGET_PW=${ERRBUDGET_PW:-pw.x}
export ERRBUDGET_QE_SHIM=0
export ERRBUDGET_QE_PSEUDO_DIR=${ERRBUDGET_QE_PSEUDO_DIR:-$ERRBUDGET_ROOT/hamilton/J3_qe/pseudo}
export ERRBUDGET_LMP=${ERRBUDGET_LMP:-lmp}
export GMX_DOUBLE_PREFIX=${GMX_DOUBLE_PREFIX:-$WORK/gromacs-2024-double}

# J4 (November 2026): the memory limit of each LAMMPS run (GB; 0 = none) and its MPI launcher
# ({R} = ranks; ERRBUDGET_MPIRUN, if set, is the mpirun to use)
export ERRBUDGET_LMP_MEM_GB=${ERRBUDGET_LMP_MEM_GB:-6}
export ERRBUDGET_LMP_LAUNCH=${ERRBUDGET_LMP_LAUNCH:-"${ERRBUDGET_MPIRUN:-mpirun} --allow-run-as-root --oversubscribe -np {R}"}
