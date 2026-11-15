# Source me: loads OpenFOAM and Python for the airfoil study.
#   source slurm/env.sh            # both
#   source slurm/env.sh python     # Python only (login node: listing, submitting, analysis)
# Not `set -u`: module functions read unset variables.
_here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=/dev/null
source "$_here/config.sh"
if [ -r "$_here/config.local.sh" ]; then
    # shellcheck source=/dev/null
    source "$_here/config.local.sh"
fi
if ! type module >/dev/null 2>&1; then
    for _f in /etc/profile.d/lmod.sh /etc/profile.d/modules.sh; do
        # shellcheck source=/dev/null
        if [ -r "$_f" ]; then source "$_f"; break; fi
    done
fi

airfoil_python() {
    if [ -n "${AIRFOIL_PYTHON:-}" ]; then
        export PATH="$(dirname "$AIRFOIL_PYTHON"):$PATH"
    else
        type module >/dev/null 2>&1 && { module load "$AIRFOIL_PY_MODULE" >/dev/null 2>&1 || true; }
        if [ -f "$AIRFOIL_VENV/bin/activate" ]; then
            # shellcheck source=/dev/null
            source "$AIRFOIL_VENV/bin/activate"
        else
            echo "env.sh: no virtual environment at $AIRFOIL_VENV: run slurm/setup.sh first" >&2
            return 1
        fi
    fi
}

airfoil_openfoam() {
    # Some OpenFOAM modules set up the environment themselves; others (Hamilton's
    # openfoam/v2606) only set OPENFOAM_HOME and ask for OpenFOAM's bashrc to be
    # sourced.  Try, in order: AIRFOIL_OF_BASHRC (config.local.sh), FOAM_BASH, and the
    # usual places under OPENFOAM_HOME and WM_PROJECT_DIR.
    local m f
    for m in $AIRFOIL_PRE_MODULES; do module load "$m" || return 1; done
    module load "$AIRFOIL_OF_MODULE" || { echo "env.sh: cannot load $AIRFOIL_OF_MODULE" >&2; return 1; }
    for f in "${AIRFOIL_OF_BASHRC:-}" "${FOAM_BASH:-}" \
             "${OPENFOAM_HOME:-/nonexistent}/OpenFOAM/etc/bashrc" "${OPENFOAM_HOME:-/nonexistent}/etc/bashrc" \
             "${WM_PROJECT_DIR:-/nonexistent}/etc/bashrc"; do
        command -v simpleFoam >/dev/null 2>&1 && break
        if [ -n "$f" ] && [ -r "$f" ]; then
            # shellcheck source=/dev/null
            source "$f" >/dev/null 2>&1 || true
            export AIRFOIL_OF_SOURCED="$f"
        fi
    done
    command -v simpleFoam >/dev/null 2>&1 || {
        echo "env.sh: simpleFoam not on PATH after loading $AIRFOIL_OF_MODULE" >&2
        echo "env.sh: OPENFOAM_HOME=${OPENFOAM_HOME:-unset} WM_PROJECT_DIR=${WM_PROJECT_DIR:-unset};" \
             "set AIRFOIL_OF_BASHRC=<path to OpenFOAM's etc/bashrc> in slurm/config.local.sh" >&2
        return 1
    }
}

# one thread per allocated core for anything threaded (NumPy's BLAS in the driver)
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS MKL_NUM_THREADS=$OMP_NUM_THREADS
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$AIRFOIL_WORK/logs" "$AIRFOIL_WORK/runs" "$AIRFOIL_WORK/packs"

case "${1:-all}" in
    python) airfoil_python ;;
    *) airfoil_openfoam && airfoil_python ;;
esac
