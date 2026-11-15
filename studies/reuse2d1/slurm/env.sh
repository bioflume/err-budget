# Source me: loads OpenFOAM and Python for the reuse study.
#   source slurm/env.sh            # both
#   source slurm/env.sh python     # Python only (login node: listing, submitting)
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

reuse_python() {
    type module >/dev/null 2>&1 && { module load "$REUSE_PY_MODULE" >/dev/null 2>&1 || true; }
    if [ -f "$REUSE_VENV/bin/activate" ]; then
        # shellcheck source=/dev/null
        source "$REUSE_VENV/bin/activate"
    else
        echo "env.sh: no virtual environment at $REUSE_VENV (the airfoil study's; RUNBOOK.md, step 2)" >&2
        return 1
    fi
}

reuse_openfoam() {
    local m f
    for m in $REUSE_PRE_MODULES; do module load "$m" || return 1; done
    module load "$REUSE_OF_MODULE" || { echo "env.sh: cannot load $REUSE_OF_MODULE" >&2; return 1; }
    for f in "${REUSE_OF_BASHRC:-}" "${FOAM_BASH:-}" \
             "${OPENFOAM_HOME:-/nonexistent}/OpenFOAM/etc/bashrc" "${OPENFOAM_HOME:-/nonexistent}/etc/bashrc" \
             "${WM_PROJECT_DIR:-/nonexistent}/etc/bashrc"; do
        command -v simpleFoam >/dev/null 2>&1 && break
        if [ -n "$f" ] && [ -r "$f" ]; then
            # shellcheck source=/dev/null
            source "$f" >/dev/null 2>&1 || true
        fi
    done
    command -v simpleFoam >/dev/null 2>&1 || {
        echo "env.sh: simpleFoam not on PATH after loading $REUSE_OF_MODULE;" \
             "set REUSE_OF_BASHRC=<OpenFOAM's etc/bashrc> in slurm/config.local.sh" >&2
        return 1
    }
}

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$REUSE_WORK/logs" "$REUSE_WORK/runs" "$REUSE_WORK/tasks" "$REUSE_WORK/packs"

case "${1:-all}" in
    python) reuse_python ;;
    *) reuse_openfoam && reuse_python ;;
esac
