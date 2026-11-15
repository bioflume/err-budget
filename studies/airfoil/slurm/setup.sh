#!/bin/bash
# One-time setup on a Hamilton login node (a few minutes; safe to re-run):
#   bash /nobackup/$USER/airfoil/code/slurm/setup.sh
# 1. lists the OpenFOAM modules (choose the ESI v2606 one in config.sh),
# 2. makes the Python virtual environment (NumPy for the driver; SciPy and PyYAML for
#    the live errbudget run),
# 3. loads OpenFOAM and prints its version, and copies the airFoil2D tutorial's
#    fvSolution for the record of setting (b),
# 4. lists the sweeps.
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=/dev/null
source "$here/config.sh"
# shellcheck source=/dev/null
if [ -r "$here/config.local.sh" ]; then source "$here/config.local.sh"; fi
if ! type module >/dev/null 2>&1; then
    # shellcheck source=/dev/null
    if [ -r /etc/profile.d/lmod.sh ]; then source /etc/profile.d/lmod.sh; fi
fi
set -o pipefail
mkdir -p "$AIRFOIL_WORK"/{runs,logs,packs,cases}

echo "== 1. OpenFOAM modules (config.sh has AIRFOIL_OF_MODULE=$AIRFOIL_OF_MODULE)"
module avail openfoam 2>&1 | tee "$AIRFOIL_WORK/module_avail.txt"
module avail OpenFOAM 2>&1 | tee -a "$AIRFOIL_WORK/module_avail.txt" >/dev/null

echo "== 2. Python environment $AIRFOIL_VENV"
if [ -z "${AIRFOIL_PYTHON:-}" ]; then
    module load "$AIRFOIL_PY_MODULE" || exit 1
    python3 --version
    if [ ! -f "$AIRFOIL_VENV/bin/activate" ]; then
        python3 -m venv "$AIRFOIL_VENV" || exit 1
    fi
    # shellcheck source=/dev/null
    source "$AIRFOIL_VENV/bin/activate"
    python3 -c "import numpy, scipy, yaml" 2>/dev/null || {
        pip install --quiet --upgrade pip
        pip install --quiet "numpy>=1.22" "scipy>=1.12" pyyaml || exit 1
    }
fi
python3 -c "import sys, numpy, scipy, yaml; print('python', sys.version.split()[0], 'numpy', numpy.__version__, 'scipy', scipy.__version__)"

echo "== 3. OpenFOAM"
# shellcheck source=/dev/null
source "$here/env.sh" || { echo "OpenFOAM did not load: check AIRFOIL_OF_MODULE (step 1 above) and AIRFOIL_OF_BASHRC in $here/config.local.sh"; exit 1; }
echo "simpleFoam: $(command -v simpleFoam)   WM_PROJECT_VERSION=${WM_PROJECT_VERSION:-unset}   sourced: ${AIRFOIL_OF_SOURCED:-nothing}"
simpleFoam -help 2>&1 | head -5 || true
tut="${FOAM_TUTORIALS:-$WM_PROJECT_DIR/tutorials}/incompressible/simpleFoam/airFoil2D/system/fvSolution"
if [ -f "$tut" ]; then
    cp "$tut" "$AIRFOIL_WORK/tutorial_airFoil2D_fvSolution"
    echo "copied $tut"
else
    echo "WARNING: airFoil2D tutorial not found at $tut"
fi

echo "== 4. Sweeps"
cd "$AIRFOIL_CODE"
for s in check stage0; do echo "$s: $(python3 driver.py --sweep "$s" --count) runs"; done
echo "Setup finished.  Next: bash $AIRFOIL_CODE/slurm/submit.sh check"
