#!/bin/bash
# The 15-minute check on the test partition (RUNBOOK.md, step 4), before any other run:
#  1. checkMesh on the coarsest and finest meshes of the study, the reference and the
#     n_surf 512 mesh (PROTOCOL.md 14, A2); a mesh that fails a check fails the job;
#  2. two tiny cases end to end (mesh, checkMesh, simpleFoam, forces, y+);
#  3. the airFoil2D tutorial file.
# Submitted by:  bash slurm/submit.sh check
#SBATCH --job-name=af-check
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
set -o pipefail
# shellcheck source=/dev/null
source "${AIRFOIL_CODE:?AIRFOIL_CODE not set}/slurm/env.sh" || { echo "CHECK FAILED: environment"; exit 1; }
cd "$AIRFOIL_CODE" || exit 1
{
    date; hostname
    echo "module $AIRFOIL_OF_MODULE   WM_PROJECT_VERSION=${WM_PROJECT_VERSION:-unset}"
    module list 2>&1
    command -v simpleFoam checkMesh python3
    python3 --version
} 2>&1 | tee "$AIRFOIL_WORK/check_environment.txt"
n=$(python3 driver.py --sweep meshcheck --count)
for i in $(seq 0 $((n - 1))); do
    python3 driver.py --sweep meshcheck --index "$i" --mesh-only
done
for i in 0 1; do
    python3 driver.py --sweep check --index "$i"
done
python3 analysis/check_smoke.py --work "$AIRFOIL_WORK"
