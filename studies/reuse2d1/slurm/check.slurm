#!/bin/bash
# The 15-minute check on the test partition (RUNBOOK.md, step 3), before any other run:
# three grid runs end to end (PROTOCOL.md 11).  Submitted by:  bash slurm/submit.sh check
#SBATCH --job-name=re-check
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
set -o pipefail
# shellcheck source=/dev/null
source "${REUSE_CODE:?REUSE_CODE not set}/slurm/env.sh" || { echo "CHECK FAILED: environment"; exit 1; }
cd "$REUSE_CODE" || exit 1
{
    date; hostname
    echo "module $REUSE_OF_MODULE   WM_PROJECT_VERSION=${WM_PROJECT_VERSION:-unset}"
    module list 2>&1
    command -v blockMesh checkMesh simpleFoam python3
    python3 --version
    python3 -c "import numpy; print('numpy', numpy.__version__)"
} 2>&1 | tee "$REUSE_WORK/check_environment.txt"
n=$(python3 driver.py --sweep check --count)
for i in $(seq 0 $((n - 1))); do
    python3 driver.py --sweep check --index "$i"
done
python3 driver.py --check-report
