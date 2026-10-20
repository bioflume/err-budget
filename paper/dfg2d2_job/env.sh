# Edit this file for your machine, then every script below sources it.
#
# OpenFOAM: an ESI-OpenCFD version (v1912 or later) with pimpleFoam, blockMesh,
# checkMesh on PATH.  On a cluster, load or source it, for example
#     module load OpenFOAM/v2312
#     source $FOAM_BASH
# The Ubuntu 24.04 package (apt install openfoam) works without sourcing.
#
# Python 3.9+ with numpy and scipy, and the errbudget package
# (pip install -e <path to errbudget-paper>), e.g.
#     module load Python/3.11
#     source ~/venvs/errbudget/bin/activate
export ERRBUDGET_ROOT=${ERRBUDGET_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}
export ERRBUDGET_DFG2_WORK=${ERRBUDGET_DFG2_WORK:-$ERRBUDGET_ROOT/paper/work_dfg2d2}
export PYTHON=${PYTHON:-python3}
# one OpenFOAM process per task is serial; keep BLAS single-threaded
export OMP_NUM_THREADS=1
