#!/bin/bash
# One-off setup on a Hamilton login node (about 30-60 minutes, mostly downloads
# and the GROMACS build).  Re-running it is safe: finished steps are skipped.
#
#   cd knobs-1.2.1/hamilton && bash setup_env.sh
set -e
cd "$(dirname "$0")"
source ./config.sh
mkdir -p "$WORK" "$HOME/bin"

# 1. micromamba (user space, no root)
if [ ! -x "$HOME/bin/micromamba" ]; then
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C "$HOME" bin/micromamba
fi
eval "$("$HOME/bin/micromamba" shell hook -s bash)"

# 2. Python, Quantum ESPRESSO (current release), LAMMPS with MPI, scikit-fem, pyamg
if ! micromamba env list | grep -q "^ *$KNOBS_ENV "; then
    micromamba create -y -n "$KNOBS_ENV" -c conda-forge \
        python=3.11 numpy scipy matplotlib pyamg scikit-fem \
        "qe>=7.3" "lammps=*=*openmpi*" openmpi cmake make
fi
micromamba activate "$KNOBS_ENV"
# revision 4: SymPy (exact values of the manufactured problems L3-L5, used by the
# pooled analyses) and pytest (the test suite), also into an existing environment
python -c "import sympy, pytest" 2>/dev/null || micromamba install -y -n "$KNOBS_ENV" -c conda-forge sympy pytest
pip install -e "$KNOBS_ROOT" >/dev/null
echo "pw.x:  $(command -v pw.x)";  pw.x -h 2>&1 | head -1 || true
LMP=$(command -v lmp || command -v lmp_mpi || true)
echo "lmp:   $LMP"
[ -z "$LMP" ] && echo "WARNING: no LAMMPS executable found; set KNOBS_LMP in config.sh"

# 3. OpenFOAM container (skip if OF_MODE=module)
if [ "$OF_MODE" = "container" ] && [ ! -f "$OF_SIF" ]; then
    singularity pull "$OF_SIF" "docker://opencfd/openfoam-default:$OF_VERSION"
fi
load_openfoam
echo "OpenFOAM check:"; ${KNOBS_OF_WRAP} simpleFoam -help 2>&1 | head -3 || true

# 4. GROMACS in double precision with thread-MPI (about 20 minutes)
if [ ! -x "$GMX_DOUBLE_PREFIX/bin/gmx_d" ]; then
    module load gcc 2>/dev/null || true
    cd "$WORK"
    [ -f gromacs-2024.4.tar.gz ] || curl -LO https://ftp.gromacs.org/gromacs/gromacs-2024.4.tar.gz
    rm -rf gromacs-2024.4 && tar xzf gromacs-2024.4.tar.gz && cd gromacs-2024.4
    mkdir -p build && cd build
    cmake .. -DGMX_DOUBLE=ON -DGMX_BUILD_OWN_FFTW=ON -DGMX_MPI=OFF -DGMX_GPU=OFF \
             -DCMAKE_INSTALL_PREFIX="$GMX_DOUBLE_PREFIX" >/dev/null
    make -j 8 >/dev/null && make install >/dev/null
fi
echo "gmx_d: $GMX_DOUBLE_PREFIX/bin/gmx_d"
echo
echo "Setup finished.  Next: bash run_all.sh smoke   (a 15-minute check on the test partition)"
