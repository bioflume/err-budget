#!/usr/bin/env bash
# Rebuilds the stored runs with the solvers, runs the analyses that need a solver, and then
# replays everything (replay_all.sh).
#
#   bash reproduce.sh            # delete the stored runs in data/ and rebuild them
#   bash reproduce.sh --cached   # keep what is already in data/
#
# Needs Python 3.11 with requirements-article.txt, LAMMPS with its Python module, GROMACS
# (gmx and gmx_d), ABINIT with its pseudopotentials in /usr/share/abinit/psp, OpenFOAM v1912
# (simpleFoam), Quantum ESPRESSO (pw.x) and gcc; on Ubuntu 24.04: apt install lammps
# python3-lammps gromacs abinit openfoam quantum-espresso gcc.  About 16 CPU-hours on two
# cores.  The coarser 2D-2 study (dfg2d2_job/, about 40 CPU-hours) and the cluster jobs
# (../hamilton/) are built separately.  Wall-clock costs depend on the machine; every
# other number does not.
set -e
cd "$(dirname "$0")" || exit 1
mkdir -p data figures
if [ "${1:-}" != "--cached" ]; then
    rm -f data/*.pkl data/*.json
fi
step() { name=$1; shift; echo "== $name"; "$@" 2>&1 | tee "figures/$name.txt"; }

# The stored runs
step build_L2        python exp00_build_L2.py
step build_L3        python exp07_2d.py
step build_L4        python exp20_reactdiff.py
step sec5x_cavity    python exp48_cavity.py
for N in 20 40 64; do step sec81_lammps_n$N env ERRBUDGET_LMP_N=$N python exp22_lammps.py; done
step sec82_gromacs_d python exp24_gromacs.py
step sec82_gromacs_mixed env ERRBUDGET_GMX=gmx python exp24_gromacs.py
step sec83_abinit    python exp27_abinit.py
step build_dfg       python exp46_dfg_build.py
step sec8x_qe        python exp53_qe.py
step rev3_flags      python exp57_flags.py
[ -f data/heldout_H1.pkl ] || step rev3_build_H1 python exp61_heldout_build.py H1
[ -f data/heldout_H2.pkl ] || step rev3_build_H2 python exp61_heldout_build.py H2

# Analyses that run a solver
step sec82_floor     python exp44_gromacs_floor.py
step sec83_reference python exp45_abinit_reference.py
step rev3_abinit_scf python exp68_abinit_scf.py
step rev3_fine       python exp69_dfg_fine.py
step dfg_plateau     python exp59_dfg_plateau.py
step pppm_mechanism  python exp74_pppm_mechanism.py
step abinit_reference python exp75_abinit_reference.py

bash replay_all.sh
