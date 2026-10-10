#!/usr/bin/env bash
# Regenerate every number, table and figure in the article FROM SCRATCH.
#
#   ./reproduce.sh            delete all cached solver output and rebuild it
#   ./reproduce.sh --cached   reuse whatever is already in data/
#
# Requirements: Python 3.11, numpy 1.26.4, scipy 1.13.1 (>= 1.12 is required),
# matplotlib, sympy, pytest; LAMMPS with its Python module, GROMACS (gmx and
# gmx_d), ABINIT with its pseudopotential directory /usr/share/abinit/psp.
# Revision 2 adds OpenFOAM (simpleFoam) and Quantum ESPRESSO (pw.x) and gcc.
# On Ubuntu 24.04: apt install lammps python3-lammps gromacs abinit openfoam
# quantum-espresso gcc.
#
# Wall-clock costs (LAMMPS, GROMACS and ABINIT timings) depend on the machine;
# OpenFOAM, Quantum ESPRESSO and the 8-parameter L3 study (CPU time) also
# depend on it; every other number is machine-independent.  About 16 CPU-hours
# on 2 cores.
set -e
cd "$(dirname "$0")"
mkdir -p data figures logs
if [ "$1" != "--cached" ]; then
  echo "== removing cached solver output"
  rm -f data/*.pkl data/*.json
fi
step(){ name=$1; shift; echo "== $name"; "$@" 2>&1 | tee "figures/$name.txt"; }

echo "######## 1. Build the exhaustive grids used to grade the method"
step build_L2        python exp00_build_L2.py
step build_L3        python exp07_2d.py
step build_L4        python exp20_reactdiff.py
step build_fine      python exp33_fineladder.py
step tab07_spacing   python exp34_spacing.py
step sec77_hardest   python exp37_calhardest.py
for N in 20 40 64; do step sec81_lammps_n$N env KNOBS_LMP_N=$N python exp22_lammps.py; done
step sec82_gromacs_d python exp24_gromacs.py
step sec82_gromacs_mixed env KNOBS_GMX=gmx python exp24_gromacs.py
step sec83_abinit    python exp27_abinit.py
step build_dfg       python exp46_dfg_build.py
step sec5x_cavity    python exp48_cavity.py
step sec8x_qe        python exp53_qe.py

echo "######## 2. Results"
step sec31_classes   python exp17_classify.py
step sec5_masking    python exp40_section5.py
step sec5_floor      python exp25_referee.py
step tab03_transfer  python exp32_conditions.py
step sec52_mechanism python exp31_mechanism.py
step tab04_ablations python exp41_ablations.py
step tab05_outcomes  python exp26_outcomes.py
step tab06_rho       python exp28_separability.py
step tab07_whyfine   python exp35_whyfine.py
step tab07_gated     python exp36_gated.py
step tab08_baselines python exp38_baselines.py
step sec76_margin    python exp29_margin.py
step tab09_transfer  python exp30_transfer.py
step tab10_calibration python exp42_calibration.py
step tab10_screened  python exp19_screened2d.py
step sec34_screening python exp16_screen8.py
step sec81_baselines python exp39_lmp_baselines.py
step sec81_misses    python exp43_lammps_misses.py
step sec82_floor     python exp44_gromacs_floor.py
step sec83_reference python exp45_abinit_reference.py
step sec84_learned   python exp23_mlhybrid.py
step sec8x_dfg       python exp46_dfg.py
step tab03_all       python exp47_transfer.py
step tab_misdeclare  python exp49_misdeclare.py
step app_margin_held python exp50_margin_heldout.py
step sec9x_adaptive  python exp51_adaptive.py
step tab05_cputime   python exp52_wallclock.py
step sec7x_eight     python exp54_eight.py
step sec8x_parallel  python exp55_parallel.py
# DFG 2D-2: the grid takes about 24 CPU-hours; build it with dfg2d2_job/ (see its README)
if [ -f data/dfg2d2_tensor.pkl ]; then step sec8x_dfg2d2 python exp56_dfg2d2.py; fi

echo "######## 3. Figures"
KNOBS_FIG_DIR=figures/ python make_solution_figure.py
KNOBS_FIG_DIR=figures/ python make_lammps_figure.py

echo "######## 4. Behavioural claims"
python -m pytest ../tests -q
ls -la figures/
