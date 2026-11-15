"""Brings the results of the Hamilton jobs J1-J5 into the article.

    python exp60_hamilton.py hamilton_results.tar.gz

Unpacks the archive written by hamilton/collect.sh into hamilton/results, replays H1-H3
with the frozen and the released program (exp62_heldout.py), prints the summaries of J1-J4
and regenerates the tables (exp70_tables.py).  J1: 2D-1 with OpenFOAM v2406 and its function
object (S-5.1, S-5.6); J2: the refined 2D-2 study (Section 7.2, S-5.2); J3: Quantum
ESPRESSO 7.6 (S-5.6); J4: MPI ranks (S-5.6); J5: H3 (Section 6.2, Table 5).
"""
import os
import subprocess
import sys
import tarfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
if len(sys.argv) < 2:
    raise SystemExit(__doc__)
with tarfile.open(sys.argv[1]) as tf:
    tf.extractall(os.path.join(ROOT, "hamilton"), filter="data")
res = os.path.join(ROOT, "hamilton", "results")
for j in ("j1", "j2", "j3", "j4", "j5"):
    p = os.path.join(res, f"{j}.txt")
    print(f"\n===== {j.upper()} =====")
    print(open(p).read() if os.path.exists(p) else "(missing)")
py = sys.executable
for mode in ([], ["--current"]):
    print(f"\n===== held-out replay {'(released 1.2.1)' if mode else '(frozen 1.2)'} =====", flush=True)
    subprocess.run([py, "-W", "ignore", os.path.join(HERE, "exp62_heldout.py")] + mode, check=True)
subprocess.run([py, os.path.join(HERE, "exp70_tables.py")], check=True, stdout=subprocess.DEVNULL)
print("\nTables regenerated in paper/figures (tab_heldout.tex now includes H3).")
