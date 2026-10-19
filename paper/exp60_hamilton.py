"""Bring the cluster results into the article.

    python exp60_hamilton.py hamilton_results.tar.gz

1. unpacks the archive written by hamilton/collect.sh into hamilton/results;
2. replays the held-out problem H3 (J5) with the frozen program and with the released
   version, as registered (exp62_heldout.py, both modes; H1 and H2 are replayed again
   alongside, unchanged);
3. prints the summaries of J1-J4 that the placeholders of the article refer to:
   J1  Section 7.1, forces function object of OpenFOAM v2406 against the article's
       post-processing, and the 2D-1 grid with a current release;
   J2  Section 7.2, the refined 2D-2 study;
   J3  Sections 7.4 and 9, Quantum ESPRESSO with a current release;
   J4  Supplement S11, parallel widths 1-64;
   J5  Section 6.2 and Table 6, H3.
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
