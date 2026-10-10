"""Did any ABINIT run of the tensor stop at nstep without meeting toldfe?

The tensor of Section 7 was built before the driver reported SCF convergence.  The
settings most likely to stop at the limit (nstep = 100) are those with the tightest
toldfe and the fewest line minimisations; they are re-run here with the driver's
convergence flag, and their energies are compared with the stored ones.
    python exp68_abinit_scf.py
"""
import os, pickle, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from abinit_silicon import AbinitSilicon
import cases as C

case = C._abinit()
LAD = case.LAD
p = AbinitSilicon()
bad = 0
for a in range(len(LAD[0])):
    for b in range(len(LAD[1])):
        idx = (a, b, len(LAD[2]) - 1, 0)
        th = dict(ecut=LAD[0][a], ngkpt=LAD[1][b], toldfe=LAD[2][-1], nline=LAD[3][0])
        obs, cost, wall, info = p.run(th)
        stored = case.tensor[idx][0]["energy"]
        bad += not info.get("converged", True)
        print(th, "converged" if info.get("converged", True) else "NOT CONVERGED",
              f"energy {obs['energy']:.10f} stored {stored:.10f}", flush=True)
print(f"{bad} of {len(LAD[0]) * len(LAD[1])} runs at the tightest toldfe and nline = {LAD[3][0]} "
      "stopped at nstep")
