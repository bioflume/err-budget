"""J4: does communication at 1-64 ranks couple the parameters or change the
recommendation? (referee point 11)

  lmp_R<R>   LAMMPS, 262,144 ions, the 64-setting grid of Section 8, on R MPI ranks
  gmx_R<R>   GROMACS in double precision, 13,824 SPC waters (nbox 4), the
             225-setting grid of Section 8, on R thread-MPI ranks
  gmx_ref    the converged GROMACS reference for that box (beyond every value)
R = 1, 2, 4, 8, 16, 32, 64.  Each task runs its whole grid on one whole node,
so that the timings of one width are taken under the same conditions.
    python -m farm parallel_hpc array K NK
"""
import itertools
import json
import os
import pickle
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper"))
import exp55_parallel as E

RANKS = [1, 2, 4, 8, 16, 32, 64]
NBOX = 4
SMOKE = os.environ.get("ERRBUDGET_SMOKE") == "1"
if SMOKE:
    RANKS, NBOX = [1, 2], 2
NG = ["rcoulomb", "fourierspacing", "pme_order", "ewald_rtol"]
LADG = [[0.9, 1.0, 1.1, 1.2, 1.4], [0.20, 0.16, 0.12, 0.10, 0.08], [4, 5, 6], [1e-4, 1e-5, 1e-6]]
GMX = os.path.join(os.environ.get("GMX_DOUBLE_PREFIX", ""), "bin", "gmx_d")


def TASKS():
    T = [dict(id="gmx_ref", kind="gmx_ref", R=8, nprocs=8)]
    T += [dict(id=f"lmp_R{R}", kind="lmp", R=R, nprocs=R) for R in RANKS]
    T += [dict(id=f"gmx_R{R}", kind="gmx", R=R, nprocs=R) for R in RANKS]
    return T


def run(task, out):
    if task["kind"] == "lmp":
        t = {}
        grid = list(itertools.product(*[range(len(l)) for l in E.LADL]))
        if SMOKE:
            grid = grid[:3]
        for idx in grid:
            th = tuple(E.LADL[k][j] for k, j in enumerate(idx))
            e, c = E.lmp_run(th, task["R"])
            t[idx] = ({"energy_per_ion": e}, c)
        if not any(np.isfinite(c) for _, c in t.values()):
            # no setting ran at all: an infrastructure failure, not a result (error.txt, retried)
            raise RuntimeError(f"no LAMMPS run succeeded on {task['R']} ranks; see this task's log")
        pickle.dump(t, open(os.path.join(out, "tensor.pkl"), "wb"))
    else:
        os.environ["GMXDATA"] = os.path.join(os.environ.get("GMX_DOUBLE_PREFIX", ""), "share", "gromacs")
        from gromacs_water import GromacsWater
        p = GromacsWater(nbox=NBOX, nsteps=20, binary=GMX, threads=1, ntmpi=task["R"])
        if task["kind"] == "gmx_ref":
            vals = []
            for th in [p.reference(), dict(p.reference(), fourierspacing=0.06, pme_order=10)]:
                obs, c, w, info = p.run(th)
                if not np.isfinite(obs["potential"]):
                    raise RuntimeError(f"GROMACS reference run failed: {info}")
                vals.append(obs["potential"])
            json.dump(dict(QS=vals[0], spread=abs(vals[1] - vals[0])),
                      open(os.path.join(out, "ref.json"), "w"))
        else:
            t = {}
            grid = list(itertools.product(*[range(len(l)) for l in LADG]))
            if SMOKE:
                grid = grid[:3]
            for idx in grid:
                obs, c, w, info = p.run({NG[k]: LADG[k][j] for k, j in enumerate(idx)})
                t[idx] = (obs, c)
            if not any(np.isfinite(c) for _, c in t.values()):
                raise RuntimeError(f"no GROMACS run succeeded on {task['R']} ranks; last: {info}")
            pickle.dump(t, open(os.path.join(out, "tensor.pkl"), "wb"))
    json.dump(dict(task), open(os.path.join(out, "result.json"), "w"))
