"""Run one LAMMPS configuration in a separate process.

At production system sizes the binding constraint is memory, not time: a tight
PPPM accuracy at a short cutoff demands a grid that does not fit, and an
out-of-memory kill takes down the whole interpreter rather than raising. Running
each configuration in its own process turns that into an ordinary unavailable
result, which is the failure path the method already handles.
"""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def run_one(n_cells, theta, n_eval=3, timeout=600, mem_gb=5):
    code = (
        "import sys, json, resource\n"
        f"resource.setrlimit(resource.RLIMIT_AS, ({int(mem_gb*1e9)}, {int(mem_gb*1e9)}))\n"
        f"sys.path.insert(0, {HERE!r})\n"
        "from lammps_madelung import LammpsMadelung\n"
        f"p = LammpsMadelung(n_cells={n_cells}, n_eval={n_eval})\n"
        f"obs, cost, wall, info = p.run({theta!r})\n"
        "print('RESULT'+json.dumps([obs, cost, wall]))\n")
    try:
        r = subprocess.run([sys.executable, "-c", code], capture_output=True,
                           timeout=timeout)
    except subprocess.TimeoutExpired:
        return ({"coulomb_energy": float("nan"),
                 "energy_per_ion": float("nan")}, float("inf"),
                {"unavailable": "timeout"})
    out = r.stdout.decode(errors="ignore")
    if "RESULT" in out:
        obs, cost, wall = json.loads(out.split("RESULT", 1)[1].split("\n")[0])
        return obs, cost, {}
    why = ("out of memory" if r.returncode in (-9, 137, 1)
           and ("MemoryError" in r.stderr.decode(errors="ignore")
                or r.returncode in (-9, 137))
           else f"exit {r.returncode}")
    return ({"coulomb_energy": float("nan"), "energy_per_ion": float("nan")},
            float("inf"), {"unavailable": why})
