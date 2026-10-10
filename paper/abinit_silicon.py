"""A third production code: ABINIT, bulk silicon with a displaced atom.

Parameters are ABINIT's own, and map onto the taxonomy exactly as Table 5
predicts for a plane-wave DFT code:

    ecut      plane-wave cutoff             discretisation
    ngkpt     Brillouin-zone sampling       discretisation
    toldfe    SCF convergence threshold     nested solver
    nline     line minimisations per step   cost-only

Two observables are computed from the same runs: the total energy and the
magnitude of the force on the displaced atom.  This is the configuration in
which a Quantum ESPRESSO developer advised a user on the mailing list in 2009
that a threshold is "acceptable if you just want the energy but not if you want
forces".  If that is a general property rather than a quirk of one code, the
method should recover it here without being told.

Ground truth is a heavily converged reference, set strictly beyond the ladder in
every parameter.
"""
from __future__ import annotations

import os
import re
import resource
import shutil
import subprocess
import tempfile
import time

import numpy as np

TEMPLATE = """ecut       {ecut}
ngkpt      {nk} {nk} {nk}
nshiftk    1
shiftk     0.5 0.5 0.5
toldfe     {toldfe}
nline      {nline}
nstep      100
optforces  1
nsym       1
tolsym     1.0d-8
acell      3*10.26
rprim      0.0 0.5 0.5  0.5 0.0 0.5  0.5 0.5 0.0
ntypat     1
znucl      14
natom      2
typat      1 1
xred       0.00 0.00 0.00
           0.26 0.25 0.25
pp_dirpath "/usr/share/abinit/psp/"
pseudos    "14si.psp"
"""


class AbinitSilicon:
    def __init__(self):
        self.dir = tempfile.mkdtemp(prefix="knobs_abinit_")

    def __del__(self):
        shutil.rmtree(getattr(self, "dir", ""), ignore_errors=True)

    def run(self, theta: dict):
        inp = TEMPLATE.format(ecut=float(theta["ecut"]), nk=int(theta["ngkpt"]),
                              toldfe=float(theta["toldfe"]),
                              nline=int(theta["nline"]))
        path = os.path.join(self.dir, "run.abi")
        with open(path, "w") as f:
            f.write(inp)
        t0 = time.perf_counter()
        ru0 = resource.getrusage(resource.RUSAGE_CHILDREN)
        r = subprocess.run(["abinit", "run.abi"], cwd=self.dir,
                           capture_output=True)
        wall = time.perf_counter() - t0
        ru1 = resource.getrusage(resource.RUSAGE_CHILDREN)
        # CPU time of the ABINIT process: unlike wall time it does not depend
        # on what else the machine is running
        cpu = (ru1.ru_utime - ru0.ru_utime) + (ru1.ru_stime - ru0.ru_stime)
        out = r.stdout.decode(errors="ignore")
        if r.returncode != 0:
            return ({"energy": float("nan"), "force": float("nan")},
                    float("inf"), wall, {"unavailable": out[-200:]})
        e = re.search(r"etotal\s+(-?[\d.]+E?[+-]?\d*)", out)
        f = re.search(r"fcart\s+(-?[\d.]+E?[+-]?\d*)\s+(-?[\d.]+E?[+-]?\d*)"
                      r"\s+(-?[\d.]+E?[+-]?\d*)", out)
        if not e or not f:
            return ({"energy": float("nan"), "force": float("nan")},
                    float("inf"), wall, {"unavailable": "could not parse output"})
        fx, fy, fz = (float(f.group(i)) for i in (1, 2, 3))
        # ABINIT warns when the SCF cycle stops at nstep without meeting toldfe
        abo = [os.path.join(self.dir, n) for n in os.listdir(self.dir) if ".abo" in n]
        log = out + (open(max(abo, key=os.path.getmtime), errors="ignore").read() if abo else "")
        converged = "not enough SCF cycles to converge" not in log
        return ({"energy": float(e.group(1)),
                 "force": float(np.sqrt(fx * fx + fy * fy + fz * fz)),
                 "converged": converged},
                float(cpu), wall, {"converged": converged})

    def reference(self):
        return dict(ecut=48.0, ngkpt=8, toldfe=1e-14, nline=8)

    def runner(self, observable: str):
        def g(theta):
            obs, cost, _, _ = self.run(theta)
            return obs[observable], cost
        return g


if __name__ == "__main__":
    p = AbinitSilicon()
    print("computing the converged reference ...", flush=True)
    ref, c, w, info = p.run(p.reference())
    print(f"  energy = {ref['energy']:.10f} Ha   |force| = {ref['force']:.10e} "
          f"Ha/Bohr   ({w:.1f} s)")
    print()
    print("ecut  nk  toldfe    energy rel.err   force rel.err    s")
    for ecut in [6.0, 12.0, 20.0]:
        for tol in [1e-6, 1e-10]:
            obs, c, w, i = p.run(dict(ecut=ecut, ngkpt=4, toldfe=tol, nline=4))
            if i.get("unavailable"):
                print(f"{ecut:5.1f}   4  {tol:.0e}   UNAVAILABLE")
                continue
            de = abs(obs["energy"] - ref["energy"]) / abs(ref["energy"])
            df = abs(obs["force"] - ref["force"]) / abs(ref["force"])
            print(f"{ecut:5.1f}   4  {tol:.0e}   {de:.3e}        {df:.3e}   {w:5.1f}")
