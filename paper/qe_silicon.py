"""Quantum ESPRESSO (pw.x): bulk silicon with one displaced atom.

The same system as the ABINIT study (Section 8.3): diamond Si, a = 10.26 bohr,
the second atom displaced from (1/4, 1/4, 1/4) to (0.26, 0.25, 0.25), symmetry
off, PseudoDojo norm-conserving PBE pseudopotential (the one shipped with
ABINIT).  Parameters are pw.x's own:

    ecutwfc      plane-wave cutoff (Ry)                discretisation
    nk           k-point mesh nk x nk x nk, shifted     discretisation
    conv_thr     SCF convergence threshold (Ry)         nested solver
    mixing_beta  density mixing factor                  cost-only

Observables: the total energy (Ry) and the magnitude of the force on the
displaced atom (Ry/bohr).  A QE developer's advice on the mailing list, quoted
in the introduction, is that a conv_thr adequate for the energy may not be
adequate for forces; this study measures it.

Ubuntu 24.04's pw.x 6.7 aborts in a fortified snprintf when it hashes the
pseudopotential file; snfix.c, preloaded, restores the unfortified behaviour.
Cost is the CPU time of the pw.x process.
"""
import os
import re
import resource
import shutil
import subprocess
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PSEUDO_DIR = os.environ.get("KNOBS_QE_PSEUDO_DIR",
                            "/usr/share/abinit/psp/Pseudodojo_nc_sr_04_pbe_standard_psp8/")
# The shim is needed only for the Ubuntu 24.04 package of pw.x 6.7.
USE_SHIM = os.environ.get("KNOBS_QE_SHIM", "1") == "1"
PW = os.environ.get("KNOBS_PW", "pw.x")
TEMPLATE = """&control
  calculation='scf', prefix='si', pseudo_dir='{pd}', outdir='./tmp', tprnfor=.true.
/
&system
  ibrav=2, celldm(1)=10.26, nat=2, ntyp=1, ecutwfc={ecut}, nosym=.true.
/
&electrons
  conv_thr={thr:.1e}, mixing_beta={beta}, electron_maxstep=200
/
ATOMIC_SPECIES
Si 28.086 Si.upf
ATOMIC_POSITIONS crystal
Si 0.00 0.00 0.00
Si 0.26 0.25 0.25
K_POINTS automatic
{nk} {nk} {nk} 1 1 1
"""


def _shim():
    so = os.path.join(tempfile.gettempdir(), "knobs_snfix.so")
    if not os.path.exists(so):
        subprocess.run(["gcc", "-shared", "-fPIC", "-O2", "-o", so,
                        os.path.join(HERE, "snfix.c")], check=True)
    return so


class QESilicon:
    def __init__(self):
        self.dir = tempfile.mkdtemp(prefix="knobs_qe_")
        self.env = dict(os.environ, OMP_NUM_THREADS="1")
        if USE_SHIM:
            self.env["LD_PRELOAD"] = _shim()

    def __del__(self):
        shutil.rmtree(getattr(self, "dir", ""), ignore_errors=True)

    def reference(self):
        return dict(ecutwfc=90.0, nk=8, conv_thr=1e-14, mixing_beta=0.7)

    def run(self, theta):
        inp = TEMPLATE.format(pd=PSEUDO_DIR, ecut=float(theta["ecutwfc"]),
                              nk=int(theta["nk"]), thr=float(theta["conv_thr"]),
                              beta=float(theta["mixing_beta"]))
        shutil.rmtree(os.path.join(self.dir, "tmp"), ignore_errors=True)
        with open(os.path.join(self.dir, "pw.in"), "w") as f:
            f.write(inp)
        ru0 = resource.getrusage(resource.RUSAGE_CHILDREN)
        t0 = time.perf_counter()
        r = subprocess.run([PW, "-in", "pw.in"], cwd=self.dir, env=self.env,
                           capture_output=True, text=True)
        wall = time.perf_counter() - t0
        ru1 = resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu = (ru1.ru_utime - ru0.ru_utime) + (ru1.ru_stime - ru0.ru_stime)
        out = r.stdout
        nan = {"energy": float("nan"), "force": float("nan")}
        e = re.search(r"^!\s+total energy\s+=\s+(-?[\d.]+) Ry", out, flags=re.M)
        f = re.search(r"atom\s+2 type\s+1\s+force =\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)", out)
        if r.returncode != 0 or not e or not f or "convergence NOT achieved" in out:
            return nan, float("inf"), wall, {"unavailable": out[-300:]}
        F = np.array([float(f.group(i)) for i in (1, 2, 3)])
        nit = re.search(r"convergence has been achieved in\s+(\d+) iterations", out)
        return ({"energy": float(e.group(1)), "force": float(np.linalg.norm(F))},
                float(cpu), wall, dict(scf_iters=int(nit.group(1)) if nit else -1))
