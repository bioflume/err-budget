"""A second third-party solver: GROMACS electrostatics on a water box.

The observable is the potential energy of one fixed configuration of 1728 SPC
water molecules (5184 atoms) -- the value at step 0, which is a deterministic
function of the numerical parameters alone.  Cost is the per-step wall time of a
short run from that configuration.  The parameters are GROMACS's own documented
.mdp options:

    rcoulomb            real-space electrostatic cutoff   discretisation
    fourierspacing      PME grid spacing                  discretisation
    ewald-rtol          real/reciprocal splitting         discretisation
    pme-order           PME interpolation order           cost-only

GROMACS ships defaults for three of these (fourierspacing 0.12, pme-order 4,
ewald-rtol 1e-5) and its documentation describes `verlet-buffer-tolerance` as
setting "the maximum allowed error for pair interactions per particle", i.e. an
accuracy-targeted mechanism of the same kind this paper automates.  The cutoff
is left to the user.

Ground truth is a heavily converged reference rather than an analytic result,
which is the usual situation for a production code.  The reference is set
strictly beyond the ladder in every parameter, so that no tested configuration
inherits an artificially small error by sharing a value with it.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import time

TOP = """#include "oplsaa.ff/forcefield.itp"
#include "oplsaa.ff/spc.itp"

[ system ]
water box

[ molecules ]
SOL {nmol}
"""

MDP = """integrator              = md
nsteps                  = {nsteps}
dt                      = 0.002
cutoff-scheme           = Verlet
coulombtype             = PME
vdwtype                 = cutoff
rvdw                    = 0.9
rlist                   = {rlist}
rcoulomb                = {rcoulomb}
fourierspacing          = {fourierspacing}
pme-order               = {pme_order}
ewald-rtol              = {ewald_rtol}
verlet-buffer-tolerance = -1
nstlist                 = 10
freezegrps              = System
freezedim               = Y Y Y
nstcalcenergy           = 1
nstenergy               = 1
continuation            = yes
constraints             = none
"""


class GromacsWater:
    def __init__(self, nbox: int = 2, nsteps: int = 40, threads: int = 2,
                 binary: str = "gmx", ntmpi: int = 1):
        self.nsteps, self.threads, self.gmx = nsteps, threads, binary
        self.ntmpi = ntmpi              # thread-MPI ranks (parallel timing study)
        self.dir = tempfile.mkdtemp(prefix="errbudget_gmx_")
        self.nmol = 216 * nbox ** 3
        subprocess.run(
            [binary, "genconf", "-f", os.path.join(os.environ.get("GMXDATA", "/usr/share/gromacs"), "top", "spc216.gro"),
             "-nbox", str(nbox), str(nbox), str(nbox), "-o", "conf.gro"],
            cwd=self.dir, capture_output=True, check=True)
        with open(os.path.join(self.dir, "topol.top"), "w") as f:
            f.write(TOP.format(nmol=self.nmol))

    def __del__(self):
        shutil.rmtree(getattr(self, "dir", ""), ignore_errors=True)

    def run(self, theta: dict):
        rc = float(theta["rcoulomb"])
        fs = float(theta["fourierspacing"])
        order = int(theta["pme_order"])
        rtol = float(theta["ewald_rtol"])
        d = self.dir

        with open(os.path.join(d, "md.mdp"), "w") as f:
            f.write(MDP.format(nsteps=self.nsteps, rcoulomb=rc, rlist=rc,
                               fourierspacing=fs, pme_order=order,
                               ewald_rtol=rtol))
        # GROMACS renames existing output to #file.N# and refuses past 99
        # backups, which silently kills a long sweep partway through.
        env = dict(os.environ, GMX_MAXBACKUP="-1")
        g = subprocess.run(
            [self.gmx, "grompp", "-f", "md.mdp", "-c", "conf.gro", "-p", "topol.top",
             "-o", "t.tpr", "-maxwarn", "10"], cwd=d, capture_output=True, env=env)
        if g.returncode != 0:
            return ({"potential": float("nan")}, float("inf"), 0.0,
                    {"unavailable": g.stderr.decode()[-200:]})

        t0 = time.perf_counter()
        m = subprocess.run(
            [self.gmx, "mdrun", "-s", "t.tpr", "-deffnm", "run", "-ntmpi", str(self.ntmpi),
             "-ntomp", str(self.threads), "-nb", "cpu", "-pme", "cpu",
             "-notunepme"], cwd=d, capture_output=True, env=env)
        wall = time.perf_counter() - t0
        if m.returncode != 0:
            return ({"potential": float("nan")}, float("inf"), wall,
                    {"unavailable": m.stderr.decode()[-200:]})

        # The summary line gmx energy prints carries too few significant
        # figures to resolve the differences being measured, so the per-frame
        # values are read from the data file instead.
        subprocess.run([self.gmx, "energy", "-f", "run.edr", "-o", "e.xvg",
                        "-xvg", "none"], cwd=d, input=b"Potential\n",
                       capture_output=True, env=env)
        vals = []
        try:
            with open(os.path.join(d, "e.xvg")) as fh:
                for line in fh:
                    if line[0] not in "#@":
                        vals.append(float(line.split()[1]))
        except OSError:
            vals = []
        val = vals[0] if vals else None
        spread = (max(vals) - min(vals)) if vals else 0.0
        if val is None:
            return ({"potential": float("nan")}, float("inf"), wall,
                    {"unavailable": "energy term not found"})
        # cost is the per-force-evaluation wall time
        return ({"potential": val}, float(wall / max(self.nsteps, 1)), wall,
                {"frames": len(vals), "frame_spread": spread})

    def reference(self):
        """A heavily converged setting, used as ground truth."""
        return dict(rcoulomb=1.6, fourierspacing=0.05, pme_order=8,
                    ewald_rtol=1e-9)


if __name__ == "__main__":
    p = GromacsWater(nbox=2, nsteps=40)
    print(f"{p.nmol} SPC waters ({p.nmol*3} atoms), all atoms frozen")
    ref, cost, wall, info = p.run(p.reference())
    print(f"converged reference: potential = {ref['potential']:.6f} kJ/mol "
          f"({cost*1000:.1f} ms/step)")
    print()
    print("rcoulomb  fourier  order  rtol     potential        rel.err     ms/eval")
    for rc in [0.9, 1.1, 1.4]:
        for fs in [0.16, 0.12, 0.08]:
            obs, c, w, i = p.run(dict(rcoulomb=rc, fourierspacing=fs,
                                      pme_order=4, ewald_rtol=1e-5))
            if i.get("unavailable"):
                print(f"{rc:8.1f} {fs:8.2f}  4     1e-05   UNAVAILABLE")
                continue
            rel = abs(obs["potential"] - ref["potential"]) / abs(ref["potential"])
            print(f"{rc:8.1f} {fs:8.2f}  4     1e-05  {obs['potential']:14.6f}"
                  f"   {rel:.3e}   {c*1000:7.1f}")
