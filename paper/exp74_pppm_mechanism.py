"""Revision 4 (referee point on the LAMMPS cost ordering): what kspace_modify order does
to the PPPM grid and to the cost, at settings of the study.

For each system size (n_cells 20, 40, 64: 8,000, 64,000 and 262,144 ions) and each of two
settings of the study, LAMMPS is run at PPPM orders 4-7, every run in its own process
(as lammps_subproc.py does, so that a grid that does not fit in memory is an ordinary
failed run).  From LAMMPS's own log it records:
  * the PPPM grid LAMMPS chose for the requested accuracy, and its number of points;
  * the G vector (Ewald splitting parameter) and the estimated relative force accuracy;
  * the "3d grid and FFT values/proc";
  * the time of one force evaluation, measured as in the study (run 0 post no, the
    median of n_eval repetitions after a first one that includes the setup);
  * from a 10-step run without integration (forces only), the split of the loop time
    between Pair, Kspace and Neigh in LAMMPS's timing breakdown;
  * the energy error against the Madelung energy.
The lattice, charges and pair style are those of lammps_madelung.py.

    python exp74_pppm_mechanism.py [n_cells ...]        # default 20 40 64
    ERRBUDGET_LMP_MEM_GB (default 6), ERRBUDGET_LMP_TIMEOUT (s, default 1800), ERRBUDGET_LMP_NEVAL (3)
Writes data/pppm_mechanism_v6.json.  About 30 minutes on the article's machine.
"""
import json
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lammps_madelung import MADELUNG, KE

SETTINGS = [dict(cutoff=10.0, kspace_accuracy=1e-5), dict(cutoff=8.0, kspace_accuracy=1e-6)]
ORDERS = [4, 5, 6, 7]
R0 = 2.82
MEM_GB = float(os.environ.get("ERRBUDGET_LMP_MEM_GB", 6))
TIMEOUT = float(os.environ.get("ERRBUDGET_LMP_TIMEOUT", 1800))
NEVAL = int(os.environ.get("ERRBUDGET_LMP_NEVAL", 3))

CHILD = r'''
import json, resource, sys, time
import numpy as np
resource.setrlimit(resource.RLIMIT_AS, ({mem}, {mem}))
import lammps
n, r0, cut, acc, order, neval, logf = {n}, {r0}, {cut}, {acc}, {order}, {neval}, {logf!r}
L = lammps.lammps(cmdargs=["-log", logf, "-screen", "none", "-nocite"])
c = L.command
c("units metal"); c("atom_style charge"); c("boundary p p p")
c(f"lattice sc {{r0}}"); c(f"region box block 0 {{n}} 0 {{n}} 0 {{n}}")
c("create_box 2 box"); c("create_atoms 1 box"); c("mass 1 22.99"); c("mass 2 35.45")
x = L.numpy.extract_atom("x")
idx = np.rint(x / r0).astype(int)
parity = (idx.sum(axis=1) % 2).astype(bool)
q = L.numpy.extract_atom("q"); q[:] = np.where(parity, -1.0, 1.0)
t = L.numpy.extract_atom("type"); t[:] = np.where(parity, 2, 1)
c(f"pair_style lj/cut/coul/long 1.0 {{cut}}"); c("pair_coeff * * 0.0 1.0")
c(f"kspace_style pppm {{acc}}"); c(f"kspace_modify order {{order}}")
c("neigh_modify every 1 delay 0 check no")
t0 = time.perf_counter(); c("run 0 post no"); setup = time.perf_counter() - t0
ts = []
for _ in range(neval):
    t1 = time.perf_counter(); c("run 0 post no"); ts.append(time.perf_counter() - t1)
e = float(L.get_thermo("ecoul") + L.get_thermo("elong"))
c("thermo 10"); c("run 10")
L.close()
print("RESULT" + json.dumps(dict(setup=setup, evals=ts, energy=e)))
'''


def parse_log(txt):
    out = {}
    pat = {
        "g_ewald": r"G vector \(1/distance\)\s*=\s*([-+\d.eE]+)",
        "stencil_order": r"stencil order\s*=\s*(\d+)",
        "acc_abs": r"estimated absolute RMS force accuracy\s*=\s*([-+\d.eE]+)",
        "acc_rel": r"estimated relative force accuracy\s*=\s*([-+\d.eE]+)",
    }
    for k, p in pat.items():
        m = re.findall(p, txt)
        if m:
            out[k] = float(m[-1])
    m = re.findall(r"^\s*grid\s*=\s*(\d+)\s+(\d+)\s+(\d+)", txt, flags=re.M)
    if m:
        out["grid"] = [int(v) for v in m[-1]]
        out["grid_points"] = out["grid"][0] * out["grid"][1] * out["grid"][2]
    m = re.findall(r"3d grid and FFT values/proc\s*=\s*(\d+)\s+(\d+)", txt)
    if m:
        out["grid_values_proc"], out["fft_values_proc"] = int(m[-1][0]), int(m[-1][1])
    m = re.findall(r"using (single|double) precision (\S+)", txt)
    if m:
        out["fft"] = " ".join(m[-1])
    m = re.findall(r"Loop time of ([-+\d.eE]+) on (\d+) procs for (\d+) steps", txt)
    if m:
        out["loop_time_10"] = float(m[-1][0])
        out["procs"] = int(m[-1][1])
    sec = {}
    for name, avg in re.findall(r"^(Pair|Bond|Kspace|Neigh|Comm|Output|Modify|Other)\s*\|"
                                r"\s*[-+\d.eE]*\s*\|\s*([-+\d.eE]+)", txt, flags=re.M):
        sec[name] = float(avg)
    if sec:
        out["breakdown_10"] = sec
    return out


def run_one(n, theta, order):
    fd, logf = tempfile.mkstemp(prefix="errbudget_pppm_", suffix=".log")
    os.close(fd)
    code = CHILD.format(mem=int(MEM_GB * 1e9), n=n, r0=R0, cut=theta["cutoff"],
                        acc=theta["kspace_accuracy"], order=order, neval=NEVAL, logf=logf)
    rec = dict(n_cells=n, ions=n ** 3, order=order, **theta)
    try:
        r = subprocess.run([sys.executable, "-c", code], capture_output=True, timeout=TIMEOUT)
        out = r.stdout.decode(errors="ignore")
        if "RESULT" in out:
            res = json.loads(out.split("RESULT", 1)[1].split("\n")[0])
            evals = sorted(res["evals"])
            rec.update(setup_s=res["setup"], per_eval_s=evals[len(evals) // 2], evals_s=res["evals"],
                       energy=res["energy"])
            ex = -(n ** 3 / 2.0) * MADELUNG * KE / R0
            rec["rel_error"] = abs(res["energy"] - ex) / abs(ex)
        else:
            rec["unavailable"] = f"exit {r.returncode}: {r.stderr.decode(errors='ignore')[-300:]}"
    except subprocess.TimeoutExpired:
        rec["unavailable"] = "timeout"
    if os.path.exists(logf):
        rec.update(parse_log(open(logf, errors="ignore").read()))
        os.remove(logf)
    return rec


def main():
    sizes = [int(x) for x in sys.argv[1:]] or [20, 40, 64]
    rows = []
    print(f"{'ions':>7} {'cut':>4} {'acc':>6} {'ord':>3} {'grid':>12} {'points':>9} {'G':>7} "
          f"{'est.acc':>8} {'t/eval s':>9} {'Pair':>7} {'Kspace':>7} {'Neigh':>7} {'rel.err':>8}")
    for n in sizes:
        for th in SETTINGS:
            for order in ORDERS:
                r = run_one(n, th, order)
                rows.append(r)
                b = r.get("breakdown_10", {})
                g = "x".join(map(str, r.get("grid", []))) or "-"
                print(f"{n**3:>7} {th['cutoff']:>4g} {th['kspace_accuracy']:>6.0e} {order:>3} {g:>12} "
                      f"{r.get('grid_points', 0):>9} {r.get('g_ewald', float('nan')):7.4f} "
                      f"{r.get('acc_rel', float('nan')):8.1e} {r.get('per_eval_s', float('nan')):9.4f} "
                      f"{b.get('Pair', float('nan')):7.3f} {b.get('Kspace', float('nan')):7.3f} "
                      f"{b.get('Neigh', float('nan')):7.3f} {r.get('rel_error', float('nan')):8.1e}"
                      + (f"  UNAVAILABLE: {r['unavailable'][:60]}" if "unavailable" in r else ""),
                      flush=True)
    path = os.path.join(HERE, "data", "pppm_mechanism_v6.json")
    json.dump(rows, open(path, "w"), indent=1)
    print("written", path)
    # the cheapest order per (size, setting), and how its grid compares with order 5
    for n in sizes:
        for th in SETTINGS:
            R = [r for r in rows if r["n_cells"] == n and r["cutoff"] == th["cutoff"]
                 and "per_eval_s" in r]
            if not R:
                continue
            best = min(R, key=lambda r: r["per_eval_s"])
            o5 = [r for r in R if r["order"] == 5]
            msg = ""
            if o5 and "grid_points" in best and "grid_points" in o5[0]:
                msg = (f"; grid points {best['grid_points']} against {o5[0]['grid_points']} at order 5, "
                       f"time {best['per_eval_s'] / o5[0]['per_eval_s']:.2f} of order 5")
            print(f"  {n**3} ions, cutoff {th['cutoff']:g}, accuracy {th['kspace_accuracy']:.0e}: "
                  f"cheapest order {best['order']}{msg}")


if __name__ == "__main__":
    main()
