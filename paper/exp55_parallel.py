"""Item 4: does parallel execution change the cost ordering or the recommendation?

LAMMPS at 262,144 ions (the 64-setting grid of Section 8.1) is re-run with the
`lmp` executable under MPI on R ranks, and GROMACS (gmx_d, 1,728 waters) with T
OpenMP threads.  For each R or T: the correlation of the per-setting costs with
the serial (or 1-thread) costs, the cost ratio of each parameter's steps, and
the method's recommendation from every start.  On the 2-core machine used for
the article only R, T = 1, 2 are possible; set ERRBUDGET_RANKS="1 2 8 16" and
ERRBUDGET_THREADS="1 2 8" to repeat it on a cluster node.
"""
import itertools
import os
import pickle
import re
import subprocess
import sys
import tempfile
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from errbudget import run_protocol
from lammps_madelung import LammpsMadelung, MADELUNG, KE

D = os.path.join(os.path.dirname(__file__), "data")
RANKS = [int(r) for r in os.environ.get("ERRBUDGET_RANKS", "1 2").split()]
THREADS = [int(t) for t in os.environ.get("ERRBUDGET_THREADS", "1 2").split()]

# ---------------------------------------------------------------- LAMMPS
N = 64
LADL = [[8.0, 10.0, 12.0, 14.0], [1e-3, 1e-4, 1e-5, 1e-6], [4, 5, 6, 7]]
HSL = [[1 / c for c in LADL[0]], list(LADL[1]), [1 / o for o in LADL[2]]]
KL = ["discretisation", "discretisation", "cost-only"]
R0 = 2.82
QS = -0.5 * MADELUNG * KE / R0
INPUT = """units metal
atom_style charge
boundary p p p
lattice sc {r0}
region box block 0 {n} 0 {n} 0 {n}
create_box 2 box
create_atoms 1 box
mass * 22.99
variable par atom (round(x/{r0})+round(y/{r0})+round(z/{r0}))-2*floor((round(x/{r0})+round(y/{r0})+round(z/{r0}))/2)
group odd variable par
set group all charge 1.0
set group odd charge -1.0
set group odd type 2
pair_style lj/cut/coul/long 1.0 {cut}
pair_coeff * * 0.0 1.0
kspace_style pppm {acc}
kspace_modify order {order}
neigh_modify every 1 delay 0 check no
thermo_style custom step ecoul elong
thermo_modify format float %.15g
run 0 post no
run {neval} pre no post no
"""


LAST_ERROR = ""          # why the last failed lmp_run failed (for the cluster job J4)


def lmp_run(th, ranks, neval=3, mem_gb=6):
    # On a cluster the memory limit and the MPI launcher may need changing (hamilton/config.sh):
    # ERRBUDGET_LMP_MEM_GB (0: no limit) and ERRBUDGET_LMP_LAUNCH (with {R} for the ranks).
    global LAST_ERROR
    mem_gb = float(os.environ.get("ERRBUDGET_LMP_MEM_GB", mem_gb))
    launch = os.environ.get("ERRBUDGET_LMP_LAUNCH", "mpirun --allow-run-as-root --oversubscribe -np {R}")
    d = tempfile.mkdtemp(prefix="errbudget_lmp_")
    open(os.path.join(d, "in.lmp"), "w").write(INPUT.format(
        r0=R0, n=N, cut=th[0], acc=th[1], order=th[2], neval=neval))
    cmd = launch.format(R=ranks).split() + [
        os.environ.get("ERRBUDGET_LMP", "lmp"), "-in", "in.lmp", "-log", "log.lammps", "-nocite"]
    limit = f"ulimit -v {int(mem_gb*1e6)}; " if mem_gb > 0 else ""
    try:
        r = subprocess.run(["bash", "-c", limit + " ".join(cmd)],
                           cwd=d, capture_output=True, text=True, timeout=900)
    except subprocess.TimeoutExpired:
        LAST_ERROR = f"timed out after 900 s: {limit}{' '.join(cmd)}"
        return np.nan, np.inf
    log = open(os.path.join(d, "log.lammps")).read() if os.path.exists(os.path.join(d, "log.lammps")) else ""
    m = re.findall(r"^\s+0\s+(-?[\d.eE+-]+)\s+(-?[\d.eE+-]+)\s*$", log, flags=re.M)
    lt = re.findall(r"Loop time of ([\d.eE+-]+) on (\d+) procs for (\d+) steps", log)
    if lt and int(lt[-1][1]) != ranks:
        # an lmp built for another MPI (or without MPI) runs once per rank, each copy
        # alone: the energy is right but the time is that of one rank
        LAST_ERROR = (f"LAMMPS ran on {lt[-1][1]} process(es) instead of {ranks}: the mpirun "
                      f"({launch}) does not match the MPI of {os.environ.get('ERRBUDGET_LMP', 'lmp')}")
        sys.stderr.write(f"lmp_run {th} on {ranks} ranks: {LAST_ERROR}\n")
        return np.nan, np.inf
    lt = [(a, c) for a, b, c in lt]
    if r.returncode != 0 or not m or len(lt) < 2:
        LAST_ERROR = (f"exit {r.returncode}: {limit}{' '.join(cmd)}\nstdout: {r.stdout[-800:]}\n"
                      f"stderr: {r.stderr[-1200:]}\nlog.lammps: {log[-800:]}")
        sys.stderr.write(f"lmp_run {th} on {ranks} ranks: {LAST_ERROR}\n")
        return np.nan, np.inf
    e = (float(m[0][0]) + float(m[0][1])) / N**3
    return e, float(lt[-1][0]) / int(lt[-1][1])


def compare(name, tensors, LAD, HS, KINDS, look_ob, QS, targets):
    keys = sorted(tensors)
    base = tensors[keys[0]]
    ALL = [i for i in itertools.product(*[range(len(l)) for l in LAD])
           if all(np.isfinite(tensors[k][i][1]) for k in keys)]
    print(f"\n{name}: {len(ALL)} settings usable at every parallel width")
    for k in keys:
        t = tensors[k]
        c0 = np.log([base[i][1] for i in ALL]); c1 = np.log([t[i][1] for i in ALL])
        rk = lambda v: np.argsort(np.argsort(v))
        rho = np.corrcoef(rk(c0), rk(c1))[0, 1]
        steps = []
        for ax in range(len(LAD)):
            r = [t[i[:ax] + (i[ax] + 1,) + i[ax + 1:]][1] / t[i][1] for i in ALL
                 if i[ax] + 1 < len(LAD[ax]) and i[:ax] + (i[ax] + 1,) + i[ax + 1:] in ALL]
            steps.append(np.median(r))
        print(f"  width {k:>2}: median cost {np.median([t[i][1] for i in ALL])*1e3:8.1f} ms, "
              f"rank correlation with width {keys[0]}: {rho:.3f}; median cost ratio of one "
              f"step finer: " + ", ".join(f"{s:.2f}" for s in steps))
    for tg in targets:
        eps = tg * abs(QS)
        recs = {}
        for k in keys:
            t = tensors[k]
            lk = lambda i, t=t: (look_ob(t[tuple(i)][0]), t[tuple(i)][1])
            out = {}
            for a in ALL:
                r = run_protocol(lk, LAD, a, eps, hs=HS, kinds=KINDS)
                out[a] = (r["proposal"], r["certified"])
            recs[k] = out
            feas = [i for i in ALL if abs(look_ob(t[i][0]) - QS) <= eps]
            opt = min(feas, key=lambda i: t[i][1]) if feas else None
            pen = [t[q][1] / t[opt][1] for q, c in out.values() if c and opt]
            print(f"  target {tg:.0e}, width {k:>2}: certified {sum(c for _, c in out.values())}"
                  f"/{len(ALL)}, cheapest attaining {opt}, median penalty "
                  f"{np.median(pen) if pen else float('nan'):.2f}")
        for k in keys[1:]:
            same = sum(recs[k][a] == recs[keys[0]][a] for a in ALL)
            print(f"    same recommendation at width {k} as at width {keys[0]}: {same}/{len(ALL)}")


if __name__ == "__main__":
    tl = {}
    for R in RANKS:
        path = os.path.join(D, f"lammps_mpi_n{N}_r{R}.pkl")
        if os.path.exists(path):
            tl[R] = pickle.load(open(path, "rb")); continue
        t, t0 = {}, time.time()
        for idx in itertools.product(*[range(len(l)) for l in LADL]):
            th = tuple(LADL[k][j] for k, j in enumerate(idx))
            e, c = lmp_run(th, R)
            t[idx] = ({"energy_per_ion": e}, c)
        pickle.dump(t, open(path, "wb")); tl[R] = t
        print(f"LAMMPS {R} ranks: {len(t)} runs in {time.time()-t0:.0f}s", flush=True)
    ser = pickle.load(open(os.path.join(D, f"lammps_tensor_n{N}.pkl"), "rb"))
    dev = max(abs(tl[RANKS[0]][i][0]["energy_per_ion"] - ser[i][0]["energy_per_ion"])
              for i in ser if np.isfinite(ser[i][1]) and np.isfinite(tl[RANKS[0]][i][1]))
    print(f"largest difference from the in-process tensor of Section 8.1: {dev:.2e} eV/ion")
    compare("LAMMPS, 262,144 ions, MPI ranks", tl, LADL, HSL, KL,
            lambda o: o["energy_per_ion"], QS, [1e-5, 1e-6])

    from gromacs_water import GromacsWater
    NG = ["rcoulomb", "fourierspacing", "pme_order", "ewald_rtol"]
    LADG = [[0.9, 1.0, 1.1, 1.2, 1.4], [0.20, 0.16, 0.12, 0.10, 0.08], [4, 5, 6],
            [1e-4, 1e-5, 1e-6]]
    KG = ["discretisation", "discretisation", "cost-only", "discretisation"]
    HSG = [[1.0 / c for c in LADG[0]], list(LADG[1]), [1.0 / o for o in LADG[2]],
           list(LADG[3])]
    tg, QG = {}, None
    for T in THREADS:
        path = os.path.join(D, f"gromacs_gmx_d_threads{T}.pkl")
        if os.path.exists(path):
            tg[T], QG = pickle.load(open(path, "rb")); continue
        p = GromacsWater(nbox=2, nsteps=20, binary="gmx_d", threads=T)
        QG = pickle.load(open(os.path.join(D, "gromacs_tensor_gmx_d.pkl"), "rb"))[1]
        t = {}
        for idx in itertools.product(*[range(len(l)) for l in LADG]):
            obs, c, w, info = p.run({NG[k]: LADG[k][j] for k, j in enumerate(idx)})
            t[idx] = (obs, c)
        pickle.dump((t, QG), open(path, "wb")); tg[T] = t
        print(f"GROMACS gmx_d {T} threads: {len(t)} runs", flush=True)
    compare("GROMACS gmx_d, OpenMP threads", tg, LADG, HSG, KG,
            lambda o: o["potential"], QG, [1e-5, 1e-6])
