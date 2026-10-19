"""Section 8.2: the arithmetic floor of the mixed- and double-precision GROMACS builds.

Six settings, each converged beyond the ladder in every parameter, are run with
each binary.  If the binary resolved the energy to better than a target, these
six values would agree to that accuracy; their spread is the floor below which
no parameter setting can certify the target.  The method is then run on the
mixed-precision tensor twice: without a floor probe (as in version 1.0) and
with the probe of version 1.1, which evaluates three of these settings, takes
their spread as the round-off floor and declines targets below ten times it.
"""
import sys, os, pickle, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from gromacs_water import GromacsWater
from knobs import run_protocol, measure_floor

CONVERGED = [dict(rcoulomb=rc, fourierspacing=fs, pme_order=po, ewald_rtol=rt)
             for rc, fs, po, rt in [(1.4, 0.06, 8, 1e-8), (1.4, 0.05, 8, 1e-9),
                                    (1.5, 0.06, 8, 1e-9), (1.6, 0.05, 8, 1e-9),
                                    (1.6, 0.06, 10, 1e-8), (1.5, 0.05, 10, 1e-9)]]
floor, probe = {}, {}
for b in ["gmx", "gmx_d"]:
    p = GromacsWater(nbox=2, nsteps=20, binary=b)
    v = np.array([p.run(th)[0]["potential"] for th in CONVERGED])
    floor[b] = (v.max() - v.min()) / abs(v.mean())
    probe[b] = measure_floor(v[:3])     # what the method itself measures: 3 runs
    print(f"{b:<6} floor probe (first three settings): {probe[b]:.4f} kJ/mol, "
          f"relative {probe[b]/abs(v.mean()):.2e}; declines targets below "
          f"{10*probe[b]/abs(v.mean()):.1e}")
    print(f"{b:<6} six converged settings: spread {v.max()-v.min():.4f} kJ/mol, "
          f"relative {floor[b]:.2e}   values {np.round(v, 4).tolist()}")

path = os.path.join(os.path.dirname(__file__), "data", "gromacs_tensor_gmx.pkl")
if os.path.exists(path):
    tensor, QS = pickle.load(open(path, "rb"))
    LAD = [[0.9, 1.0, 1.1, 1.2, 1.4], [0.20, 0.16, 0.12, 0.10, 0.08], [4, 5, 6], [1e-4, 1e-5, 1e-6]]
    HS = [[1/c for c in LAD[0]], list(LAD[1]), [1/o for o in LAD[2]], list(LAD[3])]
    KINDS = ["discretisation", "discretisation", "cost-only", "discretisation"]
    look = lambda i: (tensor[tuple(i)][0]["potential"], tensor[tuple(i)][1])
    ALL = [i for i in itertools.product(*[range(len(l)) for l in LAD]) if np.isfinite(tensor[i][1])]
    dpath = os.path.join(os.path.dirname(__file__), "data", "gromacs_tensor_gmx_d.pkl")
    QD = pickle.load(open(dpath, "rb"))[1] if os.path.exists(dpath) else None
    for use_floor in (False, True):
        print("without floor probe (v1.0)" if not use_floor else "with floor probe (v1.1)")
        for tg in [1e-4, 3e-5, 1e-5, 3e-6, 1e-6]:
            below = tg < floor["gmx"]
            cert = miss = dec = 0; worst = 0.0
            for a in ALL:
                r = run_protocol(look, LAD, a, tg * abs(QS), hs=HS, kinds=KINDS,
                                 use_screen=False,
                                 noise_floor=probe["gmx"] if use_floor else None)
                dec += r.get("declined") == "below_floor"
                if r["certified"]:
                    cert += 1
                    v = tensor[r["proposal"]][0]["potential"]
                    e = abs(v - (QD if QD is not None else QS)) / abs(QD if QD is not None else QS)
                    miss += e > tg; worst = max(worst, e)
            print(f"  mixed precision, target {tg:.0e} ({'below' if below else 'above'} the floor): "
                  f"certified from {cert}/{len(ALL)} starts, {miss} of them miss the target "
                  f"against the double-precision reference (worst error {worst:.2e}); "
                  f"{dec} declined at the floor check")
    # the probe must not get in the way where the arithmetic is adequate
    if QD is not None:
        td, QSd = pickle.load(open(dpath, "rb"))
        lookd = lambda i: (td[tuple(i)][0]["potential"], td[tuple(i)][1])
        for tg in [1e-5, 1e-6]:
            dec = sum(run_protocol(lookd, LAD, a, tg * abs(QSd), hs=HS, kinds=KINDS,
                                   noise_floor=probe["gmx_d"]).get("declined") == "below_floor"
                      for a in ALL)
            print(f"double precision, target {tg:.0e}: {dec}/{len(ALL)} starts declined at the floor check")
