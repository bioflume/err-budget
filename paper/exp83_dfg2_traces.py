"""OpenFOAM 2D-2 (refined study, case 2D-2h): the reference remark of S-2.5 and the misses
at the two ungraded targets of S-5.2, from hamilton/results/j2_tensor.pkl.

1. Observed orders of the three observables through n_arc = 64, 96 and 128 (dt = 2.5e-4,
   tightest tolerances, R = 2, W = 4), the limits that Richardson extrapolation puts beyond
   the reference at these orders and at second order, and the lift at 3e-2 graded against
   them.
2. The released program at the ungraded targets (drag 3e-3, lift 1e-2): accepted searches,
   misses, the settings returned.
3. The traced search for the drag from n_arc = 24, dt = 2e-3, residualControl 1e-2, GAMG
   1e-6, R = W = 1: the budget at the returned setting and at the other missed setting; the
   changes of the drag and the lift along the mesh at the two time steps at which every
   finer mesh runs.

    python exp83_dfg2_traces.py        # about a minute; writes data/dfg2_traces.json
"""
import json
import os

import numpy as np

import cases as C
import pooled as P
from errbudget import run_protocol
from paths import OUT

CASE = C._dfg2_refined()
CFG = P.CONFIGS["v1.2.1"]
T = CASE.tensor
REFS = C.pickle.load(open(os.path.join(C.D, "..", "..", "hamilton", "results", "j2_tensor.pkl"),
                          "rb"))["refs"]
OBS = ("cd_max", "cl_max", "st")
RES = {}


def look_for(ob):
    def look(i):
        v = T[tuple(i)]
        if not np.isfinite(v[1]):
            return np.nan, np.inf
        if not CASE.converged(tuple(i)):
            return np.nan, v[1]
        return v[0][ob], v[1]
    return look


def observed_order(f1, f2, f3, r21, r32):
    """Celik et al. (2008), Eq. (3): order from three meshes with ratios r21 = h2/h1 and
    r32 = h3/h2 (1 finest), by fixed-point iteration."""
    e21, e32 = f2 - f1, f3 - f2
    s = np.sign(e32 / e21)
    p = 2.0
    for _ in range(200):
        p = abs(np.log(abs(e32 / e21)) + np.log((r21**p - s) / (r32**p - s))) / np.log(r21)
    return p


# 1. Limits beyond the reference (S-2.5, remark)
fine = {64: T[(3, 3, 2, 2, 1, 2)][0], 96: T[(4, 3, 2, 2, 1, 2)][0],
        128: REFS["ref_na128_dt0.00025_o1e-05_p1e-08"][(2.0, 4.0)][0]}
r = 128 / 96
print("1. Limits beyond the reference (n_arc 64, 96, 128 at dt 2.5e-4; relative to the reference)")
limits = {}
for ob in OBS:
    p = observed_order(fine[128][ob], fine[96][ob], fine[64][ob], r, 96 / 64)
    step = fine[128][ob] - fine[96][ob]
    row = dict(order=p)
    for name, q in (("observed", p), ("second", 2.0)):
        lim = CASE.exact[ob] + step / (r**q - 1)
        row[name] = abs(lim - CASE.exact[ob]) / abs(CASE.exact[ob])
        limits[(ob, name)] = lim
    RES.setdefault("limits", {})[ob] = row
    print(f"  {ob:<7} order {p:.2f}   beyond the reference: {row['observed']:.1e} (observed order), "
          f"{row['second']:.1e} (second order)")
ob, tg = "cl_max", 3e-2
eps = tg * abs(CASE.exact[ob])
o = P.replay(CASE, ob, tg, CFG)
for name in ("observed", "second"):
    errs = [abs(T[q][0][ob] - limits[(ob, name)]) / eps
            for q, n in o["proposals"].items() for _ in range(n)]
    miss = [e for e in errs if e > 1]
    RES.setdefault("lift_3e-2_regraded", {})[name] = dict(accepted=len(errs), missed=len(miss),
                                                           worst=max(errs))
    print(f"  lift at 3e-2 against the limit at the {name} order: {len(miss)} of {len(errs)} "
          f"accepted searches miss, worst {max(errs):.2f}")
print(f"  lift at 3e-2: {tg / CASE.unc[ob]:.1f} times the measured uncertainty")

# 2. The ungraded targets (S-5.2)
print("\n2. The released program at the two ungraded targets")
for ob, tg in (("cd_max", 3e-3), ("cl_max", 1e-2)):
    o = P.replay(CASE, ob, tg, CFG)
    es = [e for _, _, e in o["misses"]]
    RES.setdefault("ungraded", {})[f"{ob} {tg:g}"] = dict(
        accepted=o["accepted"], missed=o["missed"], min=min(es), max=max(es),
        unc_ratio=tg / CASE.unc[ob], settings={str(q): n for q, n in o["proposals"].items()})
    print(f"  {ob} {tg:g}: {tg / CASE.unc[ob]:.1f} times the uncertainty "
          f"(uncertainty {CASE.unc[ob] / tg:.2f} of the target); accepted {o['accepted']}, "
          f"missed {o['missed']}, error {min(es):.3f}-{max(es):.3f} times the target")
    for q, n in sorted(o["proposals"].items()):
        print(f"    {n:3d} return n_arc {CASE.LAD[0][q[0]]}, dt {CASE.LAD[1][q[1]]:g}, "
              f"residualControl {CASE.LAD[2][q[2]]:g}, GAMG {CASE.LAD[3][q[3]]:g}, "
              f"R {CASE.LAD[4][q[4]]:g}, W {CASE.LAD[5][q[5]]:g}")

# 3. The traced search for the drag at 3e-3
ob, tg = "cd_max", 3e-3
eps = tg * abs(CASE.exact[ob])
look = look_for(ob)
kw = dict(CFG["kw"])
NAMES = ["mesh", "time step", "residualControl", "GAMG", "R", "W"]
print("\n3. The traced search for the drag at 3e-3")
for start in ((0, 0, 0, 1, 0, 0), (1, 1, 0, 2, 0, 0)):
    res = run_protocol(look, CASE.LAD, start, eps, hs=CASE.HS, kinds=CASE.KINDS,
                       models=CASE.models, **kw)
    b, q = res["budget"], res["proposal"]
    true = (T[q][0][ob] - CASE.exact[ob]) / eps
    rk = [x / eps for x in b.rk]
    RES.setdefault("traces", {})[str(start)] = dict(
        trace=[list(s) for s in res["trace"]], returned=list(q), accepted=bool(res["certified"]),
        budget=dict(zip(NAMES, map(float, rk))), signed_sum=float(np.nansum(rk)), Fs=float(b.Fs),
        estimate=float(b.E0 / eps), true=float(true), fine_failed=list(map(bool, b.fine_failed)),
        interaction_skipped=bool((b.joint or {}).get("skipped")))
    print(f"  from {start}: path {[tuple(s) for s in res['trace']]}, accepted {res['certified']}")
    print("    budget / eps: " + ", ".join(f"{n} {v:+.3f}" for n, v in zip(NAMES, rk)))
    print(f"    signed sum {np.nansum(rk):.3f}, safety factor {b.Fs:.3f}, estimate {b.E0 / eps:.3f}, "
          f"true error {true:.3f}; finer run failed: "
          f"{[n for n, f in zip(NAMES, b.fine_failed) if f]}; interaction check skipped: "
          f"{bool((b.joint or {}).get('skipped'))} (axes {(b.joint or {}).get('axes')})")
print("  at the two time steps where every finer mesh runs, with the tolerances and windows of each "
      "missed setting: changes along n_arc 24-32-48-64-96 and the error at n_arc 32, in units of the target")
for key, row in RES["ungraded"].items():
    ob, tg = key.split()[0], float(key.split()[1])
    eps = tg * abs(CASE.exact[ob])
    for q in sorted({tuple(int(x) for x in k.strip("()").split(",")) for k in row["settings"]}):
        for dt in (2, 3):
            idx = [(m, dt) + q[2:] for m in range(5)]
            v = [T[i][0][ob] if np.isfinite(T[i][1]) else np.nan for i in idx]
            d = np.diff(v) / eps
            e32 = (v[1] - CASE.exact[ob]) / eps
            RES.setdefault("along_mesh", {})[f"{ob} {tg:g} {q[2:]} dt {CASE.LAD[1][dt]:g}"] = \
                dict(changes=[float(x) for x in d], error_at_32=float(e32),
                     change_32_to_96=float(np.sum(d[1:])))
            print(f"    {ob} {tg:g}, rc/GAMG/R/W {q[2:]}, dt {CASE.LAD[1][dt]:g}: "
                  f"{' '.join(f'{x:+.3f}' for x in d)}; error at 32 {e32:+.3f}; "
                  f"change 32 to 96 {np.sum(d[1:]):+.3f}")

json.dump(RES, open(os.path.join(OUT, "dfg2_traces.json"), "w"), indent=1)
