"""Grading the method against the baselines of baselines.py."""
import itertools, json, os, pickle, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from knobs import run_protocol
from knobs.problems import Burgers1D, Burgers2D, ReactionDiffusion2D
from baselines import seq, coord, ref_check

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

# ------------------------------------------------------------------ grading
def grade(tensor, LAD, HS, KIN, ob, QS, target, rel):
    ALL = list(itertools.product(*[range(len(l)) for l in LAD]))
    eps = target * abs(QS) if rel else target
    err = {i: abs(tensor[i][0][ob] - QS) for i in ALL}
    feas = [i for i in ALL if err[i] <= eps]
    if not feas or len(feas) == len(ALL):
        return None
    c_opt = min(tensor[i][1] for i in feas)
    look = lambda i: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
    out = {}

    # the method
    safe = cert = 0
    pens, runs = [], []
    for a in ALL:
        r = run_protocol(look, LAD, a, eps, hs=HS, kinds=KIN, use_screen=False)
        if not r["certified"]:
            continue
        cert += 1
        safe += err[r["proposal"]] <= eps
        pens.append(tensor[r["proposal"]][1] / c_opt)
        runs.append(r["runs"])
    out["method"] = dict(n=len(ALL), certified=cert, safe=safe,
                         pen=float(np.median(pens)) if pens else np.nan,
                         runs=float(np.median(runs)) if runs else np.nan)

    # SEQ and COORD at three thresholds
    for name, fn in [("seq", seq), ("coord", coord)]:
        for frac, tag in [(1.0, "eps"), (1 / 3, "eps/3"), (0.1, "eps/10")]:
            s, pn, rn = 0, [], []
            for a in ALL:
                q, r = fn(look, LAD, a, eps * frac)
                s += err[q] <= eps
                pn.append(tensor[q][1] / c_opt)
                rn.append(r)
            out[f"{name}@{tag}"] = dict(n=len(ALL), certified=len(ALL), safe=s,
                                        pen=float(np.median(pn)),
                                        runs=float(np.median(rn)))

    # reference check (start-independent, so one outcome)
    q, r = ref_check(look, LAD, ALL, eps)
    out["ref"] = dict(n=len(ALL), certified=len(ALL),
                      safe=len(ALL) * int(err[q] <= eps),
                      pen=tensor[q][1] / c_opt, runs=r)
    return out


# ------------------------------------------------------------------ problems
b1 = Burgers1D(nu=0.05, T=1.0)
t1 = pickle.load(open(os.path.join(D, "burgers_tensor.pkl"), "rb"))
LAD1 = [[32, 64, 128, 256, 512], [25, 50, 100, 200, 400],
        [1e-2, 1e-4, 1e-6, 1e-8, 1e-10], [1e-2, 1e-4, 1e-6, 1e-8, 1e-10]]
HS1 = [[2 * np.pi / n for n in LAD1[0]], [1.0 / s for s in LAD1[1]],
       list(LAD1[2]), list(LAD1[3])]
KIN1 = ["discretisation", "discretisation", "nested-solver", "nested-solver"]

b2 = Burgers2D(nu=0.05, T=1.0)
t2 = pickle.load(open(os.path.join(D, "burgers2d_tensor.pkl"), "rb"))
LAD2 = [[16, 32, 64, 128], [10, 20, 40, 80], [1e-2, 1e-5, 1e-8],
        [1e-2, 1e-5, 1e-8], [1e-1, 1e-2, 1e-3]]
HS2 = [[2 * np.pi / n for n in LAD2[0]], [1.0 / s for s in LAD2[1]],
       list(LAD2[2]), list(LAD2[3]), list(LAD2[4])]
KIN2 = ["discretisation", "discretisation", "nested-solver", "nested-solver",
        "cost-only"]

p4 = ReactionDiffusion2D()
t4 = pickle.load(open(os.path.join(D, "reactdiff_tensor.pkl"), "rb"))
LAD4 = [[16, 32, 64, 128], [5, 10, 20, 40], [1, 2, 4, 8],
        [1e-2, 1e-5, 1e-8], [1e-1, 1e-2, 1e-3]]
HS4 = [[2 * np.pi / n for n in LAD4[0]], [1.0 / s for s in LAD4[1]],
       [1.0 / r for r in LAD4[2]], list(LAD4[3]), list(LAD4[4])]
KIN4 = ["discretisation", "discretisation", "discretisation",
        "nested-solver", "cost-only"]

CELLS = [
    ("L2 energy", t1, LAD1, HS1, KIN1, "energy", b1.exact_observables(),
     [1e-2, 3e-3, 1e-3], False),
    ("L2 dissipation", t1, LAD1, HS1, KIN1, "dissipation",
     b1.exact_observables(), [1e-2, 3e-3, 1e-3], False),
    ("L3 energy", t2, LAD2, HS2, KIN2, "energy", b2.exact_observables(),
     [1e-2, 3e-3], True),
    ("L3 gradient", t2, LAD2, HS2, KIN2, "gradient", b2.exact_observables(),
     [1e-2, 3e-3], True),
    ("L4 mass", t4, LAD4, HS4, KIN4, "inhibitor_mass", p4.exact_observables(),
     [1e-2, 3e-3], True),
    ("L4 sharpness", t4, LAD4, HS4, KIN4, "pattern_sharpness",
     p4.exact_observables(), [1e-2, 3e-3], True),
]

METHODS = ["method", "seq@eps", "seq@eps/3", "seq@eps/10",
           "coord@eps", "coord@eps/3", "coord@eps/10", "ref"]

rows = []
print("=" * 108)
print("the method against three baselines, same tensors, same starting "
      "configurations")
print("=" * 108)
for name, T, LAD, HS, KIN, ob, EX, targets, rel in CELLS:
    for t in targets:
        r = grade(T, LAD, HS, KIN, ob, EX[ob], t, rel)
        if r is None:
            continue
        rows.append(dict(cell=f"{name} {t:.0e}", **r))
        print(f"\n{name}, target {t:.0e}{' rel' if rel else ''}")
        print(f"  {'procedure':<14} {'certified':>10} {'attained':>9} "
              f"{'safety':>8} {'med penalty':>12} {'med runs':>9}")
        for m in METHODS:
            d = r[m]
            print(f"  {m:<14} {d['certified']:>10} {d['safe']:>9} "
                  f"{100*d['safe']/d['n']:>7.1f}% {d['pen']:>11.2f}x "
                  f"{d['runs']:>9.0f}")

json.dump(rows, open(os.path.join(D, "baselines.json"), "w"), indent=1)

print()
print("=" * 108)
print(f"  {'procedure':<14} {'cells':>6} {'total attained':>15} {'safety':>9} "
      f"{'med penalty':>12} {'med runs':>9}")
for m in METHODS:
    n = sum(r[m]["n"] for r in rows)
    s = sum(r[m]["safe"] for r in rows)
    pen = float(np.median([r[m]["pen"] for r in rows]))
    run = float(np.median([r[m]["runs"] for r in rows]))
    print(f"  {m:<14} {len(rows):>6} {s:>8}/{n:<6} {100*s/n:>8.1f}% "
          f"{pen:>11.2f}x {run:>9.0f}")
print()
print("  Safety is over ALL starts, not only certified ones, so that the")
print("  method's refusals count against it on the same footing as a")
print("  baseline's wrong answers.")
