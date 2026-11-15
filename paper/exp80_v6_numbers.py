"""Numbers of the text that no table holds.

    python exp80_v6_numbers.py

1. Section 7.4, S-5.5: the error budgets of ABINIT (12 Ha, 3^3 k-points) and Quantum
   ESPRESSO (16 Ry, 3^3 k-points, default threshold and mixing), with the released program,
   with its fits alone and with version 1.1; the cheapest ABINIT setting that meets 1e-2 in
   the energy and its force error; the energy and force errors along the ABINIT cutoffs at
   6^3 k-points.
2. The Phase 1 gains where the cost-only coarsest level of 2D-1 moves in the search of
   Fig. S-5.2.
3. L5: the setting at which version 1.2 with cancellation allowed missed, and what the
   interaction check measures there.
4. S-3.2: the L5 streamfunction minimum, its errors along N at the tightest tolerances and
   the settings that meet 1e-2 and 3e-3.
5. S-3.3: the y-coordinate of the L5 minimum at 3e-4 at the settings that miss with a
   constant safety factor (from data/pooled_v1.2.1-constant.json).

L3-L5 take their exact values from SymPy; without it, from data/accounting_v6.json
(exp76_accounting.py).  The printed headings keep the section numbers of revision 4.
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
try:
    import sympy  # noqa: F401
except ImportError:                                   # fall back to the stored references
    import grids
    _acc = {a["case"]: a for a in json.load(open(os.path.join(HERE, "data", "accounting_v6.json")))}
    _orig = grids.exact
    grids.exact = lambda P: dict(_acc[P["name"]]["reference"]) if P["name"] in _acc else _orig(P)
import cases as C
from errbudget import V12, V121, LocalBudget, run_protocol


def lookup(case, ob):
    T = case.tensor

    def look(i):
        v = T[tuple(i)]
        if not np.isfinite(v[1]):
            return np.nan, np.inf
        if not case.converged(tuple(i)):
            return np.nan, v[1]
        return v[0][ob], v[1]
    return look


def budget(case, start, ob, tg, version):
    eps = tg * abs(case.exact[ob])
    look = lookup(case, ob)
    if version in ("v1.2.1", "fits"):
        lb = LocalBudget(look, case.LAD, start, hs=case.HS, fs_cap=1.25, resid_margin=5.0,
                         models=case.models, kinds=case.KINDS, local_order=(version == "v1.2.1"),
                         resolution=0.1 * eps, cancel_rule="none", joint_res=1.0)
    else:                                             # version 1.1, as in pooled.py
        lb = LocalBudget(look, case.LAD, start, hs=case.HS, models=case.models)
    lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
    fin = [abs(x) for x in lb.rk if np.isfinite(x)]
    tot = sum(fin) or 1.0
    return ["unbounded" if not np.isfinite(x) else f"{100 * abs(x) / tot:.1f}%" for x in lb.rk]


print("1. Section 7.4: error budgets (shares of the parameters, in the order of Table 3)")
ab, qe = C._abinit(), C._qe()
for name, case, start, cells in (("ABINIT 12 Ha, 3^3", ab, (1, 1, 1, 1), [("energy", 1e-2), ("force", 1e-1)]),
                                 ("QE 16 Ry, 3^3, defaults", qe, (1, 1, 1, 2), [("energy", 1e-3), ("force", 1e-1)])):
    for ob, tg in cells:
        for ver in ("v1.2.1", "fits", "v1.1"):
            print(f"   {name:<24} {ob:<7} {ver:<7} {budget(case, start, ob, tg, ver)}")
T = ab.tensor
rel = lambda i, ob: abs(T[i][0][ob] - ab.exact[ob]) / abs(ab.exact[ob])
ok = [i for i in T if ab.converged(i) and rel(i, "energy") <= 1e-2]
best = min(ok, key=lambda i: T[i][1])
print(f"   cheapest ABINIT setting with energy error <= 1e-2: "
      f"{[ab.LAD[k][j] for k, j in enumerate(best)]}, energy {rel(best, 'energy'):.2e}, "
      f"force {rel(best, 'force'):.2e}")
k6 = ab.LAD[1].index(6) if 6 in ab.LAD[1] else len(ab.LAD[1]) - 1
row = [(j, k6, len(ab.LAD[2]) - 1, 1) for j in range(len(ab.LAD[0]))]
e = [rel(i, "energy") for i in row]
print("   ABINIT at 6^3: energy errors along the cutoff", [f"{x:.3g}" for x in e])
print("   ratios per step", [f"{b / a:.2f}" for a, b in zip(e, e[1:])])
print("   force errors along the cutoff", [f"{rel(i, 'force'):.2g}" for i in row])

print("\n2. Section 4.6, Fig. 3: Phase 1 gains where the coarsest level moves")
c = C._dfg1(cost="deterministic")
ob = "cl"; eps = 3e-2 * abs(c.exact[ob]); look = lookup(c, ob)
for step, a in ((3, (6, 0, 2, 0)), (7, (6, 1, 4, 1))):
    lb = LocalBudget(look, c.LAD, a, hs=c.HS, kinds=c.KINDS, models=c.models, local_order=True,
                     resid_margin=5.0, resolution=0.1 * eps, cancel_rule="none", joint_res=1.0,
                     fs_cap=1.25)
    lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
    print(f"   before step {step}: {[c.LAD[k][j] for k, j in enumerate(a)]}, estimate/eps {lb.E0 / eps:.3g}")
    for k in range(len(a)):
        j = a[k] + 1
        if j >= len(c.LAD[k]) or j not in lb.sampled[k]:
            continue
        d, g = lb.e[k][j], lb.g[k][j]
        red = abs(lb.E0) - abs(lb.E0 + d)
        print(f"      parameter {k}: one-step change/eps {d / eps:+.4f}, d ln C {g:+.4f}, "
              f"predicted reduction/eps {red / eps:+.4f}, G {red / max(g, 1e-12):.3g}")

print("\n3. Section 6.1: L5, crediting cancellation in versions 1.2 and 1.2.1")
L5 = [c for c in C._manufactured() if c.key == "L5"][0]
ob, tg = "y_min", 3e-4
eps = tg * abs(L5.exact[ob]); look = lookup(L5, ob)
for name, jr in (("v1.2, credit", 3.0), ("v1.2.1, credit", 1.0)):
    lb = LocalBudget(look, L5.LAD, (5, 1, 1, 1), hs=L5.HS, fs_cap=1.25, resid_margin=5.0,
                     models=L5.models, kinds=L5.KINDS, local_order=True, resolution=0.1 * eps,
                     cancel_rule="credit", joint_res=jr)
    lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
    E0 = lb.E0
    lb.joint_correct(look, eps, cancel_rule="credit")
    steps = [lb.e[k][lb.anchor[k] + 1] / eps if lb.anchor[k] + 1 < lb.m[k] else None
             for k in range(len(L5.LAD))]
    print(f"   {name}: setting {[L5.LAD[k][j] for k, j in enumerate((5, 1, 1, 1))]}, one-step "
          f"changes/eps {[None if s is None else round(float(s), 3) for s in steps]}, estimate/eps "
          f"{E0 / eps:.3f} before and {lb.E0 / eps:.3f} after the joint check")
    J = lb.joint
    if J.get("skipped"):
        print("      joint check skipped (fewer than two parameters above the threshold)")
    else:
        print(f"      joint check over parameters {J['axes']}: interaction I/eps {float(J['I']) / eps:.3f}, "
              f"strong {bool(J['strong'])}")
for name, kw in (("v1.2, credit", dict(V12, cancel_rule="credit")),
                 ("v1.2.1, credit", dict(V121, cancel_rule="credit"))):
    out = run_protocol(look, L5.LAD, (5, 1, 1, 1), eps, hs=L5.HS, kinds=L5.KINDS, models=L5.models, **kw)
    q = tuple(out["proposal"])
    print(f"   {name}: from (5, 1, 1, 1) accepted {out['certified']}, returns {q}, true error/eps "
          f"{abs(L5.tensor[q][0][ob] - L5.exact[ob]) / eps:.3f}")

print("\n4. Section 6.1: the L5 streamfunction minimum")
T = L5.tensor; ex = L5.exact["psi_min"]
r5 = lambda i: abs(T[i][0]["psi_min"] - ex) / abs(ex)
fin = tuple(len(l) - 1 for l in L5.LAD)
err = [r5((a,) + fin[1:]) for a in range(len(L5.LAD[0]))]
print("   N", L5.LAD[0])
print("   errors at the tightest tolerances", [f"{x:.2e}" for x in err])
print("   observed orders", [f"{np.log(a / b) / np.log(L5.LAD[0][i + 1] / L5.LAD[0][i]):.2f}"
                             for i, (a, b) in enumerate(zip(err, err[1:]))])
for tg in (1e-2, 3e-3):
    ok = sorted(i for i in T if L5.converged(i) and r5(i) <= tg)
    print(f"   target {tg:g}: {len(ok)} converged settings meet it; coarsest N among them "
          f"{min(L5.LAD[0][i[0]] for i in ok)}" + (f"; settings {ok}" if len(ok) <= 5 else ""))

print("\n5. Section 4.8, Remark: L5 y-coordinate at 3e-4, the settings missed with a constant safety factor")
ob, tg = "y_min", 3e-4
eps = tg * abs(L5.exact[ob]); look = lookup(L5, ob)
cell = [c for c in json.load(open(os.path.join(HERE, "data", "pooled_v1.2.1-constant.json")))
        if c["case"] == "L5" and c["ob"] == ob and abs(c["target"] - tg) < 1e-12][0]
T = L5.tensor
for q in sorted({tuple(m[1]) for m in cell["miss_list"]}):
    E = T[q][0][ob] - L5.exact[ob]
    print(f"   setting {[L5.LAD[k][j] for k, j in enumerate(q)]}: true error/eps {E / eps:+.3f}")
    for name, rm in (("constant factor", 0.0), ("gated factor", 5.0)):
        lb = LocalBudget(look, L5.LAD, q, hs=L5.HS, fs_cap=1.25, resid_margin=rm, models=L5.models,
                         kinds=L5.KINDS, local_order=True, resolution=0.1 * eps, cancel_rule="none",
                         joint_res=1.0)
        lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
        lb.joint_correct(look, eps, cancel_rule="none")
        fin = [abs(x) for x in lb.rk if np.isfinite(x)]
        print(f"      {name}: r_k/eps {[round(float(x) / eps, 3) for x in lb.rk]}, F_s {lb.Fs:.3f}, "
              f"estimate/eps {lb.E0 / eps:+.3f}, |E|/sum|r_k| {abs(E) / sum(fin):.2f}, "
              f"local-order check {'none' if not lb.local else sorted(lb.local)}")
    ch = []
    for j in range(len(L5.LAD[0]) - 1):
        a_ = tuple(j if k == 0 else q[k] for k in range(len(q)))
        b_ = tuple(j + 1 if k == 0 else q[k] for k in range(len(q)))
        ch.append((T[a_][0][ob] - T[b_][0][ob]) / eps)
    print("      changes D along N (/eps)", [round(float(x), 3) for x in ch], "N", L5.LAD[0])
