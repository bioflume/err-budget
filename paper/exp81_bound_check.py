"""Conditions of Proposition 1 tested on synthetic problems (S-1.9).  Synthetic errors on
candidate lists of constant ratio; the released program (errbudget.V121) searches from every
setting.

    python exp81_bound_check.py         # about three minutes

1. Lemma 1: the bound of the observed-order check against the geometric extrapolation
   |D_j0| rho^(a - j0) / (1 - rho), rho = D_j1 / D_j0, at each choice of positions.
2. Corollary 1 at targets 1% above and below 1.25 |E(theta_f)|.
3. (A1): higher-order terms of the sign of the leading term, and of the opposite sign.
4. (A3): a bilinear interaction b u_k u_l, with and without the interaction check.
5. An interaction that reverses the sign of the error at a parameter's finest value.
"""
import itertools
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
from errbudget import LocalBudget, V121, run_protocol  # noqa: E402

PG = np.linspace(0.25, 6.0, 121)                     # order grid of the fit (Section 4.2)

# 1. Lemma 1
rng = np.random.default_rng(1)
LAD1 = [[8, 16, 32, 64, 128, 256], [1e-2, 1e-4, 1e-6]]
HS1 = [[1.0 / v for v in LAD1[0]], list(LAD1[1])]
K1 = ["discretisation", "nested-solver"]
maxdev, n_checked, aligned = 0.0, 0, set()
for trial in range(300):
    rat = rng.uniform(0.05, 0.9, size=8)              # changes of one sign, arbitrary ratios
    D = np.cumprod(np.r_[1.0, rat[:-1]]) * rng.uniform(0.5, 2.0)
    tail = D[5] * rat[5] / (1 - rat[5])
    err = [np.sum(D[i:6]) + tail for i in range(6)]   # Q(h_i) - Q*
    for a, fail in [(a, False) for a in range(6)] + [(a, True) for a in range(2, 5)]:
        # fail: the run one position finer than the current one fails (Section 4.2)
        look = lambda i, a=a, fail=fail: ((np.nan, 1.0) if fail and i[0] == a + 1
                                          else (err[i[0]], 1.0 + i[0] + i[1]))
        lb = LocalBudget(look, LAD1, (a, 0), hs=HS1, kinds=K1, local_order=True, resolution=1e-30,
                         cancel_rule="none", joint_res=1.0, fs_cap=1.25, resid_margin=5.0)
        lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
        L = lb.local.get(0)
        if L is None or L["unbounded"] or not np.isfinite(L["p"]):
            continue
        j0, j1, j2 = L["pts"]
        rho = (err[j1] - err[j2]) / (err[j0] - err[j1])
        geo = abs(err[j0] - err[j1]) * rho ** (a - j0) / (1 - rho)
        maxdev = max(maxdev, abs(L["E"] - geo) / geo)
        n_checked += 1; aligned.add((j0 - a, j1 - a, j2 - a))
print("1. Lemma 1: local-order bound against the geometric extrapolation")
print(f"   {n_checked} checks at the alignments {sorted(aligned)}: largest relative difference {maxdev:.1e}")

# the synthetic problems of 2-4
LAD = [[8, 16, 32, 64, 128], [10, 20, 40, 80, 160], [1e-2, 1e-4, 1e-6, 1e-8], [1, 2, 3, 4]]
HS = [[1.0 / v for v in LAD[0]], [1.0 / v for v in LAD[1]], list(LAD[2]), [1.0, 0.5, 0.25, 0.125]]
KINDS = ["discretisation", "discretisation", "nested-solver", "cost-only"]
W = [1.3, 1.0, 1.1, 1.5]                              # cost-only parameter: interior minimum
GRID = list(itertools.product(*[range(len(l)) for l in LAD]))
FIN = tuple(len(l) - 1 for l in LAD)


def cost(i):
    return LAD[0][i[0]] ** 2 * LAD[1][i[1]] * (-np.log10(LAD[2][i[2]])) * W[i[3]]


def replay(Qfun, eps):
    """Replay from every start; Q* = 1.  Returns per-start records."""
    T = {i: Qfun(i) for i in GRID}
    look = lambda i: (T[tuple(i)], cost(tuple(i)))
    out = []
    for a in GRID:
        r = run_protocol(look, LAD, a, eps, hs=HS, kinds=KINDS, **V121)
        q = tuple(r["proposal"])
        lb = r["budget"]
        out.append(dict(start=a, acc=bool(r["certified"]), reason=r.get("declined"), q=q,
                        E=T[q] - 1.0, Ehat=float(lb.E0), Fs=float(lb.Fs),
                        acc_setting=tuple(lb.anchor), E_acc=T[tuple(lb.anchor)] - 1.0,
                        joint=lb.joint))
    return T, out


def summary(T, out, eps):
    acc = [o for o in out if o["acc"]]
    dec = {}
    for o in out:
        if not o["acc"]:
            dec[o["reason"]] = dec.get(o["reason"], 0) + 1
    miss = [o for o in acc if abs(o["E"]) > eps]
    # Proposition 1 at the accepted setting (before the cost-only pass): |E| <= |Ehat|,
    # and |E| <= |Ehat| / Fs when the model is additive
    over = [o for o in acc if abs(o["E_acc"]) > abs(o["Ehat"]) * (1 + 1e-9)]
    ratio = max((abs(o["E_acc"]) * o["Fs"] / abs(o["Ehat"]) for o in acc), default=0.0)
    worst = max((abs(o["E"]) / eps for o in acc), default=0.0)
    return acc, dec, miss, over, ratio, worst


u = lambda i, k, p: HS[k][i[k]] ** p[k]

# 2. Corollary 1
p = [PG[37], PG[16], PG[16], 1.0]
c = [1.0, 0.5, 2.0, 0.0]
QA = lambda i: 1.0 + sum(c[k] * u(i, k, p) for k in range(4))
Efin = QA(FIN) - 1.0
print("\n2. Corollary 1: exact separable model, one sign, orders on the grid "
      f"(p = {p[0]:.4f}, {p[1]:.4f}, {p[2]:.4f}); |E(finest)| = {Efin:.4e}")
for eps in (1e-2, 5e-3, Efin / 0.8 * 1.01, Efin / 0.8 * 0.99, 3e-3):
    T, out = replay(QA, eps)
    acc, dec, miss, over, ratio, worst = summary(T, out, eps)
    viol = sum(1 for o in out if o["acc"] != (abs(Efin) <= 0.8 * eps))
    exact = max((abs(abs(o["Ehat"]) - 1.25 * abs(o["E_acc"])) / abs(o["Ehat"]) for o in acc), default=0.0)
    print(f"   eps {eps:.4e} (finest meets 0.8 eps: {abs(Efin) <= 0.8 * eps}): accepted {len(acc)}, "
          f"declined {dec}, missed {len(miss)}, worst |E|/eps {worst:.3f}; corollary violated at {viol} "
          f"starts; largest relative deviation of |Ehat| from 1.25 |E| {exact:.1e}")

# 3. Proposition 1, condition (A1)
pB = [2.0, 1.0, PG[16], 1.0]                          # leading orders, off the grid of the fit
QB = lambda i: 1.0 + 1.0 * (u(i, 0, pB) - 0.4 * u(i, 0, pB) ** 1.5) \
    + 0.5 * (u(i, 1, pB) - 0.4 * u(i, 1, pB) ** 2) + 2.0 * u(i, 2, pB)
QC = lambda i: 1.0 + 1.0 * (u(i, 0, pB) + 3.0 * u(i, 0, pB) ** 1.5) \
    + 0.5 * (u(i, 1, pB) + 3.0 * u(i, 1, pB) ** 2) + 2.0 * u(i, 2, pB)
for name, Q in (("opposite-sign higher-order terms (convergence speeds up)", QB),
                ("same-sign higher-order terms (convergence slows down)", QC)):
    # ratios of successive changes along the mesh, at the finest values of the others
    D = [Q((j,) + FIN[1:]) - Q((j + 1,) + FIN[1:]) for j in range(len(LAD[0]) - 1)]
    print(f"\n3. Proposition 1, (A1): {name}; |E(finest)| = {Q(FIN) - 1.0:.4e}")
    print("   ratios of successive changes along the first parameter",
          [f"{b / a:.3f}" for a, b in zip(D, D[1:])])
    for eps in (1e-2, 5e-3, 4e-3, 3e-3):
        T, out = replay(Q, eps)
        acc, dec, miss, over, ratio, worst = summary(T, out, eps)
        print(f"   eps {eps:.0e}: accepted {len(acc)}, declined {sum(dec.values())}, missed {len(miss)}, "
              f"worst |E|/eps {worst:.3f}; |E| > |Ehat| at {len(over)} accepted settings; "
              f"max |E| Fs/|Ehat| {ratio:.3f}")

# 4. Proposition 1, condition (A3), joint check
print("\n4. Proposition 1, (A3): bilinear interaction b u0 u1 (u_k = h_k^p_k)")


def estimate(look, idx, eps, joint):
    """The estimate of run_protocol at one setting (errbudget.V121), with or without the joint check."""
    lb = LocalBudget(look, LAD, idx, hs=HS, fs_cap=1.25, resid_margin=5.0, kinds=KINDS,
                     local_order=True, resolution=0.1 * eps, cancel_rule="none", joint_res=1.0)
    lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
    if joint:
        lb.joint_correct(look, eps, cancel_rule="none")
    return lb


p4 = [PG[37], PG[37], PG[16], 1.0]
rho4 = [0.5 ** p4[0], 0.5 ** p4[1]]
print(f"   (a) every setting, orders {p4[0]:.4f} and {p4[1]:.4f} on lists halved at each step: "
      f"rho = {rho4[0]:.4f}, (1 - rho0)(1 - rho1) = {(1 - rho4[0]) * (1 - rho4[1]):.3f}; eps 5e-3")
c4 = [1.0, 1.0, 0.05, 0.0]
for b in (-60.0, -40.0, -20.0, 20.0, 60.0):
    Q4 = lambda i: 1.0 + sum(c4[k] * u(i, k, p4) for k in range(4)) + b * u(i, 0, p4) * u(i, 1, p4)
    T = {i: (Q4(i), cost(i)) for i in GRID}
    look = lambda i: T[tuple(i)]
    eps = 5e-3
    rows = {}
    for q in GRID:
        E = T[q][0] - 1.0
        Eb = [(c4[0] + b * u(q, 1, p4)) * u(q, 0, p4), (c4[1] + b * u(q, 0, p4)) * u(q, 1, p4),
              c4[2] * u(q, 2, p4)]
        lbj = estimate(look, q, eps, True)
        lb0 = estimate(look, q, eps, False)
        one_sign = len({np.sign(x) for x in Eb if x}) <= 1
        A3 = abs(E - sum(Eb)) <= (lb0.Fs - 1) * sum(abs(x) for x in Eb)
        J = lbj.joint or {}
        both = (not J.get("skipped", True)) and {0, 1} <= set(J.get("axes", []))
        rows.setdefault((one_sign, A3, both), []).append(
            (abs(E) > abs(lbj.E0) * (1 + 1e-9), abs(E) > abs(lb0.E0) * (1 + 1e-9)))
    print(f"   b = {b:+.0f}:")
    for (s1, a3, both), v in sorted(rows.items()):
        print(f"      one sign {s1!s:<5} (A3) {a3!s:<5} both in the joint check {both!s:<5}: {len(v):3d} "
              f"settings, |E| > |Ehat| at {sum(x for x, _ in v)} with the joint check and at "
              f"{sum(y for _, y in v)} without it")

print("   (b) the search, orders 2.02 and 1.02, b = -150, eps 3e-3")
pp = [PG[37], PG[16], PG[16], 1.0]
c5 = [1.0, 0.5, 0.05, 0.0]
b = -150.0
Q5 = lambda i: 1.0 + sum(c5[k] * u(i, k, pp) for k in range(4)) + b * u(i, 0, pp) * u(i, 1, pp)
eps = 3e-3
T, out = replay(Q5, eps)
acc, dec, miss, over, ratio, worst = summary(T, out, eps)
print(f"      accepted {len(acc)}, missed {len(miss)}, worst |E|/eps {worst:.3f}")
rows = {}
for o in acc:
    q = o["acc_setting"]
    Eb = [(c5[0] + b * u(q, 1, pp)) * u(q, 0, pp), (c5[1] + b * u(q, 0, pp)) * u(q, 1, pp),
          c5[2] * u(q, 2, pp)]
    Ex = o["E_acc"] - sum(Eb)
    J = o["joint"] or {}
    key = (q, round(Ex / eps, 3), tuple(float(round(x / eps, 3)) for x in Eb), round(o["Ehat"] / eps, 3),
           round(o["E_acc"] / eps, 4), bool(J.get("skipped", True)), tuple(J.get("axes", [])))
    rows[key] = rows.get(key, 0) + 1
for (q, ex, eb, eh, ea, sk, ax), n in sorted(rows.items()):
    fin = [k for k in range(2) if q[k] == FIN[k]]
    print(f"      {n:3d} starts accept {q} (finest in parameter {fin}): E/eps {ea:+.4f}, Ehat/eps {eh:+.3f}, "
          f"Ebar/eps {list(eb)}, E_x/eps {ex:+.3f}, joint check skipped {sk}, axes {list(ax)}")
