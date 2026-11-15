"""Revision 4 (referee point on relTol): when does the GAMG tolerance of DFG 2D-1 act?

Both linear solvers keep relTol 0.1, as in the OpenFOAM tutorials, so a solve stops at
max(tolerance, 0.1 x its initial residual), and does no iteration at all if its initial
residual is already below the tolerance.  The outer SIMPLE iteration stops when the
initial residual of p falls below residualControl and that of U below a tenth of it.
The study's "GAMG tolerance" t sets the pressure tolerance t and the momentum tolerance
10 t, so relative to its outer threshold the momentum tolerance is 100 t / residualControl
and the pressure tolerance t / residualControl.  The table shows which of the two ratios
decides whether the tolerance acts: the observables change with t until 100 t /
residualControl falls below about 0.1 (the momentum solve then always stops by relTol),
and the outer iteration fails to converge once t / residualControl reaches 1 (the
pressure solve then does no iteration before the outer residual reaches its threshold).
This script reads the 2D-1 tensor (no runs) and reports, for every pair of
residualControl and GAMG tolerance: the fraction of settings that did not converge,
and, for converged settings, the change of each observable from the tightest GAMG
tolerance at the same mesh, residualControl and coarsest level, relative to the
exact value (median and largest over meshes and coarsest levels).
    python exp77_reltol_2d1.py         # writes data/reltol_2d1.json
The Hamilton job J1 (tasks rt_*) measures the stopping criterion of every pressure
solve directly from the solver log, with relTol 0.1 and 0.
"""
import json
import os

import numpy as np

import cases as C


def main():
    c = C._dfg1()
    LAD, T = c.LAD, c.tensor
    RC, TOL = LAD[1], LAD[2]
    out = []
    print("columns: residualControl, GAMG tolerance, ratio t/rc (pressure), 100 t/rc (momentum), settings, not converged, "
          "and for converged ones the change from the tightest GAMG tolerance relative to the "
          "exact value: median / largest, for c_D, Delta p, c_L")
    for b, rc in enumerate(RC):
        for g, tol in enumerate(TOL):
            n = nc = 0
            dev = {o: [] for o in ("cd", "dp", "cl")}
            for a in range(len(LAD[0])):
                for d in range(len(LAD[3])):
                    i = (a, b, g, d)
                    ref = (a, b, len(TOL) - 1, d)
                    if not c.usable(i):
                        continue
                    n += 1
                    if not c.converged(i):
                        nc += 1
                        continue
                    if not (c.usable(ref) and c.converged(ref)):
                        continue
                    for o in dev:
                        dev[o].append(abs(T[i][0][o] - T[ref][0][o]) / abs(c.exact[o]))
            row = dict(res_ctrl=rc, gamg_tol=tol, ratio=tol / rc, ratio_u=100 * tol / rc,
                       settings=n, not_converged=nc,
                       **{f"{o}_med": float(np.median(v)) if v else None for o, v in dev.items()},
                       **{f"{o}_max": float(np.max(v)) if v else None for o, v in dev.items()})
            out.append(row)
            f = lambda x: "    ---  " if x is None else f"{x:9.1e}"
            print(f"{rc:7.0e} {tol:7.0e} {tol/rc:8.0e} {100*tol/rc:8.0e} {n:4d} {nc:4d}   "
                  + "  ".join(f"{f(row[o+'_med'])}/{f(row[o+'_max'])}" for o in ("cd", "dp", "cl")))
    # summary by regime
    for name, sel in [("100 t/rc < 0.1 (momentum tolerance inactive)", lambda r: r["ratio_u"] < 0.1 - 1e-12),
                      ("100 t/rc >= 0.1 and t/rc < 1", lambda r: r["ratio_u"] >= 0.1 - 1e-12 and r["ratio"] < 1 - 1e-12),
                      ("t/rc >= 1 (pressure tolerance at or above the outer threshold)", lambda r: r["ratio"] >= 1 - 1e-12),
                      ("t/rc < 1/10 (the regime in which the referee expected no effect)", lambda r: r["ratio"] < 0.1 - 1e-12)]:
        R = [r for r in out if sel(r) and r["gamg_tol"] != TOL[-1]]   # the tightest is the reference
        n, nc = sum(r["settings"] for r in R), sum(r["not_converged"] for r in R)
        mx = {o: max([r[o + "_max"] for r in R if r[o + "_max"] is not None], default=None) for o in ("cd", "dp", "cl")}
        print(f"{name}: {n} settings, {nc} not converged; largest change from the tightest "
              f"tolerance among converged runs: " + ", ".join(f"{o} {v:.1e}" if v is not None else f"{o} ---"
                                                               for o, v in mx.items()))
    json.dump(out, open(os.path.join(C.D, "reltol_2d1.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
