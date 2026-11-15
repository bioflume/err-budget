"""When the GAMG tolerance t of OpenFOAM 2D-1 acts (Section 7.1, S-5.1).  Each linear solve
stops at max(tolerance, 0.1 x its initial residual) and does no iteration if its initial
residual is already below the tolerance; the outer iteration stops when the pressure
residual falls below residualControl r_c and the momentum residual below r_c/10.  The runs
are sorted into the regimes 100 t/r_c < 0.1, 100 t/r_c >= 0.1 with t < r_c, and t >= r_c,
and each observable is compared with the run at t = 1e-8.

    python exp77_reltol_2d1.py          # writes data/reltol_2d1.json
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
