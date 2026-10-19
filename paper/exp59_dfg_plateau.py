"""DFG 2D-1: is the drag plateau at n_arc >= 48 in the solver or in the post-processing?
(referee point 7 of revision 3)

For n_arc = 8 ... 128 the steady problem is solved to residualControl 1e-8
(iteration limit 20,000), and the drag, lift and pressure difference are
evaluated four ways from the same fields:

  surface-1   the post-processing of the study (as OpenFOAM's forces function
              object: wall-cell pressure, wall shear nu u_t / d);
  surface-2   second-order one-sided wall values from the two wall-normal cells
              (pressure extrapolated to the wall, shear from a quadratic profile);
  cv          a control-volume momentum balance through a closed surface of
              internal faces at radius 0.08 (0.03 from the wall), using the
              solver's own face fluxes phi: no wall values at all;
  dp-1 / dp-2 pressure difference from wall-cell pressures (study) and from
              pressures extrapolated to the wall (second order).

If the cv drag converges where surface-1 plateaus, the plateau is in the
post-processing; if all agree, it is in the discrete solution.
    python exp59_dfg_plateau.py [n_arc ...]
"""
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import openfoam_cylinder as OC
from openfoam_cylinder import (DFGCylinder, REF, NU, UBAR, D, DEPTH, _w, _foam_list,
                               _internal_vectors, _internal_scalars, _force_coeffs,
                               _interp_angle, fv_solution, control_dict, ENV)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "dfg_plateau.json")
WORK = os.environ.get("KNOBS_PLATEAU_WORK", "/tmp/knobs_plateau")


def read_mesh(case):
    pm = os.path.join(case, "constant", "polyMesh")
    n, body = _foam_list(os.path.join(pm, "points"))
    pts = np.array(re.findall(r"\(([-\d.eE+]+) ([-\d.eE+]+) ([-\d.eE+]+)\)", body)[:n], float)
    n, body = _foam_list(os.path.join(pm, "faces"))
    faces = [list(map(int, f.split())) for f in re.findall(r"\d+\(([\d ]+)\)", body)[:n]]
    n, body = _foam_list(os.path.join(pm, "owner"))
    owner = np.array(body.split(")")[0].split()[:n], int)
    n2, body2 = _foam_list(os.path.join(pm, "neighbour"))
    neigh = np.array(body2.split(")")[0].split()[:n2], int)
    ncell = owner.max() + 1
    # face centres and area vectors (quadrilaterals)
    q = np.array([pts[f[:4]] for f in faces])
    cf = q.mean(axis=1)
    sf = 0.5 * np.cross(q[:, 2] - q[:, 0], q[:, 3] - q[:, 1])
    # cell centres: mean of the face centres of each cell, area-weighted
    acc = np.zeros((ncell, 3)); w = np.zeros(ncell)
    mag = np.linalg.norm(sf, axis=1)
    np.add.at(acc, owner, cf * mag[:, None]); np.add.at(w, owner, mag)
    np.add.at(acc, neigh, cf[:len(neigh)] * mag[:len(neigh), None]); np.add.at(w, neigh, mag[:len(neigh)])
    cc = acc / w[:, None]
    btxt = open(os.path.join(pm, "boundary")).read()
    blk = re.search(r"\bcylinder\b\s*\{([^}]*)\}", btxt).group(1)
    nf = int(re.search(r"nFaces\s+(\d+)", blk).group(1))
    st = int(re.search(r"startFace\s+(\d+)", blk).group(1))
    return dict(pts=pts, faces=faces, owner=owner, neigh=neigh, cf=cf, sf=sf, cc=cc,
                wall=np.arange(st, st + nf))


def internal_face_scalars(path, n):
    txt = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", txt)
    return np.array(txt[m.end():].split(")", 1)[0].split()[:n], float)


def evaluate(case, tdir):
    M = read_mesh(case)
    U = _internal_vectors(os.path.join(tdir, "U"))
    p = _internal_scalars(os.path.join(tdir, "p"))
    nI = len(M["neigh"])
    phi = internal_face_scalars(os.path.join(tdir, "phi"), nI)
    q = 0.5 * UBAR**2 * D * DEPTH
    owner, neigh, cf, sf, cc = M["owner"], M["neigh"], M["cf"], M["sf"], M["cc"]

    # ---- control volume: cells with centre radius < R
    R = 0.08
    r = np.hypot(cc[:, 0] - 0.2, cc[:, 1] - 0.2)
    inside = r < R
    oi, ni = inside[owner[:nI]], inside[neigh]
    cut = np.where(oi != ni)[0]
    F = np.zeros(3)
    for f in cut:
        P, N = owner[f], neigh[f]
        S = sf[f].copy(); ph = phi[f]
        if not inside[P]:               # orient S and phi outward from the volume
            S, ph, P, N = -S, -ph, N, P
        dP = np.linalg.norm(cf[f] - cc[P]); dN = np.linalg.norm(cf[f] - cc[N])
        wP = dN / (dP + dN)
        uf = wP * U[P] + (1 - wP) * U[N]
        pf = wP * p[P] + (1 - wP) * p[N]
        nhat = S / np.linalg.norm(S)
        dn = abs(np.dot(cc[N] - cc[P], nhat))
        sng = (U[N] - U[P]) / dn
        F -= uf * ph + pf * S - NU * sng * np.linalg.norm(S)
    cd_cv, cl_cv = F[0] / q, F[1] / q

    # ---- wall: first-order (study) and second-order values
    geo = dict(sf=sf[M["wall"]], cf=cf[M["wall"]], cc=cc[owner[M["wall"]]], cells=owner[M["wall"]])
    cd1, cl1 = _force_coeffs(tdir, geo)
    # second cell along the wall normal: the neighbour across the face of the
    # wall cell that lies opposite the wall face
    cells_faces = {}
    for f in range(nI):
        cells_faces.setdefault(owner[f], []).append(f)
        cells_faces.setdefault(neigh[f], []).append(f)
    Fw = np.zeros(3); pw2 = []; ang = []
    for f in M["wall"]:
        P = owner[f]
        nhat = sf[f] / np.linalg.norm(sf[f])          # into the body
        best, bd = None, -1
        for g in cells_faces[P]:
            dd = np.dot(cf[g] - cf[f], -nhat)
            if dd > bd:
                best, bd = g, dd
        N = neigh[best] if owner[best] == P else owner[best]
        dP = abs(np.dot(cc[P] - cf[f], nhat)); dN = abs(np.dot(cc[N] - cf[f], nhat))
        p_w = p[P] + (p[P] - p[N]) * dP / (dN - dP)
        uP = U[P] - np.dot(U[P], nhat) * nhat
        uN = U[N] - np.dot(U[N], nhat) * nhat
        # u(d) = a d + b d^2 with u(0) = 0: wall derivative a
        a = (uP * dN**2 - uN * dP**2) / (dP * dN * (dN - dP))
        mag = np.linalg.norm(sf[f])
        Fw += p_w * sf[f] + NU * mag * a
        pw2.append(p_w); ang.append(math.atan2(cf[f][1] - 0.2, cf[f][0] - 0.2))
    cd2, cl2 = Fw[0] / q, Fw[1] / q
    ang = np.array(ang); pw2 = np.array(pw2)
    pcell = p[owner[M["wall"]]]
    dp1 = _interp_angle(ang, pcell, math.pi) - _interp_angle(ang, pcell, 0.0)
    dp2 = _interp_angle(ang, pw2, math.pi) - _interp_angle(ang, pw2, 0.0)
    return dict(cd_surface1=cd1, cd_surface2=cd2, cd_cv=cd_cv, cl_surface1=cl1,
                cl_surface2=cl2, cl_cv=cl_cv, dp1=dp1, dp2=dp2, cells=int(len(cc)))


def solve(na, res=1e-8, lin=1e-10, end=20000):
    p = DFGCylinder(workdir=os.path.join(WORK, f"na{na}"), keep=True)
    md = p.mesh(na)[0]
    d = os.path.join(WORK, f"run_{na}")
    shutil.rmtree(d, ignore_errors=True)
    shutil.copytree(md, d)
    _w(f"{d}/system/fvSolution", "dictionary", "fvSolution", fv_solution(res, lin, 0.9))
    _w(f"{d}/system/controlDict", "dictionary", "controlDict", control_dict(end))
    t0 = time.time()
    r = subprocess.run(["simpleFoam"], cwd=d, env=ENV, capture_output=True, text=True)
    log = r.stdout
    it = len(re.findall(r"^Time = ", log, flags=re.M))
    conv = "SIMPLE solution converged" in log
    times = [t for t in os.listdir(d) if re.fullmatch(r"[\d.eE+-]+", t) and t != "0"]
    last = max(times, key=float)
    out = evaluate(d, os.path.join(d, last))
    out.update(n_arc=na, iters=it, converged=conv, wall=time.time() - t0)
    shutil.rmtree(d, ignore_errors=True)
    return out


if __name__ == "__main__":
    levels = [int(x) for x in sys.argv[1:]] or [8, 12, 16, 24, 32, 48, 64, 96, 128]
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for na in levels:
        if str(na) in res:
            continue
        o = solve(na)
        res[str(na)] = o
        json.dump(res, open(OUT, "w"), indent=1)
        print(json.dumps(o), flush=True)
    print(f"\n{'n_arc':>6} {'cells':>7} " + " ".join(f"{k:>12}" for k in
          ["cd_surface1", "cd_surface2", "cd_cv", "dp1", "dp2", "cl_surface1", "cl_cv"]))
    for na in sorted(res, key=int):
        o = res[na]
        row = []
        for k in ["cd_surface1", "cd_surface2", "cd_cv"]:
            row.append(f"{abs(o[k]-REF['cd'])/REF['cd']:12.2e}")
        for k in ["dp1", "dp2"]:
            row.append(f"{abs(o[k]-REF['dp'])/REF['dp']:12.2e}")
        for k in ["cl_surface1", "cl_cv"]:
            row.append(f"{abs(o[k]-REF['cl'])/REF['cl']:12.2e}")
        print(f"{na:>6} {o['cells']:>7} " + " ".join(row)
              + f"   iters {o['iters']}{'' if o['converged'] else ' (limit)'}")
