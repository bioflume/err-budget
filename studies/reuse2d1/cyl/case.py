"""The article's 2D-1 case (DFG cylinder in a channel, steady, simpleFoam) at any Reynolds
number.

Ported from `paper/openfoam_cylinder.py` of the article repository (commit 4d32148): the
block mesh, the case files, the iteration cap and the post-processing are the article's,
unchanged; the inflow maximum U_m = 0.015 Re (Re = Ubar D / nu, Ubar = 2 U_m / 3) is the
only new input (PROTOCOL.md 2).  OpenFOAM runs through its input files only; an optional
command prefix (REUSE_OF_WRAP, shell syntax) puts it in a container or a wrapper.
"""
from __future__ import annotations

import math
import os
import re
import shlex
import shutil
import subprocess
import time

import numpy as np

NU, H, D, DEPTH = 1e-3, 0.41, 0.1, 0.01
MAX_ITER = 2000                 # the article's cap (grid runs)
BENCHMARK_CD_RE20 = 5.57953523384
_WRAP = shlex.split(os.environ.get("REUSE_OF_WRAP", ""))


def u_max(re_):
    """Inflow maximum of the parabolic profile for Re = Ubar D / nu, Ubar = 2 U_m / 3."""
    return 1.5 * re_ * NU / D


def of_cmd(cmd):
    return _WRAP + list(cmd)


HEADER = """FoamFile
{{
    version 2.0; format ascii; class {cls}; object {obj};
}}
"""


def _w(path, cls, obj, body):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(HEADER.format(cls=cls, obj=obj) + body)


def block_mesh_dict(na: int) -> str:
    """Twelve-block O-grid mesh, all sizes tied to na (cells per quarter arc)."""
    c, r = (0.2, 0.2), 0.05
    s = r / math.sqrt(2.0)
    nr = max(4, na // 2)                       # radial, graded to the wall
    nl = max(3, na // 2)                       # left of the square
    nb = max(3, na // 2)                       # below the square
    nt = max(3, (na * 11) // 20)               # above the square (0.11 tall)
    nx = 4 * na                                # downstream, expanding
    P = {
        "c0": (c[0] - s, c[1] - s), "c1": (c[0] + s, c[1] - s),
        "c2": (c[0] + s, c[1] + s), "c3": (c[0] - s, c[1] + s),
        "q0": (0.1, 0.1), "q1": (0.3, 0.1), "q2": (0.3, 0.3), "q3": (0.1, 0.3),
        "b0": (0.0, 0.0), "b1": (0.1, 0.0), "b2": (0.3, 0.0), "b3": (2.2, 0.0),
        "l1": (0.0, 0.1), "l2": (0.0, 0.3),
        "t0": (0.0, H), "t1": (0.1, H), "t2": (0.3, H), "t3": (2.2, H),
        "r1": (2.2, 0.1), "r2": (2.2, 0.3),
    }
    names = list(P)
    vid = {}
    verts = []
    for z in (0.0, DEPTH):
        for n in names:
            vid[(n, z)] = len(verts)
            verts.append((P[n][0], P[n][1], z))

    def hexa(a, b, cc, d, nx_, ny_, gx=1.0, gy=1.0):
        ids = [vid[(n, 0.0)] for n in (a, b, cc, d)] + \
              [vid[(n, DEPTH)] for n in (a, b, cc, d)]
        return (f"    hex ({' '.join(map(str, ids))}) ({nx_} {ny_} 1) "
                f"simpleGrading ({gx} {gy} 1)")

    gr = 0.25                                  # radial cells shrink 4x to the wall
    blocks = [
        hexa("q0", "q1", "c1", "c0", na, nr, 1, gr),
        hexa("q1", "q2", "c2", "c1", na, nr, 1, gr),
        hexa("q2", "q3", "c3", "c2", na, nr, 1, gr),
        hexa("q3", "q0", "c0", "c3", na, nr, 1, gr),
        hexa("b0", "b1", "q0", "l1", nl, nb),
        hexa("l1", "q0", "q3", "l2", nl, na),
        hexa("l2", "q3", "t1", "t0", nl, nt),
        hexa("b1", "b2", "q1", "q0", na, nb),
        hexa("q3", "q2", "t2", "t1", na, nt),
        hexa("b2", "b3", "r1", "q1", nx, nb, 6, 1),
        hexa("q1", "r1", "r2", "q2", nx, na, 6, 1),
        hexa("q2", "r2", "t3", "t2", nx, nt, 6, 1),
    ]

    def arc(a, b, ang):
        pt = (c[0] + r * math.cos(ang), c[1] + r * math.sin(ang))
        out = []
        for z in (0.0, DEPTH):
            out.append(f"    arc {vid[(a, z)]} {vid[(b, z)]} "
                       f"({pt[0]:.12f} {pt[1]:.12f} {z})")
        return out

    edges = (arc("c0", "c1", -math.pi / 2) + arc("c1", "c2", 0.0)
             + arc("c2", "c3", math.pi / 2) + arc("c3", "c0", math.pi))

    def face(a, b):          # side face between vertices a, b (both z levels)
        return (f"            ({vid[(a, 0.0)]} {vid[(b, 0.0)]} "
                f"{vid[(b, DEPTH)]} {vid[(a, DEPTH)]})")

    patches = {
        "inlet patch": [face("b0", "l1"), face("l1", "l2"), face("l2", "t0")],
        "outlet patch": [face("b3", "r1"), face("r1", "r2"), face("r2", "t3")],
        "walls wall": [face("b0", "b1"), face("b1", "b2"), face("b2", "b3"),
                       face("t0", "t1"), face("t1", "t2"), face("t2", "t3")],
        "cylinder wall": [face("c0", "c1"), face("c1", "c2"), face("c2", "c3"),
                          face("c3", "c0")],
    }
    body = "convertToMeters 1;\nvertices\n(\n"
    body += "\n".join(f"    ({x:.12f} {y:.12f} {z})" for x, y, z in verts)
    body += "\n);\nblocks\n(\n" + "\n".join(blocks) + "\n);\n"
    body += "edges\n(\n" + "\n".join(edges) + "\n);\n"
    body += "boundary\n(\n"
    for key, faces in patches.items():
        nm, typ = key.split()
        body += f"    {nm}\n    {{\n        type {typ};\n        faces\n        (\n"
        body += "\n".join(faces) + "\n        );\n    }\n"
    body += ("    frontAndBack\n    {\n        type empty;\n        faces\n        (\n")
    fb = []
    for blk in blocks:
        ids = list(map(int, re.search(r"hex \(([^)]*)\)", blk).group(1).split()))
        fb.append(f"            ({ids[0]} {ids[3]} {ids[2]} {ids[1]})")
        fb.append(f"            ({ids[4]} {ids[5]} {ids[6]} {ids[7]})")
    body += "\n".join(fb) + "\n        );\n    }\n);\nmergePatchPairs ();\n"
    return body


def control_dict(end=MAX_ITER):
    return f"""application simpleFoam;
startFrom latestTime; startTime 0; stopAt endTime; endTime {end};
deltaT 1; writeControl timeStep; writeInterval {end}; purgeWrite 0;
writeFormat ascii; writePrecision 16; writeCompression off;
timeFormat general; timePrecision 8; runTimeModifiable false;
"""


FV_SCHEMES = """ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes
{
    default none;
    div(phi,U) bounded Gauss linearUpwind grad(U);
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
"""


def fv_solution(res_ctrl, lin_tol, relax_u=0.9, n_coarsest=10):
    return f"""solvers
{{
    p {{ solver GAMG; smoother GaussSeidel; tolerance {lin_tol:g}; relTol 0.1;
        nCellsInCoarsestLevel {int(n_coarsest)}; }}
    U {{ solver smoothSolver; smoother symGaussSeidel; tolerance {10*lin_tol:g}; relTol 0.1; }}
}}
SIMPLE
{{
    nNonOrthogonalCorrectors 0; consistent yes;
    residualControl {{ p {res_ctrl:g}; U {res_ctrl/10:g}; }}
}}
relaxationFactors {{ equations {{ U {relax_u:g}; }} }}
"""


class OpenFOAMError(RuntimeError):
    """An OpenFOAM utility failed (an infrastructure error, not a failed setting)."""


def _sh(cmd, cwd):
    r = subprocess.run(of_cmd(cmd), cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        raise OpenFOAMError(f"{' '.join(cmd)} failed:\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
    return r.stdout


class Mesh:
    """One mesh level at one Reynolds number: blockMesh, checkMesh, the inlet profile and
    the wall geometry of the post-processing, built once and copied for every run."""

    def __init__(self, root, re_, na):
        self.re, self.na = float(re_), int(na)
        self.um = u_max(re_)
        self.ubar = 2.0 * self.um / 3.0
        d = self.dir = os.path.join(root, f"mesh_Re{re_:g}_n{na}")
        shutil.rmtree(d, ignore_errors=True)
        _w(f"{d}/system/blockMeshDict", "dictionary", "blockMeshDict", block_mesh_dict(na))
        _w(f"{d}/system/controlDict", "dictionary", "controlDict", control_dict())
        _w(f"{d}/system/fvSchemes", "dictionary", "fvSchemes", FV_SCHEMES)
        _w(f"{d}/system/fvSolution", "dictionary", "fvSolution", fv_solution(1e-3, 1e-6))
        _w(f"{d}/constant/transportProperties", "dictionary", "transportProperties",
           f"transportModel Newtonian;\nnu {NU};\n")
        _w(f"{d}/constant/turbulenceProperties", "dictionary", "turbulenceProperties",
           "simulationType laminar;\n")
        _sh(["blockMesh"], d)
        out = _sh(["checkMesh"], d)
        self.check_mesh = out
        self.cells = int(re.search(r"cells:\s+(\d+)", out).group(1))
        self.mesh_ok = "Mesh OK." in out
        _, cy = _patch_centres(d, "inlet")
        geo = self.geo = _wall_geometry(d, "cylinder")
        prof = 4.0 * self.um * cy * (H - cy) / H**2
        inlet = "nonuniform List<vector>\n%d\n(\n%s\n)" % (
            len(prof), "\n".join(f"({u:.15g} 0 0)" for u in prof))
        _w(f"{d}/0/U", "volVectorField", "U", f"""dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{{
    inlet {{ type fixedValue; value {inlet}; }}
    outlet {{ type zeroGradient; }}
    walls {{ type noSlip; }}
    cylinder {{ type noSlip; }}
    frontAndBack {{ type empty; }}
}}
""")
        _w(f"{d}/0/p", "volScalarField", "p", """dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
    inlet { type zeroGradient; }
    outlet { type fixedValue; value uniform 0; }
    walls { type zeroGradient; }
    cylinder { type zeroGradient; }
    frontAndBack { type empty; }
}
""")
        self.ang = np.arctan2(geo["cf"][:, 1] - 0.2, geo["cf"][:, 0] - 0.2)
        self.flow_rate = float(np.sum(prof * _patch_heights(d, "inlet")) * DEPTH)

    def remove(self):
        shutil.rmtree(self.dir, ignore_errors=True)


def run(mesh: Mesh, run_dir, res_ctrl, lin_tol, n_coarsest, max_iter=MAX_ITER):
    """One simpleFoam run on a built mesh.  Returns a dict: status 'ok' (finished,
    converged or at the cap) or 'failed' (simpleFoam stopped with an error or the
    result is not finite), the observables, the cost and the log facts."""
    d = run_dir
    shutil.rmtree(d, ignore_errors=True)
    shutil.copytree(mesh.dir, d)
    _w(f"{d}/system/controlDict", "dictionary", "controlDict", control_dict(max_iter))
    _w(f"{d}/system/fvSolution", "dictionary", "fvSolution",
       fv_solution(float(res_ctrl), float(lin_tol), 0.9, int(n_coarsest)))
    t0 = time.perf_counter()
    r = subprocess.run(of_cmd(["simpleFoam"]), cwd=d, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    log = r.stdout
    cpu = [float(x) for x in re.findall(r"ExecutionTime = ([\d.eE+-]+) s", log)]
    iters = len(re.findall(r"^Time = ", log, flags=re.M))
    out = dict(cpu=(cpu[-1] if cpu else wall), wall=wall, iterations=iters,
               cells=mesh.cells, cells_iters=mesh.cells * iters,
               converged="SIMPLE solution converged" in log, returncode=r.returncode,
               residuals_final=_last_residuals(log))
    nan = dict(cd=math.nan, cl=math.nan, dp=math.nan)
    if r.returncode != 0:
        out.update(status="failed", converged=False, error=(log + r.stderr)[-1500:], **nan)
        return out
    times = [t for t in os.listdir(d) if re.fullmatch(r"[\d.eE+-]+", t) and t != "0"]
    if not times:
        out.update(status="failed", converged=False, error="no time directory written", **nan)
        return out
    last = max(times, key=float)
    cd, cl = _force_coeffs(os.path.join(d, last), mesh.geo, mesh.ubar)
    pall = _internal_scalars(os.path.join(d, last, "p"))
    geo = mesh.geo
    pc = pall[geo["cells"]]
    # wall pressure extrapolated linearly along the wall normal from the two innermost
    # cells (second order), as in the article
    nh = geo["sf"] / np.linalg.norm(geo["sf"], axis=1)[:, None]
    dP = np.abs(np.einsum("ij,ij->i", geo["cc"] - geo["cf"], nh))
    dN = np.abs(np.einsum("ij,ij->i", geo["cc2"] - geo["cf"], nh))
    pw = pc + (pc - pall[geo["second"]]) * dP / (dN - dP)
    dp = _interp_angle(mesh.ang, pw, math.pi) - _interp_angle(mesh.ang, pw, 0.0)
    obs = dict(cd=cd, cl=cl, dp=dp)
    if not all(np.isfinite(list(obs.values()))):
        out.update(status="failed", converged=False, error="non-finite result", **nan)
        return out
    out.update(status="ok", end_time=float(last), **obs)
    return out


RESID = re.compile(r"Solving for (\w+), Initial residual = ([-+\d.eE]+)")


def _last_residuals(log):
    """Initial residual of each field in the last iteration of a log."""
    blocks = re.split(r"^Time = ", log, flags=re.M)
    if len(blocks) < 2:
        return {}
    return {m.group(1): float(m.group(2)) for m in RESID.finditer(blocks[-1])}


# ------------------------------------------------- the article's post-processing
def _foam_list(path):
    """Body of an ascii polyMesh list file: the lines between the outer ( )."""
    txt = open(path).read()
    txt = txt[txt.index("}", txt.index("FoamFile")) + 1:]
    m = re.search(r"^\s*(\d+)\s*\n\s*\(", txt, flags=re.M)
    n = int(m.group(1))
    return n, txt[m.end():]


def _points_faces(case):
    pm = os.path.join(case, "constant", "polyMesh")
    n, body = _foam_list(os.path.join(pm, "points"))
    pts = np.array(re.findall(r"\(([-\d.eE+]+) ([-\d.eE+]+) ([-\d.eE+]+)\)", body)[:n], float)
    n, body = _foam_list(os.path.join(pm, "faces"))
    faces = [list(map(int, f.split())) for f in re.findall(r"\d+\(([\d ]+)\)", body)[:n]]
    return pm, pts, faces


def _patch_range(pm, patch):
    btxt = open(os.path.join(pm, "boundary")).read()
    blk = re.search(r"\b%s\b\s*\{([^}]*)\}" % patch, btxt).group(1)
    nf = int(re.search(r"nFaces\s+(\d+)", blk).group(1))
    st = int(re.search(r"startFace\s+(\d+)", blk).group(1))
    return range(st, st + nf)


def _patch_centres(case, patch):
    """Face centres (x, y) of a boundary patch, from constant/polyMesh."""
    pm, pts, faces = _points_faces(case)
    c = np.array([pts[faces[i]].mean(axis=0) for i in _patch_range(pm, patch)])
    return c[:, 0], c[:, 1]


def _patch_heights(case, patch):
    """Height (y extent) of each face of a boundary patch in the x = const plane."""
    pm, pts, faces = _points_faces(case)
    return np.array([np.ptp(pts[faces[i]][:, 1]) for i in _patch_range(pm, patch)])


def _wall_geometry(case, patch):
    """Area vectors, centres and owner-cell centres of a wall patch's faces."""
    pm, pts, faces = _points_faces(case)
    n, body = _foam_list(os.path.join(pm, "owner"))
    owner = np.array(body.split(")")[0].split()[:n], int)
    fids = _patch_range(pm, patch)
    cells = owner[list(fids)]
    n2, body2 = _foam_list(os.path.join(pm, "neighbour"))
    neigh = np.array(body2.split(")")[0].split()[:n2], int)
    want = set(cells.tolist())
    cpts = {c: set() for c in want}
    for fi, f in enumerate(faces):
        o = owner[fi]
        if o in want:
            cpts[o].update(f)
        if fi < len(neigh) and neigh[fi] in want:
            cpts[neigh[fi]].update(f)
    cc = np.array([_prism_centroid(pts[sorted(cpts[c])]) for c in cells])
    sf, cf = [], []
    for i in fids:
        q = pts[faces[i]]
        sf.append(0.5 * np.cross(q[2] - q[0], q[3] - q[1]))
        cf.append(q.mean(axis=0))                  # rectangles: centroid = mean
    sf, cf = np.array(sf), np.array(cf)
    cell_faces = {c: [] for c in want}
    for fi, f in enumerate(faces[:len(neigh)]):
        if owner[fi] in want:
            cell_faces[owner[fi]].append(fi)
        if neigh[fi] in want:
            cell_faces[neigh[fi]].append(fi)
    second = []
    for k, i in enumerate(fids):
        c = cells[k]
        nh = sf[k] / np.linalg.norm(sf[k])
        best = max(cell_faces[c], key=lambda g: np.dot(pts[faces[g]].mean(axis=0) - cf[k], -nh))
        second.append(neigh[best] if owner[best] == c else owner[best])
    second = np.array(second)
    want2 = set(second.tolist()) - want
    cpts2 = {c: set() for c in want2}
    for fi, f in enumerate(faces):
        o = owner[fi]
        if o in want2:
            cpts2[o].update(f)
        if fi < len(neigh) and neigh[fi] in want2:
            cpts2[neigh[fi]].update(f)
    allc = dict(zip(cells.tolist(), cc))
    for c in want2:
        allc[c] = _prism_centroid(pts[sorted(cpts2[c])])
    cc2 = np.array([allc[c] for c in second])
    return dict(sf=sf, cf=cf, cc=cc, cells=cells, second=second, cc2=cc2)


def _prism_centroid(P):
    """Centroid of an extruded polygonal cell from its vertices (both z levels)."""
    z0 = P[:, 2].min()
    Q = P[np.abs(P[:, 2] - z0) < 1e-12][:, :2]
    m = Q.mean(axis=0)
    Q = Q[np.argsort(np.arctan2(Q[:, 1] - m[1], Q[:, 0] - m[0]))]
    x, y = Q[:, 0], Q[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cr = x * y1 - x1 * y
    A = cr.sum() / 2.0
    return np.array([((x + x1) * cr).sum() / (6 * A), ((y + y1) * cr).sum() / (6 * A),
                     0.5 * (P[:, 2].min() + P[:, 2].max())])


def _internal_vectors(path):
    txt = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<vector>\s*(\d+)\s*\(", txt)
    n = int(m.group(1))
    v = re.findall(r"\(([-\d.eE+]+) ([-\d.eE+]+) ([-\d.eE+]+)\)", txt[m.end():])[:n]
    return np.array(v, float)


def _internal_scalars(path):
    txt = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", txt)
    n = int(m.group(1))
    return np.array(txt[m.end():].split(")", 1)[0].split()[:n], float)


def _force_coeffs(tdir, geo, ubar):
    """Drag and lift coefficients from the written fields, as OpenFOAM's forces
    function object computes them: pressure times the face area vector plus the wall
    shear nu * U_t / d from the wall-adjacent cell (d from the cell centroid)."""
    U = _internal_vectors(os.path.join(tdir, "U"))[geo["cells"]]
    p = _internal_scalars(os.path.join(tdir, "p"))[geo["cells"]]   # zeroGradient wall
    sf, cf, cc = geo["sf"], geo["cf"], geo["cc"]
    mag = np.linalg.norm(sf, axis=1)
    nhat = sf / mag[:, None]                       # out of the fluid, into the body
    d = np.abs(np.einsum("ij,ij->i", cc - cf, nhat))
    ut = U - np.einsum("ij,ij->i", U, nhat)[:, None] * nhat
    F = (p[:, None] * sf).sum(axis=0) + (NU * mag[:, None] * ut / d[:, None]).sum(axis=0)
    q = 0.5 * ubar**2 * D * DEPTH
    return float(F[0] / q), float(F[1] / q)


def _interp_angle(ang, vals, target):
    """Linear interpolation in angle of patch face values at `target` (radians)."""
    d = np.angle(np.exp(1j * (ang - target)))
    lo = np.where(d <= 0, d, -np.inf).argmax()
    hi = np.where(d >= 0, d, np.inf).argmin()
    if lo == hi:
        return float(vals[lo])
    w = -d[lo] / (d[hi] - d[lo])
    return float(vals[lo] * (1 - w) + vals[hi] * w)


def openfoam_version():
    """The OpenFOAM release the wrapper or PATH gives (from simpleFoam -help's banner)."""
    try:
        r = subprocess.run(of_cmd(["simpleFoam", "-help"]), capture_output=True, text=True, timeout=120)
    except Exception as e:                       # noqa: BLE001 - recorded, not fatal
        return f"unknown ({e})"
    m = re.search(r"(?:Version|version|OpenFOAM-)\s*:?\s*v?(\d{4})", r.stdout + r.stderr)
    return f"v{m.group(1)}" if m else "unknown"
