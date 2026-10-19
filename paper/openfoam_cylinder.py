"""DFG benchmark "flow around a cylinder", case 2D-1 (steady, Re = 20), in OpenFOAM.

Schaefer & Turek (1996) define the problem; the reference values
    c_D = 5.57953523384,  c_L = 0.010618948146,  dp = 0.11752016697
are from the converged finite-element computations of John & Matthies (2001)
as tabulated in the FEATFLOW benchmark collection.  Because they are values of
the PDE solution, not of one discretisation, no converged OpenFOAM run is
needed: any consistent discretisation converges to them.

The solver is OpenFOAM's simpleFoam (SIMPLEC), driven only through its input
files, as a user would.  Four parameters are exposed:

    n_arc       cells per quarter of the cylinder; the whole block mesh scales
                with it                                        discretisation
    res_ctrl    SIMPLE residualControl for p (U: one tenth); the outer nonlinear
                iteration stops when the initial residuals fall below it
                                                               nested solver
    lin_tol     absolute tolerance of the GAMG pressure solver (U solver: ten
                times it); the inner linear solves stop there  nested solver
    n_coarsest  nCellsInCoarsestLevel of the GAMG agglomeration; changes
                the cost of each linear solve, not its answer  cost-only

A fifth input, relax_U (momentum under-relaxation, default 0.9), is documented
as affecting only the convergence rate; in OpenFOAM 1912 the converged answer
depends on it through the Rhie-Chow interpolation.  It is held at the tutorial
value in the main study and used in exp46 to test that a false "cost-only"
declaration is detected.

The OpenFOAM pitzDaily tutorial values are res_ctrl = 1e-2 (p) and 1e-3 (U),
GAMG tolerance 1e-6 with relTol 0.1 and the default nCellsInCoarsestLevel
of 10, smoothSolver tolerance 1e-5, relax_U = 0.9 with SIMPLEC, and
endTime 2000; the "tutorial" setting uses exactly these.

Cost is OpenFOAM's own ExecutionTime (CPU seconds) for the solve, excluding
meshing and post-processing.  A run that does not meet res_ctrl stops at
MAX_ITER iterations and reports what it has.
"""
from __future__ import annotations

import math
import os
import shlex
import re
import shutil
import subprocess
import tempfile
import time

import numpy as np

REF = dict(cd=5.57953523384, cl=0.010618948146, dp=0.11752016697)
MAX_ITER = 2000
NU, UM, H, D, DEPTH = 1e-3, 0.3, 0.41, 0.1, 0.01
UBAR = 2.0 * UM / 3.0
# The Ubuntu package needs WM_PROJECT_DIR; a sourced OpenFOAM environment
# (e.g. a cluster module) already sets it and is left alone.
ENV = dict(os.environ)
ENV.setdefault("WM_PROJECT_DIR", "/usr/share/openfoam")
# On a cluster OpenFOAM may live in a container: KNOBS_OF_WRAP is a command
# prefix (shell syntax) put in front of every OpenFOAM executable, e.g.
#   singularity exec of.sif bash -c 'source /usr/lib/openfoam/openfoam2406/etc/bashrc; exec "$@"' bash
_WRAP = shlex.split(os.environ.get("KNOBS_OF_WRAP", ""))


def of_cmd(cmd):
    """An OpenFOAM command line with the optional container prefix."""
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


def fv_solution(res_ctrl, lin_tol, relax_u, n_coarsest=10):
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


class DFGCylinder:
    def __init__(self, workdir=None, keep=False):
        self.root = workdir or tempfile.mkdtemp(prefix="knobs_dfg_")
        self.keep = keep
        self.meshes = {}

    def __del__(self):
        if not getattr(self, "keep", True):
            shutil.rmtree(getattr(self, "root", ""), ignore_errors=True)

    def _sh(self, cmd, cwd):
        r = subprocess.run(of_cmd(cmd), cwd=cwd, env=ENV, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"{' '.join(cmd)} failed:\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
        return r.stdout

    def mesh(self, na):
        """Build the mesh for one level once, with the inlet profile."""
        if na in self.meshes:
            return self.meshes[na]
        d = os.path.join(self.root, f"mesh_{na}")
        shutil.rmtree(d, ignore_errors=True)
        _w(f"{d}/system/blockMeshDict", "dictionary", "blockMeshDict", block_mesh_dict(na))
        _w(f"{d}/system/controlDict", "dictionary", "controlDict", control_dict())
        _w(f"{d}/system/fvSchemes", "dictionary", "fvSchemes", FV_SCHEMES)
        _w(f"{d}/system/fvSolution", "dictionary", "fvSolution", fv_solution(1e-3, 1e-6, 0.9))
        _w(f"{d}/constant/transportProperties", "dictionary", "transportProperties",
           f"transportModel Newtonian;\nnu {NU};\n")
        _w(f"{d}/constant/turbulenceProperties", "dictionary", "turbulenceProperties",
           "simulationType laminar;\n")
        self._sh(["blockMesh"], d)
        out = self._sh(["checkMesh"], d)
        ncell = int(re.search(r"cells:\s+(\d+)", out).group(1))
        # inlet profile from the inlet face centres
        cx_in, cy = _patch_centres(d, "inlet")
        geo = _wall_geometry(d, "cylinder")
        cyl_x, cyl_y = geo["cf"][:, 0], geo["cf"][:, 1]
        prof = 4.0 * UM * cy * (H - cy) / H**2
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
        ang = np.arctan2(cyl_y - 0.2, cyl_x - 0.2)
        self.meshes[na] = (d, ncell, ang, geo)
        return self.meshes[na]

    def run(self, theta: dict):
        """theta: n_arc, res_ctrl, lin_tol, n_coarsest[, relax_U] -> (obs, cost, wall, info)."""
        md, ncell, ang, geo = self.mesh(int(theta["n_arc"]))
        d = os.path.join(self.root, "run")
        shutil.rmtree(d, ignore_errors=True)
        shutil.copytree(md, d)
        _w(f"{d}/system/fvSolution", "dictionary", "fvSolution",
           fv_solution(float(theta["res_ctrl"]), float(theta["lin_tol"]),
                       float(theta.get("relax_U", 0.9)),
                       int(theta.get("n_coarsest", 10))))
        t0 = time.perf_counter()
        r = subprocess.run(of_cmd(["simpleFoam"]), cwd=d, env=ENV, capture_output=True, text=True)
        wall = time.perf_counter() - t0
        log = r.stdout
        nan = {k: float("nan") for k in ("cd", "cl", "dp", "dp_cell", "cd_vm", "cl_vm")}
        if r.returncode != 0:
            return nan, float("inf"), wall, {"unavailable": (r.stdout + r.stderr)[-300:]}
        cpu = [float(x) for x in re.findall(r"ExecutionTime = ([\d.eE+-]+) s", log)]
        iters = len(re.findall(r"^Time = ", log, flags=re.M))
        converged = "SIMPLE solution converged" in log
        times = [t for t in os.listdir(d) if re.fullmatch(r"[\d.eE+-]+", t) and t != "0"]
        last = max(times, key=float)
        cd, cl = _force_coeffs(os.path.join(d, last), geo)
        cd_vm, cl_vm = _force_coeffs(os.path.join(d, last), geo, centres="cc_vm")
        pall = _internal_scalars(os.path.join(d, last, "p"))
        pc = pall[geo["cells"]]
        dp_cell = _interp_angle(ang, pc, math.pi) - _interp_angle(ang, pc, 0.0)
        # wall pressure extrapolated linearly along the wall normal from the two
        # innermost cells (second order); the cell value alone is first order
        nh = geo["sf"] / np.linalg.norm(geo["sf"], axis=1)[:, None]
        dP = np.abs(np.einsum("ij,ij->i", geo["cc"] - geo["cf"], nh))
        dN = np.abs(np.einsum("ij,ij->i", geo["cc2"] - geo["cf"], nh))
        pw = pc + (pc - pall[geo["second"]]) * dP / (dN - dP)
        dp = _interp_angle(ang, pw, math.pi) - _interp_angle(ang, pw, 0.0)
        obs = {"cd": cd, "cl": cl, "dp": dp, "dp_cell": dp_cell, "cd_vm": cd_vm,
               "cl_vm": cl_vm}
        if not all(np.isfinite(list(obs.values()))):
            return nan, float("inf"), wall, {"unavailable": "non-finite result"}
        return obs, (cpu[-1] if cpu else wall), wall, dict(iters=iters, converged=converged,
                                                          cells=ncell)


def _foam_list(path):
    """Body of an ascii polyMesh list file: the lines between the outer ( )."""
    txt = open(path).read()
    txt = txt[txt.index("}", txt.index("FoamFile")) + 1:]
    m = re.search(r"^\s*(\d+)\s*\n\s*\(", txt, flags=re.M)
    n = int(m.group(1))
    return n, txt[m.end():]


def _patch_centres(case, patch):
    """Face centres (x, y) of a boundary patch, from constant/polyMesh."""
    pm = os.path.join(case, "constant", "polyMesh")
    n, body = _foam_list(os.path.join(pm, "points"))
    pts = np.array(re.findall(r"\(([-\d.eE+]+) ([-\d.eE+]+) ([-\d.eE+]+)\)", body)[:n], float)
    n, body = _foam_list(os.path.join(pm, "faces"))
    faces = [list(map(int, f.split())) for f in re.findall(r"\d+\(([\d ]+)\)", body)[:n]]
    btxt = open(os.path.join(pm, "boundary")).read()
    blk = re.search(r"\b%s\b\s*\{([^}]*)\}" % patch, btxt).group(1)
    nf = int(re.search(r"nFaces\s+(\d+)", blk).group(1))
    st = int(re.search(r"startFace\s+(\d+)", blk).group(1))
    c = np.array([pts[faces[i]].mean(axis=0) for i in range(st, st + nf)])
    return c[:, 0], c[:, 1]


def _wall_geometry(case, patch):
    """Area vectors, centres and owner-cell centres of a wall patch's faces."""
    pm = os.path.join(case, "constant", "polyMesh")
    n, body = _foam_list(os.path.join(pm, "points"))
    pts = np.array(re.findall(r"\(([-\d.eE+]+) ([-\d.eE+]+) ([-\d.eE+]+)\)", body)[:n], float)
    n, body = _foam_list(os.path.join(pm, "faces"))
    faces = [list(map(int, f.split())) for f in re.findall(r"\d+\(([\d ]+)\)", body)[:n]]
    n, body = _foam_list(os.path.join(pm, "owner"))
    owner = np.array(body.split(")")[0].split()[:n], int)
    btxt = open(os.path.join(pm, "boundary")).read()
    blk = re.search(r"\b%s\b\s*\{([^}]*)\}" % patch, btxt).group(1)
    nf = int(re.search(r"nFaces\s+(\d+)", blk).group(1))
    st = int(re.search(r"startFace\s+(\d+)", blk).group(1))
    fids = range(st, st + nf)
    cells = owner[list(fids)]
    # points of each wall cell: every face owned by or neighbouring it
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
    # Cell centres.  OpenFOAM uses the volume centroid; for these prismatic
    # (extruded) cells that is the centroid of the cell's polygon in the x-y
    # plane.  Version 1.1 of this script used the mean of the cell's vertices,
    # which differs from the centroid by a fixed fraction of the height of the
    # graded wall cells and biased the wall shear by a first-order term that
    # does not vanish relative to the discretisation error (it produced the
    # drag "plateau" at n_arc >= 48 of revision 2; exp59_dfg_plateau.py).
    cc_vm = np.array([pts[sorted(cpts[c])].mean(axis=0) for c in cells])
    cc = np.array([_prism_centroid(pts[sorted(cpts[c])]) for c in cells])
    sf, cf = [], []
    for i in fids:
        q = pts[faces[i]]
        sf.append(0.5 * np.cross(q[2] - q[0], q[3] - q[1]))
        cf.append(q.mean(axis=0))                  # rectangles: centroid = mean
    sf, cf = np.array(sf), np.array(cf)
    # second cell along the wall normal (for wall values of second order): the
    # neighbour of the wall cell across its face opposite the wall face
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
    return dict(sf=sf, cf=cf, cc=cc, cc_vm=cc_vm, cells=cells, second=second, cc2=cc2)


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
    v = re.findall(r"\(([-\d.eE+]+) ([-\d.eE+]+) ([-\d.eE+]+)\)", txt[m.end():], )[:n]
    return np.array(v, float)


def _internal_scalars(path):
    txt = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", txt)
    n = int(m.group(1))
    return np.array(txt[m.end():].split(")", 1)[0].split()[:n], float)


def _force_coeffs(tdir, geo, centres="cc"):
    """Drag and lift coefficients from the written fields, as OpenFOAM's forces
    function object computes them: pressure times the face area vector plus the
    wall shear nu * U_t / d from the wall-adjacent cell (d from the cell
    centroid; centres="cc_vm" reproduces the vertex-mean centres of v1.1)."""
    U = _internal_vectors(os.path.join(tdir, "U"))[geo["cells"]]
    p = _internal_scalars(os.path.join(tdir, "p"))[geo["cells"]]   # zeroGradient wall
    sf, cf, cc = geo["sf"], geo["cf"], geo[centres]
    mag = np.linalg.norm(sf, axis=1)
    nhat = sf / mag[:, None]                       # out of the fluid, into the body
    d = np.abs(np.einsum("ij,ij->i", cc - cf, nhat))
    ut = U - np.einsum("ij,ij->i", U, nhat)[:, None] * nhat
    F = (p[:, None] * sf).sum(axis=0) + (NU * mag[:, None] * ut / d[:, None]).sum(axis=0)
    q = 0.5 * UBAR**2 * D * DEPTH
    return float(F[0] / q), float(F[1] / q)


def _patch_list(path, patch):
    """Read a nonuniform (or uniform) patch value list from an ascii field."""
    txt = open(path).read()
    i = txt.index("boundaryField")
    j = re.search(r"\b%s\b\s*\{" % patch, txt[i:]).end() + i
    blk = txt[j:]
    m = re.search(r"value\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", blk)
    if m:
        n = int(m.group(1))
        vals = blk[m.end():].split(")", 1)[0].split()
        return np.array(vals[:n], float)
    m = re.search(r"value\s+uniform\s+([-\d.eE+]+)", blk)
    return np.array([float(m.group(1))])


def _interp_angle(ang, vals, target):
    """Linear interpolation in angle of patch face values at `target` (radians)."""
    d = np.angle(np.exp(1j * (ang - target)))
    lo = np.where(d <= 0, d, -np.inf).argmax()
    hi = np.where(d >= 0, d, np.inf).argmin()
    if lo == hi:
        return float(vals[lo])
    w = -d[lo] / (d[hi] - d[lo])
    return float(vals[lo] * (1 - w) + vals[hi] * w)


if __name__ == "__main__":
    p = DFGCylinder(keep=True)
    print("root", p.root)
    for na in [8, 16, 32]:
        th = dict(n_arc=na, res_ctrl=1e-6, lin_tol=1e-8, relax_U=0.9)
        obs, cost, wall, info = p.run(th)
        print(na, info, {k: f"{v:.6f} ({abs(v-REF[k])/abs(REF[k]):.1e})" for k, v in obs.items()},
              f"cpu {cost:.1f}s wall {wall:.1f}s")
