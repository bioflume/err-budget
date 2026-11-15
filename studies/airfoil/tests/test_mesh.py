"""The C-grid (computed point by point) and its polyMesh writer (PROTOCOL.md 14, A2).

1-D distributions; the quality of every design's coarsest and finest meshes as checkMesh
measures it; the inner layer's independence of the far-field radius; the writer:
face count, owner < neighbour in OpenFOAM's order, outward normals on the boundary, closed
cells (the face area vectors of each cell sum to zero), positive volumes equal to the
grid's areas times the span, and the files read back; and the measure itself, which gives
checkMesh's own figures on the first generator's meshes from Hamilton.
"""
import importlib.util
import os
import re
import subprocess
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from af import study  # noqa: E402
from af.cgrid import D_MID, SPAN, MeshSpec, build, quality  # noqa: E402
from af.grading import (Grading, Section, fit_grading, geometric, geometric_grading,  # noqa: E402
                        solve_ratio, vinokur)
from af.naca import Aerofoil  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ---------------------------------------------------------------- 1-D
@pytest.mark.parametrize("n,L,d", [(40, 1.5, 1e-5), (163, 1.5, 1e-5), (10, 1.0, 0.1), (10, 1.0, 0.2)])
def test_geometric_progression(n, L, d):
    x = geometric(n, L, d)
    assert len(x) == n + 1 and x[0] == 0 and x[-1] == pytest.approx(L)
    assert x[1] == pytest.approx(d, rel=1e-9)
    r = np.diff(x)[1:] / np.diff(x)[:-1]
    assert np.allclose(r, solve_ratio(n, L, d))


def test_vinokur_end_spacings():
    s = vinokur(128, 1.03, 2e-3, 1e-3)
    assert s[0] == 0 and s[-1] == pytest.approx(1.03) and np.all(np.diff(s) > 0)
    assert np.diff(s)[0] == pytest.approx(2e-3, rel=0.1)
    assert np.diff(s)[-1] == pytest.approx(1e-3, rel=0.1)


def test_grading_fractions_and_fit():
    g = geometric_grading(20, 2.0, 0.01)
    assert g.fractions()[1] * 2.0 == pytest.approx(0.01, rel=1e-9)
    multi = Grading([Section(0.2, 3, 4.0), Section(0.8, 5, 0.5)])
    assert np.allclose(multi.reversed().reversed().fractions(), multi.fractions())
    t = np.sqrt(vinokur(181, 1.0, 1e-3, 5e-4))
    f = fit_grading(t).fractions()
    spacing = np.diff(t / t[-1])
    assert np.max(np.abs(f - t / t[-1])[1:-1] / np.minimum(spacing[:-1], spacing[1:])) < 0.05


# ---------------------------------------------------------------- the grid
@pytest.mark.parametrize("m,t", study.DESIGNS)
def test_every_design_coarsest_and_finest(m, t):
    f = Aerofoil(m, t)
    for spec in (MeshSpec(64, 4e-5, 50.0), MeshSpec(256, 1e-5, 5000.0)):
        g = build(f, spec)
        c = g.check()
        assert c["inverted"] == 0 and c["closed_te"]
        assert c["faces_over_70_deg"] == 0 and c["max_non_orth_deg"] < 60
        assert c["max_size_jump"] < 8                       # no size jump between neighbours
        assert c["wall_non_orth_deg"] < 18
        # rays are straight and graded from dy1: the mean height normal to the wall is dy1
        # cos(tilt), a hair more where rays fan out round the leading edge
        assert c["first_cell_max"] < 1.001 * spec.dy1
        assert c["first_cell_min"] > 0.95 * spec.dy1
        assert c["cells"] == g.n_cells() == g.n_a * (g.n_in + g.n_out)
        assert c["wall_faces"] == 2 * spec.n_surf


def test_reference_meshes():
    for m, t in ((0.06, 0.16), (0.0, 0.08)):
        for spec in (MeshSpec(362, 5e-6, 50000.0), MeshSpec(256, 1e-5, 500000.0)):
            c = build(Aerofoil(m, t), spec).check()
            assert c["inverted"] == 0 and c["faces_over_70_deg"] == 0 and c["max_size_jump"] < 8


def test_cells_double_per_level_and_far_field_is_logarithmic():
    f = Aerofoil(0.03, 0.12)
    n = [build(f, MeshSpec(ns, 2e-5, 500.0)).n_cells() for ns in study.LADDERS[0]]
    ratios = np.array(n[1:]) / np.array(n[:-1])
    assert np.all((ratios > 1.7) & (ratios < 2.2))
    r = [build(f, MeshSpec(128, 2e-5, R)).n_cells() for R in (50.0, 500.0, 5000.0, 50000.0)]
    inc = np.diff(r)
    assert np.all(inc > 0) and inc.max() < 0.45 * r[0]       # a tenfold radius adds rows, not x10


def test_inner_layer_does_not_depend_on_radius():
    f = Aerofoil(0.06, 0.16)
    g1, g2 = build(f, MeshSpec(128, 2e-5, 50.0)), build(f, MeshSpec(128, 2e-5, 5000.0))
    a, b = g1.n_lw - g1.info["n_nw"], g1.n_a - (g1.n_lw - g1.info["n_nw"])   # surface and near wake
    a2 = g2.n_lw - g2.info["n_nw"]
    assert g1.n_in == g2.n_in
    assert np.allclose(g1.G[a:b + 1, :g1.n_in + 1], g2.G[a2:a2 + (b - a) + 1, :g2.n_in + 1], atol=1e-12)


def test_first_cell_exactly_dy1_along_every_wall_ray():
    g = build(Aerofoil(0.03, 0.12), MeshSpec(128, 2e-5, 500.0))
    W = g.G[g.n_lw:g.n_lw + g.n_wall + 1]
    h = np.linalg.norm(W[:, 1] - W[:, 0], axis=1)
    assert np.allclose(h, 2e-5, rtol=1e-9)


def test_mid_boundary_convex_and_far_field_shape():
    for m, t in ((0.06, 0.08), (0.06, 0.16), (0.0, 0.08)):
        g = build(Aerofoil(m, t), MeshSpec(128, 2e-5, 500.0))
        C0 = g.G[:, g.n_in]
        T = np.diff(C0, axis=0)
        th = np.unwrap(np.arctan2(T[:, 1], T[:, 0]))
        assert np.all(np.diff(th) <= 1e-9)                    # turns one way only: convex
        far = g.G[:, -1]
        assert np.allclose(far[0], [g.info["x_out"], -500.0]) and np.allclose(far[-1], [g.info["x_out"], 500.0])
        assert abs(far[:, 1]).max() == pytest.approx(500.0) and far[:, 0].min() < -400


# ---------------------------------------------------------------- the writer
def face_vectors(pts, faces):
    p = pts[faces]                                            # (nf, 4, 3)
    S = 0.5 * np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 1])
    return S, p.mean(axis=1)


@pytest.mark.parametrize("m,t,spec", [(0.06, 0.16, MeshSpec(64, 2e-5, 50.0)),
                                      (0.0, 0.08, MeshSpec(91, 1e-5, 5000.0))])
def test_polymesh_topology_and_geometry(m, t, spec):
    g = build(Aerofoil(m, t), spec)
    pts, faces, owner, neigh, patches = g.polymesh()
    na, nn = g.n_a, g.n_in + g.n_out
    ncell = na * nn
    n_int = (na - 1) * nn + na * (nn - 1) + g.n_lw
    assert len(neigh) == n_int
    assert {p[0]: p[3] for p in patches} == {"airfoil": 2 * spec.n_surf, "farfield": na + 2 * nn,
                                             "frontAndBack": 2 * ncell}
    assert len(faces) == len(owner) == n_int + 2 * spec.n_surf + na + 2 * nn + 2 * ncell
    assert patches[0][2] == n_int and all(p[2] + p[3] == q[2] for p, q in zip(patches, patches[1:]))
    # points: unique, the wake cut merged
    assert len(np.unique(np.round(pts, 12), axis=0)) == len(pts)
    # owner < neighbour, upper-triangular order
    assert np.all(owner[:n_int] < neigh)
    key = owner[:n_int].astype(np.int64) * (ncell + 1) + neigh
    assert np.all(np.diff(key) > 0)
    # every cell has six faces
    count = np.bincount(owner, minlength=ncell) + np.bincount(neigh, minlength=ncell)
    assert np.all(count == 6)
    # closed cells: the area vectors (outward) of each cell's faces sum to zero
    S, cf = face_vectors(pts, faces)
    tot = np.zeros((ncell, 3))
    np.add.at(tot, owner, S)
    np.add.at(tot, neigh, -S[:n_int])
    scale = np.zeros(ncell)
    np.add.at(scale, owner, np.linalg.norm(S, axis=1))
    assert np.all(np.linalg.norm(tot, axis=1) <= 1e-9 * scale)
    # volumes (divergence theorem) positive and equal to the grid's areas times the span
    V = np.zeros(ncell)
    np.add.at(V, owner, np.einsum("ij,ij->i", S, cf) / 3)
    np.add.at(V, neigh, -np.einsum("ij,ij->i", S[:n_int], cf[:n_int]) / 3)
    G = g.G
    a, b, c, d = G[:-1, :-1], G[1:, :-1], G[1:, 1:], G[:-1, 1:]
    area = 0.5 * ((c - a)[..., 0] * (d - b)[..., 1] - (c - a)[..., 1] * (d - b)[..., 0])
    assert np.all(V > 0) and np.allclose(V, area.ravel() * SPAN, rtol=1e-7, atol=1e-18)
    # boundary faces point out of their cells
    cc = np.zeros((ncell, 3))
    cen = g.centroids()
    cc[:, :2] = cen.reshape(-1, 2)
    cc[:, 2] = SPAN / 2
    bf = slice(n_int, None)
    assert np.all(np.einsum("ij,ij->i", S[bf], cf[bf] - cc[owner[bf]]) > 0)
    # internal faces point from owner to neighbour
    assert np.all(np.einsum("ij,ij->i", S[:n_int], cc[neigh] - cc[owner[:n_int]]) > 0)
    # the aerofoil patch is the wall, and the wake cut is internal
    af = patches[0]
    wall = faces[af[2]:af[2] + af[3]]
    assert np.allclose(np.abs(pts[wall][:, :, 1]).max(), np.abs(g.G[g.n_lw:g.n_a - g.n_lw + 1, 0, 1]).max())


def read_foam_list(path):
    txt = open(path).read()
    body = txt[txt.index("}") + 1:]
    m = re.search(r"^(\d+)\s*\n\(", body, flags=re.M)
    return int(m.group(1)), body[m.end():body.rindex(")")]


def test_polymesh_files_read_back(tmp_path):
    g = build(Aerofoil(0.03, 0.12), MeshSpec(64, 2e-5, 50.0))
    g.write_polymesh(str(tmp_path))
    pts, faces, owner, neigh, patches = g.polymesh()
    d = tmp_path / "constant" / "polyMesh"
    n, body = read_foam_list(d / "points")
    P = np.array(re.findall(r"\(([-0-9.e+]+) ([-0-9.e+]+) ([-0-9.e+]+)\)", body), float)
    assert n == len(pts) and np.array_equal(P, pts)
    n, body = read_foam_list(d / "faces")
    F = np.array(re.findall(r"4\((\d+) (\d+) (\d+) (\d+)\)", body), int)
    assert n == len(faces) and np.array_equal(F, faces)
    for name, arr in (("owner", owner), ("neighbour", neigh)):
        n, body = read_foam_list(d / name)
        assert n == len(arr) and np.array_equal(np.array(body.split(), int), arr)
        head = open(d / name).read(400)
        assert f"nCells:{g.n_cells()}" in head and f"nInternalFaces:{len(neigh)}" in head
    b = open(d / "boundary").read()
    for name, typ, start, nf in patches:
        blk = re.search(r"%s\s*\{([^}]*)\}" % name, b).group(1)
        assert f"type            {typ};" in blk and f"nFaces          {nf};" in blk and f"startFace       {start};" in blk
    assert not any(f.startswith(".") for f in os.listdir(d))           # no temporary files left


# ---------------------------------------------------------------- the measure
def _module_at(commit, name):
    src = subprocess.run(["git", "-C", ROOT, "show", f"{commit}:airfoil/af/cgrid.py"],
                         capture_output=True, text=True, check=True).stdout
    spec_ = importlib.util.spec_from_loader(name, loader=None)
    mod = importlib.util.module_from_spec(spec_)
    mod.__package__ = "af"
    sys.modules[name] = mod                                   # dataclasses look their module up there
    try:
        exec(compile(src, f"{name}.py", "exec"), mod.__dict__)
    finally:
        sys.modules.pop(name, None)
    return mod


def test_measure_reproduces_check_mesh_on_hamilton():
    """checkMesh (v2606, Hamilton, 10 Nov 2026) on the first generator's two check meshes
    (commit 1062837): non-orthogonality 88.7021922011 and 88.4648135676, minimum volume
    4.17398418547e-08 and 4.15543014105e-08.  Rebuilt from that commit and interpolated as
    blockMesh does (commit 83e08a9 reproduces block::createPoints), `quality` gives them."""
    try:
        first = _module_at("1062837", "af.cgrid_first")
        emul = _module_at("83e08a9", "af.cgrid_emul")
    except Exception:
        pytest.skip("git history not available")
    for (m, t), non_orth, vmin in (((0.03, 0.12), 88.7021922011, 4.17398418547e-08),
                                   ((0.06, 0.16), 88.4648135676, 4.15543014105e-08)):
        grid = first.CGrid(Aerofoil(m, t), first.MeshSpec(64, 2e-5, 50.0))
        grid.__class__ = emul.CGrid
        G = grid.assemble()
        n_lw = sum(b.n1 for b in grid.blocks if b.layer == "inner" and b.name in ("fw_l_in", "nw_l_in"))
        c = quality(G, n_lw, sum(b.n1 for b in grid.blocks if b.wall))
        assert c["max_non_orth_deg"] == pytest.approx(non_orth, abs=0.01)
        assert c["min_area"] == pytest.approx(vmin, rel=1e-6)
