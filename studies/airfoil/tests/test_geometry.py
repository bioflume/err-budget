"""NACA 4-digit geometry: closed trailing edge, frame, normals, arc length."""
import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from af import study  # noqa: E402
from af.naca import A_THICK, A_THICK_OPEN, Aerofoil, camber, thickness  # noqa: E402


def test_closed_trailing_edge_coefficients():
    assert thickness(1.0, 0.12) == pytest.approx(0.0, abs=1e-15)
    assert sum(A_THICK) == pytest.approx(0.0, abs=1e-15)
    # the standard section is blunt: half thickness 0.0105 t at x = 1
    assert thickness(1.0, 0.12, A_THICK_OPEN) == pytest.approx(0.0105 * 0.12, rel=1e-9)


def test_camber_line():
    yc, dyc = camber(np.array([0.0, 0.4, 1.0]), 0.06, 0.4)
    assert yc == pytest.approx([0.0, 0.06, 0.0], abs=1e-15)
    assert dyc[1] == pytest.approx(0.0, abs=1e-15)
    yc0, _ = camber(np.linspace(0, 1, 11), 0.0, 0.4)
    assert np.all(yc0 == 0.0)


@pytest.mark.parametrize("m,t", study.DESIGNS)
def test_contour_closed_and_chord(m, t):
    f = Aerofoil(m, t)
    lo, hi = f.point(-math.pi)[0], f.point(math.pi)[0]
    assert np.allclose(lo, [1.0, 0.0], atol=1e-14) and np.allclose(hi, [1.0, 0.0], atol=1e-14)
    a = math.radians(4.0)
    assert np.allclose(f.le, [1 - math.cos(a), math.sin(a)], atol=1e-14)
    assert np.linalg.norm(f.te - f.le) == pytest.approx(1.0, abs=1e-14)
    # area of a NACA 4-digit section ~ 0.685 t c (camber changes it by < 0.5 %)
    assert abs(f.area()) == pytest.approx(0.6851 * t, rel=0.01)
    # the contour runs clockwise (negative signed area) and never crosses itself:
    # upper side above lower side at every chord station
    x = np.linspace(0.001, 0.999, 200)
    b = np.arccos(1 - 2 * x)
    up = f.unrotated(b)
    lo_ = f.unrotated(-b)
    assert f.area() < 0 and np.all(up[:, 1] > lo_[:, 1])


def test_symmetric_section_mirror_and_max_thickness():
    f = Aerofoil(0.0, 0.12, alpha_deg=0.0)
    b = np.linspace(0.01, math.pi - 0.01, 50)
    up, lo = f.point(b), f.point(-b)
    assert np.allclose(up[:, 0], lo[:, 0]) and np.allclose(up[:, 1], -lo[:, 1])
    bb = np.linspace(0, math.pi, 20001)
    assert 2 * f.point(bb)[:, 1].max() == pytest.approx(0.12, rel=2e-3)


@pytest.mark.parametrize("m,t", [(0.0, 0.08), (0.06, 0.16), (0.03, 0.12), (0.06, 0.08)])
def test_normals_unit_outward_and_smooth(m, t):
    f = Aerofoil(m, t)
    b = np.linspace(-math.pi + 0.05, math.pi - 0.05, 400)
    n = f.normal(b)
    assert np.allclose(np.linalg.norm(n, axis=1), 1.0)
    # a small step along the normal leaves the section, a step against it enters it
    from matplotlib.path import Path
    poly = Path(f.contour(2000))
    P = f.point(b)
    assert not poly.contains_points(P + 1e-5 * n).any()
    assert poly.contains_points(P - 1e-5 * n).all()
    # smooth through the leading edge: neighbouring normals differ by a small angle
    nb = f.normal(np.linspace(-0.02, 0.02, 41))
    ang = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", nb[1:], nb[:-1]), -1, 1)))
    assert ang.max() < 10.0
    # leading-edge normal points upstream
    assert f.normal(0.0)[0][0] < -0.9


@pytest.mark.parametrize("side", ["upper", "lower"])
def test_arc_length_inverse(side):
    f = Aerofoil(0.045, 0.10)
    S = f.side_length(side)
    assert 1.0 < S < 1.1
    s = np.linspace(0, S, 37)
    P = f.at_s(side, s)
    d = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])
    assert np.allclose(d, s, rtol=0, atol=2e-3)          # chords of 37 points vs arc length
    assert f.s_of_x(side, 0.3) == pytest.approx(
        float(np.interp(0.3, (1 - np.cos(np.abs(f.b_of_s(side, s)))) / 2, s)), abs=2e-3)
