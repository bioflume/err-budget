"""NACA 4-digit aerofoils with a closed (sharp) trailing edge.

The surface follows the standard NACA 4-digit formulas (Abbott & von Doenhoff,
"Theory of Wing Sections", 1959, Sec. 6.4), with the thickness applied normal to
the camber line.  The last thickness coefficient is -0.1036 instead of -0.1015, the
usual modification that closes the trailing edge: the five coefficients then sum to
zero, so y_t(1) = 0 and the trailing edge is a sharp corner at (1, 0).

The aerofoil is pitched nose-up by the angle of attack about its trailing edge, so
that the free stream is along +x, the wake cut of the C-grid is the x-axis behind
the trailing edge, and the trailing edge stays at (1, 0).

The surface is parametrised by one angle b in [-pi, pi]: x = (1 - cos|b|)/2, the
lower surface for b < 0 and the upper surface for b > 0.  b = -pi is the trailing
edge on the lower side, b = 0 the leading edge and b = pi the trailing edge on the
upper side.  The surface is smooth in b through the leading edge, so tangents and
normals there need no special case.
"""
from __future__ import annotations

import math

import numpy as np

# closed trailing edge: the standard last coefficient -0.1015 replaced by -0.1036
A_THICK = (0.2969, -0.1260, -0.3516, 0.2843, -0.1036)
A_THICK_OPEN = (0.2969, -0.1260, -0.3516, 0.2843, -0.1015)


def thickness(x, t, coeffs=A_THICK):
    """Half thickness y_t(x) of a NACA 4-digit section of thickness t (fraction of chord)."""
    x = np.asarray(x, float)
    a0, a1, a2, a3, a4 = coeffs
    return 5.0 * t * (a0 * np.sqrt(np.maximum(x, 0.0)) + a1 * x + a2 * x**2 + a3 * x**3 + a4 * x**4)


def camber(x, m, p):
    """Camber line y_c(x) and its slope for maximum camber m at chord fraction p."""
    x = np.asarray(x, float)
    if m == 0.0:
        return np.zeros_like(x), np.zeros_like(x)
    front = x < p
    yc = np.where(front, m / p**2 * (2 * p * x - x**2),
                  m / (1 - p)**2 * ((1 - 2 * p) + 2 * p * x - x**2))
    dyc = np.where(front, 2 * m / p**2 * (p - x), 2 * m / (1 - p)**2 * (p - x))
    return yc, dyc


class Aerofoil:
    """A NACA 4-digit section (m, p, t), closed trailing edge, pitched by alpha_deg."""

    def __init__(self, m: float, t: float, p: float = 0.4, alpha_deg: float = 4.0,
                 coeffs=A_THICK, nfine: int = 4001):
        if not (0.0 <= m < 0.2 and 0.0 < t < 0.4 and 0.05 < p < 0.95):
            raise ValueError(f"unsupported section m={m} p={p} t={t}")
        self.m, self.t, self.p, self.alpha = float(m), float(t), float(p), float(alpha_deg)
        self.coeffs = coeffs
        a = math.radians(self.alpha)
        # rotation by -alpha (clockwise) about the trailing edge: nose up
        self._rot = np.array([[math.cos(a), math.sin(a)], [-math.sin(a), math.cos(a)]])
        self.te = np.array([1.0, 0.0])
        # dense tables of the two sides (LE -> TE), used for arc length and inversion
        self._b = {s: np.linspace(0.0, math.pi, nfine) for s in ("upper", "lower")}
        self._pts = {s: self.point(self._b[s] * (1 if s == "upper" else -1)) for s in self._b}
        self._s = {}
        for s, P in self._pts.items():
            d = np.linalg.norm(np.diff(P, axis=0), axis=1)
            self._s[s] = np.concatenate([[0.0], np.cumsum(d)])

    # ------------------------------------------------------------------ points
    def unrotated(self, b):
        """Surface point(s) before pitching, for parameter b in [-pi, pi]."""
        b = np.atleast_1d(np.asarray(b, float))
        x = 0.5 * (1.0 - np.cos(np.abs(b)))
        yt = thickness(x, self.t, self.coeffs)
        yc, dyc = camber(x, self.m, self.p)
        th = np.arctan(dyc)
        sgn = np.where(b >= 0, 1.0, -1.0)
        X = x - sgn * yt * np.sin(th)
        Y = yc + sgn * yt * np.cos(th)
        return np.stack([X, Y], axis=-1)

    def rotate(self, P):
        P = np.asarray(P, float)
        return self.te + (P - self.te) @ self._rot.T

    def point(self, b):
        """Surface point(s) in the computational frame (pitched, free stream along +x)."""
        return self.rotate(self.unrotated(b))

    def tangent(self, b, h=1e-6):
        """Unit tangent in the direction of increasing b (TE lower -> LE -> TE upper)."""
        b = np.atleast_1d(np.asarray(b, float))
        lo = np.clip(b - h, -math.pi, math.pi)
        hi = np.clip(b + h, -math.pi, math.pi)
        T = self.point(hi) - self.point(lo)
        return T / np.linalg.norm(T, axis=-1, keepdims=True)

    def normal(self, b):
        """Unit outward normal.  Increasing b runs clockwise round the section, so the
        outward normal is the tangent turned 90 degrees anticlockwise."""
        T = self.tangent(b)
        return np.stack([-T[..., 1], T[..., 0]], axis=-1)

    @property
    def le(self):
        return self.point(0.0)[0]

    def quarter_chord(self):
        return self.rotate(np.array([0.25, 0.0]))

    # ------------------------------------------------------------- arc length
    def side_length(self, side):
        return float(self._s[side][-1])

    def b_of_s(self, side, s):
        """Parameter b (signed: negative on the lower side) at arc length s from the LE."""
        s = np.asarray(s, float)
        b = np.interp(s, self._s[side], self._b[side])
        return b if side == "upper" else -b

    def s_of_x(self, side, x):
        """Arc length from the LE to the surface point with chord parameter x."""
        b = math.acos(1.0 - 2.0 * x)
        return float(np.interp(b, self._b[side], self._s[side]))

    def at_s(self, side, s):
        return self.point(self.b_of_s(side, s))

    def normal_at_s(self, side, s):
        return self.normal(self.b_of_s(side, s))

    def te_normal(self, side):
        """Outward normal of one side at the trailing edge (the corner has two)."""
        b = math.pi - 1e-7 if side == "upper" else -math.pi + 1e-7
        return self.normal(b)[0]

    def contour(self, n=400):
        """Closed contour TE(lower) -> LE -> TE(upper), cosine spaced, for checks/plots."""
        b = np.linspace(-math.pi, math.pi, 2 * n + 1)
        return self.point(b)

    def area(self):
        P = self.contour(4000)
        x, y = P[:, 0], P[:, 1]
        return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))
