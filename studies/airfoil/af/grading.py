"""One-dimensional point distributions and blockMesh gradings.

blockMesh places the points of an edge from a grading descriptor: a list of
sections (length fraction, cell fraction, expansion ratio), each section a
geometric progression whose last cell is `expansion` times its first.  The
positions are fractions of the edge's arc length.  `Grading.fractions()`
reproduces that rule exactly, so the geometric check of cgrid.py sees the
points blockMesh will make, and `fit_grading()` turns any target distribution
into a descriptor that hits the targets at every section boundary and to 2 % of a
cell in between.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


# ---------------------------------------------------------------- geometric
def solve_ratio(n: int, length: float, first: float) -> float:
    """Ratio g of a geometric progression of n cells, first cell `first`, total `length`."""
    if n < 1 or length <= 0 or first <= 0:
        raise ValueError(f"bad geometric progression n={n} length={length} first={first}")
    if n == 1:
        return 1.0
    target = length / first
    if abs(target - n) < 1e-12 * n:
        return 1.0

    def total(g):
        return n if abs(g - 1.0) < 1e-14 else (g**n - 1.0) / (g - 1.0)

    lo, hi = (1.0, 2.0) if target > n else (1e-6, 1.0)
    while target > n and total(hi) < target:
        hi *= 2.0
    for _ in range(200):                  # total(g) increases with g
        mid = 0.5 * (lo + hi)
        if total(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def geometric(n: int, length: float, first: float) -> np.ndarray:
    """Positions 0..length of n cells growing geometrically from `first`."""
    g = solve_ratio(n, length, first)
    d = first * g ** np.arange(n)
    x = np.concatenate([[0.0], np.cumsum(d)])
    return x * (length / x[-1])


def cells_for_ratio(length: float, first: float, g: float) -> int:
    """Smallest number of cells of ratio at most g that spans `length` from `first`."""
    if g <= 1.0:
        return max(1, int(math.ceil(length / first)))
    return max(1, int(math.ceil(math.log(1.0 + length * (g - 1.0) / first) / math.log(g))))


def vinokur(n: int, length: float, d0: float, d1: float) -> np.ndarray:
    """Two-sided stretching of Vinokur (J. Comput. Phys. 50, 1983): n cells on
    [0, length] with end spacings close to d0 and d1."""
    e0, e1 = d0 / length, d1 / length
    B = 1.0 / (n * math.sqrt(e0 * e1))
    A = math.sqrt(e1 / e0)
    xi = np.linspace(0.0, 1.0, n + 1)
    if abs(B - 1.0) < 1e-9:
        u = xi
    elif B > 1.0:
        lo, hi = 1e-9, 1.0
        while math.sinh(hi) / hi < B:
            hi *= 2.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if math.sinh(mid) / mid < B else (lo, mid)
        dy = 0.5 * (lo + hi)
        u = 0.5 * (1.0 + np.tanh(dy * (xi - 0.5)) / math.tanh(0.5 * dy))
    else:
        lo, hi = 1e-9, math.pi - 1e-9
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if math.sin(mid) / mid > B else (lo, mid)
        dy = 0.5 * (lo + hi)
        u = 0.5 * (1.0 + np.tan(dy * (xi - 0.5)) / math.tan(0.5 * dy))
    s = u / (A + (1.0 - A) * u)
    s[0], s[-1] = 0.0, 1.0
    return s * length


# ---------------------------------------------------------------- gradings
@dataclass
class Section:
    length: float      # fraction of the edge length
    cells: int         # number of cells
    expansion: float   # last cell / first cell


class Grading:
    """A blockMesh grading descriptor for one edge with n cells."""

    def __init__(self, sections):
        self.sections = [Section(float(s.length), int(s.cells), float(s.expansion)) for s in sections]
        tot = sum(s.length for s in self.sections)
        for s in self.sections:
            s.length /= tot
        if any(s.cells < 1 for s in self.sections):
            raise ValueError("every section needs at least one cell")

    @property
    def n(self) -> int:
        return sum(s.cells for s in self.sections)

    @classmethod
    def uniform(cls, n):
        return cls([Section(1.0, n, 1.0)])

    @classmethod
    def simple(cls, n, expansion):
        return cls([Section(1.0, n, expansion)])

    def reversed(self) -> "Grading":
        return Grading([Section(s.length, s.cells, 1.0 / s.expansion) for s in self.sections[::-1]])

    def fractions(self) -> np.ndarray:
        """Point positions in [0, 1] as blockMesh computes them."""
        out = [0.0]
        start = 0.0
        for s in self.sections:
            if s.cells == 1 or abs(s.expansion - 1.0) < 1e-15:
                d = np.full(s.cells, 1.0)
            else:
                q = s.expansion ** (1.0 / (s.cells - 1))
                d = q ** np.arange(s.cells)
            x = np.cumsum(d) / d.sum() * s.length
            out.extend((start + x).tolist())
            start += s.length
        f = np.array(out)
        f[-1] = 1.0
        return f

    def foam(self) -> str:
        """The grading as blockMesh syntax: a number, or a list of sections."""
        if len(self.sections) == 1:
            return repr(float(self.sections[0].expansion))
        n = self.n
        parts = [f"({float(s.length)!r} {s.cells / n!r} {float(s.expansion)!r})" for s in self.sections]
        return "(" + " ".join(parts) + ")"


def _section_points(f, b, e):
    """Interior points blockMesh places in one geometric section from f[b] to f[e]."""
    k = e - b
    exp = (f[e] - f[e - 1]) / (f[b + 1] - f[b]) if k > 1 else 1.0
    d = np.full(k, 1.0) if k == 1 or abs(exp - 1.0) < 1e-15 else (exp ** (1.0 / (k - 1))) ** np.arange(k)
    return f[b] + np.cumsum(d)[:-1] / d.sum() * (f[e] - f[b]), exp


def fit_grading(targets, tol: float = 0.02, max_cells: int = 8) -> Grading:
    """A multi-section grading through target positions (n+1 values, increasing).
    Sections are geometric, as blockMesh's are; each is extended while every point
    inside it stays within `tol` of the local spacing of its target, up to `max_cells`.
    So the targets are met exactly at section boundaries and to 2 % of a cell inside."""
    f = np.asarray(targets, float)
    f = (f - f[0]) / (f[-1] - f[0])
    n = len(f) - 1
    if n < 1 or np.any(np.diff(f) <= 0):
        raise ValueError("targets must increase strictly")
    secs, b = [], 0
    while b < n:
        best = b + 1
        for e in range(b + 2, min(n, b + max_cells) + 1):
            x, _ = _section_points(f, b, e)
            loc = np.minimum(f[b + 1:e] - f[b:e - 1], f[b + 2:e + 1] - f[b + 1:e])
            if np.max(np.abs(x - f[b + 1:e]) / loc) <= tol:
                best = e
            else:
                break
        _, exp = _section_points(f, b, best)
        secs.append(Section(f[best] - f[b], best - b, exp))
        b = best
    return Grading(secs)


def geometric_grading(n: int, length: float, first: float) -> Grading:
    """A one-section grading of n cells with first cell `first` on an edge of `length`."""
    g = solve_ratio(n, length, first)
    return Grading.simple(n, g ** (n - 1))
