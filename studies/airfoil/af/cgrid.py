"""Structured C-grid around a NACA section, computed point by point in Python and
written as an OpenFOAM polyMesh (PROTOCOL.md 3.1, amendments A1 and A2).

Indices.  i runs round the C: from the lower outlet along the lower side of the wake
cut to the trailing edge (TE), round the lower surface to the leading edge (LE), back
along the upper surface to the TE and along the upper side of the cut to the upper
outlet.  j runs outwards: 0 on the wall or the cut, n_in on the mid boundary, n_in +
n_out on the far field.

Inner layer (j <= n_in).  Straight rays leave every surface and cut point.  Ahead of
the split points (where the surfaces cross x = X_CENTRE) the rays start as the
leading-edge normal and turn smoothly to the vertical; behind them they turn to the
bisector of the TE and wake normals and stay parallel over the last TE_PARALLEL of the
surface and the first TE_PARALLEL of the wake, then turn to the vertical.  Around the
leading edge they are D_MID long; behind the split points they end on the lines
y = +-D_MID.  Along each ray the points grow geometrically from the first cell: dy1 on
the wall, and on the cut a height that grows from dy1 at the TE to the wake spacing.

Outer layer (j >= n_in).  Its rows are offset curves of the mid boundary at distances
growing geometrically to R - D_MID, so the far field is the line y = +-R behind the
split points and, ahead, a curve about R from the leading edge.  Its columns start
along the normals of the mid boundary and lean smoothly, row by row, to a far-field
distribution that is spread evenly along the far boundary.  So no column keeps the
trailing-edge or leading-edge clustering out to the far field: the cell size changes
smoothly everywhere (amendment A2: with straight-line block interpolation, narrow
columns over the rear met the fan of the front with a size jump growing with R, 364
at R = 500, and every run at R = 500 blew up).

Mesh level.  n_surf, the cells along each side of the section, sets the level;
r = n_surf / 64.  Every spacing scales with 1/r and every growth ratio g with g**(1/r),
except the first-cell height dy1 and the radius R, which are parameters of their own.

`check()` measures what checkMesh measures (cell centroids, every internal face,
the wake cut included) plus the size jump between neighbours; `write_polymesh()`
writes points, faces, owner, neighbour and boundary, one cell deep in z.
"""
from __future__ import annotations

import math
import os
from dataclasses import dataclass

import numpy as np

from .grading import cells_for_ratio, geometric, solve_ratio, vinokur
from .naca import Aerofoil

# ------------------------------------------------ constants of the mesh family
N_SURF0 = 64            # cells per side at level 0
D_LE0 = 2.0e-3          # surface spacing at the LE, level 0 (chords)
D_TE0 = 2.0e-3          # surface spacing at the TE, level 0
G_NORMAL0 = 1.25        # wall-normal growth ratio, level 0 (inner and outer layer)
G_AXIAL0 = 1.20         # streamwise growth ratio of the wake, level 0
D_WAKE0 = 4.0e-3        # first cell on the wake cut, away from the TE, level 0
D_MID = 1.5             # the mid boundary: rays D_MID long, lines y = +-D_MID
X_CENTRE = 0.5          # split points at x = X_CENTRE; the outlet at x = X_CENTRE + R
SPAN = 1.0              # depth in z (Aref = chord x span = 1)
TE_PARALLEL = 0.15      # parallel rays near the TE, on the surface and in the wake
FAR_EVEN = 0.85         # share of the even far-field distribution in the far columns
NORMAL_SMOOTH = 0.15    # smoothing of the mid-boundary normals (arc length, chords)
# Bumped whenever the mesh a setting produces changes; stored in every result, so that a
# result made with an older mesh is set aside and run again (PROTOCOL.md 14).
MESH_VERSION = 3


@dataclass(frozen=True)
class MeshSpec:
    n_surf: int         # cells per side of the section (mesh level)
    dy1: float          # first-cell height at the wall (chords)
    radius: float       # far-field radius (chords)

    @property
    def r(self) -> float:
        return self.n_surf / N_SURF0


def _unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def _wrap(d):
    return (d + math.pi) % (2 * math.pi) - math.pi


def _smooth(w):
    w = np.clip(w, 0.0, 1.0)
    return w * w * (3.0 - 2.0 * w)


def _arc(P):
    """Cumulative arc length along a polyline."""
    return np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])


def _at(P, cum, s):
    """Points at arc lengths s along the polyline P (cumulative lengths cum)."""
    return np.stack([np.interp(s, cum, P[:, 0]), np.interp(s, cum, P[:, 1])], axis=-1)


def _hull(P):
    """Convex hull of 2-D points, anticlockwise (Andrew's monotone chain)."""
    pts = sorted(set(map(tuple, np.round(P, 15))))

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower, upper = [], []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return np.array(lower[:-1] + upper[:-1])


def _ray_exit(P, D, H):
    """Distance along each ray P + t D (P inside the convex polygon H) to its boundary."""
    A, B = H, np.roll(H, -1, axis=0)
    E = B - A
    t_best = np.full(len(P), np.inf)
    for k in range(len(H)):
        den = D[:, 0] * E[k, 1] - D[:, 1] * E[k, 0]
        w = A[k] - P
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (w[:, 0] * E[k, 1] - w[:, 1] * E[k, 0]) / den
            u = (w[:, 0] * D[:, 1] - w[:, 1] * D[:, 0]) / den
        ok = (np.abs(den) > 1e-14) & (t > 0) & (u >= -1e-12) & (u <= 1 + 1e-12)
        t_best = np.where(ok & (t < t_best), t, t_best)
    return t_best


class CGrid:
    """The points of one mesh, G[i, j] (around x across x 2), and what it is made of."""

    def __init__(self, foil: Aerofoil, spec: MeshSpec):
        self.foil, self.spec = foil, spec
        self._build()

    # --------------------------------------------------------------- helpers
    def _split_s(self, side):
        """Arc length of the surface point at x = X_CENTRE (the first from the leading edge)."""
        f = self.foil
        S = f.side_length(side)
        s = np.linspace(0.0, S, 4001)
        g = f.at_s(side, s)[:, 0] - X_CENTRE
        k = int(np.argmax(g > 0))
        if k == 0 or g[k - 1] > 0:
            raise ValueError(f"no split point on the {side} surface")
        lo, hi = s[k - 1], s[k]
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            v = f.at_s(side, np.array([mid]))[0, 0] - X_CENTRE
            lo, hi = (mid, hi) if v < 0 else (lo, mid)
            if hi - lo < 1e-14:
                break
        return 0.5 * (lo + hi)

    # --------------------------------------------------------------- build
    def _build(self):
        f, sp = self.foil, self.spec
        r = sp.r
        g_n, g_ax = G_NORMAL0 ** (1.0 / r), G_AXIAL0 ** (1.0 / r)
        d_le, d_te = D_LE0 / r, D_TE0 / r
        R = float(sp.radius)
        x_m, x_out = 1.0 + D_MID, X_CENTRE + R
        if x_out < x_m + 2.0:
            raise ValueError(f"far-field radius {R} too small for the near wake")
        sgn = {"upper": 1.0, "lower": -1.0}
        vert = {"upper": math.pi / 2, "lower": -math.pi / 2}

        # ---- surface distributions, LE -> TE on each side, snapped to the split point
        surf = {}
        for side in ("upper", "lower"):
            S = f.side_length(side)
            s = vinokur(sp.n_surf, S, d_le, d_te)
            s_split = self._split_s(side)
            k = int(np.clip(np.argmin(np.abs(s - s_split)), 2, sp.n_surf - 2))
            front = s[:k + 1] * (s_split / s[k])
            rear = s_split + (s[k:] - s[k]) * (S - s_split) / (S - s[k])
            surf[side] = dict(S=S, s=np.concatenate([front, rear[1:]]), k=k, s_split=s_split)

        # ---- ray directions
        bis = {"upper": _unit(f.te_normal("upper") + np.array([0.0, 1.0])),
               "lower": _unit(f.te_normal("lower") + np.array([0.0, -1.0]))}
        front_table = {}
        for side in ("upper", "lower"):
            # the surface-normal angle made monotone (flat over concave stretches),
            # blended to the vertical by a smoothstep in arc length
            sd = surf[side]
            ss = f._s[side][f._s[side] < sd["s_split"]]
            ss = np.concatenate([ss, [sd["s_split"]]])
            th_n = np.unwrap(np.arctan2(*f.normal_at_s(side, ss)[:, ::-1].T))
            th_v = th_n[0] + _wrap(vert[side] - th_n[0])
            sg = 1.0 if th_v > th_n[0] else -1.0
            env = sg * np.maximum.accumulate(sg * th_n)
            front_table[side] = (ss, env + (th_v - env) * _smooth(ss / sd["s_split"]))

        def surf_dir(side, s):
            sd = surf[side]
            ss, tab = front_table[side]
            th_v = tab[-1]
            th_b = th_v + _wrap(math.atan2(bis[side][1], bis[side][0]) - th_v)
            rear = th_v + (th_b - th_v) * _smooth((s - sd["s_split"]) / (sd["S"] - TE_PARALLEL - sd["s_split"]))
            th = np.where(s <= sd["s_split"], np.interp(s, ss, tab), rear)
            return np.stack([np.cos(th), np.sin(th)], axis=-1)

        def wake_dir(side, x):
            th0 = math.atan2(bis[side][1], bis[side][0])
            th1 = th0 + _wrap(vert[side] - th0)
            th = th0 + (th1 - th0) * _smooth((x - 1.0 - TE_PARALLEL) / (D_MID - TE_PARALLEL))
            return np.stack([np.cos(th), np.sin(th)], axis=-1)

        def surf_len(side, s, d):
            """D_MID around the LE, blended to the split ray; the front mid curve meets the
            line y = +-D_MID with a horizontal tangent; behind: to that line."""
            sd = surf[side]
            P = f.at_s(side, s)
            yt = sgn[side] * D_MID
            Psp = f.at_s(side, np.array([sd["s_split"]]))[0]
            L_sp = abs(yt - Psp[1])
            w = _smooth(s / sd["s_split"])
            h = 1e-6
            slope = (f.at_s(side, np.array([sd["s_split"] + h]))[0, 1]
                     - f.at_s(side, np.array([sd["s_split"] - h]))[0, 1]) / (2 * h)
            u = s - sd["s_split"]
            fix = -sgn[side] * slope * u * np.exp(-(u / 0.1) ** 2)       # zero slope of the mid curve at the split
            to_line = np.where(sgn[side] * d[:, 1] > 1e-9, (yt - P[:, 1]) / np.where(d[:, 1] == 0, 1, d[:, 1]), np.inf)
            # never past the line y = +-D_MID: the mid boundary stays convex
            front = np.minimum(D_MID * (1 - w) + L_sp * w + fix, to_line)
            rear = (yt - P[:, 1]) / d[:, 1]
            return np.where(s <= sd["s_split"], front, rear)

        # ---- wake distribution along the cut
        d_te_eff = 0.5 * sum(np.diff(surf[s]["s"])[-1] for s in surf)
        n_nw = cells_for_ratio(D_MID, d_te_eff, g_ax)
        x_nw = 1.0 + geometric(n_nw, D_MID, d_te_eff)
        g_nw = solve_ratio(n_nw, D_MID, d_te_eff)
        first_fw = (x_nw[-1] - x_nw[-2]) * g_nw
        n_fw = cells_for_ratio(x_out - x_m, first_fw, g_ax)
        x_fw = x_m + geometric(n_fw, x_out - x_m, first_fw)
        x_cut = np.concatenate([x_nw, x_fw[1:]])                       # TE -> outlet

        # ---- columns in C order: base point, direction, length, first cell
        n_in = cells_for_ratio(D_MID, sp.dy1, g_n)
        d_wm = D_WAKE0 / r
        if d_wm * n_in >= D_MID:
            raise ValueError("wake-cut spacing too coarse for the inner layer")
        # the first cell on the cut grows from dy1 at the TE to d_wm, smoothly in log space
        # over the near wake (column by column it would tilt the rows above it)
        h_cut = sp.dy1 * (d_wm / sp.dy1) ** _smooth((x_cut - 1.0) / D_MID)
        cols = []

        def add(P, D, L, H):
            cols.append((np.asarray(P, float), np.asarray(D, float), np.asarray(L, float), np.asarray(H, float)))

        for side, order in (("lower", -1), ("upper", 1)):
            sd = surf[side]
            s = sd["s"] if order > 0 else sd["s"][::-1]
            d = surf_dir(side, s)
            surf_cols = (f.at_s(side, s), d, surf_len(side, s, d), np.full(len(s), sp.dy1))
            xc = x_cut if order > 0 else x_cut[::-1]
            hc = h_cut if order > 0 else h_cut[::-1]
            dw = wake_dir(side, xc)
            wake_cols = (np.stack([xc, np.zeros_like(xc)], axis=-1), dw, (sgn[side] * D_MID) / dw[:, 1], hc)
            if side == "lower":
                add(*wake_cols)
                add(*[a[1:] for a in surf_cols])
            else:
                add(*[a[1:] for a in surf_cols])
                add(*[a[1:] for a in wake_cols])
        P0 = np.concatenate([c[0] for c in cols])
        D0 = np.concatenate([c[1] for c in cols])
        L0 = np.concatenate([c[2] for c in cols])
        H0 = np.concatenate([c[3] for c in cols])
        n_lw = n_fw + n_nw
        n_a = len(P0) - 1
        if n_a != 2 * (n_fw + n_nw + sp.n_surf):
            raise AssertionError("column count")
        if np.any(L0 <= 0):
            raise ValueError("a ray of non-positive length")
        # The mid boundary must be convex (offset rows of a concave stretch cross a few
        # chords out): the front rays end on the convex hull of their ends and the lines
        # y = +-D_MID.  Where the ends are convex already nothing moves.
        i_sl = n_lw + (sp.n_surf - surf["lower"]["k"])
        i_su = n_lw + sp.n_surf + surf["upper"]["k"]
        fr = np.arange(i_sl, i_su + 1)
        ends = P0[fr] + L0[fr, None] * D0[fr]
        anchors = np.array([[X_CENTRE, D_MID], [X_CENTRE + 10.0, D_MID],
                            [X_CENTRE, -D_MID], [X_CENTRE + 10.0, -D_MID]])
        L0[fr] = _ray_exit(P0[fr], D0[fr], _hull(np.vstack([ends, anchors])))
        if not np.all(np.isfinite(L0[fr])):
            raise ValueError("a front ray misses the hull of the mid boundary")

        # ---- inner layer: geometric points along every ray
        G_in = np.empty((n_a + 1, n_in + 1, 2))
        last = np.empty(n_a + 1)
        for i in range(n_a + 1):
            x = geometric(n_in, L0[i], H0[i])
            G_in[i] = P0[i] + x[:, None] * D0[i]
            last[i] = x[-1] - x[-2]
        G_in[:, 0] = P0                                               # exactly on the wall and cut

        # ---- outer layer
        C0 = G_in[:, n_in]
        cum0 = _arc(C0)
        T = np.empty_like(C0)
        T[1:-1] = C0[2:] - C0[:-2]
        T[0], T[-1] = C0[1] - C0[0], C0[-1] - C0[-2]
        T = _unit(T)
        th = np.unwrap(np.arctan2(T[:, 1], T[:, 0]))
        # smooth the tangent angle in arc length (Gaussian); the straight ends stay straight
        w = np.exp(-0.5 * ((cum0[:, None] - cum0[None, :]) / NORMAL_SMOOTH) ** 2)
        th = (w @ th) / w.sum(axis=1)
        N = np.stack([-np.sin(th), np.cos(th)], axis=-1)              # outward: T turned anticlockwise
        N[0], N[-1] = (0.0, -1.0), (0.0, 1.0)
        i_c = int(np.argmin(np.abs(np.arange(n_a + 1) - (n_lw + sp.n_surf))))   # the LE column
        first_out = float(np.median(last[[n_lw, n_a - n_lw]])) * g_n
        n_out = cells_for_ratio(R - D_MID, first_out, g_n)
        delta = geometric(n_out, R - D_MID, first_out)
        Qn = C0[:, None, :] + delta[None, :, None] * N[:, None, :]    # normal lines
        # the even far distribution (rear regions even in x, front even in arc length)
        far = Qn[:, -1]
        far_even = far.copy()
        # rear regions: even in x between the outlet and the split point
        far_even[:i_sl + 1, 0] = np.linspace(far[0, 0], far[i_sl, 0], i_sl + 1)
        far_even[i_su:, 0] = np.linspace(far[i_su, 0], far[-1, 0], n_a - i_su + 1)
        far_even[:i_sl + 1, 1], far_even[i_su:, 1] = far[:i_sl + 1, 1], far[i_su:, 1]
        # front region: even in arc length along the far curve
        cf = _arc(far[i_sl:i_su + 1])
        far_even[i_sl:i_su + 1] = _at(far[i_sl:i_su + 1], cf, np.linspace(0, cf[-1], i_su - i_sl + 1))
        cum_far = _arc(far)
        # fractions of the far row's length: normal lines, and the target
        f_far_norm = cum_far / cum_far[-1]
        # position of the even points along the far row (they lie on it): project by x on
        # the straight parts and by arc length on the front
        s_even = np.empty(n_a + 1)
        s_even[:i_sl + 1] = np.interp(-far_even[:i_sl + 1, 0], -far[:i_sl + 1, 0], cum_far[:i_sl + 1])
        s_even[i_su:] = np.interp(far_even[i_su:, 0], far[i_su:, 0], cum_far[i_su:])
        s_even[i_sl:i_su + 1] = cum_far[i_sl] + np.linspace(0, cf[-1], i_su - i_sl + 1)
        f_far = (1 - FAR_EVEN) * f_far_norm + FAR_EVEN * s_even / cum_far[-1]
        if np.any(np.diff(f_far) <= 0):
            raise ValueError("far-field distribution not monotone")
        # the lean grows in proportion to the distance from the mid boundary: shifts of
        # order R along the far rows are spread over a height of order R
        beta = delta / delta[-1]
        G_out = np.empty((n_a + 1, n_out + 1, 2))
        for jj in range(n_out + 1):
            row = Qn[:, jj]
            cr = _arc(row)
            f_row = (1 - beta[jj]) * (cr / cr[-1]) + beta[jj] * f_far
            G_out[:, jj] = _at(row, cr, f_row * cr[-1])
        G_out[:, 0] = C0

        self.G = np.concatenate([G_in, G_out[:, 1:]], axis=1)
        self.n_in, self.n_out, self.n_lw, self.n_a = n_in, n_out, n_lw, n_a
        self.n_wall = 2 * sp.n_surf
        self.info = dict(n_surf=sp.n_surf, dy1=sp.dy1, radius=R, r=r, g_normal=g_n, g_axial=g_ax,
                         n_in=n_in, n_out=n_out, n_nw=n_nw, n_fw=n_fw,
                         n_front_upper=surf["upper"]["k"], n_front_lower=surf["lower"]["k"],
                         x_out=x_out, d_le=d_le, d_te=d_te, outer_first=first_out,
                         mesh_version=MESH_VERSION, cells=self.n_cells())

    # --------------------------------------------------------------- output
    def n_cells(self) -> int:
        return int(self.n_a * (self.n_in + self.n_out))

    def point_ids(self):
        """Point label of every (i, j), with the lower side of the wake cut merged into the
        upper side (the cut is internal)."""
        ni, nj = self.G.shape[:2]
        pid = np.arange(ni * nj).reshape(ni, nj)
        for i in range(self.n_lw + 1):
            pid[i, 0] = pid[self.n_a - i, 0]
        used = np.unique(pid)
        remap = -np.ones(ni * nj, dtype=np.int64)
        remap[used] = np.arange(len(used))
        return remap[pid], used

    def polymesh(self):
        """(points, faces, owner, neighbour, patches) of the mesh extruded one cell in z,
        in OpenFOAM's order: internal faces upper-triangular (by owner, then neighbour),
        then the patches.  Cell (i, j) has corners a = (i, j), b = (i+1, j), c = (i+1, j+1),
        d = (i, j+1) anticlockwise, so the outward normal of each edge lies to its right;
        a side face listed as (p, q, q', p') over the edge p -> q points that way."""
        na, nn = self.n_a, self.n_in + self.n_out
        pid, used = self.point_ids()
        npt = len(used)
        flat = self.G.reshape(-1, 2)[used]
        pts = np.concatenate([np.column_stack([flat, np.zeros(npt)]),
                              np.column_stack([flat, np.full(npt, SPAN)])])
        cell = np.arange(na)[:, None] * nn + np.arange(nn)[None, :]

        def side(p, q):
            p, q = p.ravel(), q.ravel()
            return np.column_stack([p, q, q + npt, p + npt])

        # internal: (i-1, j) | (i, j) over b -> c of the left cell; (i, j-1) | (i, j) over
        # c -> d of the lower cell; the wake cut: lower (i, 0) | upper (na-1-i, 0) over a -> b
        I, J = np.meshgrid(np.arange(1, na), np.arange(nn), indexing="ij")
        f1, o1, n1 = side(pid[I, J], pid[I, J + 1]), cell[I - 1, J].ravel(), cell[I, J].ravel()
        I, J = np.meshgrid(np.arange(na), np.arange(1, nn), indexing="ij")
        f2, o2, n2 = side(pid[I + 1, J], pid[I, J]), cell[I, J - 1].ravel(), cell[I, J].ravel()
        ic = np.arange(self.n_lw)
        f3, o3, n3 = side(pid[ic, 0], pid[ic + 1, 0]), cell[ic, 0], cell[na - 1 - ic, 0]
        fi = np.concatenate([f1, f2, f3])
        oi, ni = np.concatenate([o1, o2, o3]), np.concatenate([n1, n2, n3])
        if np.any(oi >= ni):
            raise AssertionError("owner must be the lower cell label")
        order = np.lexsort((ni, oi))
        faces, owner, neigh = [fi[order]], [oi[order]], ni[order]
        patches = []
        start = len(oi)
        iw = np.arange(self.n_lw, na - self.n_lw)                      # aerofoil: a -> b
        faces.append(side(pid[iw, 0], pid[iw + 1, 0]))
        owner.append(cell[iw, 0])
        patches.append(("airfoil", "wall", start, len(iw)))
        start += len(iw)
        ia, jn = np.arange(na), np.arange(nn)
        far = [side(pid[ia + 1, nn], pid[ia, nn]),                     # far field: c -> d
               side(pid[0, jn + 1], pid[0, jn]),                       # lower outlet: d -> a
               side(pid[na, jn], pid[na, jn + 1])]                     # upper outlet: b -> c
        faces += far
        owner += [cell[ia, nn - 1], cell[0, jn], cell[na - 1, jn]]
        nfar = sum(len(x) for x in far)
        patches.append(("farfield", "patch", start, nfar))
        start += nfar
        I, J = np.meshgrid(np.arange(na), np.arange(nn), indexing="ij")
        a, b, c, d = (pid[I, J].ravel(), pid[I + 1, J].ravel(), pid[I + 1, J + 1].ravel(),
                      pid[I, J + 1].ravel())
        fb = np.empty((2 * a.size, 4), dtype=np.int64)
        fb[0::2] = np.column_stack([a, d, c, b])                       # z = 0, normal -z
        fb[1::2] = np.column_stack([a, b, c, d]) + npt                 # z = SPAN, normal +z
        faces.append(fb)
        owner.append(np.repeat(cell.ravel(), 2))
        patches.append(("frontAndBack", "empty", start, len(fb)))
        return pts, np.concatenate(faces), np.concatenate(owner), neigh, patches

    def write_polymesh(self, case):
        """Write constant/polyMesh/{points,faces,owner,neighbour,boundary} (ascii)."""
        pts, faces, owner, neigh, patches = self.polymesh()
        d = os.path.join(case, "constant", "polyMesh")
        os.makedirs(d, exist_ok=True)
        ncell = self.n_cells()
        note = (f'"nPoints:{len(pts)}  nCells:{ncell}  nFaces:{len(faces)}  '
                f'nInternalFaces:{len(neigh)}"')

        def header(cls, obj, extra=""):
            return ("FoamFile\n{\n    version     2.0;\n    format      ascii;\n"
                    f"    class       {cls};\n{extra}    location    \"constant/polyMesh\";\n"
                    f"    object      {obj};\n}}\n\n")

        def write(name, cls, body_lines, n, extra=""):
            tmp = os.path.join(d, f".{name}.tmp")
            with open(tmp, "w") as fh:
                fh.write(header(cls, name, extra))
                fh.write(f"{n}\n(\n")
                fh.write(body_lines)
                fh.write(")\n")
            os.replace(tmp, os.path.join(d, name))

        write("points", "vectorField", "".join(f"({p[0]!r} {p[1]!r} {p[2]!r})\n"
                                               for p in pts.astype(float).tolist()), len(pts))
        write("faces", "faceList", "".join(f"4({a} {b} {c} {e})\n" for a, b, c, e in faces.tolist()), len(faces))
        write("owner", "labelList", "\n".join(map(str, owner.tolist())) + "\n", len(owner),
              f"    note        {note};\n")
        write("neighbour", "labelList", "\n".join(map(str, neigh.tolist())) + "\n", len(neigh),
              f"    note        {note};\n")
        body = ""
        for name, typ, start, n in patches:
            body += (f"    {name}\n    {{\n        type            {typ};\n"
                     + ("        inGroups        1(wall);\n" if typ == "wall" else "")
                     + f"        nFaces          {n};\n        startFace       {start};\n    }}\n")
        write("boundary", "polyBoundaryMesh", body, len(patches))

    # ---------------------------------------------------------- the check
    def centroids(self):
        G = self.G
        a, b, c, d = G[:-1, :-1], G[1:, :-1], G[1:, 1:], G[:-1, 1:]
        P = np.stack([a, b, c, d], axis=2)
        x, y = P[..., 0], P[..., 1]
        x1, y1 = np.roll(x, -1, axis=2), np.roll(y, -1, axis=2)
        cr = x * y1 - x1 * y
        A6 = 3.0 * cr.sum(axis=2)
        return np.stack([((x + x1) * cr).sum(axis=2) / A6, ((y + y1) * cr).sum(axis=2) / A6], axis=-1)

    def check(self) -> dict:
        return quality(self.G, self.n_lw, self.n_wall, self.centroids())


def quality(G, n_lw, n_wall, cen=None) -> dict:
    """What checkMesh measures on a structured C-grid G (cell centroids, every internal
    face including the wake cut), plus the largest size ratio of neighbouring cells, the
    first-cell height and the wall non-orthogonality."""
    a, b_, c, d = G[:-1, :-1], G[1:, :-1], G[1:, 1:], G[:-1, 1:]

    def cross(p, q):
        return p[..., 0] * q[..., 1] - p[..., 1] * q[..., 0]

    area = 0.5 * cross(c - a, d - b_)
    if cen is None:
        P = np.stack([a, b_, c, d], axis=2)
        x, y = P[..., 0], P[..., 1]
        x1, y1 = np.roll(x, -1, axis=2), np.roll(y, -1, axis=2)
        cr = x * y1 - x1 * y
        A6 = 3.0 * cr.sum(axis=2)
        cen = np.stack([((x + x1) * cr).sum(axis=2) / A6, ((y + y1) * cr).sum(axis=2) / A6], axis=-1)
    out = dict(cells=int(area.size), min_area=float(area.min()), inverted=int(np.sum(area <= 0)))
    ang, jumps = [], []
    for axis in (0, 1):
        if axis == 0:
            dc = cen[1:, :] - cen[:-1, :]
            f0, f1 = G[1:-1, :-1], G[1:-1, 1:]
            rr = area[1:] / area[:-1]
        else:
            dc = cen[:, 1:] - cen[:, :-1]
            f0, f1 = G[:-1, 1:-1], G[1:, 1:-1]
            rr = area[:, 1:] / area[:, :-1]
        t = f1 - f0
        nrm = np.stack([t[..., 1], -t[..., 0]], axis=-1)
        cosang = np.abs(np.sum(nrm * dc, axis=-1)) / (np.linalg.norm(nrm, axis=-1) * np.linalg.norm(dc, axis=-1))
        ang.append(np.degrees(np.arccos(np.clip(cosang, -1, 1))))
        jumps.append(np.maximum(rr, 1.0 / rr))
    n_a = G.shape[0] - 1
    low, up = cen[:n_lw, 0][::-1], cen[n_a - n_lw:, 0]
    t = G[n_a - n_lw + 1:, 0] - G[n_a - n_lw:-1, 0]
    nrm = np.stack([t[..., 1], -t[..., 0]], axis=-1)
    dc = up - low
    cosang = np.abs(np.sum(nrm * dc, axis=-1)) / (np.linalg.norm(nrm, axis=-1) * np.linalg.norm(dc, axis=-1))
    ang.append(np.degrees(np.arccos(np.clip(cosang, -1, 1))))
    al, au = area[:n_lw, 0][::-1], area[n_a - n_lw:, 0]
    jumps.append(np.maximum(al / au, au / al))
    out["max_non_orth_deg"] = float(max(x.max() for x in ang))
    out["faces_over_70_deg"] = int(sum((x > 70).sum() for x in ang))
    out["max_size_jump"] = float(max(x.max() for x in jumps))
    W0, W1 = G[n_lw:n_lw + n_wall + 1, 0], G[n_lw:n_lw + n_wall + 1, 1]
    face_t = W0[1:] - W0[:-1]
    face_len = np.linalg.norm(face_t, axis=-1)
    h = np.abs(area[n_lw:n_lw + n_wall, 0]) / face_len
    out["first_cell_min"], out["first_cell_max"] = float(h.min()), float(h.max())
    tn = face_t / face_len[:, None]
    rays = W1 - W0
    cosw = np.abs(cross(tn, rays[:-1])) / np.linalg.norm(rays[:-1], axis=-1)
    out["wall_non_orth_deg"] = float(np.degrees(np.arccos(np.clip(cosw, -1, 1))).max())
    ab, ad = np.linalg.norm(b_ - a, axis=-1), np.linalg.norm(d - a, axis=-1)
    out["aspect_max"] = float(np.max(np.maximum(ab, ad) / np.minimum(ab, ad)))
    out["wall_faces"] = int(n_wall)
    out["closed_te"] = bool(np.allclose(W0[0], W0[-1], atol=1e-12))
    return out


def build(foil: Aerofoil, spec: MeshSpec) -> CGrid:
    return CGrid(foil, spec)
