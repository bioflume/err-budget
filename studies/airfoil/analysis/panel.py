"""Inviscid lift of the study's sections: linear-strength vortex panels (Kuethe & Chow,
Foundations of Aerodynamics, 5th ed., 1998, §5.10), as a reference for the RANS lift.

At 4 degrees the method gives C_l = 0.468 for m 0, t 8 % (thin-aerofoil theory with
the thickness correction: 0.465), 0.871 for the centre design and 1.304 for m 6 %,
t 16 %.  A turbulent boundary layer at Re 3e6 lowers these by roughly a tenth while it
stays attached; much more than that means the flow has separated."""
from __future__ import annotations

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from af.naca import Aerofoil  # noqa: E402


def cl_panel(P, alpha=0.0):
    """C_l of the closed contour P (trailing edge, lower side, leading edge, upper side,
    trailing edge), free stream at angle alpha (radians) to x, chord 1."""
    X, Y = P[:, 0], P[:, 1]
    N = len(X) - 1
    xc, yc = 0.5 * (X[:-1] + X[1:]), 0.5 * (Y[:-1] + Y[1:])
    S = np.hypot(np.diff(X), np.diff(Y))
    th = np.arctan2(np.diff(Y), np.diff(X))
    CN1, CN2 = np.zeros((N, N)), np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i == j:
                CN1[i, j], CN2[i, j] = -1.0, 1.0
                continue
            A = -(xc[i] - X[j]) * math.cos(th[j]) - (yc[i] - Y[j]) * math.sin(th[j])
            B = (xc[i] - X[j]) ** 2 + (yc[i] - Y[j]) ** 2
            C = math.sin(th[i] - th[j])
            D = math.cos(th[i] - th[j])
            E = (xc[i] - X[j]) * math.sin(th[j]) - (yc[i] - Y[j]) * math.cos(th[j])
            F = math.log(1 + S[j] * (S[j] + 2 * A) / B)
            G = math.atan2(E * S[j], B + A * S[j])
            Q = (xc[i] - X[j]) * math.cos(th[i] - 2 * th[j]) - (yc[i] - Y[j]) * math.sin(th[i] - 2 * th[j])
            CN2[i, j] = D + 0.5 * Q * F / S[j] - (A * C + D * E) * G / S[j]
            CN1[i, j] = 0.5 * D * F + C * G - CN2[i, j]
    AN = np.zeros((N + 1, N + 1))
    rhs = np.zeros(N + 1)
    AN[:N, 0], AN[:N, N] = CN1[:, 0], CN2[:, N - 1]
    AN[:N, 1:N] = CN1[:, 1:] + CN2[:, :-1]
    rhs[:N] = np.sin(th - alpha)
    AN[N, 0] = AN[N, N] = 1.0                      # Kutta condition
    gam = 2 * math.pi * np.linalg.solve(AN, rhs)
    return 2.0 * float(np.sum(0.5 * (gam[:-1] + gam[1:]) * S))


def inviscid_cl(m, t, panels=160):
    """C_l of the section as the study pitches it (alpha = 4 degrees, free stream along x)."""
    return cl_panel(Aerofoil(m, t).contour(panels))


if __name__ == "__main__":
    for m, t in [(0.0, 0.08), (0.0, 0.12), (0.03, 0.12), (0.06, 0.08), (0.06, 0.16)]:
        print(f"m {m * 100:g} %, t {t * 100:g} %: inviscid C_l = {inviscid_cl(m, t):.3f}")
