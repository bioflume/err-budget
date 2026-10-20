"""L1: 1D advection-diffusion, FD2 + Crank-Nicolson + GMRES."""
from __future__ import annotations
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


class AdvectionDiffusion:
    """u_t + a u_x = nu u_xx on [0,2pi) periodic.

    Space: 2nd-order central differences (p=2 in h=2pi/N)
    Time:  Crank-Nicolson (p=2 in dt)
    Solve: GMRES to relative tolerance tol (p=1 in tol, cost ~ log(1/tol))

    Observable: Q = int_0^{2pi} u(x,T)^2 dx, known exactly.
    Cost: a deterministic flop proxy (noise-free, reproducible), plus wall time.
    """

    def __init__(self, a=1.0, nu=0.02, T=1.0, modes=((1, 1.0), (3, 0.5)),
                 warm_start=True, precond=True):
        self.a, self.nu, self.T, self.modes = a, nu, T, modes
        self.warm_start, self.precond = warm_start, precond

    def exact_Q(self) -> float:
        return float(sum(np.pi * A**2 * np.exp(-2 * self.nu * m**2 * self.T)
                         for m, A in self.modes))

    def u0(self, x):
        return sum(A * np.sin(m * x) for m, A in self.modes)

    def operator(self, N):
        h = 2 * np.pi / N
        off = np.ones(N)
        # central advection
        Dx = sp.diags([off, -off], [1, -1], shape=(N, N), format="lil") / (2 * h)
        Dx[0, -1] = -1 / (2 * h)
        Dx[-1, 0] = 1 / (2 * h)
        # laplacian
        Lx = sp.diags([off, -2 * off, off], [-1, 0, 1], shape=(N, N), format="lil") / h**2
        Lx[0, -1] = 1 / h**2
        Lx[-1, 0] = 1 / h**2
        return (-self.a * Dx.tocsr() + self.nu * Lx.tocsr())

    def run(self, theta):
        import time as _time

        N, steps, tol = int(theta["N"]), int(theta["steps"]), float(theta["tol"])
        dt = self.T / steps
        x = np.arange(N) * (2 * np.pi / N)
        u = self.u0(x)
        A = self.operator(N)
        I = sp.identity(N, format="csr")
        M = (I - 0.5 * dt * A).tocsr()
        R = (I + 0.5 * dt * A).tocsr()

        total_iters = 0
        t0 = _time.perf_counter()
        for _ in range(steps):
            b = R @ u
            it = [0]

            def cb(_):
                it[0] += 1

            x0 = u if self.warm_start else np.zeros(N)
            try:
                u, info = spla.gmres(M, b, x0=x0, rtol=tol, atol=0.0,
                                     restart=min(50, N), maxiter=400,
                                     callback=cb, callback_type="pr_norm")
            except TypeError:  # older scipy
                u, info = spla.gmres(M, b, x0=x0, tol=tol, atol=0.0,
                                     restart=min(50, N), maxiter=400,
                                     callback=cb, callback_type="pr_norm")
            total_iters += max(it[0], 1)
        wall = _time.perf_counter() - t0

        Q = float(np.sum(u**2) * (2 * np.pi / N))       # spectrally accurate on periodic grid
        cost = float(N * (total_iters + 2 * steps))      # flop proxy
        return Q, cost, wall

    def runner(self):
        def f(theta):
            Q, cost, _ = self.run(theta)
            return Q, cost

        return f
