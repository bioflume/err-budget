"""Bring your own solver.

Everything below the line marked YOUR SOLVER is ordinary code with no knowledge
of `knobs`.  Replace it with your own and the rest of the file still works.

The solver here integrates a stiff scalar ODE with an implicit trapezoidal rule
and a fixed-point inner iteration, so it has the three parameter classes in
miniature: a step count (discretisation), an inner-iteration tolerance and cap
(nested solver), and how often the iteration matrix is refreshed (cost-only).
The exact solution is known, so you can check the advice.
"""
import numpy as np
from knobs import Knob, tune

LAM, T, U0 = -50.0, 1.0, 1.0          # u' = LAM*(u - cos t) - sin t


def exact(t):
    return np.cos(t) + (U0 - 1.0) * np.exp(LAM * t)


# ------------------------------------------------------------- YOUR SOLVER ---
def solver(steps, inner_tol, inner_max, jac_every):
    dt = T / steps
    u, t = U0, 0.0
    work = 0
    gamma = 1.0 / (1.0 - 0.5 * dt * LAM)      # the "expensive" iteration matrix
    for n in range(steps):
        f_old = LAM * (u - np.cos(t)) - np.sin(t)
        if n % jac_every == 0:
            gamma = 1.0 / (1.0 - 0.5 * dt * LAM)
            work += 8                          # refreshing it costs something
        v, tn = u, t + dt
        for _ in range(inner_max):
            f_new = LAM * (v - np.cos(tn)) - np.sin(tn)
            res = v - u - 0.5 * dt * (f_new + f_old)
            work += 1
            if abs(res) <= inner_tol * (abs(u) + 1e-30):
                break
            v -= gamma * res
        u, t = v, tn
    return {"u_final": u}, float(work)
# -----------------------------------------------------------------------------


knobs = [
    Knob("steps",     [25, 50, 100, 200, 400, 800],      kind="discretisation"),
    Knob("inner_tol", [1e-2, 1e-4, 1e-6, 1e-8, 1e-10],   kind="nested-solver"),
    Knob("inner_max", [2, 4, 8, 16],                     kind="nested-solver"),
    Knob("jac_every", [16, 8, 4, 2, 1],                  kind="cost-only"),
]

start = dict(steps=100, inner_tol=1e-10, inner_max=16, jac_every=1)
result = tune(solver, knobs, target=1e-4, observable="u_final", start=start)
print(result)

truth = exact(T)
for label, cfg in [("your setting", start), ("recommended", result.settings)]:
    out, cost = solver(**cfg)
    print(f"{label:>14}: error {abs(out['u_final'] - truth):.3e},  work {cost:.0f}")
