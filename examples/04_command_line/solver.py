#!/usr/bin/env python3
"""A stand-in for your simulation code: a program that reads its numerical
parameters from the command line and prints a result.

It integrates the stiff ODE  u' = -50 (u - cos t) - sin t,  u(0) = 1,  with an
implicit trapezoidal rule and a fixed-point inner iteration, so it has the three
kinds of parameter in miniature:

    --steps     number of time steps                 (discretisation)
    --tol       inner-iteration tolerance            (nested solver)
    --refresh   how often the iteration matrix is rebuilt (cost-only)

errbudget never looks inside this file.  It only runs it and reads the two lines it
prints:  u_final = ...   and   work = ...
"""
import argparse
import math

ap = argparse.ArgumentParser()
ap.add_argument("--steps", type=int, required=True)
ap.add_argument("--tol", type=float, required=True)
ap.add_argument("--refresh", type=int, required=True)
ap.add_argument("--inner-max", type=int, default=16)
a = ap.parse_args()

LAM, T, U0 = -50.0, 1.0, 1.0
dt = T / a.steps
u, t, work = U0, 0.0, 0
gamma = 1.0 / (1.0 - 0.5 * dt * LAM)
for n in range(a.steps):
    f_old = LAM * (u - math.cos(t)) - math.sin(t)
    if n % a.refresh == 0:
        gamma = 1.0 / (1.0 - 0.5 * dt * LAM)
        work += 8
    v, tn = u, t + dt
    for _ in range(a.inner_max):
        f_new = LAM * (v - math.cos(tn)) - math.sin(tn)
        res = v - u - 0.5 * dt * (f_new + f_old)
        work += 1
        if abs(res) <= a.tol * (abs(u) + 1e-30):
            break
        v -= gamma * res
    u, t = v, tn

print(f"steps={a.steps} tol={a.tol} refresh={a.refresh}")
print(f"u_final = {u:.16g}")
print(f"work = {work}")
