"""Cells of a mesh and the predicted cost of a run (COST.md; unmeasured until Stage 0)."""
from __future__ import annotations

import functools

from . import study
from .cgrid import MeshSpec, build
from .naca import Aerofoil


@functools.lru_cache(maxsize=None)
def cells(n_surf: int, dy1: float, radius: float, design=study.CENTRE) -> int:
    """Cells of the C-grid; independent of the design to within a few cells."""
    return build(Aerofoil(*design, p=study.P_CAMBER, alpha_deg=study.ALPHA_DEG),
                 MeshSpec(int(n_surf), float(dy1), float(radius))).n_cells()


def predicted_core_seconds(run, t_cell_iter=None, iters_mid=None, growth=None) -> float:
    it = study.predicted_iters(run, iters_mid, growth)
    return cells(run.n_surf, run.dy1, run.radius) * it * (t_cell_iter or study.T_CELL_ITER)
