"""The design of the reuse study (PROTOCOL.md): members, parameter lists, runs, tasks."""
from __future__ import annotations

import itertools
from dataclasses import asdict, dataclass

# PROTOCOL.md 2: Re = 10, 12.5, ..., 40
MEMBERS = [10.0 + 2.5 * i for i in range(13)]

# PROTOCOL.md 3: the article's 2D-1 lists and classes
NAMES = ["n_arc", "res_ctrl", "lin_tol", "n_coarsest"]
LADDERS = [[8, 12, 16, 24, 32, 48, 64],
           [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7],
           [1e-4, 1e-5, 1e-6, 1e-7, 1e-8],
           [10, 50, 250, 1000]]
KINDS = ["discretisation", "nested-solver", "nested-solver", "cost-only"]
H_OF = {"n_arc": "1/value", "n_coarsest": "1/value"}       # the others: h = value
START_IDX = tuple(len(v) // 2 for v in LADDERS)            # (3, 3, 2, 2)
GRID = list(itertools.product(*[range(len(v)) for v in LADDERS]))     # 840

# PROTOCOL.md 4, 5: the references
REF_LEVELS = [96, 128, 192]
REF_SOLVER = dict(res_ctrl=1e-8, lin_tol=1e-10, n_coarsest=10)
MAX_ITER, MAX_ITER_REF = 2000, 20000

# PROTOCOL.md 6, 8
TARGETS = [0.01, 0.003]
TUTORIAL = dict(n_arc=32, res_ctrl=1e-2, lin_tol=1e-6, n_coarsest=10)
CAUTIOUS = dict(n_arc=64, res_ctrl=1e-7, lin_tol=1e-8, n_coarsest=10)
BENCHMARK = {20.0: 5.57953523384}


@dataclass(frozen=True)
class Run:
    re: float
    n_arc: int
    res_ctrl: float
    lin_tol: float
    n_coarsest: int
    ref: bool = False                 # a reference run (20,000 iterations at most)

    @property
    def tag(self) -> str:
        s = (f"Re{self.re:g}_n{self.n_arc}_rc{self.res_ctrl:.0e}_lt{self.lin_tol:.0e}"
             f"_nc{self.n_coarsest}").replace("+", "")
        return s + ("_ref" if self.ref else "")

    @property
    def max_iter(self) -> int:
        return MAX_ITER_REF if self.ref else MAX_ITER

    def setting(self) -> dict:
        return {n: getattr(self, n) for n in NAMES}

    def to_dict(self) -> dict:
        d = asdict(self)
        d["tag"] = self.tag
        return d


def setting_of(idx) -> dict:
    return {n: LADDERS[k][j] for k, (n, j) in enumerate(zip(NAMES, idx))}


def index_of(setting) -> tuple:
    return tuple(LADDERS[k].index(setting[n]) for k, n in enumerate(NAMES))


def run_at(re_, setting, ref=False) -> Run:
    s = dict(setting)
    return Run(re=float(re_), n_arc=int(s["n_arc"]), res_ctrl=float(s["res_ctrl"]),
               lin_tol=float(s["lin_tol"]), n_coarsest=int(s["n_coarsest"]), ref=ref)


def default_start() -> dict:
    return setting_of(START_IDX)


def grid_runs(re_):
    return [run_at(re_, setting_of(idx)) for idx in GRID]


def reference_runs(re_):
    return [run_at(re_, dict(REF_SOLVER, n_arc=n), ref=True) for n in REF_LEVELS]


# ------------------------------------------------------------------ tasks
@dataclass(frozen=True)
class Task:
    """One Slurm array element: the runs of one mesh at one member."""
    name: str
    re: float
    n_arc: int
    runs: tuple


def _grid_task(re_, na):
    runs = tuple(r for r in grid_runs(re_) if r.n_arc == na)
    return Task(f"grid_Re{re_:g}_n{na}", re_, na, runs)


def tasks_grid():
    """91 tasks: every member x every mesh level, the 120 solver settings of the mesh."""
    return [_grid_task(re_, na) for re_ in MEMBERS for na in LADDERS[0]]


def tasks_ref():
    """39 tasks: every member x every reference level, one run each."""
    return [Task(f"ref_Re{re_:g}_n{r.n_arc}", re_, r.n_arc, (r,))
            for re_ in MEMBERS for r in reference_runs(re_)]


CHECK_RUNS = (run_at(20.0, dict(n_arc=8, res_ctrl=1e-4, lin_tol=1e-6, n_coarsest=10)),
              run_at(20.0, dict(n_arc=12, res_ctrl=1e-4, lin_tol=1e-6, n_coarsest=10)),
              run_at(40.0, dict(n_arc=8, res_ctrl=1e-4, lin_tol=1e-6, n_coarsest=10)))


def tasks_check():
    """PROTOCOL.md 11: three grid runs, one task per mesh."""
    return [Task(f"check_Re{r.re:g}_n{r.n_arc}", r.re, r.n_arc, (r,)) for r in CHECK_RUNS]


SWEEPS = {"check": tasks_check, "grid": tasks_grid, "ref": tasks_ref}


def sweep(name):
    if name not in SWEEPS:
        raise KeyError(f"unknown sweep {name!r}; known: {', '.join(SWEEPS)}")
    return SWEEPS[name]()


def resource_class(task):
    """(time limit, memory) of a task's Slurm array."""
    if task.name.startswith("ref_"):
        return ("12:00:00", "6G") if task.n_arc >= 192 else ("06:00:00", "4G")
    if task.n_arc >= 48:
        return "12:00:00", "3G"
    return "04:00:00", "2G"
