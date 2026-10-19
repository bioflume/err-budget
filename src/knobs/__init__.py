"""knobs -- choose a solver's numerical parameters for a target accuracy.

Give it a function (or a shell command) that runs your solver, the numerical
parameters with their candidate values, the quantity you care about and how
accurate it must be.  You get a setting, an error estimate, and a statement of
which parameter limits the accuracy and which ones you are refining for nothing.

Python API:
    from knobs import Knob, tune
    result = tune(solver, [Knob("N", [32, 64, 128, 256]), ...], target=1e-3, observable="drag")

Command line (shell-command solvers, no Python needed):
    knobs init          # writes an annotated job file
    knobs check job.yaml
    knobs run job.yaml

The engine is version 1.2.1 of the article (Computer Physics Communications,
2026); ENGINE_VERSION records it.  The package version counts the tooling
around it.
"""
from .api import Knob, TuneResult, tune, KINDS, RunFailed
from ._core import LocalBudget, measure_floor, run_protocol, screen, V11, V12, V121
from .jobfile import JobSpec, load_job, run_job

__version__ = "1.3.0"
ENGINE_VERSION = "1.2.1"
CONFIGURATIONS = {"v1.2.1": V121, "v1.2": V12, "v1.1": V11}

__all__ = ["Knob", "tune", "TuneResult", "KINDS", "RunFailed",
           "JobSpec", "load_job", "run_job", "CONFIGURATIONS",
           "LocalBudget", "measure_floor", "run_protocol", "screen", "V11", "V12", "V121",
           "__version__", "ENGINE_VERSION"]
