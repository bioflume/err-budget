"""errbudget -- choose a solver's numerical parameters for a target accuracy.

Give it a function (or a shell command) that runs your solver, the numerical
parameters with their candidate values, the quantity you care about and how
accurate it must be.  You get a setting, an error estimate, and a statement of
which parameter limits the accuracy and which ones you are refining for nothing.

Python API:
    from errbudget import Parameter, tune
    result = tune(solver, [Parameter("N", [32, 64, 128, 256]), ...], target=1e-3, observable="drag")

Command line (shell-command solvers, no Python needed):
    errbudget init      # writes an annotated job file
    errbudget check job.yaml
    errbudget run job.yaml

The engine is version 1.2.2; ENGINE_VERSION records it.  The article (Computer
Physics Communications, 2026) describes version 1.2.1, which the configuration
"v1.2.1" (errbudget.V121) reproduces exactly.  The package version counts the
tooling around the engine.
"""
from .api import (Parameter, TuneResult, tune, KINDS, RunFailed, CONFIGURATIONS,
                  full_configuration)
from ._core import LocalBudget, measure_floor, run_protocol, screen, V11, V12, V121, V122
from .jobfile import JobSpec, load_job, run_job

__version__ = "2.1.0"
ENGINE_VERSION = "1.2.2"

__all__ = ["Parameter", "tune", "TuneResult", "KINDS", "RunFailed",
           "JobSpec", "load_job", "run_job", "CONFIGURATIONS", "full_configuration",
           "LocalBudget", "measure_floor", "run_protocol", "screen", "V11", "V12", "V121",
           "V122",
           "__version__", "ENGINE_VERSION"]
