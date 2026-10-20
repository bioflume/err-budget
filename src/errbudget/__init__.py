"""errbudget -- choose a solver's numerical parameters for a target accuracy.

Public API:
    Parameter, tune, TuneResult
"""
from .api import Parameter, TuneResult, tune, KINDS, RunFailed
from ._core import LocalBudget, measure_floor, run_protocol, screen, V11, V12, V121

__version__ = "1.2.1"
__all__ = ["Parameter", "tune", "TuneResult", "KINDS", "RunFailed",
           "LocalBudget", "measure_floor", "run_protocol", "screen", "V11", "V12", "V121"]
