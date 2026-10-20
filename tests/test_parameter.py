import pytest
from errbudget import Parameter

def test_ordering_inferred_for_counts_and_tolerances():
    assert Parameter("N", [32, 64, 128]).h[0] > Parameter("N", [32, 64, 128]).h[-1]
    assert Parameter("tol", [1e-2, 1e-6]).h == [1e-2, 1e-6]

def test_rejects_wrong_direction():
    with pytest.raises(ValueError, match="cannot infer the refinement direction"):
        Parameter("tol", [1e-8, 1e-2])

def test_rejects_bad_kind_and_short_ladder():
    with pytest.raises(ValueError, match="kind must be"):
        Parameter("x", [1, 2], kind="nonsense")
    with pytest.raises(ValueError, match="at least two"):
        Parameter("x", [1])

def test_custom_h():
    k = Parameter("cells", [10, 100], h=lambda v: 1.0 / v)
    assert k.h == [0.1, 0.01]
