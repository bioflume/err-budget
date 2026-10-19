"""Make the tests runnable from a fresh clone without installing the package:
if `knobs` is not importable, add ../src to the path."""
import os
import sys

try:
    import knobs  # noqa: F401
except ModuleNotFoundError:
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
