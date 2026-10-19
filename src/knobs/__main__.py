"""`python -m knobs` is the same as the `knobs` command."""
import sys

from .cli import main

sys.exit(main())
