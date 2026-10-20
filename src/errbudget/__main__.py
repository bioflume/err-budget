"""`python -m errbudget` is the same as the `errbudget` command."""
import sys

from .cli import main

sys.exit(main())
