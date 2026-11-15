"""Where the replay scripts read the stored runs (data/) and write their outputs."""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = DATA
