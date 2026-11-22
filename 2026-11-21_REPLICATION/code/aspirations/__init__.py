"""Shared machinery for the empirical work of the paper: the data root, the output roots,
the data loader and the regression tools.

The paths are resolved from this file, so they are the same wherever a script is run from.
"""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")      # input data, distributed with the repository
OUT = os.path.join(ROOT, "out")        # text reports
ODATA = os.path.join(OUT, "data")      # results in machine-readable form
os.makedirs(ODATA, exist_ok=True)
