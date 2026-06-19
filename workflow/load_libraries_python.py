# =============================================================================
# load_libraries_python.py — Python dependencies for pseudotime analysis
# =============================================================================
# Loaded automatically by capitulo3_pseudotime.ipynb.
# Do not modify unless adding new packages.

import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.preprocessing import scale

import scanpy as sc
import scFates as scf
import palantir
from IPython import get_ipython

import matplotlib

shell = get_ipython()
if shell is not None and "IPKernelApp" in shell.config:
    shell.run_line_magic("matplotlib", "inline")
else:
    matplotlib.use("Agg")  # non-interactive backend for script execution

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

warnings.filterwarnings("ignore")
sc.settings.verbosity = 3
sc.settings.logfile = sys.stdout
