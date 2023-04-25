"""Utility functions."""

import numpy as np
from pathlib import Path


def get_path(path): return str(Path(path).resolve())


def normalize_cols(array):
    return array*(1.0/array.sum(1)[:, np.newaxis])
