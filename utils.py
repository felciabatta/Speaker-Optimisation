"""Utility functions."""

import numpy as np
from pathlib import Path


def get_path(path): return str(Path(path).resolve())


def normalize_rows(array):
    return array*(1.0/array.sum(1)[:, np.newaxis])


def distance_mat(f):
    f = np.array(f, dtype=np.float64)
    return abs(f-f[:, np.newaxis])
