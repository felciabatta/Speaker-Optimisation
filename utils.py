"""Utility functions."""

from pathlib import Path


def get_path(path): return str(Path(path).resolve())
