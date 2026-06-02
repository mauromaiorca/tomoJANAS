#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/metadata/imod_xf.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
Parser and helpers for IMOD ``.xf`` 2D transform files.

Each ``.xf`` line has six numbers::

    A11 A12 A21 A22 DX DY

mapping raw -> aligned by default::

    X_ali = A11 * X_raw + A12 * Y_raw + DX
    Y_ali = A21 * X_raw + A22 * Y_raw + DY

The inverse maps aligned -> raw and is used to place a particle (picked in the
aligned tomogram space) back into the raw tilt frames.
"""

from __future__ import annotations

from typing import Mapping, Tuple, Union

import numpy as np
import pandas as pd

__all__ = [
    "read_xf",
    "xf_matrix",
    "xf_inverse_matrix",
    "apply_xf_to_point",
    "invert_xf_for_point",
]

_FIELDS = ["a11", "a12", "a21", "a22", "dx", "dy"]


def read_xf(path: str) -> pd.DataFrame:
    """Parse an ``.xf`` file into a DataFrame.

    Columns: ``tilt_index`` (1-based), ``a11, a12, a21, a22, dx, dy``.
    """
    rows = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            toks = line.split()
            if len(toks) < 6:
                continue
            vals = [float(t) for t in toks[:6]]
            rows.append(vals)
    df = pd.DataFrame(rows, columns=_FIELDS)
    df.insert(0, "tilt_index", list(range(1, len(df) + 1)))
    return df


def _row_vals(row: Union[Mapping, pd.Series]) -> Tuple[float, float, float, float, float, float]:
    g = (lambda k: float(row[k]))
    return g("a11"), g("a12"), g("a21"), g("a22"), g("dx"), g("dy")


def xf_matrix(row: Union[Mapping, pd.Series]) -> np.ndarray:
    """Return the 3x3 homogeneous affine (raw -> aligned) for one xf row."""
    a11, a12, a21, a22, dx, dy = _row_vals(row)
    return np.array(
        [[a11, a12, dx],
         [a21, a22, dy],
         [0.0, 0.0, 1.0]],
        dtype=np.float64,
    )


def xf_inverse_matrix(row: Union[Mapping, pd.Series]) -> np.ndarray:
    """Return the 3x3 homogeneous affine (aligned -> raw) for one xf row."""
    return np.linalg.inv(xf_matrix(row))


def apply_xf_to_point(x: float, y: float, row: Union[Mapping, pd.Series]) -> Tuple[float, float]:
    """Map a raw-stack point (x, y) to aligned-stack coordinates."""
    a11, a12, a21, a22, dx, dy = _row_vals(row)
    return (a11 * x + a12 * y + dx, a21 * x + a22 * y + dy)


def invert_xf_for_point(x_ali: float, y_ali: float, row: Union[Mapping, pd.Series]) -> Tuple[float, float]:
    """Map an aligned-stack point (x_ali, y_ali) back to raw-stack coordinates."""
    inv = xf_inverse_matrix(row)
    v = inv @ np.array([x_ali, y_ali, 1.0])
    return (float(v[0]), float(v[1]))
