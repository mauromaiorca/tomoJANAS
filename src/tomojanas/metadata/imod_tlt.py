#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/metadata/imod_tlt.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""Parser for IMOD ``.tlt`` / ``.rawtlt`` tilt-angle files (one angle per line)."""

from __future__ import annotations

from typing import List

import pandas as pd

__all__ = ["read_tlt", "read_tlt_angles"]


def read_tlt_angles(path: str) -> List[float]:
    """Return the list of tilt angles (degrees), one per non-empty line."""
    angles: List[float] = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            # a .tlt line is usually a single float; take the first token
            tok = line.split()[0]
            angles.append(float(tok))
    return angles


def read_tlt(path: str) -> pd.DataFrame:
    """Return a DataFrame with ``tilt_index`` (1-based) and ``tilt_angle``."""
    angles = read_tlt_angles(path)
    return pd.DataFrame(
        {
            "tilt_index": list(range(1, len(angles) + 1)),
            "tilt_angle": angles,
        }
    )
