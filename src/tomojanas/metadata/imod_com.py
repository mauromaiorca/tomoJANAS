#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/metadata/imod_com.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
Parser for IMOD ``.com`` command files (newst.com, tilt.com, ctfplotter.com).

A ``.com`` file is a sequence of lines. Lines starting with ``$`` are program
invocations and are ignored for metadata purposes. Other non-empty lines are
``KEY value...`` pairs. Values may be space- or comma-separated.

We keep the parse permissive: every key seen is stored (raw remainder string),
and helpers extract typed values for the keys the importer cares about.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Set, Tuple

__all__ = [
    "parse_com_file",
    "com_get",
    "com_get_float",
    "com_get_int",
    "com_get_numbers",
    "parse_exclude_list",
]


def parse_com_file(path: str) -> Dict[str, str]:
    """Parse a ``.com`` file into ``{KEY: "remainder string"}``.

    Lines beginning with ``$`` (program calls) and ``#`` (comments) are
    skipped. The first whitespace-delimited token is the key; the rest of the
    line (stripped) is the value. Later occurrences overwrite earlier ones.
    """
    out: Dict[str, str] = {}
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("$") or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            key = parts[0]
            value = parts[1].strip() if len(parts) > 1 else ""
            out[key] = value
    return out


def com_get(com: Dict[str, str], key: str, default=None) -> Optional[str]:
    return com.get(key, default)


def com_get_numbers(value: Optional[str]) -> List[float]:
    """Split a value on commas/whitespace and parse floats (best effort)."""
    if value is None:
        return []
    toks = value.replace(",", " ").split()
    nums: List[float] = []
    for t in toks:
        try:
            nums.append(float(t))
        except ValueError:
            break  # stop at first non-numeric token (e.g. trailing comment)
    return nums


def com_get_float(com: Dict[str, str], key: str, default=None) -> Optional[float]:
    nums = com_get_numbers(com.get(key))
    return nums[0] if nums else default


def com_get_int(com: Dict[str, str], key: str, default=None) -> Optional[int]:
    v = com_get_float(com, key, None)
    return int(round(v)) if v is not None else default


def parse_exclude_list(value: Optional[str], n_frames: int) -> Set[int]:
    """Parse an IMOD EXCLUDELIST value into a set of 0-based frame indices.

    IMOD exclusion indices are 1-based and ranges are inclusive, e.g.
    ``"1,2,5-7"`` -> {0,1,4,5,6}. Faithful to RELION's readTiltCom range
    handling. Out-of-range indices (``< 1`` or ``> n_frames``) raise ValueError.
    """
    excluded: Set[int] = set()
    if not value:
        return excluded
    # keep only the first token group (IMOD writes a single comma list)
    token = value.split()[0] if value.split() else ""
    for part in token.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a_str, b_str = part.split("-", 1)
            a, b = int(a_str), int(b_str)
            if a < 1 or a > n_frames or b < 1 or b > n_frames:
                raise ValueError(f"bad exclusion range '{part}' (n_frames={n_frames})")
            lo, hi = (a, b) if a <= b else (b, a)
            for x in range(lo, hi + 1):
                excluded.add(x - 1)
        else:
            x = int(part)
            if x < 1 or x > n_frames:
                raise ValueError(f"bad exclusion index '{x}' (n_frames={n_frames})")
            excluded.add(x - 1)
    return excluded
