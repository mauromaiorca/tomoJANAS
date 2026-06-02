#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/io/star_reader.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
Multi-block STAR reader.

Parses a STAR file into an ordered mapping of block-name -> parsed block.
Each block is one of:

* ``{"type": "loop", "columns": [...], "df": DataFrame}`` — a loop_ table
  (values kept as strings to avoid lossy coercion; coerce downstream).
* ``{"type": "pairs", "pairs": {key: value}}`` — bare ``_key value`` lines.

Whitespace-tolerant. Lines beginning with ``#`` are comments. Handles the
files written by :mod:`tomojanas.io.star_writer` as well as typical RELION
multi-block files (``data_optics`` + ``data_particles``).
"""

from __future__ import annotations

from collections import OrderedDict
from typing import Dict, List, Optional

import pandas as pd

__all__ = [
    "read_star",
    "read_star_loop",
    "read_star_pairs",
    "star_block_names",
]


def _strip_index_suffix(token: str) -> str:
    """``_rlnTomoName #1`` -> ``_rlnTomoName`` (drop the ``#N`` if present)."""
    return token.split()[0].strip()


def read_star(path: str) -> "OrderedDict[str, dict]":
    """Parse a STAR file into ``OrderedDict[name -> block]``."""
    blocks: "OrderedDict[str, dict]" = OrderedDict()
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    current_name: Optional[str] = None
    state = "idle"  # idle | loop_header | loop_rows | pairs
    columns: List[str] = []
    rows: List[List[str]] = []
    pairs: "OrderedDict[str, str]" = OrderedDict()

    def flush():
        nonlocal current_name, state, columns, rows, pairs
        if current_name is None:
            return
        if state in ("loop_header", "loop_rows"):
            df = pd.DataFrame(rows, columns=columns, dtype=object)
            blocks[current_name] = {"type": "loop", "columns": list(columns), "df": df}
        elif state == "pairs":
            blocks[current_name] = {"type": "pairs", "pairs": dict(pairs)}
        else:
            # empty block
            blocks[current_name] = {"type": "pairs", "pairs": {}}
        columns = []
        rows = []
        pairs = OrderedDict()

    for raw in lines:
        line = raw.strip()
        if line == "" or line.startswith("#"):
            continue
        if line.startswith("data_"):
            flush()
            current_name = line[len("data_"):].strip()
            state = "idle"
            columns, rows, pairs = [], [], OrderedDict()
            continue
        if current_name is None:
            # stray content before any data_ block: ignore
            continue
        if line == "loop_":
            state = "loop_header"
            columns, rows = [], []
            continue

        if state == "loop_header":
            if line.startswith("_"):
                columns.append(_strip_index_suffix(line))
                continue
            # first data row
            state = "loop_rows"
            rows.append(line.split())
            continue

        if state == "loop_rows":
            rows.append(line.split())
            continue

        # not in a loop: a bare "_key value" line -> pairs block
        if line.startswith("_"):
            state = "pairs"
            parts = line.split(None, 1)
            key = parts[0]
            value = parts[1].strip() if len(parts) > 1 else ""
            pairs[key] = value
            continue
        # unknown stray line inside a block: ignore

    flush()
    return blocks


def star_block_names(path: str) -> List[str]:
    return list(read_star(path).keys())


def read_star_loop(path: str, block_name: str) -> pd.DataFrame:
    """Return the loop table of ``block_name`` as a string DataFrame."""
    blocks = read_star(path)
    if block_name not in blocks:
        raise KeyError(f"block 'data_{block_name}' not found in {path}")
    blk = blocks[block_name]
    if blk["type"] != "loop":
        raise TypeError(f"block 'data_{block_name}' is not a loop block")
    return blk["df"]


def read_star_pairs(path: str, block_name: str) -> Dict[str, str]:
    blocks = read_star(path)
    if block_name not in blocks:
        raise KeyError(f"block 'data_{block_name}' not found in {path}")
    blk = blocks[block_name]
    if blk["type"] != "pairs":
        raise TypeError(f"block 'data_{block_name}' is not a pairs block")
    return blk["pairs"]
