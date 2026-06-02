#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/io/star_writer.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
Multi-block STAR writer.

A STAR file is a sequence of named blocks. tomoJANAS uses two kinds:

* **loop blocks** — a ``loop_`` table with ordered columns and N rows. Used
  for RELION-compatible particle/tomogram tables and for tomoJANAS tables.
* **pair blocks** — a set of ``_key value`` lines (no ``loop_``). Used for
  small key/value metadata (e.g. manifest summaries, settings).

Design goals:
* preserve column order exactly;
* keep RELION-compatible blocks free of tomoJANAS-specific tags (the caller
  decides which columns go in which block);
* represent missing values (``None``/NaN/empty) consistently as ``na_rep``.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Sequence, Union

__all__ = [
    "LoopBlock",
    "PairBlock",
    "write_star",
    "loop_from_records",
    "loop_from_dataframe",
    "format_value",
]


@dataclass
class LoopBlock:
    name: str
    columns: List[str]
    rows: List[Sequence[Any]] = field(default_factory=list)


@dataclass
class PairBlock:
    name: str
    pairs: "Dict[str, Any]" = field(default_factory=dict)


Block = Union[LoopBlock, PairBlock]


def format_value(v: Any, na_rep: str = "?") -> str:
    """Render a single STAR cell value as a token (no whitespace)."""
    if v is None:
        return na_rep
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, float):
        if math.isnan(v):
            return na_rep
        if math.isinf(v):
            return na_rep
        s = f"{v:.6f}"
        if "." in s:
            s = s.rstrip("0").rstrip(".")
        return s if s not in ("", "-0") else "0"
    if isinstance(v, (int,)):
        return str(v)
    s = str(v).strip()
    if s == "":
        return na_rep
    # tokens must not contain whitespace; collapse internal whitespace
    return "_".join(s.split()) if any(c.isspace() for c in s) else s


def _write_loop(handle, block: LoopBlock, na_rep: str) -> None:
    handle.write(f"data_{block.name}\n\n")
    handle.write("loop_\n")
    for i, col in enumerate(block.columns, start=1):
        handle.write(f"{col} #{i}\n")
    # Pre-format all cells and compute per-column widths for readable alignment.
    formatted: List[List[str]] = []
    for row in block.rows:
        if isinstance(row, Mapping):
            cells = [format_value(row.get(c), na_rep) for c in block.columns]
        else:
            if len(row) != len(block.columns):
                raise ValueError(
                    f"row has {len(row)} values but block '{block.name}' has "
                    f"{len(block.columns)} columns"
                )
            cells = [format_value(v, na_rep) for v in row]
        formatted.append(cells)
    widths = [len(c) for c in block.columns]
    for cells in formatted:
        for j, cell in enumerate(cells):
            if len(cell) > widths[j]:
                widths[j] = len(cell)
    for cells in formatted:
        line = "  ".join(cell.ljust(widths[j]) for j, cell in enumerate(cells))
        handle.write(line.rstrip() + "\n")
    handle.write("\n")


def _write_pairs(handle, block: PairBlock, na_rep: str) -> None:
    handle.write(f"data_{block.name}\n\n")
    keys = list(block.pairs.keys())
    kw = max((len(k) for k in keys), default=0)
    for k, v in block.pairs.items():
        handle.write(f"{k.ljust(kw)}  {format_value(v, na_rep)}\n")
    handle.write("\n")


def write_star(
    path: str,
    blocks: Sequence[Block],
    na_rep: str = "?",
    header_comment: str | None = None,
) -> None:
    """Write a list of blocks (in order) to a STAR file at ``path``."""
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        if header_comment:
            for line in str(header_comment).splitlines():
                f.write(f"# {line}\n")
            f.write("\n")
        for block in blocks:
            if isinstance(block, LoopBlock):
                _write_loop(f, block, na_rep)
            elif isinstance(block, PairBlock):
                _write_pairs(f, block, na_rep)
            else:  # pragma: no cover - defensive
                raise TypeError(f"Unknown block type: {type(block)!r}")


def loop_from_records(
    name: str, columns: Sequence[str], records: Sequence[Mapping[str, Any]]
) -> LoopBlock:
    """Build a LoopBlock from a list of dict records (missing keys -> None)."""
    cols = list(columns)
    rows = [[rec.get(c) for c in cols] for rec in records]
    return LoopBlock(name=name, columns=cols, rows=rows)


def loop_from_dataframe(name: str, df, columns: Sequence[str] | None = None) -> LoopBlock:
    """Build a LoopBlock from a pandas DataFrame, preserving column order."""
    cols = list(columns) if columns is not None else list(df.columns)
    rows = df[cols].values.tolist()
    return LoopBlock(name=name, columns=cols, rows=rows)
