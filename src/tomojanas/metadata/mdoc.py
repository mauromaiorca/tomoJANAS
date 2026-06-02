#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/metadata/mdoc.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
Permissive parser for SerialEM ``.mdoc`` files.

An mdoc has a global header section followed by per-image ``[ZValue = N]``
sections, each a block of ``Key = value`` lines. We never fail on missing
fields; absent values are simply omitted.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional

__all__ = ["read_mdoc"]

_ZVALUE_RE = re.compile(r"^\[\s*ZValue\s*=\s*(\d+)\s*\]", re.IGNORECASE)
_SECTION_RE = re.compile(r"^\[(.+)\]\s*$")


def _coerce(value: str) -> Any:
    """Best-effort scalar/vector coercion of an mdoc value string."""
    toks = value.split()
    out = []
    for t in toks:
        try:
            out.append(int(t))
            continue
        except ValueError:
            pass
        try:
            out.append(float(t))
            continue
        except ValueError:
            out.append(t)
    if len(out) == 1:
        return out[0]
    return out


def read_mdoc(path: str) -> Dict[str, Any]:
    """Parse an mdoc into ``{"global": {...}, "sections": [...], ...}``.

    Convenience top-level keys (when available): ``pixel_spacing`` (Å) and
    ``tilt_angles`` (list, ordered by ZValue).
    """
    global_block: Dict[str, Any] = {}
    sections: List[Dict[str, Any]] = []
    current: Optional[Dict[str, Any]] = None

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            zmatch = _ZVALUE_RE.match(line)
            if zmatch:
                current = {"ZValue": int(zmatch.group(1))}
                sections.append(current)
                continue
            smatch = _SECTION_RE.match(line)
            if smatch and "=" not in line:
                # a non-ZValue bracket section (e.g. [T = ...]); start a generic block
                current = {"_section": smatch.group(1).strip()}
                sections.append(current)
                continue
            if "=" in line:
                key, _, val = line.partition("=")
                key = key.strip()
                val = val.strip()
                target = current if current is not None else global_block
                target[key] = _coerce(val)

    # convenience extractions
    pixel_spacing = global_block.get("PixelSpacing")
    if pixel_spacing is None:
        # sometimes only present per-section
        for s in sections:
            if "PixelSpacing" in s:
                pixel_spacing = s["PixelSpacing"]
                break

    tilt_angles: List[float] = []
    for s in sections:
        if "TiltAngle" in s:
            try:
                tilt_angles.append(float(s["TiltAngle"]))
            except (TypeError, ValueError):
                pass

    return {
        "global": global_block,
        "sections": sections,
        "pixel_spacing": float(pixel_spacing) if isinstance(pixel_spacing, (int, float)) else None,
        "tilt_angles": tilt_angles,
        "num_sections": len(sections),
    }
