#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/metadata/ctf.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
CTF metadata parsing (IMOD CtfPlotter ``.defocus``, CTFFind STAR) and unit
handling. This module parses metadata ONLY — it never touches image pixels.

RELION expects defocus in Ångström. CtfPlotter ``.defocus`` files store
defocus in nanometres (RELION multiplies by 10 on import).
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional, Tuple

__all__ = [
    "defocus_to_angstrom",
    "infer_defocus_unit",
    "read_ctfplotter_defocus",
    "detect_ctf_source",
]

_NM_TO_A = 10.0
_UM_TO_A = 10000.0


def defocus_to_angstrom(value: float, unit: str) -> float:
    """Convert a defocus value to Ångström. ``unit`` in {angstrom, nm, micrometer}."""
    u = (unit or "").lower()
    if u in ("a", "angstrom", "angstroms", "ang"):
        return float(value)
    if u in ("nm", "nanometer", "nanometre", "nanometers"):
        return float(value) * _NM_TO_A
    if u in ("um", "µm", "micrometer", "micrometre", "micron", "microns"):
        return float(value) * _UM_TO_A
    raise ValueError(f"unknown defocus unit: {unit!r}")


def infer_defocus_unit(values: List[float]) -> str:
    """Heuristic unit inference from the magnitude of defocus values.

    Typical cryo-ET defocus ~1–8 µm. Expressed as:
      Ångström  -> ~10000–80000
      nanometre -> ~1000–8000
      micrometre-> ~1–8
    Returns one of {"angstrom", "nm", "micrometer"}; the caller MUST record
    the decision in metadata/logs (never silently trust it).
    """
    vals = [abs(v) for v in values if v is not None]
    if not vals:
        return "nm"
    med = sorted(vals)[len(vals) // 2]
    if med < 100.0:
        return "micrometer"
    if med < 1000.0:
        # ambiguous low-nm / could be nm; default nm
        return "nm"
    if med < 9000.0:
        return "nm"
    return "angstrom"


def read_ctfplotter_defocus(path: str) -> Dict[str, object]:
    """Parse an IMOD CtfPlotter ``.defocus`` file (best effort, permissive).

    Returns ``{"entries": [...], "has_astigmatism": bool, "raw_unit": "nm"}``
    where each entry is a dict with keys: ``view_start, view_end,
    angle_start, angle_end, defocus_u, defocus_v, astig_angle, phase_shift``.
    Defocus values are returned in the file's native unit (nm); conversion to
    Å is the caller's responsibility (so the unit decision is explicit).
    """
    entries: List[Dict[str, float]] = []
    has_astig = False
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = [ln.strip() for ln in f if ln.strip() and not ln.strip().startswith("#")]

    if not lines:
        return {"entries": [], "has_astigmatism": False, "raw_unit": "nm"}

    # Optional header/flags line: 6 tokens where the "defocus" slot (idx 4) is ~0.
    start = 0
    first = lines[0].split()
    if len(first) >= 6:
        try:
            if abs(float(first[4])) < 1e-6:
                has_astig = int(float(first[5])) in (1, 3)
                start = 1
        except ValueError:
            pass

    for ln in lines[start:]:
        toks = ln.split()
        if len(toks) < 5:
            continue
        try:
            vs, ve = int(float(toks[0])), int(float(toks[1]))
            a_start, a_end = float(toks[2]), float(toks[3])
            du = float(toks[4])
        except ValueError:
            continue
        dv = du
        astig_angle = 0.0
        phase = 0.0
        if len(toks) >= 7:
            # extended: defocusU defocusV astigAngle [phaseShift]
            try:
                dv = float(toks[5])
                astig_angle = float(toks[6])
                has_astig = True
            except ValueError:
                dv = du
            if len(toks) >= 8:
                try:
                    phase = float(toks[7])
                except ValueError:
                    phase = 0.0
        entries.append(
            {
                "view_start": vs,
                "view_end": ve,
                "angle_start": a_start,
                "angle_end": a_end,
                "defocus_u": du,
                "defocus_v": dv,
                "astig_angle": astig_angle,
                "phase_shift": phase,
            }
        )
    return {"entries": entries, "has_astigmatism": has_astig, "raw_unit": "nm"}


def detect_ctf_source(
    imod_dir: str,
    basename: str,
    *,
    explicit_ctfplotter: Optional[str] = None,
    explicit_defocus: Optional[str] = None,
    explicit_ctffind: Optional[str] = None,
    ctfplotter_info: Optional[str] = None,
    ctfplotter_log: Optional[str] = None,
) -> Tuple[str, Optional[str]]:
    """Resolve the CTF source following the documented priority.

    Returns ``(source, path)`` where source is one of
    {"ctfplotter", "ctffind", "defocus", "none"} and path may be None.
    """
    # 1. explicit files
    if explicit_ctfplotter and os.path.isfile(explicit_ctfplotter):
        return ("ctfplotter", explicit_ctfplotter)
    if explicit_defocus and os.path.isfile(explicit_defocus):
        return ("ctfplotter", explicit_defocus)
    if explicit_ctffind and os.path.isfile(explicit_ctffind):
        return ("ctffind", explicit_ctffind)

    # 2. file referenced in ctfplotter.com
    com_path = os.path.join(imod_dir, "ctfplotter.com")
    if os.path.isfile(com_path):
        from .imod_com import parse_com_file
        com = parse_com_file(com_path)
        ref = com.get("DefocusFile") or com.get("SaveAndExit") or com.get("OutputFile")
        if ref:
            cand = ref.split()[0]
            cand_path = cand if os.path.isabs(cand) else os.path.join(imod_dir, cand)
            if os.path.isfile(cand_path):
                return ("ctfplotter", cand_path)

    # 3. <basename>.defocus
    cand = os.path.join(imod_dir, f"{basename}.defocus")
    if os.path.isfile(cand):
        return ("ctfplotter", cand)

    # 4. ctfplotter.info / explicit info
    for cand in [ctfplotter_info, os.path.join(imod_dir, "ctfplotter.info")]:
        if cand and os.path.isfile(cand):
            return ("ctfplotter", cand)

    # 5. ctfplotter.log / explicit log
    for cand in [ctfplotter_log, os.path.join(imod_dir, "ctfplotter.log")]:
        if cand and os.path.isfile(cand):
            return ("ctfplotter", cand)

    # 6. ctf3d*.log
    try:
        for name in sorted(os.listdir(imod_dir)):
            if name.startswith("ctf3d") and name.endswith(".log"):
                return ("ctfplotter", os.path.join(imod_dir, name))
    except OSError:
        pass

    # 7. nothing found
    return ("none", None)
