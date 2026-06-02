#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/importers/ctf_importer.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
``tomojanas-import ctf`` — import CTF metadata without modifying image pixels.

Updates ``tilt_series/<tomo>.star`` and writes ``imod_settings/ctf_settings.star``.
Never touches image data; never creates premultiplied stacks.
"""
from __future__ import annotations

import os
import sys
from typing import Dict, List, Optional

from tomojanas import get_version
from tomojanas.io.star_writer import LoopBlock, PairBlock, write_star
from tomojanas.io.star_reader import read_star
from tomojanas.io.logs import ImportLogger, now_iso
from tomojanas.io import project_writer as pw
from tomojanas.models.project import Project
from tomojanas.metadata.ctf import detect_ctf_source, read_ctfplotter_defocus, defocus_to_angstrom, infer_defocus_unit

__all__ = ["import_ctf"]


def import_ctf(args) -> int:
    project_root = os.path.abspath(args.project)
    tomo_name = args.tomo_name

    logger = ImportLogger(project_root, get_version())
    logger.invocation("tomojanas-import ctf", argv=sys.argv)

    try:
        proj = Project.open_or_create(project_root, fail_if_existing=False, append=True)
    except Exception as exc:
        logger.error(str(exc))
        return 1

    # Load tomogram metadata for basename and imod_dir
    tomo_meta = _load_tomo_meta(proj.tomograms_star, tomo_name)
    if tomo_meta is None:
        logger.error(f"Tomogram '{tomo_name}' not found in tomograms.star. Import imod first.")
        return 1

    imod_dir = tomo_meta.get("_tomoJANASImodDir", project_root) or project_root
    # Resolve relative imod_dir
    if not os.path.isabs(imod_dir):
        imod_dir = os.path.join(project_root, imod_dir)

    # Detect CTF source
    ctf_src = getattr(args, "ctf_source", "auto") or "auto"
    if ctf_src == "none":
        logger.info("CTF source set to none; marking CTF as missing")
        _update_ctf_settings(pw.imod_settings_dir(project_root, tomo_name),
                             tomo_name, "none", None, "none", now_iso())
        logger.summary()
        return 0

    # Get basename for auto-detection
    basename = os.path.basename(imod_dir)  # fallback

    ctf_source, ctf_path = detect_ctf_source(
        imod_dir, basename,
        explicit_ctfplotter=getattr(args, "ctfplotter_file", None),
        explicit_defocus=getattr(args, "defocus_file", None),
        explicit_ctffind=getattr(args, "ctffind_star", None),
        ctfplotter_info=getattr(args, "ctfplotter_info", None),
        ctfplotter_log=getattr(args, "ctfplotter_log", None),
    )

    if ctf_source == "none" or ctf_path is None:
        logger.warning(f"No CTF file found for {tomo_name}; marking CTF as missing")
        _update_ctf_settings(pw.imod_settings_dir(project_root, tomo_name),
                             tomo_name, "none", None, "missing", now_iso())
        logger.summary()
        return 0

    logger.info(f"CTF source: {ctf_source} @ {ctf_path}")

    # Parse defocus values
    defocus_entries = []
    if ctf_source in ("ctfplotter",):
        try:
            data = read_ctfplotter_defocus(ctf_path)
            defocus_entries = data.get("entries", [])
        except Exception as exc:
            logger.warning(f"Cannot parse CTF file {ctf_path}: {exc}")

    # Determine defocus unit
    defocus_unit = getattr(args, "defocus_unit", "auto") or "auto"
    if defocus_unit == "auto" and defocus_entries:
        raw_vals = [e["defocus_u"] for e in defocus_entries]
        inferred_unit = infer_defocus_unit(raw_vals)
        logger.info(f"Inferred defocus unit: {inferred_unit} (median ~{sorted(raw_vals)[len(raw_vals)//2]:.1f})")
        defocus_unit_resolved = inferred_unit
    elif defocus_unit == "auto":
        defocus_unit_resolved = "nm"
    else:
        defocus_unit_resolved = defocus_unit

    # Convert to Å and build per-tilt CTF table
    # ctfplotter .defocus has view_start..view_end ranges; we expand to per-view rows
    tilt_ctf_rows = _expand_defocus_entries(defocus_entries, defocus_unit_resolved,
                                             getattr(args, "phase_shift", None))

    # Load existing tilt STAR and update CTF columns
    tilt_star = pw.tilt_series_star_path(project_root, tomo_name)
    if os.path.isfile(tilt_star):
        _update_tilt_star_ctf(tilt_star, tomo_name, tilt_ctf_rows, logger)
    else:
        logger.warning(f"tilt_series STAR not found: {tilt_star}")

    # Write ctf_settings.star
    _update_ctf_settings(
        pw.imod_settings_dir(project_root, tomo_name),
        tomo_name, ctf_source, ctf_path, defocus_unit_resolved, now_iso()
    )

    logger.info(f"CTF updated for {tomo_name}: {len(tilt_ctf_rows)} tilt entries")
    logger.summary()
    return 0


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def _load_tomo_meta(tomograms_star: str, tomo_name: str) -> Optional[Dict]:
    if not os.path.isfile(tomograms_star):
        return None
    try:
        blks = read_star(tomograms_star)
        # Check both RELION global and tomoJANAS source blocks
        for bname in ("global", "tomoJANAS_tomogram_sources"):
            if bname in blks and blks[bname]["type"] == "loop":
                df = blks[bname]["df"]
                for _, r in df.iterrows():
                    name_col = "_rlnTomoName" if bname == "global" else "_tomoJANASTomoName"
                    if str(r.get(name_col, "")) == tomo_name:
                        return dict(r)
    except Exception:
        pass
    return None


def _expand_defocus_entries(entries: List[Dict], unit: str, phase_shift_override) -> Dict[int, Dict]:
    """Expand view-range defocus entries to {view_index: ctf_dict}."""
    result: Dict[int, Dict] = {}
    for entry in entries:
        vs = int(entry.get("view_start", 1))
        ve = int(entry.get("view_end", vs))
        du_A = defocus_to_angstrom(float(entry.get("defocus_u", 0.0)), unit)
        dv_A = defocus_to_angstrom(float(entry.get("defocus_v", entry.get("defocus_u", 0.0))), unit)
        astig_angle = float(entry.get("astig_angle", 0.0))
        astigmatism = abs(du_A - dv_A)
        phase = float(phase_shift_override) if phase_shift_override is not None else float(entry.get("phase_shift", 0.0))
        for v in range(vs, ve + 1):
            result[v] = {
                "_rlnDefocusU": f"{du_A:.2f}",
                "_rlnDefocusV": f"{dv_A:.2f}",
                "_rlnDefocusAngle": f"{astig_angle:.2f}",
                "_rlnCtfAstigmatism": f"{astigmatism:.2f}",
                "_rlnCtfFigureOfMerit": "?",
                "_rlnCtfMaxResolution": "?",
                "_rlnCtfScalefactor": "?",
                "_rlnPhaseShift": f"{phase:.2f}",
            }
    return result


def _update_tilt_star_ctf(tilt_star: str, tomo_name: str,
                           tilt_ctf: Dict[int, Dict], logger) -> None:
    """Write updated CTF columns into the RELION per-tilt block."""
    from tomojanas.metadata.relion_labels import TILT_SERIES_COLUMNS
    try:
        blks = read_star(tilt_star)
    except Exception as exc:
        logger.warning(f"Cannot read {tilt_star}: {exc}")
        return

    if tomo_name not in blks or blks[tomo_name]["type"] != "loop":
        logger.warning(f"Block 'data_{tomo_name}' not found in {tilt_star}")
        return

    df = blks[tomo_name]["df"]
    new_rows = []
    for _, r in df.iterrows():
        row = dict(r)
        tilt_idx = None
        try:
            tilt_idx = int(float(row.get("_rlnTomoTiltMovieIndex", 0)))
        except (TypeError, ValueError):
            pass
        if tilt_idx and tilt_idx in tilt_ctf:
            row.update(tilt_ctf[tilt_idx])
        new_rows.append(row)

    # Rebuild file: replace RELION block, keep other blocks
    out_blocks = []
    for bname, bdata in blks.items():
        if bname == tomo_name and bdata["type"] == "loop":
            out_blocks.append(LoopBlock(name=tomo_name, columns=TILT_SERIES_COLUMNS, rows=new_rows))
        elif bdata["type"] == "loop":
            out_blocks.append(LoopBlock(name=bname, columns=bdata["columns"], rows=bdata["df"].values.tolist()))
        else:
            out_blocks.append(PairBlock(name=bname, pairs=bdata.get("pairs", {})))
    write_star(tilt_star, out_blocks)


def _update_ctf_settings(imod_settings_dir: str, tomo_name: str,
                          ctf_source: str, ctf_path: Optional[str],
                          unit: str, timestamp: str) -> None:
    os.makedirs(imod_settings_dir, exist_ok=True)
    write_star(os.path.join(imod_settings_dir, "ctf_settings.star"), [
        PairBlock(name="tomoJANAS_ctf_settings", pairs={
            "_tomoJANASTomoName": tomo_name,
            "_tomoJANASCtfSource": ctf_source,
            "_tomoJANASCtfFile": ctf_path or "?",
            "_tomoJANASCtfDefocusUnit": unit,
            "_tomoJANASCtfPremultiplied": "0",
            "_tomoJANASCtfImportTimestamp": timestamp,
        })
    ])
