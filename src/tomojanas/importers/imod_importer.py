#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/importers/imod_importer.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
``tomojanas-import imod`` — import an IMOD reconstruction directory into a
tomoJANAS project.

Creates or updates:
    <project>/tomograms.star
    <project>/optimisation_set.star
    <project>/tilt_series/<tomo>.star
    <project>/tilt_series/<tomo>/imod_settings/*.star
    <project>/project_manifest.{json,star}
    <project>/logs/tomojanas_import.log
"""
from __future__ import annotations

import os
import sys
from typing import Optional

import numpy as np

from tomojanas import get_version
from tomojanas.io.mrc import read_mrc_header
from tomojanas.io.star_writer import LoopBlock, PairBlock, write_star
from tomojanas.io.star_reader import read_star
from tomojanas.io.logs import ImportLogger, now_iso
from tomojanas.io import project_writer as pw
from tomojanas.models.project import Project, ProjectExistsError, ProjectMissingError
from tomojanas.metadata.relion_labels import (
    TOMOGRAMS_GLOBAL_COLUMNS, TOMOGRAMS_SOURCE_COLUMNS, TILT_SERIES_COLUMNS,
    TILT_MAPPING_COLUMNS, OPTIMISATION_SET_COLUMNS, PROJECTION_MATRIX_COLUMNS,
    RelionGeometryStatus,
)
from tomojanas.metadata.imod_com import parse_com_file
from tomojanas.metadata.imod_tlt import read_tlt
from tomojanas.metadata.imod_xf import read_xf
from tomojanas.metadata.ctf import detect_ctf_source
from tomojanas.geometry.imod_mapping import build_relion_imod_mapping

__all__ = ["import_imod_project"]


def _resolve_file(explicit: Optional[str], imod_dir: str, basename: str, suffixes) -> Optional[str]:
    """Return explicit path if it exists, else try basename+suffix in imod_dir."""
    if explicit and os.path.isfile(explicit):
        return os.path.abspath(explicit)
    if explicit and os.path.isfile(os.path.join(imod_dir, explicit)):
        return os.path.abspath(os.path.join(imod_dir, explicit))
    for suf in suffixes:
        cand = os.path.join(imod_dir, basename + suf)
        if os.path.isfile(cand):
            return os.path.abspath(cand)
    return None


def import_imod_project(args) -> int:
    """Main entry point for ``tomojanas-import imod``."""
    project_root = os.path.abspath(args.project)
    imod_dir = os.path.abspath(args.imod_dir)
    tomo_name = args.tomo_name.rstrip("/").rstrip("\\")
    basename = getattr(args, "basename", None) or tomo_name

    logger = ImportLogger(project_root, get_version())
    logger.invocation("tomojanas-import imod", argv=sys.argv)

    # ------------------------------------------------------------------ #
    # Project open/create
    # ------------------------------------------------------------------ #
    try:
        proj = Project.open_or_create(
            project_root,
            create_if_missing=getattr(args, "create_if_missing", False),
            overwrite=getattr(args, "overwrite", False),
            append=getattr(args, "append", False),
            update_existing=getattr(args, "update_existing", False),
            fail_if_existing=getattr(args, "fail_if_existing", True),
        )
    except (ProjectExistsError, ProjectMissingError) as exc:
        logger.error(str(exc))
        return 1

    pw.ensure_tomogram_dirs(project_root, tomo_name)
    relative = getattr(args, "relative_paths", True)

    # ------------------------------------------------------------------ #
    # Resolve IMOD files
    # ------------------------------------------------------------------ #
    raw_stack = _resolve_file(getattr(args, "raw_stack", None), imod_dir, basename,
                               [".mrc", ".st", ".rawst"])
    mdoc_file = _resolve_file(getattr(args, "mdoc", None), imod_dir, basename,
                               [".mrc.mdoc", ".mdoc"])
    ali_stack = _resolve_file(getattr(args, "ali_stack", None), imod_dir, basename,
                               ["_ali.mrc", "_ali.st", ".ali"])
    rec_tomo  = _resolve_file(getattr(args, "rec_tomo", None), imod_dir, basename,
                               ["_rec.mrc", "_full.mrc"])
    xf_file   = _resolve_file(getattr(args, "xf", None), imod_dir, basename, [".xf"])
    tlt_file  = _resolve_file(getattr(args, "tlt", None), imod_dir, basename, [".tlt"])
    newst_com = _resolve_file(getattr(args, "newst_com", None), imod_dir, "", [""])
    if newst_com is None:
        newst_com = os.path.join(imod_dir, "newst.com") if os.path.isfile(os.path.join(imod_dir, "newst.com")) else None
    tilt_com  = _resolve_file(getattr(args, "tilt_com", None), imod_dir, "", [""])
    if tilt_com is None:
        tilt_com = os.path.join(imod_dir, "tilt.com") if os.path.isfile(os.path.join(imod_dir, "tilt.com")) else None

    logger.info(f"imod_dir     = {imod_dir}")
    logger.info(f"tomo_name    = {tomo_name}")
    logger.info(f"ali_stack    = {ali_stack}")
    logger.info(f"rec_tomo     = {rec_tomo}")
    logger.info(f"raw_stack    = {raw_stack}")

    # ------------------------------------------------------------------ #
    # MRC headers
    # ------------------------------------------------------------------ #
    ali_hdr = None
    rec_hdr = None
    raw_hdr = None
    if ali_stack and os.path.isfile(ali_stack):
        try:
            ali_hdr = read_mrc_header(ali_stack)
        except Exception as exc:
            logger.warning(f"Cannot read ali stack header: {exc}")
    if rec_tomo and os.path.isfile(rec_tomo):
        try:
            rec_hdr = read_mrc_header(rec_tomo)
        except Exception as exc:
            logger.warning(f"Cannot read rec tomo header: {exc}")
    if raw_stack and os.path.isfile(raw_stack):
        try:
            raw_hdr = read_mrc_header(raw_stack)
        except Exception as exc:
            logger.warning(f"Cannot read raw stack header: {exc}")

    # ------------------------------------------------------------------ #
    # Pixel sizes and binning
    # ------------------------------------------------------------------ #
    a_raw = getattr(args, "raw_pixel_size", None)
    a_ali = getattr(args, "ali_pixel_size", None)
    a_rec = getattr(args, "rec_pixel_size", None)

    if getattr(args, "infer_pixel_size_from_mrc", True):
        if a_raw is None and raw_hdr:
            a_raw = raw_hdr.pixel_x
        if a_ali is None and ali_hdr:
            a_ali = ali_hdr.pixel_x
        if a_rec is None and rec_hdr:
            a_rec = rec_hdr.pixel_x

    # mdoc fallback
    if getattr(args, "infer_pixel_size_from_mdoc", False) and mdoc_file:
        from tomojanas.metadata.mdoc import read_mdoc
        try:
            mdoc_data = read_mdoc(mdoc_file)
            ps = mdoc_data.get("pixel_spacing")
            if ps and a_raw is None:
                a_raw = float(ps)
        except Exception:
            pass

    if a_ali is None:
        a_ali = a_raw or 1.0
        logger.warning("Ali pixel size unknown; using raw pixel size or 1.0 Å")
    if a_rec is None:
        a_rec = a_ali
        logger.warning("Rec pixel size unknown; assuming same as ali pixel size")
    if a_raw is None:
        a_raw = a_ali

    B_rec_ali = float(a_rec) / float(a_ali) if float(a_ali) > 0 else 1.0
    if getattr(args, "require_ali_rec_same_bin", False):
        tol = 1e-4
        if abs(B_rec_ali - 1.0) > tol:
            logger.error(f"--require-ali-rec-same-bin set but B_rec_ali={B_rec_ali:.6f} != 1.0")
            return 1

    # ------------------------------------------------------------------ #
    # RELION IMOD mapping
    # ------------------------------------------------------------------ #
    mapping = None
    geom_status = RelionGeometryStatus.MATRICES_MISSING

    if ali_stack and os.path.isfile(ali_stack) and (newst_com or tilt_com):
        try:
            # We pass the ali stack as the tilt-series stack (it has the right nz)
            mapping = build_relion_imod_mapping(
                ts_path=ali_stack,
                imod_dir=imod_dir,
                newst_com=os.path.basename(newst_com) if newst_com else "newst.com",
                tilt_com=os.path.basename(tilt_com) if tilt_com else "tilt.com",
                ali=True,
                ali_size=True,
                flip_yz=getattr(args, "flip_yz", False),
                flip_z=getattr(args, "flip_z", False),
                flip_angles=getattr(args, "flip_angles", False),
                offset_x=getattr(args, "import_offset_x", 0.0) or 0.0,
                offset_y=getattr(args, "import_offset_y", 0.0) or 0.0,
                offset_z=getattr(args, "import_offset_z", 0.0) or 0.0,
                thickness_override=getattr(args, "thickness_override", None),
            )
            for w in mapping.warnings:
                logger.warning(w)
            geom_status = mapping.status
        except Exception as exc:
            logger.warning(f"IMOD mapping failed: {exc}; continuing without projection matrices")
    else:
        if not ali_stack:
            logger.warning("No aligned stack provided; cannot build projection matrices")

    # Validate frame counts
    if ali_hdr and tlt_file and os.path.isfile(tlt_file):
        tlt_df = read_tlt(tlt_file)
        if ali_hdr.nz != len(tlt_df):
            logger.warning(
                f"ali stack nz={ali_hdr.nz} != .tlt rows={len(tlt_df)}"
            )
    if ali_hdr and xf_file and os.path.isfile(xf_file):
        xf_df = read_xf(xf_file)
        if ali_hdr.nz != len(xf_df):
            logger.warning(
                f"ali stack nz={ali_hdr.nz} != .xf rows={len(xf_df)}"
            )

    # ------------------------------------------------------------------ #
    # Tomogram size in bin-1 ali pixels (rlnTomoSizeX/Y/Z)
    # ------------------------------------------------------------------ #
    if rec_hdr:
        tomo_size_x = int(round(rec_hdr.nx * B_rec_ali))
        tomo_size_y = int(round(rec_hdr.ny * B_rec_ali))
        tomo_size_z = int(round(rec_hdr.nz * B_rec_ali))
    elif ali_hdr:
        tomo_size_x = ali_hdr.nx
        tomo_size_y = ali_hdr.ny
        tomo_size_z = ali_hdr.nz
        logger.warning("rec_tomo missing; using ali stack dimensions as tomo size")
    else:
        tomo_size_x = tomo_size_y = tomo_size_z = 0
        logger.warning("No ali or rec MRC available; tomo size set to 0")

    n_frames = ali_hdr.nz if ali_hdr else 0
    micrograph_ref = getattr(args, "micrograph_reference", "ali")

    # ------------------------------------------------------------------ #
    # CTF detection
    # ------------------------------------------------------------------ #
    ctf_source_opt = getattr(args, "ctf_source", "auto")
    ctf_source, ctf_path = "none", None
    if ctf_source_opt != "none":
        ctf_source, ctf_path = detect_ctf_source(
            imod_dir, basename,
            explicit_ctfplotter=getattr(args, "ctfplotter_file", None),
            explicit_defocus=None,
            explicit_ctffind=getattr(args, "ctffind_star", None),
            ctfplotter_info=getattr(args, "ctfplotter_info", None),
            ctfplotter_log=getattr(args, "ctfplotter_log", None),
        )
    ctf_premultiplied = not getattr(args, "do_not_premultiply", True)
    if ctf_premultiplied:
        logger.warning("CTF premultiplication requested; pixels will NOT be modified during import")

    # ------------------------------------------------------------------ #
    # Write tilt_series/<tomo>.star
    # ------------------------------------------------------------------ #
    tilt_star_path = pw.tilt_series_star_path(project_root, tomo_name)
    _write_tilt_series_star(
        tilt_star_path, tomo_name, ali_stack, raw_stack, mapping,
        tlt_file, xf_file, n_frames, micrograph_ref, project_root,
        relative, B_rec_ali, a_ali,
    )

    # ------------------------------------------------------------------ #
    # Write imod_settings/*.star
    # ------------------------------------------------------------------ #
    imod_dir_out = pw.imod_settings_dir(project_root, tomo_name)
    _write_imod_settings(
        imod_dir_out, tomo_name, imod_dir, ali_stack, raw_stack, rec_tomo,
        xf_file, tlt_file, newst_com, tilt_com, mapping, ctf_path,
    )

    # ------------------------------------------------------------------ #
    # Write / update tomograms.star
    # ------------------------------------------------------------------ #
    tomograms_star = pw.tomograms_star_path(project_root)
    tilt_star_rel = pw.store_path(tilt_star_path, project_root, relative)
    rec_tomo_stor = pw.store_path(rec_tomo, project_root, relative) if rec_tomo else "?"

    global_row = {
        "_rlnTomoName": tomo_name,
        "_rlnTomoTiltSeriesStarFile": tilt_star_rel,
        "_rlnTomoSizeX": tomo_size_x,
        "_rlnTomoSizeY": tomo_size_y,
        "_rlnTomoSizeZ": tomo_size_z,
        "_rlnTomoHand": "?",
        "_rlnMicrographOriginalPixelSize": f"{a_raw:.4f}",
        "_rlnTomoTiltSeriesPixelSize": f"{a_ali:.4f}",
        "_rlnTomoTomogramBinning": f"{B_rec_ali:.6f}",
        "_rlnTomoReconstructedTomogram": rec_tomo_stor,
        "_rlnTomoDenoisedTomogram": "?",
        "_rlnTomoFrameCount": n_frames,
        "_rlnEtomoDirectiveFile": "?",
        "_rlnOpticsGroupName": f"opticsGroup{tomo_name}",
        "_rlnVoltage": "?",
        "_rlnSphericalAberration": "?",
        "_rlnAmplitudeContrast": "?",
    }
    source_row = {
        "_tomoJANASTomoName": tomo_name,
        "_tomoJANASImodDir": imod_dir,
        "_tomoJANASRawStack": pw.store_path(raw_stack, project_root, relative) if raw_stack else "?",
        "_tomoJANASMdocFile": pw.store_path(mdoc_file, project_root, relative) if mdoc_file else "?",
        "_tomoJANASAliStack": pw.store_path(ali_stack, project_root, relative) if ali_stack else "?",
        "_tomoJANASRecTomogram": rec_tomo_stor,
        "_tomoJANASXfFile": pw.store_path(xf_file, project_root, relative) if xf_file else "?",
        "_tomoJANASTltFile": pw.store_path(tlt_file, project_root, relative) if tlt_file else "?",
        "_tomoJANASNewstCom": pw.store_path(newst_com, project_root, relative) if newst_com else "?",
        "_tomoJANASTiltCom": pw.store_path(tilt_com, project_root, relative) if tilt_com else "?",
        "_tomoJANASMicrographReference": micrograph_ref,
        "_tomoJANASCtfSource": ctf_source,
        "_tomoJANASCtfPremultiplied": "1" if ctf_premultiplied else "0",
        "_tomoJANASImportOffsetX": str(getattr(args, "import_offset_x", 0.0) or 0.0),
        "_tomoJANASImportOffsetY": str(getattr(args, "import_offset_y", 0.0) or 0.0),
        "_tomoJANASImportOffsetZ": str(getattr(args, "import_offset_z", 0.0) or 0.0),
        "_tomoJANASRelionGeometryStatus": geom_status,
    }

    _append_or_create_tomograms_star(tomograms_star, tomo_name, global_row, source_row)

    # ------------------------------------------------------------------ #
    # Write optimisation_set.star
    # ------------------------------------------------------------------ #
    opt_star = pw.optimisation_set_star_path(project_root)
    particles_star_rel = pw.store_path(pw.particles_all_star_path(project_root), project_root, relative)
    tomograms_star_rel = pw.store_path(tomograms_star, project_root, relative)
    write_star(opt_star, [LoopBlock(
        name="optimisation_set",
        columns=OPTIMISATION_SET_COLUMNS,
        rows=[[tomograms_star_rel, particles_star_rel]],
    )])

    # ------------------------------------------------------------------ #
    # Empty particles_all.star if it doesn't exist
    # ------------------------------------------------------------------ #
    if not os.path.isfile(proj.particles_all_star):
        _write_empty_particles_all(proj.particles_all_star)

    # ------------------------------------------------------------------ #
    # Manifest
    # ------------------------------------------------------------------ #
    manifest_tomos = []
    manifest_json_path = os.path.join(project_root, "project_manifest.json")
    import json
    existing_manifest = {}
    if os.path.isfile(manifest_json_path):
        try:
            with open(manifest_json_path) as f:
                existing_manifest = json.load(f)
        except Exception:
            pass
    existing_tomos = existing_manifest.get("tomograms", []) or []
    existing_tomos = [t for t in existing_tomos if t.get("tomo_name") != tomo_name]
    existing_tomos.append({
        "tomo_name": tomo_name,
        "imod_dir": imod_dir,
        "rec_tomo": rec_tomo,
        "ali_stack": ali_stack,
        "geometry_status": geom_status,
    })
    from tomojanas.io.project_writer import write_manifest
    write_manifest(project_root, {
        "tomojanas_version": get_version(),
        "created": existing_manifest.get("created", now_iso()),
        "updated": now_iso(),
        "tomograms": existing_tomos,
    })

    generated = [tilt_star_path, tomograms_star, opt_star]
    logger.info(f"geometry_status = {geom_status}")
    logger.invocation("tomojanas-import imod [complete]",
                      generated_files=generated)
    logger.summary()

    if getattr(args, "validate", False):
        from tomojanas.importers.validators import validate_project
        report = validate_project(proj, tomo_name=tomo_name,
                                  strict=getattr(args, "strict", False))
        from tomojanas.io.logs import write_validation_logs
        write_validation_logs(project_root, report)
        if not report["ok"] and getattr(args, "strict", False):
            return 1

    print(f"[tomojanas-import imod] project: {project_root}")
    print(f"  tomogram: {tomo_name}")
    print(f"  geometry: {geom_status}")
    print(f"  tilts:    {n_frames}")
    if mapping:
        print(f"  excluded: {sum(mapping.excluded)}")
    return 0


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #

def _write_tilt_series_star(
    path, tomo_name, ali_stack, raw_stack, mapping,
    tlt_file, xf_file, n_frames, micrograph_ref, project_root,
    relative, B_rec_ali, a_ali,
):
    ali_rel = pw.store_path(ali_stack, project_root, relative) if ali_stack else "?"
    raw_rel = pw.store_path(raw_stack, project_root, relative) if raw_stack else "?"

    # Load tilt angles and xf if available
    tilt_angles = []
    if tlt_file and os.path.isfile(tlt_file):
        from tomojanas.metadata.imod_tlt import read_tlt_angles
        tilt_angles = read_tlt_angles(tlt_file)

    xf_rows = {}
    if xf_file and os.path.isfile(xf_file):
        from tomojanas.metadata.imod_xf import read_xf
        xf_df = read_xf(xf_file)
        for _, r in xf_df.iterrows():
            xf_rows[int(r["tilt_index"])] = r

    n_tot = n_frames or (len(tilt_angles) if tilt_angles else 0)
    if n_tot == 0 and mapping:
        n_tot = len(mapping.tilt_angles_all)

    relion_rows = []
    for i in range(n_tot):
        angle = tilt_angles[i] if i < len(tilt_angles) else "?"
        # 1-based STAR index
        mic_name = f"{i + 1:06d}@{ali_rel}" if micrograph_ref == "ali" else f"{i + 1:06d}@{raw_rel}"
        row = {
            "_rlnMicrographName": mic_name,
            "_rlnMicrographMovieName": "?",
            "_rlnMicrographMetadata": "?",
            "_rlnTomoTiltMovieFrameCount": "?",
            "_rlnTomoTiltMovieIndex": i + 1,
            "_rlnTomoNominalStageTiltAngle": f"{angle:.4f}" if isinstance(angle, float) else "?",
            "_rlnTomoNominalTiltAxisAngle": "?",
            "_rlnTomoNominalDefocus": "?",
            "_rlnMicrographPreExposure": "?",
            "_rlnDefocusU": "?",
            "_rlnDefocusV": "?",
            "_rlnDefocusAngle": "?",
            "_rlnCtfAstigmatism": "?",
            "_rlnCtfFigureOfMerit": "?",
            "_rlnCtfMaxResolution": "?",
            "_rlnCtfScalefactor": "?",
            "_rlnPhaseShift": "?",
            "_rlnTomoXTilt": "?",
            "_rlnTomoYTilt": f"{angle:.4f}" if isinstance(angle, float) else "?",
            "_rlnTomoZRot": "?",
            "_rlnTomoXShiftAngst": "?",
            "_rlnTomoYShiftAngst": "?",
            # Projection matrices stored in tomoJANAS block; RELION slots left empty
            "_rlnTomoProjX": "?",
            "_rlnTomoProjY": "?",
            "_rlnTomoProjZ": "?",
            "_rlnTomoProjW": "?",
        }
        relion_rows.append(row)

    relion_block = LoopBlock(
        name=tomo_name,
        columns=TILT_SERIES_COLUMNS,
        rows=relion_rows,
    )

    # tomoJANAS tilt mapping block
    tilt_map_rows = []
    excluded_set = set(
        mapping.excluded.index(True) for _ in range(sum(mapping.excluded))
    ) if mapping else set()
    if mapping:
        excluded_indices = {i for i, e in enumerate(mapping.excluded) if e}
    else:
        excluded_indices = set()

    for i in range(n_tot):
        xr = xf_rows.get(i + 1, {})
        tilt_map_rows.append({
            "_tomoJANASTiltIndex": i + 1,
            "_tomoJANASAliStack": ali_rel,
            "_tomoJANASAliSlice": i + 1,
            "_tomoJANASRawStack": raw_rel,
            "_tomoJANASRawSlice": i + 1,
            "_tomoJANASXfA11": f"{float(xr.get('a11', 1.0)):.6f}",
            "_tomoJANASXfA12": f"{float(xr.get('a12', 0.0)):.6f}",
            "_tomoJANASXfA21": f"{float(xr.get('a21', 0.0)):.6f}",
            "_tomoJANASXfA22": f"{float(xr.get('a22', 1.0)):.6f}",
            "_tomoJANASXfDX": f"{float(xr.get('dx', 0.0)):.6f}",
            "_tomoJANASXfDY": f"{float(xr.get('dy', 0.0)):.6f}",
            "_tomoJANASXfDirection": "raw_to_ali",
            "_tomoJANASXfDirectionStatus": "inferred_from_newst",
        })
    tilt_map_block = LoopBlock(name="tomoJANAS_tilt_mapping",
                                columns=TILT_MAPPING_COLUMNS,
                                rows=tilt_map_rows)

    blocks = [relion_block, tilt_map_block]

    # If mapping built, store matrices in tomoJANAS scalar block
    if mapping and mapping.projections:
        mat_rows = []
        for f_idx, M in enumerate(mapping.projections):
            orig_idx = mapping.old_frame_index[f_idx]
            row = {"_tomoJANASTiltIndex": orig_idx + 1}
            for r in range(4):
                for c in range(4):
                    row[f"_tomoJANASProj{r}{c}"] = f"{M[r, c]:.8f}"
            mat_rows.append(row)
        blocks.append(LoopBlock(
            name="tomoJANAS_projection_matrices",
            columns=PROJECTION_MATRIX_COLUMNS,
            rows=mat_rows,
        ))

    write_star(path, blocks)


def _write_imod_settings(
    imod_dir_out, tomo_name, imod_dir, ali_stack, raw_stack, rec_tomo,
    xf_file, tlt_file, newst_com, tilt_com, mapping, ctf_path,
):
    os.makedirs(imod_dir_out, exist_ok=True)

    # imod_files.star
    write_star(os.path.join(imod_dir_out, "imod_files.star"), [
        LoopBlock(name="tomoJANAS_imod_files", columns=[
            "_tomoJANASTomoName", "_tomoJANASImodDir",
            "_tomoJANASRawStack", "_tomoJANASAliStack", "_tomoJANASRecTomogram",
            "_tomoJANASXfFile", "_tomoJANASTltFile", "_tomoJANASNewstCom", "_tomoJANASTiltCom",
        ], rows=[[
            tomo_name, imod_dir,
            raw_stack or "?", ali_stack or "?", rec_tomo or "?",
            xf_file or "?", tlt_file or "?", newst_com or "?", tilt_com or "?",
        ]])
    ])

    # imod_geometry.star
    geom_pairs: dict = {
        "_tomoJANASTomoName": tomo_name,
        "_tomoJANASRelionGeometryStatus": mapping.status if mapping else RelionGeometryStatus.MATRICES_MISSING,
    }
    if mapping:
        geom_pairs["_tomoJANASMappingW"] = mapping.w
        geom_pairs["_tomoJANASMappingH"] = mapping.h
        geom_pairs["_tomoJANASMappingD"] = mapping.d
        geom_pairs["_tomoJANASFramesMissing"] = 1 if mapping.frames_missing else 0
        if mapping.warnings:
            geom_pairs["_tomoJANASMappingWarnings"] = "; ".join(mapping.warnings)
    write_star(os.path.join(imod_dir_out, "imod_geometry.star"), [
        PairBlock(name="tomoJANAS_imod_geometry", pairs=geom_pairs)
    ])

    # relion_geometry_status.star
    write_star(os.path.join(imod_dir_out, "relion_geometry_status.star"), [
        PairBlock(name="tomoJANAS_relion_status", pairs={
            "_tomoJANASTomoName": tomo_name,
            "_tomoJANASRelionGeometryStatus": mapping.status if mapping else RelionGeometryStatus.MATRICES_MISSING,
            "_tomoJANASRelionOracleStatus": "not_run",
        })
    ])

    # xf_matrices.star (if xf available)
    if xf_file and os.path.isfile(xf_file):
        xf_df = read_xf(xf_file)
        cols = ["_tomoJANASTiltIndex", "_tomoJANASXfA11", "_tomoJANASXfA12",
                "_tomoJANASXfA21", "_tomoJANASXfA22", "_tomoJANASXfDX", "_tomoJANASXfDY"]
        rows = [[int(r.tilt_index), float(r.a11), float(r.a12),
                  float(r.a21), float(r.a22), float(r.dx), float(r.dy)]
                 for _, r in xf_df.iterrows()]
        write_star(os.path.join(imod_dir_out, "xf_matrices.star"), [
            LoopBlock(name="tomoJANAS_xf_matrices", columns=cols, rows=rows)
        ])

    # ctf_settings.star
    write_star(os.path.join(imod_dir_out, "ctf_settings.star"), [
        PairBlock(name="tomoJANAS_ctf_settings", pairs={
            "_tomoJANASCtfSource": "none" if not ctf_path else "ctfplotter",
            "_tomoJANASCtfFile": ctf_path or "?",
            "_tomoJANASCtfPremultiplied": "0",
        })
    ])

    # coordinate_inference.star (empty placeholder for particle import to fill)
    write_star(os.path.join(imod_dir_out, "coordinate_inference.star"), [
        PairBlock(name="tomoJANAS_coordinate_inference", pairs={
            "_tomoJANASInferenceStatus": "pending",
        })
    ])

    # imod_import.star (full provenance)
    write_star(os.path.join(imod_dir_out, "imod_import.star"), [
        PairBlock(name="tomoJANAS_imod_import", pairs={
            "_tomoJANASImportTimestamp": now_iso(),
            "_tomoJANASTomoName": tomo_name,
            "_tomoJANASImodDir": imod_dir,
        })
    ])


def _append_or_create_tomograms_star(path, tomo_name, global_row, source_row):
    """Append (or create) a tomogram row in tomograms.star."""
    global_rows = []
    source_rows = []

    if os.path.isfile(path):
        try:
            blks = read_star(path)
            if "global" in blks and blks["global"]["type"] == "loop":
                df = blks["global"]["df"]
                for _, r in df.iterrows():
                    if r.get("_rlnTomoName") == tomo_name:
                        continue  # replace existing
                    global_rows.append(dict(r))
            if "tomoJANAS_tomogram_sources" in blks and blks["tomoJANAS_tomogram_sources"]["type"] == "loop":
                df2 = blks["tomoJANAS_tomogram_sources"]["df"]
                for _, r in df2.iterrows():
                    if r.get("_tomoJANASTomoName") == tomo_name:
                        continue
                    source_rows.append(dict(r))
        except Exception:
            pass

    global_rows.append(global_row)
    source_rows.append(source_row)

    write_star(path, [
        LoopBlock(name="global", columns=TOMOGRAMS_GLOBAL_COLUMNS, rows=global_rows),
        LoopBlock(name="tomoJANAS_tomogram_sources", columns=TOMOGRAMS_SOURCE_COLUMNS, rows=source_rows),
    ])


def _write_empty_particles_all(path):
    from tomojanas.metadata.relion_labels import PARTICLE_OPTICS_COLUMNS, PARTICLES_COLUMNS
    write_star(path, [
        LoopBlock(name="optics", columns=PARTICLE_OPTICS_COLUMNS, rows=[]),
        LoopBlock(name="particles", columns=PARTICLES_COLUMNS, rows=[]),
    ])
