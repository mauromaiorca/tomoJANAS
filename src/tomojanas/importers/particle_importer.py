#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/importers/particle_importer.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
``tomojanas-import particles`` — import picked coordinates into a project.

Writes:
    <project>/particles_all.star          (data_optics + data_particles)
    <project>/tilt_series/<tomo>/individual_particles/P*.star
"""
from __future__ import annotations

import math
import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from tomojanas import get_version
from tomojanas.io.mrc import read_mrc_header
from tomojanas.io.star_writer import LoopBlock, write_star
from tomojanas.io.star_reader import read_star
from tomojanas.io.logs import ImportLogger, now_iso
from tomojanas.io import project_writer as pw
from tomojanas.models.project import Project, ProjectMissingError
from tomojanas.metadata.relion_labels import (
    PARTICLE_OPTICS_COLUMNS, PARTICLES_COLUMNS,
    PARTICLE_SOURCE_COLUMNS, PARTICLE_ROI_COLUMNS, PARTICLE_PROJECTION_COLUMNS,
    ProjectionStatus,
)
from tomojanas.geometry.coordinates import (
    rec_voxel_to_relion_centered_angst,
    normalise_coordinate_order,
    coordinate_roundtrip_error,
)
from tomojanas.geometry.roi import (
    roi_radius_angst_from_args,
    projection_radius_px,
    storage_box_from_radius,
    circle_inside_frame,
    sphere_inside_volume,
)
from tomojanas.geometry.imod_mapping import project_particle_to_ali, invert_xf_to_raw

__all__ = ["import_particles"]


# --------------------------------------------------------------------------- #
# Coordinate readers
# --------------------------------------------------------------------------- #

def _read_coords_csv(path: str, axis_order: str) -> List[Tuple[float, float, float]]:
    """CSV with x,y,z or X,Y,Z or (if axis_order==zyx) z,y,x columns."""
    df = pd.read_csv(path, sep=None, engine="python")
    df.columns = [c.strip() for c in df.columns]
    lower_map = {c.lower(): c for c in df.columns}
    ao = axis_order.lower()
    if ao == "zyx" and "z" in lower_map:
        z = df[lower_map["z"]].values.astype(float)
        y = df[lower_map["y"]].values.astype(float)
        x = df[lower_map["x"]].values.astype(float)
        return list(zip(x.tolist(), y.tolist(), z.tolist()))
    if "x" in lower_map:
        x = df[lower_map["x"]].values.astype(float)
        y = df[lower_map["y"]].values.astype(float)
        z = df[lower_map["z"]].values.astype(float)
        return list(zip(x.tolist(), y.tolist(), z.tolist()))
    # Try first three columns
    cols = list(df.columns)[:3]
    x = df[cols[0]].values.astype(float)
    y = df[cols[1]].values.astype(float)
    z = df[cols[2]].values.astype(float)
    return list(zip(x.tolist(), y.tolist(), z.tolist()))


def _read_coords_napari(path: str) -> List[Tuple[float, float, float]]:
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    if "axis-0" in df.columns:
        z = df["axis-0"].values.astype(float)
        y = df["axis-1"].values.astype(float)
        x = df["axis-2"].values.astype(float)
        return list(zip(x.tolist(), y.tolist(), z.tolist()))
    lower_map = {c.lower(): c for c in df.columns}
    if "z" in lower_map:
        z = df[lower_map["z"]].values.astype(float)
        y = df[lower_map["y"]].values.astype(float)
        x = df[lower_map["x"]].values.astype(float)
        return list(zip(x.tolist(), y.tolist(), z.tolist()))
    raise ValueError(f"Cannot parse napari CSV columns: {list(df.columns)}")


def _read_single_point(s: str, axis_order: str) -> List[Tuple[float, float, float]]:
    """Parse "X,Y,Z" or "X Y Z" string and apply axis_order normalisation."""
    s = s.replace(",", " ")
    vals = [float(v) for v in s.split()]
    if len(vals) < 3:
        raise ValueError(f"--input-single-point needs 3 values, got: {s!r}")
    x, y, z = normalise_coordinate_order(vals, axis_order)
    return [(x, y, z)]


# --------------------------------------------------------------------------- #
# Main function
# --------------------------------------------------------------------------- #

def import_particles(args) -> int:
    project_root = os.path.abspath(args.project)
    tomo_name = args.tomo_name

    logger = ImportLogger(project_root, get_version())
    logger.invocation("tomojanas-import particles", argv=sys.argv)

    # Open project (must exist already)
    try:
        proj = Project.open_or_create(
            project_root,
            create_if_missing=False,
            fail_if_existing=False,
            append=True,
        )
    except (ProjectMissingError, Exception) as exc:
        logger.error(str(exc))
        return 1

    # ------------------------------------------------------------------ #
    # Load tomogram metadata from tomograms.star
    # ------------------------------------------------------------------ #
    tomo_meta = _load_tomo_metadata(proj.tomograms_star, tomo_name)
    if tomo_meta is None:
        logger.error(f"Tomogram '{tomo_name}' not found in {proj.tomograms_star}. "
                     f"Run tomojanas-import imod first.")
        return 1

    a_ali = float(tomo_meta.get("_rlnTomoTiltSeriesPixelSize", 1.0) or 1.0)
    B_rec_ali = float(tomo_meta.get("_rlnTomoTomogramBinning", 1.0) or 1.0)
    tomo_size_x = int(float(tomo_meta.get("_rlnTomoSizeX", 0) or 0))
    tomo_size_y = int(float(tomo_meta.get("_rlnTomoSizeY", 0) or 0))
    tomo_size_z = int(float(tomo_meta.get("_rlnTomoSizeZ", 0) or 0))
    voltage = float(tomo_meta.get("_rlnVoltage", 300.0) or 300.0) if tomo_meta.get("_rlnVoltage", "?") != "?" else 300.0
    cs = float(tomo_meta.get("_rlnSphericalAberration", 2.7) or 2.7) if tomo_meta.get("_rlnSphericalAberration", "?") != "?" else 2.7
    amp_contrast = float(tomo_meta.get("_rlnAmplitudeContrast", 0.1) or 0.1) if tomo_meta.get("_rlnAmplitudeContrast", "?") != "?" else 0.1

    # ------------------------------------------------------------------ #
    # Load tilt series STAR (projection matrices + xf)
    # ------------------------------------------------------------------ #
    tilt_star = pw.tilt_series_star_path(project_root, tomo_name)
    proj_matrices, xf_df_tomo, tilt_angles, ali_stack_path = _load_tilt_series_data(tilt_star, tomo_name)

    # ------------------------------------------------------------------ #
    # Load MRC headers
    # ------------------------------------------------------------------ #
    rec_tomo_path_str = tomo_meta.get("_rlnTomoReconstructedTomogram", "?") or "?"
    rec_hdr = None
    ali_hdr = None
    if rec_tomo_path_str not in ("?", "", None) and os.path.isfile(rec_tomo_path_str):
        try:
            rec_hdr = read_mrc_header(rec_tomo_path_str)
        except Exception:
            pass
    # Also resolve relative to project root
    if rec_hdr is None and rec_tomo_path_str not in ("?", "", None):
        cand = os.path.join(project_root, rec_tomo_path_str)
        if os.path.isfile(cand):
            try:
                rec_hdr = read_mrc_header(cand)
            except Exception:
                pass
    if ali_stack_path and os.path.isfile(ali_stack_path):
        try:
            ali_hdr = read_mrc_header(ali_stack_path)
        except Exception:
            pass

    ali_nx = ali_hdr.nx if ali_hdr else tomo_size_x
    ali_ny = ali_hdr.ny if ali_hdr else tomo_size_y

    # ------------------------------------------------------------------ #
    # Read coordinates
    # ------------------------------------------------------------------ #
    input_path = getattr(args, "input", None)
    input_single = getattr(args, "input_single_point", None)
    fmt = getattr(args, "format", "csv")
    axis_order = getattr(args, "axis_order", "auto") or "auto"
    if axis_order == "auto":
        axis_order = "xyz"

    coords_raw: List[Tuple[float, float, float]] = []
    if input_single:
        coords_raw = _read_single_point(input_single, axis_order)
    elif input_path:
        if fmt in ("napari",):
            coords_raw = _read_coords_napari(input_path)
        else:
            coords_raw = _read_coords_csv(input_path, axis_order)
    else:
        logger.error("No --input or --input-single-point provided")
        return 1

    coord_system = getattr(args, "coordinate_system", "rec-voxel") or "rec-voxel"
    if coord_system == "auto":
        coord_system = "rec-voxel"
    indexing = getattr(args, "indexing", "zero-based") or "zero-based"
    if indexing == "auto":
        indexing = "zero-based"

    logger.info(f"Read {len(coords_raw)} coordinates (system={coord_system}, indexing={indexing})")

    # ------------------------------------------------------------------ #
    # ROI
    # ------------------------------------------------------------------ #
    try:
        roi_radius_angst = roi_radius_angst_from_args(args, a_ali)
    except ValueError:
        roi_radius_angst = 0.0
        logger.warning("No ROI radius specified; ROI will be zero")
    roi_radius_ali_px = projection_radius_px(roi_radius_angst, a_ali)
    roi_storage_box_px = storage_box_from_radius(roi_radius_ali_px, padding=float(getattr(args, "roi_padding_angst", 0.0) or 0.0) / a_ali)
    roi_radius_vox = roi_radius_angst / a_ali / B_rec_ali if (a_ali > 0 and B_rec_ali > 0) else 0.0

    # ------------------------------------------------------------------ #
    # Particle naming
    # ------------------------------------------------------------------ #
    particle_prefix = getattr(args, "particle_prefix", "P") or "P"
    particle_name_override = getattr(args, "particle_name", None)
    start_id = int(getattr(args, "start_particle_id", 1) or 1)

    # Find next free ID
    ip_dir = pw.individual_particles_dir(project_root, tomo_name)
    os.makedirs(ip_dir, exist_ok=True)
    next_id = start_id
    if len(coords_raw) == 1 and particle_name_override:
        particle_names = [particle_name_override]
    else:
        particle_names = []
        for i in range(len(coords_raw)):
            particle_names.append(f"{particle_prefix}{next_id + i:06d}")

    # ------------------------------------------------------------------ #
    # Optics group
    # ------------------------------------------------------------------ #
    optics_group_id = 1
    optics_group_name = f"opticsGroup{tomo_name}"
    optics_row = {
        "_rlnOpticsGroup": optics_group_id,
        "_rlnOpticsGroupName": optics_group_name,
        "_rlnTomoTiltSeriesPixelSize": f"{a_ali:.4f}",
        "_rlnVoltage": f"{voltage:.1f}",
        "_rlnSphericalAberration": f"{cs:.2f}",
        "_rlnAmplitudeContrast": f"{amp_contrast:.3f}",
    }

    # ------------------------------------------------------------------ #
    # Process each particle
    # ------------------------------------------------------------------ #
    all_particle_rows: List[Dict] = []
    generated = []

    for idx, (x_in, y_in, z_in) in enumerate(coords_raw):
        pname = particle_names[idx]
        particle_id = next_id + idx

        # --- coordinate conversion ---
        if coord_system == "rec-voxel":
            x_ang, y_ang, z_ang = rec_voxel_to_relion_centered_angst(
                x_in, y_in, z_in,
                tomo_size_x, tomo_size_y, tomo_size_z,
                a_ali, B_rec_ali, indexing,
            )
        elif coord_system == "relion-centered-angst":
            x_ang, y_ang, z_ang = x_in, y_in, z_in
        else:
            # fallback: treat as rec-voxel
            x_ang, y_ang, z_ang = rec_voxel_to_relion_centered_angst(
                x_in, y_in, z_in,
                tomo_size_x, tomo_size_y, tomo_size_z,
                a_ali, B_rec_ali, indexing,
            )

        particle_row = {
            "_rlnTomoName": tomo_name,
            "_rlnTomoParticleName": pname,
            "_rlnTomoParticleId": particle_id,
            "_rlnCenteredCoordinateXAngst": f"{x_ang:.4f}",
            "_rlnCenteredCoordinateYAngst": f"{y_ang:.4f}",
            "_rlnCenteredCoordinateZAngst": f"{z_ang:.4f}",
            "_rlnOriginXAngst": "0.0000",
            "_rlnOriginYAngst": "0.0000",
            "_rlnOriginZAngst": "0.0000",
            "_rlnTomoSubtomogramRot": "0.0000",
            "_rlnTomoSubtomogramTilt": "0.0000",
            "_rlnTomoSubtomogramPsi": "0.0000",
            "_rlnAngleRot": "0.0000",
            "_rlnAngleTilt": "0.0000",
            "_rlnAnglePsi": "0.0000",
            "_rlnOpticsGroup": optics_group_id,
        }
        all_particle_rows.append(particle_row)

        # --- per-tilt projections ---
        proj_rows = _compute_projections(
            x_in, y_in, z_in, B_rec_ali,
            proj_matrices, xf_df_tomo, tilt_angles,
            ali_stack_path, ali_nx, ali_ny, roi_radius_ali_px,
        )

        # --- individual P*.star ---
        if getattr(args, "write_individual_particles", True):
            p_star_path = os.path.join(ip_dir, f"{pname}.star")
            _write_particle_star(
                p_star_path,
                optics_row, particle_row,
                x_in, y_in, z_in, coord_system, indexing, axis_order,
                roi_radius_angst, roi_storage_box_px, roi_radius_ali_px,
                proj_rows,
            )
            generated.append(p_star_path)

    # ------------------------------------------------------------------ #
    # Update particles_all.star
    # ------------------------------------------------------------------ #
    if getattr(args, "update_particles_all", True):
        _update_particles_all(
            proj.particles_all_star, optics_row, all_particle_rows, tomo_name
        )
        generated.append(proj.particles_all_star)

    logger.info(f"Wrote {len(coords_raw)} particles for {tomo_name}")
    logger.invocation("tomojanas-import particles [complete]", generated_files=generated)
    logger.summary()

    if getattr(args, "validate", False):
        from tomojanas.importers.validators import validate_project
        from tomojanas.io.logs import write_validation_logs
        report = validate_project(proj, tomo_name=tomo_name,
                                  strict=getattr(args, "strict", False))
        write_validation_logs(project_root, report)
        if not report["ok"] and getattr(args, "strict", False):
            return 1

    print(f"[tomojanas-import particles] {len(coords_raw)} particle(s) added to '{tomo_name}'")
    return 0


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def _load_tomo_metadata(tomograms_star: str, tomo_name: str) -> Optional[Dict]:
    if not os.path.isfile(tomograms_star):
        return None
    try:
        blks = read_star(tomograms_star)
        if "global" not in blks:
            return None
        df = blks["global"]["df"]
        for _, r in df.iterrows():
            if str(r.get("_rlnTomoName", "")) == tomo_name:
                return dict(r)
    except Exception:
        pass
    return None


def _load_tilt_series_data(tilt_star: str, tomo_name: str):
    """Returns (proj_matrices, xf_df, tilt_angles, ali_stack_path)."""
    proj_matrices = []   # list of (4,4) ndarray, good frames in order
    xf_df_tomo = None
    tilt_angles = []
    ali_stack_path = None

    if not os.path.isfile(tilt_star):
        return proj_matrices, xf_df_tomo, tilt_angles, ali_stack_path

    try:
        blks = read_star(tilt_star)

        # Tilt angles from RELION block
        if tomo_name in blks and blks[tomo_name]["type"] == "loop":
            df = blks[tomo_name]["df"]
            for _, r in df.iterrows():
                ang = r.get("_rlnTomoNominalStageTiltAngle", "?")
                try:
                    tilt_angles.append(float(ang))
                except (TypeError, ValueError):
                    tilt_angles.append(0.0)
                # Get ali stack path from first row
                if ali_stack_path is None:
                    mic = str(r.get("_rlnMicrographName", "?"))
                    if "@" in mic:
                        ali_stack_path = mic.split("@", 1)[1]

        # Projection matrices from tomoJANAS block
        if "tomoJANAS_projection_matrices" in blks and blks["tomoJANAS_projection_matrices"]["type"] == "loop":
            mat_df = blks["tomoJANAS_projection_matrices"]["df"]
            for _, r in mat_df.iterrows():
                M = np.eye(4, dtype=np.float64)
                for row_i in range(4):
                    for col_i in range(4):
                        key = f"_tomoJANASProj{row_i}{col_i}"
                        try:
                            M[row_i, col_i] = float(r.get(key, 1.0 if row_i == col_i else 0.0))
                        except (TypeError, ValueError):
                            pass
                proj_matrices.append(M)

        # XF data from tilt mapping
        if "tomoJANAS_tilt_mapping" in blks and blks["tomoJANAS_tilt_mapping"]["type"] == "loop":
            tm_df = blks["tomoJANAS_tilt_mapping"]["df"]
            from tomojanas.metadata.imod_xf import _FIELDS
            rows_for_xf = []
            for _, r in tm_df.iterrows():
                rows_for_xf.append({
                    "tilt_index": int(float(r.get("_tomoJANASTiltIndex", 0))),
                    "a11": float(r.get("_tomoJANASXfA11", 1.0)),
                    "a12": float(r.get("_tomoJANASXfA12", 0.0)),
                    "a21": float(r.get("_tomoJANASXfA21", 0.0)),
                    "a22": float(r.get("_tomoJANASXfA22", 1.0)),
                    "dx": float(r.get("_tomoJANASXfDX", 0.0)),
                    "dy": float(r.get("_tomoJANASXfDY", 0.0)),
                })
            import pandas
            xf_df_tomo = pandas.DataFrame(rows_for_xf)
    except Exception:
        pass

    return proj_matrices, xf_df_tomo, tilt_angles, ali_stack_path


def _compute_projections(
    x_rec, y_rec, z_rec, B_rec_ali,
    proj_matrices, xf_df, tilt_angles,
    ali_stack, ali_nx, ali_ny, roi_radius_px,
) -> List[Dict]:
    """Compute per-tilt projection info for one particle."""
    rows = []
    n_tilts = max(len(proj_matrices), len(tilt_angles))

    for f in range(n_tilts):
        angle = tilt_angles[f] if f < len(tilt_angles) else 0.0
        ali_cx, ali_cy = None, None
        proj_status = ProjectionStatus.MISSING
        ali_inside = False
        raw_cx, raw_cy = None, None
        raw_inside = False
        visible = False

        if f < len(proj_matrices):
            M = proj_matrices[f]
            try:
                ali_cx, ali_cy = project_particle_to_ali(
                    float(x_rec), float(y_rec), float(z_rec), float(B_rec_ali), M
                )
                proj_status = ProjectionStatus.OK
                ali_inside = circle_inside_frame(ali_cx, ali_cy, roi_radius_px, ali_nx, ali_ny)
                visible = (
                    0 <= ali_cx <= ali_nx and 0 <= ali_cy <= ali_ny
                )
                # Inverse xf to raw
                if xf_df is not None and f < len(xf_df):
                    raw_cx, raw_cy = invert_xf_to_raw(ali_cx, ali_cy, xf_df, f)
                    raw_inside = circle_inside_frame(raw_cx, raw_cy, roi_radius_px, ali_nx, ali_ny)
            except Exception:
                proj_status = ProjectionStatus.MISSING

        rows.append({
            "_tomoJANASTiltIndex": f + 1,
            "_tomoJANASStageTiltAngle": f"{angle:.4f}",
            "_tomoJANASAlignedStack": ali_stack or "?",
            "_tomoJANASAlignedSlice": f + 1,
            "_tomoJANASAlignedCenterX": f"{ali_cx:.4f}" if ali_cx is not None else "?",
            "_tomoJANASAlignedCenterY": f"{ali_cy:.4f}" if ali_cy is not None else "?",
            "_tomoJANASAlignedRadiusPixel": f"{roi_radius_px:.4f}",
            "_tomoJANASAlignedCircleInsideFrame": "1" if ali_inside else "0",
            "_tomoJANASRawStack": "?",
            "_tomoJANASRawSlice": f + 1,
            "_tomoJANASRawCenterX": f"{raw_cx:.4f}" if raw_cx is not None else "?",
            "_tomoJANASRawCenterY": f"{raw_cy:.4f}" if raw_cy is not None else "?",
            "_tomoJANASRawRadiusPixel": f"{roi_radius_px:.4f}",
            "_tomoJANASRawCircleInsideFrame": "1" if raw_inside else "0",
            "_tomoJANASVisibleInTilt": "1" if visible else "0",
            "_tomoJANASProjectionStatus": proj_status,
            "_tomoJANASDefocusU": "?",
            "_tomoJANASDefocusV": "?",
            "_tomoJANASDefocusAngle": "?",
        })
    return rows


def _write_particle_star(
    path, optics_row, particle_row,
    x_in, y_in, z_in, coord_system, indexing, axis_order,
    roi_radius_angst, storage_box_px, proj_radius_px,
    proj_rows,
):
    pname = particle_row["_rlnTomoParticleName"]

    source_row = {
        "_tomoJANASParticleName": pname,
        "_tomoJANASPickedVolume": particle_row.get("_rlnTomoName", "?"),
        "_tomoJANASPickedCoordinateX": f"{x_in:.4f}",
        "_tomoJANASPickedCoordinateY": f"{y_in:.4f}",
        "_tomoJANASPickedCoordinateZ": f"{z_in:.4f}",
        "_tomoJANASPickedCoordinateSystem": coord_system,
        "_tomoJANASPickedIndexing": indexing,
        "_tomoJANASPickedAxisOrder": axis_order,
        "_tomoJANASPickedSoftware": "unknown",
    }
    roi_row = {
        "_tomoJANASParticleName": pname,
        "_tomoJANASRoiShape3D": "sphere",
        "_tomoJANASRoiShape2D": "circle",
        "_tomoJANASRoiRadiusAngst": f"{roi_radius_angst:.4f}",
        "_tomoJANASRoiDiameterAngst": f"{2 * roi_radius_angst:.4f}",
        "_tomoJANASRoiPaddingAngst": "0.0000",
        "_tomoJANASStorageBoxShape3D": "cube",
        "_tomoJANASStorageBoxSizeVoxel": str(storage_box_px),
        "_tomoJANASProjectionRadiusPixel": f"{proj_radius_px:.4f}",
        "_tomoJANASProjectionStorageBoxSizePixel": str(storage_box_px),
    }

    blocks = [
        LoopBlock(name="optics", columns=PARTICLE_OPTICS_COLUMNS, rows=[optics_row]),
        LoopBlock(name="particles", columns=PARTICLES_COLUMNS, rows=[particle_row]),
        LoopBlock(name="tomoJANAS_particle_source", columns=PARTICLE_SOURCE_COLUMNS, rows=[source_row]),
        LoopBlock(name="tomoJANAS_particle_roi", columns=PARTICLE_ROI_COLUMNS, rows=[roi_row]),
    ]
    if proj_rows:
        blocks.append(LoopBlock(name="tomoJANAS_particle_projections",
                                columns=PARTICLE_PROJECTION_COLUMNS, rows=proj_rows))
    write_star(path, blocks)


def _update_particles_all(path, optics_row, new_particle_rows, tomo_name):
    existing_optics: List[Dict] = []
    existing_particles: List[Dict] = []

    if os.path.isfile(path):
        try:
            blks = read_star(path)
            if "optics" in blks and blks["optics"]["type"] == "loop":
                df = blks["optics"]["df"]
                for _, r in df.iterrows():
                    existing_optics.append(dict(r))
            if "particles" in blks and blks["particles"]["type"] == "loop":
                df2 = blks["particles"]["df"]
                for _, r in df2.iterrows():
                    existing_particles.append(dict(r))
        except Exception:
            pass

    # Merge optics (add if not present for this tomo)
    og_name = optics_row["_rlnOpticsGroupName"]
    if not any(str(r.get("_rlnOpticsGroupName", "")) == og_name for r in existing_optics):
        existing_optics.append(optics_row)

    existing_particles.extend(new_particle_rows)

    write_star(path, [
        LoopBlock(name="optics", columns=PARTICLE_OPTICS_COLUMNS, rows=existing_optics),
        LoopBlock(name="particles", columns=PARTICLES_COLUMNS, rows=existing_particles),
    ])
