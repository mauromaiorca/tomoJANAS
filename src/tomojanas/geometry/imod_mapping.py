#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/geometry/imod_mapping.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
Python port of RELION's ``src/jaz/tomography/imod_import.cpp``.

Constructs ``worldToImage`` 4×4 matrices for each (non-excluded) tilt exactly
as RELION does, enabling tomoJANAS to compute per-tilt particle projections
with the same geometry as RELION.

Coordinate convention (no flip_yz)
-----------------------------------
The "world" coordinate is in bin-1 aligned-stack pixels:

    world_x  ↔  ali-image X (horizontal)
    world_y  ↔  reconstruction thickness / depth direction
    world_z  ↔  ali-image Y (vertical)

A particle at reconstruction voxel (x_rec, y_rec, z_rec) maps to world by:

    world_x = x_rec * B_rec_ali
    world_y = z_rec * B_rec_ali   (rec MRC Z = depth = IMOD Y)
    world_z = y_rec * B_rec_ali   (rec MRC Y = ali-image Y = world Z)

See project_particle_to_ali() for the full projection helper.

Geometry status
---------------
Matrices built by this algorithm are marked ``relion_imod_algorithm_ported``.
They are NOT marked ``relion_extract_ready`` unless external validation confirms
they match RELION's output within tolerance.
"""
from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from tomojanas.metadata.imod_com import parse_com_file, com_get, com_get_float, com_get_int, com_get_numbers, parse_exclude_list
from tomojanas.metadata.imod_tlt import read_tlt_angles
from tomojanas.metadata.imod_xf import read_xf
from tomojanas.io.mrc import read_mrc_header
from tomojanas.metadata.relion_labels import RelionGeometryStatus

__all__ = [
    "ImodMapping",
    "build_relion_imod_mapping",
    "project_particle_to_ali",
    "invert_xf_to_raw",
]


@dataclass
class ImodMapping:
    """Output of build_relion_imod_mapping: one 4×4 worldToImage per good tilt."""
    projections: List[np.ndarray]   # list of (4,4) float64 worldToImage matrices
    w: int                          # output tomogram width (bin-1 ali pixels)
    h: int                          # output tomogram height
    d: int                          # output tomogram depth/thickness
    frames_missing: bool
    old_frame_index: List[int]      # maps good-frame index → original frame index
    excluded: List[bool]            # one per original frame
    xf_df: object                   # pandas DataFrame from read_xf (all frames)
    tilt_angles_all: List[float]    # all tilt angles (including excluded)
    tilt_angles_good: List[float]   # tilt angles for good frames
    newst_metadata: Dict            # parsed newst.com contents
    tilt_metadata: Dict             # parsed tilt.com contents
    status: str = RelionGeometryStatus.ALGORITHM_PORTED
    warnings: List[str] = field(default_factory=list)


# --------------------------------------------------------------------------- #
# Internal helpers (ports of ImodImport private methods)
# --------------------------------------------------------------------------- #

def _identity4() -> np.ndarray:
    return np.eye(4, dtype=np.float64)


def _load_tilt_projections(
    tlt_path: str,
    center_x: float,
    center_y: float,
    flip_angles: bool,
) -> List[np.ndarray]:
    """Port of RELION ImodImport::loadTiltProjections."""
    angles = read_tlt_angles(tlt_path)
    # w2i0: translates world by (-centerX, 0, -centerY) before rotating
    w2i0 = _identity4()
    w2i0[0, 3] = -center_x
    w2i0[2, 3] = -center_y

    projections: List[np.ndarray] = []
    for a_deg in angles:
        a = math.radians(float(a_deg))
        if flip_angles:
            a *= -1.0

        ca, sa = math.cos(a), math.sin(a)

        w2i = _identity4()
        # Row 0: projects to image-X
        w2i[0, 0] = ca
        w2i[0, 1] = -sa
        w2i[0, 2] = 0.0
        w2i[0, 3] = center_x
        # Row 1: projects to image-Y  (maps world_z → image_y)
        w2i[1, 0] = 0.0
        w2i[1, 1] = 0.0
        w2i[1, 2] = 1.0
        w2i[1, 3] = center_y
        # Row 2: "depth" in the image plane (not used for 2D projection)
        w2i[2, 0] = -sa
        w2i[2, 1] = -ca
        w2i[2, 2] = 0.0
        w2i[2, 3] = 0.0

        projections.append(w2i @ w2i0)
    return projections


def _load_inv_affine_transforms(
    xf_path: str,
    center_orig: Tuple[float, float],
    center_ali: Tuple[float, float],
    binning: float,
    for_ali: bool,
) -> Tuple[List[np.ndarray], object]:
    """Port of RELION ImodImport::loadInvAffineTransforms.

    Returns (xforms, xf_df) where xforms[i] is the 4×4 B matrix and
    xf_df is the raw DataFrame from read_xf (for metadata storage).
    """
    xf_df = read_xf(xf_path)

    P = _identity4()
    P[0, 3] = center_orig[0]
    P[1, 3] = center_orig[1]

    Q = _identity4()
    Q[0, 3] = -center_ali[0]
    Q[1, 3] = -center_ali[1]

    xforms: List[np.ndarray] = []
    for _, row in xf_df.iterrows():
        A = _identity4()
        A[0, 0] = float(row["a11"])
        A[0, 1] = float(row["a12"])
        A[1, 0] = float(row["a21"])
        A[1, 1] = float(row["a22"])
        A[0, 3] = float(row["dx"])
        A[1, 3] = float(row["dy"])

        Ai = np.linalg.inv(A)

        if for_ali:
            B = P @ Q
        else:
            B = P @ Ai @ Q

        # Scale upper 3×3 by binning (matches RELION exactly)
        B_scaled = B.copy()
        for r in range(3):
            for c in range(3):
                B_scaled[r, c] *= binning

        xforms.append(B_scaled)

    return xforms, xf_df


# --------------------------------------------------------------------------- #
# Public function
# --------------------------------------------------------------------------- #

def build_relion_imod_mapping(
    ts_path: str,
    imod_dir: str,
    newst_com: str = "newst.com",
    tilt_com: str = "tilt.com",
    ali: bool = True,
    ali_size: bool = True,
    flip_yz: bool = False,
    flip_z: bool = False,
    flip_angles: bool = False,
    offset_x: float = 0.0,
    offset_y: float = 0.0,
    offset_z: float = 0.0,
    thickness_override: Optional[float] = None,
) -> ImodMapping:
    """
    Build RELION-compatible worldToImage matrices from an IMOD directory.

    Faithful port of RELION ``ImodImport::import()``.  The returned
    ImodMapping.projections list contains one (4×4) matrix per *non-excluded*
    tilt, in acquisition order.  Status is always at least
    ``relion_imod_algorithm_ported``; external RELION oracle validation
    is required to upgrade to ``relion_extract_ready``.
    """
    warnings: List[str] = []

    # --- Step 1: read raw tilt-series stack size ----------------------------
    ts_hdr = read_mrc_header(ts_path)
    w_ts = ts_hdr.nx
    h_ts = ts_hdr.ny
    fc_ts = ts_hdr.nz

    # --- Step 2: parse newst.com --------------------------------------------
    newst_path = os.path.join(imod_dir, newst_com)
    newst_meta: Dict = {}
    orig_stack_fn = ""
    ali_fn = ""
    xf_fn = ""
    w_ali_bin = -1
    h_ali_bin = -1
    binning_ali = 1.0

    if os.path.isfile(newst_path):
        newst_meta = parse_com_file(newst_path)
        orig_stack_fn = com_get(newst_meta, "InputFile", "")
        ali_fn = com_get(newst_meta, "OutputFile", "")
        xf_fn = com_get(newst_meta, "TransformFile", "")
        size_str = com_get(newst_meta, "SizeToOutputInXandY")
        if size_str:
            nums = com_get_numbers(size_str)
            if len(nums) >= 2:
                w_ali_bin, h_ali_bin = int(nums[0]), int(nums[1])
        binning_ali_v = com_get_float(newst_meta, "BinByFactor")
        if binning_ali_v is not None:
            binning_ali = binning_ali_v
    else:
        warnings.append(f"newst.com not found at {newst_path}")

    # --- Step 3: parse tilt.com ---------------------------------------------
    tilt_path = os.path.join(imod_dir, tilt_com)
    tilt_meta: Dict = {}
    ali_fn_tilt = ""
    tlt_fn = ""
    w0: float = -1.0
    h0: float = -1.0
    file_thickness: float = -1.0
    file_binning: float = -1.0
    exclude_set: set = set()
    imod_shift_x: float = 0.0
    imod_shift_y: float = 0.0
    x_axis_tilt: Optional[float] = None

    if os.path.isfile(tilt_path):
        tilt_meta = parse_com_file(tilt_path)
        ali_fn_tilt = com_get(tilt_meta, "InputProjections", "")
        size_str2 = com_get(tilt_meta, "FULLIMAGE")
        if size_str2:
            nums2 = com_get_numbers(size_str2)
            if len(nums2) >= 2:
                w0, h0 = nums2[0], nums2[1]
        tlt_fn = com_get(tilt_meta, "TILTFILE", "")
        thickness_v = com_get_float(tilt_meta, "THICKNESS")
        if thickness_v is not None:
            file_thickness = thickness_v
        binning_v = com_get_float(tilt_meta, "IMAGEBINNED")
        if binning_v is not None:
            file_binning = binning_v
        shift_str = com_get(tilt_meta, "SHIFT")
        if shift_str:
            snums = com_get_numbers(shift_str)
            if len(snums) >= 2:
                imod_shift_x, imod_shift_y = snums[0], snums[1]
        x_axis_tilt_v = com_get_float(tilt_meta, "XAXISTILT")
        if x_axis_tilt_v is not None and abs(x_axis_tilt_v) > 1e-6:
            x_axis_tilt = x_axis_tilt_v
            warnings.append(
                f"XAXISTILT={x_axis_tilt:.3f} detected in tilt.com. "
                "XAXISTILT correction is stored in metadata but is NOT applied "
                "to projection matrices in this port. RELION compatibility requires "
                "external validation."
            )
        # EXCLUDELIST / EXCLUDELIST2 / EXCLUDE
        for excl_key in ("EXCLUDELIST", "EXCLUDELIST2", "EXCLUDE"):
            excl_val = com_get(tilt_meta, excl_key)
            if excl_val:
                try:
                    exclude_set |= parse_exclude_list(excl_val, fc_ts)
                except ValueError as e:
                    warnings.append(f"EXCLUDELIST parse warning: {e}")
    else:
        warnings.append(f"tilt.com not found at {tilt_path}")

    # Warn if newst/tilt reference different aligned stacks
    if ali_fn and ali_fn_tilt and (ali_fn != ali_fn_tilt):
        warnings.append(
            f"Different aligned stacks in newst.com ({ali_fn}) "
            f"and tilt.com ({ali_fn_tilt})"
        )
    # Warn if binning levels differ
    if binning_ali >= 0 and file_binning >= 0 and binning_ali != file_binning:
        warnings.append(
            f"Different binning in newst.com ({binning_ali}) and tilt.com ({file_binning})"
        )

    # --- Step 4: handle excluded frames ------------------------------------
    excluded = [i in exclude_set for i in range(fc_ts)]
    old_frame_index = [i for i in range(fc_ts) if not excluded[i]]
    fc_good = len(old_frame_index)

    if w0 < 0.0:
        warnings.append(f"Original image size (FULLIMAGE) not found in {tilt_com}; using ts stack size")
        w0 = float(w_ts)
        h0 = float(h_ts)

    # --- Step 5: determine aligned image dimensions -------------------------
    if ali_size and w_ali_bin > 0 and h_ali_bin > 0 and binning_ali > 0:
        w_ali = int(round(w_ali_bin * binning_ali))
        h_ali = int(round(h_ali_bin * binning_ali))
    else:
        w_ali = int(round(w0))
        h_ali = int(round(h0))

    # --- Step 6: determine thickness ----------------------------------------
    thickness: float
    if thickness_override is not None and thickness_override > 0:
        thickness = float(thickness_override)
    elif file_thickness > 0:
        thickness = file_thickness
    else:
        warnings.append("Reconstruction thickness not found; using ali-stack height as fallback")
        thickness = float(h_ts)

    # --- Step 7: compute image centres (RELION formula) ---------------------
    orig_center = ((w_ts - 1.0) / 2.0, (h_ts - 1.0) / 2.0)
    ali_center  = ((w_ali - 1.0) / 2.0, (h_ali - 1.0) / 2.0)

    # --- Step 8: find .tlt file ---------------------------------------------
    tlt_path_full: Optional[str] = None
    if tlt_fn:
        cand = tlt_fn if os.path.isabs(tlt_fn) else os.path.join(imod_dir, tlt_fn)
        if os.path.isfile(cand):
            tlt_path_full = cand
    if tlt_path_full is None:
        warnings.append(f"Tilt file '{tlt_fn}' from tilt.com not found; searching imod_dir")
        for name in os.listdir(imod_dir):
            if name.endswith(".tlt") or name.endswith(".rawtlt"):
                tlt_path_full = os.path.join(imod_dir, name)
                break
    if tlt_path_full is None:
        raise FileNotFoundError(f"Cannot locate a .tlt file in {imod_dir}")

    tilt_projs = _load_tilt_projections(
        tlt_path_full, ali_center[0], ali_center[1], flip_angles
    )
    tilt_angles_all = read_tlt_angles(tlt_path_full)

    # --- Step 9: find .xf file and load inverse affine transforms -----------
    xf_path_full: Optional[str] = None
    for cand_fn in [xf_fn] + [f for f in os.listdir(imod_dir) if f.endswith(".xf")]:
        if not cand_fn:
            continue
        cand = cand_fn if os.path.isabs(cand_fn) else os.path.join(imod_dir, cand_fn)
        if os.path.isfile(cand):
            xf_path_full = cand
            break

    if xf_path_full is None:
        # No .xf: use identity transforms
        warnings.append("No .xf file found; using identity 2D transforms")
        n = len(tilt_angles_all)
        aff_xforms = [_identity4() for _ in range(n)]
        import pandas as pd
        xf_df = pd.DataFrame(
            {"tilt_index": range(1, n + 1), "a11": 1.0, "a12": 0.0, "a21": 0.0, "a22": 1.0, "dx": 0.0, "dy": 0.0}
        )
    else:
        center_for_xf = ali_center if ali else orig_center
        aff_xforms, xf_df = _load_inv_affine_transforms(
            xf_path_full, center_for_xf, ali_center, 1.0, ali
        )

    if len(aff_xforms) != fc_ts:
        raise ValueError(
            f"Mismatched frame counts between {tlt_path_full} and .xf: "
            f"tlt={len(tilt_angles_all)}, xf={len(aff_xforms)}, stack={fc_ts}"
        )

    # --- Step 10: toImodOrigin3D --------------------------------------------
    to_imod_origin = _identity4()
    to_imod_origin[0, 3] = -1.0
    to_imod_origin[1, 3] = -thickness / 2.0
    to_imod_origin[2, 3] = -1.0

    # --- Step 11: YzFlip and output dimensions ------------------------------
    w_out = w_ali
    if not flip_yz:
        yz_flip = _identity4()
        h_out = int(round(thickness))
        d_out = h_ali
    elif flip_yz and flip_z:
        yz_flip = np.array(
            [[1, 0, 0, 0],
             [0, 0, -1, thickness - 1],
             [0, 1, 0, 0],
             [0, 0, 0, 1]], dtype=np.float64
        )
        h_out = h_ali
        d_out = int(round(thickness))
    else:  # flip_yz and not flip_z
        yz_flip = np.array(
            [[1, 0, 0, 0],
             [0, 0, 1, 0],
             [0, 1, 0, 0],
             [0, 0, 0, 1]], dtype=np.float64
        )
        h_out = h_ali
        d_out = int(round(thickness))

    # --- Step 12: offset matrix ---------------------------------------------
    Off = _identity4()
    Off[0, 3] += float(offset_x) - imod_shift_x
    Off[1, 3] += float(offset_y) - imod_shift_y
    Off[2, 3] += float(offset_z)

    # --- Step 13: compose worldToImage matrices ----------------------------
    world_to_image: List[np.ndarray] = []
    tilt_angles_good: List[float] = []

    use_old_index = (len(tilt_projs) != fc_good and len(tilt_projs) == fc_ts)

    for f in range(fc_good):
        src = old_frame_index[f] if use_old_index else f
        M = (
            aff_xforms[src]
            @ tilt_projs[src]
            @ to_imod_origin
            @ Off
            @ yz_flip
        )
        world_to_image.append(M)
        tilt_angles_good.append(tilt_angles_all[old_frame_index[f]])

    return ImodMapping(
        projections=world_to_image,
        w=w_out,
        h=h_out,
        d=d_out,
        frames_missing=bool(exclude_set),
        old_frame_index=old_frame_index,
        excluded=excluded,
        xf_df=xf_df,
        tilt_angles_all=tilt_angles_all,
        tilt_angles_good=tilt_angles_good,
        newst_metadata=newst_meta,
        tilt_metadata=tilt_meta,
        status=RelionGeometryStatus.ALGORITHM_PORTED,
        warnings=warnings,
    )


# --------------------------------------------------------------------------- #
# Projection helpers
# --------------------------------------------------------------------------- #

def project_particle_to_ali(
    x_rec: float,
    y_rec: float,
    z_rec: float,
    binning_rec_ali: float,
    wti_matrix: np.ndarray,
) -> Tuple[float, float]:
    """Project a particle from reconstruction-voxel space to aligned-stack image coords.

    Convention (no flip_yz):
        world_x = x_rec * B_rec_ali   (ali-image X)
        world_y = z_rec * B_rec_ali   (reconstruction depth / IMOD thickness direction)
        world_z = y_rec * B_rec_ali   (ali-image Y)

    Returns (x_img, y_img) in bin-1 ali-stack pixel coordinates.
    """
    B = float(binning_rec_ali)
    world = np.array([
        float(x_rec) * B,
        float(z_rec) * B,   # rec MRC-Z (depth) → IMOD world Y
        float(y_rec) * B,   # rec MRC-Y → IMOD world Z
        1.0,
    ], dtype=np.float64)
    img = wti_matrix @ world
    return float(img[0]), float(img[1])


def invert_xf_to_raw(
    x_ali: float,
    y_ali: float,
    xf_df,
    frame_index0: int,
) -> Tuple[float, float]:
    """Map aligned-stack coordinates back to raw-stack coordinates using .xf inverse.

    ``frame_index0`` is the 0-based original frame index (before exclusions).
    Returns (x_raw, y_raw) in raw-stack pixel coordinates.
    """
    from tomojanas.metadata.imod_xf import invert_xf_for_point
    row = xf_df.iloc[frame_index0]
    return invert_xf_for_point(x_ali, y_ali, row)
