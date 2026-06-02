#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/geometry/coordinates.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
Coordinate conversions between reconstruction voxels and RELION centered Angstrom.

Formula (RELION convention):
    x_ali_px = x_rec * rlnTomoTomogramBinning
    rlnCenteredCoordinateXAngst = (x_ali_px - rlnTomoSizeX / 2) * rlnTomoTiltSeriesPixelSize

Indexing note:
    "zero-based" : rec voxel coordinates start at 0 (Python/NumPy convention).
    "one-based"  : IMOD model convention — subtract 1 before applying formula.
"""
from __future__ import annotations
from typing import Tuple, Sequence

__all__ = [
    "rec_voxel_to_relion_centered_angst",
    "relion_centered_angst_to_rec_voxel",
    "normalise_coordinate_order",
    "coordinate_roundtrip_error",
]


def rec_voxel_to_relion_centered_angst(
    x_rec: float,
    y_rec: float,
    z_rec: float,
    tomo_size_x: float,
    tomo_size_y: float,
    tomo_size_z: float,
    tomo_tilt_series_pixel_size: float,
    tomo_tomogram_binning: float,
    indexing: str = "zero-based",
) -> Tuple[float, float, float]:
    """Convert 0- or 1-based reconstruction voxel coords to RELION centered Å."""
    if indexing == "one-based":
        x_rec = float(x_rec) - 1.0
        y_rec = float(y_rec) - 1.0
        z_rec = float(z_rec) - 1.0
    else:
        x_rec, y_rec, z_rec = float(x_rec), float(y_rec), float(z_rec)

    B = float(tomo_tomogram_binning)
    A = float(tomo_tilt_series_pixel_size)

    x_ang = (x_rec * B - float(tomo_size_x) / 2.0) * A
    y_ang = (y_rec * B - float(tomo_size_y) / 2.0) * A
    z_ang = (z_rec * B - float(tomo_size_z) / 2.0) * A
    return x_ang, y_ang, z_ang


def relion_centered_angst_to_rec_voxel(
    x_ang: float,
    y_ang: float,
    z_ang: float,
    tomo_size_x: float,
    tomo_size_y: float,
    tomo_size_z: float,
    tomo_tilt_series_pixel_size: float,
    tomo_tomogram_binning: float,
    indexing: str = "zero-based",
) -> Tuple[float, float, float]:
    """Inverse: RELION centered Å → reconstruction voxel coords."""
    A = float(tomo_tilt_series_pixel_size)
    B = float(tomo_tomogram_binning)

    x_ali = float(x_ang) / A + float(tomo_size_x) / 2.0
    y_ali = float(y_ang) / A + float(tomo_size_y) / 2.0
    z_ali = float(z_ang) / A + float(tomo_size_z) / 2.0

    x_rec = x_ali / B
    y_rec = y_ali / B
    z_rec = z_ali / B

    if indexing == "one-based":
        x_rec += 1.0
        y_rec += 1.0
        z_rec += 1.0

    return x_rec, y_rec, z_rec


def normalise_coordinate_order(
    values: Sequence[float], axis_order: str
) -> Tuple[float, float, float]:
    """Reorder coordinates from ``axis_order`` convention to ``(x, y, z)``.

    ``axis_order`` is any permutation of the letters x, y, z describing the
    order in which the three input values are given. All 6 permutations are
    supported (``xyz``, ``xzy``, ``yxz``, ``yzx``, ``zxy``, ``zyx``).

    Example: ``axis_order="xzy"`` means the inputs are (X, Z, Y) — the common
    case for IMOD 3dmod readouts of a "flipped" tomogram, where Y and Z are
    swapped relative to the file's axes.
    """
    v = list(values)
    if len(v) < 3:
        raise ValueError(f"need at least 3 values, got {len(v)}")
    ao = str(axis_order).lower()
    if sorted(ao) != ["x", "y", "z"]:
        raise ValueError(
            f"unsupported axis_order: {axis_order!r}; must be a permutation of "
            f"'x', 'y', 'z' (e.g. xyz, xzy, zyx)"
        )
    m = {ao[i]: float(v[i]) for i in range(3)}
    return m["x"], m["y"], m["z"]


def coordinate_roundtrip_error(
    x_rec: float, y_rec: float, z_rec: float,
    tomo_size_x: float, tomo_size_y: float, tomo_size_z: float,
    tomo_tilt_series_pixel_size: float,
    tomo_tomogram_binning: float,
    indexing: str = "zero-based",
) -> float:
    """Max absolute voxel error for a rec→angst→rec round-trip."""
    ang = rec_voxel_to_relion_centered_angst(
        x_rec, y_rec, z_rec,
        tomo_size_x, tomo_size_y, tomo_size_z,
        tomo_tilt_series_pixel_size, tomo_tomogram_binning, indexing,
    )
    back = relion_centered_angst_to_rec_voxel(
        *ang,
        tomo_size_x, tomo_size_y, tomo_size_z,
        tomo_tilt_series_pixel_size, tomo_tomogram_binning, indexing,
    )
    return max(abs(back[i] - [x_rec, y_rec, z_rec][i]) for i in range(3))
