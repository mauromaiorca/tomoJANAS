#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/geometry/roi.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
ROI (Region of Interest) geometry helpers.

The scientific ROI is:
  3D: sphere     storage: cube
  2D: circle     storage: square

Cubes/squares are MRC storage containers. The biological mask is sphere/circle.
"""
from __future__ import annotations
import math
from typing import Optional

__all__ = [
    "roi_radius_angst_from_args",
    "projection_radius_px",
    "storage_box_from_radius",
    "circle_inside_frame",
    "sphere_inside_volume",
]


def roi_radius_angst_from_args(args, pixel_size: Optional[float] = None) -> float:
    """Resolve ROI radius in Å from argparse namespace.

    Priority: --roi-radius-angst > --roi-diameter-angst/2 > --roi-radius-voxel*pixel_size.
    Raises ValueError if none are provided.
    """
    r = getattr(args, "roi_radius_angst", None)
    if r is not None and r > 0:
        return float(r)
    d = getattr(args, "roi_diameter_angst", None)
    if d is not None and d > 0:
        return float(d) / 2.0
    rv = getattr(args, "roi_radius_voxel", None)
    if rv is not None and rv > 0:
        if pixel_size is None:
            raise ValueError("--roi-radius-voxel requires --apix or a pixel size from the MRC header")
        return float(rv) * float(pixel_size)
    raise ValueError(
        "Specify at least one of --roi-radius-angst, --roi-diameter-angst, --roi-radius-voxel"
    )


def projection_radius_px(radius_angst: float, pixel_size: float) -> float:
    """Circular ROI radius in pixels for a given sphere radius in Å."""
    return float(radius_angst) / float(pixel_size)


def storage_box_from_radius(radius_px_or_vox: float, padding: float = 0.0) -> int:
    """Smallest even integer side length enclosing a sphere/circle of given radius.

    The cube/square is the storage container; the sphere/circle is the actual ROI.
    """
    side = 2.0 * (float(radius_px_or_vox) + float(padding))
    box = int(math.ceil(side))
    if box % 2 == 1:
        box += 1  # keep even for FFT-friendly sizes
    return max(box, 2)


def circle_inside_frame(
    cx: float, cy: float, radius: float, nx: int, ny: int
) -> bool:
    """Return True if the entire circle fits within the frame (no clipping)."""
    return (
        cx - radius >= 0
        and cy - radius >= 0
        and cx + radius <= float(nx)
        and cy + radius <= float(ny)
    )


def sphere_inside_volume(
    cx: float, cy: float, cz: float,
    radius: float,
    nx: int, ny: int, nz: int,
) -> bool:
    """Return True if the entire sphere fits within the volume."""
    return (
        cx - radius >= 0 and cx + radius <= float(nx)
        and cy - radius >= 0 and cy + radius <= float(ny)
        and cz - radius >= 0 and cz + radius <= float(nz)
    )
