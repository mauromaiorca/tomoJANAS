#!/usr/bin/env python3
"""
MRC I/O utilities (little-endian, MRC2014-compatible).

BACKWARD-COMPATIBILITY SHIM
---------------------------
The authoritative MRC implementation now lives in :mod:`tomojanas.io.mrc`.
This module re-exports that functionality under the historical ``janas``
API so existing reconstruction/locres code keeps working unchanged.

Key compatibility guarantees preserved here:
* ``read_mrc_header(path)`` returns a plain ``dict`` (not a dataclass), with
  the same keys callers already index (``hdr["nx"]``, ``hdr.get("nsymbt",0)`` …).
* ``read_mrc_data(path)`` returns ``(volume, header_dict)``.
* ``memmap_mrc_stack(path)`` returns ``(memmap, header_dict)``.

New code should import from ``tomojanas.io.mrc`` directly (it returns the
structured :class:`tomojanas.io.mrc.MRCHeader`).

Notes
-----
* Voxel order is (Z, Y, X).
* ORIGIN is stored in Ångströms (origin_A), not angles.
* We write mode=2 (float32) on output.
* Data starts at ``1024 + nsymbt`` (extended headers handled).
"""

from __future__ import annotations

from typing import Dict, Tuple

import numpy as np

from tomojanas.io import mrc as _mrc

__all__ = [
    "write_mrc",
    "read_mrc_header",
    "read_mrc_data",
    "write_mrc_like",
    "write_mrc_header",
    "precreate_mrcs",
    "memmap_mrc_stack",
    "mrc_dtype_from_mode",
]


def mrc_dtype_from_mode(mode: int) -> np.dtype:
    """Map MRC mode to numpy dtype (see tomojanas.io.mrc)."""
    return _mrc.mrc_dtype_from_mode(mode)


def read_mrc_header(path: str) -> Dict[str, object]:
    """Parse the 1024-byte MRC header and return a legacy dict."""
    return _mrc.read_mrc_header(path).as_dict()


def read_mrc_data(path: str, astype=np.float32) -> Tuple[np.ndarray, Dict[str, object]]:
    """Read a full MRC volume; returns (volume, header_dict)."""
    vol, hdr = _mrc.read_mrc_data(path, astype=astype)
    return vol, hdr.as_dict()


def write_mrc(
    path: str,
    volume: np.ndarray,
    pixel_size: float,
    origin_angs: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    update_stats: bool = True,
) -> None:
    _mrc.write_mrc(path, volume, pixel_size, origin_angs=origin_angs,
                   update_stats=update_stats)


def write_mrc_like(
    path: str,
    volume: np.ndarray,
    template_hdr: Dict[str, object],
    update_stats: bool = True,
) -> None:
    """Write an MRC using a template header (legacy dict or MRCHeader)."""
    _mrc.write_mrc_like(path, volume, template_hdr, update_stats=update_stats)


def memmap_mrc_stack(path: str) -> Tuple[np.memmap, Dict[str, object]]:
    """Read-only memmap of an MRC stack as (nz, ny, nx); returns (mm, dict)."""
    mm, hdr = _mrc.open_mrc_memmap(path, mode="r")
    return mm, hdr.as_dict()


def write_mrc_header(
    fh,
    *,
    nx: int,
    ny: int,
    nz: int,
    apix: float,
    vmin: float = 0.0,
    vmax: float = 0.0,
    vmean: float = 0.0,
    vrms: float = 0.0,
    mode: int = 2,
) -> None:
    """Write a minimal 1024-byte MRC header at the start of ``fh``."""
    header = _mrc._pack_header(
        nx=nx, ny=ny, nz=nz, mode=mode,
        pixel_size_xyz=(apix, apix, apix),
        vmin=vmin, vmax=vmax, vmean=vmean, vrms=vrms,
    )
    fh.seek(0)
    fh.write(header)


def precreate_mrcs(path_out: str, *, nz: int, ny: int, nx: int, apix: float) -> None:
    """Pre-allocate a float32 .mrcs file on disk (header + zero body)."""
    _mrc.precreate_mrc_stack(path_out, nz=nz, ny=ny, nx=nx, pixel_size=apix, mode=2)
