#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/io/mrc.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
Consolidated MRC / MRCS / ST I/O for tomoJANAS (little-endian, MRC2014).

This is the authoritative MRC module for the whole project. The legacy
``janas.IO_utils`` re-exports from here for backward compatibility.

Conventions
-----------
* NumPy array order is ``(Z, Y, X)``. For an aligned tilt-series ``.ali``
  stack: ``Z`` = tilt index, ``Y`` = image row, ``X`` = image column.
* ORIGIN is stored in Ångström (origin_x/y/z).
* Data starts at byte ``1024 + nsymbt`` (NOT a hard-coded 1024) so files
  with extended headers are handled correctly.
* Output volumes are written as mode 2 (float32) by default.

Supported read modes: 0 (int8), 1 (int16), 2 (float32), 6 (uint16),
12 (float16). Complex modes 3 and 4 are rejected with a clear error.
"""

from __future__ import annotations

import math
import struct
from dataclasses import dataclass, asdict
from typing import Dict, Optional, Sequence, Tuple, Union

import numpy as np

__all__ = [
    "MRCHeader",
    "mrc_dtype_from_mode",
    "read_mrc_header",
    "read_mrc_data",
    "write_mrc",
    "write_mrc_like",
    "open_mrc_memmap",
    "read_mrc_slice",
    "precreate_mrc_stack",
    "write_mrc_slice",
    "write_cropped_mrc_like",
    "box_corner",
    "crop_volume_box",
    "crop_volume_sphere",
    "crop_image_square",
    "crop_image_circle",
    "make_spherical_mask",
    "make_circular_mask",
    "validate_mrc_geometry",
    "validate_mrc_pixel_size",
]

# Word offsets are in 4-byte words from the start of the 1024-byte header.
# This layout matches what the legacy janas.IO_utils reads/writes so files
# produced by either path are mutually compatible.
_MAP_STAMP = b"MAP "
_MACHINE_STAMP = b"DA\x00\x00"  # little-endian


# --------------------------------------------------------------------------- #
# Header dataclass
# --------------------------------------------------------------------------- #
@dataclass
class MRCHeader:
    nx: int
    ny: int
    nz: int
    mode: int
    nxstart: int
    nystart: int
    nzstart: int
    mx: int
    my: int
    mz: int
    xlen: float
    ylen: float
    zlen: float
    alpha: float
    beta: float
    gamma: float
    mapc: int
    mapr: int
    maps: int
    amin: float
    amax: float
    amean: float
    ispg: int
    nsymbt: int
    origin_x: float
    origin_y: float
    origin_z: float
    pixel_x: float
    pixel_y: float
    pixel_z: float
    stamp: str = "MAP "

    # --- derived helpers --------------------------------------------------- #
    @property
    def data_offset(self) -> int:
        """Byte offset to the first voxel: 1024 + nsymbt (extended header)."""
        return 1024 + int(self.nsymbt)

    @property
    def shape(self) -> Tuple[int, int, int]:
        """NumPy (Z, Y, X) shape."""
        return (int(self.nz), int(self.ny), int(self.nx))

    @property
    def dtype(self) -> np.dtype:
        return mrc_dtype_from_mode(self.mode)

    @property
    def pixel_size(self) -> float:
        """Representative (X) pixel size in Ångström."""
        return float(self.pixel_x)

    def as_dict(self) -> Dict[str, object]:
        """Return the legacy dict representation (same keys as the old
        ``janas.IO_utils.read_mrc_header``)."""
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, object]) -> "MRCHeader":
        fields = cls.__dataclass_fields__  # type: ignore[attr-defined]
        kwargs = {k: d[k] for k in fields if k in d}
        return cls(**kwargs)  # type: ignore[arg-type]


# --------------------------------------------------------------------------- #
# Mode handling
# --------------------------------------------------------------------------- #
def mrc_dtype_from_mode(mode: int) -> np.dtype:
    """Map an MRC mode to a numpy dtype.

    Supported: 0->int8, 1->int16, 2->float32, 6->uint16, 12->float16.
    Complex modes 3 and 4 are not supported for import/crop and raise.
    """
    mapping = {
        0: np.int8,
        1: np.int16,
        2: np.float32,
        6: np.uint16,
        12: np.float16,
    }
    if mode in mapping:
        return np.dtype(mapping[mode])
    if mode in (3, 4):
        raise NotImplementedError(
            f"MRC mode {mode} is a complex/transform format and is not "
            "supported for tomoJANAS import/crop operations."
        )
    raise NotImplementedError(f"Unsupported MRC mode: {mode}")


# --------------------------------------------------------------------------- #
# Header read / pack
# --------------------------------------------------------------------------- #
def read_mrc_header(path: str) -> MRCHeader:
    with open(path, "rb") as f:
        h = f.read(1024)
    if len(h) < 1024:
        raise ValueError(f"File too short to contain an MRC header: {path}")

    nx, ny, nz, mode = struct.unpack("<4i", h[0:16])
    nxstart, nystart, nzstart = struct.unpack("<3i", h[16:28])
    mx, my, mz = struct.unpack("<3i", h[28:40])
    xlen, ylen, zlen = struct.unpack("<3f", h[40:52])
    alpha, beta, gamma = struct.unpack("<3f", h[52:64])
    mapc, mapr, maps = struct.unpack("<3i", h[64:76])
    amin, amax, amean = struct.unpack("<3f", h[76:88])
    ispg, nsymbt = struct.unpack("<2i", h[88:96])
    try:
        ox, oy, oz = struct.unpack("<3f", h[196:208])  # ORIGIN (Å)
    except struct.error:  # pragma: no cover - defensive
        ox, oy, oz = 0.0, 0.0, 0.0
    stamp = h[208:212]

    if mx <= 0:
        mx = nx
    if my <= 0:
        my = ny
    if mz <= 0:
        mz = nz

    pix_x = xlen / mx if mx else 1.0
    pix_y = ylen / my if my else pix_x
    pix_z = zlen / mz if mz else pix_x

    return MRCHeader(
        nx=nx, ny=ny, nz=nz, mode=mode,
        nxstart=nxstart, nystart=nystart, nzstart=nzstart,
        mx=mx, my=my, mz=mz,
        xlen=float(xlen), ylen=float(ylen), zlen=float(zlen),
        alpha=float(alpha), beta=float(beta), gamma=float(gamma),
        mapc=mapc, mapr=mapr, maps=maps,
        amin=float(amin), amax=float(amax), amean=float(amean),
        ispg=ispg, nsymbt=nsymbt,
        origin_x=float(ox), origin_y=float(oy), origin_z=float(oz),
        pixel_x=float(pix_x), pixel_y=float(pix_y), pixel_z=float(pix_z),
        stamp=stamp.decode(errors="ignore"),
    )


def _pack_header(
    *,
    nx: int, ny: int, nz: int, mode: int,
    pixel_size_xyz: Tuple[float, float, float],
    origin_angs: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    nxstart: int = 0, nystart: int = 0, nzstart: int = 0,
    mapc: int = 1, mapr: int = 2, maps: int = 3,
    ispg: int = 0, nsymbt: int = 0,
    vmin: float = 0.0, vmax: float = 0.0, vmean: float = 0.0, vrms: float = 0.0,
    angles: Tuple[float, float, float] = (90.0, 90.0, 90.0),
) -> bytearray:
    header = bytearray(1024)

    def pw(word_off: int, fmt: str, *vals):
        struct.pack_into("<" + fmt, header, word_off * 4, *vals)

    px, py, pz = pixel_size_xyz
    pw(0, "i", nx); pw(1, "i", ny); pw(2, "i", nz); pw(3, "i", mode)
    pw(4, "i", int(nxstart)); pw(5, "i", int(nystart)); pw(6, "i", int(nzstart))
    pw(7, "i", nx); pw(8, "i", ny); pw(9, "i", nz)
    pw(10, "f", nx * px); pw(11, "f", ny * py); pw(12, "f", nz * pz)
    pw(13, "f", angles[0]); pw(14, "f", angles[1]); pw(15, "f", angles[2])
    pw(16, "i", int(mapc)); pw(17, "i", int(mapr)); pw(18, "i", int(maps))
    pw(19, "f", float(vmin)); pw(20, "f", float(vmax)); pw(21, "f", float(vmean))
    pw(22, "i", int(ispg)); pw(23, "i", int(nsymbt))
    ox, oy, oz = origin_angs
    pw(49, "f", float(ox)); pw(50, "f", float(oy)); pw(51, "f", float(oz))
    header[208:212] = _MAP_STAMP
    header[212:216] = _MACHINE_STAMP
    pw(54, "f", float(vrms))
    pw(55, "i", 0)
    return header


def _stats(vol: np.ndarray, update_stats: bool) -> Tuple[float, float, float, float]:
    if not update_stats or vol.size == 0:
        return 0.0, 0.0, 0.0, 0.0
    return (
        float(np.nanmin(vol)),
        float(np.nanmax(vol)),
        float(np.nanmean(vol)),
        float(np.nanstd(vol)),
    )


# --------------------------------------------------------------------------- #
# Full-volume read / write
# --------------------------------------------------------------------------- #
def read_mrc_data(path: str, astype=np.float32) -> Tuple[np.ndarray, MRCHeader]:
    hdr = read_mrc_header(path)
    dt = mrc_dtype_from_mode(hdr.mode)
    nvox = hdr.nx * hdr.ny * hdr.nz
    with open(path, "rb") as f:
        f.seek(hdr.data_offset)
        arr = np.fromfile(f, dtype=dt, count=nvox)
    if arr.size != nvox:
        raise ValueError(
            f"MRC data size ({arr.size}) inconsistent with header "
            f"({nvox}) for {path}"
        )
    vol = arr.reshape(hdr.shape)
    if dt != np.float32:
        vol = vol.astype(np.float32)
    if astype is not None and np.dtype(astype) != np.float32:
        vol = vol.astype(astype)
    return vol, hdr


def write_mrc(
    path: str,
    volume: np.ndarray,
    pixel_size: float,
    origin_angs: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    update_stats: bool = True,
) -> None:
    vol = np.asarray(volume, dtype=np.float32)
    if vol.ndim != 3:
        raise ValueError(f"write_mrc expects a 3D volume, got ndim={vol.ndim}")
    nz, ny, nx = vol.shape
    vmin, vmax, vmean, vrms = _stats(vol, update_stats)
    header = _pack_header(
        nx=nx, ny=ny, nz=nz, mode=2,
        pixel_size_xyz=(pixel_size, pixel_size, pixel_size),
        origin_angs=origin_angs,
        vmin=vmin, vmax=vmax, vmean=vmean, vrms=vrms,
    )
    with open(path, "wb") as f:
        f.write(header)
        f.write(vol.tobytes())


def write_mrc_like(
    path: str,
    volume: np.ndarray,
    template_header: Union[MRCHeader, Dict[str, object]],
    update_stats: bool = True,
) -> None:
    """Write a float32 MRC using the geometry of an existing header.

    ``template_header`` may be an :class:`MRCHeader` or a legacy dict.
    The volume shape must match the template's (nz, ny, nx).
    """
    hdr = template_header if isinstance(template_header, MRCHeader) \
        else MRCHeader.from_dict(dict(template_header))
    vol = np.asarray(volume, dtype=np.float32)
    nz, ny, nx = vol.shape
    if (nx, ny, nz) != (hdr.nx, hdr.ny, hdr.nz):
        raise ValueError(
            "write_mrc_like: volume shape (z,y,x)="
            f"{vol.shape} does not match template (z,y,x)="
            f"{(hdr.nz, hdr.ny, hdr.nx)}"
        )
    vmin, vmax, vmean, vrms = _stats(vol, update_stats)
    header = _pack_header(
        nx=nx, ny=ny, nz=nz, mode=2,
        pixel_size_xyz=(hdr.pixel_x, hdr.pixel_y, hdr.pixel_z),
        origin_angs=(hdr.origin_x, hdr.origin_y, hdr.origin_z),
        nxstart=hdr.nxstart, nystart=hdr.nystart, nzstart=hdr.nzstart,
        mapc=hdr.mapc, mapr=hdr.mapr, maps=hdr.maps,
        ispg=hdr.ispg, nsymbt=0,  # crop/derived files drop the extended header
        vmin=vmin, vmax=vmax, vmean=vmean, vrms=vrms,
        angles=(hdr.alpha, hdr.beta, hdr.gamma),
    )
    with open(path, "wb") as f:
        f.write(header)
        f.write(vol.tobytes())


# --------------------------------------------------------------------------- #
# Memmap / slice access
# --------------------------------------------------------------------------- #
def open_mrc_memmap(path: str, mode: str = "r") -> Tuple[np.memmap, MRCHeader]:
    """Memory-map an MRC/MRCS/ST file as a ``(nz, ny, nx)`` array.

    Honours ``1024 + nsymbt`` as the data offset. ``mode`` is passed to
    :class:`numpy.memmap` (``"r"`` read-only, ``"r+"`` read/write).
    """
    hdr = read_mrc_header(path)
    dt = mrc_dtype_from_mode(hdr.mode)
    mm = np.memmap(path, mode=mode, dtype=dt, offset=hdr.data_offset, shape=hdr.shape)
    return mm, hdr


def read_mrc_slice(path: str, slice_index0: int) -> Tuple[np.ndarray, MRCHeader]:
    """Read a single Z-slice (image) ``slice_index0`` (0-based) as ``(ny, nx)``."""
    hdr = read_mrc_header(path)
    if slice_index0 < 0 or slice_index0 >= hdr.nz:
        raise IndexError(
            f"slice index {slice_index0} out of range [0, {hdr.nz - 1}] for {path}"
        )
    dt = mrc_dtype_from_mode(hdr.mode)
    itemsize = np.dtype(dt).itemsize
    plane = hdr.ny * hdr.nx
    offset = hdr.data_offset + slice_index0 * plane * itemsize
    with open(path, "rb") as f:
        f.seek(offset)
        arr = np.fromfile(f, dtype=dt, count=plane)
    if arr.size != plane:
        raise ValueError(f"Short read for slice {slice_index0} of {path}")
    img = arr.reshape((hdr.ny, hdr.nx))
    if dt != np.float32:
        img = img.astype(np.float32)
    return img, hdr


def precreate_mrc_stack(
    path: str, nz: int, ny: int, nx: int, pixel_size: float, mode: int = 2
) -> None:
    """Pre-allocate an MRC stack on disk (header + zero-filled body) without
    holding the data in RAM. Useful before streaming per-slice writes."""
    dt = mrc_dtype_from_mode(mode)
    itemsize = np.dtype(dt).itemsize
    header = _pack_header(
        nx=nx, ny=ny, nz=nz, mode=mode,
        pixel_size_xyz=(pixel_size, pixel_size, pixel_size),
    )
    nbytes = nz * ny * nx * itemsize
    with open(path, "wb") as f:
        f.write(header)
        if nbytes > 0:
            f.seek(1024 + nbytes - 1)
            f.write(b"\0")


def write_mrc_slice(path: str, slice_index0: int, image: np.ndarray) -> None:
    """Write a single Z-slice into an existing (pre-created) MRC stack."""
    hdr = read_mrc_header(path)
    if slice_index0 < 0 or slice_index0 >= hdr.nz:
        raise IndexError(
            f"slice index {slice_index0} out of range [0, {hdr.nz - 1}] for {path}"
        )
    img = np.asarray(image)
    if img.shape != (hdr.ny, hdr.nx):
        raise ValueError(
            f"slice shape {img.shape} != expected {(hdr.ny, hdr.nx)} for {path}"
        )
    dt = mrc_dtype_from_mode(hdr.mode)
    itemsize = np.dtype(dt).itemsize
    plane = hdr.ny * hdr.nx
    offset = hdr.data_offset + slice_index0 * plane * itemsize
    with open(path, "r+b") as f:
        f.seek(offset)
        f.write(np.ascontiguousarray(img, dtype=dt).tobytes())


# --------------------------------------------------------------------------- #
# Cropping primitives (ROI is sphere/circle; storage is cube/square)
# --------------------------------------------------------------------------- #
def box_corner(centre_zyx: Sequence[float], box_size: int) -> Tuple[int, int, int]:
    """Integer (z0, y0, x0) source corner for a centred cubic box."""
    half = box_size // 2
    cz, cy, cx = centre_zyx
    return (
        int(round(cz)) - half,
        int(round(cy)) - half,
        int(round(cx)) - half,
    )


def crop_volume_box(
    volume: np.ndarray,
    centre_zyx: Sequence[float],
    box_size: int,
    pad_value: float = 0.0,
) -> np.ndarray:
    """Extract a cubic ``box_size`` sub-volume centred at ``centre_zyx``.

    Out-of-bounds regions are filled with ``pad_value`` (the storage box is
    only a container; biological masking is applied separately).
    """
    if box_size <= 0:
        raise ValueError("box_size must be positive")
    nz, ny, nx = volume.shape
    z0, y0, x0 = box_corner(centre_zyx, box_size)
    z1, y1, x1 = z0 + box_size, y0 + box_size, x0 + box_size

    out = np.full((box_size, box_size, box_size), pad_value, dtype=np.float32)
    sz0, sz1 = max(z0, 0), min(z1, nz)
    sy0, sy1 = max(y0, 0), min(y1, ny)
    sx0, sx1 = max(x0, 0), min(x1, nx)
    if sz0 < sz1 and sy0 < sy1 and sx0 < sx1:
        out[sz0 - z0:sz1 - z0, sy0 - y0:sy1 - y0, sx0 - x0:sx1 - x0] = \
            volume[sz0:sz1, sy0:sy1, sx0:sx1]
    return out


def _enclosing_box(radius: float, padding: float) -> int:
    side = int(math.ceil(2.0 * (float(radius) + float(padding))))
    return max(side, 1)


def crop_volume_sphere(
    volume: np.ndarray,
    centre_zyx: Sequence[float],
    radius_vox: float,
    padding_vox: float = 0,
) -> np.ndarray:
    """Extract the cubic storage container that encloses a spherical ROI of
    ``radius_vox`` (plus optional ``padding_vox``). The returned array is the
    cube; use :func:`make_spherical_mask` to obtain the sphere mask."""
    box = _enclosing_box(radius_vox, padding_vox)
    return crop_volume_box(volume, centre_zyx, box)


def crop_image_square(
    image: np.ndarray,
    centre_yx: Sequence[float],
    box_size: int,
    pad_value: float = 0.0,
) -> np.ndarray:
    """Extract a square ``box_size`` patch centred at ``centre_yx`` (2D)."""
    if box_size <= 0:
        raise ValueError("box_size must be positive")
    ny, nx = image.shape
    half = box_size // 2
    cy, cx = centre_yx
    y0 = int(round(cy)) - half
    x0 = int(round(cx)) - half
    y1, x1 = y0 + box_size, x0 + box_size
    out = np.full((box_size, box_size), pad_value, dtype=np.float32)
    sy0, sy1 = max(y0, 0), min(y1, ny)
    sx0, sx1 = max(x0, 0), min(x1, nx)
    if sy0 < sy1 and sx0 < sx1:
        out[sy0 - y0:sy1 - y0, sx0 - x0:sx1 - x0] = image[sy0:sy1, sx0:sx1]
    return out


def crop_image_circle(
    image: np.ndarray,
    centre_yx: Sequence[float],
    radius_px: float,
    padding_px: float = 0,
) -> np.ndarray:
    """Extract the square storage patch enclosing a circular ROI of
    ``radius_px`` (plus optional ``padding_px``)."""
    box = _enclosing_box(radius_px, padding_px)
    return crop_image_square(image, centre_yx, box)


def make_spherical_mask(
    box_size: int, radius_vox: float, centre_local: Optional[Sequence[float]] = None
) -> np.ndarray:
    """Boolean ``(box, box, box)`` mask: True inside the sphere of radius
    ``radius_vox`` centred at ``centre_local`` (defaults to the box centre)."""
    if centre_local is None:
        c = (box_size - 1) / 2.0
        centre_local = (c, c, c)
    zz, yy, xx = np.ogrid[:box_size, :box_size, :box_size]
    cz, cy, cx = centre_local
    d2 = (zz - cz) ** 2 + (yy - cy) ** 2 + (xx - cx) ** 2
    return d2 <= float(radius_vox) ** 2


def make_circular_mask(
    box_size: int, radius_px: float, centre_local: Optional[Sequence[float]] = None
) -> np.ndarray:
    """Boolean ``(box, box)`` mask: True inside the circle of radius
    ``radius_px`` centred at ``centre_local`` (defaults to the box centre)."""
    if centre_local is None:
        c = (box_size - 1) / 2.0
        centre_local = (c, c)
    yy, xx = np.ogrid[:box_size, :box_size]
    cy, cx = centre_local
    d2 = (yy - cy) ** 2 + (xx - cx) ** 2
    return d2 <= float(radius_px) ** 2


# --------------------------------------------------------------------------- #
# Crop writer (keeps geometry, updates origin)
# --------------------------------------------------------------------------- #
def write_cropped_mrc_like(
    path: str,
    cropped_volume: np.ndarray,
    source_header: Union[MRCHeader, Dict[str, object]],
    crop_origin_zyx: Sequence[float],
    pixel_size: Optional[float] = None,
    update_origin: bool = True,
    update_stats: bool = True,
) -> None:
    """Write a cropped volume as float32, deriving geometry from the source.

    1. writes float32 (mode 2);
    2. updates nx, ny, nz (and mx, my, mz) from the crop;
    3. updates xlen, ylen, zlen accordingly;
    4. preserves pixel_x/y/z unless ``pixel_size`` overrides them;
    5. updates ORIGIN (Å) by the crop corner if ``update_origin``;
    6. writes amin, amax, amean and rms statistics.
    """
    hdr = source_header if isinstance(source_header, MRCHeader) \
        else MRCHeader.from_dict(dict(source_header))
    vol = np.asarray(cropped_volume, dtype=np.float32)
    if vol.ndim != 3:
        raise ValueError(f"write_cropped_mrc_like expects 3D, got ndim={vol.ndim}")
    nz, ny, nx = vol.shape

    if pixel_size is not None:
        px = py = pz = float(pixel_size)
    else:
        px, py, pz = hdr.pixel_x, hdr.pixel_y, hdr.pixel_z

    if update_origin:
        cz0, cy0, cx0 = crop_origin_zyx
        ox = hdr.origin_x + float(cx0) * px
        oy = hdr.origin_y + float(cy0) * py
        oz = hdr.origin_z + float(cz0) * pz
    else:
        ox, oy, oz = hdr.origin_x, hdr.origin_y, hdr.origin_z

    vmin, vmax, vmean, vrms = _stats(vol, update_stats)
    header = _pack_header(
        nx=nx, ny=ny, nz=nz, mode=2,
        pixel_size_xyz=(px, py, pz),
        origin_angs=(ox, oy, oz),
        mapc=hdr.mapc, mapr=hdr.mapr, maps=hdr.maps,
        ispg=hdr.ispg, nsymbt=0,
        vmin=vmin, vmax=vmax, vmean=vmean, vrms=vrms,
        angles=(hdr.alpha, hdr.beta, hdr.gamma),
    )
    with open(path, "wb") as f:
        f.write(header)
        f.write(vol.tobytes())


# --------------------------------------------------------------------------- #
# Geometry / pixel-size validation
# --------------------------------------------------------------------------- #
def validate_mrc_geometry(path_a: str, path_b: str) -> Dict[str, object]:
    """Compare the (nx, ny, nz) dimensions of two MRC files.

    Returns a dict with ``ok`` plus the two shapes.
    """
    a = read_mrc_header(path_a)
    b = read_mrc_header(path_b)
    ok = (a.nx, a.ny, a.nz) == (b.nx, b.ny, b.nz)
    return {
        "ok": ok,
        "shape_a": (a.nx, a.ny, a.nz),
        "shape_b": (b.nx, b.ny, b.nz),
        "path_a": path_a,
        "path_b": path_b,
    }


def validate_mrc_pixel_size(
    path: str, expected_apix: float, tolerance: float = 1e-3
) -> Dict[str, object]:
    """Check that a file's pixel size matches ``expected_apix`` within an
    absolute ``tolerance`` (Å). Returns a dict with ``ok`` and the delta."""
    hdr = read_mrc_header(path)
    delta = abs(hdr.pixel_x - float(expected_apix))
    return {
        "ok": delta <= float(tolerance),
        "pixel_x": hdr.pixel_x,
        "expected": float(expected_apix),
        "delta": delta,
        "tolerance": float(tolerance),
        "path": path,
    }
