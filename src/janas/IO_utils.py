#!/usr/bin/env python3
"""
MRC I/O utilities (little-endian, MRC2014-compatible).

Public API:
- write_mrc(path, volume, pixel_size, origin_A=(0,0,0), update_stats=True)
- read_mrc_header(path) -> dict
- read_mrc_data(path, astype=np.float32) -> (volume, header_dict)
- write_mrc_like(path, volume, template_hdr, update_stats=True)

Notes
-----
* Voxel order is (Z, Y, X).
* ORIGIN is stored in Ångströms (origin_A), not angles.
* We write mode=2 (float32) on output.
"""

from __future__ import annotations

import struct
from typing import Dict, Tuple, Optional

import numpy as np

__all__ = [
    "write_mrc",
    "read_mrc_header",
    "read_mrc_data",
    "write_mrc_like",
    "mrc_dtype_from_mode",
]

# ----------------------- MRC I/O -----------------------
def memmap_mrc_stack(path: str) -> Tuple[np.memmap, dict]:
    hdr = IO_utils.read_mrc_header(path)
    nx, ny, nz, mode = hdr["nx"], hdr["ny"], hdr["nz"], hdr["mode"]
    if mode == 2:
        dtype = np.float32
    elif mode == 1:
        dtype = np.int16
    elif mode == 0:
        dtype = np.int8
    elif mode == 6:
        dtype = np.uint16
    else:
        raise RuntimeError(f"Unsupported MRC mode={mode} for stack {path}")
    mm = np.memmap(path, mode="r", dtype=dtype, offset=1024, shape=(nz, ny, nx))
    return mm, hdr

def mrc_dtype_from_mode(mode: int) -> np.dtype:
    """
    Map MRC mode to numpy dtype.

    Supported modes:
      0: int8
      1: int16
      2: float32
      6: uint16

    Other modes raise NotImplementedError.
    """
    if mode == 0:
        return np.int8
    if mode == 1:
        return np.int16
    if mode == 2:
        return np.float32
    if mode == 6:
        return np.uint16
    raise NotImplementedError(f"MRC mode {mode} is not supported by this reader.")

def write_mrc(path: str,
              volume: np.ndarray,
              pixel_size: float,
              origin_angs: Tuple[float, float, float] = (0.0, 0.0, 0.0),
              update_stats: bool = True):
    vol = np.asarray(volume, dtype=np.float32)
    nz, ny, nx = vol.shape
    mode = 2
    mx, my, mz = nx, ny, nz

    if update_stats:
        vmin = float(np.nanmin(vol))
        vmax = float(np.nanmax(vol))
        vmean = float(np.nanmean(vol))
        vstd = float(np.nanstd(vol))
    else:
        vmin = vmax = vmean = vstd = 0.0

    header = bytearray(1024)
    def pw(word_off: int, fmt: str, *vals):
        struct.pack_into("<" + fmt, header, word_off * 4, *vals)

    pw(0, "i", nx); pw(1, "i", ny); pw(2, "i", nz); pw(3, "i", mode)
    pw(4, "i", 0);  pw(5, "i", 0);  pw(6, "i", 0)
    pw(7, "i", mx); pw(8, "i", my); pw(9, "i", mz)
    pw(10, "f", mx * pixel_size); pw(11, "f", my * pixel_size); pw(12, "f", mz * pixel_size)
    pw(13, "f", 90.0); pw(14, "f", 90.0); pw(15, "f", 90.0)
    pw(16, "i", 1); pw(17, "i", 2); pw(18, "i", 3)
    pw(19, "f", vmin); pw(20, "f", vmax); pw(21, "f", vmean)
    pw(22, "i", 0); pw(23, "i", 0)
    ox, oy, oz = origin_angs
    pw(49, "f", float(ox)); pw(50, "f", float(oy)); pw(51, "f", float(oz))
    header[208:212] = b"MAP "; header[212:216] = b"DA\x00\x00"
    pw(54, "f", vstd); pw(55, "i", 0)

    with open(path, "wb") as f:
        f.write(header)
        f.write(vol.tobytes())
        
        
def read_mrc_header(path: str) -> Dict[str, float]:
    with open(path, "rb") as f:
        h = f.read(1024)
    if len(h) < 1024:
        raise ValueError("File too short to contain an MRC header.")
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
    except struct.error:
        ox, oy, oz = 0.0, 0.0, 0.0
    stamp = h[208:212]

    if mx <= 0: mx = nx
    if my <= 0: my = ny
    if mz <= 0: mz = nz

    pix_x = xlen / mx if mx else 1.0
    pix_y = ylen / my if my else pix_x
    pix_z = zlen / mz if mz else pix_x

    return {
        "nx": nx, "ny": ny, "nz": nz, "mode": mode,
        "nxstart": nxstart, "nystart": nystart, "nzstart": nzstart,
        "mx": mx, "my": my, "mz": mz,
        "xlen": float(xlen), "ylen": float(ylen), "zlen": float(zlen),
        "alpha": float(alpha), "beta": float(beta), "gamma": float(gamma),
        "mapc": mapc, "mapr": mapr, "maps": maps,
        "amin": float(amin), "amax": float(amax), "amean": float(amean),
        "ispg": ispg, "nsymbt": nsymbt,
        "origin_x": float(ox), "origin_y": float(oy), "origin_z": float(oz),
        "pixel_x": float(pix_x), "pixel_y": float(pix_y), "pixel_z": float(pix_z),
        "stamp": stamp.decode(errors="ignore")
    }


def write_mrc_header(fh, *, nx, ny, nz, apix, vmin=0.0, vmax=0.0, vmean=0.0, vrms=0.0, mode=2):
    # Minimal MRC header (1024 bytes). mode=2 => float32
    # Only fields we actually need are filled.
    import struct
    header = bytearray(1024)
    struct.pack_into("<4i", header, 0, nx, ny, nz, mode)
    # nxstart,nystart,nzstart
    struct.pack_into("<3i", header, 16, 0, 0, 0)
    # mx,my,mz (grid)
    struct.pack_into("<3i", header, 28, nx, ny, nz)
    # cell dimensions (Å) — isotropic pixels
    struct.pack_into("<3f", header, 40, nx*apix, ny*apix, nz*apix)
    # cell angles
    struct.pack_into("<3f", header, 52, 90.0, 90.0, 90.0)
    # axis mapping
    struct.pack_into("<3i", header, 64, 1, 2, 3)
    # min, max, mean
    struct.pack_into("<3f", header, 76, float(vmin), float(vmax), float(vmean))
    # rms (non-standard location; many tools ignore — leave 0 is fine)
    # spacegroup, extra...
    fh.seek(0)
    fh.write(header)

def precreate_mrcs(path_out: str, *, nz: int, ny: int, nx: int, apix: float) -> None:
    # Create header + allocate body (float32), without loading to RAM
    with open(path_out, "wb") as f:
        write_mrc_header(f, nx=nx, ny=ny, nz=nz, apix=apix, mode=2)
        f.seek(1024 + (nz * ny * nx * 4) - 1)
        f.write(b"\0")



def read_mrc_data(path: str, astype=np.float32) -> Tuple[np.ndarray, Dict[str, float]]:
    hdr = read_mrc_header(path)
    nx, ny, nz, mode = hdr["nx"], hdr["ny"], hdr["nz"], hdr["mode"]
    dt = mrc_dtype_from_mode(mode)
    nvox = nx * ny * nz
    with open(path, "rb") as f:
        f.seek(1024 + hdr["nsymbt"])
        arr = np.fromfile(f, dtype=dt, count=nvox)
    if arr.size != nvox:
        raise ValueError("MRC data size inconsistent with header.")
    vol = arr.reshape((nz, ny, nx))
    if dt != np.float32:
        vol = vol.astype(np.float32)
    if (astype is not None) and (astype is not np.float32):
        vol = vol.astype(astype)
    return vol, hdr
    
def write_mrc_like(path: str,
                   volume: np.ndarray,
                   template_hdr: Dict[str, float],
                   update_stats: bool = True):
    vol = np.asarray(volume, dtype=np.float32)
    nz, ny, nx = vol.shape
    if (nx, ny, nz) != (template_hdr["nx"], template_hdr["ny"], template_hdr["nz"]):
        raise ValueError("write_mrc_like: volume shape does not match template header shape.")

    if update_stats:
        vmin = float(np.nanmin(vol))
        vmax = float(np.nanmax(vol))
        vmean = float(np.nanmean(vol))
        vstd = float(np.nanstd(vol))
    else:
        vmin = vmax = vmean = vstd = 0.0

    header = bytearray(1024)

    def pw(word_off: int, fmt: str, *vals):
        struct.pack_into("<" + fmt, header, word_off * 4, *vals)

    pw(0, "i", nx); pw(1, "i", ny); pw(2, "i", nz); pw(3, "i", 2)
    pw(4, "i", int(template_hdr["nxstart"]))
    pw(5, "i", int(template_hdr["nystart"]))
    pw(6, "i", int(template_hdr["nzstart"]))
    pw(7, "i", int(template_hdr["mx"]))
    pw(8, "i", int(template_hdr["my"]))
    pw(9, "i", int(template_hdr["mz"]))
    pw(10, "f", float(template_hdr["xlen"]))
    pw(11, "f", float(template_hdr["ylen"]))
    pw(12, "f", float(template_hdr["zlen"]))
    pw(13, "f", float(template_hdr["alpha"]))
    pw(14, "f", float(template_hdr["beta"]))
    pw(15, "f", float(template_hdr["gamma"]))
    pw(16, "i", int(template_hdr["mapc"]))
    pw(17, "i", int(template_hdr["mapr"]))
    pw(18, "i", int(template_hdr["maps"]))
    pw(19, "f", vmin); pw(20, "f", vmax); pw(21, "f", vmean)
    pw(22, "i", int(template_hdr["ispg"]))
    pw(23, "i", int(template_hdr["nsymbt"]))
    pw(49, "f", float(template_hdr["origin_x"]))
    pw(50, "f", float(template_hdr["origin_y"]))
    pw(51, "f", float(template_hdr["origin_z"]))
    header[208:212] = b"MAP "
    header[212:216] = b"DA\x00\x00"
    pw(54, "f", vstd)
    pw(55, "i", 0)

    with open(path, "wb") as f:
        f.write(header)
        f.write(vol.tobytes())
