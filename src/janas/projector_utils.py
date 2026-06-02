# File: projector_utils.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology


# Standard library
import os
from os import PathLike, path
import time

# Third-party
import numpy as np
import pandas as pd
import json
import logging
from typing import List, Optional, Tuple, Dict, Any
from scipy.spatial.transform import Rotation
from numpy.fft import fftn, ifftn, fftfreq
import struct
from dataclasses import dataclass
from scipy.ndimage import gaussian_filter1d,gaussian_filter, distance_transform_edt

import csv

# Local
import janas.janas_core as janas_core
from janas import starHandler
from janas import IO_utils
from janas import utils

import time

import sys
import time
import numpy as np  # se non già importato



def addSmoothEdge_2D_fast(MI_2d, edge_pixels):
    """
    MI_2d: 2D float32 array with values 0 or 1 (projected mask).
    edge_pixels: radius of soft edge in pixels (small integer, e.g. 2–5).

    Returns: 2D float32 array with values in [0,1].
    """
    M = np.asarray(MI_2d, dtype=np.float32)

    if edge_pixels <= 0:
        return M

    # Build a simple 1D triangular kernel of length 2*edge_pixels+1
    r = int(edge_pixels)
    x = np.arange(-r, r+1, dtype=np.float32)
    # triangle: 1 at centre, linearly down to 0 at edges
    k = (r + 1 - np.abs(x)) / (r + 1)
    k /= k.sum()

    # Convolve rows then columns (separable)
    # Convolve along axis 1 (x direction)
    tmp = np.apply_along_axis(lambda v: np.convolve(v, k, mode='same'), axis=1, arr=M)
    # Convolve along axis 0 (y direction)
    soft = np.apply_along_axis(lambda v: np.convolve(v, k, mode='same'), axis=0, arr=tmp)

    # Normalise to [0,1]
    maxval = soft.max()
    if maxval > 0:
        soft /= maxval
    soft = np.clip(soft, 0.0, 1.0)

    return soft.astype(np.float32)



class ProgressTimer:
    """
    Tracker con ETA. Può aggiornare in-place su una sola riga (TTY)
    o fare print multiline se non è un TTY o se disabilitato.
    """
    def __init__(self, total: int, log_every: int = 1000, progress_secs: float = 5.0,
                 label: str = "", single_line: bool = True):
        self.total = max(int(total), 1)
        self.log_every = max(int(log_every), 1)
        self.progress_secs = float(progress_secs) if progress_secs is not None else None
        self.label = label
        self.start = time.time()
        self.last_log_t = self.start
        self.count = 0
        # abilita single-line solo se stdout è un TTY e l’utente lo desidera
        self.single_line = bool(single_line) and hasattr(sys.stdout, "isatty") and sys.stdout.isatty()
        self._printed_anything = False

    def start_banner(self):
        # Con single-line non stampiamo un banner separato (evita righe in più)
        if not self.single_line:
            print(f"[{self.label}] started; total items = {self.total}")

    def should_log(self) -> bool:
        if self.count == 0:
            return True
        if self.count % self.log_every == 0:
            return True
        if self.progress_secs is not None:
            now = time.time()
            if now - self.last_log_t >= self.progress_secs:
                self.last_log_t = now
                return True
        return False

    def _fmt_seconds(self, s: float) -> str:
        if not np.isfinite(s) or s >= 9e6:
            return "∞"
        s = int(round(s))
        h, r = divmod(s, 3600); m, s = divmod(r, 60)
        if h: return f"{h}h{m:02d}m{s:02d}s"
        if m: return f"{m}m{s:02d}s"
        return f"{s}s"

    def _render_line(self) -> str:
        now = time.time()
        dt = now - self.start
        rate = (self.count / dt) if dt > 0 else 0.0
        remaining = max(self.total - self.count, 0)
        eta = (remaining / rate) if rate > 0 else float("inf")
        pct = 100.0 * self.count / self.total
        return (f"[{self.label}] {self.count}/{self.total} ({pct:5.1f}%) "
                f"elapsed={self._fmt_seconds(dt)} rate={rate:,.1f}/s ETA~{self._fmt_seconds(eta)}")

    def _print_line(self, final: bool = False):
        line = self._render_line()
        if self.single_line:
            # \r per riscrivere la stessa riga; pad con spazi per cancellare residui
            width = max(getattr(self, "_last_width", 0), len(line))
            sys.stdout.write("\r" + line.ljust(width))
            sys.stdout.flush()
            self._last_width = width
            if final:
                # chiudi con newline solo alla fine
                sys.stdout.write("\n")
                sys.stdout.flush()
        else:
            print(line)
        self._printed_anything = True

    def tick(self, n: int = 1):
        self.count += n
        if self.should_log():
            self._print_line(final=False)

    def summary(self):
        # stampa riga finale “definitiva”
        self._print_line(final=True)



# ---------- Minimal MRC header/data helpers (little-endian CCP4/MRC) ----------

def _read_mrc_header(fp) -> Dict[str, Any]:
    """Read 1024-byte CCP4/MRC header at current file position (assumes LE)."""
    h = fp.read(1024)
    if len(h) < 1024:
        raise ValueError("File too short for an MRC header.")
    nx, ny, nz, mode = struct.unpack("<4i", h[0:16])
    nxstart, nystart, nzstart = struct.unpack("<3i", h[16:28])
    mx, my, mz = struct.unpack("<3i", h[28:40])
    xlen, ylen, zlen = struct.unpack("<3f", h[40:52])
    alpha, beta, gamma = struct.unpack("<3f", h[52:64])
    mapc, mapr, maps = struct.unpack("<3i", h[64:76])
    amin, amax, amean = struct.unpack("<3f", h[76:88])
    ispg, nsymbt = struct.unpack("<2i", h[88:96])

    # ORIGIN floats: standard offsets 196–208 (words 49–51)
    try:
        ox, oy, oz = struct.unpack("<3f", h[196:208])
    except struct.error:
        ox = oy = oz = 0.0

    # Stamp + machine: bytes 208–216
    stamp = h[208:212]
    mach  = h[212:216]

    # guard against bogus or missing unit cell entries
    if mx <= 0: mx = nx
    if my <= 0: my = ny
    if mz <= 0: mz = nz

    pix_x = (xlen / mx) if mx else 1.0
    pix_y = (ylen / my) if my else pix_x
    pix_z = (zlen / mz) if mz else pix_x

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
        "stamp": stamp, "mach": mach,
        "_raw": h  # retain original if ever needed
    }

def _dtype_from_mode(mode: int):
    # Standard modes: 0=int8, 1=int16, 2=float32, 3=complex int16, 4=complex float32, 6=uint16
    if mode == 2: return np.float32
    if mode == 1: return np.int16
    if mode == 0: return np.int8
    if mode == 6: return np.uint16
    # Non-standard but seen in practice
    if mode == 12: return np.float16  # some stacks use 12 to mean float16
    # Complex types not supported by this copier:
    if mode in (3, 4):
        raise ValueError(f"Unsupported complex MRC mode {mode}.")
    raise ValueError(f"Unsupported MRC mode {mode}.")

def _write_mrc_header(fp, nx: int, ny: int, nz: int, apix: float,
                      vmin: float, vmax: float, vmean: float, vrms: float,
                      mode: int = 2,
                      origin_angs: Tuple[float, float, float] = (0.0, 0.0, 0.0)):
    """Write a 1024-byte header (LE)."""
    header = bytearray(1024)
    def pw(word_off: int, fmt: str, *vals):
        struct.pack_into("<" + fmt, header, word_off * 4, *vals)
    mx, my, mz = nx, ny, nz
    # words 0–3
    pw(0, "i", nx); pw(1, "i", ny); pw(2, "i", nz); pw(3, "i", mode)
    # starts
    pw(4, "i", 0); pw(5, "i", 0); pw(6, "i", 0)
    # mx,my,mz
    pw(7, "i", mx); pw(8, "i", my); pw(9, "i", mz)
    # cell dims
    pw(10, "f", mx * apix); pw(11, "f", my * apix); pw(12, "f", mz * apix)
    # angles
    pw(13, "f", 90.0); pw(14, "f", 90.0); pw(15, "f", 90.0)
    # axis order
    pw(16, "i", 1); pw(17, "i", 2); pw(18, "i", 3)
    # stats
    pw(19, "f", vmin); pw(20, "f", vmax); pw(21, "f", vmean)
    # spacegroup & nsymbt
    pw(22, "i", 0); pw(23, "i", 0)
    # origin (words 49–51)
    ox, oy, oz = origin_angs
    pw(49, "f", float(ox)); pw(50, "f", float(oy)); pw(51, "f", float(oz))
    # 'MAP ' + machine stamp 'DA\x00\x00' (little-endian)
    header[208:212] = b"MAP "
    header[212:216] = b"DA\x00\x00"
    # rms and labels count
    pw(54, "f", vrms); pw(55, "i", 0)
    # write out
    fp.seek(0)
    fp.write(header)

def _slice_offset_bytes(hdr: Dict[str, Any], slice_index0: int) -> int:
    """Byte offset of a given z-slice plane (0-based) in a stack file."""
    nx, ny = hdr["nx"], hdr["ny"]
    nsymbt  = hdr["nsymbt"]
    elsz    = np.dtype(_dtype_from_mode(hdr["mode"])).itemsize
    return 1024 + nsymbt + slice_index0 * nx * ny * elsz


# get_pixel_spacing
def get_MRC_map_pixel_spacing(filename):
    """
    Parse an MRC header to extract pixel dimensions and compute angstrom-per-voxel.

    Args:
        filename: Path to the .mrc file.
    Returns:
        (apix_x, apix_y, apix_z): Pixel spacing along each axis.
    """
    import struct

    with open(filename, "rb") as f:
        # Read the number of columns, rows, and sections (bytes 0-11)
        nx = struct.unpack("i", f.read(4))[0]
        ny = struct.unpack("i", f.read(4))[0]
        nz = struct.unpack("i", f.read(4))[0]

        # Seek to the cell dimensions (bytes 40-51)
        f.seek(40)
        x_dim = struct.unpack("f", f.read(4))[0]
        y_dim = struct.unpack("f", f.read(4))[0]
        z_dim = struct.unpack("f", f.read(4))[0]

    # Calculate the actual pixel spacing by dividing dimensions by the number
    # of voxels
    apix_x = x_dim / nx
    apix_y = y_dim / ny
    apix_z = z_dim / nz
    return apix_x, apix_y, apix_z



# ===== Shared helpers for stack creation =====
def infer_cs_project_root_from_path(cs_path: str) -> str:
    import re
    d = path.dirname(path.abspath(cs_path))
    while True:
        base = path.basename(d)
        if not re.fullmatch(r"J\d+", base):
            return d
        d = path.dirname(d)

def resolve_stack_path_from_image_name(image_name: str, project_root: Optional[str]) -> Tuple[int, str]:
    at = image_name.find("@")
    if at < 0:
        raise ValueError(f"Invalid _rlnImageName entry (no '@'): {image_name}")
    idx_str = image_name[:at].strip()
    stack_rel = image_name[at + 1:].strip()
    try:
        image_no = int(idx_str)
    except Exception:
        raise ValueError(f"Invalid image index in _rlnImageName: {image_name}")
    if path.isabs(stack_rel):
        stack_path = stack_rel
    else:
        stack_path = path.normpath(path.join(project_root, stack_rel)) if project_root else path.normpath(stack_rel)
    return image_no - 1, stack_path



def _auto_particles_section_name(star_path: str) -> str:
    version = starHandler.infoStarFile(star_path)[2]
    return "" if version == "relion_v30" else "particles"





def _open_src_mm(src_stack: str, hdr: Dict[str, Any]) -> np.memmap:
    """Memmap the image payload as (nz, ny, nx)."""
    nx, ny, nz = int(hdr["nx"]), int(hdr["ny"]), int(hdr["nz"])
    dt = _dtype_from_mode(int(hdr["mode"]))
    off = 1024 + int(hdr["nsymbt"])
    return np.memmap(src_stack, mode="r", dtype=dt, offset=off, shape=(nz, ny, nx), order="C")

def _probe_first_stack(image_names: List[str], project_root: Optional[str]) -> Tuple[int, int, float]:
    idx0, first_stack = resolve_stack_path_from_image_name(image_names[0], project_root)
    if not path.exists(first_stack):
        raise FileNotFoundError(f"Cannot locate source stack: {first_stack}")
    with open(first_stack, "rb") as f0:
        hdr0 = _read_mrc_header(f0)
    nx, ny = int(hdr0["nx"]), int(hdr0["ny"])
    apix   = float(hdr0["pixel_x"])
    return nx, ny, apix

def create_stack_from_star(
    star_in: str,
    out_root: str,
    project_root: Optional[str] = None,
    provenance_tag: Optional[str] = None,
    zfill_width: int = 6,
    override_section_name: Optional[str] = None,
    chunk_size: int = 128,  # number of slices to move per batch
) -> Tuple[str, str]:
    """
    Build a consolidated .mrcs with minimal open files:
      - group by source stack
      - open one memmap per group, copy slices (in chunks), close it
      - keep a single output memmap open
    """
    image_tag = "_rlnImageName"
    names_df = starHandler.readColumns(star_in, [image_tag])
    image_names = list(names_df[image_tag])
    if not image_names:
        raise ValueError("No _rlnImageName entries found in STAR.")

    nx, ny, apix = _probe_first_stack(image_names, project_root)
    N = len(image_names)
    out_stack = f"{out_root}.mrcs"
    out_star  = f"{out_root}.star"

    # Pre-create output file (mode=2 float32) and map as (N, ny, nx)
    with open(out_stack, "wb") as f:
        _write_mrc_header(f, nx=nx, ny=ny, nz=N, apix=apix,
                          vmin=0.0, vmax=0.0, vmean=0.0, vrms=0.0, mode=2)
        f.seek(1024 + (N * ny * nx * 4) - 1)
        f.write(b"\0")
    out_mm = np.memmap(out_stack, mode="r+", dtype=np.float32,
                       offset=1024, shape=(N, ny, nx), order="C")

    # Build work lists grouped by source stack → [(out_idx, slice_idx0), ...]
    by_stack: Dict[str, List[Tuple[int, int]]] = {}
    for ii, img_name in enumerate(image_names):
        slice_idx0, src_stack = resolve_stack_path_from_image_name(img_name, project_root)
        if not path.exists(src_stack):
            raise FileNotFoundError(f"Missing source stack '{src_stack}' for particle '{img_name}'.")
        by_stack.setdefault(src_stack, []).append((ii, slice_idx0))

    # Running stats (vectorised)
    total_count = 0
    sum_vals    = 0.0
    sum_squares = 0.0
    vmin        = np.inf
    vmax        = -np.inf

    # Process each source stack independently: open → copy → close
    processed = 0
    for src_stack, pairs in by_stack.items():
        # Header + geometry/type validation
        with open(src_stack, "rb") as fs:
            hdr = _read_mrc_header(fs)
        if int(hdr["nx"]) != nx or int(hdr["ny"]) != ny:
            raise ValueError(
                f"Geometry mismatch: first stack is {nx}x{ny}, but {src_stack} is {hdr['nx']}x{hdr['ny']}"
            )

        mm = _open_src_mm(src_stack, hdr)  # open only for this group

        # Sort by slice index for better locality (optional but cheap)
        pairs.sort(key=lambda p: p[1])
        out_idx_all = np.fromiter((p[0] for p in pairs), dtype=np.int64)
        slice_idx_all = np.fromiter((p[1] for p in pairs), dtype=np.int64)

        # Copy in chunks to cap peak RAM
        for start in range(0, len(pairs), chunk_size):
            stop = min(start + chunk_size, len(pairs))
            out_idx_chunk   = out_idx_all[start:stop]
            slice_idx_chunk = slice_idx_all[start:stop]

            # Read a (k, ny, nx) view from the source, cast once
            batch = np.asarray(mm[slice_idx_chunk, :, :], dtype=np.float32, order="C")

            # Write directly to the correct output planes (random access)
            out_mm[out_idx_chunk, :, :] = batch

            # Stats for this chunk
            # Use nan-aware reductions just in case
            vmin = min(vmin, float(np.nanmin(batch)))
            vmax = max(vmax, float(np.nanmax(batch)))
            sum_vals    += float(np.nansum(batch))
            sum_squares += float(np.nansum(batch * batch))
            total_count += int(batch.size)

            processed += int(len(out_idx_chunk))
            if (processed % 1000) == 0 or processed == N:
                print(f"Copied {processed}/{N} particles", end="\r")

            # free batch early
            del batch

        # Close this source memmap immediately
        del mm

    # Flush and close output map
    out_mm.flush()
    del out_mm

    # Final stats
    if total_count > 0:
        vmean = sum_vals / total_count
        var = (sum_squares / total_count) - (vmean * vmean)
        var = max(var, 0.0)
        vrms = float(np.sqrt(var))
    else:
        vmean = 0.0
        vrms = 0.0

    # Rewrite header with final stats
    with open(out_stack, "r+b") as f:
        _write_mrc_header(f, nx=nx, ny=ny, nz=N, apix=apix,
                          vmin=float(vmin), vmax=float(vmax),
                          vmean=float(vmean), vrms=float(vrms), mode=2)

    # STAR rewrite
    star_df = starHandler.readStar(star_in)
    if len(star_df) != N:
        raise ValueError(f"STAR row count ({len(star_df)}) does not match image list length ({N}).")

    if provenance_tag:
        star_df[provenance_tag] = image_names
    star_df["_rlnImageName"] = [f"{str(i+1).zfill(zfill_width)}@{out_stack}" for i in range(N)]
    if provenance_tag and provenance_tag in star_df.columns:
        ordered_cols = [c for c in star_df.columns if c != provenance_tag] + [provenance_tag]
        star_df = star_df[ordered_cols]

    section_name = override_section_name if override_section_name is not None else _auto_particles_section_name(star_in)
    starHandler.update_star_columns_from_sections(
        filenameIn=star_in,
        filenameOut=out_star,
        section_name=section_name,
        df=star_df
    )

    print(f"Wrote stack: {out_stack}")
    print(f"Wrote STAR:  {out_star}")
    return out_stack, out_star
# ———————————————————————————————————————————————————————————————
# project_mask_from_star
def _precreate_and_memmap_out_stack(out_stack: str, N: int, ny: int, nx: int, apix: float) -> np.memmap:
    """
    Create float32 .mrcs with MRC2014 header via IO_utils and return an r+ memmap view (N, ny, nx).
    """
    IO_utils.precreate_mrcs(out_stack, nz=N, ny=ny, nx=nx, apix=apix)
    mm = np.memmap(out_stack, mode="r+", dtype=np.float32, offset=1024, shape=(N, ny, nx), order="C")
    return mm


#################################################
##### project_map_from_star

def project_map_from_star(
    star_in: str,
    out_root: str,
    map_3d: str,
    mask_3d: str,
    project_root: Optional[str] = None,
    zfill_width: int = 6
) -> Tuple[str, str]:
    """
    Project a 3D reference map to each particle pose and build a 2D reprojection stack (.mrcs).
    STAR is rewritten to point _rlnImageName to the new stack.

    The projector call matches assessParticles.createDiffStack:
        janas_core.projectMap(mapI, nx, ny, nz, phi, theta, psi, tx, ty, 0)
    """
    if not path.exists(map_3d):
        raise FileNotFoundError(f"Reference map not found: {map_3d}")

    # 1) Read image list and probe geometry/apix from the first source stack
    image_tag = "_rlnImageName"
    names_df = starHandler.readColumns(star_in, [image_tag])
    image_names = list(names_df[image_tag])
    if not image_names:
        raise ValueError("No _rlnImageName entries found in STAR.")

    # geometry from first referenced stack
    _, first_stack = resolve_stack_path_from_image_name(image_names[0], project_root)
    if not path.exists(first_stack):
        raise FileNotFoundError(f"Cannot locate source stack: {first_stack}")
    hdr0 = IO_utils.read_mrc_header(first_stack)
    nx, ny = int(hdr0["nx"]), int(hdr0["ny"])
    apix   = float(hdr0["pixel_x"])

    # 2) Load reference map and mask
    mapI = janas_core.ReadMRC(map_3d)
    nx_m, ny_m, nz_m = janas_core.sizeMRC(map_3d)
    map_list_flat = list(mapI)

    if not path.exists(mask_3d):
        raise FileNotFoundError(f"Mask not found: {mask_3d}")

    maskI = janas_core.ReadMRC(mask_3d)
    mask_list_flat = list(maskI)

    # 3) Pre-create output stack
    N = len(image_names)
    out_stack = f"{out_root}.mrcs"
    out_star  = f"{out_root}.star"
    out_mm = _precreate_and_memmap_out_stack(out_stack, N=N, ny=ny, nx=nx, apix=apix)

    # 4) Build pose table (Relion 3.1 Å→px conversion identical to assessParticles)
    version = starHandler.infoStarFile(star_in)[2]
    if version == "relion_v31":
        coords_full = starHandler.readColumns(
            star_in,
            [
                "_rlnImageName",
                "_rlnAngleRot",
                "_rlnAngleTilt",
                "_rlnAnglePsi",
                "_rlnOriginXAngst",
                "_rlnOriginYAngst",
                "_rlnOpticsGroup",
            ],
        )
        idx = list(range(len(coords_full)))
        coords_full["idx"] = idx
        optics = starHandler.dataOptics(star_in)[["_rlnOpticsGroup", "_rlnImagePixelSize"]]
        coords = (coords_full
                  .merge(optics, on=["_rlnOpticsGroup"])
                  .sort_values(["idx"])
                  .set_index("idx"))
        coords["_rlnOriginX"] = coords["_rlnOriginXAngst"] / coords["_rlnImagePixelSize"]
        coords["_rlnOriginY"] = coords["_rlnOriginYAngst"] / coords["_rlnImagePixelSize"]
        coords = coords.drop(columns=["_rlnOpticsGroup", "_rlnOriginXAngst", "_rlnOriginYAngst"])
    else:
        coords = starHandler.readColumns(
            star_in,
            [
                "_rlnImageName",
                "_rlnAngleRot",
                "_rlnAngleTilt",
                "_rlnAnglePsi",
                "_rlnOriginX",
                "_rlnOriginY",
            ],
        )

    # 5) Project per-particle reprojection and write planes
    out_names = []
    timer = ProgressTimer(
        total=N, log_every=log_every if 'log_every' in locals() else 1000,
        progress_secs=progress_secs if 'progress_secs' in locals() else 5.0,
        label="project_map"
    )
    timer.start_banner()
    extraPixelsEdge = 10
    for ii in range(N):
        phi   = float(coords.at[ii, "_rlnAngleRot"])
        theta = float(coords.at[ii, "_rlnAngleTilt"])
        psi   = float(coords.at[ii, "_rlnAnglePsi"])
        tx    = float(coords.at[ii, "_rlnOriginX"])
        ty    = float(coords.at[ii, "_rlnOriginY"])

        # 1) Project the 3D mask to 2D
        MI_flat = janas_core.projectMask(
            mask_list_flat, nx_m, ny_m, nz_m,
            phi, theta, psi, tx, ty, 0, 0.5
        )

        RI_flat = janas_core.projectMap_with2DMask(map_list_flat, nx_m, ny_m, nz_m,phi, theta, psi, tx, ty, 0, MI_flat, extraPixelsEdge )
#        RI_flat = janas_core.projectMap_with2DMask(map_list_flat, nx_m, ny_m, nz_m,phi, theta, psi, tx, ty, 0)
        RI_2d = np.asarray(RI_flat, dtype=np.float32).reshape(ny, nx)

        # 4) Apply soft 2D mask in Python
        out_mm[ii, :, :] = RI_2d 

        out_names.append(f"{str(ii+1).zfill(zfill_width)}@{out_stack}")
        timer.tick()

    timer.summary()


    # 6) Flush, update header stats, and write STAR
    out_mm.flush()
    del out_mm

    # final stats (optional; keep minimal)
    with open(out_stack, "r+b") as f:
        # recompute stats cheaply
        mm_ro = np.memmap(out_stack, mode="r", dtype=np.float32, offset=1024, shape=(N, ny, nx))
        vmin = float(np.nanmin(mm_ro)); vmax = float(np.nanmax(mm_ro))
        vmean = float(np.nanmean(mm_ro)); vrms = float(np.sqrt(max(0.0, np.nanvar(mm_ro))))
        del mm_ro
        IO_utils.write_mrc_header(f, nx=nx, ny=ny, nz=N, apix=apix,
                                  vmin=vmin, vmax=vmax, vmean=vmean, vrms=vrms, mode=2)

    df_all = starHandler.readStar(star_in)
    if len(df_all) != N:
        raise ValueError(f"STAR row count ({len(df_all)}) does not match image list length ({N}).")
    df_all["_rlnImageName"] = out_names
    section_name = _auto_particles_section_name(star_in)
    starHandler.update_star_columns_from_sections(
        filenameIn=star_in, filenameOut=out_star, section_name=section_name, df=df_all
    )

    print(f"Wrote stack: {out_stack}")
    print(f"Wrote STAR:  {out_star}")
    return out_stack, out_star



#################################################
import numpy as np
import pandas as pd
import timeit
from os import PathLike
from multiprocessing import Pool

import janas.janas_core as janas_core
from janas import starHandler, utils
from .assessParticles import transformCtfImage

# ---------------------------------------------------------------------
# Scalar fit in real space inside the 2D mask
# ---------------------------------------------------------------------
def _fit_scalar_in_mask(
    I2d: np.ndarray,
    P2d: np.ndarray,
    mask2d: np.ndarray,
    mask_threshold: float = 0.1,
    alpha_clip: Optional[float] = None,
) -> float:
    """
    Fit a single scalar alpha in real space inside the 2D mask region:

        I(x) ≈ alpha * P(x)   for x in mask region

    alpha = sum(I P mask) / sum(P^2 mask)
    """
    m = (mask2d > mask_threshold).astype(np.float32)
    if np.sum(m) == 0:
        return 1.0

    I_roi = I2d * m
    P_roi = P2d * m

    num = float(np.sum(I_roi * P_roi))
    den = float(np.sum(P_roi * P_roi))

    if den <= 0.0:
        alpha = 1.0
    else:
        alpha = num / den

    if alpha_clip is not None and alpha_clip > 0.0:
        alpha = max(-alpha_clip, min(alpha_clip, alpha))

    return alpha


# ---------------------------------------------------------------------
# Globals for worker processes (map + mask shared per worker)
# ---------------------------------------------------------------------
_global_mapI = None
_global_maskI = None
_global_sizeM = None


def _init_subtraction_worker(referenceMap: str, referenceMask: str):
    """
    Initialiser for each worker: load reference map and mask once per process.
    """
    global _global_mapI, _global_maskI, _global_sizeM
    _global_mapI = janas_core.ReadMRC(referenceMap)
    _global_maskI = janas_core.ReadMRC(referenceMask)
    _global_sizeM = janas_core.sizeMRC(referenceMap)


def _process_single_particle(args):
    """
    Worker function: process a single particle and return (ii, Idiff_flat).

    All heavy work (projection, optional CTF, scalar fit, subtraction)
    happens here, using the worker-local map/mask in _global_mapI/_global_maskI.
    """
    (
        ii,
        line,
        phi,
        theta,
        psi,
        tx,
        ty,
        nx,
        ny,
        useCTF,
        ctf_dict,
        angpix,
        ctfMode,
        mask_threshold,
        alpha_clip,
    ) = args

    # Parse image number and stack name
    at_pos = line.find("@")
    imageNo = int(line[:at_pos])
    stackName = line[at_pos + 1 :]

    # Original particle
    I_list = janas_core.ReadMrcSlice(stackName, imageNo - 1)
    I2d = np.array(I_list, dtype=np.float32).reshape(ny, nx)

    # 3D map/mask sizes (shared per worker)
    nx_m, ny_m, nz_m = _global_sizeM

    # 2D mask in this pose
    MI_flat = janas_core.projectMask(
        _global_maskI, nx_m, ny_m, nz_m, phi, theta, psi, tx, ty, 0, 0.5
    )
    mask2d = np.array(MI_flat, dtype=np.float32).reshape(ny, nx)
    mask2d = np.clip(mask2d, 0.0, 1.0)

    # Projection of masked map (ROI)
    RI_flat = janas_core.projectMap_with2DMask(
        _global_mapI,
        nx_m,
        ny_m,
        nz_m,
        phi,
        theta,
        psi,
        tx,
        ty,
        0,
        MI_flat,
        10,
    )
    P2d = np.array(RI_flat, dtype=np.float32).reshape(ny, nx)

    # Optional CTF on the projection
    if useCTF and ctf_dict is not None:
        P_ctf_flat = transformCtfImage(
            P2d.flatten().tolist(),
            nx,
            ny,
            angpix,
            ctf_dict["Voltage"],
            ctf_dict["DefocusU"],
            ctf_dict["DefocusV"],
            ctf_dict["DefocusAngle"],
            ctf_dict["SphericalAberration"],
            ctf_dict["CtfBfactor"],
            ctf_dict["PhaseShift"],
            ctf_dict["AmplitudeContrast"],
            ctf_dict["DetectorPixelSize"],
            ctfMode=ctfMode,
        )
        P2d = np.array(P_ctf_flat, dtype=np.float32).reshape(ny, nx)

    # Scalar fit in mask
    alpha = _fit_scalar_in_mask(
        I2d,
        P2d,
        mask2d,
        mask_threshold=mask_threshold,
        alpha_clip=alpha_clip,
    )

    # Subtraction
    Idiff = (I2d - alpha * P2d).astype(np.float32)

    return ii, Idiff.flatten().tolist()


# ---------------------------------------------------------------------
# Main API
# ---------------------------------------------------------------------
def subtractParticles_weighted(
    particlesStarFile: PathLike,
    outputBasename: str,
    referenceMap: PathLike,
    referenceMask: PathLike,
    useCTF: bool = True,
    ctfMode: str = "modulate",
    mask_threshold: float = 0.1,
    alpha_clip: Optional[float] = None,
    n_procs: int = 180,
    chunksize: int = 20,
) -> None:
    """
    Projection-weighted particle subtraction (order-0, real-space) with
    optional multiprocessing.

    New stack:  outputBasename.mrcs
    New STAR:   image indices are (ii+1)@outputBasename.mrcs, independent
                of original image numbers or stacks.

    Parameters
    ----------
    n_procs : int
        Number of worker processes to use. 1 = no parallelism (serial).
    chunksize : int
        Chunk size for Pool.imap_unordered; a value like 80 matches your
        “blocks of 80 particles at a time” idea.
    """

    start_time = timeit.default_timer()

    # Pixel size from map header (assume isotropic)
    apix_x, apix_y, apix_z = utils.get_MRC_map_pixel_spacing(str(referenceMap))
    angpix = float(apix_x)

    # 3D map dimensions (map itself will be loaded in each worker)
    nx_m, ny_m, nz_m = janas_core.sizeMRC(str(referenceMap))

    # Inspect first particle to get 2D box size and allocate output stack
    imageNameTag = "_rlnImageName"
    imageNames = starHandler.readColumns(particlesStarFile, [imageNameTag])
    firstLine = imageNames[imageNameTag][0]
    first_at_pos = firstLine.find("@")
    stackName0 = firstLine[first_at_pos + 1 :]
    nx, ny, nz_stack = janas_core.sizeMRC(stackName0)

    # Create empty output stack in STAR row order
    N = len(imageNames[imageNameTag])
    janas_core.WriteEmptyMRC(
        outputBasename + ".mrcs", nx, ny, N
    )

    # Read orientations and shifts
    version = starHandler.infoStarFile(particlesStarFile)[2]
    if version == "relion_v31":
        coords_full = starHandler.readColumns(
            particlesStarFile,
            [
                "_rlnRandomSubset",
                "_rlnImageName",
                "_rlnAngleRot",
                "_rlnAngleTilt",
                "_rlnAnglePsi",
                "_rlnOriginXAngst",
                "_rlnOriginYAngst",
                "_rlnOpticsGroup",
            ],
        )
        coords_full["idx"] = range(len(coords_full))
        optics = starHandler.dataOptics(particlesStarFile)[
            ["_rlnOpticsGroup", "_rlnImagePixelSize"]
        ]
        coordinates = (
            coords_full.merge(optics, on=["_rlnOpticsGroup"])
            .sort_values(["idx"])
            .set_index("idx")
        )
        coordinates["_rlnOriginX"] = (
            coordinates["_rlnOriginXAngst"] / coordinates["_rlnImagePixelSize"]
        )
        coordinates["_rlnOriginY"] = (
            coordinates["_rlnOriginYAngst"] / coordinates["_rlnImagePixelSize"]
        )
        coordinates = coordinates.drop(
            ["_rlnOpticsGroup", "_rlnOriginXAngst", "_rlnOriginYAngst"], axis=1
        )
    else:
        coordinates = starHandler.readColumns(
            particlesStarFile,
            [
                "_rlnRandomSubset",
                "_rlnImageName",
                "_rlnAngleRot",
                "_rlnAngleTilt",
                "_rlnAnglePsi",
                "_rlnOriginX",
                "_rlnOriginY",
            ],
        )

    # CTF parameters (if used)
    ctfParameters = pd.DataFrame([])
    if useCTF:
        columns = starHandler.header_columns(particlesStarFile)
        if "_rlnPhaseShift" not in columns:
            PhaseShift = pd.DataFrame(np.zeros(len(coordinates)))
        else:
            PhaseShift = starHandler.readColumns(
                particlesStarFile, ["_rlnPhaseShift"]
            )

        if version == "relion_v31":
            params_full = starHandler.readColumns(
                particlesStarFile,
                [
                    "_rlnImageName",
                    "_rlnDefocusU",
                    "_rlnDefocusV",
                    "_rlnDefocusAngle",
                    "_rlnOpticsGroup",
                    "_rlnCtfBfactor",
                ],
            )
            params_full["idx"] = range(len(params_full))
            optics_ctf = starHandler.dataOptics(particlesStarFile)[
                [
                    "_rlnImagePixelSize",
                    "_rlnVoltage",
                    "_rlnAmplitudeContrast",
                    "_rlnSphericalAberration",
                    "_rlnOpticsGroup",
                ]
            ]
            ctfParameters = (
                params_full.merge(optics_ctf, on=["_rlnOpticsGroup"])
                .sort_values(["idx"])
                .set_index("idx")
            )
            ctfParameters = ctfParameters.drop(["_rlnOpticsGroup"], axis=1)
            ctfParameters.rename(
                columns={"_rlnImagePixelSize": "_rlnDetectorPixelSize"}, inplace=True
            )
        else:
            ctfParameters = starHandler.readColumns(
                particlesStarFile,
                [
                    "_rlnImageName",
                    "_rlnDefocusU",
                    "_rlnDefocusV",
                    "_rlnDefocusAngle",
                    "_rlnDetectorPixelSize",
                    "_rlnVoltage",
                    "_rlnAmplitudeContrast",
                    "_rlnSphericalAberration",
                    "_rlnCtfBfactor",
                ],
            )

        ctfParameters["_rlnPhaseShift"] = PhaseShift

    numParticles = len(coordinates["_rlnImageName"])
    print(
        f"[subtractParticles_weighted] num particles = {numParticles}, "
        f"useCTF = {useCTF}, ctfMode = {ctfMode}, "
        f"n_procs = {n_procs}, chunksize = {chunksize}"
    )

    # Pre-build task list
    tasks = []
    for ii in range(numParticles):
        line = coordinates["_rlnImageName"][ii]
        phi = float(coordinates.at[ii, "_rlnAngleRot"])
        theta = float(coordinates.at[ii, "_rlnAngleTilt"])
        psi = float(coordinates.at[ii, "_rlnAnglePsi"])
        tx = float(coordinates.at[ii, "_rlnOriginX"])
        ty = float(coordinates.at[ii, "_rlnOriginY"])

        if useCTF and not ctfParameters.empty:
            row_ctf = ctfParameters.iloc[ii]
            ctf_dict = {
                "Voltage": float(row_ctf["_rlnVoltage"]),
                "DefocusU": float(row_ctf["_rlnDefocusU"]),
                "DefocusV": float(row_ctf["_rlnDefocusV"]),
                "DefocusAngle": float(row_ctf["_rlnDefocusAngle"]),
                "SphericalAberration": float(row_ctf["_rlnSphericalAberration"]),
                "CtfBfactor": float(row_ctf["_rlnCtfBfactor"]),
                "PhaseShift": float(row_ctf["_rlnPhaseShift"]),
                "AmplitudeContrast": float(row_ctf["_rlnAmplitudeContrast"]),
                "DetectorPixelSize": float(row_ctf["_rlnDetectorPixelSize"]),
            }
        else:
            ctf_dict = None

        tasks.append(
            (
                ii,
                line,
                phi,
                theta,
                psi,
                tx,
                ty,
                nx,
                ny,
                useCTF,
                ctf_dict,
                angpix,
                ctfMode,
                mask_threshold,
                alpha_clip,
            )
        )

    outImageNames = [None] * numParticles
    width = first_at_pos  # zero-padding width from original label

    # Serial version (n_procs == 1): just call worker directly
    if n_procs == 1:
        for ii, args in enumerate(tasks):
            # progress + ETA
            if ii % 10 == 0 or ii == numParticles - 1:
                elapsed = timeit.default_timer() - start_time
                done = ii + 1
                if elapsed > 0:
                    speed = done / elapsed
                    remaining = max(0, numParticles - done)
                    eta_sec = remaining / speed if speed > 0 else 0.0
                else:
                    speed = 0.0
                    eta_sec = 0.0
                eta_h = int(eta_sec // 3600)
                eta_m = int((eta_sec % 3600) // 60)
                eta_s = int(eta_sec % 60)
                print(
                    f"  processing particle {done}/{numParticles} | "
                    f"elapsed {elapsed:7.1f}s | "
                    f"speed {speed:6.2f} part/s | "
                    f"ETA {eta_h:02d}:{eta_m:02d}:{eta_s:02d}",
                    end="\r",
                    flush=True,
                )

            # For serial, we still need map/mask loaded
            _init_subtraction_worker(str(referenceMap), str(referenceMask))
            ii_res, Idiff_flat = _process_single_particle(args)

            # Write to new stack
            janas_core.ReplaceMrcSlice(
                Idiff_flat,
                outputBasename + ".mrcs",
                nx,
                ny,
                ii_res,
            )
            new_img_str = str(ii_res + 1).zfill(width)
            outImageNames[ii_res] = f"{new_img_str}@{outputBasename}.mrcs"

    else:
        # Parallel version with Pool
        with Pool(
            processes=n_procs,
            initializer=_init_subtraction_worker,
            initargs=(str(referenceMap), str(referenceMask)),
        ) as pool:
            results_iter = pool.imap_unordered(
                _process_single_particle, tasks, chunksize=chunksize
            )

            processed = 0
            for ii_res, Idiff_flat in results_iter:
                processed += 1

                # Write slice to new stack
                janas_core.ReplaceMrcSlice(
                    Idiff_flat,
                    outputBasename + ".mrcs",
                    nx,
                    ny,
                    ii_res,
                )

                # New name in STAR
                new_img_str = str(ii_res + 1).zfill(width)
                outImageNames[ii_res] = f"{new_img_str}@{outputBasename}.mrcs"

                # Progress + ETA every ~chunksize particles
                if processed % chunksize == 0 or processed == numParticles:
                    elapsed = timeit.default_timer() - start_time
                    done = processed
                    if elapsed > 0:
                        speed = done / elapsed
                        remaining = max(0, numParticles - done)
                        eta_sec = remaining / speed if speed > 0 else 0.0
                    else:
                        speed = 0.0
                        eta_sec = 0.0
                    eta_h = int(eta_sec // 3600)
                    eta_m = int((eta_sec % 3600) // 60)
                    eta_s = int(eta_sec % 60)
                    print(
                        f"  processed {done}/{numParticles} | "
                        f"elapsed {elapsed:7.1f}s | "
                        f"speed {speed:6.2f} part/s | "
                        f"ETA {eta_h:02d}:{eta_m:02d}:{eta_s:02d}",
                        end="\r",
                        flush=True,
                    )

    # Update STAR
    inputStar = starHandler.readStar(particlesStarFile)
    inputStar["_rlnImageName"] = outImageNames
    starHandler.writeDataframeToStar(
        particlesStarFile, outputBasename + ".star", inputStar
    )

    # Final timing
    elapsed = timeit.default_timer() - start_time
    hours = int(elapsed // 3600)
    minutes = int((elapsed % 3600) // 60)
    seconds = elapsed % 60
    if elapsed > 0:
        speed = numParticles / elapsed
        print(
            f"\n[subtractParticles_weighted] Elapsed "
            f"{hours:02d}h:{minutes:02d}m:{seconds:05.2f}s "
            f"({speed:.2f} particles/s)"
        )
    else:
        print("[subtractParticles_weighted] Elapsed time: < 1e-6 s")


# ---------------------------------------------------------------------
# RELION-like subtraction: map * mask in 3D, project, CTF, subtract
# ---------------------------------------------------------------------

_relion_mapMasked = None
_relion_sizeM = None


def _init_subtraction_worker_relionLike():
    """
    Initialiser for each worker (RELION-like).

    On Linux with 'fork', the parent already loaded _relion_mapMasked and
    _relion_sizeM; workers just inherit them. This function is mainly here
    to make it explicit and to fail fast if globals are missing.
    """
    global _relion_mapMasked, _relion_sizeM
    if _relion_mapMasked is None or _relion_sizeM is None:
        raise RuntimeError(
            "RELION-like worker initialised without _relion_mapMasked/_relion_sizeM. "
            "Make sure subtractParticles_weighted_relionLike() initialised them in the parent."
        )


def _process_single_particle_relionLike(args):
    """
    Worker function: RELION-like subtraction for a single particle.

    Steps:
      1. Read particle image I2d.
      2. Project masked 3D map (global _relion_mapMasked) with projectMap.
      3. Optionally apply CTF to the projection (transformCtfImage).
      4. Subtract projection from I2d: Idiff = I2d - P2d.
    """
    (
        ii,
        line,
        phi,
        theta,
        psi,
        tx,
        ty,
        nx,
        ny,
        useCTF,
        ctf_dict,
        angpix,
        ctfMode,
    ) = args

    # Parse image number and stack name
    at_pos = line.find("@")
    imageNo = int(line[:at_pos])
    stackName = line[at_pos + 1 :]

    # Original particle (raw, contains experimental CTF)
    I_list = janas_core.ReadMrcSlice(stackName, imageNo - 1)
    I2d = np.array(I_list, dtype=np.float32).reshape(ny, nx)

    # 3D masked map size
    nx_m, ny_m, nz_m = _relion_sizeM

    # Project 3D masked map directly (RELION-style)
    RI_flat = janas_core.projectMap(
        _relion_mapMasked,
        nx_m,
        ny_m,
        nz_m,
        phi,
        theta,
        psi,
        tx,
        ty,
        0.0,   # centreZ
    )
    P2d = np.array(RI_flat, dtype=np.float32).reshape(ny, nx)

    # Optional CTF on the projection only
    if useCTF and ctf_dict is not None:
        P_ctf_flat = transformCtfImage(
            P2d.flatten().tolist(),
            nx,
            ny,
            angpix,
            ctf_dict["Voltage"],
            ctf_dict["DefocusU"],
            ctf_dict["DefocusV"],
            ctf_dict["DefocusAngle"],
            ctf_dict["SphericalAberration"],
            ctf_dict["CtfBfactor"],
            ctf_dict["PhaseShift"],
            ctf_dict["AmplitudeContrast"],
            ctf_dict["DetectorPixelSize"],
            ctfMode=ctfMode,
        )
        P2d = np.array(P_ctf_flat, dtype=np.float32).reshape(ny, nx)

    # Direct subtraction (no scalar fitting)
    Idiff = (I2d - P2d).astype(np.float32)

    return ii, Idiff.flatten().tolist()


def subtractParticles_weighted_relionLike(
    particlesStarFile: PathLike,
    outputBasename: str,
    referenceMap: PathLike,
    referenceMask: PathLike,
    useCTF: bool = True,
    ctfMode: str = "modulate",
    mask_threshold: float = 0.1,      # kept for API compatibility, unused
    alpha_clip: Optional[float] = None,  # kept for API compatibility, unused
    n_procs: int = 8,
    chunksize: int = 20,
    **kwargs,
) -> None:
    """
    RELION-like projection-weighted particle subtraction.

    For each particle:
      1. Project the 3D masked map (map * mask in 3D).
      2. Optionally apply CTF to the projection.
      3. Subtract the projection from the raw particle.

    New stack:  outputBasename.mrcs
    New STAR:   entries are (ii+1)@outputBasename.mrcs in STAR row order.
    """
    useCTF=True
    n_procs=80
    chunksize=80

    start_time = timeit.default_timer()

    # -------------------------
    # 0) Pre-load masked map ONCE in the parent
    # -------------------------
    global _relion_mapMasked, _relion_sizeM

    # Read map and mask and build masked map in parent
    mapI = janas_core.ReadMRC(str(referenceMap))
    maskI = janas_core.ReadMRC(str(referenceMask))
    nx_m, ny_m, nz_m = janas_core.sizeMRC(str(referenceMap))
    nxyz = nx_m * ny_m * nz_m

    if len(mapI) != nxyz or len(maskI) != nxyz:
        raise ValueError(
            f"Map ({len(mapI)}) or mask ({len(maskI)}) size mismatch with "
            f"nx*ny*nz = {nxyz}"
        )

    map_arr = np.array(mapI, dtype=np.float32)
    mask_arr = np.array(maskI, dtype=np.float32)
    map_masked = (map_arr * mask_arr).astype(np.float32)

    _relion_mapMasked = map_masked.tolist()
    _relion_sizeM = (nx_m, ny_m, nz_m)

    # Pixel size from map header (assume isotropic)
    apix_x, apix_y, apix_z = utils.get_MRC_map_pixel_spacing(str(referenceMap))
    angpix = float(apix_x)

    # Inspect first particle to get 2D box size and allocate output stack
    imageNameTag = "_rlnImageName"
    imageNames = starHandler.readColumns(particlesStarFile, [imageNameTag])
    firstLine = imageNames[imageNameTag][0]
    first_at_pos = firstLine.find("@")
    stackName0 = firstLine[first_at_pos + 1 :]
    nx, ny, nz_stack = janas_core.sizeMRC(stackName0)

    # Create empty output stack in STAR row order
    N = len(imageNames[imageNameTag])
    janas_core.WriteEmptyMRC(
        outputBasename + ".mrcs", nx, ny, N
    )

    # Read orientations and shifts
    version = starHandler.infoStarFile(particlesStarFile)[2]
    if version == "relion_v31":
        coords_full = starHandler.readColumns(
            particlesStarFile,
            [
                "_rlnRandomSubset",
                "_rlnImageName",
                "_rlnAngleRot",
                "_rlnAngleTilt",
                "_rlnAnglePsi",
                "_rlnOriginXAngst",
                "_rlnOriginYAngst",
                "_rlnOpticsGroup",
            ],
        )
        coords_full["idx"] = range(len(coords_full))
        optics = starHandler.dataOptics(particlesStarFile)[
            ["_rlnOpticsGroup", "_rlnImagePixelSize"]
        ]
        coordinates = (
            coords_full.merge(optics, on=["_rlnOpticsGroup"])
            .sort_values(["idx"])
            .set_index("idx")
        )
        coordinates["_rlnOriginX"] = (
            coordinates["_rlnOriginXAngst"] / coordinates["_rlnImagePixelSize"]
        )
        coordinates["_rlnOriginY"] = (
            coordinates["_rlnOriginYAngst"] / coordinates["_rlnImagePixelSize"]
        )
        coordinates = coordinates.drop(
            ["_rlnOpticsGroup", "_rlnOriginXAngst", "_rlnOriginYAngst"], axis=1
        )
    else:
        coordinates = starHandler.readColumns(
            particlesStarFile,
            [
                "_rlnRandomSubset",
                "_rlnImageName",
                "_rlnAngleRot",
                "_rlnAngleTilt",
                "_rlnAnglePsi",
                "_rlnOriginX",
                "_rlnOriginY",
            ],
        )

    # CTF parameters (if used)
    ctfParameters = pd.DataFrame([])
    if useCTF:
        columns = starHandler.header_columns(particlesStarFile)
        if "_rlnPhaseShift" not in columns:
            PhaseShift = pd.DataFrame(np.zeros(len(coordinates)))
        else:
            PhaseShift = starHandler.readColumns(
                particlesStarFile, ["_rlnPhaseShift"]
            )

        if version == "relion_v31":
            params_full = starHandler.readColumns(
                particlesStarFile,
                [
                    "_rlnImageName",
                    "_rlnDefocusU",
                    "_rlnDefocusV",
                    "_rlnDefocusAngle",
                    "_rlnOpticsGroup",
                    "_rlnCtfBfactor",
                ],
            )
            params_full["idx"] = range(len(params_full))
            optics_ctf = starHandler.dataOptics(particlesStarFile)[
                [
                    "_rlnImagePixelSize",
                    "_rlnVoltage",
                    "_rlnAmplitudeContrast",
                    "_rlnSphericalAberration",
                    "_rlnOpticsGroup",
                ]
            ]
            ctfParameters = (
                params_full.merge(optics_ctf, on=["_rlnOpticsGroup"])
                .sort_values(["idx"])
                .set_index("idx")
            )
            ctfParameters = ctfParameters.drop(["_rlnOpticsGroup"], axis=1)
            ctfParameters.rename(
                columns={"_rlnImagePixelSize": "_rlnDetectorPixelSize"}, inplace=True
            )
        else:
            ctfParameters = starHandler.readColumns(
                particlesStarFile,
                [
                    "_rlnImageName",
                    "_rlnDefocusU",
                    "_rlnDefocusV",
                    "_rlnDefocusAngle",
                    "_rlnDetectorPixelSize",
                    "_rlnVoltage",
                    "_rlnAmplitudeContrast",
                    "_rlnSphericalAberration",
                    "_rlnCtfBfactor",
                ],
            )

        ctfParameters["_rlnPhaseShift"] = PhaseShift

    numParticles = len(coordinates["_rlnImageName"])
    print(
        f"[subtractParticles_weighted_relionLike] num particles = {numParticles}, "
        f"useCTF = {useCTF}, ctfMode = {ctfMode}, "
        f"n_procs = {n_procs}, chunksize = {chunksize}"
    )

    # Build tasks list
    tasks = []
    for ii in range(numParticles):
        line = coordinates["_rlnImageName"][ii]
        phi = float(coordinates.at[ii, "_rlnAngleRot"])
        theta = float(coordinates.at[ii, "_rlnAngleTilt"])
        psi = float(coordinates.at[ii, "_rlnAnglePsi"])
        tx = float(coordinates.at[ii, "_rlnOriginX"])
        ty = float(coordinates.at[ii, "_rlnOriginY"])

        if useCTF and not ctfParameters.empty:
            row_ctf = ctfParameters.iloc[ii]
            ctf_dict = {
                "Voltage": float(row_ctf["_rlnVoltage"]),
                "DefocusU": float(row_ctf["_rlnDefocusU"]),
                "DefocusV": float(row_ctf["_rlnDefocusV"]),
                "DefocusAngle": float(row_ctf["_rlnDefocusAngle"]),
                "SphericalAberration": float(row_ctf["_rlnSphericalAberration"]),
                "CtfBfactor": float(row_ctf["_rlnCtfBfactor"]),
                "PhaseShift": float(row_ctf["_rlnPhaseShift"]),
                "AmplitudeContrast": float(row_ctf["_rlnAmplitudeContrast"]),
                "DetectorPixelSize": float(row_ctf["_rlnDetectorPixelSize"]),
            }
        else:
            ctf_dict = None

        tasks.append(
            (
                ii,
                line,
                phi,
                theta,
                psi,
                tx,
                ty,
                nx,
                ny,
                useCTF,
                ctf_dict,
                angpix,
                ctfMode,
            )
        )

    outImageNames = [None] * numParticles
    width = first_at_pos  # zero-padding width from original label

    # Reasonable default: do not oversubscribe cores too aggressively
    from multiprocessing import cpu_count
    if n_procs is None or n_procs <= 0:
        n_procs = min(cpu_count(), 16)

    # Parallel version with Pool
    from multiprocessing import Pool
    with Pool(
        processes=n_procs,
        initializer=_init_subtraction_worker_relionLike,
    ) as pool:
        results_iter = pool.imap_unordered(
            _process_single_particle_relionLike, tasks, chunksize=chunksize
        )

        processed = 0
        for ii_res, Idiff_flat in results_iter:
            processed += 1

            # Write slice to new stack
            janas_core.ReplaceMrcSlice(
                Idiff_flat,
                outputBasename + ".mrcs",
                nx,
                ny,
                ii_res,
            )

            # New name in STAR: (ii+1)@outputBasename.mrcs
            new_img_str = str(ii_res + 1).zfill(width)
            outImageNames[ii_res] = f"{new_img_str}@{outputBasename}.mrcs"

            # Progress + ETA
            if processed % chunksize == 0 or processed == numParticles:
                elapsed = timeit.default_timer() - start_time
                done = processed
                if elapsed > 0:
                    speed = done / elapsed
                    remaining = max(0, numParticles - done)
                    eta_sec = remaining / speed if speed > 0 else 0.0
                else:
                    speed = 0.0
                    eta_sec = 0.0
                eta_h = int(eta_sec // 3600)
                eta_m = int((eta_sec % 3600) // 60)
                eta_s = int(eta_sec % 60)
                print(
                    f"  processed {done}/{numParticles} | "
                    f"elapsed {elapsed:7.1f}s | "
                    f"speed {speed:6.2f} part/s | "
                    f"ETA {eta_h:02d}:{eta_m:02d}:{eta_s:02d}",
                    end="\r",
                    flush=True,
                )

    # Update STAR
    inputStar = starHandler.readStar(particlesStarFile)
    inputStar["_rlnImageName"] = outImageNames
    starHandler.writeDataframeToStar(
        particlesStarFile, outputBasename + ".star", inputStar
    )

    # Final timing
    elapsed = timeit.default_timer() - start_time
    hours = int(elapsed // 3600)
    minutes = int((elapsed % 3600) // 60)
    seconds = elapsed % 60
    if elapsed > 0:
        speed = numParticles / elapsed
        print(
            f"\n[subtractParticles_weighted_relionLike] Elapsed "
            f"{hours:02d}h:{minutes:02d}m:{seconds:05.2f}s "
            f"({speed:.2f} particles/s)"
        )
    else:
        print("[subtractParticles_weighted_relionLike] Elapsed time: < 1e-6 s")

