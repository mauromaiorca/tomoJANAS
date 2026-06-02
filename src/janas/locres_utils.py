#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# locres_utils.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
Local-resolution utilities for JANAS.

This module factors out the local FSC / local-resolution machinery from
janas_reconstructor.py so it can be reused by both:
  - janas_reconstructor.py  (post-reconstruction half-map analysis)
  - janas_cmd_utils.py locres (standalone CLI)

Public API:
  - LocResConfig (dataclass)
  - locres_map(cfg: LocResConfig) -> Dict[str, Any]
  - run_locres_for_all_pairs(results, preferred_dir, args) -> None
      (drop-in replacement for the reconstructor’s previous helper)

"""

from __future__ import annotations

# stdlib
import os
import re
import math
import time
from pathlib import Path
from dataclasses import dataclass, replace
from typing import Any, Dict, List, Optional, Tuple, Union

# third-party
import csv
import numpy as np
from numpy.fft import fftn, ifftn, fftfreq
from scipy.ndimage import binary_dilation,binary_erosion,  zoom as ndi_zoom, distance_transform_edt
from scipy.interpolate import interp1d
import multiprocessing as mp
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, uniform_filter
from scipy.ndimage import distance_transform_edt, uniform_filter, gaussian_filter, median_filter

# local
from janas import IO_utils
from janas import janas_mapProcess


# ----------------------- Constants -----------------------
AUX_THR = 0.143  # fixed auxiliary threshold for auto-radius and tail stats


# -----------------------------
# Globals for forked workers
# -----------------------------
_LOCINFO_F1B = None  # list[np.ndarray], band-filtered halfmap1 in real space
_LOCINFO_F2B = None  # list[np.ndarray], band-filtered halfmap2 in real space
_LOCINFO_KC  = None  # np.ndarray of band centre frequencies (Å^-1)
_LOCINFO_W   = None  # 3D window weights (float32)
_LOCINFO_R   = None  # int radius vox (window half-size)
_LOCINFO_K2B = None # np.ndarray, per-band mean of k^2 (cycles^2 / Å^2)
_LOCINFO_EPS = 1e-12


def _build_soft_sphere_window(R_vox: int, edge_vox: int) -> np.ndarray:
    R = int(R_vox)
    ew = max(int(edge_vox), 1)
    size = 2 * R + 1
    cz = cy = cx = R

    zz, yy, xx = np.ogrid[:size, :size, :size]
    rr = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2 + (zz - cz) ** 2)

    w = np.zeros((size, size, size), dtype=np.float32)
    r0 = max(R - ew, 0)
    inside = rr <= r0
    w[inside] = 1.0

    shell = (rr > r0) & (rr <= R)
    if np.any(shell):
        t = (rr[shell] - r0) / max(R - r0, 1e-6)
        w[shell] = (0.5 * (1.0 + np.cos(np.pi * t))).astype(np.float32)

    return w


def _fft_freq_grid(shape: Tuple[int, int, int], pixel_A: float) -> np.ndarray:
    nz, ny, nx = shape
    fx = np.fft.fftfreq(nx, d=pixel_A)
    fy = np.fft.fftfreq(ny, d=pixel_A)
    fz = np.fft.fftfreq(nz, d=pixel_A)
    kz, ky, kx = np.meshgrid(fz, fy, fx, indexing="ij")
    kmag = np.sqrt(kx * kx + ky * ky + kz * kz).astype(np.float32)
    return kmag


def _make_radial_bands(kmag: np.ndarray,
                      nyquist: float,
                      n_bands: int = 6,
                      k_min: float = 0.0,
                      k_max: Optional[float] = None) -> Tuple[List[np.ndarray], np.ndarray]:
    if k_max is None:
        k_max = float(nyquist)

    k_min = max(float(k_min), 0.0)
    k_max = min(float(k_max), float(nyquist))
    if k_max <= k_min:
        raise ValueError(f"Invalid band range: k_min={k_min} k_max={k_max} (nyq={nyquist})")

    edges = np.linspace(k_min, k_max, int(n_bands) + 1, dtype=np.float64)
    masks = []
    centres = []
    for i in range(len(edges) - 1):
        a, b = float(edges[i]), float(edges[i + 1])
        if b <= a:
            continue
        m = ((kmag >= a) & (kmag < b)).astype(np.float32)
        if m.sum() == 0:
            continue
        masks.append(m)
        centres.append(0.5 * (a + b))
    return masks, np.asarray(centres, dtype=np.float32)


def _prefilter_bands_realspace(vol: np.ndarray,
                              band_masks: List[np.ndarray]) -> List[np.ndarray]:
    F = np.fft.fftn(vol)
    out = []
    for B in band_masks:
        fb = np.fft.ifftn(F * B).real.astype(np.float32, copy=False)
        out.append(fb)
    return out


def _extract_patch_with_window(vol: np.ndarray,
                               W: np.ndarray,
                               cx: int, cy: int, cz: int,
                               R: int) -> Tuple[np.ndarray, np.ndarray]:
    nz, ny, nx = vol.shape
    x0 = max(cx - R, 0); x1 = min(cx + R + 1, nx)
    y0 = max(cy - R, 0); y1 = min(cy + R + 1, ny)
    z0 = max(cz - R, 0); z1 = min(cz + R + 1, nz)

    wx0 = x0 - (cx - R); wx1 = wx0 + (x1 - x0)
    wy0 = y0 - (cy - R); wy1 = wy0 + (y1 - y0)
    wz0 = z0 - (cz - R); wz1 = wz0 + (z1 - z0)

    patch = vol[z0:z1, y0:y1, x0:x1]
    wpatch = W[wz0:wz1, wy0:wy1, wx0:wx1]
    return patch, wpatch


def _compute_local_info_radius(cx: int, cy: int, cz: int) -> Optional[Dict[str, float]]:
    global _LOCINFO_F1B, _LOCINFO_F2B, _LOCINFO_K2B, _LOCINFO_W, _LOCINFO_R, _LOCINFO_EPS

    R = _LOCINFO_R
    W = _LOCINFO_W
    k2b = _LOCINFO_K2B

    # Collect per-band stats first
    fsc_list = []
    A_list = []
    k2_list = []

    for i, (f1, f2) in enumerate(zip(_LOCINFO_F1B, _LOCINFO_F2B)):
        p1, w1 = _extract_patch_with_window(f1, W, cx, cy, cz, R)
        p2, _  = _extract_patch_with_window(f2, W, cx, cy, cz, R)

        p1f = p1.astype(np.float64, copy=False)
        p2f = p2.astype(np.float64, copy=False)
        wf  = w1.astype(np.float64, copy=False)

        e1 = float(np.sum(wf * (p1f * p1f)))
        e2 = float(np.sum(wf * (p2f * p2f)))
        c  = float(np.sum(wf * (p1f * p2f)))

        denom = math.sqrt(max(e1 * e2, 0.0)) + _LOCINFO_EPS
        fsc = c / denom
        if not np.isfinite(fsc):
            continue
        fsc = max(min(fsc, 0.999999), -0.999999)

        k2 = float(k2b[i])
        if not (np.isfinite(k2) and k2 > 0.0):
            continue

        # Cross-amplitude proxy (signal-like component)
        A = max(c, _LOCINFO_EPS)

        fsc_list.append(float(fsc))
        A_list.append(float(A))
        k2_list.append(float(k2))

    if len(k2_list) < 3:
        return None

    k2_arr = np.asarray(k2_list, dtype=np.float64)
    fsc_arr = np.asarray(fsc_list, dtype=np.float64)
    A_arr  = np.asarray(A_list, dtype=np.float64)

    # ---- Soft detectability weight from FSC (no hard chopping) ----
    # Power mapping (simple, stable)
    p = 2.0
    fpos = np.maximum(fsc_arr, 0.0)
    w_fsc = np.power(fpos, p)

    # ---- Estimate blur/heterogeneity slope from log amplitude falloff ----
    # Fit only on upper-k bands to avoid low-frequency dominance
    hi = k2_arr >= np.percentile(k2_arr, 50.0)
    if np.count_nonzero(hi) >= 3:
        x = k2_arr[hi]
        y = np.log(A_arr[hi])
        # simple least-squares slope: y ~ a + s*x
        s = float(np.polyfit(x, y, 1)[0])
        alpha = max(-s, 0.0)  # alpha >= 0 means amplitude decays with k^2
    else:
        alpha = 0.0

    # ---- Combine weights: FSC detectability × blur penalty ----
    # Blur penalty suppresses high-k gradually in blurred/heterogeneous regions.
    beta = 1.0
    w_blur = np.exp(-beta * alpha * k2_arr)

    w = w_fsc * w_blur
    w_sum = float(np.sum(w))
    if w_sum <= 0.0:
        return None

    k2_mean = float(np.sum(w * k2_arr) / w_sum)

    # Raw width-like measure (cycles→radians via 2π)
    r_raw = math.sqrt(3.0 / k2_mean) / (2.0 * math.pi)

    return {
        "r_info_A": float(r_raw),
        "k2_mean": float(k2_mean),
        "w_sum": float(w_sum),
        "alpha_falloff": float(alpha),
        "fsc_mean": float(np.mean(fsc_arr)),
    }



def locinfo_map(cfg) -> Dict[str, Any]:
    """
    Drop-in compatible replacement for locres_map(cfg), but computes a fast
    local "information radius" map (Å) using band-filtered half-maps and
    SNR-derived weights.

    Compatibility contract:
      - accepts LocResConfig fields used by locres_map (including precomputed_* and rand_*; rand_* ignored)
      - returns dict with keys: maps, coarse_grids, global_fsc, params, stats_per_region
      - params contains sampling_A and radius_A
      - global_fsc contains crossings[AUX_THR] (finite if computed or injected if precomputed provided)
    """
    print("[JANAS locinfo >")

    # ---- Load inputs ----
    if isinstance(cfg.halfmap1, str):
        vol1, hdr1 = IO_utils.read_mrc_data(cfg.halfmap1)
    else:
        vol1 = np.asarray(cfg.halfmap1)
        hdr1 = {"pixel_x": 1.0, "pixel_y": 1.0, "pixel_z": 1.0,
                "origin_x": 0.0, "origin_y": 0.0, "origin_z": 0.0}

    if isinstance(cfg.halfmap2, str):
        vol2, _hdr2 = IO_utils.read_mrc_data(cfg.halfmap2)
    else:
        vol2 = np.asarray(cfg.halfmap2)

    if vol1.shape != vol2.shape:
        raise ValueError(f"Shape mismatch: {vol1.shape} vs {vol2.shape}")

    pix_in = float(hdr1["pixel_x"])
    if (
        not np.isclose(hdr1.get("pixel_x", pix_in), hdr1.get("pixel_y", pix_in))
        or not np.isclose(hdr1.get("pixel_x", pix_in), hdr1.get("pixel_z", pix_in))
    ):
        print("Warning: anisotropic pixel size in header; using pixel_x.")

    # Optional mask
    mask_bool = None
    if getattr(cfg, "mask", None) is not None:
        if isinstance(cfg.mask, str):
            mask_vol, _hdrm = IO_utils.read_mrc_data(cfg.mask)
        else:
            mask_vol = np.asarray(cfg.mask)
        if mask_vol.shape != vol1.shape:
            raise ValueError(f"Mask shape mismatch: {mask_vol.shape} vs maps {vol1.shape}")
        thr = 0.5
        mask_bool = binary_dilation(mask_vol >= thr, iterations=2)
        print(f"Mask loaded: thresholded at {thr}, dilated by 2 voxels; covered voxels: {int(mask_bool.sum()):,}")

    # ---- Reuse-mode compatibility: honour precomputed_sampling/radius if present ----
    pre_samp = getattr(cfg, "precomputed_sampling", None)
    pre_rad  = getattr(cfg, "precomputed_radius", None)
    pre_kcr  = getattr(cfg, "precomputed_first_fsc_crossing", None)

    have_precomputed = (pre_samp is not None) and (pre_rad is not None)

    # ---- Global FSC only if needed (sampling/radius auto or no precomputed) ----
    freqs_g = None
    fsc_g = None
    res_g = {}
    kcross_g = {}

    # resolve sampling/radius inputs (same semantics as locres_map usage)
    sampling_in = getattr(cfg, "sampling", "auto")
    radius_in   = getattr(cfg, "radius", "auto")

    need_global = False
    if not have_precomputed:
        # If either sampling or radius is 'auto', we need global FSC to set them.
        if isinstance(sampling_in, str) and sampling_in.strip().lower() == "auto":
            need_global = True
        if isinstance(radius_in, str) and radius_in.strip().lower() == "auto":
            need_global = True
        # If neither is auto, we can skip global FSC entirely.
    else:
        # reuse mode: skip global FSC computation
        need_global = False

    if need_global:
        freqs_g, fsc_g = compute_fsc(vol1, vol2, pix_in, n_bins=getattr(cfg, "bins", None))
        thresholds_all = [float(getattr(cfg, "local_fsc_threshold", 0.143)), AUX_THR]
        res_g, kcross_g, _, _, _, _ = find_FSC_resolutions_and_stats(freqs_g, fsc_g, thresholds_all)
        if bool(getattr(cfg, "plot", False)):
            crossings_for_plot = {t: k for t, k in kcross_g.items() if (k is not None and np.isfinite(k))}
            plot_fsc(freqs_g, fsc_g, thresholds_all, crossings_for_plot)

    # Sampling selection (compatible with your locamp/locres pattern)
    if have_precomputed:
        sampling_A = float(pre_samp)
        print(f"[locinfo] Reuse sampling from precomputed: {sampling_A:.2f} Å")
    else:
        if isinstance(sampling_in, str) and sampling_in.strip().lower() == "auto":
            # follow your locres behaviour: based on AUX_THR resolution when available
            r_aux = float(res_g.get(AUX_THR, np.nan))
            sampling_A = 0.5 * r_aux if np.isfinite(r_aux) else max(pix_in, 25.0)
            print(f"[locinfo] Auto sampling from global {AUX_THR} resolution: {sampling_A:.2f} Å")
        else:
            sampling_A = float(sampling_in)
            print(f"[locinfo] Sampling (from user): {sampling_A:.2f} Å")

    # Radius selection
    if have_precomputed:
        rad = float(pre_rad)
        R_vox = rad / pix_in
        print(f"[locinfo] Reuse radius from precomputed: {rad:.2f} Å ({R_vox:.1f} vox)")
    else:
        nyq = 0.5 / pix_in
        if isinstance(radius_in, str) and radius_in.strip().lower() == "auto":
            rmax_vox = min(vol1.shape) // 2 - 2
            R_A, R_vox = auto_radius_from_global_fsc(
                freqs_g, fsc_g,
                pixel_A=pix_in,
                cycles_min=float(getattr(cfg, "cycles", 10.0)),
                gamma=float(getattr(cfg, "gamma", 1.8)),
                rmin_vox=6,
                rmax_vox=rmax_vox,
            )
            rad = float(R_A)
            print(f"[locinfo] Auto radius: {rad:.2f} Å ({R_vox:.1f} vox)")
        else:
            rad = float(radius_in)
            R_vox = rad / pix_in
            print(f"[locinfo] Radius (from user): {rad:.2f} Å ({R_vox:.1f} vox)")

    # Ensure crossing field exists for caller caching logic:
    # - If global FSC computed, use actual.
    # - If reuse-mode, inject precomputed value (may be NaN).
    if freqs_g is None or fsc_g is None:
        freqs_g = np.asarray([], dtype=np.float32)
        fsc_g = np.asarray([], dtype=np.float32)
    if AUX_THR not in kcross_g:
        if pre_kcr is not None:
            kcross_g[AUX_THR] = float(pre_kcr)
        else:
            kcross_g[AUX_THR] = float("nan")

    edge_A = (0.3 * rad) if (getattr(cfg, "edgewidth", None) is None) else float(cfg.edgewidth)
    edge_vox = max(int(round(edge_A / pix_in)), 1)
    R_int = int(round(rad / pix_in))

    # ---- Fast band setup ----
    nyq = 0.5 / pix_in
    n_bands = int(getattr(cfg, "info_bands", 6))
    kmin = float(getattr(cfg, "info_kmin", 0.0))
    kmax = float(getattr(cfg, "info_kmax", nyq))

    print(f"[locinfo] Bands={n_bands}, k in [{kmin:.4f},{kmax:.4f}] Å^-1")

    kmag = _fft_freq_grid(vol1.shape, pix_in)
    band_masks, k_centres = _make_radial_bands(kmag, nyquist=nyq, n_bands=n_bands, k_min=kmin, k_max=kmax)
    if len(band_masks) < 2:
        raise RuntimeError("Too few usable bands; adjust info_bands/info_kmin/info_kmax.")
# ---- precompute mean k^2 in each band from the mask support ----
    kmag2 = (kmag.astype(np.float64) ** 2)
    k2_band = []
    for B in band_masks:
        s = float(B.sum())
        if s <= 0:
            k2_band.append(np.nan)
        else:
            k2_band.append(float((kmag2 * B).sum() / s))
    k2_band = np.asarray(k2_band, dtype=np.float64)


    print("[locinfo] Prefiltering bands (one-off)...")
    f1_bands = _prefilter_bands_realspace(vol1.astype(np.float32, copy=False), band_masks)
    f2_bands = _prefilter_bands_realspace(vol2.astype(np.float32, copy=False), band_masks)

    W = _build_soft_sphere_window(R_int, edge_vox)

    # ---- Centres ----
    centres_all = get_region_centres(vol1.shape, pix_in, sampling_A, rad)
    if mask_bool is not None:
        centres = [(x, y, z) for (x, y, z) in centres_all if mask_bool[z, y, x]]
        print(f"[locinfo] Centres within mask: {len(centres)}/{len(centres_all)}")
    else:
        centres = centres_all

    total = len(centres)
    if total == 0:
        print("[locinfo] No region centres to process. Exiting.")
        return {
            "maps": {},
            "coarse_grids": {},
            "global_fsc": {"freqs": freqs_g, "fsc": fsc_g, "resolutions": res_g, "crossings": kcross_g},
            "params": {
                "pixel_size_A": pix_in,
                "sampling_A": sampling_A,
                "radius_A": rad,
                "edge_A": edge_A,
                "centres_count": 0,
                "interp_mode": getattr(cfg, "interp", "cubic"),
                "resample": bool(getattr(cfg, "resample", False)),
            },
            "stats_per_region": {},
        }

    print(f"[locinfo] Total regions: {total}")

    # ---- Parallel (fork) ----
    cpu = int(getattr(cfg, "cpu", 1))
    cpu_max = mp.cpu_count()
    cpu = min(max(cpu, 1), max(1, min(cpu_max, total)))

    try:
        mp.set_start_method("fork")
    except RuntimeError:
        pass

    # Set globals for workers
    global _LOCINFO_F1B, _LOCINFO_F2B, _LOCINFO_KC, _LOCINFO_K2B, _LOCINFO_W, _LOCINFO_R
    _LOCINFO_F1B = f1_bands
    _LOCINFO_F2B = f2_bands
    _LOCINFO_KC = k_centres
    _LOCINFO_K2B = k2_band
    _LOCINFO_W = W
    _LOCINFO_R = R_int

    block = int(math.floor(total / cpu)) if cpu > 0 else total
    idx_blocks = []
    for i in range(cpu):
        b0 = i * block
        b1 = (i + 1) * block if i < cpu - 1 else total
        idx_blocks.append((b0, b1))

    def _worker(b0: int, b1: int, procnum: int, ret_dict, progress):
        results = []
        stats_per = {}
        FLUSH_EVERY = 64
        acc = 0
        for idx in range(b0, b1):
            cx, cy, cz = centres[idx]
            st = _compute_local_info_radius(cx, cy, cz)
            if st is None:
                continue
            # store centre + r_info + fsc_mean
            results.append((cx, cy, cz, st["r_info_A"], st["fsc_mean"]))
            stats_per[(cx, cy, cz)] = st

            acc += 1
            if (acc & (FLUSH_EVERY - 1)) == 0:
                with progress.get_lock():
                    progress.value += acc
                acc = 0

        if acc:
            with progress.get_lock():
                progress.value += acc
        ret_dict[procnum] = (b0, b1, results, stats_per)

    progress = mp.Value("i", 0)
    manager = mp.Manager()
    ret_dict = manager.dict()

    procs = []
    for pi, (b0, b1) in enumerate(idx_blocks):
        procs.append(mp.Process(target=_worker, args=(b0, b1, pi, ret_dict, progress)))

    start = time.time()
    for p in procs:
        p.start()

    last = -1
    while any(p.is_alive() for p in procs):
        done = progress.value
        if done != last:
            elapsed = time.time() - start
            pct = 100.0 * done / max(1, total)
            rate = done / elapsed if elapsed > 0 else 0.0
            remain = (total - done) / rate if rate > 0 else 0.0
            print(
                f"\r[locinfo] Processed {done:,}/{total:,} ({pct:5.1f}%) — "
                f"{rate:,.1f} regions/s — elapsed {time.strftime('%H:%M:%S', time.gmtime(elapsed))} — "
                f"ETA {time.strftime('%H:%M:%S', time.gmtime(remain))}",
                end="",
            )
            last = done
        time.sleep(0.2)
    for p in procs:
        p.join()
    print("")

    results: List[Tuple[int, int, int, float, float]] = []
    stats_per_region: Dict[Tuple[int, int, int], dict] = {}
    for pi in range(len(idx_blocks)):
        _b0, _b1, res_block, stats_block = ret_dict[pi]
        results.extend(res_block)
        stats_per_region.update(stats_block)

    if not results:
        print("[locinfo] No valid regions produced r_info.")
        return {
            "maps": {},
            "coarse_grids": {},
            "global_fsc": {"freqs": freqs_g, "fsc": fsc_g, "resolutions": res_g, "crossings": kcross_g},
            "params": {
                "pixel_size_A": pix_in,
                "sampling_A": sampling_A,
                "radius_A": rad,
                "edge_A": edge_A,
                "centres_count": total,
                "interp_mode": getattr(cfg, "interp", "cubic"),
                "resample": bool(getattr(cfg, "resample", False)),
            },
            "stats_per_region": stats_per_region,
        }

    # ---- Coarse grids ----
    xs = sorted({r[0] for r in results})
    ys = sorted({r[1] for r in results})
    zs = sorted({r[2] for r in results})
    n_x, n_y, n_z = len(xs), len(ys), len(zs)
    xi = {x: i for i, x in enumerate(xs)}
    yi = {y: j for j, y in enumerate(ys)}
    zi = {z: k for k, z in enumerate(zs)}

    loc_coarse = np.zeros((n_z, n_y, n_x), dtype=np.float32)   # r_info (Å)
    fscm_coarse = np.zeros_like(loc_coarse)                    # mean band FSC (diagnostic)
    for cx, cy, cz, rinfo, fscm in results:
        i, j, k = xi[cx], yi[cy], zi[cz]
        loc_coarse[k, j, i] = float(rinfo)
        fscm_coarse[k, j, i] = float(fscm)

    # ---- Upsample / embed (same pattern) ----
    origin_in = np.array(
        [float(hdr1.get("origin_x", 0.0)),
         float(hdr1.get("origin_y", 0.0)),
         float(hdr1.get("origin_z", 0.0))],
        dtype=np.float64,
    )

    if bool(getattr(cfg, "resample", False)):
        out_pix = float(sampling_A)
        x0, y0, z0 = xs[0], ys[0], zs[0]
        centre_in = origin_in + np.array(
            [(x0 + 0.5) * pix_in, (y0 + 0.5) * pix_in, (z0 + 0.5) * pix_in],
            dtype=np.float64,
        )
        origin_out = centre_in + out_pix

        loc_grid = loc_coarse
        write_specs = dict(pixel_size=out_pix, origin_angs=tuple(origin_out))
    else:
        step_vox = max(1, int(round(float(sampling_A) / pix_in)))
        z_start, y_start, x_start = zs[0], ys[0], xs[0]
        z_len = (len(zs) - 1) * step_vox + 1
        y_len = (len(ys) - 1) * step_vox + 1
        x_len = (len(xs) - 1) * step_vox + 1

        order = 1 if getattr(cfg, "interp", "cubic") == "linear" else 3

        if mask_bool is not None:
            up_loc = _upsample_with_nearest_fill(loc_coarse, step_vox, (z_len, y_len, x_len), order=order)
        else:
            up_loc = _upsample_interp_ignore_zeros(loc_coarse, step_vox, (z_len, y_len, x_len), order=order)

        loc_full = embed_into_full_grid(up_loc, (z_start, y_start, x_start), vol1.shape)
        if mask_bool is not None:
            loc_full *= mask_bool.astype(np.float32)

        nz_full = loc_full != 0
        holes = (mask_bool & (~nz_full)) if (mask_bool is not None) else (~nz_full)
        if np.any(holes):
            _, idx = distance_transform_edt(~nz_full, return_distances=True, return_indices=True)
            filled = loc_full[idx[0], idx[1], idx[2]]
            loc_full[holes] = filled[holes]

            if not np.any(nz_full):
                vals = loc_coarse[np.isfinite(loc_coarse) & (loc_coarse > 0)]
                plateau = float(np.median(vals)) if vals.size else float(rad)
                inside = (mask_bool & (loc_full == 0)) if (mask_bool is not None) else (loc_full == 0)
                loc_full[inside] = plateau

        loc_grid = loc_full.astype(np.float32, copy=False)
        write_specs = dict(pixel_size=pix_in, origin_angs=tuple(origin_in))

    # ---- Output files (naming consistent with locres style) ----
    if getattr(cfg, "output_basename", None):
        base = cfg.output_basename
        fn = f"{base}_locinfo.mrc" if not bool(getattr(cfg, "resample", False)) else f"{base}_resampledLocinfo.mrc"
        IO_utils.write_mrc(fn, loc_grid, **write_specs)
        print("[locinfo] Wrote", fn)

        fn_params = f"{base}_locinfo_params.txt"
        with open(fn_params, "w") as pf:
            pf.write(f"pixel_size_A\t{pix_in:.6f}\n")
            pf.write(f"sampling_A\t{sampling_A:.6f}\n")
            pf.write(f"radius_A\t{rad:.6f}\n")
            pf.write(f"edge_A\t{edge_A:.6f}\n")
            pf.write(f"info_bands\t{len(band_masks)}\n")
            pf.write(f"info_kmin_Ainv\t{kmin:.6f}\n")
            pf.write(f"info_kmax_Ainv\t{kmax:.6f}\n")
        print("[locinfo] Wrote", fn_params)

    out = {
        "maps": {
            "locinfo": loc_grid,  # Å
        },
        "coarse_grids": {
            "locinfo": loc_coarse,
            "fsc_mean": fscm_coarse,
        },
        "global_fsc": {
            "freqs": freqs_g,
            "fsc": fsc_g,
            "resolutions": res_g,
            "crossings": kcross_g,
        },
        "params": {
            "pixel_size_A": pix_in,
            "sampling_A": float(sampling_A),
            "radius_A": float(rad),
            "edge_A": float(edge_A),
            "centres_count": int(total),
            "interp_mode": getattr(cfg, "interp", "cubic"),
            "resample": bool(getattr(cfg, "resample", False)),
            "info_bands": int(len(band_masks)),
            "info_kmin_Ainv": float(kmin),
            "info_kmax_Ainv": float(kmax),
            "band_centres_Ainv": k_centres,
        },
        "stats_per_region": stats_per_region,
    }
    print("]")

    return out

# ----------------------- Generic FSC + helpers -----------------------
# ----------------------- Generic FSC + helpers -----------------------
def compute_fsc(vol1: np.ndarray,
                vol2: np.ndarray,
                pixel_spacing: float,
                n_bins: Optional[int] = None,
                eps: float = 1e-8) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute global FSC between two real-space half-maps by FFT'ing them.
    Returns (freqs [1/Å], fsc).
    """
    from numpy.fft import fftn, fftfreq
    if vol1.shape != vol2.shape:
        raise ValueError("Half-maps must have identical shapes.")
    if not np.isfinite(vol1).all() or not np.isfinite(vol2).all():
        raise ValueError("NaNs/Infs in input volumes.")

    F1 = fftn(vol1); F2 = fftn(vol2)
    nz, ny, nx = vol1.shape
    fx = fftfreq(nx, d=pixel_spacing)[None, None, :]
    fy = fftfreq(ny, d=pixel_spacing)[None, :, None]
    fz = fftfreq(nz, d=pixel_spacing)[:, None, None]
    radii = np.sqrt(fx * fx + fy * fy + fz * fz)

    max_freq = float(np.max(radii))
    if n_bins is None:
        n_bins = int(min(vol1.shape) // 2)
        n_bins = max(n_bins, 16)

    bins  = np.linspace(0.0, max_freq, n_bins + 1, dtype=np.float64)
    freqs = 0.5 * (bins[:-1] + bins[1:])

    idx = np.digitize(radii.ravel(), bins) - 1
    S1r = F1.ravel(); S2r = F2.ravel()
    fsc = np.zeros(n_bins, dtype=np.float32)
    for i in range(n_bins):
        mask = (idx == i)
        if not np.any(mask): 
            continue
        num = np.sum(S1r[mask] * np.conj(S2r[mask]))
        den = np.sqrt(np.sum(np.abs(S1r[mask])**2) * np.sum(np.abs(S2r[mask])**2)) + eps
        fsc[i] = float(np.real(num / den))
    return freqs.astype(np.float32), fsc


def compute_ssnr_from_halfmaps(vol1: np.ndarray,
                               vol2: np.ndarray,
                               pixel_spacing: float,
                               n_bins: Optional[int] = None,
                               eps: float = 1e-8) -> Tuple[np.ndarray, np.ndarray]:
    """
    Estimate per-shell SSNR from the gold-standard FSC between two half-maps:
        SSNR(k) = FSC(k) / (1 - FSC(k)), with clipping for stability.
    Returns (freqs [1/Å], ssnr).
    """
    freqs, fsc = compute_fsc(vol1, vol2, pixel_spacing, n_bins=n_bins, eps=eps)
    fsc = np.clip(fsc.astype(np.float32), 0.0, 0.99)
    denom = np.maximum(1.0 - fsc, 1e-6)
    ssnr  = fsc / denom
    return freqs.astype(np.float32), ssnr.astype(np.float32)

def compute_ssnr_from_paths(path1: str,
                            path2: str,
                            n_bins: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Load two half-maps from disk and compute (freqs, SSNR, angpix).
    """
    import janas.janas_core as janas_core
    s1 = janas_core.sizeMRC(path1); s2 = janas_core.sizeMRC(path2)
    if s1 != s2:
        raise ValueError(f"Half-maps must share shape: {s1} vs {s2}")
    v1 = np.reshape(np.array(janas_core.ReadMRC(path1)), s1)
    v2 = np.reshape(np.array(janas_core.ReadMRC(path2)), s2)
    angpix = float(janas_core.spacingMRC(path1))
    freqs, ssnr = compute_ssnr_from_halfmaps(v1, v2, angpix, n_bins=n_bins)
    return freqs, ssnr, angpix


def compute_fsc_from_kspace(F1: np.ndarray,
                            F2: np.ndarray,
                            pixel_spacing: float,
                            n_bins: Optional[int] = None,
                            eps: float = 1e-8) -> Tuple[np.ndarray, np.ndarray]:
    """
    FSC computed directly in Fourier space. F1 and F2 are complex64 k-space volumes
    representing the reconstructed half-maps BEFORE iFFT (i.e., after dividing by weights).
    Returns (freqs [1/Å], fsc).
    """
    from numpy.fft import fftfreq
    if F1.shape != F2.shape:
        raise ValueError("Half-maps must have identical shapes.")
    if not np.isfinite(F1.real).all() or not np.isfinite(F1.imag).all():
        raise ValueError("NaNs/Infs in input k-space volumes.")
    if not np.isfinite(F2.real).all() or not np.isfinite(F2.imag).all():
        raise ValueError("NaNs/Infs in input k-space volumes.")

    nz, ny, nx = F1.shape
    fx = fftfreq(nx, d=pixel_spacing)[None, None, :]
    fy = fftfreq(ny, d=pixel_spacing)[None, :, None]
    fz = fftfreq(nz, d=pixel_spacing)[:, None, None]
    radii = np.sqrt(fx * fx + fy * fy + fz * fz)

    max_freq = float(np.max(radii))
    if n_bins is None:
        n_bins = int(min(nx, ny, nz) // 2)
        n_bins = max(n_bins, 16)

    bins  = np.linspace(0.0, max_freq, n_bins + 1, dtype=np.float64)
    freqs = 0.5 * (bins[:-1] + bins[1:])

    idx = np.digitize(radii.ravel(), bins) - 1
    S1r = F1.ravel()
    S2r = F2.ravel()
    fsc = np.zeros(n_bins, dtype=np.float32)
    for i in range(n_bins):
        mask = (idx == i)
        if not np.any(mask): 
            continue
        num = np.sum(S1r[mask] * np.conj(S2r[mask]))
        den = np.sqrt(np.sum(np.abs(S1r[mask])**2) * np.sum(np.abs(S2r[mask])**2)) + eps
        fsc[i] = float(np.real(num / den))
    return freqs.astype(np.float32), fsc


def first_FSC_crossing(freqs: np.ndarray, y: np.ndarray, thr: float) -> Optional[float]:
    """
    Return the first frequency (1/Å) where FSC crosses from >=thr to <thr.
    Linear-interpolates between the bracketing samples. Returns None if no crossing.
    """
    if len(freqs) < 2:
        return None
    above = y >= thr
    idx = np.where(above[:-1] & (~above[1:]))[0]
    if idx.size == 0:
        return None
    i = int(idx[0])
    f1, f2 = float(freqs[i]), float(freqs[i+1])
    y1, y2 = float(y[i]), float(y[i+1])
    frac = (y1 - thr) / (y1 - y2 + 1e-12)
    return f1 + frac * (f2 - f1)


def find_FSC_resolutions_and_stats(freqs: np.ndarray,
                                   fsc: np.ndarray,
                                   thresholds: List[float],
                                   interp_factor: int = 10,
                                   kind: str = 'linear'):
    """
    (Optional utility) Interpolate FSC and report resolution (Å) for each threshold.
    Returns (res_dict, cross_freq_dict, freq_grid, fsc_interp, freqs, fsc).
    """
    from scipy.interpolate import interp1d
    N = len(freqs)
    if N < 2:
        return {}, {}, np.nan, np.nan, freqs, fsc

    fd = np.linspace(freqs[0], freqs[-1], N * interp_factor)
    interp = interp1d(freqs, fsc, kind=kind, bounds_error=False, fill_value=np.nan)
    yd = interp(fd)

    res = {}
    cross = {}
    for thr in thresholds:
        fc = first_FSC_crossing(fd, yd, float(thr))
        cross[thr] = fc
        res[thr] = (1.0 / fc) if (fc and fc > 0.0) else np.nan
    return res, cross, fd, yd, freqs, fsc



def plot_fsc(freqs, fsc_vals, thresholds, f_crossings=None):
    plt.figure()
    plt.plot(freqs, fsc_vals, lw=1, label='FSC')
    for t in thresholds:
        plt.hlines(t, float(freqs.min()), float(freqs.max()), linestyles='--', label=f'{t:.3f}')
    if f_crossings:
        for t, f in f_crossings.items():
            plt.vlines(f, 0, t, linestyles=':', label=f'cross@{t:.3f}: {f:.4f}')
    plt.xlabel('Spatial frequency (Å⁻¹)')
    plt.ylabel('FSC')
    plt.legend(loc='lower left')
    plt.grid(True)
    plt.tight_layout()
    plt.show()


# ----------------------- Windowing + randomisation -----------------------
def auto_radius_from_global_fsc(freqs, fsc, pixel_A,
                                cycles_min=10.0, gamma=1.8,
                                rmin_vox=6, rmax_vox=None):
    nyq = 0.5 / pixel_A
    k0143 = first_FSC_crossing(freqs, fsc, thr=AUX_THR)
    if not (k0143 and np.isfinite(k0143) and k0143 > 0):
        k_tgt = 0.6 * nyq
    else:
        k_tgt = min(float(gamma) * float(k0143), 0.9 * nyq)

    # Apply smooth frequency floor (plateau ~55–60 Å for cycles_min=10)
    k0 = 0.09  # Å^-1  → R_A_max ≈ cycles_min / (2*k0) = 55.6 Å
    k_eff = math.hypot(float(k_tgt), k0)  # sqrt(k_tgt^2 + k0^2)
    R_A = float(cycles_min) / (2.0 * k_eff)

    R_vox = R_A / pixel_A
    R_vox = max(R_vox, float(rmin_vox))
    if rmax_vox is not None:
        R_vox = min(R_vox, float(rmax_vox))
    return R_A, R_vox


def spherical_raised_cosine_window(shape, cx, cy, cz, R_vox, edge_vox):
    nz, ny, nx = shape
    z = np.arange(nz)[:, None, None]
    y = np.arange(ny)[None, :, None]
    x = np.arange(nx)[None, None, :]
    r = np.sqrt((x - cx)**2 + (y - cy)**2 + (z - cz)**2).astype(np.float32)
    R1 = max(float(R_vox) - float(edge_vox), 0.0)
    w = np.zeros_like(r, dtype=np.float32)
    core = (r <= R1); rim = (r > R1) & (r <= float(R_vox))
    w[core] = 1.0
    if np.any(rim):
        t = (r[rim] - R1) / max(float(edge_vox), 1e-6)
        w[rim] = 0.5 * (1.0 + np.cos(np.pi * t))
    return w


def phase_randomize_beyond(vol, k_cut, pixel_A, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    F = fftn(vol)
    Rnd = rng.standard_normal(vol.shape, dtype=np.float32)
    FR = fftn(Rnd)

    nz, ny, nx = vol.shape
    fx = fftfreq(nx, d=pixel_A)[None, None, :]
    fy = fftfreq(ny, d=pixel_A)[None, :, None]
    fz = fftfreq(nz, d=pixel_A)[:, None, None]
    rad = np.sqrt(fx * fx + fy * fy + fz * fz).astype(np.float32)

    mask = (rad > float(k_cut))
    F[mask] = np.abs(F[mask]) * np.exp(1j * np.angle(FR[mask]))
    out = np.fft.ifftn(F).real.astype(np.float32)
    return out


def resolve_kcut(option,
                 pixel_A,
                 k0143,
                 nyq,
                 frac_nyq=None,
                 auto_mult=1.25,
                 fallback_frac=0.75,
                 min_frac=0.60,
                 max_frac=0.85):
    """
    Parse randomisation option into a k-cut (1/Å).

    - option:
        'none'
        'auto'
        'auto:<mult>'  e.g., 'auto:1.35'
        numeric Å value (interpreted as resolution)
        numeric k value (if you pass a float and treat it as k explicitly)
    - frac_nyq:
        if provided, overrides option and returns frac_nyq * nyq.

    Notes:
      This returns a frequency in 1/Å.
    """

    # Explicit nyquist fraction override
    if (frac_nyq is not None) and (0.0 < float(frac_nyq) < 1.0):
        return float(frac_nyq) * float(nyq)

    # String options
    if isinstance(option, str):
        opt = option.strip().lower()

        if opt == 'none':
            return None

        if opt.startswith('auto'):
            # parse 'auto:1.35' if provided
            mult = auto_mult
            if ':' in opt:
                try:
                    mult = float(opt.split(':', 1)[1])
                except Exception:
                    mult = auto_mult

            if (k0143 is not None) and (k0143 > 0):
                base = mult * float(k0143)
            else:
                base = float(fallback_frac) * float(nyq)

            kcut = np.clip(
                base,
                float(min_frac) * float(nyq),
                float(max_frac) * float(nyq)
            )
            return float(kcut)

        # Numeric string interpreted as resolution in Å
        try:
            res_A = float(opt)
            if res_A > 0:
                return 1.0 / res_A
        except Exception:
            pass

        raise ValueError(f"Cannot parse randomisation value: {option}")

    # Non-string numeric
    # Interpret as resolution Å if it looks like Å-scale (> 1),
    # otherwise assume it's already k (1/Å).
    val = float(option)
    if val <= 0:
        raise ValueError(f"Cannot parse randomisation value: {option}")
    if val > 1.0:
        return 1.0 / val
    return val





# ----------------------- Region logic -----------------------
def get_region_centres(vol_shape: Tuple[int, int, int],
                       pixel_size_A: float,
                       sampling_A: float,
                       region_radius_A: Optional[float] = None) -> List[Tuple[int, int, int]]:
    nz, ny, nx = vol_shape
    if region_radius_A is None:
        region_radius_A = 2.0 * sampling_A
    rvox = int(round(region_radius_A / pixel_size_A))
    step = max(1, int(round(sampling_A / pixel_size_A)))

    z_range = range(rvox, nz - rvox + 1, step)
    y_range = range(rvox, ny - rvox + 1, step)
    x_range = range(rvox, nx - rvox + 1, step)

    coords = []
    for z in z_range:
        for y in y_range:
            for x in x_range:
                coords.append((x, y, z))
    return coords

def compute_local_fsc_region(vol1: np.ndarray, vol2: np.ndarray,
                             cx: int, cy: int, cz: int,
                             pixel_A: float,
                             radius_A: float,
                             thresholds: List[float],
                             n_bins: Optional[int] = None,
                             eps: float = 1e-8,
                             interp_factor: int = 10,
                             interp_kind: str = 'cubic',
                             edge_vox: int = 2,
                             kcut_local: Optional[float] = None,
                             rng: Optional[np.random.Generator] = None,
                             mask_vol: Optional[np.ndarray] = None,
                             mask_thr: float = 0.90,
                             min_coverage: float = 0.90):
    """
    Local FSC in a spherical raised-cosine window.

    If mask_vol is provided (float mask 0..1), we:
      - hard-threshold the mask patch at mask_thr to define a *support* (not weighting),
      - require that the support covers at least min_coverage of the spherical window,
      - compute FSC only on that supported region.
    Centres in the soft-edge band will return None and should be inpainted later.
    """
    R_vox = int(round(radius_A / pixel_A))
    z0, z1 = cz - R_vox, cz + R_vox
    y0, y1 = cy - R_vox, cy + R_vox
    x0, x1 = cx - R_vox, cx + R_vox

    nz, ny, nx = vol1.shape
    if z0 < 0 or y0 < 0 or x0 < 0 or z1 > nz or y1 > ny or x1 > nx:
        return None

    H1 = vol1[z0:z1, y0:y1, x0:x1].copy()
    H2 = vol2[z0:z1, y0:y1, x0:x1].copy()
    if H1.size == 0 or H2.size == 0 or min(H1.shape) < 4:
        return None

    lz, ly, lx = H1.shape
    cz_loc, cy_loc, cx_loc = (lz // 2, ly // 2, lx // 2)
    W = spherical_raised_cosine_window(
        (lz, ly, lx),
        cx_loc, cy_loc, cz_loc,
        R_vox,
        edge_vox=edge_vox
    ).astype(np.float32, copy=False)

    # ---- NEW: restrict computation to mask-supported voxels (avoid soft-edge artefacts) ----
    if mask_vol is not None:
        Mp = mask_vol[z0:z1, y0:y1, x0:x1].astype(np.float32, copy=False)
        S = (Mp >= float(mask_thr)).astype(np.float32, copy=False)  # hard support
        Weff = W * S
        cov = float(Weff.sum() / (W.sum() + 1e-12))
        if cov < float(min_coverage):
            return None  # do not compute in the soft edge band

        # mild renormalisation to keep energies away from eps-dominated regimes
        s = float(Weff.sum())
        if s > 0:
            Weff *= float(W.sum() / (s + 1e-12))
        H1 *= Weff
        H2 *= Weff
    else:
        H1 *= W
        H2 *= W

    if (kcut_local is not None) and np.isfinite(kcut_local) and (kcut_local > 0):
        if rng is None:
            rng = np.random.default_rng()
        H1 = phase_randomize_beyond(H1, kcut_local, pixel_A, rng=rng)
        H2 = phase_randomize_beyond(H2, kcut_local, pixel_A, rng=rng)

    if n_bins is None:
        n_bins = max(16, R_vox)

    freqs, fsc_vals = compute_fsc(H1, H2, pixel_A, n_bins, eps)
    N = len(freqs)
    if N < 2:
        return None

    fd = np.linspace(freqs[0], freqs[-1], N * interp_factor)
    interp = interp1d(freqs, fsc_vals, kind=interp_kind, bounds_error=False, fill_value="extrapolate")
    fsd = interp(fd)

    fcross = {}
    resol = {}
    for t in thresholds:
        above = fsd >= t
        cr = np.where(above[:-1] & (~above[1:]))[0]
        if cr.size > 0:
            i = int(cr[0])
            f1, f2 = fd[i], fd[i + 1]
            y1, y2 = fsd[i], fsd[i + 1]
            frac = (y1 - t) / (y1 - y2 + 1e-12)
            fc = f1 + frac * (f2 - f1)
        else:
            mi = int(np.nanargmin(np.abs(fsd - t)))
            fc = fd[mi]
        fcross[t] = float(fc)
        resol[t] = (1.0 / float(fc)) if fc > 0 else np.inf

    t_hi, t_lo = max(thresholds), min(thresholds)
    fh, fl = fcross[t_hi], fcross[t_lo]
    if fh > fl:
        fh, fl = fl, fh
    mask = (fd >= fh) & (fd <= fl)
    seg = fsd[mask]
    mean_between = float(np.nanmean(seg)) if seg.size else np.nan
    var_between = float(np.nanvar(seg)) if seg.size else np.nan

    mean_val = float(((H1 + H2) / 2.0).mean())

    return {
        'f_crossings': fcross,
        'resolutions': resol,
        'mean_between': mean_between,
        'var_between': var_between,
        'freqs': freqs,
        'fsc_vals': fsc_vals,
        'freqs_dense': fd,
        'fsc_dense': fsd,
        'mean_val': mean_val,
    }




# ----------------------- Coarse-grid upsampling helpers -----------------------
def _upsample_interp_ignore_zeros(coarse: np.ndarray,
                                  step_vox: int,
                                  out_len_zyx: Tuple[int, int, int],
                                  order: int = 3) -> np.ndarray:
    if step_vox <= 1:
        return coarse.astype(np.float32, copy=False)

    mask = (coarse != 0).astype(np.float32)
    vals = coarse.astype(np.float32) * mask

    zoom_factors = (float(step_vox), float(step_vox), float(step_vox))
    up_vals = ndi_zoom(vals, zoom=zoom_factors, order=order, prefilter=True, mode='nearest')
    up_mask = ndi_zoom(mask, zoom=zoom_factors, order=order, prefilter=True, mode='nearest')

    eps = 1e-6
    up = np.divide(up_vals, up_mask, out=np.zeros_like(up_vals, dtype=np.float32), where=(up_mask > eps))

    z_len, y_len, x_len = out_len_zyx
    up = up[:z_len, :y_len, :x_len]
    return up.astype(np.float32, copy=False)


# ----------------------- Coarse-grid upsampling helpers -----------------------
def _upsample_interp(coarse: np.ndarray,
                                  step_vox: int,
                                  out_len_zyx: Tuple[int, int, int],
                                  order: int = 3) -> np.ndarray:
    if step_vox <= 1:
        return coarse.astype(np.float32, copy=False)

    vals = coarse.astype(np.float32, copy=False)

    zoom_factors = (float(step_vox), float(step_vox), float(step_vox))
    up = ndi_zoom(vals, zoom=zoom_factors, order=order,
                  prefilter=(order > 1), mode='nearest')

    z_len, y_len, x_len = out_len_zyx
    up = up[:z_len, :y_len, :x_len]
    return up.astype(np.float32, copy=False)


def _upsample_with_nearest_fill(coarse: np.ndarray,
                                step_vox: int,
                                out_len_zyx: Tuple[int, int, int],
                                order: int = 3) -> np.ndarray:
    if step_vox <= 1:
        return coarse.astype(np.float32, copy=False)

    coarse = coarse.astype(np.float32, copy=False)
    nonzero = (coarse != 0)

    if np.any(nonzero):
        _, indices = distance_transform_edt(~nonzero, return_distances=True, return_indices=True)
        filled = coarse[indices[0], indices[1], indices[2]]
    else:
        filled = coarse

    zoom_factors = (float(step_vox), float(step_vox), float(step_vox))
    up = ndi_zoom(filled, zoom=zoom_factors, order=order, prefilter=True, mode='nearest')

    z_len, y_len, x_len = out_len_zyx
    up = up[:z_len, :y_len, :x_len]
    return up.astype(np.float32, copy=False)


def embed_into_full_grid(block: np.ndarray,
                         start_zyx: Tuple[int, int, int],
                         full_shape: Tuple[int, int, int]) -> np.ndarray:
    full = np.zeros(full_shape, dtype=np.float32)
    z0, y0, x0 = start_zyx
    bz, by, bx = block.shape
    full[z0:z0 + bz, y0:y0 + by, x0:x0 + bx] = block
    return full


def _bbox_from_mask_bool(mask_bool: np.ndarray) -> Optional[Tuple[Tuple[int, int], Tuple[int, int], Tuple[int, int]]]:
    idx = np.where(mask_bool)
    if idx[0].size == 0:
        return None
    return (
        (int(idx[0].min()), int(idx[0].max()) + 1),
        (int(idx[1].min()), int(idx[1].max()) + 1),
        (int(idx[2].min()), int(idx[2].max()) + 1),
    )


def _bbox_mask(shape: Tuple[int, int, int],
               bbox: Optional[Tuple[Tuple[int, int], Tuple[int, int], Tuple[int, int]]]) -> Optional[np.ndarray]:
    if bbox is None:
        return None
    out = np.zeros(shape, dtype=bool)
    (z0, z1), (y0, y1), (x0, x1) = bbox
    out[z0:z1, y0:y1, x0:x1] = True
    return out


def _bbox_contains_point(bbox: Tuple[Tuple[int, int], Tuple[int, int], Tuple[int, int]],
                         x: int, y: int, z: int) -> bool:
    (z0, z1), (y0, y1), (x0, x1) = bbox
    return (z0 <= z < z1) and (y0 <= y < y1) and (x0 <= x < x1)



def _masked_stats_for_locres(locres_vol: np.ndarray, mask_bool: Optional[np.ndarray]) -> Tuple[float,float,float,float,float]:
    """
    Returns (best, q25, mean, q75, worst) over masked voxels.
    'best' is numerically smallest value (max resolution in colloquial cryo-EM).
    """
    if mask_bool is not None:
        vals = locres_vol[mask_bool]
    else:
        vals = locres_vol.ravel()
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return (np.nan, np.nan, np.nan, np.nan, np.nan)
    best = float(np.nanmin(vals))
    q25  = float(np.nanpercentile(vals, 25))
    mean = float(np.nanmean(vals))
    q75  = float(np.nanpercentile(vals, 75))
    worst= float(np.nanmax(vals))
    return (best, q25, mean, q75, worst)


# ----------------------- Data class -----------------------
@dataclass
class LocResConfig:
    # Inputs: either paths or arrays
    halfmap1: Union[str, np.ndarray]
    halfmap2: Union[str, np.ndarray]
    mask: Optional[Union[str, np.ndarray]] = None
    output_mask: Optional[Union[str, np.ndarray]] = None

    # Parameters
    local_fsc_threshold: float = 0.143
    sampling: Union[str, float] = "auto"  # Å or 'auto'
    radius: Union[str, float] = "auto"     # Å or 'auto'
    edgewidth: Optional[float] = None       # Å, default 0.3*radius
    cycles: float = 10.0
    gamma: float = 1.8
    bins: Optional[int] = None
    rand_global: str = "auto"               # 'auto' | 'none' | Å
    rand_local: str = "auto"                # 'auto' | 'none' | Å (default auto -> global)
    rand_frac_nyq: Optional[float] = None   # overrides global if set (0..1)
    mask_support_threshold: float = 0.90
    min_mask_coverage: float = 0.90
    mask_dilation_vox: int = 6
    compute_only_mask: bool = False
    outside_mask_factor: int = 8
    fast: bool = True
    outside_box_fill_A: float = 100.0
    local_interp_factor: int = 10

    # Execution
    cpu: int = 1

    # Output behaviour
    resample: bool = False
    interp: str = 'cubic'  # 'linear' | 'cubic'
    output_basename: Optional[str] = None
    plot: bool = False
    precomputed_sampling: Optional[float] = None
    precomputed_radius: Optional[float] = None
    precomputed_first_fsc_crossing: Optional[float] = None
    reference_fsc_freqs: Optional[np.ndarray] = None
    reference_fsc_vals: Optional[np.ndarray] = None


# ----------------------- Core API -----------------------
def locres_map(cfg: LocResConfig) -> Dict[str, Any]:
    """Run local FSC resolution mapping.

    Returns a dict with keys:
      - 'maps': {'locres'} as np.ndarray (full-size if not resampled; coarse if resampled)
      - 'coarse_grids': {'locmap', 'avmap', 'invvarmap'} on sampling grid
      - 'global_fsc': {'freqs', 'fsc', 'resolutions', 'crossings'}
      - 'params': dict with pixel_size_A, sampling_A, radius_A, edge_A, centres_count
      - 'stats_per_region': dict keyed by (x,y,z) with per-region FSC details

    If cfg.output_basename is set, writes MRC/TXT files to disk mirroring
    the previous reconstructor behaviour.
    """
    print("[JANAS locres >")
    # ---- Load inputs ----
    if isinstance(cfg.halfmap1, str):
        vol1, hdr1 = IO_utils.read_mrc_data(cfg.halfmap1)
    else:
        vol1 = np.asarray(cfg.halfmap1)
        hdr1 = {
            "pixel_x": 1.0, "pixel_y": 1.0, "pixel_z": 1.0,
            "origin_x": 0.0, "origin_y": 0.0, "origin_z": 0.0,
        }
    if isinstance(cfg.halfmap2, str):
        vol2, _hdr2 = IO_utils.read_mrc_data(cfg.halfmap2)
    else:
        vol2 = np.asarray(cfg.halfmap2)

    if vol1.shape != vol2.shape:
        raise ValueError(f"Shape mismatch: {vol1.shape} vs {vol2.shape}")

    # ---- Optional: apply reference radial amplitude profile to both half-maps ----
    ref_freqs = getattr(cfg, "reference_fsc_freqs", None)
    ref_vals  = getattr(cfg, "reference_fsc_vals", None)

    if False and ref_freqs is not None and ref_vals is not None:
        angpix = float(hdr1.get("pixel_x", 1.0))
        shape = vol1.shape
        Fx, Fy, Fz = janas_mapProcess._fourier_freq_grid(shape, angpix)
        s_grid = np.sqrt(Fx * Fx + Fy * Fy + Fz * Fz)

        FT1 = np.fft.fftn(vol1.astype(np.float64, copy=False))
        FT2 = np.fft.fftn(vol2.astype(np.float64, copy=False))

        janas_mapProcess._apply_fsc_weighting_inplace(FT1, s_grid, ref_freqs, ref_vals)
        janas_mapProcess._apply_fsc_weighting_inplace(FT2, s_grid, ref_freqs, ref_vals)

        vol1 = np.real(np.fft.ifftn(FT1)).astype(np.float32, copy=False)
        vol2 = np.real(np.fft.ifftn(FT2)).astype(np.float32, copy=False)

    pix_in = float(hdr1["pixel_x"])
    if (not np.isclose(hdr1.get("pixel_x", pix_in), hdr1.get("pixel_y", pix_in)) or
        not np.isclose(hdr1.get("pixel_x", pix_in), hdr1.get("pixel_z", pix_in))):
        print("    Warning: anisotropic pixel size in header; using pixel_x for FSC.")


    compute_only_mask = bool(getattr(cfg, "compute_only_mask", False))

    # Optional masks
    mask_bool = None
    mask_bool_dilated = None
    mask_float = None
    mask_float_dilated = None
    output_mask_bool = None
    if cfg.mask is not None:
        if isinstance(cfg.mask, str):
            mask_vol, hdr_m = IO_utils.read_mrc_data(cfg.mask)
        else:
            mask_vol = np.asarray(cfg.mask)

        if mask_vol.shape != vol1.shape:
            raise ValueError(f"    Mask shape mismatch: {mask_vol.shape} vs maps {vol1.shape}")

        mask_float = np.clip(mask_vol.astype(np.float32, copy=False), 0.0, 1.0)
        thr = 0.5
        mask_bool = (mask_vol >= thr)
        mask_dilate_iters = max(int(getattr(cfg, "mask_dilation_vox", 6)), 0)
        if mask_dilate_iters > 0:
            mask_bool_dilated = binary_dilation(mask_bool, iterations=mask_dilate_iters)
        else:
            mask_bool_dilated = mask_bool.copy()
        mask_float_dilated = mask_bool_dilated.astype(np.float32, copy=False)
        print(f"    Compute mask loaded (dilation={mask_dilate_iters} vox)")

    if getattr(cfg, "output_mask", None) is not None:
        if isinstance(cfg.output_mask, str):
            output_mask_vol, _hdr_out = IO_utils.read_mrc_data(cfg.output_mask)
        else:
            output_mask_vol = np.asarray(cfg.output_mask)
        if output_mask_vol.shape != vol1.shape:
            raise ValueError(f"    Output-mask shape mismatch: {output_mask_vol.shape} vs maps {vol1.shape}")
        output_mask_bool = (np.asarray(output_mask_vol) >= 0.5)
        print("    Output mask loaded")

    compute_bbox = None
    compute_box_bool = None
    output_bbox = None
    output_box_bool = None
    if mask_bool is not None and compute_only_mask:
        compute_bbox = _bbox_from_mask_bool(mask_bool)
        if compute_bbox is None:
            raise ValueError("    Compute mask is empty.")
        compute_box_bool = _bbox_mask(vol1.shape, compute_bbox)

        if output_mask_bool is not None:
            output_bbox = _bbox_from_mask_bool(output_mask_bool)
            if output_bbox is None:
                raise ValueError("    Output mask is empty.")
        else:
            output_bbox = compute_bbox
        output_box_bool = _bbox_mask(vol1.shape, output_bbox)
        print(f"    Boxed-mask mode: compute box z={compute_bbox[0]}, y={compute_bbox[1]}, x={compute_bbox[2]}")
        print(f"    Boxed-mask mode: output  box z={output_bbox[0]}, y={output_bbox[1]}, x={output_bbox[2]}")

    # ---- Mask-support settings for local FSC (avoid smooth-edge artefacts) ----
    # These thresholds are used only for centres that are considered "inside mask"
    # Outside the mask we compute normally.
    mask_thr = float(getattr(cfg, "mask_support_threshold", 0.90))
    min_cov  = float(getattr(cfg, "min_mask_coverage", 0.90))

    # Clamp for safety
    mask_thr = min(max(mask_thr, 0.0), 1.0)
    min_cov  = min(max(min_cov, 0.0), 1.0)


    # ---- Global FSC (unrandomised) OR reuse precomputed ----
    if (cfg.precomputed_sampling is not None and
        cfg.precomputed_radius is not None and
        cfg.precomputed_first_fsc_crossing is not None):
        freqs_g, fsc_g = None, None
        thresholds_all = [cfg.local_fsc_threshold, AUX_THR]
        fcg = {t: (cfg.precomputed_first_fsc_crossing if t == AUX_THR else np.nan) for t in thresholds_all}
        rg  = {t: (1.0/fcg[t] if np.isfinite(fcg[t]) and fcg[t] > 0 else np.inf) for t in thresholds_all}
        sampling_A = float(cfg.precomputed_sampling)
        rad        = float(cfg.precomputed_radius)
        k0143      = float(cfg.precomputed_first_fsc_crossing)
        print(f"    [locres] Reusing precomputed: sampling={sampling_A:.3f} Å, radius={rad:.6f} Å, k0143={k0143:.6f} 1/Å")
    else:
        freqs_g, fsc_g = compute_fsc(vol1, vol2, pix_in, n_bins=cfg.bins)
        thresholds_all = [cfg.local_fsc_threshold, AUX_THR]
        res_g, kcross_g, _, _, _, _ = find_FSC_resolutions_and_stats(freqs_g, fsc_g, thresholds_all)
        if cfg.plot:
            crossings_for_plot = {t: k for t, k in kcross_g.items() if (k is not None and np.isfinite(k))}
            plot_fsc(freqs_g, fsc_g, thresholds_all, crossings_for_plot)

        if isinstance(cfg.sampling, str) and cfg.sampling.strip().lower() == "auto":
            sampling_A = 0.5 * float(res_g[AUX_THR]) if np.isfinite(res_g[AUX_THR]) else max(pix_in, 25.0)
            print(f"    Auto sampling from global 0.143 resolution: {sampling_A:.2f} Å")
        else:
            sampling_A = float(cfg.sampling)
            print(f"    Sampling (from user): {sampling_A:.2f} Å")

        k0143 = first_FSC_crossing(freqs_g, fsc_g, AUX_THR)  # Å^-1
        nyq   = 0.5 / pix_in
        if isinstance(cfg.radius, str) and cfg.radius.lower() == "auto":
            rmax_vox = min(vol1.shape) // 2 - 2
            R_A, R_vox = auto_radius_from_global_fsc(freqs_g, fsc_g, pixel_A=pix_in,
                                                     cycles_min=cfg.cycles, gamma=cfg.gamma,
                                                     rmin_vox=6, rmax_vox=rmax_vox)
            rad = float(R_A)
            print(f"    Auto radius: {rad:.2f} Å  ({R_vox:.1f} vox)")
        else:
            rad = float(cfg.radius)
            R_vox = rad / pix_in
            print(f"    Radius (from user): {rad:.2f} Å  ({R_vox:.1f} vox)")

        rg  = res_g
        fcg = kcross_g

    edge_A = (0.3 * rad) if (cfg.edgewidth is None) else float(cfg.edgewidth)
    edge_vox = max(int(round(edge_A / pix_in)), 1)
    nyq = 0.5 / pix_in

    # ---- Speed-oriented presets ----
    fast_mode = bool(getattr(cfg, "fast", False))
    local_bins = cfg.bins
    local_interp_factor = max(int(getattr(cfg, "local_interp_factor", 10)), 2)
    local_interp_kind = 'cubic'
    full_interp_mode = cfg.interp
    have_precomputed = (
        cfg.precomputed_sampling is not None and
        cfg.precomputed_radius is not None and
        cfg.precomputed_first_fsc_crossing is not None
    )
    if fast_mode:
        if not have_precomputed:
            if isinstance(cfg.sampling, str) and cfg.sampling.strip().lower() == "auto":
                sampling_A = max(float(sampling_A) * 1.5, 3.0 * pix_in)
            if isinstance(cfg.radius, str) and cfg.radius.strip().lower() == "auto":
                rad = max(float(rad) * 0.85, 6.0 * pix_in)
                edge_A = (0.3 * rad) if (cfg.edgewidth is None) else float(cfg.edgewidth)
                edge_vox = max(int(round(edge_A / pix_in)), 1)
        if local_bins is None:
            local_bins = max(12, int(round((rad / pix_in) * 0.5)))
        local_interp_factor = min(local_interp_factor, 4)
        local_interp_kind = 'linear'
        if full_interp_mode == 'cubic':
            full_interp_mode = 'linear'
        print(
            f"    Fast/default mode: sampling={sampling_A:.2f} Å, radius={rad:.2f} Å, "
            f"bins={local_bins if local_bins is not None else 'auto'}, "
            f"local_interp={local_interp_kind}, full_interp={full_interp_mode}"
        )

    # ---- Phase randomisation cutoffs ----
    kcut_global = resolve_kcut(cfg.rand_global, pix_in, k0143, nyq, frac_nyq=cfg.rand_frac_nyq)
    kcut_local  = resolve_kcut(cfg.rand_local,  pix_in, k0143, nyq, frac_nyq=cfg.rand_frac_nyq)
    if isinstance(cfg.rand_local, str) and cfg.rand_local.strip().lower().startswith("auto"):
        if kcut_global is not None:
            kcut_local = kcut_global

    #print(f"[kcut debug] pix={pix_in:.4f} Å  nyq={nyq:.6f}  k0143={k0143:.6f}")
    #print(f"[kcut debug] rand_global={cfg.rand_global}  rand_local={cfg.rand_local}")
    #print(f"[kcut debug] kcut_global={kcut_global}  kcut_local={kcut_local}")

    # ---- Apply global randomisation for analysis ----
    vol1_work = vol1
    vol2_work = vol2
    if kcut_global is not None:
        # Push randomisation slightly towards worse resolution
        K_SHIFT = 1.15  # 1.2–1.3 is sensible; higher = more aggressive
        kcut_eff = kcut_global / K_SHIFT

        print(
            f"    Global phase randomisation beyond {1.0/kcut_eff:.2f} Å (k={kcut_eff:.3f} Å^-1)"
        )
        vol1_work = phase_randomize_beyond(vol1, kcut_eff, pix_in)
        vol2_work = phase_randomize_beyond(vol2, kcut_eff, pix_in)

    # ---- Centres ----
    centres_all = get_region_centres(vol1.shape, pix_in, sampling_A, rad)

    if compute_box_bool is not None:
        centres = [(x, y, z) for (x, y, z) in centres_all if _bbox_contains_point(compute_bbox, x, y, z)]
        xs = sorted({c[0] for c in centres})
        ys = sorted({c[1] for c in centres})
        zs = sorted({c[2] for c in centres})
        n_x, n_y, n_z = len(xs), len(ys), len(zs)
        xi = {x: i for i, x in enumerate(xs)}
        yi = {y: j for j, y in enumerate(ys)}
        zi = {z: k for k, z in enumerate(zs)}
        print(f"    Centres within compute box only: {len(centres)}/{len(centres_all)}")
    else:
        # Build the sampling grid axes from the full lattice (so coarse maps cover the full lattice)
        xs = sorted({c[0] for c in centres_all})
        ys = sorted({c[1] for c in centres_all})
        zs = sorted({c[2] for c in centres_all})
        n_x, n_y, n_z = len(xs), len(ys), len(zs)
        xi = {x: i for i, x in enumerate(xs)}
        yi = {y: j for j, y in enumerate(ys)}
        zi = {z: k for k, z in enumerate(zs)}

        if mask_bool is not None:
            outside_factor = max(1, int(getattr(cfg, "outside_mask_factor", 8)))

            centres_in: List[Tuple[int, int, int]] = []
            centres_out: List[Tuple[int, int, int]] = []

            for (x, y, z) in centres_all:
                if mask_bool[z, y, x]:
                    centres_in.append((x, y, z))
                else:
                    if outside_factor == 1:
                        centres_out.append((x, y, z))
                    else:
                        # stride on the *sampling lattice* indices so outside points remain commensurate
                        if ((xi[x] % outside_factor) == 0 and
                            (yi[y] % outside_factor) == 0 and
                            (zi[z] % outside_factor) == 0):
                            centres_out.append((x, y, z))

            centres = centres_in + centres_out
            print(f"    Centres inside mask: {len(centres_in)}/{len(centres_all)}")
            print(f"    Centres outside mask (sparse, factor={outside_factor}): {len(centres_out)}/{len(centres_all)}")
        else:
            centres = centres_all

    total = len(centres)
    if total == 0:
        print("    No region centres to process (mask too small or sampling too coarse). Exiting.")
        out = {
            'maps': {}, 'coarse_grids': {},
            'global_fsc': {'freqs': freqs_g, 'fsc': fsc_g, 'resolutions': rg, 'crossings': fcg},
            'params': {
                'pixel_size_A': pix_in,
                'sampling_A': sampling_A,
                'radius_A': rad,
                'edge_A': edge_A,
                'step_vox': None,
                'centres_count': 0,
                'compute_only_mask': bool(compute_only_mask),
                'fast': bool(fast_mode),
            },
            'stats_per_region': {}
        }
        if cfg.output_basename:
            params_fn = f"{cfg.output_basename}_params.txt"
            with open(params_fn, "w") as pf:
                pf.write(f"pixel_size_A\t{pix_in:.6f}\n")
                pf.write(f"sampling_A\t{sampling_A:.6f}\n")
                pf.write(f"radius_A\t{rad:.6f}\n")
            print(f"    Wrote {params_fn}")
        return out

    print(f"    Total regions to compute: {total}")

    # ---- Parallel execution (CPU only) ----
    cpu = int(cfg.cpu)
    cpu_max = mp.cpu_count()
    cpu = min(max(cpu, 1), max(1, min(cpu_max, total)))

    try:
        mp.set_start_method("fork")
    except RuntimeError:
        pass

    block = int(math.floor(total / cpu)) if cpu > 0 else total
    idx_blocks = []
    for i in range(cpu):
        b0 = i * block
        b1 = (i + 1) * block if i < cpu - 1 else total
        idx_blocks.append((b0, b1))

    def _worker(b0: int, b1: int, procnum: int, ret_dict, progress):
        results = []
        stats_per = {}
        rng = np.random.default_rng(12345 + procnum)
        FLUSH_EVERY = 64
        acc = 0
        for idx in range(b0, b1):
            cx, cy, cz = centres[idx]
            # Apply hard-support masking only for "inside-mask" centres.
            # Outside-mask centres should be computed without mask gating,
            # otherwise they will almost always fail min_coverage and return None.
            mvol = None
            if (mask_bool_dilated is not None) and mask_bool_dilated[cz, cy, cx]:
                mvol = mask_float_dilated

            stats = compute_local_fsc_region(
                vol1_work, vol2_work, cx, cy, cz, pix_in, rad,
                thresholds_all, n_bins=local_bins,
                interp_factor=local_interp_factor, interp_kind=local_interp_kind,
                edge_vox=edge_vox, kcut_local=kcut_local, rng=rng,
                mask_vol=mvol,
                mask_thr=mask_thr,
                min_coverage=min_cov,
            )

            if stats is None:
                continue
            results.append((cx, cy, cz, stats['resolutions'], stats['mean_val']))
            stats_per[(cx, cy, cz)] = stats

            acc += 1
            if (acc & (FLUSH_EVERY - 1)) == 0:
                with progress.get_lock():
                    progress.value += acc
                acc = 0

        if acc:
            with progress.get_lock():
                progress.value += acc

        ret_dict[procnum] = (b0, b1, results, stats_per)

    progress = mp.Value('i', 0)
    manager = mp.Manager()
    ret_dict = manager.dict()

    procs = []
    for pi, (b0, b1) in enumerate(idx_blocks):
        p = mp.Process(target=_worker, args=(b0, b1, pi, ret_dict, progress))
        procs.append(p)

    start = time.time()
    for p in procs:
        p.start()

    last = -1
    while any(p.is_alive() for p in procs):
        done = progress.value
        if done != last:
            elapsed = time.time() - start
            pct = 100.0 * done / max(1, total)
            rate = done / elapsed if elapsed > 0 else 0.0
            remain = (total - done) / rate if rate > 0 else 0.0
            print(
                f"\rProcessed {done:,}/{total:,} ({pct:5.1f}%) — "
                f"{rate:,.1f} regions/s — elapsed {time.strftime('%H:%M:%S', time.gmtime(elapsed))} — "
                f"ETA {time.strftime('%H:%M:%S', time.gmtime(remain))}",
                end=""
            )
            last = done
        time.sleep(0.2)
    for p in procs:
        p.join()
    print("")

    results: List[Tuple[int, int, int, Dict[float, float], float]] = []
    stats_per_region: Dict[Tuple[int, int, int], dict] = {}
    for pi in range(len(idx_blocks)):
        b0, b1, res_block, stats_block = ret_dict[pi]
        results.extend(res_block)
        stats_per_region.update(stats_block)

    # ---- Logging table (always build in-memory; write if requested) ----
    log_rows = []
    hdr_row = ['x', 'y', 'z']
    for t in thresholds_all:
        hdr_row += [f'f_cross@{t}', f'res@{t}']
    hdr_row += ['mean_between', 'var_between', 'raw_fsc']
    log_rows.append('\t'.join(hdr_row))
    for (cx, cy, cz, _resdict, _) in results:
        stats = stats_per_region[(cx, cy, cz)]
        row = [str(cx), str(cy), str(cz)]
        for t in thresholds_all:
            row += [f"{stats['f_crossings'][t]:.6f}", f"{stats['resolutions'][t]:.3f}"]
        row += [f"{stats['mean_between']:.4f}", f"{stats['var_between']:.4f}"]
        row += [":".join(f"{v:.6f}" for v in stats['fsc_vals'])]
        log_rows.append('\t'.join(row))

    # ---- Build coarse maps on the FULL sampling lattice ----
    locmap = np.zeros((n_z, n_y, n_x), dtype=np.float32)
    avmap  = np.zeros_like(locmap)
    invvarmap = np.zeros_like(locmap)

    thr_min = min(cfg.local_fsc_threshold, AUX_THR)

    # Track which lattice nodes have valid locres samples (>0)
    sampled = np.zeros((n_z, n_y, n_x), dtype=bool)

    for cx, cy, cz, resdict, meanval in results:
        i, j, k = xi[cx], yi[cy], zi[cz]
        val = resdict.get(cfg.local_fsc_threshold, None)
        if val is None:
            val_f = 0.0
        else:
            val_f = float(val)
        locmap[k, j, i] = val_f
        avmap[k, j, i]  = float(meanval)

        raw_fsc = stats_per_region[(cx, cy, cz)]['fsc_vals']
        below = np.where(raw_fsc < thr_min)[0]
        end = int(below[0]) if below.size else (len(raw_fsc) - 1)
        seg = raw_fsc[:end + 1]
        var_r = float(np.var(seg)) if seg.size else 0.0
        invvarmap[k, j, i] = (1.0 / var_r) if (var_r > 0) else -1.0

        # Mark sampled if we got a non-zero resolution estimate
        if np.isfinite(val_f) and (val_f > 0.0):
            sampled[k, j, i] = True
            
    # ---- Fill unsampled coarse nodes (prevents holes after upsampling) ----
    sampled = np.isfinite(locmap) & (locmap > 0)

    if np.any(sampled) and np.any(~sampled):
        # distance_transform_edt returns, for each hole voxel, the indices of the nearest sampled voxel
        _, idx = distance_transform_edt(~sampled, return_distances=True, return_indices=True)

        nn_loc = locmap[idx[0], idx[1], idx[2]]
        nn_av  = avmap[idx[0], idx[1], idx[2]]
        nn_inv = invvarmap[idx[0], idx[1], idx[2]]

        locmap[~sampled] = nn_loc[~sampled]
        avmap[~sampled]  = nn_av[~sampled]
        invvarmap[~sampled] = nn_inv[~sampled]

    elif not np.any(sampled):
        # Degenerate case: no valid samples at all -> plateau to something deterministic
        vals = locmap[np.isfinite(locmap) & (locmap > 0)]
        plateau = float(np.median(vals)) if vals.size else float("inf")
        locmap[:] = plateau
        # avmap/invvarmap can stay as-is (zeros), they are auxiliary


    # ---- Upsampling or resampled outputs ----
    origin_in = np.array(
        [
            float(hdr1.get("origin_x", 0.0)),
            float(hdr1.get("origin_y", 0.0)),
            float(hdr1.get("origin_z", 0.0)),
        ],
        dtype=np.float64,
    )

    if cfg.resample:
        out_pix = float(sampling_A)

        x0, y0, z0 = xs[0], ys[0], zs[0]
        centre_in = origin_in + np.array(
            [(x0 + 0.5) * pix_in, (y0 + 0.5) * pix_in, (z0 + 0.5) * pix_in],
            dtype=np.float64,
        )
        origin_out = centre_in + out_pix

        locres_grid = locmap
        write_specs = dict(pixel_size=out_pix, origin_angs=tuple(origin_out))
    else:
        # Full-grid output: upsample FILLED coarse grid with LINEAR interpolation only
        step_vox = max(1, int(round(float(sampling_A) / pix_in)))
        z_start, y_start, x_start = zs[0], ys[0], xs[0]
        z_len = (len(zs) - 1) * step_vox + 1
        y_len = (len(ys) - 1) * step_vox + 1
        x_len = (len(xs) - 1) * step_vox + 1

        # Interpolate the filled coarse grid back to full voxel spacing.
        zoom_factors = (float(step_vox), float(step_vox), float(step_vox))
        interp_order = 1 if full_interp_mode == "linear" else 3
        up_loc = ndi_zoom(
            locmap.astype(np.float32, copy=False),
            zoom=zoom_factors,
            order=interp_order,
            prefilter=(interp_order > 1),
            mode="nearest"
        )

        up_loc = up_loc[:z_len, :y_len, :x_len]

        full_shape = vol1.shape
        loc_full = embed_into_full_grid(
            up_loc,
            (z_start, y_start, x_start),
            full_shape
        )

        # ---- Nearest-neighbour fill only where output is actually kept ----
        nz = np.isfinite(loc_full) & (loc_full > 0.0)
        if output_box_bool is not None:
            fill_region = output_box_bool
        elif output_mask_bool is not None:
            fill_region = output_mask_bool
        elif compute_box_bool is not None:
            fill_region = compute_box_bool
        elif mask_bool is not None:
            fill_region = mask_bool
        else:
            fill_region = np.ones_like(loc_full, dtype=bool)
        holes = fill_region & (~nz)

        if np.any(holes):
            if np.any(nz):
                _, idx = distance_transform_edt(
                    ~nz, return_distances=True, return_indices=True
                )
                loc_full[holes] = loc_full[idx[0], idx[1], idx[2]][holes]
            else:
                vals = locmap[np.isfinite(locmap) & (locmap > 0.0)]
                plateau = float(np.median(vals)) if vals.size else float("inf")
                loc_full[fill_region] = plateau

        locres_grid = loc_full
        write_specs = dict(pixel_size=pix_in, origin_angs=tuple(origin_in))

    # ---- For boxed-mask runs: force outside the final output box to a high locres value ----
    box_mask_src = output_box_bool if output_box_bool is not None else output_mask_bool
    if box_mask_src is not None:
        if cfg.resample:
            box_mask_out = box_mask_src[np.ix_(zs, ys, xs)]
        else:
            box_mask_out = box_mask_src

        inside_vals = locres_grid[box_mask_out]
        inside_vals = inside_vals[np.isfinite(inside_vals)]

        if inside_vals.size > 0:
            high_locres = max(float(np.max(inside_vals)), float(getattr(cfg, "outside_box_fill_A", 100.0)))
        else:
            high_locres = float(getattr(cfg, "outside_box_fill_A", 100.0))
        locres_grid[~box_mask_out] = high_locres


    # ---- Write outputs if requested ----
    if cfg.output_basename:
        base = cfg.output_basename

        logfn = f"{base}_log.txt"
        with open(logfn, "w") as log:
            log.write("\n".join(log_rows) + "\n")
        print(f"Wrote {logfn}")

        if cfg.resample:
            fn1 = f"{base}_resampledLocres.mrc"
        else:
            fn1 = f"{base}_locres.mrc"
        IO_utils.write_mrc(fn1, locres_grid, **write_specs)
        print("Wrote", fn1)

        if (freqs_g is not None) and (fsc_g is not None):
            fn4 = f"{base}_globalFSC.txt"
            with open(fn4, "w") as f:
                f.write("freq_A^-1\tFSC\n")
                for fr, fv in zip(freqs_g, fsc_g):
                    f.write(f"{float(fr):.6f}\t{float(fv):.6f}\n")
            print("Wrote", fn4)
        else:
            print("[locres] Skipped globalFSC.txt (reused precomputed global parameters).")

        params_fn = f"{base}_params.txt"
        with open(params_fn, "w") as pf:
            pf.write(f"pixel_size_A\t{pix_in:.6f}\n")
            pf.write(f"sampling_A\t{sampling_A:.6f}\n")
            pf.write(f"radius_A\t{rad:.6f}\n")
        print(f"Wrote {params_fn} (sampling_A={sampling_A:.3f} Å, radius_A={rad:.3f} Å)")

    # ---- Assemble return object ----
    out = {
        "maps": {
            "locres": locres_grid,
            "mean": avmap,
            "inverse_variance": invvarmap,
        },
        "coarse_grids": {
            "locmap": locmap,
            "avmap": avmap,
            "invvarmap": invvarmap,
        },
        "global_fsc": {
            "freqs": freqs_g,
            "fsc": fsc_g,
            "resolutions": rg,
            "crossings": fcg,
        },
        "params": {
            "pixel_size_A": pix_in,
            "sampling_A": sampling_A,
            "radius_A": rad,
            "edge_A": edge_A,
            "centres_count": total,
            "interp_mode": full_interp_mode,
            "resample": cfg.resample,
            "outside_mask_factor": int(getattr(cfg, "outside_mask_factor", 8)) if (mask_bool is not None and not compute_only_mask) else None,
            "compute_only_mask": bool(compute_only_mask),
            "fast": bool(fast_mode),
        },
        "stats_per_region": stats_per_region,
        "log_rows": log_rows,
    }
    print("]")

    return out




# ----------------------- Recon-orchestration helpers -----------------------
def _parse_filename_halfmap(path: str):
    """
    From '..._<COUNT>_recH{1|2}.mrc' or '..._recH{1|2}.mrc' extract:
      parent_dir (Path), base_stem (without _recH# or _<COUNT>_recH#), count (int|None), subset (1|2|None)
    """
    p = Path(path)
    m = re.match(r"^(.*)_(\d+)_recH([12])\.mrc$", p.name)
    if m:
        return p.parent, m.group(1), int(m.group(2)), int(m.group(3))
    m = re.match(r"^(.*)_recH([12])\.mrc$", p.name)
    if m:
        return p.parent, m.group(1), None, int(m.group(2))
    return p.parent, p.stem, None, None



def _aggregate_locres_set(scan_dir: Path,
                          base_prefix: str,
                          labels_sorted: List[int],
                          mask_path: Optional[str]) -> None:
    """
    For a set of checkpoint labels (particle counts) that produced:
       <base_prefix>_<COUNT>_locres.mrc
       <base_prefix>_<COUNT>_recH1.mrc
       <base_prefix>_<COUNT>_recH2.mrc
    produce:
       stats_locres.mrc
       stats_minParticles.mrc
       stats_intensity_recH1.mrc
       stats_intensity_recH2.mrc
       bestRanked_locres_values.csv
    """
    # Locate a template header from any half-map
    tmpl_hdr = None
    tmpl_shape = None

    # Optional mask
    mask_bool = None
    if mask_path:
        mvol, mh = IO_utils.read_mrc_data(mask_path)
        mask_bool = (mvol >= 0.5)
        mask_bool = binary_dilation(mask_bool, iterations=2)

    # Pre-scan to determine shape/header
    for lbl in labels_sorted:
        h1p = scan_dir / f"{base_prefix}_{lbl}_recH1.mrc"
        if h1p.exists():
            _, tmpl_hdr = IO_utils.read_mrc_data(str(h1p))
            tmpl_shape = (int(tmpl_hdr["nz"]), int(tmpl_hdr["ny"]), int(tmpl_hdr["nx"]))
            break
    if tmpl_hdr is None or tmpl_shape is None:
        print(f"[stats] No half-maps found for prefix {base_prefix}. Skipping aggregation.")
        return

    # Accumulators
    best_map = np.full(tmpl_shape, np.inf, dtype=np.float32)  # best (=min) resolution per voxel
    best_count = np.zeros(tmpl_shape, dtype=np.float32)       # particle count achieving best
    best_h1 = np.zeros(tmpl_shape, dtype=np.float32)          # intensity from H1 at best
    best_h2 = np.zeros(tmpl_shape, dtype=np.float32)          # intensity from H2 at best

    # CSV rows
    csv_rows = ["numParticles,max,highQuartile,mean,lowQuartile,min"]

    # Sweep labels in ascending order of particle count
    for lbl in sorted(labels_sorted):
        loc_path = scan_dir / f"{base_prefix}_{lbl}_locres.mrc"
        h1_path  = scan_dir / f"{base_prefix}_{lbl}_recH1.mrc"
        h2_path  = scan_dir / f"{base_prefix}_{lbl}_recH2.mrc"
        if not (loc_path.exists() and h1_path.exists() and h2_path.exists()):
            continue

        locres, _ = IO_utils.read_mrc_data(str(loc_path))
        h1, _     = IO_utils.read_mrc_data(str(h1_path))
        h2, _     = IO_utils.read_mrc_data(str(h2_path))

        # Safety: enforce shapes
        if locres.shape != tmpl_shape or h1.shape != tmpl_shape or h2.shape != tmpl_shape:
            raise ValueError(f"[stats] Shape mismatch at count {lbl}: "
                             f"locres={locres.shape}, h1={h1.shape}, h2={h2.shape}, expected={tmpl_shape}")

        # Per-voxel update where current locres is better (numerically smaller)
        better = locres < best_map
        if np.any(better):
            best_map[better]   = locres[better]
            best_count[better] = float(lbl)
            best_h1[better]    = h1[better]
            best_h2[better]    = h2[better]

        # Per-checkpoint masked stats line (note: 'max' = best = min())
        best_v, q25, mean_v, q75, worst_v = _masked_stats_for_locres(locres, mask_bool)
        csv_rows.append(f"{int(lbl)},{best_v:.5f},{q25:.5f},{mean_v:.5f},{q75:.5f},{worst_v:.5f}")

    # If nothing was updated, bail politely
    if not np.isfinite(best_map).any():
        print(f"[stats] No locres maps found for prefix {base_prefix}.")
        return

    # Write four MRCs using the template header
    IO_utils.write_mrc_like(str(scan_dir / "stats_locres.mrc"),           best_map.astype(np.float32),   tmpl_hdr, update_stats=True)
    IO_utils.write_mrc_like(str(scan_dir / "stats_minParticles.mrc"),     best_count.astype(np.float32), tmpl_hdr, update_stats=True)
    IO_utils.write_mrc_like(str(scan_dir / "stats_intensity_recH1.mrc"),  best_h1.astype(np.float32),    tmpl_hdr, update_stats=True)
    IO_utils.write_mrc_like(str(scan_dir / "stats_intensity_recH2.mrc"),  best_h2.astype(np.float32),    tmpl_hdr, update_stats=True)

    # Write CSV
    csv_path = scan_dir / "bestRanked_locres_values.csv"
    with open(csv_path, "w") as f:
        f.write("\n".join(csv_rows) + "\n")
    print(f"[stats] Wrote {csv_path.name}, stats_locres.mrc, stats_minParticles.mrc, stats_intensity_recH1.mrc, stats_intensity_recH2.mrc in {scan_dir}")


def run_locres_for_all_pairs(results: List[Tuple[str, bool, str, Optional[int], Optional[int]]],
                             preferred_dir: Optional[Path],
                             args) -> None:
    """
    Drop-in for the reconstructor: find every H1/H2 checkpoint pair sharing the SAME <COUNT>
    label and run locres on each pair. Uses options provided in `args`:
      - args.locres_mask
      - args.locres_cpu
      - args.locres_stats  (if True, aggregate a best-of stats set)

    The function scans next to the final outputs so that intermediate checkpoint files are picked up.
    """
    # Gather candidate (class_id -> list of final half-map paths we can scan around)
    finals_by_class: Dict[Optional[int], List[str]] = {}
    for outp, ok, _, cid, sid in results:
        if ok and sid in (1, 2):
            finals_by_class.setdefault(cid, []).append(outp)

    if not finals_by_class:
        print("[locres] Skipped: no half-maps reported by reconstruction.")
        return

    # For each class (or None), look in the directory of any final half-map,
    # and collect all files matching "<base>_<COUNT>_recH{1|2}.mrc"
    for cls_key, sample_paths in finals_by_class.items():
        # Use the first sample path to derive where and what to scan
        parent, base_stem, _, _ = _parse_filename_halfmap(sample_paths[0])
        scan_dir = preferred_dir if preferred_dir else parent

        # If --basename was given, rebuild the expected stem (matches naming elsewhere)
        if getattr(args, 'basename', None):
            base_prefix = f"{args.basename}" + (f"_class_{cls_key}" if cls_key is not None else "")
        else:
            star_stem = Path(args.star_file).stem
            base_prefix = f"{star_stem}" + (f"_class_{cls_key}" if cls_key is not None else "")

        # Glob both halves with count labels
        h1 = {}
        h2 = {}
        for p in Path(scan_dir).glob(f"{base_prefix}_*_recH1.mrc"):
            _, _, count_label, subset = _parse_filename_halfmap(str(p))
            if count_label is not None and subset == 1:
                h1[int(count_label)] = str(p)
        for p in Path(scan_dir).glob(f"{base_prefix}_*_recH2.mrc"):
            _, _, count_label, subset = _parse_filename_halfmap(str(p))
            if count_label is not None and subset == 2:
                h2[int(count_label)] = str(p)

        common_labels = sorted(set(h1.keys()).intersection(h2.keys()), reverse=True)
        if not common_labels:
            plain_h1 = Path(scan_dir) / f"{base_prefix}_recH1.mrc"
            plain_h2 = Path(scan_dir) / f"{base_prefix}_recH2.mrc"
            if plain_h1.exists() and plain_h2.exists():
                out_base = str(Path(scan_dir) / f"{base_prefix}")
                print(f"[locres] Computing local resolution for "
                      f"{'class '+str(cls_key) if cls_key is not None else 'all classes'} "
                      f"using H1={plain_h1.name}, H2={plain_h2.name} → basename {out_base}")
                try:
                    cfg = LocResConfig(
                        halfmap1=str(plain_h1),
                        halfmap2=str(plain_h2),
                        mask=getattr(args, 'locres_mask', None),
                        cpu=int(getattr(args, 'locres_cpu', 1)),
                        output_basename=out_base,
                    )
                    _ = locres_map(cfg)
                except Exception as e:
                    print(f"[locres] ERROR (plain pair): {e}")
            else:
                print(f"[locres] No matching H1/H2 checkpoint pairs found for "
                      f"{'class '+str(cls_key) if cls_key is not None else 'all classes'} in {scan_dir}.")
            continue

        global_locres_params = None

        for lbl in common_labels:
            half1 = h1[lbl]
            half2 = h2[lbl]
            out_base = str(Path(scan_dir) / f"{base_prefix}_{lbl}")
            print(f"[locres] Computing local resolution for "
                  f"{'class '+str(cls_key) if cls_key is not None else 'all classes'} "
                  f"using H1={Path(half1).name}, H2={Path(half2).name} → basename {out_base}")
            try:
                if global_locres_params is None:
                    # Reference run on largest particle count: compute parameters
                    cfg = LocResConfig(
                        halfmap1=half1,
                        halfmap2=half2,
                        mask=getattr(args, 'locres_mask', None),
                        cpu=int(getattr(args, 'locres_cpu', 1)),
                        output_basename=out_base,
                    )
                    res = locres_map(cfg)
                    params = res.get('params', {})
                    global_locres_params = {
                        'sampling': float(params.get('sampling_A')),
                        'radius': float(params.get('radius_A')),
                        # Prefer the 0.143 crossing if present; else use NaN sentinel
                        'first_fsc_crossing': float(res.get('global_fsc', {}).get('crossings', {}).get(0.143, np.nan)),
                        'mask': getattr(args, 'locres_mask', None),
                    }
                    if not np.isfinite(global_locres_params['first_fsc_crossing']):
                        print('[locres] Warning: reference run did not yield a finite 0.143 crossing; reuse will still use sampling/radius.')
                else:
                    # Reuse mode: pass precomputed parameters so no global FSC is recomputed
                    cfg = LocResConfig(
                        halfmap1=half1,
                        halfmap2=half2,
                        mask=global_locres_params['mask'],
                        cpu=int(getattr(args, 'locres_cpu', 1)),
                        output_basename=out_base,
                        precomputed_sampling=global_locres_params['sampling'],
                        precomputed_radius=global_locres_params['radius'],
                        precomputed_first_fsc_crossing=global_locres_params['first_fsc_crossing'],
                    )
                    _ = locres_map(cfg)
            except Exception as e:
                print(f"[locres] ERROR (count {lbl}): {e}")


        if getattr(args, 'locres_stats', False) and common_labels:
            try:
                _aggregate_locres_set(
                    scan_dir=scan_dir,
                    base_prefix=base_prefix,
                    labels_sorted=sorted(common_labels),
                    mask_path=getattr(args, 'locres_mask', None),
                )
            except Exception as e:
                print(f"[stats] ERROR while aggregating stats for {base_prefix}: {e}")


def write_resampled_mask_from_locres_result(
    locres_result: Dict[str, Any],
    like_basename: str,
    out_mask_path: str,
) -> None:
    """
    Build a binary mask on the same coarse grid as the resampled locres map,
    and write it as an MRC file.

    Assumes locres was run with cfg.resample=True for this basename, so that
    {like_basename}_resampledLocres.mrc was written.

    Mask definition: 1 where locres has finite, non-zero values; 0 elsewhere.
    """
    maps = locres_result.get("maps", {})
    locmap = maps.get("locres", None)
    if locmap is None:
        print("[locresBulk] Cannot derive mask: 'locres' map not found in result.")
        return

    # Binary mask on same grid as locres map
    mask = np.zeros_like(locmap, dtype=np.float32)
    finite = np.isfinite(locmap)
    nonzero = np.abs(locmap) > 0.0
    mask[finite & nonzero] = 1.0

    like_mrc = f"{like_basename}_resampledLocres.mrc"
    try:
        _, hdr_like = IO_utils.read_mrc_data(like_mrc)
    except Exception as e:
        print(
            f"[locresBulk] WARNING: could not read '{like_mrc}' to copy header for mask: {e}"
        )
        # Fallback: write with pixel size from params, if available
        pix = float(locres_result.get("params", {}).get("sampling_A", 1.0))
        IO_utils.write_mrc(out_mask_path, mask, pixel_size=pix)
        print(
            f"[locresBulk] Wrote resampled mask (fallback header) to {out_mask_path}"
        )
        return

    # Write mask with same header (dimensions, pixel size, origin) as the resampled locres map
    IO_utils.write_mrc_like(
        out_mask_path, mask.astype(np.float32), hdr_like, update_stats=True
    )
    print(
        f"[locresBulk] Wrote resampled mask to {out_mask_path} "
        f"(same grid as {like_mrc})"
    )


def run_locres_bulk(
    prefix: str,
    items_csv: str,
    suffixes_csv: Optional[str],
    args,
) -> None:
    """
    Bulk local-resolution over an explicit set of half-map pairs.

    Given:
      - prefix:       common prefix for all half-maps, e.g. 'ciccio_'
      - items_csv:    comma-separated list of numeric labels, e.g. '34,455,2233,4455'
      - suffixes_csv: comma-separated pair of suffixes for the two halves,
                      e.g. '_recH1.mrc,_recH2.mrc' (if None, this default is used)

    This function looks for half-map pairs

        H1 = {prefix}{item}{suffix1}
        H2 = {prefix}{item}{suffix2}

    for every item in items_csv, and runs locres on each pair.

    Behaviour:

      * The items are interpreted as integers and sorted descending; the largest
        label for which BOTH half-maps exist is used as the reference pair.

      * For the reference pair, locres_map is run in full, yielding sampling_A,
        radius_A and the auxiliary FSC crossing at AUX_THR. Those parameters are
        cached in global_locres_params.

      * For all subsequent pairs, locres_map is called with
            precomputed_sampling
            precomputed_radius
            precomputed_first_fsc_crossing
        so the same local-resolution parameters are reused.

      * The CLI flag --resample is interpreted as:
            --resample
                → args.resample is True   → resample_flag=True, no mask output
            --resample maskResampled.mrc
                → args.resample is 'maskResampled.mrc'
                   resample_flag=True, resample_mask_out='maskResampled.mrc'
    """
    skip_locinfo=True
    # -------- Parse suffixesCSV --------
    if not suffixes_csv:
        suffixes_csv = "_recH1.mrc,_recH2.mrc"
    suffix_parts = [s.strip() for s in suffixes_csv.split(",") if s.strip()]
    if len(suffix_parts) != 2:
        raise ValueError(
            "[locresBulk] suffixesCSV must contain exactly two comma-separated entries, "
            f"e.g. '_recH1.mrc,_recH2.mrc'; got: {suffixes_csv!r}"
        )
    suffix1, suffix2 = suffix_parts

    # -------- Parse itemsCSV into numeric labels --------
    item_strs = [s.strip() for s in items_csv.split(",") if s.strip()]
    if not item_strs:
        print("[locresBulk] No items provided; nothing to do.")
        return

    labels: List[Tuple[int, str]] = []
    for s in item_strs:
        try:
            v = int(s)
        except Exception:
            raise ValueError(
                f"[locresBulk] Cannot interpret item '{s}' as an integer."
            )
        labels.append((v, s))

    # Sort by numeric label, largest first → first existing pair is the reference
    labels_sorted = sorted(labels, key=lambda t: t[0], reverse=True)

    # -------- Reference FSC profile (radial amplitude) --------
    ref_fsc_freqs = None
    ref_fsc_vals = None

    # -------- Optional: global half-map analysis CSV --------
    csv_path = getattr(args, "include_global_analysis", None)
    if csv_path:
        rows = []
        best_ref = None  # will hold best combination of amp/PSF scores

        for numeric_label, label_str in labels_sorted:
            half1 = f"{prefix}{label_str}{suffix1}"
            half2 = f"{prefix}{label_str}{suffix2}"

            p1 = Path(half1)
            p2 = Path(half2)
            if not (p1.exists() and p2.exists()):
                print(
                    f"[locresBulk] (global-analysis) Skipping label {label_str}: "
                    f"missing half-maps '{half1}' or '{half2}'."
                )
                continue

            try:
                (
                    global_score,
                    res_143,
                    fsc_score,
                    amp_score,
                    iso_score,
                    freqs_valid,
                    fsc_valid,
                ) = janas_mapProcess.score_halfmaps_from_files(
                    hmap1_fname=half1,
                    hmap2_fname=half2,
                    report_path=None,
                    verbose=False,
                )

                # numeric_label is typically the number of particles in that subset
                rows.append(
                    (
                        int(numeric_label),
                        float(res_143),
                        float(fsc_score),
                        float(amp_score),
                        float(iso_score),
                        float(global_score),
                    )
                )

                # Choose reference: best combination of amplitude + PSF scores
                # here we use the product amp_score * iso_score as the metric.
                metric = float(amp_score) * float(iso_score)
                if best_ref is None or metric > best_ref["metric"]:
                    best_ref = {
                        "metric": metric,
                        "num_particles": int(numeric_label),
                        "label": label_str,
                        "freqs": freqs_valid,
                        "fsc": fsc_valid,
                        "amp_score": float(amp_score),
                        "iso_score": float(iso_score),
                        "global_score": float(global_score),
                        "res_143": float(res_143),
                    }

            except Exception as e:
                print(
                    f"[locresBulk] (global-analysis) ERROR for label {label_str}: {e}"
                )

        if rows:
            # Sort rows by numParticles ascending (to match previous behaviour)
            rows_sorted = sorted(rows, key=lambda r: r[0], reverse=False)
            try:
                with open(csv_path, "w", newline="", encoding="utf-8") as fh:
                    writer = csv.writer(fh)
                    writer.writerow(
                        [
                            "numParticles",
                            "FSC",
                            "FSC_score",
                            "ampl_falloff_score",
                            "PSF_score",
                            "global_score",
                        ]
                    )
                    for row in rows_sorted:
                        writer.writerow(row)
                print(
                    f"[locresBulk] Global half-map analysis written to '{csv_path}'."
                )
            except Exception as e:
                print(
                    f"[locresBulk] ERROR: could not write global-analysis CSV "
                    f"to '{csv_path}': {e}"
                )
        else:
            print(
                "[locresBulk] Global half-map analysis requested, but no valid "
                "half-map pairs were found to score."
            )

        # If we found a best reference pair, store its FSC curve
        if best_ref is not None:
            ref_fsc_freqs = best_ref["freqs"]
            ref_fsc_vals = best_ref["fsc"]
            print(
                "[locresBulk] Using label {lab} (numParticles={n}) as reference "
                "for radial amplitude profile (amp_score={a:.3f}, PSF_score={p:.3f}, "
                "global_score={g:.3f}, FSC_0.143={r:.2f} Å).".format(
                    lab=best_ref["label"],
                    n=best_ref["num_particles"],
                    a=best_ref["amp_score"],
                    p=best_ref["iso_score"],
                    g=best_ref["global_score"],
                    r=best_ref["res_143"],
                )
            )
        else:
            print(
                "[locresBulk] No suitable reference pair found for radial amplitude profile; "
                "local resolution will use raw half-maps."
            )

    # -------- Interpret --resample (flag or filename) --------
    resample_arg = getattr(args, "resample", False)
    resample_flag = bool(resample_arg)
    resample_mask_out = resample_arg if isinstance(resample_arg, str) else None

    global_locres_params = None
    any_done = False

    for numeric_label, label_str in labels_sorted:
        half1 = f"{prefix}{label_str}{suffix1}"
        half2 = f"{prefix}{label_str}{suffix2}"

        p1 = Path(half1)
        p2 = Path(half2)
        if not (p1.exists() and p2.exists()):
            print(
                f"[locresBulk] Skipping label {label_str}: "
                f"missing half-maps '{half1}' or '{half2}'."
            )
            continue

        out_base = f"{prefix}{label_str}"

        print(
            f"[locresBulk] Computing local resolution for label {label_str} "
            f"using H1={p1.name}, H2={p2.name} → basename {out_base}"
        )

        try:
            if global_locres_params is None:
                # ---------------- Reference run: compute parameters ----------------
                cfg = LocResConfig(
                    halfmap1=half1,
                    halfmap2=half2,
                    mask=getattr(args, "mask", None),
                    local_fsc_threshold=float(getattr(args, "threshold", 0.143)),
                    sampling=getattr(args, "sampling", "auto"),
                    radius=getattr(args, "radius", "auto"),
                    edgewidth=getattr(args, "edgewidth", None),
                    cycles=float(getattr(args, "cycles", 10.0)),
                    gamma=float(getattr(args, "gamma", 1.8)),
                    bins=getattr(args, "bins", None),
                    rand_global=getattr(args, "rand_global", "auto"),
                    rand_local=getattr(args, "rand_local", "auto"),
                    rand_frac_nyq=getattr(args, "rand_frac_nyq", None),
                    cpu=int(getattr(args, "cpu", 1)),
                    resample=resample_flag,
                    interp=getattr(args, "interp", "cubic"),
                    output_basename=out_base,
                    plot=bool(getattr(args, "plot", False)),
                    # NEW: pass reference FSC profile (may be None)
                    reference_fsc_freqs=ref_fsc_freqs,
                    reference_fsc_vals=ref_fsc_vals,
                    output_mask=getattr(args, "output_mask", None),
                    compute_only_mask=bool(getattr(args, "output_mask", None) is not None),
                    mask_dilation_vox=(0 if getattr(args, "output_mask", None) is not None else 6),
                    fast=bool(getattr(args, "fast", False)),
                )
                res = locres_map(cfg)

                params = res.get("params", {})
                global_locres_params = {
                    "sampling": float(params.get("sampling_A")),
                    "radius": float(params.get("radius_A")),
                    "first_fsc_crossing": float(
                        res.get("global_fsc", {})
                        .get("crossings", {})
                        .get(AUX_THR, np.nan)
                    ),
                    "mask": getattr(args, "mask", None),
                }
                if not np.isfinite(global_locres_params["first_fsc_crossing"]):
                    print(
                        "[locresBulk] Warning: reference run did not yield a finite "
                        f"{AUX_THR} FSC crossing; reuse will still use sampling/radius."
                    )

                # One-shot resampled mask from the reference run, if requested
                if resample_flag and resample_mask_out:
                    write_resampled_mask_from_locres_result(
                        locres_result=res,
                        like_basename=out_base,
                        out_mask_path=resample_mask_out,
                    )
                if not skip_locinfo:
                    cfg_info = replace(cfg, output_basename=out_base)
                    _ = locinfo_map(cfg_info)

            else:
                # ---------------- Reuse mode: parameters from reference ----------------
                cfg = LocResConfig(
                    halfmap1=half1,
                    halfmap2=half2,
                    mask=global_locres_params["mask"],
                    local_fsc_threshold=float(getattr(args, "threshold", 0.143)),
                    sampling=getattr(args, "sampling", "auto"),
                    radius=getattr(args, "radius", "auto"),
                    edgewidth=getattr(args, "edgewidth", None),
                    cycles=float(getattr(args, "cycles", 10.0)),
                    gamma=float(getattr(args, "gamma", 1.8)),
                    bins=getattr(args, "bins", None),
                    rand_global=getattr(args, "rand_global", "auto"),
                    rand_local=getattr(args, "rand_local", "auto"),
                    rand_frac_nyq=getattr(args, "rand_frac_nyq", None),
                    cpu=int(getattr(args, "cpu", 1)),
                    resample=resample_flag,
                    interp=getattr(args, "interp", "cubic"),
                    output_basename=out_base,
                    plot=bool(getattr(args, "plot", False)),
                    precomputed_sampling=global_locres_params["sampling"],
                    precomputed_radius=global_locres_params["radius"],
                    precomputed_first_fsc_crossing=global_locres_params[
                        "first_fsc_crossing"
                    ],
                    reference_fsc_freqs=ref_fsc_freqs,
                    reference_fsc_vals=ref_fsc_vals,
                    output_mask=getattr(args, "output_mask", None),
                    compute_only_mask=bool(getattr(args, "output_mask", None) is not None),
                    mask_dilation_vox=(0 if getattr(args, "output_mask", None) is not None else 6),
                    fast=bool(getattr(args, "fast", False)),
                )
                _ = locres_map(cfg)
                if not skip_locinfo:
                    cfg_info = replace(cfg, output_basename=out_base)
                    _ = locinfo_map(cfg_info)

            any_done = True

        except Exception as e:
            print(f"[locresBulk] ERROR (label {label_str}): {e}")

    if not any_done:
        print("[locresBulk] No valid half-map pairs found; nothing was computed.")



def _mean_min_max_for_locres_single(
    locres_file: str,
    mask_file: str,
    mask_threshold: float = 0.2,
) -> Tuple[float, float, float, float, float]:
    """
    Compute robust min / quartiles / mean / max for a single locres map
    within a mask, following the behaviour of janas_app_meanMinMax:

      - select voxels where mask > mask_threshold (default 0.2)
      - compute mean over selected voxels
      - sort values, get Q1=floor(N/4), Q3=floor(3N/4)
      - minFinal = values[2], maxFinal = values[-3] if N>4 else true min/max
      - then override maxFinal = values[floor(N * 0.925)]
      - return (minFinal, Q1, mean, Q3, maxFinal)

    This is intended to match the logic of the standalone C++ tool.
    """
    if not os.path.isfile(locres_file):
        raise FileNotFoundError(f'map file not found: "{locres_file}"')
    if not os.path.isfile(mask_file):
        raise FileNotFoundError(f'mask file not found: "{mask_file}"')

    locres_vol, _ = IO_utils.read_mrc_data(locres_file)
    mask_vol, _   = IO_utils.read_mrc_data(mask_file)

    if locres_vol.shape != mask_vol.shape:
        raise ValueError(
            f"[locresStats] shape mismatch: map {locres_vol.shape} vs mask {mask_vol.shape}"
        )

    I = np.asarray(locres_vol, dtype=np.float64)
    M = np.asarray(mask_vol, dtype=np.float64)

    sel = I[M > mask_threshold]
    N = sel.size
    if N == 0:
        raise ValueError("[locresStats] no voxels selected with mask > 0.2")

    mean_val = float(sel.mean())

    values = np.sort(sel)  # ascending
    q1_idx = int(N / 4.0)
    q3_idx = int(N * 3.0 / 4.0)
    q1_val = float(values[q1_idx])
    q3_val = float(values[q3_idx])

    if N > 4:
        min_final = float(values[2])
        max_final = float(values[-3])
    else:
        min_final = float(values[0])
        max_final = float(values[-1])

    # override max_final to index floor(N * 0.925), clamped
    idx_925 = int(N * (3.70 / 4.0))  # 0.925
    if idx_925 >= N:
        idx_925 = N - 1
    if idx_925 < 0:
        idx_925 = 0
    max_final = float(values[idx_925])

    # Return in the same order as the standalone helper:
    # minFinal, q1, mean, q3, maxFinal
    return min_final, q1_val, mean_val, q3_val, max_final

def _median_min_max_for_locres_single(
    locres_file: str,
    mask_file: str,
    mask_threshold: float = 0.2,
) -> Tuple[float, float, float, float, float]:
    """
    Compute robust min / quartiles / median / max for a single locres map
    within a mask, mirroring _mean_min_max_for_locres_single but using the
    median instead of the mean as the central statistic.

      - select voxels where mask > mask_threshold (default 0.2)
      - compute median over selected voxels
      - sort values, get Q1=floor(N/4), Q3=floor(3N/4)
      - minFinal = values[2], maxFinal = values[-3] if N>4 else true min/max
      - then override maxFinal = values[floor(N * 0.925)]
      - return (minFinal, Q1, median, Q3, maxFinal)
    """
    if not os.path.isfile(locres_file):
        raise FileNotFoundError(f'map file not found: "{locres_file}"')
    if not os.path.isfile(mask_file):
        raise FileNotFoundError(f'mask file not found: "{mask_file}"')

    locres_vol, _ = IO_utils.read_mrc_data(locres_file)
    mask_vol, _   = IO_utils.read_mrc_data(mask_file)

    if locres_vol.shape != mask_vol.shape:
        raise ValueError(
            f"[locresStats] shape mismatch: map {locres_vol.shape} vs mask {mask_vol.shape}"
        )

    I = np.asarray(locres_vol, dtype=np.float64)
    M = np.asarray(mask_vol, dtype=np.float64)

    sel = I[M > mask_threshold]
    N = sel.size
    if N == 0:
        raise ValueError("[locresStats] no voxels selected with mask > 0.2")

    median_val = float(np.median(sel))

    values = np.sort(sel)  # ascending
    q1_idx = int(N / 4.0)
    q3_idx = int(N * 3.0 / 4.0)
    q1_val = float(values[q1_idx])
    q3_val = float(values[q3_idx])

    if N > 4:
        min_final = float(values[2])
        max_final = float(values[-3])
    else:
        min_final = float(values[0])
        max_final = float(values[-1])

    # override max_final to index floor(N * 0.925), clamped
    idx_925 = int(N * (3.70 / 4.0))  # 0.925
    if idx_925 >= N:
        idx_925 = N - 1
    if idx_925 < 0:
        idx_925 = 0
    max_final = float(values[idx_925])

    return min_final, q1_val, median_val, q3_val, max_final


def _label_from_locres_filename(path: Path) -> str:
    """
    Extract the numeric label (usually particle count) from a map filename.

    Supports:
      - ..._<NUM>_locres.mrc
      - ..._<NUM>_locinfo.mrc
      - ..._best<NUM>_locres.mrc
      - ..._best<NUM>_locinfo.mrc
      - fallback: last run of digits before '.mrc'
    """
    name = path.name

    # Preferred pattern: digits immediately before '_locres.mrc' or '_locinfo.mrc'
    # e.g. subset_1_105340_locres.mrc → "105340"
    #      subset_1_105340_locinfo.mrc → "105340"
    m = re.search(r"(\d+)(?=_(?:locres|locinfo)\.mrc$)", name)
    if m:
        return m.group(1)

    # Backwards compatibility: ..._best<NUM>_locres.mrc or ..._best<NUM>_locinfo.mrc
    m = re.search(r"best(\d+)_(?:locres|locinfo)\.mrc$", name)
    if m:
        return m.group(1)

    # Generic fallback: last run of digits before '.mrc'
    m = re.search(r"(\d+)(?=\.mrc$)", name)
    if m:
        return m.group(1)

    # Last resort: use the stem
    return path.stem



def produce_LocalMinBestParticles_map(
    locres_files: List[str],
    mask_path: Optional[str],
    out_mrc: str,
) -> None:
    """
    Build a voxel-wise map of the particle count that achieves the best
    (numerically smallest) local resolution at each voxel, analogous to
    'stats_minParticles.mrc' from _aggregate_locres_set().

    Rules:
      - Only locres values that are finite and >= 0.5 Å are considered valid.
      - For each voxel inside the mask, take the particle count corresponding
        to the smallest valid locres across all maps.
      - Voxels outside the mask, and voxels inside the mask where no valid
        locres >= 0.5 Å exists, are set to the maximum particle count observed
        across the input maps.

    Parameters
    ----------
    locres_files : list of str
        Paths to locres MRC files. All maps must have identical shapes.
    mask_path : str or None
        Path to a mask MRC. Voxels with mask > 0.2 are considered inside the mask.
        If None, all voxels are treated as inside the mask.
    out_mrc : str
        Output filename for the min-particle-count map.
    """
    if not locres_files:
        print("[locresStats] No locres files given for LocalMinBestParticles map; nothing to do.")
        return

    # Deduplicate, check existence, and sort for reproducibility
    unique_paths: List[Path] = []
    seen = set()
    for f in locres_files:
        p = Path(f)
        if str(p) in seen:
            continue
        seen.add(str(p))
        if not p.is_file():
            print(f"[locresStats] WARNING: locres map not found for LocalMinBestParticles, skipping: '{p}'")
            continue
        unique_paths.append(p)

    if not unique_paths:
        print("[locresStats] No existing locres files for LocalMinBestParticles; nothing written.")
        return

    unique_paths = sorted(unique_paths, key=lambda p: str(p))

    # First pass: determine which files have a valid integer label and find max particle count
    labelled_entries: List[Tuple[Path, int]] = []
    max_particles: Optional[int] = None

    for p in unique_paths:
        label = _label_from_locres_filename(p)
        try:
            count_val = int(label)
        except Exception:
            print(
                f"[locresStats] WARNING: cannot interpret label '{label}' in "
                f"{p.name} as integer particle count; skipping in LocalMinBestParticles map."
            )
            continue

        labelled_entries.append((p, count_val))
        if max_particles is None:
            max_particles = count_val
        else:
            if count_val > max_particles:
                max_particles = count_val

    if not labelled_entries:
        print("[locresStats] No locres files with integer particle labels; LocalMinBestParticles map not written.")
        return

    assert max_particles is not None
    max_particles_f = float(max_particles)

    # Use the first labelled locres file as template for header and shape
    first_vol, first_hdr = IO_utils.read_mrc_data(str(labelled_entries[0][0]))
    tmpl_shape = first_vol.shape

    # Load mask and build boolean mask of voxels to consider "inside"
    if mask_path:
        mask_vol, _mh = IO_utils.read_mrc_data(mask_path)
        if mask_vol.shape != tmpl_shape:
            raise ValueError(
                f"[locresStats] Mask shape mismatch in LocalMinBestParticles: "
                f"{mask_vol.shape} vs {tmpl_shape}"
            )
        mask_bool = (mask_vol > 0.2)
    else:
        mask_bool = np.ones(tmpl_shape, dtype=bool)

    # Accumulators: best locres per voxel and corresponding particle count.
    # Initialise best_count to max_particles everywhere, so outside-mask voxels
    # and voxels with no valid locres will naturally end up with max_particles.
    best_map = np.full(tmpl_shape, np.inf, dtype=np.float32)
    best_count = np.full(tmpl_shape, max_particles_f, dtype=np.float32)

    RES_MIN_VALID = 0.5  # Å: ignore locres values below this threshold

    for p, count_val in labelled_entries:
        locres_vol, _ = IO_utils.read_mrc_data(str(p))
        if locres_vol.shape != tmpl_shape:
            raise ValueError(
                f"[locresStats] Shape mismatch in {p.name} for LocalMinBestParticles: "
                f"{locres_vol.shape} vs {tmpl_shape}"
            )

        # Valid voxels: finite and >= RES_MIN_VALID
        valid = np.isfinite(locres_vol) & (locres_vol >= RES_MIN_VALID)

        # Only update voxels inside the mask
        valid_inside = mask_bool & valid

        # Per-voxel update: where this locres is valid, inside mask, and better than the current best
        better = valid_inside & (locres_vol < best_map)
        if np.any(better):
            best_map[better] = locres_vol[better]
            best_count[better] = float(count_val)

    # At this point:
    #  - outside mask: best_count == max_particles_f (never updated)
    #  - inside mask, no valid locres >= 0.5: best_count still == max_particles_f
    #  - inside mask, at least one valid locres: best_count == particle count of best locres

    # If absolutely nothing was updated (everything stayed at +inf), we still have a map
    # filled with max_particles. That is consistent with the stated rule.
    IO_utils.write_mrc_like(
        out_mrc,
        best_count.astype(np.float32),
        first_hdr,
        update_stats=True,
    )
    print(f"[locresStats] Wrote LocalMinBestParticles map '{out_mrc}'")


def run_locres_stats(
    locres_files: List[str],
    mask_path: str,
    out_csv: str = "bestRanked_locres_values.csv",
    out_LocalMinBestParticles_map: str = "",
    assessmentMethod: str = "mean",
) -> None:
    """
    Compute masked summary statistics for a collection of local-resolution maps.

    Parameters
    ----------
    locres_files : list of str
        Paths to locres MRC maps. These should already be resolved (no wildcards).
    mask_path : str
        Path to a mask MRC; voxels with mask > 0.2 are included.
    out_csv : str
        Output CSV filename. The file will have one header line:

            numParticles,max,highQuartile,<mean|median>,lowQuartile,min

        followed by one line per locres map.
    assessmentMethod : str
        "mean" or "median". Controls which statistic is stored in the central column.
    """
    if not locres_files:
        print("[locresStats] No locres files provided; nothing to do.")
        return
    print ("assessmentMethod=", assessmentMethod)

    mask_path = str(mask_path)
    if not os.path.isfile(mask_path):
        raise FileNotFoundError(f"[locresStats] mask file not found: '{mask_path}'")

    method = str(assessmentMethod).strip().lower()
    if method not in ("mean", "median"):
        raise ValueError(
            f"[locresStats] assessmentMethod must be 'mean' or 'median', got: '{assessmentMethod}'"
        )

    # Pick the single-map summariser
    if method == "median":
        summarize_fn = _median_min_max_for_locres_single
        center_col_name = "median"
    else:
        summarize_fn = _mean_min_max_for_locres_single
        center_col_name = "mean"

    # Deduplicate and sort filenames for reproducible output
    unique_paths: List[Path] = []
    seen = set()
    for f in locres_files:
        p = Path(f)
        if str(p) in seen:
            continue
        seen.add(str(p))
        if not p.is_file():
            print(f"[locresStats] WARNING: map file not found, skipping: '{p}'")
            continue
        unique_paths.append(p)

    if not unique_paths:
        print("[locresStats] No existing locres files after filtering; nothing written.")
        return

    # Sort lexicographically first (for reproducibility of tie-breaks)
    unique_paths = sorted(unique_paths, key=lambda p: str(p))

    # Collect results as structured records
    # each: (sort_key, label_str, min_final, q1, center, q3, max_final)
    records = []

    for p in unique_paths:
        try:
            min_final, q1_val, center_val, q3_val, max_final = summarize_fn(
                locres_file=str(p),
                mask_file=mask_path,
                mask_threshold=0.2,
            )
        except Exception as e:
            print(f"[locresStats] ERROR while processing {p.name}: {e}")
            continue

        label = _label_from_locres_filename(p)

        # try to interpret the label as an integer particle count
        try:
            sort_key = int(label)
        except Exception:
            # if it is not an integer, push it to the end
            sort_key = float("inf")

        records.append(
            (
                sort_key,
                str(label),
                float(min_final),
                float(q1_val),
                float(center_val),
                float(q3_val),
                float(max_final),
            )
        )

    if not records:
        print("[locresStats] All files failed; no CSV written.")
        return

    # Sort by numeric particle count (ascending)
    records.sort(key=lambda r: r[0])

    # Build CSV rows
    header = f"numParticles,max,highQuartile,{center_col_name},lowQuartile,min"
    rows: List[str] = [header]

    for _, label_str, min_final, q1_val, center_val, q3_val, max_final in records:
        row = ",".join(
            [
                label_str,
                f"{min_final:.5f}",
                f"{q1_val:.5f}",
                f"{center_val:.5f}",
                f"{q3_val:.5f}",
                f"{max_final:.5f}",
            ]
        )
        rows.append(row)

    out_path = Path(out_csv)
    try:
        with out_path.open("w", encoding="utf-8") as f:
            f.write("\n".join(rows) + "\n")
    except Exception as e:
        raise RuntimeError(f"[locresStats] cannot write CSV '{out_csv}': {e}")

    print(f"[locresStats] Wrote {out_path}")

    # Optionally produce a stats_minParticles-style map
    if out_LocalMinBestParticles_map:
        try:
            produce_LocalMinBestParticles_map(
                locres_files=[str(p) for p in unique_paths],
                mask_path=mask_path,
                out_mrc=out_LocalMinBestParticles_map,
            )
        except Exception as e:
            print(
                f"[locresStats] ERROR while producing LocalMinBestParticles map "
                f"'{out_LocalMinBestParticles_map}': {e}"
            )
