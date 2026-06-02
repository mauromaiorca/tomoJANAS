#!/usr/bin/env python3
# reconstructor_utils.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
Utilities for JANAS reconstruction and local-resolution mapping.

Public API:
- Reconstruction: reconstruct_3D, post_process
- LocRes: LocResConfig, locres_map

- FSC helpers: compute_fsc, first_FSC_crossing, find_FSC_resolutions_and_stats
"""

from __future__ import annotations

import sys
import os
import fcntl
import re
import struct
import argparse
import time
import math
import shutil
import tempfile
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any, Union
from dataclasses import dataclass, field

# third-party
import numpy as np
from numpy.fft import fftn, ifftn, fftfreq
from scipy.spatial.transform import Rotation as R
from scipy.ndimage import gaussian_filter, binary_dilation, zoom as ndi_zoom, distance_transform_edt
from scipy.interpolate import interp1d

# Optional plotting (used only if requested by caller)
import matplotlib.pyplot as plt

# External deps used by reconstruction
import janas.janas_core as janas_core
from janas import assessParticles, starHandler, IO_utils  

__all__ = [
    # reconstruction
    "reconstruct_3D", "post_process",
]

# ----------------------- (Optional) Torch -----------------------
try:
    import torch
    import torch.nn.functional as F
    TORCH_AVAILABLE = True
except Exception:
    TORCH_AVAILABLE = False


# ----------------------- helper -----------------------
def error_log_path_from_out_mrc(out_mrc_path: str) -> str:
    # testUpdateRec_recH1.mrc -> testUpdateRec_recH1_recError.log
    base, _ = os.path.splitext(out_mrc_path)
    return base + "_recError.log"

def append_error_rows_csv(log_path: str, rows):
    """
    rows: list of tuples (imgname, "NaN")
    Writes CSV with header: _rlnImageName,ERROR
    Safe for concurrent writers via fcntl lock.
    """
    if not log_path or not rows:
        return

    parent = os.path.dirname(log_path)
    if parent:
        os.makedirs(parent, exist_ok=True)

    with open(log_path, "a+", encoding="utf-8", newline="") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.seek(0, os.SEEK_END)
        if f.tell() == 0:
            f.write("_rlnImageName,ERROR\n")
        for imgname, err in rows:
            f.write(f"{imgname},{err}\n")
        f.flush()
        os.fsync(f.fileno())
        fcntl.flock(f, fcntl.LOCK_UN)

def _error_log_path_from_checkpoint(checkpoint_out_path: str) -> str:
    """
    From e.g. 'testUpdateRec_recH1.mrc' -> 'testUpdateRec_recError.log'
          e.g. '.../testUpdateRec_recH2.mrc' -> '.../testUpdateRec_recError.log'
    """
    base, _ext = os.path.splitext(checkpoint_out_path)
    base = base.replace("_recH1", "").replace("_recH2", "")
    return base + "_recError.log"


def append_rec_error_log(log_path: str, rows):
    """
    Append CSV rows to log_path in format:
    _rlnImageName,ERROR
    """
    if not log_path or not rows:
        return

    parent = os.path.dirname(log_path)
    if parent:
        os.makedirs(parent, exist_ok=True)

    import fcntl  # Linux
    with open(log_path, "a+", encoding="utf-8", newline="") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.seek(0, os.SEEK_END)
        if f.tell() == 0:
            f.write("_rlnImageName,ERROR\n")
        for imgname, err in rows:
            f.write(f"{imgname},{err}\n")
        f.flush()
        os.fsync(f.fileno())
        fcntl.flock(f, fcntl.LOCK_UN)


def _insert_count_in_filename(path: str, count_label: int) -> str:
    """
    Insert the count just before '_recH' (if present) or before '_rec' (if present),
    otherwise before the final extension.
    """
    p = Path(path); name = p.name
    if "_recH" in name:
        name = name.replace("_recH", f"{count_label}_recH")
        return str(p.with_name(name))
    if "_rec" in name:
        name = name.replace("_rec", f"{count_label}_rec")
        return str(p.with_name(name))
    # fallback
    return str(p.with_name(f"{p.stem}{count_label}{p.suffix}"))


# ======================= GPU helpers =======================

def compute_ctf_kernel_gpu(
    nx: int, ny: int, angpix: float,
    Voltage: float, DefocusU: float, DefocusV: float,
    DefocusAngle: float, SphericalAberration: float,
    CtfBfactor: float, PhaseShift: float, AmplitudeContrast: float,
    device: str = "cuda"
):
    import torch, numpy as np
    # --- Reject physically invalid ranges (finite but dangerous) ---
    if (not math.isfinite(Voltage)) or (Voltage <= 0.0):
        return torch.ones((ny, nx), device=device, dtype=torch.float32)
    if (not math.isfinite(angpix)) or (angpix <= 0.0):
        return torch.ones((ny, nx), device=device, dtype=torch.float32)
    if (not math.isfinite(CtfBfactor)):
        return torch.ones((ny, nx), device=device, dtype=torch.float32)
    if CtfBfactor < 0.0:
        CtfBfactor = 0.0
    if (not math.isfinite(AmplitudeContrast)):
        return torch.ones((ny, nx), device=device, dtype=torch.float32)


    pi = np.pi
    local_Cs = SphericalAberration * 1e7
    local_kV = Voltage * 1e3
    lam = 12.2643247 / torch.sqrt(
        torch.tensor(local_kV, device=device) *
        (1.0 + torch.tensor(local_kV, device=device) * 0.978466e-6)
    )
    K1 = torch.tensor(pi, device=device) * lam
    K2 = (torch.tensor(pi, device=device) / 2) * torch.tensor(local_Cs, device=device) * lam**3
    Ac = torch.tensor(float(AmplitudeContrast), device=device)
    Ac = torch.nan_to_num(Ac, nan=0.1, posinf=0.1, neginf=0.1)
    Ac = torch.clamp(Ac, 0.0, 0.999999)
    K3 = torch.atan(Ac / torch.sqrt(1 - Ac*Ac))
    K4 = -CtfBfactor / 4.0
    K5 = PhaseShift * pi / 180.0

    rad_az = DefocusAngle * pi / 180.0
    def_avg = -(DefocusU + DefocusV) * 0.5
    def_dev = -(DefocusU - DefocusV) * 0.5

    ys, xs = ny*angpix, nx*angpix
    yy = torch.arange(ny, device=device) - (ny//2)
    xx = torch.arange(nx, device=device) - (nx//2)
    Y, X = torch.meshgrid(yy, xx, indexing="ij")
    fx = X / xs; fy = Y / ys
    u2 = fx**2 + fy**2
    ell_ang = torch.atan2(fy, fx) - rad_az
    deltaf = def_avg + def_dev * torch.cos(2*ell_ang)
    arg = K1 * deltaf * u2 + K2 * (u2**2) - K5 - K3
    return -torch.sin(arg) * torch.exp(K4 * u2)


# ================================================================================
# ======================= Core Reconstruction (single job) =======================
def compute_ctf_kernel_cpu(
    nx: int, ny: int, angpix: float,
    Voltage: float, DefocusU: float, DefocusV: float,
    DefocusAngle: float, SphericalAberration: float,
    CtfBfactor: float, PhaseShift: float, AmplitudeContrast: float,
):
    # --- Sanity: avoid NaN/Inf in CTF parameters (NaN would poison the whole kernel) ---
    if (not np.isfinite(Voltage)) or (Voltage <= 0.0):
        return np.ones((ny, nx), dtype=np.float32)
    if (not np.isfinite(angpix)) or (angpix <= 0.0):
        return np.ones((ny, nx), dtype=np.float32)
    if (not np.isfinite(CtfBfactor)):
        return np.ones((ny, nx), dtype=np.float32)
    if CtfBfactor < 0.0:
        CtfBfactor = 0.0
    if (not np.isfinite(AmplitudeContrast)):
        return np.ones((ny, nx), dtype=np.float32)

    # If amplitude contrast is finite but outside [0,1), clamp; if it is NaN we already returned
    AmplitudeContrast = float(np.clip(AmplitudeContrast, 0.0, 0.999999))
    PhaseShift = float(PhaseShift)
    CtfBfactor = float(CtfBfactor)

    pi = np.pi
    CsA  = SphericalAberration * 1e7
    kV   = Voltage * 1e3
    lam  = 12.2643247 / np.sqrt(kV * (1.0 + kV * 0.978466e-6))
    K1   = pi * lam
    K2   = (pi/2.0) * CsA * (lam**3)
    Ac   = np.clip(AmplitudeContrast, 0.0, 0.999999)
    K3   = np.arctan(Ac / np.sqrt(1 - Ac*Ac))
    K4   = -CtfBfactor / 4.0
    K5   = PhaseShift * pi / 180.0

    rad_az  = DefocusAngle * pi / 180.0
    def_avg = -(DefocusU + DefocusV) * 0.5
    def_dev = -(DefocusU - DefocusV) * 0.5

    ys, xs = ny*angpix, nx*angpix
    yy = np.arange(ny, dtype=np.float32) - (ny//2)
    xx = np.arange(nx, dtype=np.float32) - (nx//2)
    Y, X = np.meshgrid(yy, xx, indexing="ij")
    fx = X / xs; fy = Y / ys
    u2 = fx*fx + fy*fy
    ell_ang = np.arctan2(fy, fx) - rad_az
    deltaf = def_avg + def_dev * np.cos(2.0*ell_ang)
    arg = K1 * deltaf * u2 + K2 * (u2**2) - K5 - K3
    return -np.sin(arg) * np.exp(K4 * u2)


def reconstruct_3D(
    image_entries: List[Tuple[str, int]],
    params:       List[Tuple[float, float, float, float, float]],
    pixel_size:   float,
    box_size:     int,
    acc: Optional[np.ndarray] = None,
    weight: Optional[np.ndarray] = None,
    use_ctf: bool = False,
    ctf_params: Optional[List[Tuple[float, float, float, float, float, float, float, float, float]]] = None,
    ctf_mode: Optional[str] = None,   # "modulate" | "phaseflip" | "wiener"
    wiener_tau: float = 0.1,
    wiener_tau_csv: Optional[str] = None,
    device: str = "cpu",
    cpu_workers: int = 1,        # number of CPU workers if device=="cpu"
    gpu_index: Optional[int] = None,
    gpu_batch: int = 20,         # mini-batch on GPU

    # Progressive checkpoints
    checkpoint_counts: Optional[List[int]] = None,
    checkpoint_out_path: Optional[str] = None,
    lp_start: float = 0.9,
    lp_rolloff: Optional[int] = None,
    lp_outer_scale: float = 1.1,
    lp_add_px: Optional[int] = None,

    # NEW: labels to use in filenames (global targets)
    checkpoint_labels: Optional[List[int]] = None,
    on_checkpoint_kspace=None
) -> Tuple[np.ndarray, np.ndarray]:

    N = box_size
    halfN = N // 2
    freq_unit = 1.0 / (N * pixel_size)
    two_pi = 2.0 * np.pi

    # ---- Optional: build a per-shell τ(k) 2D image once (from SSNR CSV) ----
    def _load_ssnr_csv_to_tau_curve(csv_path: str):
        arr = np.loadtxt(csv_path, delimiter=",", comments="#")
        # Accept "f,SSNR" or arbitrary two columns; robustly pick the first two numeric cols
        if arr.ndim == 1:
            arr = arr[None, :]
        freq = arr[:, 0].astype(np.float32)
        ssnr = arr[:, 1].astype(np.float32)
        freq = np.clip(freq, 0.0, None)
        ssnr = np.clip(ssnr, 0.0, np.finfo(np.float32).max)
        # τ = 1/SSNR; guard SSNR ~ 0
        tau = np.where(ssnr > 1e-6, 1.0/ssnr, np.finfo(np.float32).max)
        # Also cap τ for extreme lows at high frequency (numerical safety)
        p99 = np.nanpercentile(tau[np.isfinite(tau)], 99.9) if np.any(np.isfinite(tau)) else 1e3
        tau = np.clip(tau, 0.0, float(p99))
        return freq, tau

    tau_image_np = None
    tau_image_t  = None
    if (ctf_mode == "wiener") and (wiener_tau_csv is not None):
        # Radial freq grid in 2D, normalised to 1/Å
        yy = np.arange(N, dtype=np.float32) - (N//2)
        xx = np.arange(N, dtype=np.float32) - (N//2)
        Y, X = np.meshgrid(yy, xx, indexing="ij")
        fx = X / (N * pixel_size)
        fy = Y / (N * pixel_size)
        fr = np.sqrt(fx*fx + fy*fy).astype(np.float32)  # [1/Å]
        fr = np.fft.ifftshift(fr)  # match our kernel shift convention
        f_curve, tau_curve = _load_ssnr_csv_to_tau_curve(wiener_tau_csv)
        # Interp τ on absolute frequency (1/Å). Clamp outside range to edge values.
        fmin, fmax = float(np.min(f_curve)), float(np.max(f_curve))
        tau_interp = np.interp(np.clip(fr, fmin, fmax), f_curve, tau_curve).astype(np.float32)
        tau_image_np = tau_interp  # store as IFFT-shifted to match our kernels below
        if TORCH_AVAILABLE and device != "cpu":
            tau_image_t = torch.from_numpy(tau_image_np).to(device)


    # Prepare progressive checkpoints
    counts = []
    if checkpoint_counts:
        counts = sorted({int(c) for c in checkpoint_counts if int(c) > 0})
    next_ci = 0

    # Map write labels (global) vs trigger counts (local per-job)
    labels = None
    if checkpoint_labels and counts and len(checkpoint_labels) >= len(counts):
        labels = [int(x) for x in checkpoint_labels[:len(counts)]]

    def _maybe_checkpoint_cpu(processed: int, acc_ref: np.ndarray, w_ref: np.ndarray):
        nonlocal next_ci
        while next_ci < len(counts) and processed >= counts[next_ci]:
            label = labels[next_ci] if labels is not None else counts[next_ci]
            outp = _insert_count_in_filename(checkpoint_out_path, label) if checkpoint_out_path else None
            if outp:
                vol = post_process(
                    acc_ref.copy(), w_ref.copy(),
                    pixel_size=pixel_size,
                    device="cpu",
                    lp_start=lp_start,
                    lp_rolloff=lp_rolloff,
                    lp_outer_scale=lp_outer_scale,
                    lp_add_px=lp_add_px,
                )
                IO_utils.write_mrc(outp, vol, pixel_size)
                print(f"[checkpoint] wrote {outp} ")
            next_ci += 1

    def _maybe_checkpoint_gpu(processed: int, acc_t, w_t):
        nonlocal next_ci
        while next_ci < len(counts) and processed >= counts[next_ci]:
            label = labels[next_ci] if labels is not None else counts[next_ci]
            outp = _insert_count_in_filename(checkpoint_out_path, label) if checkpoint_out_path else None
            if outp:
                vol = post_process(
                    acc_t.detach().cpu().numpy().copy(),
                    w_t.detach().cpu().numpy().copy(),
                    pixel_size=pixel_size,
                    device="cpu",
                    lp_start=lp_start,
                    lp_rolloff=lp_rolloff,
                    lp_outer_scale=lp_outer_scale,
                    lp_add_px=lp_add_px,
                )
                IO_utils.write_mrc(outp, vol, pixel_size)
                print(f"[checkpoint] wrote {outp} ")
            next_ci += 1

    # --- Robust slice reader (shared by CPU and GPU paths) ---
    # A single non-finite pixel (NaN/Inf) in one particle poisons the FFT and can turn the entire map into NaN.
    bad_reads = 0
    def _read_slice(stack_path: str, sl_in: int) -> np.ndarray:
        nonlocal bad_reads
        # Try a small set of candidate indices (handles -1, off-by-one, etc.)
        cand = [int(sl_in), int(sl_in) - 1, int(sl_in) + 1]
        # de-duplicate while preserving order; drop negatives
        cand = [c for c in dict.fromkeys(cand) if c >= 0]

        last = None
        for c in cand:
            try:
                arr = np.asarray(janas_core.ReadMrcSlice(stack_path, c), dtype=np.float32).reshape(N, N)
            except Exception:
                continue
            if np.isfinite(arr).all():
                return arr
            last = arr

        bad_reads += 1
        if bad_reads <= 10:
            print(f"[WARN] Non-finite or unreadable slice: {stack_path} (slice={sl_in}); replacing non-finite with 0")
        if last is None:
            return np.zeros((N, N), dtype=np.float32)
        return np.nan_to_num(last, nan=0.0, posinf=0.0, neginf=0.0)

    error_rows = []  # list of (imgname, "NaN")
    def _log_nan_particle(global_i: int):
        stack_i, sl_i = image_entries[global_i]
        imgname = f"{sl_i}@{stack_i}"  # formato _rlnImageName
        error_rows.append((imgname, "NaN"))

    # --------------- CPU PATH ---------------
    if device == "cpu":
        total = len(image_entries)
        if total == 0:
            return (np.zeros((N,N,N), np.__dict__['complex64']), np.zeros((N,N,N), np.float32))

        # Progressive checkpoints require single-process on CPU
        if counts and cpu_workers > 1:
            print("NOTE: --subrec checkpoints requested on CPU; forcing single-process to enable incremental writes.")
            cpu_workers = 1

        # --- Single-process path ---
        if cpu_workers <= 1:
            acc    = np.zeros((N, N, N), dtype=np.__dict__['complex64'])
            weight = np.zeros((N, N, N), dtype=np.float32)
            int_freq = np.arange(N, dtype=np.int32)
            int_freq[int_freq > halfN] -= N
            KX_cpu, KY_cpu = np.meshgrid(int_freq, int_freq, indexing="xy")
            FX_cpu = KX_cpu * freq_unit
            FY_cpu = KY_cpu * freq_unit

            start_time = time.time()
            for idx, (entry, par) in enumerate(zip(image_entries, params), start=1):
                stack, sl = entry
                rot, tilt, psi, ox, oy = par
                Rmat = R.from_euler("ZYZ", [rot, tilt, psi], degrees=True).as_matrix()


                #raw = np.asarray(janas_core.ReadMrcSlice(stack, sl), dtype=np.float32).reshape(N, N)
                raw = _read_slice(stack, sl)
                raw = raw - raw.mean()
                F2  = np.fft.fft2(np.fft.ifftshift(raw))
                if use_ctf and ctf_params is not None:
                    V, Du, Dv, Da, Cs, Bf, Ph, Ac, _ = ctf_params[idx-1]
                    k = compute_ctf_kernel_cpu(N, N, pixel_size, V, Du, Dv, Da, Cs, Bf, Ph, Ac)
                    k = np.fft.ifftshift(k)
                    mode = (ctf_mode or "modulate")
                    if mode == "phaseflip":
                        F2 *= np.sign(k)
                    elif mode == "wiener":
                        if tau_image_np is not None:
                            # per-shell τ image (already ifftshifted to match k)
                            denom = (k*k + tau_image_np)
                            F2 *= (k / denom)
                        else:
                            F2 *= (k / (k*k + float(wiener_tau)))
                    else:  # "modulate"
                        F2 *= k

                #dx, dy = ox / pixel_size, oy / pixel_size
                dx = ox / pixel_size
                dy = oy / pixel_size
                if not np.isfinite(dx): dx = 0.0
                if not np.isfinite(dy): dy = 0.0
                if dx or dy:
                    phase = np.exp(-1j * two_pi * (KX_cpu * dx + KY_cpu * dy) / N)
                    F2 *= phase

                for ky in range(N):
                    fy = FY_cpu[ky, 0]
                    for kx in range(N):
                        val = F2[ky, kx]
                        if val == 0:
                            continue
                        fx = FX_cpu[ky, kx]
                        xp, yp, zp = Rmat.dot([fx, fy, 0.0])
                        if xp < 0:
                            xp, yp, zp = -xp, -yp, -zp
                            val = np.conj(val)
                        kxr = int(round(xp / freq_unit))
                        kyr = int(round(yp / freq_unit))
                        kzr = int(round(zp / freq_unit))
                        if abs(kxr) > halfN or abs(kyr) > halfN or abs(kzr) > halfN:
                            continue
                        ix, iy, iz = kxr % N, kyr % N, kzr % N
                        acc[iz, iy, ix]   += val
                        weight[iz, iy, ix] += 1.0

                if counts:
                    _maybe_checkpoint_cpu(idx, acc, weight)

                elapsed = time.time() - start_time
                avg = elapsed / idx
                rem = avg * (total - idx)
                sys.stdout.write(
                    f"\rProcessed particles {idx}/{total} ({100.0*idx/total:.1f}%) — "
                    f"elapsed {time.strftime('%H:%M:%S', time.gmtime(elapsed))} — "
                    f"ETA {time.strftime('%H:%M:%S', time.gmtime(rem))}"
                )
                sys.stdout.flush()
            sys.stdout.write("\n")
            print("CPU reconstruction (single-process) done.")
            if use_ctf:
                print("CTF modulation applied (CPU).")
            return acc, weight

        # --- Multiprocessing path (no progressive checkpoints) ---
        import multiprocessing as mp
        n_workers = min(int(cpu_workers), mp.cpu_count())
        sizes = [(total // n_workers) + (1 if i < (total % n_workers) else 0) for i in range(n_workers)]
        starts = np.cumsum([0] + sizes[:-1]).tolist()
        blocks = [(s, s + sz) for s, sz in zip(starts, sizes) if sz > 0]
        n_workers = len(blocks)

        gdict = dict(
            N=N,
            pixel_size=pixel_size,
            entries=image_entries,
            params=params,
            use_ctf=use_ctf,
            ctf_params=ctf_params
        )

        tmp_dir = tempfile.mkdtemp(prefix="recon_parts_")
        try:
            try:
                mp.set_start_method("fork")
            except RuntimeError:
                pass

            progress = mp.Value('i', 0)
            procs = []
            start = time.time()
            _init_globals(gdict)
            for pi, (b0, b1) in enumerate(blocks):
                p = mp.Process(
                    target=_worker_reconstruct,
                    args=(b0, b1, pi, tmp_dir, progress),
                )
                procs.append(p)

            for p in procs: p.start()

            last = -1
            while any(p.is_alive() for p in procs):
                done = progress.value
                if done != last:
                    elapsed = time.time() - start
                    rate = done / elapsed if elapsed > 0 else 0.0
                    remain = (total - done) / rate if rate > 0 else 0.0
                    pct = 100.0 * done / max(1, total)
                    sys.stdout.write(
                        f"\rProcessed particles {done:,}/{total:,} ({pct:5.1f}%) — "
                        f"{rate:,.1f} imgs/s — elapsed {time.strftime('%H:%M:%S', time.gmtime(elapsed))} — "
                        f"ETA {time.strftime('%H:%M:%S', time.gmtime(remain))}"
                    )
                    sys.stdout.flush()
                    last = done
                time.sleep(0.2)

            for p in procs: p.join()
            sys.stdout.write("\n")

            acc = np.zeros((N, N, N), dtype=np.__dict__['complex64'])
            weight = np.zeros((N, N, N), dtype=np.float32)
            for pi in range(n_workers):
                acc += np.load(os.path.join(tmp_dir, f"acc_{pi:04d}.npy"), allow_pickle=False)
                weight += np.load(os.path.join(tmp_dir, f"w_{pi:04d}.npy"), allow_pickle=False)

            total_time = time.time() - start
            print(f"CPU reconstruction (workers={n_workers}) done in {time.strftime('%H:%M:%S', time.gmtime(total_time))}")
            if use_ctf:
                print("CTF modulation applied (CPU).")
            return acc, weight

        finally:
            shutil.rmtree(tmp_dir, ignore_errors=True)

    # --------------- GPU PATH ---------------
    if not TORCH_AVAILABLE:
        raise RuntimeError("CUDA requested but PyTorch not available.")
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but no GPU is visible.")

    if gpu_index is not None:
        try:
            torch.cuda.set_device(int(gpu_index))
        except Exception:
            print(f"WARNING: could not set CUDA device {gpu_index}; using default.", file=sys.stderr)
    device_t = torch.device("cuda")
    acc_t    = torch.zeros((N, N, N), dtype=torch.complex64, device=device_t)
    weight_t = torch.zeros((N, N, N), dtype=torch.float32,   device=device_t)

    f = torch.arange(N, device=device_t)
    f[f > halfN] -= N
    KX = f.view(1, -1).repeat(N, 1)
    KY = f.view(-1, 1).repeat(1, N)
    FX = KX * freq_unit
    FY = KY * freq_unit
    two_pi_t = 2.0 * torch.pi

    Rmats_cpu = [
        R.from_euler("ZYZ", [r, t, p], degrees=True).as_matrix()
        for (r, t, p, ox, oy) in params
    ]
    Rmats_t = torch.from_numpy(np.stack(Rmats_cpu, axis=0)).to(device_t)

    total = len(image_entries)
    processed = 0
    start_time = time.time()

    bs = max(1, int(gpu_batch))

    # --- Skip bad particles (only triggers on errors; no effect on clean data) ---
    SKIP_BAD = True
    skipped = 0
    skipped_rows = []
    def _mark_skip(global_i: int, reason: str):
        nonlocal skipped
        skipped += 1
        stack_i, sl_i = image_entries[global_i]
        imgname = f"{sl_i}@{stack_i}"
        err = "NaN"
        skipped_rows.append((imgname, err))



    for b0 in range(0, total, bs):
        b1 = min(b0 + bs, total)
        M  = b1 - b0
        keep = [True] * M

        raw_list = []
        for i in range(b0, b1):
            j = i - b0
            stack_i, sl_i = image_entries[i]
            cpu_slice = janas_core.ReadMrcSlice(stack_i, sl_i)
            raw_np = np.asarray(cpu_slice, np.float32).reshape(N, N)
            raw_t  = torch.from_numpy(raw_np).to(device_t)
            m = raw_t.mean()  # already needed for normalisation
            if SKIP_BAD and (not torch.isfinite(m)):
                keep[j] = False
                _mark_skip(i, "raw_mean_nonfinite")
                raw_t = torch.zeros((N, N), dtype=torch.float32, device=device_t)
                m = raw_t.mean()
            raw_t = raw_t - m
            raw_list.append(raw_t)
        raw = torch.stack(raw_list, dim=0)


        img_c = torch.fft.ifftshift(raw, dim=(-2,-1))
        F2    = torch.fft.fft2(img_c)


        if use_ctf and (ctf_params is not None):
            mode = (ctf_mode or "modulate")
            tau  = float(wiener_tau)
            for j in range(M):
                if SKIP_BAD and (not keep[j]):
                    continue                
                V, Du, Dv, Da, Cs, Bf, Ph, Ac, _ = ctf_params[b0 + j]

                # Reject finite-but-invalid values that can yield NaNs inside CTF
                if SKIP_BAD:
                    if (not math.isfinite(V)) or (V <= 0.0):
                        keep[j] = False; _mark_skip(b0 + j, "ctf_voltage_invalid"); continue
                    if (not math.isfinite(Ac)):
                        keep[j] = False; _mark_skip(b0 + j, "ctf_ampcontrast_nonfinite"); continue
                    if (not math.isfinite(Bf)):
                        keep[j] = False; _mark_skip(b0 + j, "ctf_bfactor_nonfinite"); continue
                    if (not math.isfinite(Ph)):
                        keep[j] = False; _mark_skip(b0 + j, "ctf_phaseshift_nonfinite"); continue

                # clamp safe ranges (does not change sane data)
                if Bf < 0.0: Bf = 0.0
                if Ac < 0.0: Ac = 0.0
                if Ac >= 1.0: Ac = 0.999999                
                
                kernel = compute_ctf_kernel_gpu(
                    N, N, pixel_size, V, Du, Dv, Da, Cs, Bf, Ph, Ac, device=device_t
                )
                kernel = torch.nan_to_num(kernel, nan=0.0, posinf=0.0, neginf=0.0)
                kernel = torch.fft.ifftshift(kernel, dim=(0,1))
                if not torch.isfinite(kernel).all():
                    _log_nan_particle(b0 + j)
                    kernel = torch.nan_to_num(kernel, nan=0.0, posinf=0.0, neginf=0.0)

                if mode == "phaseflip":
                    F2[j] = F2[j] * torch.sign(kernel)
                elif mode == "wiener":
                    if tau_image_t is not None:
                        denom = (kernel*kernel + tau_image_t)
                        F2[j] = F2[j] * (kernel / denom)
                    else:
                        F2[j] = F2[j] * (kernel / (kernel*kernel + tau))
                else:  # "modulate"
                    F2[j] = F2[j] * kernel

        dx = torch.tensor([params[i][3] for i in range(b0,b1)], device=device_t).view(M,1,1) / pixel_size
        dy = torch.tensor([params[i][4] for i in range(b0,b1)], device=device_t).view(M,1,1) / pixel_size
        dx = torch.nan_to_num(dx, nan=0.0, posinf=0.0, neginf=0.0)
        dy = torch.nan_to_num(dy, nan=0.0, posinf=0.0, neginf=0.0)
        phase = torch.exp(-1j * two_pi_t * (KX*dx + KY*dy) / N)
        F2 = F2 * phase

        finite = torch.isfinite(F2).view(M, -1).all(dim=1)
        if not bool(finite.all()):
            bad_js = (~finite).nonzero(as_tuple=False).view(-1).tolist()
            for j in bad_js:
                _log_nan_particle(b0 + j)
            F2 = torch.nan_to_num(F2, nan=0.0, posinf=0.0, neginf=0.0)



        # Accumulate one image at a time so checkpoints hit exact targets
        for j in range(M):
            if SKIP_BAD and (not keep[j]):
                processed += 1
                if counts:
                    _maybe_checkpoint_gpu(processed, acc_t, weight_t)
                continue            
            Rmat = Rmats_t[b0 + j]
            F2_j = F2[j]
            xp = Rmat[0,0]*FX + Rmat[0,1]*FY
            yp = Rmat[1,0]*FX + Rmat[1,1]*FY
            zp = Rmat[2,0]*FX + Rmat[2,1]*FY

            conj  = xp < 0
            xp    = torch.where(conj, -xp, xp)
            yp    = torch.where(conj, -yp, yp)
            zp    = torch.where(conj, -zp, zp)
            F2_j  = torch.where(conj, torch.conj(F2_j), F2_j)

            kx_r = torch.round(xp / freq_unit).to(torch.int64)
            ky_r = torch.round(yp / freq_unit).to(torch.int64)
            kz_r = torch.round(zp / freq_unit).to(torch.int64)

            valid = (kx_r.abs() <= halfN) & (ky_r.abs() <= halfN) & (kz_r.abs() <= halfN)
            if valid.any():
                ix   = torch.where(kx_r[valid]>=0, kx_r[valid], kx_r[valid]+N)
                iy   = torch.where(ky_r[valid]>=0, ky_r[valid], ky_r[valid]+N)
                iz   = torch.where(kz_r[valid]>=0, kz_r[valid], kz_r[valid]+N)
                vals = F2_j[valid]
                acc_t.index_put_((iz, iy, ix), vals,   accumulate=True)
                weight_t.index_put_((iz, iy, ix), torch.ones_like(vals.real), accumulate=True)

            processed += 1
            if counts:
                _maybe_checkpoint_gpu(processed, acc_t, weight_t)

        elapsed = time.time() - start_time
        avg     = elapsed / max(processed,1)
        rem     = avg * (total - processed)
        pct     = 100.0 * processed / max(1,total)
        sys.stdout.write(
            f"\rProcessed particles {processed}/{total} ({pct:.1f}%) — "
            f"elapsed {time.strftime('%H:%M:%S', time.gmtime(elapsed))} — ETA {time.strftime('%H:%M:%S', time.gmtime(rem))}"
        )
        sys.stdout.flush()

    if SKIP_BAD and skipped_rows and checkpoint_out_path:
        log_path = _error_log_path_from_checkpoint(checkpoint_out_path)
        append_rec_error_log(log_path, skipped_rows)

    total_time = time.time() - start_time
    print(f"\nGPU reconstruction done in {time.strftime('%H:%M:%S', time.gmtime(total_time))}")
    return acc_t.cpu().numpy(), weight_t.cpu().numpy()

# ======================= Post-processing =======================
# ======================= Low-pass helpers =======================

def _compute_R_out(N: int, outer_scale: Optional[float], add_px: Optional[int]) -> float:
    """
    Compute outer radius R_out (in Fourier index units).
    Base Nyquist per-axis radius is N/2. We allow extending up to the cube corner at sqrt(3)*N/2.
    R_out = clamp( (outer_scale * N/2) + add_px, 1, sqrt(3)*N/2 )
    """
    halfN = N // 2
    rmax = math.sqrt(3.0) * halfN
    s = 1.0 if (outer_scale is None) else float(outer_scale)
    apx = 0.0 if (add_px is None) else float(add_px)
    Rout = halfN * s + apx
    Rout = max(1.0, min(rmax, Rout))
    return float(Rout)

def _spherical_soft_mask_numpy(
    N: int,
    start_frac: Optional[float] = 0.9,
    rolloff_px: Optional[int] = None,
    outer_scale: Optional[float] = 1.1,
    add_px: Optional[int] = None,
) -> np.ndarray:
    """
    Spherical raised-cosine mask (in index units).
    Passband: r <= R_in          (mask = 1)
    Roll-off: R_in < r <= R_out  (cosine from 1 -> 0)
    Stopband: r > R_out          (mask = 0)
    """
    halfN = N // 2
    R_out = _compute_R_out(N, outer_scale, add_px)

    if rolloff_px is not None:
        w = max(0.0, float(rolloff_px))
        R_in = max(0.0, R_out - w)
    else:
        sf = 0.9 if start_frac is None else float(start_frac)
        sf = max(0.0, min(1.1, sf))
        R_in = sf * R_out
        w = R_out - R_in

    f = np.arange(N, dtype=np.float32)
    f[f > halfN] -= N
    kz = f.reshape(N, 1, 1)
    ky = f.reshape(1, N, 1)
    kx = f.reshape(1, 1, N)
    r = np.sqrt(kx * kx + ky * ky + kz * kz)

    mask = np.ones_like(r, dtype=np.float32)
    mask[r > R_out] = 0.0
    if w <= 0.0:
        return mask
    roll = (r > R_in) & (r <= R_out)
    mask[roll] = 0.5 * (1.0 + np.cos(np.pi * (r[roll] - R_in) / w))
    return mask

def _spherical_soft_mask_torch(
    N: int,
    start_frac: Optional[float],
    rolloff_px: Optional[int],
    outer_scale: Optional[float],
    add_px: Optional[int],
    device,
) -> "torch.Tensor":
    halfN = N // 2
    R_out = _compute_R_out(N, outer_scale, add_px)

    if rolloff_px is not None:
        w = max(0.0, float(rolloff_px))
        R_in = max(0.0, R_out - w)
    else:
        sf = 0.9 if start_frac is None else float(start_frac)
        sf = max(0.0, min(1.1, sf))
        R_in = sf * R_out
        w = R_out - R_in

    f = torch.arange(N, device=device, dtype=torch.float32)
    f[f > halfN] -= N
    kz = f.view(N, 1, 1)
    ky = f.view(1, N, 1)
    kx = f.view(1, 1, N)
    r = torch.sqrt(kx * kx + ky * ky + kz * kz)

    mask = torch.ones_like(r, dtype=torch.float32)
    mask[r > R_out] = 0.0
    if w <= 0.0:
        return mask
    roll = (r > R_in) & (r <= R_out)
    mask[roll] = 0.5 * (1.0 + torch.cos(math.pi * (r[roll] - R_in) / w))
    return mask


def post_process(
    acc,
    weight,
    pixel_size: float,
    device: str = "cpu",
    lp_start: float = 0.9,
    lp_rolloff: Optional[int] = None,
    lp_outer_scale: float = 1.1,
    lp_add_px: Optional[int] = None,
) -> np.ndarray:
    N = acc.shape[0]; halfN = N // 2

    # Normalise Fourier voxels by hit count
    if device == "cpu":
        int_freq = np.arange(N, dtype=np.int32)
        int_freq[int_freq > halfN] -= N
        nzmask = weight > 0
        acc[nzmask] /= weight[nzmask]
        acc[~nzmask] = 0

        # Spherical low-pass
        lpf = _spherical_soft_mask_numpy(
            N,
            start_frac=lp_start,
            rolloff_px=lp_rolloff,
            outer_scale=lp_outer_scale,
            add_px=lp_add_px,
        )
        acc = acc * lpf.astype(np.float32)

        vol = np.fft.ifftn(acc).real
        vol = np.roll(vol, shift=(-halfN, -halfN, -halfN), axis=(0, 1, 2))
        return vol.astype(np.float32)

    if not TORCH_AVAILABLE:
        raise RuntimeError("GPU post_process requested but PyTorch is not available.")
    import torch
    if not torch.cuda.is_available():
        raise RuntimeError("GPU post_process requested but no CUDA device is visible.")

    device_t = torch.device("cuda")
    acc_t = torch.tensor(acc, device=device_t)
    weight_t = torch.tensor(weight, device=device_t)

    int_freq = torch.arange(N, device=device_t)
    int_freq[int_freq > halfN] -= N
    nzmask = weight_t > 0
    acc_t[nzmask] = acc_t[nzmask] / weight_t[nzmask]
    acc_t[~nzmask] = 0

    lpf_t = _spherical_soft_mask_torch(
        N,
        start_frac=lp_start,
        rolloff_px=lp_rolloff,
        outer_scale=lp_outer_scale,
        add_px=lp_add_px,
        device=device_t,
    )
    acc_t = acc_t * lpf_t

    vol_t = torch.fft.ifftn(acc_t).real
    vol_t = torch.roll(vol_t, shifts=(-halfN, -halfN, -halfN), dims=(0, 1, 2))
    return vol_t.cpu().numpy().astype(np.float32)
