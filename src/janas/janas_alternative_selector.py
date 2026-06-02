# File: selector_utils.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
#
# Direct scoring kernels used by assessParticles’ block scorers when tags
# request FIR-based SCI (non-recursive) or LoG. Implemented in NumPy/FFT,
# mask-aware, and amplitude-equalised in Fourier domain to mirror the SCI behaviour.
#
# Terminology:
# - FIR Gaussian derivatives: non-recursive (finite impulse response) realisations
#   of Gaussian derivative operators (exact up to numerical precision).
# - Recursive/Deriche: IIR approximations (not implemented here).
#
# Public APIs:
#   - sci_fir_gaussian_derivative_2d(I, R, M, sigma_pixels=1.0) → float
#   - log_score_2d(I, R, M, sigma_pixels=1.0) → float
# Backward-compatibility alias:
#   - sci_non_deriche_2d = sci_fir_gaussian_derivative_2d

from __future__ import annotations
import numpy as np
from numpy.fft import fftn, ifftn, fftfreq
import os 

import math, numpy as np
from scipy.ndimage import distance_transform_edt  # NEW

from typing import Optional, Tuple, Dict, List
from janas import IO_utils, locres_utils

import janas.janas_core as janas_core


_DEBUG = os.environ.get("JANAS_DEBUG", "0") not in ("0", "", "false", "False")
if _DEBUG and not _PRINTED_FLAGS["sci_fir_cpu"]:
    print("[DEBUG] SCIFIR CPU sigma=...", flush=True)

# --------------------------- amplitude equalisation ---------------------------

def _amplitude_equalise_fft(I: np.ndarray, R: np.ndarray, eps: float = 1e-8):
    """
    Equalise Fourier amplitudes of I and R to their average magnitude spectrum.
    This matches the amplitude-equalisation performed in the canonical SCI.
    """
    if I.shape != R.shape:
        raise ValueError("Images must have the same shape for amplitude equalisation.")
    FI = fftn(I)
    FR = fftn(R)
    magI = np.abs(FI) + eps
    magR = np.abs(FR) + eps
    ampAvg = 0.5 * (magI + magR)
    Ieq = ifftn(FI * (ampAvg / magI)).real
    Req = ifftn(FR * (ampAvg / magR)).real
    return Ieq.astype(np.float32), Req.astype(np.float32)

# --------------------------- masked NCC (mean-subtracted) ---------------------------

def _masked_ncc(I: np.ndarray, R: np.ndarray, M: np.ndarray, min_mask_value: float = 0.0, eps: float = 1e-12) -> float:
    """
    Mean-subtracted, normalised cross-correlation under a (soft) mask with optional
    thresholding. Negative values are truncated to zero to match SCI behaviour.
    """
    if (I.shape != R.shape) or (I.shape != M.shape):
        raise ValueError("I, R, and M must share the same shape.")
    w = np.asarray(M, dtype=np.float64)
    if min_mask_value > 0.0:
        w = np.where(w > min_mask_value, w, 0.0)
    sw = w.sum()
    if sw <= 0.0:
        return 0.0

    I = np.asarray(I, dtype=np.float64)
    R = np.asarray(R, dtype=np.float64)

    muI = (w * I).sum() / sw
    muR = (w * R).sum() / sw
    I0 = I - muI
    R0 = R - muR

    num = (w * I0 * R0).sum()
    den = np.sqrt((w * I0 * I0).sum() * (w * R0 * R0).sum()) + eps
    cc = num / den
    if not np.isfinite(cc):
        return 0.0
    return float(max(0.0, cc))

# --------------------------- FIR Gaussian derivative filters (Fourier) ---------------------------

def _gaussian_derivative_filters_2d(shape, sigma_pixels: float):
    """
    Construct Fourier-domain FIR (non-recursive) Gaussian and Gaussian-derivative
    transfer functions for a 2D image of given shape.
    """
    ny, nx = shape
    fy = fftfreq(ny)
    fx = fftfreq(nx)
    FY, FX = np.meshgrid(fy, fx, indexing="xy")
    FX = FX.T
    FY = FY.T

    t = float(sigma_pixels) ** 2
    two_pi = 2.0 * np.pi
    sq = FX**2 + FY**2
    G = np.exp(-2.0 * (np.pi**2) * t * sq)

    Hx  = (1j * two_pi * FX) * G
    Hy  = (1j * two_pi * FY) * G
    Hxx = -(two_pi * FX)**2 * G
    Hyy = -(two_pi * FY)**2 * G
    return Hx, Hy, Hxx, Hyy, G

def _apply_fir_gaussian_derivatives(I: np.ndarray, sigma_pixels: float):
    """
    Return FIR Gaussian fields for 2D image I via FFT:
      L, Lx, Ly, Lxx, Lyy
    """
    Hx, Hy, Hxx, Hyy, G = _gaussian_derivative_filters_2d(I.shape, sigma_pixels)
    F = fftn(I)
    out = {
        "L":    ifftn(G   * F).real.astype(np.float32),
        "L_x":  ifftn(Hx  * F).real.astype(np.float32),
        "L_y":  ifftn(Hy  * F).real.astype(np.float32),
        "L_xx": ifftn(Hxx * F).real.astype(np.float32),
        "L_yy": ifftn(Hyy * F).real.astype(np.float32),
    }
    return out

# --------------------------- LoG filter (Fourier) ---------------------------

def _log_filter_2d(shape, sigma_pixels: float):
    """
    Fourier-domain Laplacian-of-Gaussian (LoG) kernel.
    Overall scale is irrelevant because NCC is scale-invariant.
    """
    ny, nx = shape
    fy = fftfreq(ny)
    fx = fftfreq(nx)
    FY, FX = np.meshgrid(fy, fx, indexing="xy")
    FX = FX.T
    FY = FY.T

    t = float(sigma_pixels) ** 2
    two_pi = 2.0 * np.pi
    sq = FX**2 + FY**2
    G = np.exp(-2.0 * (np.pi**2) * t * sq)
    LoG = - (two_pi**2) * (2.0 * t) * sq * G
    return LoG

def _apply_log(I: np.ndarray, sigma_pixels: float):
    F = fftn(I)
    H = _log_filter_2d(I.shape, sigma_pixels)
    return ifftn(H * F).real.astype(np.float32)

# --------------------------- public scoring APIs ---------------------------

def sci_fir_gaussian_derivative_2d(I: np.ndarray, R: np.ndarray, M: np.ndarray, sigma_pixels: float = 1.0) -> float:
    """
    SCI using FIR (non-recursive) Gaussian-derivative filters:
      SCI = CC(Ieq,Req) * CC(Ix,Rx) * CC(Iy,Ry) * CC(Ixx,Rxx) * CC(Iyy,Ryy)
    where Ieq, Req are amplitude-equalised images in Fourier domain.
    """
    Ieq, Req = _amplitude_equalise_fft(I, R)
    Id = _apply_fir_gaussian_derivatives(Ieq, sigma_pixels)
    Rd = _apply_fir_gaussian_derivatives(Req, sigma_pixels)

    c00 = _masked_ncc(Ieq,       Req,       M)
    cx  = _masked_ncc(Id["L_x"], Rd["L_x"], M)
    cy  = _masked_ncc(Id["L_y"], Rd["L_y"], M)
    cxx = _masked_ncc(Id["L_xx"],Rd["L_xx"],M)
    cyy = _masked_ncc(Id["L_yy"],Rd["L_yy"],M)
    sci = c00 * cx * cy * cxx * cyy
    return float(sci if np.isfinite(sci) else 0.0)


def CC_2d(I: np.ndarray, R: np.ndarray, M: np.ndarray, sigma_pixels: float = 1.0) -> float:
    c00 = _masked_ncc(I,       R,       M)
    return float(c00 if np.isfinite(c00) else 0.0)

def sci_2d(I, R, M, sigmaBlur="1", preprocessingMethod="unprocessed"):
    """
    Structural Cross-correlation Index (SCI), 2D version.

    Matches assessParticles.py comparisonMethod == "SCI":
      - Fourier amplitude equalisation between I and R
      - then janas_core.MaskedImageComparison on the equalised images with mask M

    Parameters
    ----------
    I, R : 2D numpy arrays (ny, nx), float32/float64
        Input image and reference image.
    M : 2D numpy array (ny, nx), float32
        Mask in image space (0/1 or soft weights).
    sigmaBlur : str
        Passed through to MaskedImageComparison to keep signature compatible.
    preprocessingMethod : str
        Passed through to MaskedImageComparison.

    Returns
    -------
    float
        SCI score.
    """
    # Ensure float arrays
    I2d = np.asarray(I, dtype=np.float32, order="C")
    R2d = np.asarray(R, dtype=np.float32, order="C")
    M2d = np.asarray(M, dtype=np.float32, order="C")

    if I2d.shape != R2d.shape or I2d.shape != M2d.shape:
        raise ValueError(f"Shape mismatch: I={I2d.shape}, R={R2d.shape}, M={M2d.shape}")

    ny, nx = I2d.shape

    # --- Fourier amplitude equalisation (as in assessParticles.py) ---
    I_fft = np.fft.fft2(I2d)
    R_fft = np.fft.fft2(R2d)

    I_abs = np.abs(I_fft) + 1e-7
    R_abs = np.abs(R_fft) + 1e-7
    ampAvg = 0.5 * (I_abs + R_abs)

    I_out = np.fft.ifft2(I_fft * (ampAvg / I_abs)).real.astype(np.float32, copy=False)
    R_out = np.fft.ifft2(R_fft * (ampAvg / R_abs)).real.astype(np.float32, copy=False)

    # Flatten to lists because janas_core expects flat vectors (as used in assessParticles.py)
    I_out_list = I_out.ravel(order="C").tolist()
    R_out_list = R_out.ravel(order="C").tolist()
    M_list     = M2d.ravel(order="C").tolist()

    # Call the same comparator used in the pipeline
    # (comparisonMethod="SCI" is what your pipeline uses in this branch)
    score = janas_core.MaskedImageComparison(
        R_out_list,
        I_out_list,
        M_list,
        nx, ny, 1,
        "SCI",
        preprocessingMethod,
        str(sigmaBlur),
    )
    return float(score)




def log_score_2d_equalized(I: np.ndarray, R: np.ndarray, M: np.ndarray, sigma_pixels: float = 1.0) -> float:
    """
    LoG variant with the same normalisation logic:
      score = CC(Ieq,Req) * CC(LoG(Ieq), LoG(Req))
    """
    Ieq, Req = _amplitude_equalise_fft(I, R)
    LI = _apply_log(Ieq, sigma_pixels)
    LR = _apply_log(Req, sigma_pixels)

    c00  = _masked_ncc(Ieq, Req, M)
    clog = _masked_ncc(LI,  LR,  M)
    score = c00 * clog
    return float(score if np.isfinite(score) else 0.0)

def log_score_2d(I: np.ndarray, R: np.ndarray, M: np.ndarray, sigma_pixels: float = 1.0) -> float:
    LI = _apply_log(I, sigma_pixels)
    LR = _apply_log(R, sigma_pixels)
    score  = _masked_ncc(LI, LR, M)
    return float(score if np.isfinite(score) else 0.0)

# --- GPU backend (optional, via CuPy) ----------------------------------------
def _get_cupy():
    try:
        import cupy as cp
        return cp
    except Exception:
        return None

def _amplitude_equalise_fft_gpu(cp, I, R, eps: float = 1e-8, device: int = 0):
    with cp.cuda.Device(device):
        I_d = cp.asarray(I, dtype=cp.float32)
        R_d = cp.asarray(R, dtype=cp.float32)
        FI = cp.fft.fftn(I_d)
        FR = cp.fft.fftn(R_d)
        magI = cp.abs(FI) + eps
        magR = cp.abs(FR) + eps
        ampAvg = 0.5 * (magI + magR)
        Ieq = cp.fft.ifftn(FI * (ampAvg / magI)).real
        Req = cp.fft.ifftn(FR * (ampAvg / magR)).real
        return Ieq.astype(cp.float32), Req.astype(cp.float32)

def _masked_ncc_gpu(cp, I, R, M, min_mask_value: float = 0.0, eps: float = 1e-12, device: int = 0):
    with cp.cuda.Device(device):
        I = cp.asarray(I, dtype=cp.float32)
        R = cp.asarray(R, dtype=cp.float32)
        w = cp.asarray(M, dtype=cp.float32)
        if min_mask_value > 0.0:
            w = cp.where(w > min_mask_value, w, cp.float32(0.0))
        sw = float(w.sum().get())
        if sw <= 0.0:
            return 0.0
        muI = float((w * I).sum().get()) / sw
        muR = float((w * R).sum().get()) / sw
        I0 = I - muI
        R0 = R - muR
        num = float((w * I0 * R0).sum().get())
        den = float(cp.sqrt((w * I0 * I0).sum() * (w * R0 * R0).sum() + eps).get())
        cc = 0.0 if den == 0.0 else (num / den)
        return float(max(0.0, cc))

def _gaussian_derivative_filters_2d_gpu(cp, shape, sigma_pixels: float, device: int = 0):
    with cp.cuda.Device(device):
        ny, nx = int(shape[0]), int(shape[1])
        fy = cp.fft.fftfreq(ny)
        fx = cp.fft.fftfreq(nx)
        FY, FX = cp.meshgrid(fy, fx, indexing="xy")
        FX = FX.T
        FY = FY.T
        t = float(sigma_pixels) ** 2
        two_pi = 2.0 * cp.pi
        sq = FX**2 + FY**2
        G = cp.exp(-2.0 * (cp.pi**2) * t * sq)
        Hx  = (1j * two_pi * FX) * G
        Hy  = (1j * two_pi * FY) * G
        Hxx = -(two_pi * FX)**2 * G
        Hyy = -(two_pi * FY)**2 * G
        return Hx, Hy, Hxx, Hyy, G

def _apply_fir_derivatives_gpu(cp, I, sigma_pixels: float, device: int = 0):
    with cp.cuda.Device(device):
        F = cp.fft.fftn(cp.asarray(I, dtype=cp.float32))
        Hx, Hy, Hxx, Hyy, G = _gaussian_derivative_filters_2d_gpu(cp, I.shape, sigma_pixels, device=device)
        outL    = cp.fft.ifftn(G   * F).real
        outLx   = cp.fft.ifftn(Hx  * F).real
        outLy   = cp.fft.ifftn(Hy  * F).real
        outLxx  = cp.fft.ifftn(Hxx * F).real
        outLyy  = cp.fft.ifftn(Hyy * F).real
        return {
            "L":    outL.astype(cp.float32),
            "L_x":  outLx.astype(cp.float32),
            "L_y":  outLy.astype(cp.float32),
            "L_xx": outLxx.astype(cp.float32),
            "L_yy": outLyy.astype(cp.float32),
        }

def _log_filter_2d_gpu(cp, shape, sigma_pixels: float, device: int = 0):
    with cp.cuda.Device(device):
        ny, nx = int(shape[0]), int(shape[1])
        fy = cp.fft.fftfreq(ny)
        fx = cp.fft.fftfreq(nx)
        FY, FX = cp.meshgrid(fy, fx, indexing="xy")
        FX = FX.T
        FY = FY.T
        t = float(sigma_pixels) ** 2
        two_pi = 2.0 * cp.pi
        sq = FX**2 + FY**2
        G = cp.exp(-2.0 * (cp.pi**2) * t * sq)
        LoG = - (two_pi**2) * (2.0 * t) * sq * G
        return LoG

def _apply_log_gpu(cp, I, sigma_pixels: float, device: int = 0):
    with cp.cuda.Device(device):
        F = cp.fft.fftn(cp.asarray(I, dtype=cp.float32))
        H = _log_filter_2d_gpu(cp, I.shape, sigma_pixels, device=device)
        return cp.fft.ifftn(H * F).real.astype(cp.float32)

def sci_fir_gaussian_derivative_2d_gpu(I, R, M, sigma_pixels: float = 1.0, device: int = 0):
    """
    Versione GPU (CuPy) del kernel.
    Ritorna float (trasferito su host).
    """
    cp = _get_cupy()
    if cp is None:
        return sci_fir_gaussian_derivative_2d(I, R, M, sigma_pixels=sigma_pixels)
    Ieq, Req = _amplitude_equalise_fft_gpu(cp, I, R, device=device)
    Id = _apply_fir_derivatives_gpu(cp, Ieq, sigma_pixels, device=device)
    Rd = _apply_fir_derivatives_gpu(cp, Req, sigma_pixels, device=device)
    c00 = _masked_ncc_gpu(cp, Ieq,    Req,    M, device=device)
    cx  = _masked_ncc_gpu(cp, Id["L_x"],  Rd["L_x"],  M, device=device)
    cy  = _masked_ncc_gpu(cp, Id["L_y"],  Rd["L_y"],  M, device=device)
    cxx = _masked_ncc_gpu(cp, Id["L_xx"], Rd["L_xx"], M, device=device)
    cyy = _masked_ncc_gpu(cp, Id["L_yy"], Rd["L_yy"], M, device=device)
    return float(c00 * cx * cy * cxx * cyy)

def log_score_2d_gpu(I, R, M, sigma_pixels: float = 1.0, device: int = 0):
    """
    Versione GPU (CuPy) del LoG score: CC(Ieq,Req) * CC(LoG(Ieq),LoG(Req)).
    """
    cp = _get_cupy()
    if cp is None:
        return log_score_2d(I, R, M, sigma_pixels=sigma_pixels)
    Ieq, Req = _amplitude_equalise_fft_gpu


# --------------------------- Automatic sigma ------------------------------------------
# --- Auto-sigma from half-maps (mask-aware FSC, uses IO_utils + locres_utils) ---



def auto_sigma_from_halfmapsOld(
    map_h1_paths: List[str],
    map_h2_paths: List[str],
    mask_path: Optional[str],
    apix: float,
    *,
    gamma_at_f0143: float = 0.5,
    sigma_scale: float = 2.0,
    sigma_minmax_px: Tuple[float, float] = (0.25, 8.0),
    mask_threshold: float = 0.5,
    mask_soft_edge_add: float = 0.0,   # NEW (pixels)
) -> Tuple[float, Dict[str, float]]:
    assert len(map_h1_paths) == len(map_h2_paths) and len(map_h1_paths) > 0, "Need equal H1/H2 lists."

    # Load + prepare mask (optional)
    mask_vol = None
    if mask_path:
        mask_vol, _ = IO_utils.read_mrc_data(mask_path)
        # binarise
        bin_mask = (mask_vol >= float(mask_threshold)).astype(np.uint8)

        if mask_soft_edge_add and mask_soft_edge_add > 0:
            # Soft rim outside the binary mask:
            # value decays linearly from 1 at boundary to 0 after 'mask_soft_edge_add' pixels.
            # 1) distance outside (to nearest '1' voxel)
            dist_out = distance_transform_edt(1 - bin_mask)
            rim = np.clip(1.0 - (dist_out / float(mask_soft_edge_add)), 0.0, 1.0)
            # 2) keep interior as 1, take max with rim for smooth extension
            mask_vol = np.maximum(bin_mask.astype(np.float32), rim.astype(np.float32))
        else:
            mask_vol = bin_mask.astype(np.float32)

    f05_list: List[float] = []
    f0143_list: List[float] = []
    sigma_px_list: List[float] = []
    sigma_A_list: List[float] = []

    nyq = 0.5 / float(apix)
    sig_min, sig_max = float(sigma_minmax_px[0]), float(sigma_minmax_px[1])
    two_pi_sq = (2.0 * np.pi) ** 2

    for p1, p2 in zip(map_h1_paths, map_h2_paths):
        V1, _ = IO_utils.read_mrc_data(p1)
        V2, _ = IO_utils.read_mrc_data(p2)
        if mask_vol is not None:
            if mask_vol.shape != V1.shape:
                raise ValueError(f"Mask shape {mask_vol.shape} != map shape {V1.shape}")
            V1 = (V1 * mask_vol).astype(np.float32, copy=False)
            V2 = (V2 * mask_vol).astype(np.float32, copy=False)

        freqs, fsc = locres_utils.compute_fsc(V1, V2, float(apix))  # cycles/Å

        _, crossings, _, _, _, _ = locres_utils.find_FSC_resolutions_and_stats(
            freqs, fsc, thresholds=[0.5, 0.143], interp_factor=10, kind="linear"
        )

        # fallbacks if interpolation fails
        f05 = crossings.get(0.5)
        if f05 is None or not np.isfinite(f05):
            f05 = locres_utils.first_FSC_crossing(freqs, fsc, 0.5)

        f0143 = crossings.get(0.143)
        if f0143 is None or not np.isfinite(f0143):
            f0143 = locres_utils.first_FSC_crossing(freqs, fsc, 0.143)

        f_tgt = min(0.95 * nyq, f0143 if (f0143 and np.isfinite(f0143) and f0143 > 0.0) else 0.95 * nyq)
        k_tgt_px = float(f_tgt) * float(apix)

        base_sigma_px = math.sqrt(
            max(1e-12, math.log(1.0 / float(gamma_at_f0143))) / (two_pi_sq * (k_tgt_px ** 2))
        )
        sigma_px = float(np.clip(base_sigma_px, sig_min, sig_max)) * float(sigma_scale)
        sigma_A  = sigma_px * float(apix)

        f05_list.append(float(f05) if (f05 is not None) else float("nan"))
        f0143_list.append(float(f0143) if (f0143 is not None) else float("nan"))
        sigma_px_list.append(sigma_px)
        sigma_A_list.append(sigma_A)

    sigma_px_avg = float(np.nanmean(sigma_px_list))
    meta = {
        "sigma_px_avg": sigma_px_avg,
        "sigma_A_avg":  sigma_px_avg * float(apix),
        "sigma_px_list": sigma_px_list,
        "sigma_A_list":  sigma_A_list,
        "f05_list":      f05_list,
        "f0143_list":    f0143_list,
        "nyquist":       nyq,
        "gamma":         float(gamma_at_f0143),
        "scale":         float(sigma_scale),
        "clamp_px":      (sig_min, sig_max),
        "mask_soft_edge_add": float(mask_soft_edge_add),
    }
    return sigma_px_avg, meta




import math
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy.ndimage import distance_transform_edt


def auto_sigma_from_halfmaps(
    map_h1_paths: List[str],
    map_h2_paths: List[str],
    mask_path: Optional[str],
    apix: float,
    *,
    gamma_at_f0143: float = 0.5,
    sigma_scale: float = 2.0,
    mask_threshold: float = 0.5,
    mask_soft_edge_add: float = 0.0,  # pixels

    # Backward-compatibility: accept but do NOT enforce bounds
    sigma_minmax_px: Tuple[float, float] = (0.25, 8.0),

    # Variability control (NO hard bounds)
    sigma0_px: float = 1.0,           # prior/default scale in pixels
    shrink_k0_px: float = 0.20,       # cycles/pixel where shrinkage is ~50%
    shrink_power: float = 1.0,        # >1 makes transition sharper, <1 smoother

    # Printing
    verbose: bool = True,
    print_fn: Callable[[str], None] = print,

    # Extra backward-compatible sink for any other legacy kwargs
    **_unused_kwargs,
) -> Tuple[float, Dict[str, float]]:
    """
    FSC-linked sigma estimate with smooth variability reduction (no clipping, no quantisation).

    Notes:
    - sigma_minmax_px is accepted for backwards compatibility but is not used to clip sigma.
    - Any additional unexpected legacy kwargs are accepted via **_unused_kwargs.
    """
    assert len(map_h1_paths) == len(map_h2_paths) and len(map_h1_paths) > 0, \
        "Need equal H1/H2 lists."

    def _p(msg: str) -> None:
        if verbose:
            print_fn(msg)

    # Inform if legacy args are being ignored
    if verbose:
        _p(f"[auto_sigma] NOTE: sigma_minmax_px={sigma_minmax_px} provided but NOT enforced (no hard bounds).")
        if _unused_kwargs:
            _p(f"[auto_sigma] NOTE: Ignoring unused kwargs: {sorted(_unused_kwargs.keys())}")

    # Load + prepare mask (optional)
    mask_vol = None
    if mask_path:
        mask_vol, _ = IO_utils.read_mrc_data(mask_path)
        bin_mask = (mask_vol >= float(mask_threshold)).astype(np.uint8)

        if mask_soft_edge_add and mask_soft_edge_add > 0:
            dist_out = distance_transform_edt(1 - bin_mask)
            rim = np.clip(1.0 - (dist_out / float(mask_soft_edge_add)), 0.0, 1.0)
            mask_vol = np.maximum(bin_mask.astype(np.float32), rim.astype(np.float32))
            _p(f"[auto_sigma] Mask loaded: {mask_path} | threshold={mask_threshold} | soft_edge_add(px)={mask_soft_edge_add}")
        else:
            mask_vol = bin_mask.astype(np.float32)
            _p(f"[auto_sigma] Mask loaded: {mask_path} | threshold={mask_threshold} | soft_edge_add(px)=0")

    # Lists for reporting
    f05_list: List[float] = []
    f0143_list: List[float] = []
    f_tgt_list: List[float] = []
    k_tgt_px_list: List[float] = []

    sigma_raw_px_list: List[float] = []
    w_list: List[float] = []
    sigma_px_list: List[float] = []
    sigma_A_list: List[float] = []

    nyq = 0.5 / float(apix)  # cycles/Å
    two_pi_sq = (2.0 * np.pi) ** 2
    eps = 1e-12  # numerical floor only

    _p("[auto_sigma] ---- Parameters ----")
    _p(f"[auto_sigma] apix={apix:.8f} Å/px | Nyquist={nyq:.8f} cycles/Å")
    _p(f"[auto_sigma] gamma_at_f0143={gamma_at_f0143} | sigma_scale={sigma_scale}")
    _p(f"[auto_sigma] shrink: sigma0_px={sigma0_px} | shrink_k0_px={shrink_k0_px} | shrink_power={shrink_power}")
    _p("[auto_sigma] -------------------\n")

    for i, (p1, p2) in enumerate(zip(map_h1_paths, map_h2_paths), start=1):
        V1, _ = IO_utils.read_mrc_data(p1)
        V2, _ = IO_utils.read_mrc_data(p2)

        if mask_vol is not None:
            if mask_vol.shape != V1.shape:
                raise ValueError(f"Mask shape {mask_vol.shape} != map shape {V1.shape}")
            V1 = (V1 * mask_vol).astype(np.float32, copy=False)
            V2 = (V2 * mask_vol).astype(np.float32, copy=False)

        freqs, fsc = locres_utils.compute_fsc(V1, V2, float(apix))  # cycles/Å

        _, crossings, _, _, _, _ = locres_utils.find_FSC_resolutions_and_stats(
            freqs, fsc, thresholds=[0.5, 0.143], interp_factor=10, kind="linear"
        )

        f05 = crossings.get(0.5)
        if f05 is None or not np.isfinite(f05):
            f05 = locres_utils.first_FSC_crossing(freqs, fsc, 0.5)

        f0143 = crossings.get(0.143)
        if f0143 is None or not np.isfinite(f0143):
            f0143 = locres_utils.first_FSC_crossing(freqs, fsc, 0.143)

        if f0143 is not None and np.isfinite(f0143) and f0143 > 0.0:
            f_tgt = min(0.95 * nyq, float(f0143))
            f_tgt_src = "FSC0.143"
        else:
            f_tgt = 0.95 * nyq
            f_tgt_src = "NyqFallback"

        k_tgt_px = float(f_tgt) * float(apix)
        if not np.isfinite(k_tgt_px) or k_tgt_px <= 0.0:
            raise ValueError(f"Invalid k_tgt_px={k_tgt_px} (from f_tgt={f_tgt}, apix={apix})")

        # Raw sigma (definitionally FSC-linked)
        sigma_raw_px = float(sigma_scale) * math.sqrt(
            max(eps, math.log(1.0 / float(gamma_at_f0143))) / (two_pi_sq * (k_tgt_px ** 2))
        )

        # Smooth shrinkage weight: w in (0,1)
        w = (k_tgt_px / (k_tgt_px + float(shrink_k0_px))) ** float(shrink_power)

        # Geometric shrinkage towards sigma0 (no bounds)
        sigma_px = float(sigma0_px) * (max(eps, sigma_raw_px) / float(sigma0_px)) ** float(w)
        sigma_A = sigma_px * float(apix)

        # Store
        f05_list.append(float(f05) if (f05 is not None) else float("nan"))
        f0143_list.append(float(f0143) if (f0143 is not None) else float("nan"))
        f_tgt_list.append(float(f_tgt))
        k_tgt_px_list.append(float(k_tgt_px))
        sigma_raw_px_list.append(float(sigma_raw_px))
        w_list.append(float(w))
        sigma_px_list.append(float(sigma_px))
        sigma_A_list.append(float(sigma_A))

        # Print per pair
        d0143 = (1.0 / float(f0143)) if (f0143 is not None and np.isfinite(f0143) and f0143 > 0) else float("nan")
        d05 = (1.0 / float(f05)) if (f05 is not None and np.isfinite(f05) and f05 > 0) else float("nan")

        _p(f"[auto_sigma] Pair {i}/{len(map_h1_paths)}")
        _p(f"  H1: {p1}")
        _p(f"  H2: {p2}")
        _p(f"  FSC: f0.5={float(f05) if f05 is not None else float('nan'):.8f} cycles/Å (d0.5={d05:.8f} Å) | "
           f"f0.143={float(f0143) if f0143 is not None else float('nan'):.8f} cycles/Å (d0.143={d0143:.8f} Å)")
        _p(f"  Target: f_tgt={f_tgt:.8f} cycles/Å ({f_tgt_src}) | k_tgt={k_tgt_px:.8f} cycles/px")
        _p(f"  Sigma raw: {sigma_raw_px:.8f} px")
        _p(f"  Shrink weight w(k): {w:.6f}")
        _p(f"  Sigma final: {sigma_px:.8f} px | {sigma_A:.8f} Å\n")

    sigma_px_avg = float(np.nanmean(sigma_px_list))
    sigma_A_avg = sigma_px_avg * float(apix)

    _p("[auto_sigma] ==== Final outcome ====")
    _p(f"[auto_sigma] sigma_px_avg={sigma_px_avg:.8f} px | sigma_A_avg={sigma_A_avg:.8f} Å")
    _p(f"[auto_sigma] per-pair sigma_px={', '.join(f'{v:.6f}' for v in sigma_px_list)}")
    _p("[auto_sigma] =======================\n")

    meta: Dict[str, float] = {
        "sigma_px_avg": sigma_px_avg,
        "sigma_A_avg": sigma_A_avg,
        "sigma_px_list": sigma_px_list,
        "sigma_A_list": sigma_A_list,
        "sigma_raw_px_list": sigma_raw_px_list,
        "w_list": w_list,
        "f05_list": f05_list,
        "f0143_list": f0143_list,
        "f_tgt_list": f_tgt_list,
        "k_tgt_px_list": k_tgt_px_list,
        "nyquist": nyq,
        "gamma": float(gamma_at_f0143),
        "scale": float(sigma_scale),
        "sigma0_px": float(sigma0_px),
        "shrink_k0_px": float(shrink_k0_px),
        "shrink_power": float(shrink_power),
        "mask_soft_edge_add": float(mask_soft_edge_add),
        "legacy_sigma_minmax_px_passed": True,
        "legacy_sigma_minmax_px_value": (float(sigma_minmax_px[0]), float(sigma_minmax_px[1])),
        "unused_kwargs": list(_unused_kwargs.keys()),
    }

    return sigma_px_avg, meta





# --------------------------- backward-compatibility aliases ---------------------------

# Keep previous name used in earlier integration steps.
sci_non_deriche_2d = sci_fir_gaussian_derivative_2d
