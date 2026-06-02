#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# File: janas_mapProcess.py
# Utilities for map processing (B-factor sharpening, FSC weighting, low-pass)

from typing import Tuple
import numpy as np
import math


# Local modules
from janas import janas_core
from janas import utils



def _safe_angpix(spacing) -> float:
    """Normalise utils.get_MRC_map_pixel_spacing() return to a float Å/px."""
    try:
        return float(spacing[0])
    except (TypeError, IndexError):
        return float(spacing)

def _load_mrc_as_array(fname: str) -> Tuple[np.ndarray, Tuple[int, int, int], float]:
    """
    Reads an MRC via janas_core.* and returns (volume[Z,Y,X], shape, angpix).
    """
    size = janas_core.sizeMRC(fname)          # (nz, ny, nx)
    vol = np.array(janas_core.ReadMRC(fname)).reshape(size)
    angpix = _safe_angpix(utils.get_MRC_map_pixel_spacing(fname))
    return vol.astype(np.float32, copy=False), size, float(angpix)


def _fourier_freq_grid(shape, angpix):
    """
    Return broadcastable frequency grids Fx, Fy, Fz in 1/Å for a volume of shape (nz, ny, nx).
    Axes match numpy's fftn order (z, y, x).
    """
    nz, ny, nx = shape
    fz = np.fft.fftfreq(nz, d=angpix)  # (nz,)
    fy = np.fft.fftfreq(ny, d=angpix)  # (ny,)
    fx = np.fft.fftfreq(nx, d=angpix)  # (nx,)
    Fz = fz[:, None, None]             # (nz,1,1)
    Fy = fy[None, :, None]             # (1,ny,1)
    Fx = fx[None, None, :]             # (1,1,nx)
    return Fx, Fy, Fz



def _radial_bin_mean(absFT, s_grid, angpix, fit_minres, fit_maxres):
    """
    Build Guinier arrays:
      x = 1/d^2 with d = 1/s; y = ln(<|F|>) per shell; w = 1 within fit range else 0.
    Excludes beyond Nyquist (d < 2*angpix).
    """
    s = s_grid.ravel()
    a = absFT.ravel()
    # Bin by s to shells (use uniform s bins)
    s_max = 0.5 / angpix
    n_bins = int(min(absFT.shape))  # simple dense binning
    edges = np.linspace(0, s_max, n_bins + 1)
    idx = np.clip(np.searchsorted(edges, s, side="right") - 1, 0, n_bins - 1)

    sums = np.bincount(idx, weights=a, minlength=n_bins)
    counts = np.bincount(idx, minlength=n_bins)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_abs = np.where(counts > 0, sums / counts, 0.0)

    s_centers = 0.5 * (edges[:-1] + edges[1:])
    # Convert to d (Å), guard s=0
    d = np.where(s_centers > 0, 1.0 / s_centers, np.inf)

    # Nyquist guard: keep only d >= 2*angpix
    valid = d >= (2.0 * angpix)
    d = d[valid]
    mean_abs = mean_abs[valid]

    # Fit window [fit_maxres, fit_minres] in Å; if fit_maxres==0 -> use all up to Nyquist
    if fit_maxres <= 0:
        in_fit = (d <= fit_minres)
    else:
        in_fit = (d <= fit_minres) & (d >= fit_maxres)

    with np.errstate(divide="ignore"):
        y = np.where(mean_abs > 0, np.log(mean_abs), np.nan)
    x = 1.0 / (d ** 2)
    w = (in_fit & np.isfinite(y))

    return x, y, w


def _apply_bfactor_inplace(FT, s_grid, B):
    """Multiply by exp(-B * s^2 / 4) with s in Å^-1."""
    if B == 0.0:
        return
    arg = -(B * (s_grid.astype(np.float64) ** 2)) / 4.0
    W = np.exp(arg)  # float64
    FT *= W.astype(FT.dtype, copy=False)


def normalize_amplitudes(
    hmap1_fname: str,
    hmap2_fname: str,
    fit_minres: float = 15.0,
    fit_maxres: float = 0.0,
):
    """
    Estimate a single B-factor from the average of two half-maps and apply
    that same B-factor to both half-maps.

    The B-factor is obtained via a Guinier fit of ln(<|F|>) vs 1/d^2 using
    `_radial_bin_mean`, analogous to the automatic mode of `bfactor_utils`.

    Parameters
    ----------
    hmap1_fname, hmap2_fname : str
        Filenames of the two half-maps (MRC). They must have identical
        shape and pixel size.
    fit_minres : float
        Upper resolution (Å) of the fit window (low-resolution end).
    fit_maxres : float
        Lower resolution (Å) of the fit window (high-resolution end).
        If <= 0, use all shells up to Nyquist.

    Returns
    -------
    vol1_out, vol2_out : np.ndarray
        Sharpened half-maps, shape (nz, ny, nx), dtype float32.
    B : float
        Estimated B-factor (Å^2) applied to both half-maps.
    """
    # --- Load half-maps as 3D arrays ---
    vol1, shape1, ang1 = _load_mrc_as_array(hmap1_fname)
    vol2, shape2, ang2 = _load_mrc_as_array(hmap2_fname)

    if shape1 != shape2:
        raise ValueError("normalize_amplitudes: half-maps must have identical shapes.")
    if not np.isclose(ang1, ang2):
        raise ValueError("normalize_amplitudes: half-maps must have identical pixel size.")

    nz, ny, nx = shape1
    angpix = float(ang1)

    # --- Combined map for B-factor fit ---
    vol_comb = 0.5 * (
        vol1.astype(np.float64, copy=False) +
        vol2.astype(np.float64, copy=False)
    )

    FT_comb = np.fft.fftn(vol_comb)  # complex128

    # Frequency grid and radial coordinate s = |k|
    Fx, Fy, Fz = _fourier_freq_grid(shape1, angpix)
    s_grid = np.sqrt(Fx * Fx + Fy * Fy + Fz * Fz).astype(np.float64, copy=False)

    absFT = np.abs(FT_comb)

    # Build Guinier arrays x = 1/d^2, y = ln(<|F|>) per shell
    x, y, w = _radial_bin_mean(absFT, s_grid, angpix, fit_minres, fit_maxres)
    sel = (w & np.isfinite(x) & np.isfinite(y))
    if np.count_nonzero(sel) < 3:
        raise ValueError(
            "normalize_amplitudes: not enough valid shells for Guinier fit."
        )

    x_fit = x[sel]
    y_fit = y[sel]

    # Linear least-squares fit: y = a + b x  → B = 4 b
    A = np.vstack([np.ones_like(x_fit), x_fit]).T
    a, b = np.linalg.lstsq(A, y_fit, rcond=None)[0]
    B = float(4.0 * b)

    # --- Apply B-factor to each half-map separately ---
    FT1 = np.fft.fftn(vol1.astype(np.float64, copy=False))
    FT2 = np.fft.fftn(vol2.astype(np.float64, copy=False))

    _apply_bfactor_inplace(FT1, s_grid, B)
    _apply_bfactor_inplace(FT2, s_grid, B)

    vol1_out = np.real(np.fft.ifftn(FT1)).astype(np.float32, copy=False)
    vol2_out = np.real(np.fft.ifftn(FT2)).astype(np.float32, copy=False)

    return vol1_out, vol2_out, B



def _apply_fsc_weighting_inplace(FT, s_grid, freqs, fsc_vals):
    """
    RELION-like FSC weighting:
      - Find the first shell where FSC < 1e-4  -> ires_max
      - Define a global hard cut at s_cut = freqs[ires_max]
      - For s <= s_cut, multiply by sqrt(2F/(1+F)) if F>0; for s > s_cut set to 0.
    """
    # 1) find ires_max (first near-zero)
    ires_max = 0
    for i, f in enumerate(fsc_vals):
        if f < 1e-4:
            break
        ires_max = i
    s_cut = float(freqs[ires_max])

    # 2) interpolate F(s) on the 3D grid
    F = np.interp(s_grid.ravel(), freqs, fsc_vals,
                  left=fsc_vals[0], right=fsc_vals[-1]).reshape(s_grid.shape).astype(np.float64, copy=False)

    # 3) apply RELION weight with a hard cut beyond s_cut
    inside = (s_grid <= s_cut)
    W = np.zeros_like(F, dtype=np.float64)
    pos = inside & (F > 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        W[pos] = np.sqrt((2.0 * F[pos]) / (1.0 + F[pos]))

    FT *= W.astype(FT.dtype, copy=False)

def score_halfmaps_from_files(
    hmap1_fname: str,
    hmap2_fname: str,
    report_path: str = None,
    verbose: bool = False,
):
    """Compute half-map quality metrics combining FSC, amplitude falloff and PSF isotropy.

    Parameters
    ----------
    hmap1_fname, hmap2_fname : str
        Paths to the two half-maps (MRC format). Both must have identical shape and angpix.
    report_path : str or None
        If not None, a short text report is written to this path.
    verbose : bool
        If True, print a text report to stdout. Default False.

    Returns
    -------
    global_score : float
        Global quality score in [0,1], where 0 = worst, 1 = best.
    res_143 : float
        FSC=0.143 resolution (Å) used as reference.
    fsc_score : float
        Information-weighted mean FSC up to the s_cut shell.
    amp_score : float
        Amplitude falloff score in [0,1] (Guinier linearity R²).
    iso_score : float
        PSF isotropy score in [0,1].
    freqs_valid : np.ndarray
        1D array of spatial frequencies (1/Å) used for the FSC curve (after masking/cleaning).
    fsc_valid : np.ndarray
        1D array of FSC values corresponding to freqs_valid.
    """
    # --- Load half-maps ---
    vol1, shape1, ang1 = _load_mrc_as_array(hmap1_fname)
    vol2, shape2, ang2 = _load_mrc_as_array(hmap2_fname)
    if shape1 != shape2:
        raise ValueError("Half-maps must have identical shapes.")
    if not np.isclose(ang1, ang2):
        raise ValueError("Half-maps must have identical pixel spacing (Å/px).")

    nz, ny, nx = shape1
    angpix = float(ang1)

    # --- FSC and resolution cutoff (FSC=0.143) ---
    freqs, fsc_vals = utils.compute_fsc_3d(vol1, vol2, angpix, n_bins=None)
    freqs = np.asarray(freqs, dtype=np.float64)
    fsc_vals = np.asarray(fsc_vals, dtype=np.float64)

    # Guard against pathological values
    fsc_vals = np.clip(fsc_vals, 0.0, 1.0)
    valid_f = np.isfinite(freqs) & np.isfinite(fsc_vals) & (freqs > 0.0)
    freqs_valid = freqs[valid_f]
    fsc_valid = fsc_vals[valid_f]

    if freqs_valid.size == 0:
        raise ValueError("No valid FSC samples; cannot score half-maps.")

    nyquist_s = 0.5 / angpix

    # Resolution from FSC=0.143 (if available)
    res_143 = None
    s_cut = None
    try:
        f_cross, res_val = utils.find_resolution_at_threshold(freqs, fsc_vals, 0.143)
        if res_val is not None and np.isfinite(res_val) and res_val > 0.0:
            res_143 = float(res_val)
            s_cut = 1.0 / res_143
    except Exception:
        res_143 = None
        s_cut = None

    if s_cut is None or not np.isfinite(s_cut) or s_cut <= 0.0:
        # Fallback: use last shell with FSC > 0.143, else max valid freq
        above = fsc_valid > 0.143
        if np.any(above):
            s_cut = float(freqs_valid[above][-1])
            res_143 = 1.0 / s_cut
        else:
            s_cut = float(freqs_valid[-1])
            res_143 = 1.0 / s_cut

    s_cut = float(min(s_cut, nyquist_s))

    # --- FSC-based score: information-weighted mean FSC up to s_cut ---
    shell_mask = (freqs_valid > 0.0) & (freqs_valid <= s_cut)
    s_shell = freqs_valid[shell_mask]
    f_shell = fsc_valid[shell_mask]
    if s_shell.size == 0:
        raise ValueError("No FSC samples within scoring range; check input half-maps.")

    weights = 4.0 * math.pi * (s_shell ** 2)
    fsc_score = float(np.sum(f_shell * weights) / np.sum(weights))

    # --- Build combined map and Fourier grid for amplitude / isotropy ---
    vol_comb = 0.5 * (vol1.astype(np.float64, copy=False) + vol2.astype(np.float64, copy=False))
    FT = np.fft.fftn(vol_comb)  # complex128

    Fx, Fy, Fz = _fourier_freq_grid(shape1, angpix)
    Fx = np.asarray(Fx, dtype=np.float64)
    Fy = np.asarray(Fy, dtype=np.float64)
    Fz = np.asarray(Fz, dtype=np.float64)

    # Broadcast to full 3D grids so boolean masks match in size
    s_grid = np.sqrt(Fx * Fx + Fy * Fy + Fz * Fz)
    Fx_full = np.broadcast_to(Fx, s_grid.shape)
    Fy_full = np.broadcast_to(Fy, s_grid.shape)
    Fz_full = np.broadcast_to(Fz, s_grid.shape)

    # --- Amplitude falloff score via Guinier linearity (R^2) ---
    absFT = np.abs(FT)
    # Choose a fit band around the FSC=0.143 resolution
    if res_143 is not None and np.isfinite(res_143) and res_143 > 0.0:
        fit_maxres = res_143
        fit_minres = max(res_143 * 2.0, res_143 + 2.0, 8.0)
    else:
        fit_maxres = max(3.0 * angpix, 3.0)
        fit_minres = 15.0

    x, y, w = _radial_bin_mean(absFT, s_grid, angpix, fit_minres, fit_maxres)
    sel = (w & np.isfinite(x) & np.isfinite(y))
    if np.count_nonzero(sel) >= 3:
        x_fit = x[sel].astype(np.float64)
        y_fit = y[sel].astype(np.float64)
        # Subtract means
        xm = x_fit - np.mean(x_fit)
        ym = y_fit - np.mean(y_fit)
        num = float(np.sum(xm * ym))
        den = float(math.sqrt(np.sum(xm * xm) * np.sum(ym * ym)))
        if den > 0.0:
            r = num / den
            amp_score = float(max(0.0, min(1.0, r * r)))
        else:
            amp_score = 0.0
    else:
        # Not enough points to fit; treat as neutral
        amp_score = 1.0

    # --- PSF isotropy score from directional second moments ---
    P = (absFT ** 2).astype(np.float64, copy=False)

    s_nz = s_grid > 0.0
    if np.any(s_nz):
        # High-frequency band near s_cut
        s_min_band = 0.7 * s_cut
        s_max_band = s_cut
        band = s_nz & (s_grid >= s_min_band) & (s_grid <= s_max_band)
        if np.count_nonzero(band) < 100:
            # Fallback: use upper third of spectrum if the band is too sparse
            s_thresh = 0.66 * nyquist_s
            band = s_nz & (s_grid >= s_thresh)
    else:
        band = None

    if band is None or np.count_nonzero(band) == 0:
        iso_score = 1.0
    else:
        s_band = s_grid[band]
        w_band = P[band]
        # Avoid zero or negative weights
        w_band = np.where(w_band > 0.0, w_band, 0.0)
        if np.sum(w_band) <= 0.0:
            iso_score = 1.0
        else:
            # Unit direction vectors n = k / |k|
            nx = Fx_full[band] / s_band
            ny_ = Fy_full[band] / s_band
            nz_ = Fz_full[band] / s_band

            # Weighted second-moment matrix C_ij = <n_i n_j>
            w_norm = w_band / np.sum(w_band)
            c_xx = float(np.sum(w_norm * nx * nx))
            c_yy = float(np.sum(w_norm * ny_ * ny_))
            c_zz = float(np.sum(w_norm * nz_ * nz_))
            c_xy = float(np.sum(w_norm * nx * ny_))
            c_xz = float(np.sum(w_norm * nx * nz_))
            c_yz = float(np.sum(w_norm * ny_ * nz_))

            C = np.array([[c_xx, c_xy, c_xz],
                          [c_xy, c_yy, c_yz],
                          [c_xz, c_yz, c_zz]], dtype=np.float64)

            # For an isotropic PSF, C ~ (1/3) * I; eigenvalues are all 1/3
            evals = np.linalg.eigvalsh(C)
            lam_min = float(np.min(evals.real))
            iso_score = float(max(0.0, min(1.0, 3.0 * lam_min)))

    # --- Combine into a single scalar via geometric mean ---
    eps = 1e-8
    fsc_c = float(max(eps, min(1.0, fsc_score)))
    amp_c = float(max(eps, min(1.0, amp_score)))
    iso_c = float(max(eps, min(1.0, iso_score)))
    global_score = float((fsc_c * amp_c * iso_c) ** (1.0 / 3.0))

    # --- Prepare a human-readable report ---
    lines = []
    lines.append("Half-map quality scoring (FSC + amplitudes + PSF isotropy)")
    lines.append(f"  Half-map 1: {hmap1_fname}")
    lines.append(f"  Half-map 2: {hmap2_fname}")
    lines.append(f"  Pixel size (Å/px): {angpix:.4f}")
    lines.append("")
    lines.append(f"  FSC=0.143 resolution (Å): {res_143:.3f}")
    lines.append(f"  FSC-based score:         {fsc_score:.4f}")
    lines.append(f"  Amplitude falloff score: {amp_score:.4f}")
    lines.append(f"  PSF isotropy score:      {iso_score:.4f}")
    lines.append("")
    lines.append(f"  Global half-map quality score: {global_score:.4f}")
    report_text = "\n".join(lines)

    if verbose:
        print(report_text)

    if report_path is not None:
        try:
            with open(report_path, "w", encoding="utf-8") as fh:
                fh.write(report_text + "\n")
        except Exception as e:
            print(f"[scoreHM] WARNING: could not write report to '{report_path}': {e}")

    return global_score, res_143, fsc_score, amp_score, iso_score, freqs_valid, fsc_valid



def _raised_cosine_lowpass_inplace(FT, s_grid, cutoff_A, angpix, edge_width_shells=2):
    """
    Raised-cosine low-pass with an edge width expressed in *shells*.
    """
    if cutoff_A <= 0:
        return
    s_c = 1.0 / cutoff_A
    s_max = 0.5 / angpix
    if s_c >= s_max:
        return

    # shell spacing ~ 1/(angpix * max_dim) in Å^-1
    max_dim = max(s_grid.shape)
    shell_ds = 1.0 / (angpix * float(max_dim))
    delta = float(edge_width_shells) * shell_ds

    s = s_grid
    W = np.ones_like(s, dtype=FT.real.dtype)
    W[s > (s_c + delta)] = 0.0
    in_edge = (s >= (s_c - delta)) & (s <= (s_c + delta))
    t = (s[in_edge] - (s_c - delta)) / (2.0 * delta)
    W[in_edge] = 0.5 * (1.0 + np.cos(np.pi * t))

    FT *= W.astype(FT.dtype, copy=False)



def bfactor_utils(args):
    # --- Load input map (file or in-memory) ---
    if getattr(args, "i_array", None) is not None:
        vol = np.asarray(args.i_array, dtype=np.float64, order="C")
        shape = tuple(args.shape)
        angpix = float(args.angpix)
        spacing = float(args.spacing)
    else:
        vol, shape, angpix = _load_mrc_as_array(args.i)
        vol = vol.astype(np.float64, copy=False)
        spacing = round(janas_core.spacingMRC(args.i), 4)

    nz, ny, nx = shape

    # --- FFT (float64 -> complex128) and frequency grid ---
    FT = np.fft.fftn(vol)  # complex128
    Fx, Fy, Fz = _fourier_freq_grid(shape, angpix)
    s_grid = np.sqrt(Fx*Fx + Fy*Fy + Fz*Fz).astype(np.float64, copy=False)  # 1/Å

    # --- Optional FSC weighting (from half-maps) ---
    freqs = fsc_vals = None
    if getattr(args, "hmaps", None):
        h1, h2 = args.hmaps
        vol1, shape1, ang1 = _load_mrc_as_array(h1)
        vol2, shape2, ang2 = _load_mrc_as_array(h2)
        if shape1 != shape2 or shape1 != shape:
            raise ValueError("Half-map shapes must match each other and the input map.")
        if not (np.isclose(ang1, ang2) and np.isclose(ang1, angpix)):
            raise ValueError("Å/px must match between input map and half-maps.")
        freqs, fsc_vals = utils.compute_fsc_3d(vol1, vol2, angpix, n_bins=None)
        _apply_fsc_weighting_inplace(FT, s_grid, freqs, fsc_vals)

    # --- B-factor (auto or ad-hoc) ---
    if getattr(args, "auto", False):
        absFT = np.abs(FT)
        x, y, w = _radial_bin_mean(absFT, s_grid, angpix, args.fit_minres, args.fit_maxres)
        sel = (w & np.isfinite(y))
        if not np.any(sel):
            raise ValueError("No valid shells in the selected Guinier fit range.")
        A = np.vstack([np.ones_like(x[sel]), x[sel]]).T
        a, b = np.linalg.lstsq(A, y[sel], rcond=None)[0]
        B = float(4.0 * b)
        print(f"[bfactor] auto-fit slope={b:.3f} -> B={B:.1f} Å^2")
    else:
        B = float(getattr(args, "B", 0.0))
        if B > 0:
            print("[bfactor] WARNING: positive B will damp high frequencies.")

    # --- Apply B ---
    _apply_bfactor_inplace(FT, s_grid, B)

    # --- Optional low-pass (RELION-style negative flag for FSC=0.143) ---
    if getattr(args, "lowpass", 0.0) != 0.0:
        cutoff_A = args.lowpass
        if cutoff_A < 0.0 and (freqs is not None and fsc_vals is not None):
            f_cross, res = utils.find_resolution_at_threshold(freqs, fsc_vals, 0.143)
            if res is not None:
                cutoff_A = float(res)
                print(f"[bfactor] Using FSC=0.143 crossing: {cutoff_A:.3f} Å as low-pass")
            else:
                cutoff_A = 0.0
        if cutoff_A > 0.0:
            _raised_cosine_lowpass_inplace(FT, s_grid, cutoff_A, angpix, edge_width_shells=2)

    # --- iFFT and write ---
    out_vol = np.real(np.fft.ifftn(FT))
    janas_core.WriteMRC(out_vol.flatten().tolist(), args.o, nx, ny, nz, spacing)
    print(f"[bfactor] Wrote: {args.o} (B={B:.1f} Å^2)")


def gaussian_blur_array(vol, angpix, sigmaA, device='cpu', gpu_index=None):
    """
    Gaussian blur a 3D volume with σ in Å on CPU or CUDA.

    vol       : np.ndarray or torch.Tensor, shape (Z, Y, X)
    angpix    : float, Å/pixel
    sigmaA    : float, Gaussian σ in Å
    device    : 'cpu' or 'cuda'
    gpu_index : optional CUDA device index
    returns   : np.ndarray (float32) on CPU
    """
    import numpy as np
    import math

    if sigmaA is None or sigmaA <= 0:
        return np.asarray(vol, dtype=np.float32, order='C')

    if device == 'cuda':
        try:
            import torch
        except Exception as e:
            raise RuntimeError('CUDA blur requested but torch is not available') from e

        dev = torch.device(f'cuda:{gpu_index}' if gpu_index is not None else 'cuda')
        with torch.cuda.device(dev):
            t = torch.as_tensor(vol, dtype=torch.float32, device=dev)  # (Z,Y,X)
            FT = torch.fft.fftn(t)

            nz, ny, nx = t.shape
            fx = torch.fft.fftfreq(nx, d=float(angpix), device=dev)
            fy = torch.fft.fftfreq(ny, d=float(angpix), device=dev)
            fz = torch.fft.fftfreq(nz, d=float(angpix), device=dev)
            Fx = fx.view(1, 1, -1)
            Fy = fy.view(1, -1, 1)
            Fz = fz.view(-1, 1, 1)
            s2 = Fx * Fx + Fy * Fy + Fz * Fz

            factor = -2.0 * (math.pi ** 2) * (float(sigmaA) ** 2)
            W = torch.exp(factor * s2)
            FT = FT * W

            out = torch.fft.ifftn(FT).real
            return out.detach().cpu().numpy().astype(np.float32, copy=False)

    # CPU path (NumPy)
    vol64 = np.asarray(vol, dtype=np.float64, order='C')
    FT = np.fft.fftn(vol64)
    # build frequency grid in Å^-1
    nx = vol64.shape[2]; ny = vol64.shape[1]; nz = vol64.shape[0]
    fx = np.fft.fftfreq(nx, d=float(angpix)).reshape(1, 1, nx)
    fy = np.fft.fftfreq(ny, d=float(angpix)).reshape(1, ny, 1)
    fz = np.fft.fftfreq(nz, d=float(angpix)).reshape(nz, 1, 1)
    s2 = fx * fx + fy * fy + fz * fz
    factor = -2.0 * (np.pi ** 2) * (float(sigmaA) ** 2)
    W = np.exp(factor * s2, dtype=np.float64)
    FT *= W.astype(FT.dtype, copy=False)
    out = np.real(np.fft.ifftn(FT)).astype(np.float32, copy=False)
    return out



def gaussian_blur_utils(args):
    """
    Gaussian blur a 3D map with σ in Å.

    Inputs (SimpleNamespace or argparse args):
      - args.i:  input .mrc
      - args.o:  output .mrc
      - args.sigmaA: Gaussian σ in Å (float > 0)

    The blur is applied in Fourier domain:
        FT_out = FT_in * exp(-2 * pi^2 * sigmaA^2 * s^2),
    where s is spatial frequency in Å^-1 computed from np.fft.fftfreq with d=angpix.
    """
    in_mrc = args.i
    out_mrc = args.o
    sigmaA = float(args.sigmaA)
    if sigmaA <= 0.0:
        raise ValueError(f"sigmaA must be > 0 (Å); got {sigmaA}")

    # Load volume and metadata
    vol, shape, angpix = _load_mrc_as_array(in_mrc)  # vol is float32 [Z,Y,X]
    vol = vol.astype(np.float64, copy=False)
    nz, ny, nx = shape
    spacing = round(janas_core.spacingMRC(in_mrc), 4)

    # FFT and frequency grid (s in Å^-1)
    FT = np.fft.fftn(vol)  # complex128
    Fx, Fy, Fz = _fourier_freq_grid(shape, angpix)
    s_grid = np.sqrt(Fx*Fx + Fy*Fy + Fz*Fz).astype(np.float64, copy=False)

    # Fourier-domain Gaussian (σ in Å): exp(-2 π^2 σ^2 s^2)
    factor = -2.0 * (np.pi ** 2) * (sigmaA ** 2)
    W = np.exp(factor * (s_grid ** 2))
    FT *= W.astype(FT.dtype, copy=False)

    # iFFT and write
    out_vol = np.real(np.fft.ifftn(FT))
    janas_core.WriteMRC(out_vol.flatten().tolist(), out_mrc, nx, ny, nz, spacing)
    print(f"[blur] σ={sigmaA:.3f} Å  ->  {out_mrc}")


def normalize_utils(args):
    """
    Intensity normalisation for a 3D MRC map.

    Inputs (SimpleNamespace or argparse args):
      - args.i:      input .mrc
      - args.o:      output .mrc
      - args.method: 'zscore' (default) or 'minmax'

    Behaviour:
      zscore:  out = (vol - mean) / std          (global mean/std over all voxels)
      minmax:  out = (vol - min) / (max - min)   (global min/max; if flat -> zeros)
    """
    in_mrc = args.i
    out_mrc = args.o
    method = getattr(args, "method", "zscore").lower()

    # Load volume and metadata
    vol, shape, _ = _load_mrc_as_array(in_mrc)  # vol: float32 (Z,Y,X)
    vol = vol.astype(np.float64, copy=False)
    nz, ny, nx = shape
    spacing = round(janas_core.spacingMRC(in_mrc), 4)

    if method == "zscore":
        mu = float(np.mean(vol))
        sigma = float(np.std(vol))
        if not np.isfinite(sigma) or sigma == 0.0:
            # Degenerate case: constant map -> write zeros
            out_vol = np.zeros_like(vol, dtype=np.float64)
        else:
            out_vol = (vol - mu) / sigma
    elif method == "minmax":
        vmin = float(np.min(vol))
        vmax = float(np.max(vol))
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax == vmin:
            out_vol = np.zeros_like(vol, dtype=np.float64)
        else:
            out_vol = (vol - vmin) / (vmax - vmin)
    else:
        raise ValueError(f"Unknown method '{method}'. Use 'zscore' or 'minmax'.")

    # Write output
    janas_core.WriteMRC(out_vol.astype(np.float32, copy=False).flatten().tolist(),
                          out_mrc, nx, ny, nz, spacing)
    print(f"[normalize] method={method} -> {out_mrc}")



