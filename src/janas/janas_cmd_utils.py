#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# File: janas_cmd_utils.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology


"""
Module: janas_utils.py

Utility commands for the JANAS pipeline:
- maskedCrop: automatically crop volumes based on a binary mask
- equalize_images: match Fourier amplitude spectra across multiple volumes
- scores_to_csv: merge multiple SCI score STAR files into CSV + assign class labels
- randomize_halves: shuffle particles into two half‐stacks
- extract_particles_from_label_value: select particles by metadata label
"""

# Standard library
import argparse
import os
import os.path
from os import PathLike, makedirs, path
import glob

# Third-party
import numpy as np
from numpy.fft import fftn, ifftn, fftfreq
import pandas as pd
import matplotlib.pyplot as plt 
import re
import tempfile
from typing import Tuple, Dict, Any, Optional
from types import SimpleNamespace
from pathlib import Path
import json
import os

import sys
import time


# Local
import janas.janas_core as janas_core
from janas import starHandler
from janas import utils
from janas.version import get_version
from janas import locres_utils
from janas import projector_utils
from janas import janas_mapProcess
from janas import janas_alternative_selector as alt_selector

janas_parser = argparse.ArgumentParser(
    prog="janas_utils",
    usage="%(prog)s [command] [arguments]",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
janas_parser.add_argument(
    "-V", "--version",
    action="version",
    version=get_version(),
    help="show program’s version number and exit"
)
command = janas_parser.add_subparsers(dest="command")


# ---------------------------------
# compare2D (direct 2D vs 2D scoring)
janas_compare2D = command.add_parser(
    "compare2D",
    description="Compare two 2D MRC images with SCIFIR/LoG (optionally SCI) under an optional mask.",
    help="compare two 2D images (no stack, no CTF)"
)
janas_compare2D.add_argument("I", type=str, help="input image I (.mrc; may be 2D or 3D with nz=1)")
janas_compare2D.add_argument("R", type=str, help="reference image R (.mrc; same size as I)")
janas_compare2D.add_argument("--mask", default=None,
                               help="mask (.mrc). Use 'none' for all-ones mask, or 'round' or 'round:RAD' for a centred round mask.")
janas_compare2D.add_argument("--score", default="SCI",
                               choices=["SCI", "SCI_FT", "LoG","LoG_Equalized","CC"],
                               help="scoring kernel")
janas_compare2D.add_argument("--sigma", type=float, default=1.0, help="sigma in pixels")
janas_compare2D.add_argument("--json", dest="json_out", default=None, help="optional JSON output path")

def _read_mrc_as_2d(fname: str) -> Tuple[np.ndarray, Tuple[int, int, int], float]:
    """
    Read an MRC using janas_core and return:
      img2d as float32 with shape (ny, nx),
      (nx, ny, nz) header tuple,
      spacing (float, Å/px).
    Accepts true-2D (if your reader supports it) or 3D with nz==1.
    """
    if not os.path.exists(fname):
        raise FileNotFoundError(f"File not found: {fname}")

    nx, ny, nz = janas_core.sizeMRC(fname)  # header order in this codebase
    spacing = float(round(janas_core.spacingMRC(fname), 6))

    flat = np.array(janas_core.ReadMRC(fname), dtype=np.float32)

    # Critical: MRC is X-fast -> reshape as (nz, ny, nx) so X is last in NumPy
    vol = flat.reshape((nz, ny, nx), order="C")

    if nz == 1:
        img2d = vol[0]  # (ny, nx)
        return img2d, (nx, ny, nz), spacing

    # If it's not nz==1, treat as an error for compare2D
    raise ValueError(f"compare2D expects a 2D image (nz==1). Got (nx,ny,nz)=({nx},{ny},{nz}) for {fname}")


def _make_round_mask_2d(nx: int, ny: int, radius: float = None) -> np.ndarray:
    """
    Return a centred circular binary mask (ny,nx).
    If radius is None, use maximum that fits.
    """
    cx = (nx - 1) / 2.0
    cy = (ny - 1) / 2.0
    if radius is None:
        radius = float(min(cx, cy))
    yy = np.arange(ny, dtype=np.float32)[:, None]
    xx = np.arange(nx, dtype=np.float32)[None, :]
    d2 = (xx - cx) ** 2 + (yy - cy) ** 2
    return (d2 <= (radius ** 2)).astype(np.float32)


def compare2D_cli(args):
    I, (nxI, nyI, nzI), apI = _read_mrc_as_2d(args.I)
    R, (nxR, nyR, nzR), apR = _read_mrc_as_2d(args.R)

    if (nxI, nyI, nzI) != (nxR, nyR, nzR):
        raise ValueError(f"I and R shape mismatch: I={nxI,nyI,nzI} vs R={nxR,nyR,nzR}")

    # Mask handling
    if args.mask is None:
        # default: all ones
        M = np.ones((nyI, nxI), dtype=np.float32)
    else:
        m = str(args.mask).strip()
        if m.lower() == "none":
            M = np.ones((nyI, nxI), dtype=np.float32)
        elif m.lower().startswith("round"):
            # allow: "round" or "round:30"
            rad = None
            if ":" in m:
                _, rad_s = m.split(":", 1)
                rad = float(rad_s)
            M = _make_round_mask_2d(nxI, nyI, radius=rad)
        else:
            M2, (nxM, nyM, nzM), _ = _read_mrc_as_2d(m)
            if (nxM, nyM, nzM) != (nxI, nyI, nzI):
                raise ValueError(f"Mask shape mismatch: mask={nxM,nyM,nzM} vs I={nxI,nyI,nzI}")
            # For masks, treat >0 as inside; keep as float32 weights
            M = (M2 > 0.0).astype(np.float32)

    sigma = float(args.sigma)

    score_name = args.score.upper()
    if score_name == "SCI_FT":
        score = alt_selector.sci_fir_gaussian_derivative_2d(I, R, M, sigma_pixels=sigma)
    elif score_name == "LOG":
        score = alt_selector.log_score_2d(I, R, M, sigma_pixels=sigma)
    elif score_name == "LOG_EQUALIZED":
        score = alt_selector.log_score_2d_equalized(I, R, M, sigma_pixels=sigma)
    elif score_name == "SCI":
        score = alt_selector.sci_2d(I, R, M, sigmaBlur=str(args.sigma))
    elif score_name == "CC":
        score = alt_selector.CC_2d(I, R, M, sigma_pixels=sigma)
    else:
        raise RuntimeError(f"Unknown score: {args.score}")

    out = {
        "I": args.I,
        "R": args.R,
        "mask": args.mask,
        "score": args.score,
        "sigma_pixels": sigma,
        "value": float(score),
        "shape_nx_ny_nz": [nxI, nyI, nzI],
        "apix": apI,  # informative; compare2D does not require it
    }

    print(json.dumps(out, indent=2))

    if args.json_out:
        with open(args.json_out, "w") as f:
            json.dump(out, f, indent=2)


# ---------------------------------
# compare_maps (direct 3D vs 3D CC scoring)
janas_compare_maps = command.add_parser(
    "compare_maps",
    description=(
        "Compare 3D MRC maps using CC, optionally within a mask. "
        "With --target_map, each target map is compared against --maps. "
        "Without --target_map, all maps passed to --maps are mutually compared and written as a square matrix."
    ),
    help="compare 3D maps using CC"
)
janas_compare_maps.add_argument(
    "--target_map",
    required=False,
    nargs="+",
    type=str,
    help=(
        "optional target 3D maps (.mrc). Wildcards such as * and ? are accepted. "
        "If omitted, compare_maps runs in mutual all-vs-all mode on --maps."
    )
)
janas_compare_maps.add_argument(
    "--maps",
    required=True,
    nargs="+",
    type=str,
    help="one or more 3D maps (.mrc). Wildcards such as * and ? are accepted."
)
janas_compare_maps.add_argument(
    "--mask",
    default=None,
    help="optional mask (.mrc). Only voxels with mask > 0 are used."
)
janas_compare_maps.add_argument(
    "--o",
    default=None,
    type=str,
    help=(
        "optional CSV output path. In mutual mode this writes the square similarity matrix; "
        "in target mode it writes a target-vs-map table."
    )
)


def _read_mrc_as_3d(fname: str) -> Tuple[np.ndarray, Tuple[int, int, int], float]:
    """
    Read an MRC using janas_core and return:
      vol as float32 with shape (nz, ny, nx),
      (nx, ny, nz) header tuple,
      spacing (float, Å/px).
    """
    if not os.path.exists(fname):
        raise FileNotFoundError(f"File not found: {fname}")

    nx, ny, nz = janas_core.sizeMRC(fname)
    spacing = float(round(janas_core.spacingMRC(fname), 6))

    flat = np.array(janas_core.ReadMRC(fname), dtype=np.float32)

    # MRC is X-fast -> reshape as (nz, ny, nx) so X is last in NumPy
    vol = flat.reshape((nz, ny, nx), order="C")
    return vol, (nx, ny, nz), spacing


def _expand_map_patterns(items):
    """
    Expand wildcard patterns (*, ?, []) in a list of filenames/patterns.
    If a token has no wildcard, keep it as-is.
    Remove duplicates while preserving order.
    """
    collected = []

    if items is None:
        return []

    for item in items:
        if glob.has_magic(item):
            matches = sorted(glob.glob(item))
            if not matches:
                raise FileNotFoundError(f"No files matched pattern: {item}")
            collected.extend(matches)
        else:
            if not os.path.exists(item):
                raise FileNotFoundError(f"File not found: {item}")
            collected.append(item)

    seen = set()
    out = []
    for f in collected:
        if f in seen:
            continue
        seen.add(f)
        out.append(f)

    return out


def _load_optional_mask_3d(mask_arg: Optional[str], expected_shape: Tuple[int, int, int]) -> Optional[np.ndarray]:
    """
    Load and binarise a 3D mask for comparisons.
    Returns a boolean array with shape (nz, ny, nx), or None.
    """
    if mask_arg is None:
        return None

    m = str(mask_arg).strip()
    if m == "" or m.lower() == "none":
        return None

    mask_vol, mask_shape, _ = _read_mrc_as_3d(m)
    if mask_shape != expected_shape:
        raise ValueError(f"Mask shape mismatch: mask={mask_shape} vs maps={expected_shape}")

    inside = mask_vol > 0.0
    if not np.any(inside):
        raise ValueError("Mask contains no voxels > 0.")

    return inside


def _normalised_cc_vector_3d(vol: np.ndarray, mask: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Convert a 3D map into a zero-mean, unit-norm 1D vector for fast repeated CC.
    The returned vector is float32 to reduce memory footprint when many maps are compared.
    """
    if mask is None:
        vec = vol.ravel().astype(np.float64, copy=False)
    else:
        vec = vol[mask].astype(np.float64, copy=False)

    vec = vec - np.mean(vec)
    norm = np.sqrt(np.dot(vec, vec))
    if norm <= 0.0:
        return np.zeros(vec.shape, dtype=np.float32)

    return (vec / norm).astype(np.float32, copy=False)


def _cc_from_normalised_vectors(vec_a: np.ndarray, vec_b: np.ndarray) -> float:
    """CC between two zero-mean, unit-norm vectors."""
    return float(np.dot(vec_a, vec_b))


def compare_maps_cli(args):
    map_list = _expand_map_patterns(args.maps)
    target_list = _expand_map_patterns(args.target_map)

    if len(map_list) == 0:
        raise ValueError("No maps were provided to --maps.")

    # Load all comparison maps once, because they are reused in either mode.
    loaded_maps = []
    reference_shape = None
    for fname in map_list:
        vol, shape, _ = _read_mrc_as_3d(fname)
        if reference_shape is None:
            reference_shape = shape
        elif shape != reference_shape:
            raise ValueError(
                f"Map shape mismatch in --maps: {fname} has shape {shape}, expected {reference_shape}"
            )
        loaded_maps.append((fname, vol))

    mask_bool = _load_optional_mask_3d(args.mask, reference_shape)

    # Mutual all-vs-all mode: compare only --maps.
    if len(target_list) == 0:
        n_maps = len(loaded_maps)
        names = [fname for fname, _ in loaded_maps]
        norm_vectors = [_normalised_cc_vector_3d(vol, mask_bool) for _, vol in loaded_maps]

        matrix = np.eye(n_maps, dtype=np.float32)
        for i in range(n_maps):
            for j in range(i + 1, n_maps):
                cc = _cc_from_normalised_vectors(norm_vectors[i], norm_vectors[j])
                matrix[i, j] = cc
                matrix[j, i] = cc

        df = pd.DataFrame(matrix, index=names, columns=names)
        df.index.name = "maps"

        if args.o:
            df.to_csv(args.o, float_format="%.6f")
        print(df.to_csv(float_format="%.6f"))
        return

    # Target-vs-map mode: keep previous behaviour, but optionally write CSV.
    records = []
    norm_map_vectors = [_normalised_cc_vector_3d(vol, mask_bool) for _, vol in loaded_maps]

    for target_fname in target_list:
        target, target_shape, _ = _read_mrc_as_3d(target_fname)
        if target_shape != reference_shape:
            raise ValueError(
                f"Target map shape mismatch: {target_fname} has shape {target_shape}, expected {reference_shape}"
            )

        target_vec = _normalised_cc_vector_3d(target, mask_bool)
        results = []
        for (fname, _), map_vec in zip(loaded_maps, norm_map_vectors):
            cc = _cc_from_normalised_vectors(target_vec, map_vec)
            results.append((fname, cc))
            records.append({"target_map": target_fname, "map": fname, "cc": cc})

        results.sort(key=lambda x: (-x[1], x[0]))

        print(f"# TARGET_MAP\t{target_fname}")
        for fname, cc in results:
            print(f"{fname}\t{cc:.6f}")
        if target_fname != target_list[-1]:
            print()

    if args.o:
        pd.DataFrame.from_records(records).to_csv(args.o, index=False, float_format="%.6f")
# ---------------------------
# Automatic sigma estimate from half-maps (mask-aware FSC)
janas_sigma_estimate = command.add_parser(
    "sigma_estimate",
    description=(
        "Estimate the Gaussian-derivative sigma (in pixels) from a half-map pair via FSC. "
        "This calls janas_alternative_selector.auto_sigma_from_halfmaps()."
    ),
    help="estimate sigma (pixels) from half-maps via FSC",
)
janas_sigma_estimate.add_argument("halfmap1", type=str, help="first half-map (.mrc)")
janas_sigma_estimate.add_argument("halfmap2", type=str, help="second half-map (.mrc)")
janas_sigma_estimate.add_argument(
    "mask",
    nargs="?",
    default=None,
    help="optional mask (.mrc). Omit or pass 'none' to disable masking."
)
janas_sigma_estimate.add_argument(
    "--gamma",
    dest="gamma_at_f0143",
    type=float,
    default=0.5,
    help="gamma value used at FSC=0.143 (see auto_sigma_from_halfmaps)."
)
janas_sigma_estimate.add_argument(
    "--scale",
    dest="sigma_scale",
    type=float,
    default=2.0,
    help="multiplicative scale factor applied to the base sigma."
)
janas_sigma_estimate.add_argument(
    "--sigma-min-px",
    dest="sigma_min_px",
    type=float,
    default=0.25,
    help="minimum base sigma in pixels before scaling (clamp)."
)
janas_sigma_estimate.add_argument(
    "--sigma-max-px",
    dest="sigma_max_px",
    type=float,
    default=8.0,
    help="maximum base sigma in pixels before scaling (clamp)."
)
janas_sigma_estimate.add_argument(
    "--mask-threshold",
    type=float,
    default=0.5,
    help="mask binarisation threshold (voxels >= threshold are inside)."
)
janas_sigma_estimate.add_argument(
    "--mask-soft-edge-add",
    type=float,
    default=0.0,
    help="add a soft rim outside the binary mask (in pixels)."
)
janas_sigma_estimate.add_argument(
    "--full",
    action="store_true",
    help="print a full JSON report (in addition to the sigma value)."
)
janas_sigma_estimate.add_argument(
    "--json",
    dest="json_out",
    default=None,
    help="optional path to write the full JSON report."
)

def sigma_estimate_utils(args):
    h1 = str(args.halfmap1)
    h2 = str(args.halfmap2)

    if not os.path.exists(h1):
        raise FileNotFoundError(f"Half-map 1 not found: {h1}")
    if not os.path.exists(h2):
        raise FileNotFoundError(f"Half-map 2 not found: {h2}")

    # Optional mask
    mask = args.mask
    if mask is not None:
        mask = str(mask).strip()
        if mask == "" or mask.lower() == "none":
            mask = None
    if mask is not None and not os.path.exists(mask):
        raise FileNotFoundError(f"Mask not found: {mask}")

    # Validate shapes from headers (nx, ny, nz)
    nx1, ny1, nz1 = janas_core.sizeMRC(h1)
    nx2, ny2, nz2 = janas_core.sizeMRC(h2)
    if (nx1, ny1, nz1) != (nx2, ny2, nz2):
        raise ValueError(f"Half-map shape mismatch: {(nx1, ny1, nz1)} vs {(nx2, ny2, nz2)}")
    if mask is not None:
        nxm, nym, nzm = janas_core.sizeMRC(mask)
        if (nxm, nym, nzm) != (nx1, ny1, nz1):
            raise ValueError(f"Mask shape mismatch: {(nxm, nym, nzm)} vs {(nx1, ny1, nz1)}")

    # Pixel size (Å/px): normalise utils.get_MRC_map_pixel_spacing() return type
    def _safe_angpix(spacing):
        try:
            return float(spacing[0])
        except (TypeError, IndexError):
            return float(spacing)

    ap1 = _safe_angpix(utils.get_MRC_map_pixel_spacing(h1))
    ap2 = _safe_angpix(utils.get_MRC_map_pixel_spacing(h2))
    if not np.isclose(ap1, ap2, rtol=0, atol=1e-6):
        raise ValueError(f"Å/px mismatch: {ap1} vs {ap2}")

    sigma_px, meta = alt_selector.auto_sigma_from_halfmaps(
        map_h1_paths=[h1],
        map_h2_paths=[h2],
        mask_path=mask,
        apix=float(ap1),
        gamma_at_f0143=float(args.gamma_at_f0143),
        sigma_scale=float(args.sigma_scale),
        sigma_minmax_px=(float(args.sigma_min_px), float(args.sigma_max_px)),
        mask_threshold=float(args.mask_threshold),
        mask_soft_edge_add=float(args.mask_soft_edge_add),
    )

    # Default: print sigma in pixels as a single number (copy/paste into --sigma).
    print(f"{float(sigma_px):.6f}")

    if args.full or args.json_out:
        report = {
            "halfmap1": h1,
            "halfmap2": h2,
            "mask": mask,
            "shape_nx_ny_nz": [int(nx1), int(ny1), int(nz1)],
            "apix": float(ap1),
            "sigma_px": float(sigma_px),
            "sigma_A": float(sigma_px) * float(ap1),
            "meta": meta,
        }
        if args.full:
            print(json.dumps(report, indent=2))
        if args.json_out:
            with open(args.json_out, "w") as f:
                json.dump(report, f, indent=2)


# ---------------------------
# SSNR from half-maps (global, no mask)
janas_ssnr = command.add_parser(
    "ssnr",
    description="Estimate per-shell SSNR from half-maps via gold-standard FSC.",
    help="Compute SSNR(f) = FSC/(1-FSC) per frequency shell (no mask).",
)
janas_ssnr.add_argument(
    "--hmaps", required=True, nargs=2, action="append",
    metavar=("HALF1.mrc", "HALF2.mrc"),
    help="Pairs of half-maps (can be repeated)."
)
janas_ssnr.add_argument(
    "--bins", type=int, default=None, help="Number of frequency bins (default auto)."
)
janas_ssnr.add_argument(
    "--o", "--output_prefix", dest="output_prefix", default=None,
    help="If set, write one CSV per pair as <prefix>_setN.csv"
)
janas_ssnr.add_argument(
    "--no-plot", action="store_true", help="Disable plotting."
)

def ssnr_utils(args):
    import os
    curves = []
    for idx_pair, pair in enumerate(args.hmaps, start=1):
        f1, f2 = pair
        if not (os.path.exists(f1) and os.path.exists(f2)):
            raise FileNotFoundError(f"Set {idx_pair}: missing half-map(s)")
        freqs, ssnr, ang = locres_utils.compute_ssnr_from_paths(f1, f2, n_bins=args.bins)
        if args.output_prefix:
            out_csv = f"{args.output_prefix}_set{idx_pair}.csv"
            np.savetxt(
                out_csv,
                np.column_stack((freqs, ssnr)),
                header="Frequency(1/Å),SSNR",
                fmt="%.6f,%.6f",
                delimiter=","
            )
            print(f"Saved SSNR curve to {out_csv}")
        curves.append({"freqs": freqs, "ssnr": ssnr,
                       "label": f"{os.path.basename(f1)} | {os.path.basename(f2)}"})
    if not args.no_plot and curves:
        plt.figure()
        for c in curves:
            plt.plot(c["freqs"], c["ssnr"], lw=2.0, label=c["label"])
        plt.xlabel("Spatial frequency (Å⁻¹)")
        plt.ylabel("SSNR")
        plt.grid(True); plt.legend(loc="upper right", fontsize=8)
        plt.xlim(left=0.0)
        plt.tight_layout(); plt.show()


#################################
# Half-map quality analysis (global scalar 0–1)
janas_quality = command.add_parser(
    "quality_analysis",
    description=(
        "Compute a global half-map quality score (0–1) "
        "combining FSC, amplitude falloff and PSF isotropy."
    ),
    help="Score half-map quality via FSC, amplitudes and PSF isotropy.",
)

janas_quality.add_argument(
    "--hmaps",
    required=True,
    nargs=2,
    metavar=("HALF1.mrc", "HALF2.mrc"),
    help="Paths to the two half-maps (MRC format).",
)

janas_quality.add_argument(
    "--report",
    "-r",
    dest="report_path",
    default=None,
    help="Optional path for a text report (same content as printed to screen).",
)
def quality_analysis(args):
    """CLI wrapper around janas_mapProcess.score_halfmaps_from_files.

    Expects:
        args.hmaps: list/tuple of two paths to half-maps
        args.report_path: optional output report filename
    """
    import os

    h1, h2 = args.hmaps
    if not os.path.exists(h1):
        raise FileNotFoundError(f"Half-map 1 not found: {h1}")
    if not os.path.exists(h2):
        raise FileNotFoundError(f"Half-map 2 not found: {h2}")

    global_score, *_ = janas_mapProcess.score_halfmaps_from_files(
        hmap1_fname=h1,
        hmap2_fname=h2,
        report_path=args.report_path,
        verbose=True,
    )
    print(
        "[quality_analysis] "
        f"FSC_0.143={res_143:.3f} Å, "
        f"FSC_score={fsc_score:.4f}, "
        f"amp_score={amp_score:.4f}, "
        f"PSF_score={iso_score:.4f}, "
        f"global_score={global_score:.4f}"
    )




#################################
## split_mask
janas_split_mask = command.add_parser(
    "split_mask",
    description=(
        "Split the values of a map inside a mask into low/high binary masks using an Otsu threshold. "
        "Optionally Gaussian-blur the map before thresholding."
    ),
    help="split a masked map into low/high masks using Otsu thresholding",
)
janas_split_mask.add_argument(
    "--map",
    required=True,
    type=str,
    help="input map (.mrc) used to compute the thresholded split",
)
janas_split_mask.add_argument(
    "--mask",
    required=True,
    type=str,
    help="input mask (.mrc) with the same size as --map; voxels > 0 define the valid region",
)
janas_split_mask.add_argument(
    "--blur",
    required=False,
    type=float,
    default=0.0,
    help="optional Gaussian sigma in voxels applied to the map before thresholding (default 0)",
)
janas_split_mask.add_argument(
    "--o_below",
    required=False,
    type=str,
    default=None,
    help="output mask (.mrc) for values <= Otsu threshold inside the input mask",
)
janas_split_mask.add_argument(
    "--o_above",
    required=False,
    type=str,
    default=None,
    help="output mask (.mrc) for values > Otsu threshold inside the input mask",
)


def _otsu_threshold_from_values(values: np.ndarray, n_bins: int = 256) -> float:
    vals = np.asarray(values, dtype=np.float64)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        raise ValueError("No finite values are available to compute the Otsu threshold.")

    vmin = float(np.min(vals))
    vmax = float(np.max(vals))
    if not np.isfinite(vmin) or not np.isfinite(vmax):
        raise ValueError("Input values are not finite.")
    if vmax <= vmin:
        return vmin

    hist, bin_edges = np.histogram(vals, bins=int(max(16, n_bins)), range=(vmin, vmax))
    hist = hist.astype(np.float64, copy=False)
    if hist.sum() <= 0:
        return 0.5 * (vmin + vmax)

    bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    weight1 = np.cumsum(hist)
    weight2 = hist.sum() - weight1

    mean1_num = np.cumsum(hist * bin_centres)
    mean_total = mean1_num[-1]

    valid = (weight1 > 0.0) & (weight2 > 0.0)
    if not np.any(valid):
        return 0.5 * (vmin + vmax)

    mean1 = np.zeros_like(bin_centres, dtype=np.float64)
    mean2 = np.zeros_like(bin_centres, dtype=np.float64)
    mean1[valid] = mean1_num[valid] / weight1[valid]
    mean2[valid] = (mean_total - mean1_num[valid]) / weight2[valid]

    sigma_b2 = np.full_like(bin_centres, -np.inf, dtype=np.float64)
    sigma_b2[valid] = weight1[valid] * weight2[valid] * (mean1[valid] - mean2[valid]) ** 2
    idx = int(np.argmax(sigma_b2))
    return float(bin_centres[idx])


def split_mask_utils(args):
    if args.o_below is None and args.o_above is None:
        raise ValueError("Provide at least one output: --o_below and/or --o_above.")

    map_vol, map_shape_header, map_spacing = _read_mrc_as_3d(args.map)
    mask_vol, mask_shape_header, _mask_spacing = _read_mrc_as_3d(args.mask)

    if map_shape_header != mask_shape_header:
        raise ValueError(
            f"Map and mask shape mismatch: map={map_shape_header} vs mask={mask_shape_header}"
        )

    inside = mask_vol > 0.0
    n_inside = int(np.count_nonzero(inside))
    if n_inside == 0:
        raise ValueError("The input mask contains no voxels > 0.")

    work_map = map_vol.astype(np.float32, copy=False)
    blur_sigma = float(args.blur)
    if blur_sigma > 0.0:
        try:
            from scipy.ndimage import gaussian_filter
        except Exception as exc:
            raise RuntimeError("scipy.ndimage.gaussian_filter is required for --blur.") from exc
        work_map = gaussian_filter(work_map, sigma=blur_sigma, mode="nearest").astype(np.float32, copy=False)

    vals = work_map[inside]
    thr = _otsu_threshold_from_values(vals)

    inside_f = inside.astype(np.float32, copy=False)
    low_mask = ((work_map <= thr) & inside).astype(np.float32, copy=False) * inside_f
    above_mask = ((work_map > thr) & inside).astype(np.float32, copy=False) * inside_f

    nx, ny, nz = map_shape_header
    if args.o_below is not None:
        janas_core.WriteMRC(
            low_mask.ravel(order="C").tolist(),
            args.o_below,
            nx,
            ny,
            nz,
            map_spacing,
        )
    if args.o_above is not None:
        janas_core.WriteMRC(
            above_mask.ravel(order="C").tolist(),
            args.o_above,
            nx,
            ny,
            nz,
            map_spacing,
        )

    n_low = int(np.count_nonzero(low_mask))
    n_above = int(np.count_nonzero(above_mask))
    print(
        "[split_mask] "
        f"threshold={thr:.6f} "
        f"inside_voxels={n_inside} "
        f"low_voxels={n_low} "
        f"above_voxels={n_above} "
        f"blur_sigma_vox={blur_sigma:.3f}"
    )


#################################
## janas_maskedCrop
janas_maskedCrop = command.add_parser(
    "maskedCrop",
    description="compute automatic crop of an image based on mask",
    help="compute automatic crop of an image based on mask",
)
janas_maskedCrop.add_argument(
    "--i", required=True, type=str, help="file with the input mrc map"
)
janas_maskedCrop.add_argument(
    "--mask",
    required=True,
    type=str,
    default="",
    help="file with the input mask mrc map",
)
janas_maskedCrop.add_argument(
    "--padding",
    required=False,
    type=int,
    default=2,
    help="padding value for the mask, default=2",
)
janas_maskedCrop.add_argument("--o", required=True, type=str, help="outputFilename")


def maskedCrop(args):
    sizeMap = janas_core.sizeMRC(args.i)
    sizeMask = janas_core.sizeMRC(args.mask)
    inputMap = np.array(janas_core.ReadMRC(args.i))
    inputMap = np.reshape(inputMap, sizeMap)
    inputMask = np.array(janas_core.ReadMRC(args.mask))
    inputMask = np.reshape(inputMask, sizeMap)

    nonzero_indices = np.where(inputMask > 0.1)
    min_z = np.min(nonzero_indices[0])
    max_z = np.max(nonzero_indices[0]) + 1
    min_y = np.min(nonzero_indices[1])
    max_y = np.max(nonzero_indices[1]) + 1
    min_x = np.min(nonzero_indices[2])
    max_x = np.max(nonzero_indices[2]) + 1

    max_size = max(max_x - min_x, max_y - min_y, max_z - min_z)

    # Adjust bounding box
    half_size = max_size // 2
    center_z = (max_z + min_z) // 2
    center_y = (max_y + min_y) // 2
    center_x = (max_x + min_x) // 2
    min_z = center_z - half_size
    max_z = center_z + half_size + max_size % 2
    min_y = center_y - half_size
    max_y = center_y + half_size + max_size % 2
    min_x = center_x - half_size
    max_x = center_x + half_size + max_size % 2

    # Introduce padding
    padding = args.padding
    min_z -= padding
    max_z += padding
    min_y -= padding
    max_y += padding
    min_x -= padding
    max_x += padding

    # Ensure the bounding box doesn't go outside the original image
    min_z = max(min_z, 0)
    max_z = min(max_z, inputMap.shape[0])
    min_y = max(min_y, 0)
    max_y = min(max_y, inputMap.shape[1])
    min_x = max(min_x, 0)
    max_x = min(max_x, inputMap.shape[2])
    croppedMap = inputMap[min_z:max_z, min_y:max_y, min_x:max_x]

    print("Size:", inputMap.size)
    print("Shape:", inputMap.shape)
    print("Data Type:", inputMap.dtype)
    janas_core.WriteMRC(
        croppedMap.flatten().tolist(),
        args.o,
        croppedMap.shape[2],
        croppedMap.shape[1],
        croppedMap.shape[0],
        1,
    )

    # print("sizeMap=",sizeMap)
    # janas_core.WriteEmptyMRC(outputStackBasename+'.mrcs',sizeI[0],sizeI[1],len(imageNames[imageNameTag]))





def get_rotation_matrix(phi, theta, psi):
    # Convert angles from degrees to radians
    alpha = np.radians(phi)
    beta = np.radians(theta)
    gamma = np.radians(psi)

    ca = np.cos(alpha)
    cb = np.cos(beta)
    cg = np.cos(gamma)
    sa = np.sin(alpha)
    sb = np.sin(beta)
    sg = np.sin(gamma)
    cc = cb * ca
    cs = cb * sa
    sc = sb * ca
    ss = sb * sa
    RMatrix = np.array(
        [
            [cg * cc - sg * sa, cg * cs + sg * ca, -cg * sb],
            [-sg * cc - cg * sa, -sg * cs + cg * ca, sg * sb],
            [sc, ss, cb],
        ]
    )
    return RMatrix


def kaiser_bessel_window(x, alpha):
    """
    x: normalized distance to the center (should be between 0 and 0.5)
    alpha: shape parameter. Higher values yield a wider window.
    """
    # Ensure x is within the valid range
    if x < 0 or x > 0.5:
        return 0
    else:
        # Compute the window value
        value = 1 - (2 * x) ** 2
        return np.i0(alpha * np.sqrt(value)) / np.i0(alpha)


def create_circle_mask(image_shape, radius):
    h, w = image_shape
    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - w / 2) ** 2 + (Y - h / 2) ** 2)

    mask = dist_from_center <= radius
    return mask


def insert_image(image2D, fourier2D_mask, outmap, outmapCounter, phi, theta, psi):
    nx, ny, nz = outmap.shape

    # Calculate the center of the volume
    center = np.array([nx // 2, ny // 2, nz // 2])

    # Calculate the rotation matrix
    rotation_matrix = get_inverse_rotation_matrix(phi, theta, psi)
    # Compute the centered Fourier transform of the 2D image
    fourier_image2D = ifftshift(fft2(image2D)) * fourier2D_mask

    # Rotate the 2D image
    # rotated_image = rotate_image(image2D, rotation_matrix, center)
    rotated_image = rotate_image(fourier_image2D, rotation_matrix, center)

    # Insert the rotated image into the 3D volume at its center
    for i in range(rotated_image.shape[0]):
        for j in range(rotated_image.shape[1]):
            x, y, z = rotated_image[j, i]
            if 0 <= x < nx and 0 <= y < ny and 0 <= z < nz:
                outmap[int(z), int(y), int(x)] += fourier_image2D[j, i]
                outmapCounter[int(z), int(y), int(x)] += 1


def get_inverse_rotation_matrix(phi, theta, psi):
    # Convert angles from degrees to radians
    alpha = np.radians(phi)
    beta = np.radians(theta)
    gamma = np.radians(psi)

    ca = np.cos(alpha)
    cb = np.cos(beta)
    cg = np.cos(gamma)
    sa = np.sin(alpha)
    sb = np.sin(beta)
    sg = np.sin(gamma)
    cc = cb * ca
    cs = cb * sa
    sc = sb * ca
    ss = sb * sa
    RMatrix = np.array(
        [
            [cg * cc - sg * sa, cg * cs + sg * ca, -cg * sb],
            [-sg * cc - cg * sa, -sg * cs + cg * ca, sg * sb],
            [sc, ss, cb],
        ]
    )

    # Return the transpose of the rotation matrix
    return RMatrix.T


def rotate_image(image2D, rotation_matrix, center):
    height, width = image2D.shape
    rotated_image = np.zeros((height, width, 3), dtype=int)
    offset = np.array([height // 2, width // 2, 0])
    for i in range(height):
        for j in range(width):
            point = np.array([i, j, 0]) - offset
            rotated_point = np.dot(rotation_matrix, point) + center
            rotated_image[i, j] = rotated_point
    return rotated_image


#################################
## equalize_images
janas_equalize_images = command.add_parser(
    "equalize_images", description="equalize images", help="equalize images"
)
janas_equalize_images.add_argument(
    "--i", required=True, nargs="+", type=str, help="files with the input MRC maps"
)
janas_equalize_images.add_argument(
    "--o_suffix",
    default="_amplEqualized",
    type=str,
    help="suffix to add before the file extension",
)
janas_equalize_images.add_argument(
    "--dir", default="./", type=str, help="output directory"
)


def equalize_images(args):
    fourier_transforms = []
    sizes = []
    sum_amplitudes = None
    count = 0

    # First pass: Compute Fourier transforms and sum their amplitudes
    for input_file in args.i:
        if not os.path.exists(input_file):
            print(f"[ERROR] Input file not found: {input_file}")
            sys.exit(1)
        sizeMap = janas_core.sizeMRC(input_file)
        inputMap = np.array(janas_core.ReadMRC(input_file))
        inputMap = np.reshape(inputMap, sizeMap)

        spacingMRC = round(janas_core.spacingMRC(input_file), 4)
        # print ("spacing=",spacingMRC)

        ft_inputMap = fftn(inputMap)
        amplitude = np.abs(ft_inputMap)

        if sum_amplitudes is None:
            sum_amplitudes = np.zeros_like(amplitude, dtype=np.complex128)

        sum_amplitudes += amplitude
        fourier_transforms.append(ft_inputMap)
        sizes.append(sizeMap)
        count += 1

    # Calculate average amplitude
    average_amplitude = sum_amplitudes / count
    if not os.path.exists(args.dir):
        os.makedirs(args.dir)

    # Second pass: Replace amplitude and inverse Fourier transform
    for i, ft_map in enumerate(fourier_transforms):
        modified_ft = average_amplitude * (ft_map / np.abs(ft_map))
        modified_map = np.real(ifftn(modified_ft))
        nx, ny, nz = sizes[i]
        input_file = args.i[i]
        base_name, ext = os.path.splitext(os.path.basename(input_file))
        output_filename = f"{base_name}{args.o_suffix}{ext}"
        output_filename = os.path.join(args.dir, output_filename)
        janas_core.WriteMRC(
            modified_map.flatten().tolist(), output_filename, nx, ny, nz, spacingMRC
        )


#################################
## scores_to_csv
def retrieveScoringTag(columnList):
    def check_sore_tag_correct_format(s):
        # Splitting the string by underscores
        parts = s.split("_")
        # Checking the structure of the string
        if (
            len(parts) == 8
            and parts[1] in ("janas", "emprove")
            and parts[2] == "SCI"
            and parts[5] == "scored"
            and parts[6] == "selection"
        ):
            return True
        return False

    matching_string = None
    for name in columnList:
        if check_sore_tag_correct_format(name):
            matching_string = name
    return matching_string


janas_scores_to_csv = command.add_parser(
    "scores_to_csv", description="scores_to_csv", help="scores_to_csv"
)
janas_scores_to_csv.add_argument(
    "--i", required=True, nargs="+", type=str, help="files with the input scores"
)
janas_scores_to_csv.add_argument(
    "--csv", default="scoreFile.csv", type=str, help="score_file_to_csv"
)
janas_scores_to_csv.add_argument(
    "--o", default="scoreFile.star", type=str, help="scoreFile.star"
)


def scores_to_csv(args):
    result_df = pd.DataFrame()
    file_class_mapping = {
        file_name: index + 1 for index, file_name in enumerate(args.i)
    }
    reference_star = None
    for input_file in args.i:
        if reference_star == None:
            reference_star = input_file
        columns = starHandler.header_columns(input_file)
        tagScore = retrieveScoringTag(columns)
        referenceColumns = starHandler.readColumns(input_file, [tagScore])
        result_df[input_file] = referenceColumns[tagScore]
    result_df["Max_Score_File"] = result_df.idxmax(axis=1)
    result_df["_rlnClassNumber"] = result_df["Max_Score_File"].map(file_class_mapping)
    result_df.to_csv(args.csv)
    # now update the class on the star file
    version = starHandler.infoStarFile(input_file)[2]
    main_section_name = "particles"
    if version == "relion_v30":
        main_section_name = ""
    starHandler.removeColumnsTagsStartingWith(reference_star, args.o, "_janas_")
    starHandler.removeColumnsTagsStartingWith(args.o, args.o, "_emprove_")
    starHandler.replace_star_columns_from_sections(
        args.o, args.o, main_section_name, "_rlnClassNumber", result_df
    )


janas_randomize_halves = command.add_parser(
    "randomize_halves", description="randomize_halves", help="randomize_halves"
)
janas_randomize_halves.add_argument(
    "--i", required=True, type=str, help="input star file"
)
janas_randomize_halves.add_argument(
    "--o", required=True, type=str, help="output starfile"
)


def randomize_halves_utils(args):
    starHandler.randomize_halves(args.i, args.o)


#################################
janas_extract_particles_from_label_value = command.add_parser(
    "extract_particles_from_label_value",
    description="extract_particles_from_label_value",
    help="extract_particles_from_label_value",
)
janas_extract_particles_from_label_value.add_argument(
    "--i", required=True, type=str, help="input star file"
)
janas_extract_particles_from_label_value.add_argument(
    "--label", required=True, type=str, help="label to extract values from"
)
janas_extract_particles_from_label_value.add_argument(
    "--value", required=True, type=str, help="value to extract"
)
janas_extract_particles_from_label_value.add_argument(
    "--o", required=True, type=str, help="output star file"
)


def extract_particles_from_label_value(args):
    version = starHandler.infoStarFile(args.i)[2]
    main_section_name = "particles"
    if version == "relion_v30":
        main_section_name = ""
    starHandler.extract_particles_from_label_from_sections(
        args.i, args.o, main_section_name, args.label, args.value
    )

######################
## convert cs to star
janas_csparc2star = command.add_parser(
    "csparc2star",
    description="convert CryoSPARC .cs or .csv to Relion .star",
    help="convert CryoSPARC .cs or .csv to Relion .star"
)
janas_csparc2star.add_argument(
    "input", nargs=1, help="input .cs (CryoSPARC2+) or .csv (CryoSPARC0.6.x) file"
)
janas_csparc2star.add_argument(
    "output", help="output Relion .star filename"
)
janas_csparc2star.add_argument(
    "--loglevel", "-l", default="WARNING",
    help="logging level (DEBUG, INFO, WARNING)"
)
janas_csparc2star.add_argument(
    "--clean_path", action="store_true",
    help=(
        "Drop the directory from each blob path and keep only the filename "
        "in the generated _rlnImageName. Overridden by --fix_path if both are set."
    ),
)
janas_csparc2star.add_argument(
    "--clean_prefix", action="store_true",
    help=(
        "Strip the leading numeric CryoSPARC prefix (matching '^\\d+_') from "
        "the filename in each generated _rlnImageName. Applies only to the "
        "filename, never to the directory."
    ),
)
janas_csparc2star.add_argument(
    "--clean_suffix", action="store_true",
    help=(
        "Strip a terminal '_particles' from the filename stem (before the "
        "extension) in each generated _rlnImageName."
    ),
)
janas_csparc2star.add_argument(
    "--fix_path", default=None,
    help=(
        "Replace the original directory in each blob path with this path. "
        "A trailing slash is tolerated. Takes precedence over --clean_path."
    ),
)
janas_csparc2star.add_argument(
    "--missing_pose_to_zero", action="store_true",
    help=(
        "Allow conversion when the input .cs has no 3D alignment metadata "
        "(typical for extraction/picking/passthrough/coordinate-only jobs). "
        "Missing alignments3D/pose -> _rlnAngleRot/Tilt/Psi = 0; missing "
        "alignments3D/shift -> _rlnOriginXAngst/_rlnOriginYAngst = 0. "
        "Without this flag, missing pose/shift raises a clear error. "
        "Note: the resulting STAR is unaligned, and the origin columns are "
        "refinement shifts (not extraction coordinates)."
    ),
)
def list_cs_fields(cs_file: str):
    arr = np.load(cs_file, max_header_size=100000)
    print(">>> available dtype.names in", cs_file)
    for name in arr.dtype.names:
        print("   ", name)

def csparc2star_utils(args):
    #list_cs_fields(args.input[0])
    utils.csparc2star(
        infile=args.input[0],
        outfile=args.output,
        clean_path=args.clean_path,
        clean_prefix=args.clean_prefix,
        clean_suffix=args.clean_suffix,
        fix_path=args.fix_path,
        missing_pose_to_zero=args.missing_pose_to_zero,
    )


######################
## convert cs to star
janas_update_from_csparc = command.add_parser(
    "update_from_csparc",
    description="convert CryoSPARC .cs or .csv to Relion .star",
    help="convert CryoSPARC .cs or .csv to Relion .star"
)
janas_update_from_csparc.add_argument(
    "input_cs", nargs=1, help="input .cs (CryoSPARC2+) or .csv (CryoSPARC0.6.x) file"
)
janas_update_from_csparc.add_argument(
    "input_star", nargs=1, help="input star file (Relion >3.1) file"
)
janas_update_from_csparc.add_argument(
    "output", help="output  .star (Relion >3.1) filename"
)
janas_update_from_csparc.add_argument(
    "--loglevel", "-l", default="WARNING",
    help="logging level (DEBUG, INFO, WARNING)"
)


def janas_update_from_csparc_utils(args):
    utils.update_star_from_csparc(
        csfile=args.input_cs[0],
        starfile_in=args.input_star[0],
        starfile_out=args.output,
        loglevel=args.loglevel
    )


######################
## csparc2star-stack
janas_csparc2star_stack = command.add_parser(
    "csparc2star-stack",
    description="convert CryoSPARC .cs to a RELION-style STAR and build a new consolidated particle stack (.mrcs)",
    help="convert CryoSPARC .cs to STAR and assemble a new .mrcs stack"
)
janas_csparc2star_stack.add_argument(
    "input", nargs=1, help="input .cs (CryoSPARC2+) particles file"
)
janas_csparc2star_stack.add_argument(
    "output_root", help="output rootname (writes output_root.star and output_root.mrcs)"
)
janas_csparc2star_stack.add_argument(
    "--root", default=None,
    help="CryoSPARC project root (directory that contains job folders J1, J2, ...). "
         "If omitted, inferred as parent of the job directory that contains the input .cs"
)
janas_csparc2star_stack.add_argument(
    "--loglevel", "-l", default="WARNING",
    help="logging level (DEBUG, INFO, WARNING)"
)

def _infer_cs_project_root(cs_path: str) -> str:
    """
    Infer CryoSPARC project root from a path like .../J44/J44_006_particles.cs
    by returning the parent directory of the first 'J<digits>' folder up the tree.
    """
    d = os.path.dirname(os.path.abspath(cs_path))
    # Walk up until we exit the job directory (Jxx)
    while True:
        base = os.path.basename(d)
        if not re.fullmatch(r"J\d+", base):
            # 'd' is now at project root OR above it; if the child we just left was Jxx,
            # then this 'd' is the project root.
            return d
        d = os.path.dirname(d)

def _resolve_stack_path(image_name: str, project_root: str) -> Tuple[int, str]:
    """
    Parse 'NNNNNN@path/to/stack.(mrc|mrcs)' and return (index_0based, absolute_stack_path).
    """
    at = image_name.find("@")
    if at < 0:
        raise ValueError(f"Invalid _rlnImageName entry (no '@'): {image_name}")
    idx_str = image_name[:at].strip()
    stack_rel = image_name[at+1:].strip()

    try:
        image_no = int(idx_str)
    except Exception:
        # Some STARs use zero-padded strings; int() handles them.
        raise ValueError(f"Invalid image index in _rlnImageName: {image_name}")

    # Absolute vs relative (CryoSPARC often stores relative like 'J2/imported/.../stack.mrc')
    if os.path.isabs(stack_rel):
        stack_path = stack_rel
    else:
        stack_path = os.path.normpath(os.path.join(project_root, stack_rel))

    return image_no - 1, stack_path  # zero-based slice index



def csparc2star_stack_utils(args):
    cs_file = args.input[0]
    out_root = args.output_root
    loglevel = args.loglevel

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_star = path.join(tmpdir, "tmp_particles.star")
        utils.csparc2star(infile=cs_file, outfile=tmp_star, loglevel=loglevel)

        # Resolve project root (as before)
        project_root = path.abspath(args.root) if args.root else utils.infer_cs_project_root_from_path(cs_file)

        # Reuse the shared builder; for csparc we keep a specific provenance tag
        utils.create_stack_from_star(
            star_in=tmp_star,
            out_root=out_root,
            project_root=project_root,
            provenance_tag="_janas_csparc_rlnImageName",
            override_section_name="particles"  # csparc2star writes Relion 3.1 layout
        )



########################
######################
## create_stack (from STAR)
janas_create_stack = command.add_parser(
    "create_stack",
    description="Create a consolidated .mrcs stack from a STAR file and repoint _rlnImageName",
    help="create a consolidated stack from STAR"
)
janas_create_stack.add_argument(
    "input_star", help="input .star file with _rlnImageName entries"
)
janas_create_stack.add_argument(
    "output_root", help="output rootname (writes output_root.star and output_root.mrcs)"
)
janas_create_stack.add_argument(
    "--root", default=None,
    help="Base directory to resolve relative stack paths found in _rlnImageName. "
         "For CryoSPARC-style relative paths you can point this to the project root."
)
janas_create_stack.add_argument(
    "--provenance-tag", default="_janas_source_rlnImageName",
    help="Optional STAR tag to store the original _rlnImageName (default: _janas_source_rlnImageName). "
         "Use '' to disable."
)
janas_create_stack.add_argument(
    "--path_mode", "--path-mode",
    dest="path_mode",
    default="auto",
    choices=["auto", "root", "as_is", "star_dir"],
    help=(
        "How to resolve relative stack paths read from _rlnImageName. "
        "'auto' (default): try --root, then the path as written (CWD-relative), "
        "then the directory of the input STAR file; the first existing match wins. "
        "'root': resolve against --root only (original CryoSPARC behaviour). "
        "'as_is': use the path exactly as written; relative paths are CWD-relative. "
        "'star_dir': resolve relative paths against the directory containing the input STAR."
    ),
)

def create_stack_utils(args):
    prov = args.provenance_tag if (args.provenance_tag and len(args.provenance_tag.strip()) > 0) else None
    utils.create_stack_from_star(
        star_in=args.input_star,
        out_root=args.output_root,
        project_root=(os.path.abspath(args.root) if args.root else None),
        provenance_tag=prov,
        path_mode=args.path_mode,
    )


######################
## backmap_stars
janas_backmap_stars = command.add_parser(
    "backmap_stars",
    description=(
        "Restore original source _rlnImageName values in a downstream STAR file by "
        "joining against the STAR produced by create_stack (which carries "
        "_janas_source_rlnImageName). Inverse companion of create_stack."
    ),
    help="restore original _rlnImageName via a stack-generation STAR mapping"
)
janas_backmap_stars.add_argument(
    "--processed", required=True, dest="processed_star",
    help="downstream STAR file to fix (its _rlnImageName points to the consolidated stack)."
)
janas_backmap_stars.add_argument(
    "--mapping", required=True, dest="mapping_star",
    help="stack-generation STAR file (contains _rlnImageName and _janas_source_rlnImageName)."
)
janas_backmap_stars.add_argument(
    "--output", required=True, dest="output_star",
    help="output STAR file with restored source _rlnImageName."
)
janas_backmap_stars.add_argument(
    "--image-tag", default="_rlnImageName",
    help="column holding the consolidated-stack reference (default: _rlnImageName)."
)
janas_backmap_stars.add_argument(
    "--source-tag", default="_janas_source_rlnImageName",
    help="column in --mapping that holds the original source reference "
         "(default: _janas_source_rlnImageName)."
)
janas_backmap_stars.add_argument(
    "--stack-reference-tag", default="_janas_stack_rlnImageName",
    help="audit column to write into the output with the previous (stack-based) "
         "_rlnImageName values. Pass '' to disable. "
         "(default: _janas_stack_rlnImageName)"
)
janas_backmap_stars.add_argument(
    "--section-name", default=None,
    help="override the particle data block name. Inferred automatically if omitted."
)
janas_backmap_stars.add_argument(
    "--no-strict", action="store_true",
    help="leave unmapped _rlnImageName values unchanged instead of failing."
)


def backmap_stars_utils(args):
    ref_tag = args.stack_reference_tag
    if ref_tag is not None and ref_tag.strip() == "":
        ref_tag = None
    report = utils.backmap_stars(
        processed_star=args.processed_star,
        mapping_star=args.mapping_star,
        output_star=args.output_star,
        image_tag=args.image_tag,
        source_tag=args.source_tag,
        stack_reference_tag=ref_tag,
        section_name=args.section_name,
        strict=not args.no_strict,
    )
    print(
        f"backmap_stars: mapped {report['n_mapped']}/{report['n_processed']} rows "
        f"(missing={report['n_missing']}); wrote {report['output_star']}"
    )
    if report["n_missing"] > 0 and report["missing_examples"]:
        print("  example missing keys:")
        for k in report["missing_examples"]:
            print(f"    - {k}")



######################
## get angpix
janas_angpix = command.add_parser(
    "angpix",
    description="print pixel spacing (Å/pixel) for an MRC map",
    help="print pixel spacing for an MRC file"
)
janas_angpix.add_argument(
    "file", type=str, help="input MRC file to query"
)

def print_angpix(args):
    spacing = utils.get_MRC_map_pixel_spacing(args.file)
    try:
        pixel_size = spacing[0]
    except (TypeError, IndexError):
        pixel_size = spacing
    formatted = f"{pixel_size:.3f}"
    print(formatted)




######################
## janas_mask
janas_mask = command.add_parser(
    "mask",
    description="JANAS mask: mask manipulation utilities, for more details use janas_utils mask --help",
    help="JANAS mask: mask manipulation utilities, for more details use janas_utils mask --help",
)
mask_subparsers = janas_mask.add_subparsers(dest="mask_command")

# mask merge
mask_merge = mask_subparsers.add_parser(
    "merge",
    help=(
        "Merge (union) multiple binary MRC masks into one. "
        "Usage: janas_utils mask merge mask1.mrc mask2.mrc ... maskN.mrc output.mrc"
    ),
)
mask_merge.add_argument(
    "inputs", nargs="+", type=str,
    help="input MRC mask files followed by output MRC filename (last argument is output)"
)


# mask countPixels
mask_countPixels = mask_subparsers.add_parser(
    "countPixels",
    help=(
        "Count voxels inside a mask MRC volume.\n\n"
        "Default mode (no --value): counts all voxels strictly greater than 0.01.\n"
        "This is intended for non-binary floating-valued masks.\n\n"
        "Value mode (--value V): counts voxels within a tolerance window around V,\n"
        "suitable for floating-point mask data where exact equality is unreliable.\n"
        "For example, --value 1 counts voxels in the interval [0.99, 1.1].\n\n"
        "Usage:\n"
        "  janas_utils mask countPixels mask.mrc\n"
        "  janas_utils mask countPixels mask.mrc --value 1\n"
    ),
)
mask_countPixels.add_argument(
    "mask_file", type=str,
    help="Input MRC mask file."
)
mask_countPixels.add_argument(
    "--value", type=float, default=None,
    help=(
        "If provided, count voxels near this value using a tolerance window "
        "rather than counting all positive voxels. "
        "E.g. --value 1 counts voxels in [0.99, 1.1]."
    ),
)


def countPixels_core(mask_path: str, value: float = None) -> int:
    """Count voxels in an MRC mask.

    Parameters
    ----------
    mask_path : str
        Path to the MRC file.
    value : float or None
        If *None* (default), count voxels with value > 0.01.
        If a float is given, count voxels within a tolerance window around
        that value.  The window is [value - 0.01*max(|value|,1),
        value + 0.1*max(|value|,1)], which for ``value=1`` yields
        exactly [0.99, 1.1] as required by specification.

    Returns
    -------
    int
        Number of matching voxels.
    """
    if not os.path.exists(mask_path):
        raise FileNotFoundError(f"File not found: {mask_path}")

    vol, _shape, _spacing = _read_mrc_as_3d(mask_path)

    if value is None:
        # Default: count all voxels > 0.01
        return int(np.count_nonzero(vol > 0.01))
    else:
        # Tolerance window: [value - 0.01*max(|value|,1), value + 0.1*max(|value|,1)]
        # For value=1 this gives exactly [0.99, 1.1] as required by specification.
        scale = max(abs(value), 1.0)
        lo = value - 0.01 * scale
        hi = value + 0.1 * scale
        return int(np.count_nonzero((vol >= lo) & (vol <= hi)))


def mask_countPixels_utils(args):
    """CLI entry point for ``janas_utils mask countPixels``."""
    mask_path = args.mask_file
    if not mask_path.lower().endswith('.mrc'):
        raise SystemExit(f"Input must be an MRC file; got: {mask_path}")
    count = countPixels_core(mask_path, value=args.value)
    print(count)


def mask_merge_utils(args):
    """Compute the voxel-wise union of one or more binary masks.

    Any voxel that is > 0 in *any* input mask will be 1.0 in the output.
    All inputs must share the same grid dimensions.

    Usage: janas_utils mask merge m1.mrc m2.mrc ... mN.mrc output.mrc
    """
    items = args.inputs
    if len(items) < 2:
        raise SystemExit(
            "Usage: janas_utils mask merge mask1.mrc [mask2.mrc ...] output.mrc\n"
            "At least one input mask and one output path are required."
        )

    input_files = items[:-1]
    output_file = items[-1]

    # Read the first mask to establish reference shape and spacing
    merged, ref_shape, ref_spacing = _read_mrc_as_3d(input_files[0])
    merged = (merged > 0.0).astype(np.float32)

    # Union with remaining masks
    for fname in input_files[1:]:
        vol, shape, _sp = _read_mrc_as_3d(fname)
        if shape != ref_shape:
            raise ValueError(
                f"Mask shape mismatch: {input_files[0]} has {ref_shape} "
                f"but {fname} has {shape}"
            )
        merged = np.maximum(merged, (vol > 0.0).astype(np.float32))

    nx, ny, nz = ref_shape
    janas_core.WriteMRC(
        merged.ravel(order="C").tolist(),
        output_file,
        nx, ny, nz,
        ref_spacing,
    )


######################
## janas_clip
janas_clip = command.add_parser(
    "clip",
    description="JANAS clip: JANAS command line image processing, for more details use janas_utils clip --help",
    help="JANAS clip: JANAS command line image processing, for more details use janas_utils clip --help",
)
clip_subparsers = janas_clip.add_subparsers(dest="clip_command")

# clip add
clip_add = clip_subparsers.add_parser(
    "add", help="add multiple MRC maps into a single volume"
)
clip_add.add_argument(
    "inputs", nargs="*", type=str, help="input MRC files to sum"
)
clip_add.add_argument(
    "output", type=str, help="output MRC filename"
)

# clip amplReplace
clip_amplReplace = clip_subparsers.add_parser(
    "amplReplace",
    help=("Replace Fourier amplitudes.\n"
          "Usage: janas_utils clip amplReplace PHASE.mrc AMPL.mrc OUT.mrc")
)
clip_amplReplace.add_argument(
    "items", nargs="+",
    help="PHASE.mrc AMPL.mrc OUT.mrc"
)

def clip_amplReplace_utils(args):
    items = args.items
    if len(items) != 3:
        raise SystemExit("Usage: janas_utils clip amplReplace PHASE.mrc AMPL.mrc OUT.mrc")

    phase_mrc, ampl_mrc, out_mrc = items
    if not (phase_mrc.lower().endswith(".mrc") and ampl_mrc.lower().endswith(".mrc") and out_mrc.lower().endswith(".mrc")):
        raise SystemExit("amplReplace expects: PHASE.mrc AMPL.mrc OUT.mrc")

    s1 = janas_core.sizeMRC(phase_mrc)     # shape used everywhere else in this file
    s2 = janas_core.sizeMRC(ampl_mrc)
    if s1 != s2:
        raise ValueError(f"Input size mismatch: {phase_mrc} {s1} vs {ampl_mrc} {s2}")

    v1 = np.array(janas_core.ReadMRC(phase_mrc), dtype=np.float32).reshape(s1)
    v2 = np.array(janas_core.ReadMRC(ampl_mrc), dtype=np.float32).reshape(s2)

    ft1 = fftn(v1)
    ft2 = fftn(v2)

    amp2 = np.abs(ft2)
    amp1 = np.abs(ft1)
    phase1 = ft1 / np.maximum(amp1, 1e-12)  # avoid division by zero

    out = np.real(ifftn(amp2 * phase1)).astype(np.float32, copy=False)

    spacing = round(janas_core.spacingMRC(phase_mrc), 4)
    nz, ny, nx = out.shape
    janas_core.WriteMRC(out.flatten().tolist(), out_mrc, nx, ny, nz, spacing)

# clip roundMask
clip_roundMask = clip_subparsers.add_parser(
    "roundMask",
    help=("Create a centred circular/spherical mask that fits in the map.\n"
          "Usage:\n"
          "  janas_utils clip roundMask IN.mrc OUT_MASK.mrc\n"
          "  janas_utils clip roundMask IN.mrc OUT_MASK.mrc RADIUS\n"
          "  janas_utils clip roundMask IN.mrc OUT_MASK.mrc RADIUS [X,Y,Z]")
)
clip_roundMask.add_argument(
    "items", nargs="+",
    help="IN.mrc OUT_MASK.mrc [RADIUS] [[X,Y,Z]]"
)


# clip average
clip_average = clip_subparsers.add_parser(
    "average", help="average multiple MRC maps into a single volume"
)
clip_average.add_argument(
    "inputs", nargs="*", type=str, help="input MRC files to average"
)
clip_average.add_argument(
    "output", type=str, help="output MRC filename"
)

# clip blur
clip_blur = clip_subparsers.add_parser(
    "blur",
    help=("Gaussian blur a volume with σ specified in Å.\n"
          "Usage: janas_utils clip blur IN.mrc OUT.mrc SIGMA_A")
)
clip_blur.add_argument(
    "items", nargs="+",
    help="Arguments: IN.mrc OUT.mrc SIGMA_A  (σ in Angstrom)"
)



# clip rot2D
clip_rot2D = clip_subparsers.add_parser(
    "rot2D",
    help=("Rotate a 2D MRC image by an angle in degrees.\n"
          "Usage: janas_utils clip rot2D IN.mrc ANGLE_DEG OUT.mrc")
)
clip_rot2D.add_argument(
    "items", nargs="+",
    help="Arguments: IN.mrc ANGLE_DEG OUT.mrc"
)

# clip Bfactor weighting
clip_bfac = clip_subparsers.add_parser(
    "bfac",
    help=("B-factor weighting.\n"
          "automatic:  clip bfac HALF1.mrc HALF2.mrc OUT.mrc\n"
          "user-driven: clip bfac HALF1.mrc HALF2.mrc OUT.mrc BVALUE")
)
clip_bfac.add_argument(
    "items", nargs="+",
    help="automatic: half1 half2 out; user-driven: half1 half2 out BVALUE"
)



def clip_add_maps(args):
    # Read all volumes
    volumes = []
    for fname in args.inputs:
        size = janas_core.sizeMRC(fname)
        data = np.array(janas_core.ReadMRC(fname)).reshape(size)
        volumes.append(data)
    if not volumes:
        raise ValueError("No input volumes provided for 'clip add'.")
    # Ensure consistent dimensions
    ref_shape = volumes[0].shape
    if any(v.shape != ref_shape for v in volumes):
        raise ValueError("All input volumes must have identical dimensions.")
    # Sum them
    result = sum(volumes)
    # Retrieve pixel spacing from the first volume
    spacing = round(janas_core.spacingMRC(args.inputs[0]), 4)
    nz, ny, nx = ref_shape
    janas_core.WriteMRC(
        result.flatten().tolist(),
        args.output,
        nx, ny, nz,
        spacing
    )

def clip_average_maps(args):
    # Read all volumes
    volumes = []
    for fname in args.inputs:
        size = janas_core.sizeMRC(fname)
        data = np.array(janas_core.ReadMRC(fname)).reshape(size)
        volumes.append(data)
    if not volumes:
        raise ValueError("No input volumes provided for 'clip average'.")
    # Ensure consistent dimensions
    ref_shape = volumes[0].shape
    if any(v.shape != ref_shape for v in volumes):
        raise ValueError("All input volumes must have identical dimensions.")
    # Compute the average
    summed = sum(volumes)
    result = summed / len(volumes)
    # Retrieve pixel spacing
    spacing = round(janas_core.spacingMRC(args.inputs[0]), 4)
    nz, ny, nx = ref_shape
    janas_core.WriteMRC(
        result.flatten().tolist(),
        args.output,
        nx, ny, nz,
        spacing
    )


def clip_blur_utils(args):
    items = args.items
    # Expect exactly: IN.mrc OUT.mrc SIGMA_A
    if len(items) != 3 or not items[0].lower().endswith(".mrc") or not items[1].lower().endswith(".mrc"):
        raise SystemExit(
            "Usage:\n"
            "  janas_utils clip blur IN.mrc OUT.mrc SIGMA_A\n"
            "Where SIGMA_A is the Gaussian σ in Å."
        )

    in_mrc, out_mrc, sigma_str = items
    try:
        sigmaA = float(sigma_str)
    except Exception:
        raise SystemExit(f"SIGMA_A must be numeric (Å); got: {sigma_str}")

    if sigmaA <= 0.0:
        raise SystemExit(f"SIGMA_A must be > 0 (Å); got: {sigmaA}")

    ns = SimpleNamespace(
        i=in_mrc,
        o=out_mrc,
        sigmaA=sigmaA
    )
    janas_mapProcess.gaussian_blur_utils(ns)





def _parse_center_list(s: str):
    """
    Parse a string like "[40,50,60]" or "40,50,60" into a list of floats.
    """
    t = s.strip()
    if t.startswith("[") and t.endswith("]"):
        t = t[1:-1]
    parts = [p.strip() for p in t.split(",") if p.strip() != ""]
    if len(parts) == 0:
        raise ValueError(f"Invalid centre specification: {s}")
    return [float(p) for p in parts]


def clip_roundMask_utils(args):
    """
    Create a circular (2D) or spherical (3D) mask in pixel units.
    Dimension(s) with size < 2 are ignored when deciding 2D vs 3D.
    """
    items = args.items
    if len(items) not in (2, 3, 4):
        raise SystemExit(
            "Usage:\n"
            "  janas_utils clip roundMask IN.mrc OUT_MASK.mrc\n"
            "  janas_utils clip roundMask IN.mrc OUT_MASK.mrc RADIUS\n"
            "  janas_utils clip roundMask IN.mrc OUT_MASK.mrc RADIUS [X,Y,Z]"
        )

    in_mrc = items[0]
    out_mrc = items[1]
    if not in_mrc.lower().endswith(".mrc") or not out_mrc.lower().endswith(".mrc"):
        raise SystemExit("roundMask expects .mrc input and .mrc output.")

    if not os.path.exists(in_mrc):
        raise FileNotFoundError(f"Input file not found: {in_mrc}")

    # Header sizes in your environment: (nx, ny, nz)
    nx, ny, nz = janas_core.sizeMRC(in_mrc)
    spacingMRC = round(janas_core.spacingMRC(in_mrc), 4)

    # Decide active dimensions (ignore dims < 2)
    active_axes = []
    sizes = []
    if nx_bot := nx >= 2:
        active_axes.append("x"); sizes.append(nx)
    if ny_bot := ny >= 2:
        active_axes.append("y"); sizes.append(ny)
    if nz_bot := nz >= 2:
        active_axes.append("z"); sizes.append(nz)

    if len(active_axes) == 0:
        raise SystemExit(f"Invalid map size (all dims < 2): nx,ny,nz = {nx},{ny},{nz}")

    # Default centre: geometric centre in pixel coordinates, in X,Y,Z indexing
    default_center = [ (nx - 1) / 2.0, (ny - 1) / 2.0, (nz - 1) / 2.0 ]

    # Default radius: maximum that fits, considering only active dimensions
    # Distance to nearest face along each active axis
    # For axis with size N, max radius from centre is min(cx, (N-1)-cx).
    cx0, cy0, cz0 = default_center
    rmax_candidates = []
    if nx >= 2:
        rmax_candidates.append(min(cx0, (nx - 1) - cx0))
    if ny >= 2:
        rmax_candidates.append(min(cy0, (ny - 1) - cy0))
    if nz >= 2:
        rmax_candidates.append(min(cz0, (nz - 1) - cz0))
    r_default = float(min(rmax_candidates))

    # Optional radius
    if len(items) >= 3:
        try:
            radius = float(items[2])
        except Exception:
            raise SystemExit(f"RADIUS must be numeric (pixels); got: {items[2]}")
        if radius <= 0:
            raise SystemExit(f"RADIUS must be > 0; got: {radius}")
    else:
        radius = r_default

    # Optional centre
    center = default_center
    if len(items) == 4:
        c_list = _parse_center_list(items[3])
        # Accept [X,Y] for 2D; [X,Y,Z] for 3D. Missing coords fall back to defaults.
        if len(c_list) > 3:
            raise SystemExit("Centre must be [X,Y] or [X,Y,Z].")
        center = default_center.copy()
        for i, v in enumerate(c_list):
            center[i] = float(v)

    cx, cy, cz = center

    # Build mask in MRC-consistent memory layout: (nz, ny, nx) with X fast
    # We write out a mask matching nx,ny,nz header.
    # If nz==1, still make (1,ny,nx).
    zz = np.arange(nz, dtype=np.float32)[:, None, None]  # (nz,1,1)
    yy = np.arange(ny, dtype=np.float32)[None, :, None]  # (1,ny,1)
    xx = np.arange(nx, dtype=np.float32)[None, None, :]  # (1,1,nx)

    # Distance-squared, ignoring any dimension < 2 by zeroing its contribution
    d2 = np.zeros((nz, ny, nx), dtype=np.float32)

    if nx >= 2:
        d2 += (xx - cx) ** 2
    if ny >= 2:
        d2 += (yy - cy) ** 2
    if nz >= 2:
        d2 += (zz - cz) ** 2

    mask = (d2 <= (radius ** 2)).astype(np.float32)

    # Write mask: flatten in C-order keeps X fastest for (nz,ny,nx)
    janas_core.WriteMRC(
        mask.flatten(order="C").tolist(),
        out_mrc,
        nx, ny, nz,
        spacingMRC
    )




def clip_rot2D_utils(args):
    items = args.items
    if len(items) != 3:
        raise SystemExit("Usage: janas_utils clip rot2D IN.mrc ANGLE_DEG OUT.mrc")

    in_mrc, angle_str, out_mrc = items
    try:
        angle_deg = float(angle_str)
    except Exception:
        raise SystemExit(f"ANGLE_DEG must be numeric (degrees); got: {angle_str}")

    if not os.path.exists(in_mrc):
        raise FileNotFoundError(f"Input file not found: {in_mrc}")

    try:
        from scipy.ndimage import rotate as ndi_rotate
    except ImportError:
        raise SystemExit("rot2D requires SciPy: pip install scipy")

    # sizeMRC (in your environment) matches header order: (nx, ny, nz)
    nx, ny, nz = janas_core.sizeMRC(in_mrc)

    flat = np.array(janas_core.ReadMRC(in_mrc), dtype=np.float32)

    # CRITICAL: reshape so that X is the FAST axis (last axis in NumPy):
    # MRC storage is X-fast, then Y, then Z -> NumPy shape must be (nz, ny, nx)
    vol = flat.reshape((nz, ny, nx), order="C")

    spacingMRC = round(janas_core.spacingMRC(in_mrc), 4)

    if nz != 1:
        # If you truly only want 2D behaviour, you can refuse here instead.
        # This rotates each Z-slice independently in XY.
        out = np.empty_like(vol, dtype=np.float32)
        for k in range(nz):
            out[k] = ndi_rotate(
                vol[k],
                angle=angle_deg,
                reshape=False,
                order=1,              # bilinear
                mode="constant",
                cval=0.0,
                prefilter=True
            ).astype(np.float32, copy=False)
    else:
        img2d = vol[0]  # shape (ny, nx) i.e., rows(Y), cols(X) as SciPy expects
        rot2d = ndi_rotate(
            img2d,
            angle=angle_deg,
            reshape=False,
            order=1,
            mode="constant",
            cval=0.0,
            prefilter=True
        ).astype(np.float32, copy=False)

        out = rot2d[None, :, :]  # back to (1, ny, nx)

    # Write back: your WriteMRC expects (nx, ny, nz) in header fields.
    # Flattening (nz, ny, nx) in C-order preserves X-fast storage.
    janas_core.WriteMRC(
        out.flatten(order="C").tolist(),
        out_mrc,
        nx, ny, nz,
        spacingMRC
    )



def clip_normalize_utils(args):
    items = args.items
    # Expect: IN.mrc OUT.mrc [METHOD]
    if len(items) < 2 or len(items) > 3 \
       or not items[0].lower().endswith(".mrc") \
       or not items[1].lower().endswith(".mrc"):
        raise SystemExit(
            "Usage:\n"
            "  janas_utils clip normalize IN.mrc OUT.mrc [METHOD]\n"
            "Where METHOD is one of: zscore, minmax (default: zscore)"
        )

    in_mrc, out_mrc = items[0], items[1]
    method = items[2].lower() if len(items) == 3 else "zscore"
    if method not in ("zscore", "minmax"):
        raise SystemExit(f"Unknown METHOD '{method}'. Use 'zscore' or 'minmax'.")

    ns = SimpleNamespace(
        i=in_mrc,
        o=out_mrc,
        method=method
    )
    janas_mapProcess.normalize_utils(ns)



def clip_bfac_utils(args):
    items = args.items

    # Helper: load two halves, validate, average in RAM
    def _load_and_avg(h1, h2):
        vol1, shape1, ap1 = janas_mapProcess._load_mrc_as_array(h1)
        vol2, shape2, ap2 = janas_mapProcess._load_mrc_as_array(h2)
        if shape1 != shape2:
            raise ValueError(f"Half-map shapes differ: {shape1} vs {shape2}")
        if not np.isfinite(ap1) or not np.isfinite(ap2) or abs(ap1 - ap2) > 1e-6:
            raise ValueError(f"Å/px mismatch: {ap1} vs {ap2}")
        avg = 0.5 * (vol1.astype(np.float64, copy=False) + vol2.astype(np.float64, copy=False))
        spacing = round(janas_core.spacingMRC(h1), 4)
        return avg, shape1, float(ap1), float(spacing)

    # AUTO: half1 half2 out
    if len(items) == 3 and items[2].lower().endswith(".mrc"):
        half1, half2, out = items
        avg, shape, angpix, spacing = _load_and_avg(half1, half2)
        ns = SimpleNamespace(
            i=None,             # not used (in-memory path)
            i_array=avg,        # in-memory averaged map
            shape=shape,
            angpix=angpix,
            spacing=spacing,
            o=out,
            B=None,
            auto=True,
            fit_minres=10.0,
            fit_maxres=0.0,
            hmaps=(half1, half2),   # FSC weighting from halves
            lowpass=0.0
        )
        janas_mapProcess.bfactor_utils(ns)
        return

    # AD-HOC: half1 half2 out BVALUE
    if len(items) == 4 and items[2].lower().endswith(".mrc"):
        half1, half2, out, bval = items
        try:
            B = float(bval)
        except Exception:
            raise ValueError(f"BVALUE must be numeric (Å^2); got: {bval}")
        avg, shape, angpix, spacing = _load_and_avg(half1, half2)
        ns = SimpleNamespace(
            i=None,
            i_array=avg,
            shape=shape,
            angpix=angpix,
            spacing=spacing,
            o=out,
            B=B,
            auto=False,
            fit_minres=10.0,
            fit_maxres=0.0,
            hmaps=(half1, half2),   # retain FSC weighting step (as in RELION)
            lowpass=0.0
        )
        janas_mapProcess.bfactor_utils(ns)
        return

    raise SystemExit(
        "Usage:\n"
        "  AUTO : janas_utils clip bfac HALF1.mrc HALF2.mrc OUT.mrc\n"
        "  ADHOC: janas_utils clip bfac HALF1.mrc HALF2.mrc OUT.mrc BVALUE"
    )


#################################
## map_histogram
janas_map_histogram = command.add_parser(
    "map_histogram",
    description="Plot intensity histogram of a 3D MRC map, optionally within a mask (mask > 0.1).",
    help="plot histogram of a map (optionally masked)"
)
janas_map_histogram.add_argument(
    "map", type=str, help="input MRC map filename"
)
janas_map_histogram.add_argument(
    "mask", nargs="?", default=None, help="optional mask MRC filename (voxels > 0.1 are included)"
)
janas_map_histogram.add_argument(
    "--bins", type=int, default=512, help="number of histogram bins (default: 512)"
)
janas_map_histogram.add_argument(
    "--density", action="store_true",
    help="normalise histogram to probability density (default: counts)"
)
janas_map_histogram.add_argument(
    "--clip", type=float, nargs=2, metavar=("LOW_PCT", "HIGH_PCT"),
    default=None,
    help="optional percentile clipping, e.g. --clip 0.5 99.5"
)


######################
## CryoSparc Utils
######################
## csparc_nurefinement
janas_csparc_nurefinement = command.add_parser(
    "csparc_nurefinement",
    description=(
        "Run CryoSPARC non-uniform refinement starting from a Relion STAR "
        "and update the STAR with refined angles and shifts."
    ),
    help="CryoSPARC non-uniform refinement + STAR update",
)

janas_csparc_nurefinement.add_argument(
    "input_star",
    help="input particles STAR file (Relion >=3.1)",
)

janas_csparc_nurefinement.add_argument(
    "output_basename",
    help=(
        "basename for outputs (without extension). "
        "The STAR will be written as <basename>.star and maps as "
        "<basename>_rec.mrc, <basename>_recH1.mrc, <basename>_recH2.mrc."
    ),
)

janas_csparc_nurefinement.add_argument(
    "--particle-dir",
    required=False,
    default=".",
    help="directory containing the particle stack(s) referenced by _rlnImageName",
)

janas_csparc_nurefinement.add_argument(
    "--user",
    required=False,
    type=str,
    help=(
        "CryoSPARC user e-mail. "
        "If omitted, it is read from CRYOSPARC_EMAIL or from "
        "~/.janas/cryosparc_config.json (created by `janas csparc_setup`)."
    ),
)
janas_csparc_nurefinement.add_argument(
    "--mask", default=None,
    help="optional static mask (.mrc) to feed into non-uniform refinement"
)
janas_csparc_nurefinement.add_argument(
    "--project",
    required=True,
    help="CryoSPARC project UID (e.g. P1)",
)

janas_csparc_nurefinement.add_argument(
    "--workspace",
    required=True,
    help="CryoSPARC workspace UID (e.g. W1)",
)

janas_csparc_nurefinement.add_argument(
    "--lane",
    required=True,
    help="CryoSPARC compute lane (as defined in CryoSPARC)",
)

janas_csparc_nurefinement.add_argument(
    "--sym",
    default="C1",
    help="symmetry to use in non-uniform refinement (default: C1)",
)

janas_csparc_nurefinement.add_argument(
    "--ref",
    default=None,
    help="optional reference map (.mrc); if omitted, ab-initio is used",
)

janas_csparc_nurefinement.add_argument(
    "--ini-high",
    dest="ini_high",
    type=float,
    default=None,
    help="initial resolution for refinement (Å)",
)

janas_csparc_nurefinement.add_argument(
    "--resplit",
    action="store_true",
    help="force re-splitting of half-sets in NU refinement",
)

janas_csparc_nurefinement.add_argument(
    "--local",
    action="store_true",
    help="run local refinement after non-uniform refinement",
)

janas_csparc_nurefinement.add_argument(
    "--min-angular-step",
    dest="min_angular_step",
    type=float,
    default=0.01,
    help="minimum angular step for local refinement (default: 0.01)",
)

janas_csparc_nurefinement.add_argument(
    "--loglevel",
    "-l",
    default="WARNING",
    help="logging level for STAR update (DEBUG, INFO, WARNING)",
)

janas_csparc_nurefinement.add_argument(
    "--precomputed",
    type=str,
    default=None,
    help=(
        "Reuse an existing CryoSPARC refinement job (e.g. J89) instead of "
        "running a new non-uniform (and optional local) refinement."
    ),
)

def _resolve_cryosparc_user(cli_user: Optional[str]) -> str:
    """
    returns e-mail of CryoSPARC user:
    """
    if cli_user:
        return cli_user

    env_email = os.getenv("CRYOSPARC_EMAIL")
    if env_email:
        return env_email

    cfg_path = Path.home() / ".janas" / "cryosparc_config.json"
    if cfg_path.is_file():
        try:
            with cfg_path.open("r") as f:
                cfg = json.load(f)
            email = cfg.get("email")
            if email:
                return email
        except Exception:
            pass  # non bloccare, si passa al fallback

    raise RuntimeError(
        "CryoSPARC user e-mail not specified. "
        "Use --user or run `janas csparc_setup` to configure it."
    )

def csparc_nurefinement_cli(args):
    """
    CLI entry point for 'csparc_nurefinement'.

    This is intentionally tolerant: if cryosparc-tools is not installed
    or a CryoSPARC instance is not available, a clear RuntimeError is raised
    and reported to the user, without breaking other JANAS commands.
    """
    from pathlib import Path
    in_star = Path(args.input_star)
    out_base = Path(args.output_basename)
    if out_base.suffix == ".star":
        out_base = out_base.with_suffix("")
    star_out = out_base.with_suffix(".star")

    try:
        # Lazy import, so that janas_cmd_utils can be imported without
        # requiring cryosparc-tools when the user does not use this command.
        from janas import utils_csparc
    except ImportError:
        print(
            "[csparc_nurefinement] CryoSPARC integration is not available: "
            "the Python package 'cryosparc-tools' is not installed.\n"
            "Install it with 'pip install cryosparc-tools' or avoid using "
            "the 'csparc_nurefinement' command."
        )
        raise SystemExit(1)

    try:
        user_email=_resolve_cryosparc_user(args.user)
        utils_csparc.run_csparc_nu_refinement_and_update_star(
            star_in=str(in_star),
            star_out=str(star_out),
            particle_dir=args.particle_dir,
            user_email=user_email,
            project_uid=args.project,
            workspace_uid=args.workspace,
            lane=args.lane,
            symmetry=args.sym,
            ref_map=args.ref,
            ini_high=args.ini_high,
            resplit=args.resplit,
            do_local=args.local,
            min_angular_step=args.min_angular_step,
            loglevel=args.loglevel,
            mask_map=args.mask,
            output_basename=str(out_base),
            precomputed_job_uid=args.precomputed,
        )
    except RuntimeError as e:
        # Qualsiasi problema lato CryoSPARC (connessione, autenticazione,
        # job falliti, ecc.) viene riportato con un messaggio leggibile.
        print(f"[csparc_nurefinement] {e}")
        raise SystemExit(1)



# ------------------------------------------------------------
# csparc_localnurefinement
# ------------------------------------------------------------
janas_csparc_localnurefinement = command.add_parser(
    "csparc_localnurefinement",
    help=(
        "Run CryoSPARC **local** refinement (new_local_refine) on a STAR file "
        "and update the STAR + export maps."
    ),
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

janas_csparc_localnurefinement.add_argument(
    "--particle-dir",
    required=False,
    default=".",
    type=str,
    help="Directory containing particle image files (RELION-style).",
)
janas_csparc_localnurefinement.add_argument(
    "--project",
    required=True,
    type=str,
    help="CryoSPARC project UID, e.g. P31.",
)
janas_csparc_localnurefinement.add_argument(
    "--workspace",
    required=True,
    type=str,
    help="CryoSPARC workspace UID, e.g. W1.",
)
janas_csparc_localnurefinement.add_argument(
    "--lane",
    required=True,
    type=str,
    help="CryoSPARC scheduler lane, e.g. 'default'.",
)
janas_csparc_localnurefinement.add_argument(
    "--sym",
    required=True,
    type=str,
    help="Symmetry string (e.g. C1, D7, etc.).",
)
janas_csparc_localnurefinement.add_argument(
    "--ref",
    required=True,
    type=str,
    help="Path to the reference map (.mrc) to import as initial volume.",
)
janas_csparc_localnurefinement.add_argument(
    "--mask",
    required=True,
    type=str,
    help="Path to the static mask (.mrc) to use in local refinement.",
)
janas_csparc_localnurefinement.add_argument(
    "--resplit",
    action="store_true",
    help=(
        "Force re-splitting particles into two GS halves "
        "(refine_gs_resplit=True)."
    ),
)
janas_csparc_localnurefinement.add_argument(
    "--min-angular-step",
    type=float,
    default=0.2,
    help="Maximum alignment resolution (degrees) for local refinement.",
)
# user is optional: if omitted, we read it from ~/.janas/cryosparc_config.json
janas_csparc_localnurefinement.add_argument(
    "--user",
    required=False,
    type=str,
    help="CryoSPARC user e-mail. If omitted, read from csparc_setup config.",
)
janas_csparc_localnurefinement.add_argument(
    "--loglevel",
    required=False,
    type=str,
    default="info",
    help="Logging level for this command.",
)
janas_csparc_localnurefinement.add_argument(
    "input_star",
    type=str,
    help="Input RELION STAR file to be updated.",
)
janas_csparc_localnurefinement.add_argument(
    "output_basename",
    type=str,
    help=(
        "Basename for outputs. STAR will be '<basename>.star'; "
        "maps will be '<basename>_rec.mrc', '_recH1.mrc', '_recH2.mrc'."
    ),
)
janas_csparc_localnurefinement.add_argument(
    "--precomputed",
    type=str,
    default=None,
    help="Reuse an existing CryoSPARC refinement job (e.g. J89) instead of running a new local refinement.",
)



def csparc_localnurefinement_cli(args):
    """
    Pure local refinement:
      import_particles -> import_ref -> import_mask -> new_local_refine
      then update STAR and export maps.
    """
    import janas.utils_csparc as utils_csparc  # or from janas import utils_csparc

    user_email = _resolve_cryosparc_user(args.user)

    utils_csparc.run_csparc_local_refinement_and_update_star(
        star_in=args.input_star,
        star_out=args.output_basename + ".star",
        output_basename=args.output_basename,
        particle_dir=args.particle_dir,
        user_email=user_email,
        project_uid=args.project,
        workspace_uid=args.workspace,
        lane=args.lane,
        symmetry=args.sym,
        ref_map=args.ref,
        mask_map=args.mask,
        resplit=args.resplit,
        min_angular_step=args.min_angular_step,
        loglevel=args.loglevel,
        precomputed_job_uid=args.precomputed,
    )


######################
## FSC Utils



# FSC HELPERS:
def _safe_angpix(spacing):
    """Normalise utils.get_MRC_map_pixel_spacing() return shape to a float Å/px."""
    try:
        return float(spacing[0])
    except (TypeError, IndexError):
        return float(spacing)



def _load_mrc_as_array(fname: str) -> Tuple[np.ndarray, Tuple[int, int, int], float]:
    """
    Read an MRC via janas_core.* and return:
      vol as float32 shaped (Z, Y, X) = (nz, ny, nx),
      shape tuple (nz, ny, nx),
      angpix (Å/px).
    """
    if not os.path.exists(fname):
        raise FileNotFoundError(f"File not found: {fname}")

    # In this codebase sizeMRC is used as header order (nx, ny, nz)
    nx, ny, nz = janas_core.sizeMRC(fname)

    flat = np.array(janas_core.ReadMRC(fname), dtype=np.float32)

    # MRC is X-fast: reshape as (nz, ny, nx) so X is the last NumPy axis
    vol = flat.reshape((nz, ny, nx), order="C")

    angpix = _safe_angpix(utils.get_MRC_map_pixel_spacing(fname))
    return vol, (nz, ny, nx), float(angpix)


# CLI PARSER
janas_fsc = command.add_parser(
    "fsc",
    description="Compute Fourier Shell Correlation (FSC) between one or more half-map pairs.",
    help="compute FSC between half-map pairs"
)
# Allow multiple --hmaps, each expecting exactly two filenames
janas_fsc.add_argument(
    "--hmaps", action="append", nargs=2, metavar=("HALF1.mrc", "HALF2.mrc"),
    required=True,
    help="pair of half-maps (.mrc). Use multiple --hmaps to overlay multiple FSC curves."
)
janas_fsc.add_argument(
    "--bins", type=int, default=None,
    help="number of frequency shells (default: half of smallest volume dimension)"
)
janas_fsc.add_argument(
    "--threshold", type=float, default=0.143,
    help="FSC threshold for resolution estimate (default: 0.143)"
)
janas_fsc.add_argument(
    "--output_prefix", default=None,
    help="if set, save each FSC curve to CSV as {output_prefix}_set{N}.csv"
)
janas_fsc.add_argument(
    "--no-plot", action="store_true",
    help="suppress display of the FSC plot"
)
def fsc_utils(args):
    curves = []
    for idx_pair, pair in enumerate(args.hmaps, start=1):
        f1, f2 = pair
        vol1, shape1, ang1 = _load_mrc_as_array(f1)
        vol2, shape2, ang2 = _load_mrc_as_array(f2)

        if shape1 != shape2:
            raise ValueError(f"Shape mismatch for pair {idx_pair}: {shape1} vs {shape2}")
        if not np.isclose(ang1, ang2, rtol=0, atol=1e-6):
            raise ValueError(f"Å/px mismatch for pair {idx_pair}: {ang1} vs {ang2}")

        #freqs, fsc_vals = utils.compute_fsc_3d(vol1, vol2, ang1, n_bins=args.bins)
        # Auto-detect 2D (one singleton dimension) vs 3D
        active_dims = sum(s > 1 for s in vol1.shape)
        if active_dims == 2:
            img1 = np.squeeze(vol1)
            img2 = np.squeeze(vol2)
            if img1.ndim != 2 or img2.ndim != 2:
                raise ValueError(f"Detected 2D but squeeze did not yield 2D arrays: {vol1.shape} -> {img1.shape}")
            freqs, fsc_vals = utils.compute_fsc_2d(img1, img2, ang1, n_bins=args.bins)
        elif active_dims == 3:
            freqs, fsc_vals = utils.compute_fsc_3d(vol1, vol2, ang1, n_bins=args.bins)
        else:
            raise ValueError(f"Unsupported dimensionality for FSC/FRC: vol shape {vol1.shape}")


        # Optional CSV per set
        if args.output_prefix:
            out_csv = f"{args.output_prefix}_set{idx_pair}.csv"
            np.savetxt(
                out_csv,
                np.column_stack((freqs, fsc_vals)),
                header="Frequency(1/Å)    FSC",
                fmt="%.6f    %.6f"
            )
            print(f"Saved FSC curve to {out_csv}")

        # Threshold crossing summary
        f_cross, res = utils.find_resolution_at_threshold(freqs, fsc_vals, args.threshold)
        if f_cross is not None:
            print(
                f"Set {idx_pair}: FSC={args.threshold:.3f} crosses at "
                f"{f_cross:.6f} Å^-1 -> {res:.3f} Å"
            )
        else:
            print(f"Set {idx_pair}: no crossing found for FSC={args.threshold:.3f}")

        label = f"{os.path.basename(f1)}\n{os.path.basename(f2)}"
        curves.append({
            "freqs": freqs, "fsc": fsc_vals, "label": label,
            "f_cross": f_cross, "res": res, "angpix": ang1
        })

    # Plot overlay unless suppressed
    if not args.no_plot and curves:
        plt.figure()

        # Plot curves, track minima for y-limits, and remember line colours
        y_mins = []
        line_info = []  # (freqs, fsc, color)
        for c in curves:
            y = c["fsc"]
            line, = plt.plot(c["freqs"], y, lw=2.2, label=c["label"])  # thicker line
            color = line.get_color()
            line_info.append((c["freqs"], y, color))
            finite = np.isfinite(y)
            if np.any(finite):
                y_mins.append(np.min(y[finite]))

        # Reference horizontal lines (no legend entries)
        plt.hlines(0.5, 0, 0.6, linestyles='--', colors='gray')
        plt.hlines(0.143, 0, 0.6, linestyles='--', colors='gray')
        # Zero baseline
        plt.hlines(0.0, 0, 0.6, linestyles='-', colors='black', linewidth=0.8)

        # Axes labels and limits
        plt.xlabel('Spatial frequency (Å⁻¹)')
        plt.ylabel('FSC')
        plt.title('Fourier Shell Correlation')
        plt.legend(loc='upper right', fontsize=8)
        plt.grid(True)
        plt.xlim(0, 0.6)

        # Y axis: set lower bound just below observed min, small margin
        if y_mins:
            min_val = min(y_mins)
        else:
            min_val = 0.0

        # Small visual margin
        margin = 0.03

        # If the minimum is non-negative, start exactly at 0.0; otherwise just a hair below the min
        if min_val >= 0.0:
            lower = 0.0
        else:
            lower = min_val - margin

        upper = 1.05
        plt.ylim(lower, upper)

        # Ensure y-ticks include 0, 0.143, 0.5 (and are within limits)
        yticks = list(plt.yticks()[0]) + [0.0, 0.143, 0.5]
        lo, hi = plt.ylim()
        yticks = sorted({t for t in yticks if lo <= t <= hi})
        plt.yticks(yticks)

        # --- Vertical dotted lines at crossings (extend to the X axis baseline, i.e., y=lower) ---
        thresholds_for_vlines = (0.5, 0.143)
        y_bottom = lower  # draw down to the plot's bottom (x-axis baseline)
        for (freqs, fsc_vals, color) in line_info:
            for thr in thresholds_for_vlines:
                f_cross, _ = utils.find_resolution_at_threshold(freqs, fsc_vals, thr)
                if f_cross is not None and 0 <= f_cross <= 0.6:
                    plt.vlines(f_cross, y_bottom, thr, linestyles=':', colors=color, linewidth=1.0)

        plt.tight_layout()
        plt.show()



def _make_box_mask_from_mrc(mask_path: str,
                            expected_shape: Tuple[int, int, int],
                            pad_vox: int = 3) -> np.ndarray:
    nx, ny, nz = janas_core.sizeMRC(mask_path)
    mask_shape = (nz, ny, nx)

    if mask_shape != tuple(expected_shape):
        raise ValueError(f"Box-mask shape mismatch: mask={mask_shape} vs half-maps={expected_shape}")

    mask_vol = np.array(janas_core.ReadMRC(mask_path), dtype=np.float32).reshape(mask_shape, order="C")
    inside = (mask_vol >= 0.5)

    if not np.any(inside):
        raise ValueError(f"Mask contains no voxels >= 0.5: {mask_path}")

    z_idx, y_idx, x_idx = np.where(inside)

    z0 = max(int(z_idx.min()) - pad_vox, 0)
    z1 = min(int(z_idx.max()) + pad_vox + 1, mask_shape[0])
    y0 = max(int(y_idx.min()) - pad_vox, 0)
    y1 = min(int(y_idx.max()) + pad_vox + 1, mask_shape[1])
    x0 = max(int(x_idx.min()) - pad_vox, 0)
    x1 = min(int(x_idx.max()) + pad_vox + 1, mask_shape[2])

    box_mask = np.zeros(mask_shape, dtype=np.float32)
    box_mask[z0:z1, y0:y1, x0:x1] = 1.0
    return box_mask


def _make_box_mask_pair_from_mrc(mask_path: str,
                                 expected_shape: Tuple[int, int, int],
                                 pad_vox: int = 3,
                                 compute_pad_vox: int = 1) -> Tuple[np.ndarray, np.ndarray]:
    output_box = _make_box_mask_from_mrc(mask_path, expected_shape, pad_vox=pad_vox)
    compute_box = _make_box_mask_from_mrc(mask_path, expected_shape, pad_vox=pad_vox + compute_pad_vox)
    return output_box, compute_box


def _make_bulk_box_mask_from_first_pair(prefix: str,
                                        items_csv: str,
                                        suffixes_csv: Optional[str],
                                        mask_path: str,
                                        pad_vox: int = 3) -> np.ndarray:
    if not suffixes_csv:
        suffixes_csv = "_recH1.mrc,_recH2.mrc"

    suffix_parts = [s.strip() for s in suffixes_csv.split(",") if s.strip()]
    if len(suffix_parts) != 2:
        raise ValueError(
            "[locresBulk] suffixesCSV must contain exactly two comma-separated entries, "
            f"e.g. '_recH1.mrc,_recH2.mrc'; got: {suffixes_csv!r}"
        )
    suffix1, suffix2 = suffix_parts

    item_strs = [s.strip() for s in items_csv.split(",") if s.strip()]
    if not item_strs:
        raise ValueError("[locresBulk] No items provided.")

    labels = []
    for s in item_strs:
        try:
            labels.append((int(s), s))
        except Exception:
            raise ValueError(f"[locresBulk] Cannot interpret item '{s}' as an integer.")

    labels_sorted = sorted(labels, key=lambda t: t[0], reverse=True)

    for _numeric_label, label_str in labels_sorted:
        half1 = f"{prefix}{label_str}{suffix1}"
        half2 = f"{prefix}{label_str}{suffix2}"
        if Path(half1).exists() and Path(half2).exists():
            nx, ny, nz = janas_core.sizeMRC(half1)
            return _make_box_mask_from_mrc(mask_path, (nz, ny, nx), pad_vox=pad_vox)

    raise ValueError(
        "[locresBulk] No existing half-map pair was found to determine the box-mask shape."
    )
    
    

# ---------------------------------
# Local resolution (half-map pair)
janas_locres = command.add_parser(
    "locres",
    description="Compute local-resolution map from two half-maps (H1/H2). Writes *_locres.mrc and auxiliary files.",
    help="local-resolution from half-maps"
)
janas_locres.add_argument("halfmap1", type=str, help="first half-map (.mrc)")
janas_locres.add_argument("halfmap2", type=str, help="second half-map (.mrc)")
janas_locres.add_argument("--mask", default=None, help="optional mask .mrc used to build a rectangular box mask (voxels >=0.5 define the box)")
janas_locres.add_argument("--tight-mask", default=None, help="optional mask .mrc used directly as a tight mask (voxels >=0.5 included)")
janas_locres.add_argument("--o", dest="output_basename", required=True, help="output basename (no extension)")
janas_locres.add_argument("--threshold", type=float, default=0.143, help="local FSC threshold (default 0.143)")
janas_locres.add_argument("--sampling", default="auto", help="'auto' or Å value for sampling step (default auto)")
janas_locres.add_argument("--radius", default="auto", help="'auto' or Å radius for local windows (default auto)")
janas_locres.add_argument("--edgewidth", type=float, default=None, help="edge width (Å); default 0.3*radius")
janas_locres.add_argument("--cycles", type=float, default=10.0, help="cycles for auto-radius (default 10)")
janas_locres.add_argument("--gamma", type=float, default=1.8, help="gamma for auto-radius (default 1.8)")
janas_locres.add_argument("--bins", type=int, default=None, help="FSC bins (default: inferred)")
janas_locres.add_argument("--rand-global", default="auto", help="'auto'|'none'|Å for global phase randomisation")
janas_locres.add_argument("--rand-local", default="auto", help="'auto'|'none'|Å for local phase randomisation")
janas_locres.add_argument("--rand-frac-nyq", type=float, default=None, help="override global k-cut as fraction of Nyquist (0..1)")
janas_locres.add_argument("--cpu", type=int, default=1, help="CPU processes (default 1)")
janas_locres.add_argument("--resample", action="store_true", help="write coarse-grid map (resampled outputs)")
janas_locres.add_argument("--interp", choices=("linear","cubic"), default="cubic", help="coarse-grid interpolation")
janas_locres.add_argument("--plot", action="store_true", help="plot global FSC")
janas_locres.add_argument("--accurate", action="store_true", help="use the previous, slower local-resolution parameters instead of the new faster defaults")
janas_locres.add_argument("--fast", action="store_true", help=argparse.SUPPRESS)
janas_locres.add_argument("--box-pad", type=int, default=0, help="padding (voxels) added to the mask-derived output box when using --mask")
janas_locres.add_argument("--box-compute-pad", type=int, default=2, help="extra padding (voxels) added only for computation before trimming back to the output box")

def locres_cli(args):
    if args.mask is not None and args.tight_mask is not None:
        raise ValueError("Use either --mask or --tight-mask, not both.")

    mask_arg = None
    output_mask_arg = None
    compute_only_mask = False
    mask_dilation_vox = 6

    if args.mask is not None:
        nx, ny, nz = janas_core.sizeMRC(args.halfmap1)
        output_mask_arg, mask_arg = _make_box_mask_pair_from_mrc(
            args.mask,
            (nz, ny, nx),
            pad_vox=int(args.box_pad),
            compute_pad_vox=int(args.box_compute_pad),
        )
        compute_only_mask = True
        mask_dilation_vox = 0
    elif args.tight_mask is not None:
        mask_arg = args.tight_mask

    use_fast_defaults = (not bool(getattr(args, "accurate", False))) or bool(getattr(args, "fast", False))

    cfg = locres_utils.LocResConfig(
        halfmap1=args.halfmap1,
        halfmap2=args.halfmap2,
        mask=mask_arg,
        output_mask=output_mask_arg,
        compute_only_mask=compute_only_mask,
        mask_dilation_vox=mask_dilation_vox,
        local_fsc_threshold=float(args.threshold),
        sampling=args.sampling,
        radius=args.radius,
        edgewidth=args.edgewidth,
        cycles=float(args.cycles),
        gamma=float(args.gamma),
        bins=args.bins,
        rand_global=args.rand_global,
        rand_local=args.rand_local,
        rand_frac_nyq=args.rand_frac_nyq,
        cpu=int(args.cpu),
        resample=bool(args.resample),
        interp=args.interp,
        output_basename=args.output_basename,
        plot=bool(args.plot),
        fast=bool(use_fast_defaults),
    )
    _ = locres_utils.locres_map(cfg)

# ---------------------------------
# Local resolution (bulk over multiple half-map pairs)
janas_locresBulk = command.add_parser(
    "locresBulk",
    description=(
        "Compute local-resolution maps for multiple half-map pairs that share a "
        "common prefix and numeric labels. The pair with the largest label is used "
        "to determine sampling/radius/global FSC parameters, which are then reused "
        "for all other pairs."
    ),
    help="local-resolution for multiple half-map pairs (shared prefix + labels)"
)

janas_locresBulk.add_argument(
    "prefix",
    type=str,
    help="common prefix for half-map filenames, e.g. 'ciccio_'"
)
janas_locresBulk.add_argument(
    "itemsCSV",
    type=str,
    help="comma-separated list of numeric labels, e.g. '34,455,2233,4455'"
)
janas_locresBulk.add_argument(
    "suffixesCSV",
    nargs="?",
    default="_recH1.mrc,_recH2.mrc",
    help=(
        "optional comma-separated pair of suffixes for half-maps "
        "(default: '_recH1.mrc,_recH2.mrc')"
    )
)


janas_locresBulk.add_argument(
    "--mask", default=None,
    help="optional mask .mrc used to build a rectangular box mask (voxels >=0.5 define the box)"
)
janas_locresBulk.add_argument(
    "--tight-mask", default=None,
    help="optional mask .mrc used directly as a tight mask (voxels >=0.5 included)"
)

janas_locresBulk.add_argument(
    "--threshold", type=float, default=0.143,
    help=(
        "Local FSC threshold used to convert FSC→resolution (default 0.143). "
        "Common choices: 0.143 (half-map), 0.5 (more conservative)."
    )
)

janas_locresBulk.add_argument(
    "--sampling", default="auto",
    help=(
        "Spacing (Å) between local-resolution evaluation centres. "
        "'auto' uses a map-dependent heuristic. "
        "Smaller values give denser sampling and slower runs. "
        "Typical range: 1.5–4.0 Å."
    )
)
janas_locresBulk.add_argument(
    "--radius", default="auto",
    help=(
        "Radius (Å) of the spherical local-FSC window. "
        "'auto' derives a radius from global resolution and the "
        "--gamma/--cycles settings. "
        "Larger radii improve statistical stability but reduce locality "
        "and increase runtime. Typical range: ~15–35 Å."
    )
)

janas_locresBulk.add_argument(
    "--edgewidth", type=float, default=None,
    help=(
        "Raised-cosine edge width (Å) of the spherical window. "
        "If not set, defaults to ~0.3 × radius. "
        "Usually keep within ~0.1–0.5 × radius."
    )
)

janas_locresBulk.add_argument(
    "--cycles", type=float, default=15.0,
    help=(
        "Controls the effective window size when --radius is 'auto'. "
        "Higher values typically increase the auto radius and smooth local FSC. "
        "Practical range: ~10–25."
    )
)

janas_locresBulk.add_argument(
    "--gamma", type=float, default=1.5,
    help=(
        "Scaling factor for --radius 'auto'. "
        "Higher gamma increases the local window radius. "
        "Practical range: ~1.0–2.0."
    )
)

janas_locresBulk.add_argument(
    "--bins", type=int, default=None,
    help=(
        "Number of radial frequency bins used to compute each local FSC curve. "
        "If omitted, inferred from the window size. "
        "More bins give finer FSC sampling but can be noisier in small windows."
    )
)

janas_locresBulk.add_argument(
    "--rand-global", default="auto",
    help=(
        "Global phase-randomisation cutoff. "
        "'auto' selects a cutoff from the reference pair's global FSC "
        "and clamps to a Nyquist fraction. "
        "'none' disables global randomisation. "
        "A numeric value is interpreted as a resolution in Å; phases "
        "beyond that resolution are randomised."
    )
)

janas_locresBulk.add_argument(
    "--rand-local", default="auto",
    help=(
        "Local phase-randomisation cutoff applied to each spherical window. "
        "'auto' follows the reference-pair-derived cutoff used for bulk, "
        "providing a steady k-cut across particle subsets. "
        "'none' disables local randomisation. "
        "A numeric value is interpreted as a resolution in Å."
    )
)

janas_locresBulk.add_argument(
    "--rand-frac-nyq", type=float, default=None,
    help=(
        "Override phase-randomisation cutoff as a fraction of Nyquist frequency "
        "(0 < value < 1). Use to enforce a fixed band across runs "
        "(e.g., 0.6–0.8). When set, this override is applied consistently "
        "in bulk."
    )
)

janas_locresBulk.add_argument(
    "--cpu", type=int, default=1,
    help=(
        "Number of CPU worker processes (default 1). "
        "Increase for faster processing of many local regions."
    )
)

janas_locresBulk.add_argument(
    "--resample", action="store_true",
    help=(
        "Write coarse-grid outputs at the chosen --sampling spacing "
        "in addition to the full-grid maps. Useful for speed/diagnostics."
    )
)

janas_locresBulk.add_argument(
    "--interp", choices=("linear", "cubic"), default="cubic",
    help=(
        "Interpolation method used to expand coarse-grid local-resolution "
        "estimates back to the full voxel grid when needed."
    )
)

janas_locresBulk.add_argument(
    "--plot", action="store_true",
    help=(
        "Plot the global FSC curve for the reference pair used to set "
        "auto parameters in bulk."
    )
)
janas_locresBulk.add_argument(
    "--accurate", action="store_true",
    help="use the previous, slower local-resolution parameters instead of the new faster defaults"
)
janas_locresBulk.add_argument(
    "--fast", action="store_true",
    help=argparse.SUPPRESS
)
janas_locresBulk.add_argument(
    "--box-pad", type=int, default=0,
    help="padding (voxels) added to the mask-derived output box when using --mask"
)
janas_locresBulk.add_argument(
    "--box-compute-pad", type=int, default=2,
    help="extra padding (voxels) added only for computation before trimming back to the output box"
)
janas_locresBulk.add_argument(
    "--include-global-analysis",
    metavar="CSV",
    dest="include_global_analysis",
    default=None,
    help=(
        "If set, compute global half-map quality metrics for each pair "
        "and write them to the given CSV file. "
        "Columns: numParticles,FSC,FSC_score,ampl_falloff_score,PSF_score,global_score."
    ),
)




def locresBulk_cli(args):
    """
    CLI entry point for 'janas_utils locresBulk'.

    Example:
        janas_utils locresBulk ciccio_ 34,455,2233,4455 _recH1.mrc,_recH2.mrc

    will look for:
        ciccio_34_recH1.mrc / ciccio_34_recH2.mrc
        ciccio_455_recH1.mrc / ciccio_455_recH2.mrc
        ciccio_2233_recH1.mrc / ciccio_2233_recH2.mrc
        ciccio_4455_recH1.mrc / ciccio_4455_recH2.mrc

    and use the largest existing label as reference for the shared locres parameters.
    """
    if args.mask is not None and args.tight_mask is not None:
        raise ValueError("Use either --mask or --tight-mask, not both.")

    args.output_mask = None
    args.fast = (not bool(getattr(args, "accurate", False))) or bool(getattr(args, "fast", False))

    if args.mask is not None:
        output_box, compute_box = _make_box_mask_pair_from_mrc(
            args.mask,
            _make_bulk_box_mask_from_first_pair(
                prefix=args.prefix,
                items_csv=args.itemsCSV,
                suffixes_csv=args.suffixesCSV,
                mask_path=args.mask,
                pad_vox=0,
            ).shape,
            pad_vox=int(args.box_pad),
            compute_pad_vox=int(args.box_compute_pad),
        )
        args.output_mask = output_box
        args.mask = compute_box
    elif args.tight_mask is not None:
        args.mask = args.tight_mask
    else:
        args.mask = None

    locres_utils.run_locres_bulk(
        prefix=args.prefix,
        items_csv=args.itemsCSV,
        suffixes_csv=args.suffixesCSV,
        args=args,
    )


# ---------------------------------
# Local-resolution statistics over existing locres maps
janas_locresStats = command.add_parser(
    "locresStats",
    description=(
        "Compute masked summary statistics (max/Q1/mean/Q3/min) for one or more "
        "local-resolution maps, using a supplied mask. The output is a CSV "
        "similar to bestRanked_locres_values.csv."
    ),
    help="summary stats for locres maps"
)

janas_locresStats.add_argument(
    "--locres-files", "-L",
    nargs="+",
    required=True,
    help=(
        "one or more locres .mrc files or shell-style patterns (quoted), "
        "e.g. 'norm__*_locres.mrc'"
    ),
)

janas_locresStats.add_argument(
    "--mask",
    required=True,
    help="mask MRC file; voxels with mask > 0.2 are included in the statistics",
)
# add near other locresStats arguments
janas_locresStats.add_argument(
    "--assessmentMethod",
    required=False,
    default="mean",
    choices=["mean", "median"],
    help="Statistic used for the central-tendency column in --out-stats CSV (mean or median).",
)
janas_locresStats.add_argument(
    "--out-stats",
    default="bestRanked_locres_values.csv",
    help="output CSV filename (default: bestRanked_locres_values.csv)",
)
janas_locresStats.add_argument(
    "--out-local-best-particles",
    default="",
    help="mrc file with the local-best-particles estimation (ignored if empty)",
)
def locresStats_cli(args):
    """
    CLI entry point for 'janas_utils locresStats'.

    Example:
        janas_utils locresStats \
            --locres-files 'norm__janas_SCI__1.0_scored_selection_1_best*_locres.mrc' \
            --mask masks/mask_partB.mrc \
            --out-stats class1_selectionFull/.../bestRanked_locres_values.csv
    """
    import glob

    patterns = args.locres_files
    collected: List[str] = []

    for pat in patterns:
        matches = sorted(glob.glob(pat))
        if not matches:
            print(f"[locresStats] WARNING: pattern '{pat}' matched no files")
            continue
        collected.extend(matches)

    if not collected:
        raise SystemExit(
            "[locresStats] No locres files found for the given --locres-files patterns."
        )

    # Remove duplicates while preserving order
    seen = set()
    locres_list: List[str] = []
    for f in collected:
        if f in seen:
            continue
        seen.add(f)
        locres_list.append(f)

    locres_utils.run_locres_stats(
        locres_files=locres_list,
        mask_path=args.mask,
        out_csv=args.out_stats,
        out_LocalMinBestParticles_map=args.out_local_best_particles,
        assessmentMethod=args.assessmentMethod,
    )


######################
######################
## project_map (from STAR)
janas_project_map = command.add_parser(
    "project_map",
    description=("Project a 3D reference map at each particle pose and write a 2D reprojection "
                 "stack (.mrcs) plus a STAR that repoints _rlnImageName to the new stack."),
    help="project 3D map to per-particle 2D re-projections and build a .mrcs + .star"
)
janas_project_map.add_argument(
    "input_star", help="input .star with particle poses (_rlnImageName, angles, origins)"
)
janas_project_map.add_argument(
    "output_root",
    help="output rootname (writes output_root.mrcs and output_root.star)"
)
janas_project_map.add_argument(
    "--map", required=True,
    help="3D reference map .mrc to project"
)
janas_project_map.add_argument(
    "--mask", required=False,
    help="3D reference map .mrc to project"
)
janas_project_map.add_argument(
    "--root", default=None,
    help=("Base directory to resolve relative stack paths in _rlnImageName "
          "(useful for CryoSPARC-style relative paths)")
)

def project_map_utils(args):
    projector_utils.project_map_from_star(
        star_in=args.input_star,
        out_root=args.output_root,
        map_3d=args.map,
        mask_3d=args.mask,
        project_root=(os.path.abspath(args.root) if args.root else None),
    )


##########################################################
##### subtract_particles
##########################################################
janas_subtract_particles = command.add_parser(
    "subtract_particles",
    description=("subtract particles."),
    help="subtract particles"
)
janas_subtract_particles.add_argument(
    "input_star", help="input .star with particle poses (_rlnImageName, angles, origins)"
)
janas_subtract_particles.add_argument(
    "output_root",
    help="output rootname (writes output_root.mrcs and output_root.star)"
)
janas_subtract_particles.add_argument(
    "--map", required=True,
    help="3D reference map .mrc for subtraction"
)
janas_subtract_particles.add_argument(
    "--map2", required=False,
    help="3D reference map .mrc second halfMap for subtraction"
)
janas_subtract_particles.add_argument(
    "--mask_to_retain", required=True,
    help="mask 3d retain"
)
janas_subtract_particles.add_argument(
    "--root", default=None,
    help=("Base directory to resolve relative stack paths in _rlnImageName "
          "(useful for CryoSPARC-style relative paths)")
)

def subtract_particles_utils(args):
    projector_utils.subtractParticles_weighted_relionLike(
        particlesStarFile=args.input_star,
        outputBasename=args.output_root,
        referenceMap=args.map,
        referenceMask=args.mask_to_retain,
    )


######################
## Main

def main(command_line=None):
    args = janas_parser.parse_args(command_line)
    if args.command == "split_mask":
        split_mask_utils(args)
    elif args.command == "maskedCrop":
        maskedCrop(args)
    elif args.command == "equalize_images":
        equalize_images(args)
    elif args.command == "randomize_halves":
        randomize_halves_utils(args)
    elif args.command == "scores_to_csv":
        scores_to_csv(args)
    elif args.command == "extract_particles_from_label_value":
        extract_particles_from_label_value(args)
    elif args.command == "csparc2star":
        csparc2star_utils(args)
    elif args.command == "csparc2star-stack":
        csparc2star_stack_utils(args)
    elif args.command == "create_stack":
        create_stack_utils(args)
    elif args.command == "backmap_stars":
        backmap_stars_utils(args)
    elif args.command == "update_from_csparc":
        janas_update_from_csparc_utils(args)
    elif args.command == "angpix":
        print_angpix(args)
    elif args.command == "mask":
        if args.mask_command == "merge":
            mask_merge_utils(args)
        elif args.mask_command == "countPixels":
            mask_countPixels_utils(args)
        else:
            janas_parser.print_help()
    elif args.command == "clip":
        if args.clip_command == "add":
            clip_add_maps(args)
        elif args.clip_command == "average":
            clip_average_maps(args)
        elif args.clip_command == "bfac":
            clip_bfac_utils(args)
        elif args.clip_command == "blur":
            clip_blur_utils(args)
        elif args.clip_command == "rot2D":
            clip_rot2D_utils(args)
        elif args.clip_command == "amplReplace":
            clip_amplReplace_utils(args)
        elif args.clip_command == "roundMask":
            clip_roundMask_utils(args)
        else:
            janas_parser.print_help()
    elif args.command == "map_histogram":
        utils.map_histogram_plot(args)
    elif args.command == "locres":
        locres_cli(args)
    elif args.command == "locresBulk":
        locresBulk_cli(args)
    elif args.command == "locresStats":
        locresStats_cli(args)
    elif args.command == "fsc":
        fsc_utils(args)
    elif args.command == "sigma_estimate":
        sigma_estimate_utils(args)        
    elif args.command == "ssnr":
        ssnr_utils(args)
    elif args.command == "quality_analysis":
        quality_analysis(args)
    elif args.command == "project_map":
        project_map_utils(args)
    elif args.command == "subtract_particles":
        subtract_particles_utils(args)
    elif args.command == "csparc_nurefinement":
        csparc_nurefinement_cli(args)        
    elif args.command == "csparc_localnurefinement":
        csparc_localnurefinement_cli(args)
    elif args.command == "compare_maps":
        compare_maps_cli(args)
    elif args.command == "compare2D":
        compare2D_cli(args)
    else:
        janas_parser.print_help()


if __name__ == "__main__":
    main()
