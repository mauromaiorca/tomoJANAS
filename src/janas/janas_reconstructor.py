#!/usr/bin/env python3
import sys
import os
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

import json
import threading
from collections import defaultdict

import numpy as np
from numpy.fft import fftn, ifftn, fftfreq
from scipy.spatial.transform import Rotation as R
from scipy.ndimage import gaussian_filter
from scipy.ndimage import binary_dilation, zoom as ndi_zoom, distance_transform_edt
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import multiprocessing as mp

import janas.janas_core as janas_core
from janas import assessParticles, starHandler, IO_utils, reconstructor_utils
from janas import locres_utils
from janas import janas_mapProcess

# ======================= STAR I/O (Relion v3.1) =======================

def write_relion31_star(out_path: str, optics_df, particles_df,
                        float_fmt: str = "{:.6f}", int_fmt: str = "{:d}"):
    def _write_loop(fh, df):
        fh.write("loop_\n")
        cols = list(df.columns)
        for i, c in enumerate(cols, start=1):
            fh.write(f"{c} #{i}\n")
        for _, row in df.iterrows():
            vals = []
            for c in cols:
                v = row[c]
                if isinstance(v, (int, np.integer)):
                    vals.append(int_fmt.format(int(v)))
                elif isinstance(v, (float, np.floating)):
                    vals.append(float_fmt.format(float(v)))
                else:
                    vals.append(str(v))
            fh.write(" ".join(vals) + "\n")

    with open(out_path, "w") as fh:
        fh.write("data_optics\n\n")
        _write_loop(fh, optics_df)
        fh.write("\n")
        fh.write("data_particles\n\n")
        _write_loop(fh, particles_df)

# ======================= STAR parsing (relion_v31) =======================

def parse_star_with_ctf_and_labels(star_file: str, use_ctf: bool,
                                   sort_label: Optional[str] = None,
                                   sort_ascending: bool = True):
    version = starHandler.infoStarFile(star_file)[2]
    if version != "relion_v31":
        sys.exit("Compatibility error: STAR file version is not 'relion_v31'.")

    optics = starHandler.dataOptics(star_file)
    if optics is None or "_rlnImagePixelSize" not in optics.columns:
        sys.exit("Could not find _rlnImagePixelSize in data_optics.")
    pixel_size = float(optics["_rlnImagePixelSize"].iloc[0])

    # Header columns to validate/sanitize sort label and to include all needed fields
    hdr_cols = starHandler.header_columns(star_file)
    norm_map = {str(c).strip(): str(c) for c in hdr_cols}  # normalized -> actual name

    # Base particle columns
    base_cols = [
        "_rlnAngleRot", "_rlnAngleTilt", "_rlnAnglePsi",
        "_rlnOriginXAngst", "_rlnOriginYAngst",
        "_rlnImageName", "_rlnOpticsGroup",
        "_rlnClassNumber", "_rlnRandomSubset"
    ]
    part_cols = base_cols.copy()

    # CTF columns
    if use_ctf:
        ctf_cols = [
            "_rlnDefocusU", "_rlnDefocusV", "_rlnDefocusAngle",
            "_rlnCtfBfactor"
        ]
        ctf_cols.append("_rlnPhaseShift" if "_rlnPhaseShift" in hdr_cols else None)
        part_cols += [c for c in ctf_cols if c is not None]

    # Optional sort column: validate against header & include in read
    sort_col_real = None
    if sort_label is not None:
        key = str(sort_label).strip()
        if key not in norm_map:
            sys.exit(f"ERROR: sort label '{sort_label}' not found in STAR.")
        sort_col_real = norm_map[key]
        if sort_col_real not in part_cols:
            part_cols.append(sort_col_real)

    # Read particles with the required columns
    df = starHandler.readColumns(star_file, part_cols)
    if df is None or df.shape[0] == 0:
        sys.exit("Failed to read required columns from particles section.")

    # Merge optics if CTF requested
    if use_ctf:
        optics_needed = [
            "_rlnImagePixelSize", "_rlnVoltage",
            "_rlnAmplitudeContrast", "_rlnSphericalAberration"
        ]
        df = df.merge(
            optics[["_rlnOpticsGroup"] + optics_needed],
            on="_rlnOpticsGroup", how="left", validate="many_to_one"
        )
        if df[optics_needed].isnull().any().any():
            sys.exit("Missing some optics parameters after merge.")

    # Apply sorting if requested
    if sort_col_real is not None:
        try:
            df[sort_col_real] = df[sort_col_real].astype(float)
        except Exception:
            sys.exit(f"ERROR: sort label '{sort_label}' cannot be parsed as float for numeric sorting.")
        df = df.sort_values(by=sort_col_real, ascending=sort_ascending, kind="mergesort")
        #print(f"SORT: particles sorted by '{sort_col_real}' ({'ascending' if sort_ascending else 'descending'}).")

    # Parse image names
    slices = df["_rlnImageName"].str.split("@", n=1, expand=True)
    slice_idx = slices[0].astype(int).to_numpy() - 1
    stack_names = slices[1].to_numpy()

    entries = list(zip(stack_names.tolist(), slice_idx.tolist()))
    params = list(zip(
        df["_rlnAngleRot"].astype(float),
        df["_rlnAngleTilt"].astype(float),
        df["_rlnAnglePsi"].astype(float),
        df["_rlnOriginXAngst"].astype(float),
        df["_rlnOriginYAngst"].astype(float),
    ))

    size_cache = {}
    for stack in np.unique(stack_names):
        size_cache[stack] = janas_core.sizeMRC(stack)

    ctf_params = None
    if use_ctf:
        c = df.copy()
        if "_rlnPhaseShift" not in c.columns:
            c["_rlnPhaseShift"] = 0.0
        ctf_params = list(zip(
            c["_rlnVoltage"].astype(float),
            c["_rlnDefocusU"].astype(float),
            c["_rlnDefocusV"].astype(float),
            c["_rlnDefocusAngle"].astype(float),
            c["_rlnSphericalAberration"].astype(float),
            c["_rlnCtfBfactor"].astype(float),
            c["_rlnPhaseShift"].astype(float),
            c["_rlnAmplitudeContrast"].astype(float),
            c["_rlnImagePixelSize"].astype(float),
        ))

    labels = {
        "class": df["_rlnClassNumber"].astype(int).to_numpy(),
        "subset": df["_rlnRandomSubset"].astype(int).to_numpy()
    }
    extra = {"optics": optics, "particles_df": df}

    return entries, params, size_cache, pixel_size, ctf_params, labels, extra

class _FSCStore:
    """
    Disk-backed rendezvous for k-space checkpoints:
    each half writes a small .npz with normalized k-space at count 'C'.
    When both halves are present, we compute FSC and append to CSV.
    """
    def __init__(self, root_dir: Path, csv_path: Path, pixel_size: float):
        self.root = Path(root_dir)
        self.root.mkdir(parents=True, exist_ok=True)
        self.csv = Path(csv_path)
        self.pixel_size = float(pixel_size)
        self._lock = threading.Lock()
        if not self.csv.exists():
            with open(self.csv, "w") as fh:
                fh.write("numparticles,FSC05,fsc0143\n")

    def _npz_path(self, cls_id: Optional[int], count: int, half: int) -> Path:
        ctag = "all" if cls_id is None else f"class{int(cls_id)}"
        return self.root / f"kspace_{ctag}_H{int(half)}_{int(count)}.npz"

    def write_half(self, cls_id: Optional[int], count: int, half: int,
                   acc: np.ndarray, weight: np.ndarray):
        # produce normalized k-space S = acc / max(weight, 1)
        w = weight.astype(np.float32)
        S = acc.astype(np.complex64) / (w + 1e-8)
        # zero out unobserved coefficients
        S[(w <= 0)] = 0.0 + 0.0j
        np.savez_compressed(self._npz_path(cls_id, count, half), S=S)

    def try_fuse_and_append(self, cls_id: Optional[int], count: int):
        p1 = self._npz_path(cls_id, count, 1)
        p2 = self._npz_path(cls_id, count, 2)
        if not (p1.exists() and p2.exists()):
            return False
        # single-writer append with lock
        with self._lock:
            # re-check under lock
            if not (p1.exists() and p2.exists()):
                return False
            d1 = np.load(p1)
            d2 = np.load(p2)
            F1 = d1["S"]; F2 = d2["S"]
            freqs, fsc = locres_utils.compute_fsc_from_kspace(F1, F2, pixel_spacing=self.pixel_size)
            # thresholds at 0.5 and 0.143, in 1/Å
            def _cross_at(thr: float) -> float:
                above = fsc >= thr
                idx = np.where(above[:-1] & (~above[1:]))[0]
                if idx.size == 0:
                    return float("nan")
                i = int(idx[0])
                f1, f2 = float(freqs[i]), float(freqs[i+1])
                y1, y2 = float(fsc[i]), float(fsc[i+1])
                frac = (y1 - thr) / (y1 - y2 + 1e-12)
                return f1 + frac * (f2 - f1)
            fc05 = _cross_at(0.5)
            fc143 = _cross_at(0.143)
            with open(self.csv, "a") as fh:
                fh.write(f"{int(count)},{fc05:.6g},{fc143:.6g}\n")
        return True


# ----------------------- (Optional) Torch -----------------------
try:
    import torch
    import torch.nn.functional as F
    TORCH_AVAILABLE = True
except Exception:
    TORCH_AVAILABLE = False


# ======================= CPU multiprocessing core =======================

_G = {}  # fork-visible globals

def _init_globals(gdict):
    _G.update(gdict)

def _worker_reconstruct(b0: int, b1: int, procnum: int, tmp_dir: str, progress):
    N          = _G["N"]
    halfN      = N // 2
    pixel_size = _G["pixel_size"]
    entries    = _G["entries"]
    params     = _G["params"]
    use_ctf    = _G["use_ctf"]
    ctf_params = _G["ctf_params"]

    freq_unit = 1.0 / (N * pixel_size)
    two_pi = 2.0 * np.pi
    int_freq = np.arange(N, dtype=np.int32)
    int_freq[int_freq > halfN] -= N
    KX_cpu, KY_cpu = np.meshgrid(int_freq, int_freq, indexing="xy")
    FX_cpu = KX_cpu * freq_unit
    FY_cpu = KY_cpu * freq_unit

    acc    = np.zeros((N, N, N), dtype=np.__dict__['complex64'])
    weight = np.zeros((N, N, N), dtype=np.float32)

    for i in range(b0, b1):
        stack, sl = entries[i]
        rot, tilt, psi, ox, oy = params[i]
        Rmat = R.from_euler("ZYZ", [rot, tilt, psi], degrees=True).as_matrix()

        raw = np.asarray(janas_core.ReadMrcSlice(stack, sl), dtype=np.float32).reshape(N, N)

        if use_ctf:
            raw = np.array(
                assessParticles.transformCtfImage(
                    raw.flatten().tolist(),
                    N, N, pixel_size,
                    *ctf_params[i]
                ),
                dtype=np.float32
            ).reshape(N, N)

        raw -= raw.mean()
        F2 = np.fft.fft2(np.fft.ifftshift(raw))

        dx, dy = ox / pixel_size, oy / pixel_size
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

        with progress.get_lock():
            progress.value += 1

    acc_path = os.path.join(tmp_dir, f"acc_{procnum:04d}.npy")
    w_path   = os.path.join(tmp_dir, f"w_{procnum:04d}.npy")
    np.save(acc_path, acc, allow_pickle=False)
    np.save(w_path, weight, allow_pickle=False)
    return (b0, b1, acc_path, w_path)

# ======================= Helpers =======================



def _ensure_suffix(parts: List[str], include_subset: bool, sub: Optional[int]) -> List[str]:
    """Append rec/recH# suffix for consistent naming."""
    if include_subset and (sub is not None):
        parts.append(f"recH{sub}")
    else:
        parts.append("rec")
    return parts

def write_checkpoint_star(out_dir: Path, basename: Optional[str],
                          count: int, optics_df, particles_df_sorted_filtered):
    out_dir.mkdir(parents=True, exist_ok=True)
    base = basename if basename else "reconstruction"
    out_path = out_dir / f"{base}{count}.star"
    drop_cols = ["_rlnImagePixelSize", "_rlnVoltage",
                 "_rlnAmplitudeContrast", "_rlnSphericalAberration"]
    particles_df_sorted_filtered = particles_df_sorted_filtered.drop(columns=drop_cols, errors="ignore")
    write_relion31_star(str(out_path), optics_df, particles_df_sorted_filtered)
    print(f"[checkpoint] wrote {out_path}")








# ======================= Selection & helpers =======================

def filter_indices_by_class_subset(labels: Dict[str, np.ndarray],
                                   classes: Optional[List[int]],
                                   subsets: Optional[List[int]]) -> np.ndarray:
    n = labels["class"].shape[0]
    mask = np.ones(n, dtype=bool)
    if classes is not None:
        mask &= np.isin(labels["class"], np.asarray(classes, dtype=int))
    if subsets is not None:
        mask &= np.isin(labels["subset"], np.asarray(subsets, dtype=int))
    return mask

def select_available(labels: Dict[str, np.ndarray]) -> Tuple[List[int], List[int]]:
    classes = sorted(list({int(x) for x in labels["class"].tolist()}))
    subsets = sorted(list({int(x) for x in labels["subset"].tolist()}))
    return classes, subsets

def build_output_name(base_out: Optional[str],
                      basename: Optional[str],
                      star_path: str,
                      include_class: bool,
                      include_subset: bool,
                      cls: Optional[int],
                      sub: Optional[int],
                      single_recon: bool,
                      bulk_dir: Optional[Path] = None) -> str:
    """
    Ensure names end in _rec.mrc (no subsets) or _recH#.mrc (with subsets).
    If bulk_dir is provided, place outputs there; otherwise use CWD.
    """
    target_dir = Path(bulk_dir) if bulk_dir else Path(".")
    if single_recon:
        if base_out:
            return str(Path(base_out))
        stem = basename if basename else Path(star_path).stem
        fname = f"{stem}_rec.mrc"
        return str(target_dir / fname)

    parts = []
    if basename:
        parts.append(basename)
    if include_class and (cls is not None):
        parts.append(f"class_{cls}")
    parts = _ensure_suffix(parts, include_subset, sub)
    fname = "_".join(parts) + ".mrc"
    return str(target_dir / fname)

# ======================= Single job runner (top-level for pickling) =======================
def run_single_reconstruction(job: Dict[str, Any]) -> Tuple[str, bool, str, Optional[int], Optional[int]]:
    """
    Executes one (class, subset) reconstruction job.
    Returns (out_path, success, message, class_id, subset_id)
    NOTE: out_path MUST be the actual file that was written, including any subrec count tag.
    """
    star_file      = job["star_file"]
    device         = job["device"]
    gpu_index      = job.get("gpu_index", None)
    cpu_workers    = job.get("cpu_workers", 1)
    gpu_batch      = job.get("gpu_batch", 20)
    use_ctf        = job["use_ctf"]
    pixel_size_opt = job["pixel_size"]
    base_out       = job["out_path"]
    include_class  = job["include_class"]
    include_subset = job["include_subset"]
    cls_filter     = job["cls_filter"]
    sub_filter     = job["sub_filter"]
    class_id       = job.get("class_id", None)
    subset_id      = job.get("subset_id", None)
    sort_label     = job.get("sort_label", None)
    sort_ascending = job.get("sort_ascending", True)
    subrec_only    = bool(job.get("subrec_only", False))


    # --- FSC rendezvous (only when requested and we are reconstructing a half) ---
    fsc_enabled = bool(job.get("fsc_enabled", False))
    fsc_csv = job.get("fsc_csv_path", None)
    base_dir_for_fsc = job.get("basename_dir", None)
    cls_here = job.get("class_id", None)
    sub_here = job.get("subset_id", None)

    fsc_store = None  # constructed once we know the effective pixel size

    # LPF controls
    lp_start       = float(job.get("lp_start", 0.9))
    lp_rolloff     = job.get("lp_rolloff", None)
    lp_rolloff     = None if lp_rolloff in (None, "None") else int(lp_rolloff)
    lp_outer_scale = float(job.get("lp_outer_scale", 1.1))
    lp_add_px      = job.get("lp_add_px", None)
    lp_add_px      = None if lp_add_px in (None, "None") else int(lp_add_px)

    # Progressive counts (per-job, already split by main for halves/classes)
    subrec_counts  = job.get("subrec_counts", []) or []
    checkpoint_labels = job.get("checkpoint_labels", None)

    # Parse STAR and selection
    entries_all, params_all, size_cache, optics_pixel, ctf_params_all, labels, _ = parse_star_with_ctf_and_labels(
        star_file, use_ctf=use_ctf, sort_label=sort_label, sort_ascending=bool(sort_ascending)
    )
    pixel_size = pixel_size_opt if (pixel_size_opt is not None) else optics_pixel

    mask = filter_indices_by_class_subset(labels, cls_filter, sub_filter)
    if not np.any(mask):
        return base_out, False, "No particles match the selection.", class_id, subset_id

    entries = [entries_all[i] for i, m in enumerate(mask) if m]
    params  = [params_all[i]  for i, m in enumerate(mask) if m]
    ctf_params = [ctf_params_all[i] for i, m in enumerate(mask) if m] if (use_ctf and ctf_params_all is not None) else None
    if subrec_only:
        sub_counts = job.get("subrec_counts", []) or []
        if len(sub_counts) > 0:
            max_count = max(int(c) for c in sub_counts if c > 0)
            if max_count < len(entries):
                entries = entries[:max_count]
                params = params[:max_count]
                if ctf_params is not None:
                    ctf_params = ctf_params[:max_count]

    first_stack = entries[0][0]
    if not os.path.exists(first_stack):
        sys.exit(
            f"ERROR: image stack '{first_stack}' not found on disk "
            f"(first stack referenced in STAR file '{star_file}')."
        )
    try:
        box_size = int(janas_core.sizeMRC(first_stack)[0])
    except Exception as e:
        sys.exit(
            f"ERROR: failed to read MRC header from '{first_stack}': {e}"
        )

    # Now that pixel_size is known, instantiate the store if needed
    if fsc_enabled and (sub_here in (1, 2)) and (fsc_csv is not None):
        kspace_root = (Path(base_dir_for_fsc) / "_kspace_checkpoints") if base_dir_for_fsc else Path(fsc_csv).parent / "_kspace_checkpoints"
        # pixel_size is in Å/pixel; FSC freqs will be in 1/Å
        fsc_store = _FSCStore(kspace_root, Path(fsc_csv), pixel_size)

    # Hook to capture k-space before iFFT at each checkpoint
    def _on_checkpoint_kspace(acc_ck: np.ndarray, weight_ck: np.ndarray, count_label: int):
        if fsc_store is None or (sub_here not in (1, 2)):
            return
        try:
            # Write this half's k-space and try to fuse with its mate
            fsc_store.write_half(cls_here, int(count_label), int(sub_here), acc_ck, weight_ck)
            fsc_store.try_fuse_and_append(cls_here, int(count_label))
        except Exception as e:
            print(f"[FSC] WARNING at count {count_label}: {e}")

    acc, weight = reconstructor_utils.reconstruct_3D(
        entries, params, pixel_size, box_size,
        acc=None, weight=None,
        use_ctf=use_ctf, ctf_params=ctf_params,
        ctf_mode=job.get("ctf_mode", None),
        wiener_tau=float(job.get("wiener_tau", 0.1)),
        device=("cuda" if device == "gpu" else "cpu"),
        cpu_workers=cpu_workers,
        gpu_index=gpu_index,
        gpu_batch=gpu_batch,
        checkpoint_counts=subrec_counts,
        checkpoint_out_path=base_out,
        lp_start=lp_start,
        lp_rolloff=lp_rolloff,
        lp_outer_scale=lp_outer_scale,
        lp_add_px=lp_add_px,
        checkpoint_labels=checkpoint_labels,  # global labels for filenames
        on_checkpoint_kspace=_on_checkpoint_kspace,
    )

    vol = reconstructor_utils.post_process(
        acc, weight,
        pixel_size=pixel_size,
        device=("cuda" if device == "gpu" else "cpu"),
        lp_start=lp_start,
        lp_rolloff=lp_rolloff,
        lp_outer_scale=lp_outer_scale,
        lp_add_px=lp_add_px,
    )
    # Optional Gaussian blur (σ in Å) on the same device used for reconstruction
    post_blur_sigmaA = job.get("post_blur_sigmaA", None)
    if post_blur_sigmaA not in (None, "None"):
        try:
            dev = "cuda" if device == "gpu" else "cpu"
            vol = janas_mapProcess.gaussian_blur_array(
                vol, angpix=pixel_size, sigmaA=float(post_blur_sigmaA),
                device=dev, gpu_index=gpu_index
            )
        except Exception as e:
            print(f"[post-blur] WARNING: failed to blur with σ={post_blur_sigmaA}: {e}")

    # Only append the total particle count to the FINAL filename
    # if a subrec workflow was requested (explicit --subrec or --regular-subrec).
    from pathlib import Path

    final_out = base_out
    if len(subrec_counts) > 0:
        global_label = job.get("final_count_label", None)
        count_label = int(global_label) if (global_label is not None) else len(entries)

        out_dir = Path(base_out).parent

        basename_root = job.get("basename_root")
        if not basename_root:
            basename_root = Path(base_out).stem

        suffix_parts = []

        if class_id is not None:
            suffix_parts.append(f"C{int(class_id)}")

        if subset_id is not None:
            suffix_parts.append(f"recH{int(subset_id)}")
        elif class_id is not None:
            suffix_parts.append("rec")
        else:
            suffix_parts.append("rec")

        suffix_str = "_".join(suffix_parts)

        filename = f"{basename_root}{count_label}_{suffix_str}.mrc"
        final_out = str(out_dir / filename)

    IO_utils.write_mrc(final_out, vol, pixel_size)

    # Return the ACTUAL filename written so downstream (locres) can find it.
    return final_out, True, f"Wrote {final_out}", class_id, subset_id




# For multiprocessing Pool on CPU bulk
def cpu_job_worker(jobdict: Dict[str, Any]) -> Tuple[str, bool, str, Optional[int], Optional[int]]:
    try:
        return run_single_reconstruction(jobdict)
    except Exception as e:
        return jobdict["out_path"], False, f"ERROR: {e}", jobdict.get("class_id", None), jobdict.get("subset_id", None)

# ======================= CLI =======================

def parse_args():
    p = argparse.ArgumentParser(
        description="3D backprojection from JANAS slices with class/subset selection, multi-GPU/CPU orchestration, and progressive checkpoints with STAR exports."
    )
    p.add_argument("star_file", help="particles.star input")

    # Reconstruction controls
    p.add_argument("--pixel_size", type=float, default=None,
                   help="Override pixel size in Å (otherwise use _rlnImagePixelSize from optics)")

    p.add_argument(
        "--ctf",
        nargs="?",
        const="wiener",
        choices=["modulate", "phaseflip", "wiener"],
        help="CTF correction mode. No argument == 'modulate'. Choices: modulate|phaseflip|wiener."
    )
    p.add_argument(
        "--wiener-tau",
        type=float,
        default=0.1,
        help="Regularization τ for --ctf wiener (default: 0.1)."
    )
    p.add_argument(
        "--wiener-tau-auto",
        type=str,
        default=None,
        help="CSV with columns: Frequency(1/Å),SSNR. Used to build a per-shell τ(k)=1/SSNR(k) image (no sharpening)."
    )

    # Selection controls
    p.add_argument("--class", dest="classes", nargs="+", default=None,
                   help="Class numbers or 'all'. Example: --class 1 2 5  or  --class all")
    p.add_argument("--subset", dest="subsets", nargs="+", default=None,
                   help="Subset numbers or 'all'. Example: --subset 1 2  or  --subset all")
    p.add_argument("--all-recs", nargs="?", const="reconstructions", default=None,
                   help="Bulk: equivalent to --class all --subset all. Optional DIR argument (default: 'reconstructions').")

    # Progressive checkpoints
    p.add_argument("--subrec", dest="subrec", nargs="+", default=None,
                   help="Progressive checkpoints by particle count (e.g., --subrec 3000,6000). Sorted ascending; writes outputs and STAR at each target (STAR uses --basename).")
    
    p.add_argument("--subrec-only", dest="subrec_only", nargs="+", default=None,
                   help="Like --subrec, but only reconstruct up to the largest checkpoint (no full-particle reconstruction). Example: --subrec-only 3000,6000.")

    p.add_argument("--regular-subrec", nargs="+", default=None,
                   help="Regularly spaced checkpoints: --regular-subrec <count> [min] [max]. min/max can be numbers or 'min'/'max'. Defaults: min=0, max=<selected N>.")

    # Sorting
    p.add_argument("--sort", nargs="+", default=None,
                   help="Sort particles by STAR column before backprojection. Usage: --sort <label> [ascending|descending]. Default ascending.")

    # Output controls
    p.add_argument("--out-basename", "--basename", dest="basename", type=str, default=None,
                   help="Basename or DIR/BASENAME for outputs and checkpoint STARs. Example: --out-basename results/EMPIAR_10667_run")
    
    p.add_argument("--o", dest="single_output", default=None,
                   help="Output filename (used only for single reconstruction).")

    # Device selection
    p.add_argument("--gpu", nargs="*", default=None,
                   help="Use GPU(s). Example: --gpu (==[0]), or --gpu 0, or --gpu 0 1")
    p.add_argument("--cpu", nargs="?", const="1", default=None,
                   help="CPU multiprocessing. Single recon: worker count. Multiple recons: number of concurrent reconstructions (each runs single-process). Example: --cpu 30")

    # GPU batch
    p.add_argument("--gpu-batch", type=int, default=20,
                   help="GPU mini-batch size per iteration (controls memory/speed trade-off).")

    # Low-pass mask controls
    p.add_argument("--lp-start", type=float, default=0.9,
                   help="Fraction of R_out where attenuation begins (0<start<=1). Default 0.9.")
    p.add_argument("--lp-rolloff", type=int, default=None,
                   help="Override: raised-cosine roll-off width (in Fourier pixels) up to R_out. If set, overrides --lp-start.")
    p.add_argument("--lp-outer-scale", type=float, default=1.1,
                   help="Scale R_out relative to per-axis Nyquist (N/2). Range (0, sqrt(3)]. E.g., 1.1 keeps up to 1.1×Nyquist radially.")
    p.add_argument("--lp-add-px", type=int, default=None,
                   help="Additive pixels to extend R_out beyond lp-outer-scale*(N/2). For example, use N/8 to extend by 1/8 of box.")

    # ---- Local-resolution options (minimal surface) ----
    p.add_argument("--locres", action="store_true",
                   help="After reconstructing half-maps H1 and H2 (requires --subset 1 2), compute local-resolution maps.")
    p.add_argument("--locres-mask", type=str, default=None,
                   help="Optional mask (.mrc) to restrict local FSC centres (used only if --locres).")
    p.add_argument("--locres-cpu", type=int, default=1,
                   help="Number of CPU processes for local-resolution computation (used only if --locres).")
    p.add_argument("--locres-stats", action="store_true", help="After computing locres at checkpoints, aggregate per-voxel best values and write stats_locres.mrc, stats_minParticles.mrc, stats_intensity_recH1.mrc, stats_intensity_recH2.mrc, and bestRanked_locres_values.csv.")
    p.add_argument("--postprocessBlur", type=float, default=0.02,
               help="Gaussian blur σ in Å applied after reconstruction to each written map. "
                    "If reconstruction uses GPU, blurring runs on the same CUDA device.")



    return p.parse_args()




# ======================= Main orchestration =======================

def _fmt_hms(sec: float) -> str:
    return time.strftime('%H:%M:%S', time.gmtime(max(0, sec)))

def _split_basename_path(basename_opt: Optional[str]) -> Tuple[Optional[Path], Optional[str]]:
    """
    Split a possibly path-like basename into (directory, stem).
    If no directory is present, returns (None, basename).
    """
    if not basename_opt:
        return None, None
    p = Path(basename_opt).expanduser()
    dir_part = p.parent if str(p.parent) not in ("", ".") else None

    known_exts = {".mrc", ".map"}
    if p.suffix.lower() in known_exts:
        stem = p.stem
    else:
        stem = p.name

    return dir_part, stem



def main():
    args = parse_args()

    # Parse --sort
    sort_label = None
    sort_ascending = True
    if args.sort is not None:
        if len(args.sort) < 1:
            sys.exit("ERROR: --sort requires at least a label.")
        sort_label = args.sort[0]
        if len(args.sort) >= 2:
            order = args.sort[1].strip().lower()
            if order not in ("ascending", "descending"):
                sys.exit("ERROR: --sort order must be 'ascending' or 'descending'.")
            sort_ascending = (order == "ascending")
        else:
            sort_ascending = True  # default ascending

    # Handle basename directory redirection
    basename_dir, basename_stem = _split_basename_path(args.basename)
    if basename_dir:
        basename_dir.mkdir(parents=True, exist_ok=True)
    # Overwrite args.basename to just the stem (used in filenames)
    if basename_stem:
        args.basename = basename_stem

    # Initial parse to enumerate classes/subsets (unsorted) and to carry full particles df for STAR export later
    _, _, _, _, _, labels0, extra0 = parse_star_with_ctf_and_labels(args.star_file, use_ctf=args.ctf)
    avail_classes, avail_subsets = select_available(labels0)

    # Interpret --all-recs DIR
    raw_bulk_dir: Optional[Path] = None
    want_all = args.all_recs is not None
    if want_all:
        raw_bulk_dir = Path(args.all_recs).expanduser().resolve()
        raw_bulk_dir.mkdir(parents=True, exist_ok=True)

    # Choose preferred output directory: basename_dir has priority over --all-recs dir
    preferred_dir: Optional[Path] = basename_dir if basename_dir else raw_bulk_dir

    # Decode --class
    classes_sel = None
    include_class_in_name = False
    if want_all or (args.classes is not None):
        include_class_in_name = True
        if want_all or (len(args.classes) == 1 and str(args.classes[0]).lower() == "all"):
            classes_sel = avail_classes
        else:
            tokens = []
            for tok in args.classes:
                tokens += [t for t in str(tok).split(",") if t]
            classes_sel = sorted({int(x) for x in tokens})

    # Decode --subset (half-map IDs only, keep original behavior)
    subsets_sel = None
    include_subset_in_name = False
    if want_all or (args.subsets is not None):
        include_subset_in_name = True
        if want_all or (len(args.subsets) == 1 and str(args.subsets[0]).lower() == "all"):
            subsets_sel = avail_subsets
        else:
            subsets_sel = sorted({int(x) for tok in args.subsets for x in str(tok).split(",") if x})

    # Decode --subrec (explicit progressive checkpoints)
    # Decode --subrec / --subrec-only (explicit progressive checkpoints)
    subrec_counts: List[int] = []
    subrec_only_mode = False

    def _parse_subrec_list(arg_list, flag_name: str) -> List[int]:
        tokens = []
        for tok in arg_list:
            tokens += [t.strip() for t in str(tok).split(",") if t.strip()]
        vals: List[int] = []
        for t in tokens:
            try:
                vals.append(int(t))
            except Exception:
                sys.exit(f"ERROR: --{flag_name} token '{t}' is not an integer.")
        return vals

    # Do not allow both at once
    if (args.subrec is not None) and (args.subrec_only is not None):
        sys.exit("ERROR: --subrec and --subrec-only cannot be used together.")

    if args.subrec is not None:
        subrec_counts = _parse_subrec_list(args.subrec, "subrec")
    elif args.subrec_only is not None:
        subrec_counts = _parse_subrec_list(args.subrec_only, "subrec-only")
        subrec_only_mode = True

    subrec_counts = sorted(set([c for c in subrec_counts if c > 0]))

    # Parse --regular-subrec (defer numeric min/max until we know selection size)
    reg_spec = None
    if args.regular_subrec is not None:
        toks = args.regular_subrec
        if len(toks) < 1:
            sys.exit("ERROR: --regular-subrec requires at least <count>.")
        try:
            reg_count = int(str(toks[0]))
            if reg_count <= 0:
                raise ValueError
        except Exception:
            sys.exit(f"ERROR: --regular-subrec count '{toks[0]}' must be a positive integer.")
        reg_min = str(toks[1]) if len(toks) >= 2 else "min"
        reg_max = str(toks[2]) if len(toks) >= 3 else "max"
        reg_spec = (reg_count, reg_min, reg_max)

    # Build job list with explicit metadata (class_id, subset_id)
    jobs: List[Dict[str, Any]] = []
    effective_dir = preferred_dir  # pass this into build_output_name (treated as bulk_dir)
    if (classes_sel is None) and (subsets_sel is None):
        # Single recon
        out_name = build_output_name(
            base_out=args.single_output,
            basename=args.basename,
            star_path=args.star_file,
            include_class=False,
            include_subset=False,
            cls=None, sub=None,
            single_recon=True,
            bulk_dir=effective_dir
        )
        jobs.append(dict(
            star_file=args.star_file,
            device=("gpu" if args.gpu is not None else "cpu"),
            gpu_index=None,
            cpu_workers=(int(args.cpu) if args.cpu is not None else 1),
            gpu_batch=int(args.gpu_batch),
            use_ctf=bool(args.ctf),
            ctf_mode=(args.ctf if args.ctf else None),
            wiener_tau=float(args.wiener_tau),
            wiener_tau_csv=args.wiener_tau_auto,
            pixel_size=args.pixel_size,
            out_path=str(out_name),
            include_class=False,
            include_subset=False,
            cls_filter=None,
            sub_filter=None,
            class_id=None,
            subset_id=None,
            lp_start=float(args.lp_start),
            lp_rolloff=(None if args.lp_rolloff is None else int(args.lp_rolloff)),
            lp_outer_scale=float(args.lp_outer_scale),
            lp_add_px=(None if args.lp_add_px is None else int(args.lp_add_px)),
            subrec_counts=subrec_counts,           # may be overwritten below by split plan
            checkpoint_labels=subrec_counts,       # labels for filenames
            sort_label=sort_label,
            sort_ascending=bool(sort_ascending),
            post_blur_sigmaA=(None if args.postprocessBlur is None else float(args.postprocessBlur)),
            subrec_only=subrec_only_mode,
            basename_root=(args.basename if args.basename else Path(args.star_file).stem),
        ))
    else:
        cls_list = classes_sel if classes_sel is not None else [None]
        sub_list = subsets_sel if subsets_sel is not None else [None]
        for cls in cls_list:
            for sub in sub_list:
                out_name = build_output_name(
                    base_out=None,
                    basename=args.basename,
                    star_path=args.star_file,
                    include_class=include_class_in_name,
                    include_subset=include_subset_in_name,
                    cls=cls, sub=sub,
                    single_recon=False,
                    bulk_dir=effective_dir
                )
                jobs.append(dict(
                    star_file=args.star_file,
                    device=("gpu" if args.gpu is not None else "cpu"),
                    gpu_index=None,
                    cpu_workers=1,   # single-process per job in bulk
                    gpu_batch=int(args.gpu_batch),
                    use_ctf=bool(args.ctf),
                    ctf_mode=(args.ctf if args.ctf else None),
                    wiener_tau=float(args.wiener_tau),
                    wiener_tau_csv=args.wiener_tau_auto,
                    pixel_size=args.pixel_size,
                    out_path=str(out_name),
                    include_class=include_class_in_name,
                    include_subset=include_subset_in_name,
                    cls_filter=(None if cls is None else [int(cls)]),
                    sub_filter=(None if sub is None else [int(sub)]),
                    class_id=(None if cls is None else int(cls)),
                    subset_id=(None if sub is None else int(sub)),
                    lp_start=float(args.lp_start),
                    lp_rolloff=(None if args.lp_rolloff is None else int(args.lp_rolloff)),
                    lp_outer_scale=float(args.lp_outer_scale),
                    lp_add_px=(None if args.lp_add_px is None else int(args.lp_add_px)),
                    subrec_counts=subrec_counts,     # will be split below
                    checkpoint_labels=subrec_counts, # will remain global
                    sort_label=sort_label,
                    sort_ascending=bool(sort_ascending),
                    post_blur_sigmaA=(None if args.postprocessBlur is None else float(args.postprocessBlur)),
                    subrec_only=subrec_only_mode,
                    basename_root=(args.basename if args.basename else Path(args.star_file).stem),
                ))


    # ---------- STAR checkpoint planning & per-job split (for halves/classes) ----------
    want_checkpoints = (len(subrec_counts) > 0) or (reg_spec is not None)
    if want_checkpoints:
        # Parse sorted once to produce the selection and split per job
        _, _, _, _, _, labels_sorted, extra_sorted = parse_star_with_ctf_and_labels(
            args.star_file, use_ctf=args.ctf, sort_label=sort_label, sort_ascending=sort_ascending
        )
        df_sorted = extra_sorted["particles_df"].copy()

        # Build master selection mask across requested classes/subsets
        master_mask = filter_indices_by_class_subset(
            {"class": df_sorted["_rlnClassNumber"].astype(int).to_numpy(),
             "subset": df_sorted["_rlnRandomSubset"].astype(int).to_numpy()},
            classes_sel, subsets_sel
        )
        master_df = df_sorted.loc[master_mask].reset_index(drop=True)
        total_sel = len(master_df)
        if subrec_only_mode and subrec_counts:
            global_final = subrec_counts[-1]
        else:
            global_final = int(total_sel)

        for j in jobs:
            j["final_count_label"] = int(global_final)


        # Compute regular checkpoints if requested
        if reg_spec is not None:
            rcount, rmin_tok, rmax_tok = reg_spec
            def _tok_to_int(tok: str, default_minmax: int) -> int:
                t = tok.strip().lower()
                if t in ("min", "start"):
                    return 0
                if t in ("max", "end"):
                    return default_minmax
                try:
                    v = int(float(tok))
                except Exception:
                    sys.exit(f"ERROR: --regular-subrec bound '{tok}' must be an integer or 'min'/'max'.")
                return v
            rmin = _tok_to_int(rmin_tok, 0)
            rmax = _tok_to_int(rmax_tok, total_sel)
            if rmin < 0: rmin = 0
            if rmax < 0: rmax = 0
            if rmax < rmin:
                rmin, rmax = rmax, rmin
            # linspace inclusive between rmin and rmax
            regs = np.linspace(rmin, rmax, num=rcount, endpoint=True, dtype=float)
            regs = np.rint(regs).astype(int).tolist()
            # Remove 0 and duplicates, clamp to [1, total_sel]
            regs = sorted({c for c in regs if 1 <= c <= total_sel})
            # Union with explicit subrec_counts
            subrec_counts = sorted(set(list(subrec_counts) + regs))

        if not subrec_counts:
            pass
        else:
            # Directory for checkpoint STARs: prefer preferred_dir, else parent of first job's output
            if preferred_dir:
                out_dir = preferred_dir
            else:
                out_dir = Path(jobs[0])["out_path"].parent  # type: ignore[index]
            base_for_star = args.basename if args.basename else Path(args.star_file).stem

            # Prepare lookup columns
            SUBCOL = "_rlnRandomSubset"
            CLASSCOL = "_rlnClassNumber"

            # Precompute job masks (on master_df) and assign per-job count list
            job_masks = []
            for j in jobs:
                jmask = np.ones(len(master_df), dtype=bool)
                if j["cls_filter"] is not None:
                    jmask &= master_df[CLASSCOL].astype(int).isin(j["cls_filter"]).to_numpy()
                if j["sub_filter"] is not None:
                    jmask &= master_df[SUBCOL].astype(int).isin(j["sub_filter"]).to_numpy()
                job_masks.append(jmask)

            # Fill per-job split counts for each checkpoint and write STARs
            per_job_counts = [ [] for _ in jobs ]
            optics_df_for_star = extra0["optics"]  # optics unchanged
            csv_rows = [ "num_particles,basename" ]


            for C in subrec_counts:
                topN = min(C, len(master_df))
                top_df = master_df.iloc[:topN].copy()
                write_checkpoint_star(out_dir, base_for_star, C, optics_df_for_star, top_df)
                try:
                    rel_dir = os.path.relpath(out_dir, Path.cwd())
                except Exception:
                    rel_dir = str(out_dir)
                basename_with_dir = (Path(rel_dir) / f"{base_for_star}_{topN}").as_posix()
                csv_rows.append(f"{topN},{basename_with_dir}")
                    
                for ji, jmask in enumerate(job_masks):
                    count_j = int(np.count_nonzero(jmask[:topN]))
                    per_job_counts[ji].append(count_j)

            # Store per-job counts back into job dicts and keep global labels
            for ji, j in enumerate(jobs):
                j["subrec_counts"] = per_job_counts[ji]
                j["checkpoint_labels"] = subrec_counts
            csv_path = out_dir / "info_checkpoints_recs.csv"
            with open(csv_path, "w") as fh:
                fh.write("\n".join(csv_rows) + "\n")
            print(f"[checkpoint] wrote {csv_path}")

    # Storage for results with metadata
    results: List[Tuple[str, bool, str, Optional[int], Optional[int]]] = []

    # -------- GPU orchestration --------
    if args.gpu is not None:
        if not TORCH_AVAILABLE:
            print(
                "\nERROR: PyTorch is required for GPU reconstruction but is not installed.\n"
                "Install it before running with --gpu. Pick the wheel matching your CUDA\n"
                "driver (check with `nvidia-smi`), for example:\n\n"
                "    python -m pip install torch --index-url https://download.pytorch.org/whl/cu128\n"
                "    python -m pip install torch --index-url https://download.pytorch.org/whl/cu121\n"
                "    python -m pip install torch --index-url https://download.pytorch.org/whl/cu118\n\n"
                "Or, for a generic install from the default PyPI index:\n\n"
                "    pip install 'janas[gpu]'   # or:   pip install torch\n\n"
                "See docs/installation.md (section 'GPU support (PyTorch)') for details.\n",
                file=sys.stderr,
            )
            sys.exit(1)
        # Resolve GPU list
        if len(args.gpu) == 0:
            gpu_list = [0]
        else:
            gpu_list = [int(x) for x in args.gpu]

        # Decide grouping dimension:
        group_by_subset = (subsets_sel is not None) or (args.subsets is not None)

        if group_by_subset:
            def key_fn(j):  # subset id or None
                return j.get("subset_id", None)
        else:
            def key_fn(j):  # class id or None
                return j.get("class_id", None)

        jobs_by_key: Dict[Optional[int], List[Dict[str, Any]]] = {}
        for j in jobs:
            jobs_by_key.setdefault(key_fn(j), []).append(j)

        def _sort_key(x):
            return (-1 if x is None else x)
        keys_sorted = sorted(jobs_by_key.keys(), key=_sort_key)

        # Assign each key-group wholesale to a GPU in round-robin order.
        gpu_queues: Dict[int, List[Dict[str, Any]]] = {g: [] for g in gpu_list}
        for i, k in enumerate(keys_sorted):
            g = gpu_list[i % len(gpu_list)]
            for j in jobs_by_key[k]:
                j["gpu_index"] = g
                gpu_queues[g].append(j)

        import multiprocessing as mp
        try:
            mp.set_start_method("fork")
        except RuntimeError:
            pass
        manager = mp.Manager()
        results_q = manager.Queue()

        def _gpu_worker(gid: int, queue_jobs: List[Dict[str, Any]], results_q):
            # up to two jobs in parallel per GPU
            import threading
            from queue import Queue, Empty
            MAX_PAR = 2
            q = Queue()
            for j in queue_jobs:
                q.put(j)

            def describe(job: Dict[str, Any]) -> Tuple[str,str]:
                c = job.get("class_id", None); s = job.get("subset_id", None)
                cs = f"class {c}" if c is not None else "all classes"
                ss = f"H{s}" if s is not None else "all subsets"
                return cs, ss

            def runner():
                while True:
                    try:
                        job = q.get_nowait()
                    except Empty:
                        return
                    cs, ss = describe(job)
                    start = time.time()
                    print(f"[GPU{gid}] START {cs}, subset {ss} — out: {Path(job['out_path']).name}")
                    try:
                        outp, ok, msg, cid, sid = run_single_reconstruction(job)
                    except Exception as e:
                        outp, ok, msg, cid, sid = job["out_path"], False, f"ERROR: {e}", job.get("class_id", None), job.get("subset_id", None)
                    elapsed = time.time() - start
                    tag = "DONE" if ok else "FAIL"
                    print(f"[GPU{gid}] {tag}  {cs}, subset {ss} — elapsed {_fmt_hms(elapsed)} — {msg}")
                    results_q.put((outp, ok, msg, cid, sid))
                    q.task_done()

            threads = [threading.Thread(target=runner, daemon=True) for _ in range(min(MAX_PAR, len(queue_jobs)))]
            for t in threads: t.start()
            for t in threads: t.join()

        procs = []
        for gid, qjobs in gpu_queues.items():
            if not qjobs:
                continue
            p = mp.Process(target=_gpu_worker, args=(gid, qjobs, results_q))
            p.start()
            procs.append(p)

        alive = True
        while alive:
            alive = any(p.is_alive() for p in procs)
            while not results_q.empty():
                results.append(results_q.get())
            time.sleep(0.2)
        for p in procs: p.join()
        while not results_q.empty():
            results.append(results_q.get())

    # -------- CPU orchestration (bulk or single) --------
    else:
        cpu_flag = args.cpu is not None
        if len(jobs) == 1:
            job = jobs[0]
            job["cpu_workers"] = int(args.cpu) if cpu_flag else 1
            outp, ok, msg, cid, sid = run_single_reconstruction(job)
            results.append((outp, ok, msg, cid, sid))
            print(msg if ok else f"FAILED: {msg}")
        else:
            max_concurrent = int(args.cpu) if cpu_flag else 1
            if max_concurrent < 1: max_concurrent = 1
            import multiprocessing as mp
            try:
                mp.set_start_method("fork")
            except RuntimeError:
                pass
            with mp.Pool(processes=max_concurrent) as pool:
                for outp, ok, msg, cid, sid in pool.imap_unordered(cpu_job_worker, jobs, chunksize=1):
                    results.append((outp, ok, msg, cid, sid))
                    print(msg if ok else f"FAILED: {msg}")

    # -------- Per-class STARs (only when bulk/all-recs) --------
    want_all = args.all_recs is not None
    if want_all:
        optics_df = extra0["optics"]
        particles_df_full = extra0["particles_df"]
        class_list = classes_sel if classes_sel is not None else avail_classes
        # Choose directory for these STARs: prefer preferred_dir
        star_out_dir = preferred_dir if preferred_dir else Path(args.all_recs).expanduser().resolve()
        star_out_dir.mkdir(parents=True, exist_ok=True)
        for cls in class_list:
            part_sub = particles_df_full[particles_df_full["_rlnClassNumber"].astype(int) == int(cls)].copy()
            if part_sub.empty:
                continue
            star_name = (args.basename + f"_class_{cls}.star") if args.basename else f"class_{cls}.star"
            star_path = star_out_dir / star_name
            drop_cols = ["_rlnImagePixelSize", "_rlnVoltage","_rlnAmplitudeContrast", "_rlnSphericalAberration"]
            part_sub = part_sub.drop(columns=drop_cols, errors="ignore")
            write_relion31_star(str(star_path), optics_df, part_sub)
            print(f"Wrote {star_path}")

    # -------- Averages H1/H2 -> rec.mrc (only when bulk/all-recs) --------
    if want_all:
        # Build class -> {subset_id: path}
        by_class: Dict[int, Dict[int, str]] = {}
        for outp, ok, _, cid, sid in results:
            if not ok:
                continue
            if (cid is None) or (sid is None):
                continue
            by_class.setdefault(int(cid), {})[int(sid)] = outp

        # Helper to read data + header
        def _read_mrc_data(path: str) -> Tuple[np.ndarray, Dict[str,float]]:
            hdr = IO_utils.read_mrc_header(path)
            nx, ny, nz = hdr["nx"], hdr["ny"], hdr["nz"]
            with open(path, "rb") as f:
                f.seek(1024 + hdr.get("nsymbt", 0))
                arr = np.fromfile(f, dtype=np.float32, count=nx * ny * nz)
            return arr.reshape((nz, ny, nx)), hdr

        avg_out_dir = preferred_dir if preferred_dir else Path(args.all_recs).expanduser().resolve()
        avg_out_dir.mkdir(parents=True, exist_ok=True)
        for cls, subs in sorted(by_class.items()):
            if 1 in subs and 2 in subs:
                p1 = subs[1]; p2 = subs[2]
                v1, h1 = _read_mrc_data(p1)
                v2, _  = _read_mrc_data(p2)
                vavg = 0.5 * (v1 + v2)

                # Optional blur on averaged map (CPU)
                if args.postprocessBlur is not None and float(args.postprocessBlur) > 0:
                    try:
                        angpix_avg = float(janas_core.spacingMRC(p1))
                        vavg = janas_mapProcess.gaussian_blur_array(
                            vavg, angpix=angpix_avg, sigmaA=float(args.postprocessBlur), device='cpu'
                        )
                    except Exception as e:
                        print(f"[post-blur] WARNING on average map: {e}")

                if args.basename:
                    avg_name = f"{args.basename}_class_{cls}_rec.mrc"
                else:
                    avg_name = f"class_{cls}_rec.mrc"
                avg_path = avg_out_dir / avg_name
                IO_utils.write_mrc_like(str(avg_path), vavg, h1, update_stats=True)




                IO_utils.write_mrc_like(str(avg_path), vavg, h1, update_stats=True)
                print(f"Wrote average map {avg_path}")

    # -------- Local resolution: only if requested and both halves were reconstructed --------
    if getattr(args, "locres", False):
        print("doing locres")
        locres_utils.run_locres_for_all_pairs(results, preferred_dir, args)



if __name__ == "__main__":
    main()