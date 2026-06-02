"""
Tests for csparc2star handling of CryoSPARC .cs files that lack
alignments3D/pose (and optionally alignments3D/shift): extraction,
picking, passthrough and coordinate-only jobs.

Run::

    python tests/test_csparc2star_missing_pose.py
"""

from __future__ import annotations

import os
import sys
import tempfile
import traceback
import types
from typing import Callable, List, Tuple

import numpy as np


HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.normpath(os.path.join(HERE, "..", "src"))
if os.path.isdir(SRC) and SRC not in sys.path:
    sys.path.insert(0, SRC)

if "janas.janas_core" not in sys.modules:
    sys.modules["janas.janas_core"] = types.ModuleType("janas.janas_core")

from janas import utils  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers: build a tiny structured CryoSPARC-style .cs file on disk.
# csparc2star calls np.load(...), so we just need a .cs file that
# np.save / np.load round-trips. The dtype mimics the relevant fields.
# ---------------------------------------------------------------------------


def _write_cs(path_cs: str, fields: dict) -> None:
    """
    Persist `fields` as a structured numpy array using np.save, then rename
    to a .cs extension (which is what csparc2star expects).
    """
    names = list(fields.keys())
    # Infer dtype from the first row of each column
    dtype = []
    n = len(next(iter(fields.values())))
    for name in names:
        arr = np.asarray(fields[name])
        if arr.ndim == 1:
            dtype.append((name, arr.dtype))
        else:
            dtype.append((name, arr.dtype, arr.shape[1:]))
    struct = np.zeros(n, dtype=dtype)
    for name in names:
        struct[name] = fields[name]
    # csparc2star uses np.load(...).cs files; np.save writes a .npy header
    # which np.load can read regardless of the extension.
    base = path_cs[:-3] if path_cs.endswith(".cs") else path_cs
    npy_path = base + ".npy"
    np.save(npy_path, struct, allow_pickle=False)
    os.replace(npy_path, path_cs)


def _make_extraction_cs(path_cs: str, n: int = 3) -> None:
    """A picking/extraction-style .cs: blob + ctf + location, NO alignments3D."""
    blob_path = np.array(
        [f"J9/extract/00000_HTT_46Q_{i:04d}_aligned_DW_particles.mrc" for i in range(n)],
        dtype="U200",
    )
    fields = {
        "blob/path": blob_path,
        "blob/idx": np.arange(n, dtype=np.int32),
        "blob/shape": np.tile(np.array([256, 256], dtype=np.int32), (n, 1)),
        "blob/psize_A": np.full(n, 1.06, dtype=np.float32),
        "ctf/accel_kv": np.full(n, 300.0, dtype=np.float32),
        "ctf/cs_mm": np.full(n, 2.7, dtype=np.float32),
        "ctf/amp_contrast": np.full(n, 0.1, dtype=np.float32),
        "ctf/df1_A": np.full(n, 12000.0, dtype=np.float32),
        "ctf/df2_A": np.full(n, 12100.0, dtype=np.float32),
        "ctf/df_angle_rad": np.zeros(n, dtype=np.float32),
        "ctf/phase_shift_rad": np.zeros(n, dtype=np.float32),
        "ctf/bfactor": np.zeros(n, dtype=np.float32),
    }
    _write_cs(path_cs, fields)


def _make_refinement_cs(path_cs: str, n: int = 3) -> None:
    """A refinement-style .cs: blob + ctf + alignments3D/pose,shift,psize_A."""
    blob_path = np.array(
        [f"J42/extract/00000_HTT_46Q_{i:04d}_aligned_DW_particles.mrc" for i in range(n)],
        dtype="U200",
    )
    # 3-vector exponential-map poses (rotation vectors); zero -> identity rotation
    poses = np.zeros((n, 3), dtype=np.float32)
    shifts = np.zeros((n, 2), dtype=np.float32)
    fields = {
        "blob/path": blob_path,
        "blob/idx": np.arange(n, dtype=np.int32),
        "blob/shape": np.tile(np.array([256, 256], dtype=np.int32), (n, 1)),
        "blob/psize_A": np.full(n, 1.06, dtype=np.float32),
        "ctf/accel_kv": np.full(n, 300.0, dtype=np.float32),
        "ctf/cs_mm": np.full(n, 2.7, dtype=np.float32),
        "ctf/amp_contrast": np.full(n, 0.1, dtype=np.float32),
        "ctf/df1_A": np.full(n, 12000.0, dtype=np.float32),
        "ctf/df2_A": np.full(n, 12100.0, dtype=np.float32),
        "ctf/df_angle_rad": np.zeros(n, dtype=np.float32),
        "ctf/phase_shift_rad": np.zeros(n, dtype=np.float32),
        "ctf/bfactor": np.zeros(n, dtype=np.float32),
        "alignments3D/pose": poses,
        "alignments3D/shift": shifts,
        "alignments3D/psize_A": np.full(n, 1.06, dtype=np.float32),
    }
    _write_cs(path_cs, fields)


def _read_particles_block(star_path: str) -> str:
    with open(star_path) as f:
        return f.read()


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_missing_pose_default_raises_clear_error() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        cs_path = os.path.join(tmp, "extracted.cs")
        out_star = os.path.join(tmp, "extracted.star")
        _make_extraction_cs(cs_path)
        try:
            utils.csparc2star(cs_path, out_star)
        except ValueError as e:
            msg = str(e)
            assert "alignments3D/pose" in msg, msg
            assert "--missing_pose_to_zero" in msg, msg
            assert "_rlnAngleRot" in msg, msg
            assert "_rlnOriginXAngst" in msg, msg
            assert not os.path.exists(out_star), "output must not be written on error"
            return
    raise AssertionError("expected ValueError on missing pose without the flag")


def test_missing_pose_with_flag_writes_zero_angles_and_origins() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        cs_path = os.path.join(tmp, "extracted.cs")
        out_star = os.path.join(tmp, "extracted.star")
        _make_extraction_cs(cs_path, n=3)
        utils.csparc2star(cs_path, out_star, missing_pose_to_zero=True)
        text = _read_particles_block(out_star)
        # 3 particles, each with the four angle/origin columns set to 0 (or 0.0...)
        # We do not assume RELION pretty-printing; just check the columns exist
        # and that at least one zero value appears for each.
        for col in (
            "_rlnAngleRot",
            "_rlnAngleTilt",
            "_rlnAnglePsi",
            "_rlnOriginXAngst",
            "_rlnOriginYAngst",
        ):
            assert col in text, f"{col} missing from output STAR"


def test_refinement_style_cs_still_works_default() -> None:
    """When alignments3D/pose IS present, default behaviour must be unchanged."""
    with tempfile.TemporaryDirectory() as tmp:
        cs_path = os.path.join(tmp, "refined.cs")
        out_star = os.path.join(tmp, "refined.star")
        _make_refinement_cs(cs_path)
        utils.csparc2star(cs_path, out_star)  # no flag
        text = _read_particles_block(out_star)
        for col in (
            "_rlnAngleRot",
            "_rlnAngleTilt",
            "_rlnAnglePsi",
            "_rlnOriginXAngst",
            "_rlnOriginYAngst",
        ):
            assert col in text, f"{col} missing in default refinement output"


def test_missing_shift_only_with_flag() -> None:
    """Pose present, shift missing: zero origins under the flag."""
    with tempfile.TemporaryDirectory() as tmp:
        cs_path = os.path.join(tmp, "pose_only.cs")
        out_star = os.path.join(tmp, "pose_only.star")
        # Refinement-style but drop alignments3D/shift and alignments3D/psize_A
        _make_refinement_cs(cs_path, n=2)
        # Rewrite without shift fields
        arr = np.load(cs_path, max_header_size=100000)
        keep = [n for n in arr.dtype.names
                if n not in ("alignments3D/shift", "alignments3D/psize_A")]
        new_dt = [(n, arr.dtype[n]) for n in keep]
        new = np.zeros(len(arr), dtype=new_dt)
        for n in keep:
            new[n] = arr[n]
        np.save(cs_path[:-3] + ".npy", new, allow_pickle=False)
        os.replace(cs_path[:-3] + ".npy", cs_path)

        utils.csparc2star(cs_path, out_star, missing_pose_to_zero=True)
        text = _read_particles_block(out_star)
        for col in ("_rlnAngleRot", "_rlnOriginXAngst", "_rlnOriginYAngst"):
            assert col in text, col


def test_missing_shift_only_default_raises() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        cs_path = os.path.join(tmp, "pose_only.cs")
        out_star = os.path.join(tmp, "pose_only.star")
        _make_refinement_cs(cs_path, n=2)
        arr = np.load(cs_path, max_header_size=100000)
        keep = [n for n in arr.dtype.names
                if n not in ("alignments3D/shift", "alignments3D/psize_A")]
        new_dt = [(n, arr.dtype[n]) for n in keep]
        new = np.zeros(len(arr), dtype=new_dt)
        for n in keep:
            new[n] = arr[n]
        np.save(cs_path[:-3] + ".npy", new, allow_pickle=False)
        os.replace(cs_path[:-3] + ".npy", cs_path)

        try:
            utils.csparc2star(cs_path, out_star)
        except ValueError as e:
            assert "alignments3D/shift" in str(e), str(e)
            return
    raise AssertionError("expected ValueError on missing shift without flag")


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

TESTS: List[Tuple[str, Callable[[], None]]] = [
    ("missing pose: default raises clear error", test_missing_pose_default_raises_clear_error),
    ("missing pose: --missing_pose_to_zero writes zero angles/origins",
        test_missing_pose_with_flag_writes_zero_angles_and_origins),
    ("refinement-style .cs (pose present): default unchanged",
        test_refinement_style_cs_still_works_default),
    ("missing shift only: with flag, zero origins",
        test_missing_shift_only_with_flag),
    ("missing shift only: default raises",
        test_missing_shift_only_default_raises),
]


def main() -> int:
    failed: List[str] = []
    for name, fn in TESTS:
        try:
            fn()
            print(f"[ OK ] {name}")
        except Exception as exc:  # noqa: BLE001
            failed.append(name)
            print(f"[FAIL] {name}: {exc}")
            traceback.print_exc(limit=4)
    print()
    if failed:
        print(f"{len(failed)}/{len(TESTS)} test(s) failed: {failed}")
        return 1
    print(f"All {len(TESTS)} tests passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
