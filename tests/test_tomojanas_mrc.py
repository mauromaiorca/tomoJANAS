#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Self-contained tests for tomojanas.io.mrc.

Run directly:  python tests/test_tomojanas_mrc.py
(or via pytest:  pytest tests/test_tomojanas_mrc.py)
"""
import os
import sys
import tempfile
from typing import Callable, List, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from tomojanas.io import mrc  # noqa: E402


def _approx(a, b, tol=1e-4):
    return abs(float(a) - float(b)) <= tol


def test_mode_mapping():
    assert mrc.mrc_dtype_from_mode(0) == np.dtype(np.int8)
    assert mrc.mrc_dtype_from_mode(1) == np.dtype(np.int16)
    assert mrc.mrc_dtype_from_mode(2) == np.dtype(np.float32)
    assert mrc.mrc_dtype_from_mode(6) == np.dtype(np.uint16)
    assert mrc.mrc_dtype_from_mode(12) == np.dtype(np.float16)
    for bad in (3, 4):
        try:
            mrc.mrc_dtype_from_mode(bad)
        except NotImplementedError:
            pass
        else:
            raise AssertionError(f"mode {bad} should raise")


def test_write_read_roundtrip():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "vol.mrc")
        vol = np.arange(4 * 5 * 6, dtype=np.float32).reshape(4, 5, 6)  # (z,y,x)
        mrc.write_mrc(p, vol, pixel_size=3.42, origin_angs=(1.0, 2.0, 3.0))
        hdr = mrc.read_mrc_header(p)
        assert (hdr.nx, hdr.ny, hdr.nz) == (6, 5, 4)
        assert _approx(hdr.pixel_x, 3.42)
        assert _approx(hdr.origin_x, 1.0) and _approx(hdr.origin_z, 3.0)
        assert hdr.data_offset == 1024
        back, hdr2 = mrc.read_mrc_data(p)
        assert back.shape == (4, 5, 6)
        assert np.allclose(back, vol)
        assert hdr2.shape == (4, 5, 6)


def test_as_dict_legacy_keys():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "vol.mrc")
        mrc.write_mrc(p, np.zeros((2, 3, 4), np.float32), 1.5)
        hdr = mrc.read_mrc_header(p)
        dd = hdr.as_dict()
        expected = {
            "nx", "ny", "nz", "mode", "nxstart", "nystart", "nzstart",
            "mx", "my", "mz", "xlen", "ylen", "zlen", "alpha", "beta", "gamma",
            "mapc", "mapr", "maps", "amin", "amax", "amean", "ispg", "nsymbt",
            "origin_x", "origin_y", "origin_z", "pixel_x", "pixel_y", "pixel_z",
            "stamp",
        }
        assert expected.issubset(set(dd.keys())), set(dd.keys()) ^ expected
        assert mrc.MRCHeader.from_dict(dd).nx == hdr.nx


def test_write_mrc_like():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "vol.mrc")
        q = os.path.join(d, "like.mrc")
        vol = np.random.rand(3, 4, 5).astype(np.float32)
        mrc.write_mrc(p, vol, 2.0, origin_angs=(5, 6, 7))
        hdr = mrc.read_mrc_header(p)
        vol2 = np.random.rand(3, 4, 5).astype(np.float32)
        mrc.write_mrc_like(q, vol2, hdr)
        h2 = mrc.read_mrc_header(q)
        assert _approx(h2.pixel_x, 2.0) and _approx(h2.origin_y, 6.0)
        back, _ = mrc.read_mrc_data(q)
        assert np.allclose(back, vol2)


def _write_with_nsymbt(path, vol, pixel_size, nsymbt):
    """Helper: write an MRC with a non-zero extended header (nsymbt bytes)."""
    nz, ny, nx = vol.shape
    header = mrc._pack_header(
        nx=nx, ny=ny, nz=nz, mode=2,
        pixel_size_xyz=(pixel_size, pixel_size, pixel_size),
        nsymbt=nsymbt,
    )
    with open(path, "wb") as f:
        f.write(header)
        f.write(b"\xAB" * nsymbt)  # fake extended header payload
        f.write(np.asarray(vol, np.float32).tobytes())


def test_nsymbt_offset_read_and_memmap():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "ext.mrc")
        vol = np.arange(2 * 3 * 4, dtype=np.float32).reshape(2, 3, 4)
        _write_with_nsymbt(p, vol, 1.0, nsymbt=128)
        hdr = mrc.read_mrc_header(p)
        assert hdr.nsymbt == 128
        assert hdr.data_offset == 1024 + 128
        back, _ = mrc.read_mrc_data(p)
        assert np.allclose(back, vol), "nsymbt offset must be honoured on read"
        mm, h2 = mrc.open_mrc_memmap(p)
        assert np.allclose(np.asarray(mm), vol)
        del mm


def test_read_slice_with_nsymbt():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "stack.mrc")
        vol = np.stack([np.full((3, 4), k, np.float32) for k in range(5)])  # (5,3,4)
        _write_with_nsymbt(p, vol, 1.0, nsymbt=80)
        for k in range(5):
            img, _ = mrc.read_mrc_slice(p, k)
            assert img.shape == (3, 4)
            assert np.allclose(img, k)
        try:
            mrc.read_mrc_slice(p, 5)
        except IndexError:
            pass
        else:
            raise AssertionError("out-of-range slice must raise")


def test_precreate_and_write_slice():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "pre.mrcs")
        mrc.precreate_mrc_stack(p, nz=3, ny=4, nx=5, pixel_size=2.0)
        for k in range(3):
            mrc.write_mrc_slice(p, k, np.full((4, 5), 10 + k, np.float32))
        back, hdr = mrc.read_mrc_data(p)
        assert back.shape == (3, 4, 5)
        for k in range(3):
            assert np.allclose(back[k], 10 + k)
        assert _approx(hdr.pixel_x, 2.0)


def test_modes_int16_uint16_float16():
    with tempfile.TemporaryDirectory() as d:
        for mode, dt, val in ((1, np.int16, 1234), (6, np.uint16, 60000), (12, np.float16, 1.5)):
            p = os.path.join(d, f"m{mode}.mrc")
            vol = np.full((2, 2, 2), val, dtype=dt)
            header = mrc._pack_header(nx=2, ny=2, nz=2, mode=mode,
                                      pixel_size_xyz=(1.0, 1.0, 1.0))
            with open(p, "wb") as f:
                f.write(header)
                f.write(vol.tobytes())
            back, hdr = mrc.read_mrc_data(p)
            assert hdr.mode == mode
            assert back.dtype == np.float32
            assert np.allclose(back, float(val))


def test_crop_box_inside_and_partial():
    vol = np.arange(10 * 10 * 10, dtype=np.float32).reshape(10, 10, 10)
    # fully inside
    c = mrc.crop_volume_box(vol, (5, 5, 5), 4)
    assert c.shape == (4, 4, 4)
    z0, y0, x0 = mrc.box_corner((5, 5, 5), 4)
    assert np.allclose(c, vol[z0:z0 + 4, y0:y0 + 4, x0:x0 + 4])
    # partially outside -> padded with pad_value
    c2 = mrc.crop_volume_box(vol, (0, 0, 0), 4, pad_value=-7.0)
    assert c2.shape == (4, 4, 4)
    assert (c2 == -7.0).any(), "out-of-bounds region must be padded"
    assert (c2 != -7.0).any(), "in-bounds region must be copied"


def test_crop_sphere_and_mask():
    vol = np.ones((20, 20, 20), dtype=np.float32)
    cube = mrc.crop_volume_sphere(vol, (10, 10, 10), radius_vox=4, padding_vox=1)
    assert cube.shape[0] == cube.shape[1] == cube.shape[2]
    assert cube.shape[0] == 10  # ceil(2*(4+1))
    m = mrc.make_spherical_mask(cube.shape[0], radius_vox=4)
    assert m.dtype == np.bool_
    assert m[cube.shape[0] // 2, cube.shape[0] // 2, cube.shape[0] // 2]
    assert not m[0, 0, 0]


def test_crop_image_circle():
    img = np.ones((30, 30), dtype=np.float32)
    patch = mrc.crop_image_circle(img, (15, 15), radius_px=5)
    assert patch.shape == (10, 10)
    cm = mrc.make_circular_mask(10, radius_px=5)
    assert cm[5, 5] and not cm[0, 0]


def test_write_cropped_mrc_like_origin():
    with tempfile.TemporaryDirectory() as d:
        src = os.path.join(d, "src.mrc")
        out = os.path.join(d, "crop.mrc")
        vol = np.random.rand(20, 20, 20).astype(np.float32)
        mrc.write_mrc(src, vol, pixel_size=2.0, origin_angs=(0, 0, 0))
        hdr = mrc.read_mrc_header(src)
        origin_zyx = (4, 5, 6)
        cube = vol[4:12, 5:13, 6:14]
        mrc.write_cropped_mrc_like(out, cube, hdr, crop_origin_zyx=origin_zyx)
        h2 = mrc.read_mrc_header(out)
        assert (h2.nx, h2.ny, h2.nz) == (8, 8, 8)
        assert _approx(h2.pixel_x, 2.0)
        # origin shifted by corner * pixel
        assert _approx(h2.origin_x, 6 * 2.0)
        assert _approx(h2.origin_y, 5 * 2.0)
        assert _approx(h2.origin_z, 4 * 2.0)
        back, _ = mrc.read_mrc_data(out)
        assert np.allclose(back, cube)


def test_validators():
    with tempfile.TemporaryDirectory() as d:
        a = os.path.join(d, "a.mrc")
        b = os.path.join(d, "b.mrc")
        c = os.path.join(d, "c.mrc")
        mrc.write_mrc(a, np.zeros((4, 4, 4), np.float32), 1.0)
        mrc.write_mrc(b, np.zeros((4, 4, 4), np.float32), 1.0)
        mrc.write_mrc(c, np.zeros((4, 4, 5), np.float32), 1.0)
        assert mrc.validate_mrc_geometry(a, b)["ok"]
        assert not mrc.validate_mrc_geometry(a, c)["ok"]
        assert mrc.validate_mrc_pixel_size(a, 1.0, 1e-3)["ok"]
        assert not mrc.validate_mrc_pixel_size(a, 2.0, 1e-3)["ok"]


TESTS: List[Tuple[str, Callable]] = [
    ("test_mode_mapping", test_mode_mapping),
    ("test_write_read_roundtrip", test_write_read_roundtrip),
    ("test_as_dict_legacy_keys", test_as_dict_legacy_keys),
    ("test_write_mrc_like", test_write_mrc_like),
    ("test_nsymbt_offset_read_and_memmap", test_nsymbt_offset_read_and_memmap),
    ("test_read_slice_with_nsymbt", test_read_slice_with_nsymbt),
    ("test_precreate_and_write_slice", test_precreate_and_write_slice),
    ("test_modes_int16_uint16_float16", test_modes_int16_uint16_float16),
    ("test_crop_box_inside_and_partial", test_crop_box_inside_and_partial),
    ("test_crop_sphere_and_mask", test_crop_sphere_and_mask),
    ("test_crop_image_circle", test_crop_image_circle),
    ("test_write_cropped_mrc_like_origin", test_write_cropped_mrc_like_origin),
    ("test_validators", test_validators),
]


def main() -> int:
    failed = []
    for name, fn in TESTS:
        try:
            fn()
            print(f"PASS {name}")
        except Exception as exc:  # noqa: BLE001
            failed.append(name)
            print(f"FAIL {name}: {exc!r}")
    print(f"\n{len(TESTS) - len(failed)}/{len(TESTS)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
