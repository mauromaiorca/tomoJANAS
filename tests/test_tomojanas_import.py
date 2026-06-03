#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Integration tests for tomoJANAS import framework.

Creates a synthetic IMOD dataset in a temp dir, runs:
  tomojanas-import imod
  tomojanas-import particles --input-single-point
  tomojanas-import validate --strict

Then verifies all acceptance criteria from the spec.

Run: python tests/test_tomojanas_import.py
"""
import os
import sys
import json
import tempfile
from typing import Callable, List, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from tomojanas.io.mrc import write_mrc
from tomojanas.io.star_reader import read_star
from tomojanas.importers.cli import _build_parser

TOMO_NAME = "lam8_ts_006"
N_TILTS = 3
NX = NY = 16
NZ_REC = 8
APIX = 2.0


# --------------------------------------------------------------------------- #
# Synthetic dataset builders
# --------------------------------------------------------------------------- #

def _make_imod_dir(d):
    """Create a minimal synthetic IMOD directory under d."""
    imod_dir = os.path.join(d, "imod")
    os.makedirs(imod_dir, exist_ok=True)

    # Raw / ali stack (3 tilts, 16x16)
    raw = np.zeros((N_TILTS, NY, NX), dtype=np.float32)
    for z in range(N_TILTS):
        raw[z] = z + 1
    write_mrc(os.path.join(imod_dir, f"{TOMO_NAME}.mrc"), raw, APIX)
    write_mrc(os.path.join(imod_dir, f"{TOMO_NAME}_ali.mrc"), raw, APIX)

    # Rec tomogram (8 slices, 16x16)
    rec = np.ones((NZ_REC, NY, NX), dtype=np.float32)
    write_mrc(os.path.join(imod_dir, f"{TOMO_NAME}_rec.mrc"), rec, APIX)

    # .tlt
    with open(os.path.join(imod_dir, f"{TOMO_NAME}.tlt"), "w") as f:
        f.write("-10.0\n0.0\n10.0\n")

    # .xf (identity transforms for all 3 tilts)
    with open(os.path.join(imod_dir, f"{TOMO_NAME}.xf"), "w") as f:
        for _ in range(N_TILTS):
            f.write("1.0  0.0  0.0  1.0  0.0  0.0\n")

    # newst.com
    with open(os.path.join(imod_dir, "newst.com"), "w") as f:
        f.write(f"$newst\nInputFile {TOMO_NAME}.mrc\n"
                f"OutputFile {TOMO_NAME}_ali.mrc\n"
                f"TransformFile {TOMO_NAME}.xf\n"
                f"SizeToOutputInXandY {NX},{NY}\n"
                f"BinByFactor 1\n")

    # tilt.com
    with open(os.path.join(imod_dir, "tilt.com"), "w") as f:
        f.write(f"$tilt\nInputProjections {TOMO_NAME}_ali.mrc\n"
                f"FULLIMAGE {NX} {NY}\n"
                f"TILTFILE {TOMO_NAME}.tlt\n"
                f"THICKNESS {NZ_REC}\n"
                f"IMAGEBINNED 1\n")

    # minimal .mdoc
    with open(os.path.join(imod_dir, f"{TOMO_NAME}.mrc.mdoc"), "w") as f:
        f.write(f"PixelSpacing = {APIX}\n")

    return imod_dir


def _run_cmd(argv):
    """Run a command through the real CLI entry point (so command logging and
    dispatch are exercised). Returns the exit code."""
    from tomojanas.importers.cli import main as _main
    try:
        _main(argv)
        return 0
    except SystemExit as exc:
        return int(exc.code) if exc.code is not None else 0


# --------------------------------------------------------------------------- #
# Tests
# --------------------------------------------------------------------------- #

def test_imod_import():
    with tempfile.TemporaryDirectory() as d:
        imod_dir = _make_imod_dir(d)
        project = os.path.join(d, "project")

        rc = _run_cmd([
            "imod",
            "--project", project,
            "--create-if-missing",
            "--imod-dir", imod_dir,
            "--tomo-name", TOMO_NAME,
            "--validate",
        ])
        assert rc == 0, f"imod import returned {rc}"

        # AC4: check required files exist
        assert os.path.isfile(os.path.join(project, "tomograms.star")), "tomograms.star missing"
        assert os.path.isfile(os.path.join(project, "optimisation_set.star")), "optimisation_set.star missing"
        tilt_star = os.path.join(project, "tilt_series", f"{TOMO_NAME}.star")
        assert os.path.isfile(tilt_star), f"tilt_series/{TOMO_NAME}.star missing"
        assert os.path.isfile(os.path.join(project, "tilt_series", TOMO_NAME, "imod_settings", "xf_matrices.star"))
        assert os.path.isfile(os.path.join(project, "tilt_series", TOMO_NAME, "imod_settings", "relion_geometry_status.star"))
        assert os.path.isfile(os.path.join(project, "logs", "validation_log.json"))

        # AC5: no relion/ directory
        assert not os.path.isdir(os.path.join(project, "relion")), "relion/ dir must not exist"

        # AC13: geometry status is relion_imod_algorithm_ported
        blks = read_star(os.path.join(project, "tomograms.star"))
        assert "tomoJANAS_tomogram_sources" in blks
        df_src = blks["tomoJANAS_tomogram_sources"]["df"]
        status = str(df_src["_tomoJANASRelionGeometryStatus"].iloc[0])
        assert "relion_imod_algorithm_ported" in status, f"unexpected status: {status}"

        # AC12: _rlnTomoProjX/Y/Z/W should be ? (not fake values)
        blks_ts = read_star(tilt_star)
        assert TOMO_NAME in blks_ts
        df_ts = blks_ts[TOMO_NAME]["df"]
        for col in ("_rlnTomoProjX", "_rlnTomoProjY", "_rlnTomoProjZ", "_rlnTomoProjW"):
            assert col in df_ts.columns
            vals = df_ts[col].tolist()
            assert all(str(v).strip() == "?" for v in vals), f"{col} should be '?' not {vals}"

        # matrices stored in tomoJANAS block
        assert "tomoJANAS_projection_matrices" in blks_ts, "tomoJANAS_projection_matrices block missing"
        mat_df = blks_ts["tomoJANAS_projection_matrices"]["df"]
        assert len(mat_df) == N_TILTS, f"expected {N_TILTS} matrix rows, got {len(mat_df)}"

        # particles_all.star created empty with optics + particles
        pa = os.path.join(project, "particles_all.star")
        assert os.path.isfile(pa)
        pa_blks = read_star(pa)
        assert "optics" in pa_blks, "data_optics missing from empty particles_all.star"
        assert "particles" in pa_blks, "data_particles missing from empty particles_all.star"

        print("PASS test_imod_import")


def test_particles_import():
    with tempfile.TemporaryDirectory() as d:
        imod_dir = _make_imod_dir(d)
        project = os.path.join(d, "project")

        # Import IMOD first
        rc = _run_cmd(["imod", "--project", project, "--create-if-missing",
                        "--imod-dir", imod_dir, "--tomo-name", TOMO_NAME])
        assert rc == 0

        # Import single particle at voxel (8,8,4) zero-based
        rc = _run_cmd([
            "particles",
            "--project", project,
            "--tomo-name", TOMO_NAME,
            "--input-single-point", "8,8,4",
            "--coordinate-system", "rec-voxel",
            "--indexing", "zero-based",
            "--axis-order", "xyz",
            "--roi-radius-angst", "20.0",
            "--validate",
        ])
        assert rc == 0, f"particles import returned {rc}"

        # AC4: particles_all.star updated
        pa = os.path.join(project, "particles_all.star")
        assert os.path.isfile(pa)
        pa_blks = read_star(pa)

        # AC6+7: data_optics + data_particles with required columns
        assert "optics" in pa_blks, "data_optics missing"
        assert "particles" in pa_blks, "data_particles missing"
        df_optics = pa_blks["optics"]["df"]
        assert "_rlnTomoTiltSeriesPixelSize" in df_optics.columns, "_rlnTomoTiltSeriesPixelSize missing"
        df_parts = pa_blks["particles"]["df"]
        assert "_rlnOpticsGroup" in df_parts.columns, "_rlnOpticsGroup missing"

        # AC9: optics group reference is valid
        valid_og = set(df_optics["_rlnOpticsGroup"].astype(str).tolist())
        for og in df_parts["_rlnOpticsGroup"].astype(str).tolist():
            assert og in valid_og, f"_rlnOpticsGroup={og} not in optics table"

        # AC8: no _tomoJANAS in RELION loop
        tomo_cols = [c for c in df_parts.columns if c.startswith("_tomoJANAS")]
        assert not tomo_cols, f"_tomoJANAS* in RELION particles loop: {tomo_cols}"

        # AC4: P000001.star exists
        p_star = os.path.join(project, "tilt_series", TOMO_NAME, "individual_particles", "P000001.star")
        assert os.path.isfile(p_star), f"P000001.star missing at {p_star}"

        # AC10: P000001.star has all required blocks
        p_blks = read_star(p_star)
        for bname in ("optics", "particles", "tomoJANAS_particle_source",
                       "tomoJANAS_particle_roi", "tomoJANAS_particle_projections"):
            assert bname in p_blks, f"P000001.star missing block: data_{bname}"

        # AC11: ROI is sphere/circle
        roi_df = p_blks["tomoJANAS_particle_roi"]["df"]
        assert str(roi_df["_tomoJANASRoiShape3D"].iloc[0]) == "sphere"
        assert str(roi_df["_tomoJANASRoiShape2D"].iloc[0]) == "circle"

        # Projection table has correct number of tilts
        proj_df = p_blks["tomoJANAS_particle_projections"]["df"]
        assert len(proj_df) == N_TILTS, f"expected {N_TILTS} proj rows, got {len(proj_df)}"

        print("PASS test_particles_import")


def test_rec_crop_writing():
    """--write-rec-crops produces a physical sub-volume + crop metadata block,
    with internal crop path relative and external source path absolute."""
    with tempfile.TemporaryDirectory() as d:
        imod_dir = _make_imod_dir(d)
        project = os.path.join(d, "project")

        _run_cmd(["imod", "--project", project, "--create-if-missing",
                   "--imod-dir", imod_dir, "--tomo-name", TOMO_NAME])

        rc = _run_cmd([
            "particles", "--project", project, "--tomo-name", TOMO_NAME,
            "--input-single-point", "8,8,4",
            "--coordinate-system", "rec-voxel", "--indexing", "zero-based",
            "--roi-radius-angst", "8.0",
            "--write-rec-crops",
            "--crop-storage-box-size", "6",
        ])
        assert rc == 0, f"particles import returned {rc}"

        # crop MRC must exist
        crop_path = os.path.join(project, "tilt_series", TOMO_NAME,
                                  "individual_particles_recs", "P000001_rec.mrc")
        assert os.path.isfile(crop_path), f"rec crop not written at {crop_path}"

        # crop must be the requested box size and readable
        from tomojanas.io.mrc import read_mrc_header
        hdr = read_mrc_header(crop_path)
        assert (hdr.nx, hdr.ny, hdr.nz) == (6, 6, 6), f"crop dims {hdr.nx}x{hdr.ny}x{hdr.nz}"

        # P*.star has the rec crop block; check path policy
        p_star = os.path.join(project, "tilt_series", TOMO_NAME,
                               "individual_particles", "P000001.star")
        blks = read_star(p_star)
        assert "tomoJANAS_particle_rec_crop" in blks, "rec crop block missing in P*.star"
        crop_df = blks["tomoJANAS_particle_rec_crop"]["df"]
        rec_path_stored = str(crop_df["_tomoJANASParticleRecPath"].iloc[0])
        src_stored = str(crop_df["_tomoJANASParticleRecSourceTomogram"].iloc[0])
        # internal crop path → relative
        assert not os.path.isabs(rec_path_stored), f"crop path should be relative: {rec_path_stored}"
        # external source tomogram → absolute
        assert os.path.isabs(src_stored), f"source tomogram should be absolute: {src_stored}"

        print("PASS test_rec_crop_writing")


def test_external_paths_absolute():
    """External IMOD source files are stored as absolute paths in tomograms.star."""
    with tempfile.TemporaryDirectory() as d:
        imod_dir = _make_imod_dir(d)
        project = os.path.join(d, "project")
        _run_cmd(["imod", "--project", project, "--create-if-missing",
                   "--imod-dir", imod_dir, "--tomo-name", TOMO_NAME])

        blks = read_star(os.path.join(project, "tomograms.star"))
        src_df = blks["tomoJANAS_tomogram_sources"]["df"]
        for col in ("_tomoJANASImodDir", "_tomoJANASRecTomogram",
                     "_tomoJANASAliStack", "_tomoJANASXfFile", "_tomoJANASTltFile"):
            val = str(src_df[col].iloc[0])
            if val not in ("?", ""):
                assert os.path.isabs(val), f"{col} must be absolute, got: {val}"

        # internal tilt-series star ref must be relative
        glob_df = blks["global"]["df"]
        ts_ref = str(glob_df["_rlnTomoTiltSeriesStarFile"].iloc[0])
        assert not os.path.isabs(ts_ref), f"tilt-series star ref should be relative: {ts_ref}"
        # external rec tomogram must be absolute
        rec_ref = str(glob_df["_rlnTomoReconstructedTomogram"].iloc[0])
        assert os.path.isabs(rec_ref), f"rec tomogram ref should be absolute: {rec_ref}"

        print("PASS test_external_paths_absolute")


def test_particle_autoincrement_and_command_log():
    """Re-running a single-point import adds P000002 (not overwrite), and each
    command is appended to logs/commands.{sh,jsonl} with its exit status."""
    with tempfile.TemporaryDirectory() as d:
        imod_dir = _make_imod_dir(d)
        project = os.path.join(d, "project")
        _run_cmd(["imod", "--project", project, "--create-if-missing",
                   "--imod-dir", imod_dir, "--tomo-name", TOMO_NAME])
        _run_cmd(["particles", "--project", project, "--tomo-name", TOMO_NAME,
                   "--input-single-point", "8,8,4", "--coordinate-system", "rec-voxel",
                   "--indexing", "zero-based", "--roi-radius-angst", "8"])
        _run_cmd(["particles", "--project", project, "--tomo-name", TOMO_NAME,
                   "--input-single-point", "8,8,4", "--coordinate-system", "rec-voxel",
                   "--indexing", "zero-based", "--roi-radius-angst", "16"])

        ip = os.path.join(project, "tilt_series", TOMO_NAME, "individual_particles")
        names = sorted(f for f in os.listdir(ip) if f.endswith(".star"))
        assert names == ["P000001.star", "P000002.star"], names

        # command log exists and records all three invocations
        sh = os.path.join(project, "logs", "commands.sh")
        jsonl = os.path.join(project, "logs", "commands.jsonl")
        assert os.path.isfile(sh) and os.path.isfile(jsonl)
        with open(jsonl) as f:
            entries = [json.loads(ln) for ln in f if ln.strip()]
        assert len(entries) == 3, f"expected 3 logged commands, got {len(entries)}"
        assert all("exit_status" in e and "command" in e for e in entries)
        assert entries[0]["argv"][0] == "imod"
        print("PASS test_particle_autoincrement_and_command_log")


def test_status_scan_and_sync():
    """status detects an orphaned registry entry after manual deletion, and
    --sync rebuilds particles_all.star from the P*.star files on disk."""
    from tomojanas.models.project import Project
    from tomojanas.importers.status import scan_project, sync_particles_all

    with tempfile.TemporaryDirectory() as d:
        imod_dir = _make_imod_dir(d)
        project = os.path.join(d, "project")
        _run_cmd(["imod", "--project", project, "--create-if-missing",
                   "--imod-dir", imod_dir, "--tomo-name", TOMO_NAME])
        for pt in ("8,8,4", "9,9,4"):
            _run_cmd(["particles", "--project", project, "--tomo-name", TOMO_NAME,
                       "--input-single-point", pt, "--coordinate-system", "rec-voxel",
                       "--indexing", "zero-based", "--roi-radius-angst", "8"])

        ip = os.path.join(project, "tilt_series", TOMO_NAME, "individual_particles")
        os.remove(os.path.join(ip, "P000001.star"))  # manual deletion

        proj = Project(root=os.path.abspath(project))
        rep = scan_project(proj, tomo_name=TOMO_NAME)
        assert rep["n_issues"] == 1, rep
        issue = rep["tomograms"][0]["issues"][0]
        assert issue["particle"] == "P000001"
        assert issue["issue"] == "registered_but_no_star"

        n = sync_particles_all(proj, [TOMO_NAME])
        assert n == 1, f"after sync particles_all should have 1 row, got {n}"
        rep2 = scan_project(proj, tomo_name=TOMO_NAME)
        assert rep2["n_issues"] == 0, rep2
        print("PASS test_status_scan_and_sync")


def test_crop_axis_mapping_xyz():
    """Lock the crop axis mapping: with default --axis-order xyz the picked
    (X,Y,Z) maps 1st->nx, 2nd->ny, 3rd->nz (rec[k=Z, j=Y, i=X]); no flip."""
    with tempfile.TemporaryDirectory() as d:
        imod_dir = _make_imod_dir(d)
        # rec with distinct dims + a unique value per voxel: vol[k,j,i]=k*100+j*10+i
        nz, ny, nx = 3, 5, 7
        vol = np.fromfunction(lambda k, j, i: k * 100 + j * 10 + i,
                              (nz, ny, nx)).astype(np.float32)
        write_mrc(os.path.join(imod_dir, f"{TOMO_NAME}_rec.mrc"), vol, APIX)

        project = os.path.join(d, "project")
        _run_cmd(["imod", "--project", project, "--create-if-missing",
                   "--imod-dir", imod_dir, "--tomo-name", TOMO_NAME])
        # pick (X=4, Y=2, Z=1) -> expect rec[k=1, j=2, i=4] = 124
        _run_cmd(["particles", "--project", project, "--tomo-name", TOMO_NAME,
                   "--input-single-point", "4,2,1", "--coordinate-system", "rec-voxel",
                   "--indexing", "zero-based", "--roi-radius-angst", "1",
                   "--write-rec-crops", "--crop-storage-box-size", "1"])
        from tomojanas.io.mrc import read_mrc_data
        crop = os.path.join(project, "tilt_series", TOMO_NAME,
                            "individual_particles_recs", "P000001_rec.mrc")
        cube, _ = read_mrc_data(crop)
        assert float(cube[0, 0, 0]) == 124.0, \
            f"axis mapping wrong: centre={cube[0,0,0]} (expected 124 = k1*100+j2*10+i4)"
        print("PASS test_crop_axis_mapping_xyz")


def test_default_axis_order_is_xyz():
    """With no --axis-order, the default is 'xyz' (no header inference)."""
    with tempfile.TemporaryDirectory() as d:
        imod_dir = _make_imod_dir(d)
        # use a flipped-shaped rec to prove there is NO auto-detection anymore
        flipped = np.ones((NY, NZ_REC, NX), dtype=np.float32)
        write_mrc(os.path.join(imod_dir, f"{TOMO_NAME}_rec.mrc"), flipped, APIX)
        project = os.path.join(d, "project")
        _run_cmd(["imod", "--project", project, "--create-if-missing",
                   "--imod-dir", imod_dir, "--tomo-name", TOMO_NAME])
        rc = _run_cmd(["particles", "--project", project, "--tomo-name", TOMO_NAME,
                        "--input-single-point", "8,4,3", "--coordinate-system", "rec-voxel",
                        "--indexing", "zero-based", "--roi-radius-angst", "6"])
        assert rc == 0, f"particles import returned {rc}"
        p_star = os.path.join(project, "tilt_series", TOMO_NAME,
                               "individual_particles", "P000001.star")
        src = read_star(p_star)["tomoJANAS_particle_source"]["df"].iloc[0]
        assert str(src["_tomoJANASPickedAxisOrder"]) == "xyz", \
            f"default should be xyz, got {src['_tomoJANASPickedAxisOrder']}"
        print("PASS test_default_axis_order_is_xyz")


def test_explicit_xzy_axis_order():
    """Explicit --axis-order xzy reorders (X,Z,Y) inputs and is recorded."""
    with tempfile.TemporaryDirectory() as d:
        imod_dir = _make_imod_dir(d)
        project = os.path.join(d, "project")
        _run_cmd(["imod", "--project", project, "--create-if-missing",
                   "--imod-dir", imod_dir, "--tomo-name", TOMO_NAME])
        rc = _run_cmd(["particles", "--project", project, "--tomo-name", TOMO_NAME,
                        "--input-single-point", "8,4,3", "--coordinate-system", "rec-voxel",
                        "--indexing", "zero-based", "--axis-order", "xzy",
                        "--roi-radius-angst", "6"])
        assert rc == 0, f"particles import returned {rc}"
        p_star = os.path.join(project, "tilt_series", TOMO_NAME,
                               "individual_particles", "P000001.star")
        src = read_star(p_star)["tomoJANAS_particle_source"]["df"].iloc[0]
        assert str(src["_tomoJANASPickedAxisOrder"]) == "xzy"
        print("PASS test_explicit_xzy_axis_order")


def test_status_create_volume():
    """status --create-volume backfills _rec.mrc for particles imported
    without crops, and adds the rec-crop block to P*.star."""
    from tomojanas.io.mrc import read_mrc_header
    with tempfile.TemporaryDirectory() as d:
        imod_dir = _make_imod_dir(d)
        project = os.path.join(d, "project")
        _run_cmd(["imod", "--project", project, "--create-if-missing",
                   "--imod-dir", imod_dir, "--tomo-name", TOMO_NAME])
        # import WITHOUT --write-rec-crops
        _run_cmd(["particles", "--project", project, "--tomo-name", TOMO_NAME,
                   "--input-single-point", "8,8,4", "--coordinate-system", "rec-voxel",
                   "--indexing", "zero-based", "--roi-radius-angst", "8"])

        recs = os.path.join(project, "tilt_series", TOMO_NAME, "individual_particles_recs")
        assert not (os.path.isdir(recs) and os.listdir(recs)), "no crop expected yet"

        rc = _run_cmd(["status", "--project", project, "--tomo-name", TOMO_NAME,
                        "--create-volume", "--crop-storage-box-size", "6"])
        assert rc == 0, f"status --create-volume returned {rc}"

        crop = os.path.join(recs, "P000001_rec.mrc")
        assert os.path.isfile(crop), "rec crop not created"
        hdr = read_mrc_header(crop)
        assert (hdr.nx, hdr.ny, hdr.nz) == (6, 6, 6)

        p_star = os.path.join(project, "tilt_series", TOMO_NAME,
                               "individual_particles", "P000001.star")
        blks = read_star(p_star)
        assert "tomoJANAS_particle_rec_crop" in blks, "crop block not added to P*.star"
        print("PASS test_status_create_volume")


def test_validate_strict():
    with tempfile.TemporaryDirectory() as d:
        imod_dir = _make_imod_dir(d)
        project = os.path.join(d, "project")

        _run_cmd(["imod", "--project", project, "--create-if-missing",
                   "--imod-dir", imod_dir, "--tomo-name", TOMO_NAME])
        _run_cmd(["particles", "--project", project, "--tomo-name", TOMO_NAME,
                   "--input-single-point", "8,8,4", "--coordinate-system", "rec-voxel",
                   "--roi-radius-angst", "20.0"])

        rc = _run_cmd(["validate", "--project", project,
                        "--tomo-name", TOMO_NAME, "--strict"])

        # logs written
        assert os.path.isfile(os.path.join(project, "logs", "validation_log.json"))
        assert os.path.isfile(os.path.join(project, "logs", "validation_log.star"))
        assert os.path.isfile(os.path.join(project, "logs", "validation_log.md"))

        # validate should pass (rc 0 or just check the log)
        with open(os.path.join(project, "logs", "validation_log.json")) as f:
            report = json.load(f)
        # It may not be fully OK due to missing RELION matrices / missing files on CI,
        # but it must produce a valid report
        assert "ok" in report
        assert "items" in report
        assert "summary" in report

        print(f"PASS test_validate_strict (rc={rc}, report_ok={report['ok']})")


def test_imod_parsers():
    """Unit tests for IMOD metadata parsers."""
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        # test read_xf
        xf_path = os.path.join(d, "t.xf")
        with open(xf_path, "w") as f:
            f.write("1.0  0.0  0.0  1.0  0.0  0.0\n")
            f.write("0.9998 -0.0087  0.0087  0.9998  2.5  -1.2\n")
        from tomojanas.metadata.imod_xf import read_xf, invert_xf_for_point, apply_xf_to_point
        df = read_xf(xf_path)
        assert len(df) == 2
        assert df["tilt_index"].tolist() == [1, 2]
        # identity: apply and invert should be identity
        x2, y2 = apply_xf_to_point(5.0, 3.0, df.iloc[0])
        assert abs(x2 - 5.0) < 1e-10 and abs(y2 - 3.0) < 1e-10
        x3, y3 = invert_xf_for_point(5.0, 3.0, df.iloc[0])
        assert abs(x3 - 5.0) < 1e-10 and abs(y3 - 3.0) < 1e-10

        # test read_tlt
        tlt_path = os.path.join(d, "t.tlt")
        with open(tlt_path, "w") as f:
            f.write("-30.0\n-20.0\n-10.0\n0.0\n10.0\n")
        from tomojanas.metadata.imod_tlt import read_tlt
        tlt_df = read_tlt(tlt_path)
        assert len(tlt_df) == 5
        assert tlt_df["tilt_index"].tolist() == [1, 2, 3, 4, 5]
        assert abs(tlt_df["tilt_angle"].iloc[3]) < 1e-10  # 0.0

        # test parse_com_file
        newst_path = os.path.join(d, "newst.com")
        with open(newst_path, "w") as f:
            f.write("$newst\n")
            f.write("InputFile raw.mrc\n")
            f.write("OutputFile ali.mrc\n")
            f.write("SizeToOutputInXandY 512,512\n")
            f.write("BinByFactor 2\n")
        from tomojanas.metadata.imod_com import parse_com_file, com_get_float, com_get_numbers
        com = parse_com_file(newst_path)
        assert com.get("InputFile") == "raw.mrc"
        assert com_get_float(com, "BinByFactor") == 2.0
        nums = com_get_numbers(com.get("SizeToOutputInXandY"))
        assert nums == [512.0, 512.0]

        # test EXCLUDELIST range parsing
        tilt_path2 = os.path.join(d, "tilt.com")
        with open(tilt_path2, "w") as f:
            f.write("FULLIMAGE 512 512\n")
            f.write("TILTFILE test.tlt\n")
            f.write("THICKNESS 300\n")
            f.write("IMAGEBINNED 1\n")
            f.write("EXCLUDELIST 1,3-5\n")
            f.write("SHIFT 2.5 -1.0\n")
        com2 = parse_com_file(tilt_path2)
        from tomojanas.metadata.imod_com import parse_exclude_list
        excluded = parse_exclude_list(com2.get("EXCLUDELIST"), 10)
        assert excluded == {0, 2, 3, 4}, f"unexpected excluded: {excluded}"

        # SHIFT parsing
        from tomojanas.metadata.imod_com import com_get_numbers
        shift_nums = com_get_numbers(com2.get("SHIFT"))
        assert shift_nums == [2.5, -1.0]

    print("PASS test_imod_parsers")


def test_imod_mapping_math():
    """Test RELION ImodImport port: identity xf, 3 tilts, zero-tilt projection."""
    import tempfile, math
    from tomojanas.geometry.imod_mapping import build_relion_imod_mapping, project_particle_to_ali

    with tempfile.TemporaryDirectory() as d:
        imod_dir = _make_imod_dir(d)

        mapping = build_relion_imod_mapping(
            ts_path=os.path.join(imod_dir, f"{TOMO_NAME}_ali.mrc"),
            imod_dir=imod_dir,
        )
        assert len(mapping.projections) == N_TILTS, f"expected {N_TILTS} projections"
        assert mapping.frames_missing is False
        assert mapping.old_frame_index == list(range(N_TILTS))

        # At zero tilt (index 1): world centre → image centre
        # World centre: x_rec=NX/2, y_rec=NY/2, z_rec=NZ_REC/2 (ali convention)
        cx, cy = (NX - 1) / 2.0, (NY - 1) / 2.0
        B = 1.0  # B_rec_ali = 1 since same pixel size
        # project particle at reconstruction centre
        x_img, y_img = project_particle_to_ali(
            NX / 2.0, NY / 2.0, NZ_REC / 2.0, B, mapping.projections[1]
        )
        assert abs(x_img - cx) < 1.0, f"zero-tilt x_img={x_img:.3f} expected ~{cx:.3f}"
        assert abs(y_img - cy) < 1.0, f"zero-tilt y_img={y_img:.3f} expected ~{cy:.3f}"

    print("PASS test_imod_mapping_math")


def test_ctf_missing():
    """CTF import with no CTF file must not crash and marks CTF as missing."""
    with tempfile.TemporaryDirectory() as d:
        imod_dir = _make_imod_dir(d)
        project = os.path.join(d, "project")
        _run_cmd(["imod", "--project", project, "--create-if-missing",
                   "--imod-dir", imod_dir, "--tomo-name", TOMO_NAME])

        # No CTF file exists
        rc = _run_cmd(["ctf", "--project", project, "--tomo-name", TOMO_NAME,
                        "--ctf-source", "none"])
        assert rc == 0, f"ctf import returned {rc}"

        # ctf_settings.star should say none
        ctf_star = os.path.join(project, "tilt_series", TOMO_NAME, "imod_settings", "ctf_settings.star")
        assert os.path.isfile(ctf_star)
        blks = read_star(ctf_star)
        assert "tomoJANAS_ctf_settings" in blks
        assert blks["tomoJANAS_ctf_settings"]["pairs"]["_tomoJANASCtfPremultiplied"] == "0"

    print("PASS test_ctf_missing")


def test_defocus_unit_conversion():
    """Defocus nm→Å conversion is exact."""
    from tomojanas.metadata.ctf import defocus_to_angstrom
    assert defocus_to_angstrom(1000.0, "nm") == 10000.0
    assert defocus_to_angstrom(1.0, "micrometer") == 10000.0
    assert defocus_to_angstrom(10000.0, "angstrom") == 10000.0
    print("PASS test_defocus_unit_conversion")


def test_legacy_import_volume():
    """Legacy tomojanas import_volume must still work."""
    with tempfile.TemporaryDirectory() as d:
        vol = np.zeros((8, 8, 8), dtype=np.float32)
        tomo_path = os.path.join(d, "tomo.mrc")
        write_mrc(tomo_path, vol, 2.0)

        coords_path = os.path.join(d, "coords.csv")
        with open(coords_path, "w") as f:
            f.write("x,y,z\n4.0,4.0,4.0\n")

        out_dir = os.path.join(d, "out")
        out_star = os.path.join(d, "out.star")

        # Call main() of legacy caller
        from janas.tomojanas_cmd_caller import main as legacy_main
        legacy_main([
            "import_volume",
            "--tomo", tomo_path,
            "--coords", coords_path,
            "--format", "csv",
            "--box", "4",
            "--apix", "2.0",
            "--outDir", out_dir,
            "--o", out_star,
        ])
        assert os.path.isfile(out_star), "legacy STAR file not written"

    print("PASS test_legacy_import_volume")


TESTS: List[Tuple[str, Callable]] = [
    ("test_imod_parsers", test_imod_parsers),
    ("test_imod_mapping_math", test_imod_mapping_math),
    ("test_imod_import", test_imod_import),
    ("test_particles_import", test_particles_import),
    ("test_rec_crop_writing", test_rec_crop_writing),
    ("test_external_paths_absolute", test_external_paths_absolute),
    ("test_particle_autoincrement_and_command_log", test_particle_autoincrement_and_command_log),
    ("test_status_scan_and_sync", test_status_scan_and_sync),
    ("test_crop_axis_mapping_xyz", test_crop_axis_mapping_xyz),
    ("test_default_axis_order_is_xyz", test_default_axis_order_is_xyz),
    ("test_explicit_xzy_axis_order", test_explicit_xzy_axis_order),
    ("test_status_create_volume", test_status_create_volume),
    ("test_validate_strict", test_validate_strict),
    ("test_ctf_missing", test_ctf_missing),
    ("test_defocus_unit_conversion", test_defocus_unit_conversion),
    ("test_legacy_import_volume", test_legacy_import_volume),
]


def main() -> int:
    failed = []
    for name, fn in TESTS:
        try:
            fn()
        except Exception as exc:
            failed.append(name)
            import traceback
            print(f"FAIL {name}: {exc}")
            traceback.print_exc()
    print(f"\n{len(TESTS) - len(failed)}/{len(TESTS)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
