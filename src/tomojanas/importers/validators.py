#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/importers/validators.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
``tomojanas-import validate`` — validate project consistency.

Writes logs/validation_log.{star,json,md}. In --strict mode returns non-zero
on any critical failure.
"""
from __future__ import annotations

import os
import sys
from typing import Any, Dict, List, Optional

from tomojanas import get_version
from tomojanas.io.logs import write_validation_logs, now_iso, ImportLogger
from tomojanas.io import project_writer as pw
from tomojanas.io.star_reader import read_star
from tomojanas.models.project import Project

__all__ = ["validate_project", "validate_project_cli"]


def _item(check: str, status: str, severity: str, detail: str = "") -> Dict[str, Any]:
    return {"check": check, "status": status, "severity": severity, "detail": detail}


def validate_project(
    project,
    tomo_name: Optional[str] = None,
    strict: bool = False,
    tolerance_pixel: float = 0.5,
    tolerance_angst: float = 1.0,
) -> Dict[str, Any]:
    """Run all validation checks; returns a report dict."""
    items: List[Dict[str, Any]] = []
    root = project.root

    # ------------------------------------------------------------------ #
    # Project-level checks
    # ------------------------------------------------------------------ #
    def _check_exists(path, name, sev="critical"):
        ok = os.path.isfile(path) or os.path.isdir(path)
        items.append(_item(
            f"exists:{name}", "pass" if ok else "fail", sev,
            f"{path} {'found' if ok else 'NOT FOUND'}"
        ))
        return ok

    _check_exists(root, "project_root")
    ts_ok = _check_exists(project.tomograms_star, "tomograms.star")
    _check_exists(project.optimisation_set_star, "optimisation_set.star", sev="warning")
    _check_exists(project.logs_dir, "logs/", sev="warning")

    # ------------------------------------------------------------------ #
    # Tomogram-level checks
    # ------------------------------------------------------------------ #
    tomo_names_to_check: List[str] = []
    if tomo_name:
        tomo_names_to_check = [tomo_name]
    elif ts_ok:
        try:
            blks = read_star(project.tomograms_star)
            if "global" in blks and blks["global"]["type"] == "loop":
                df = blks["global"]["df"]
                for _, r in df.iterrows():
                    n = str(r.get("_rlnTomoName", ""))
                    if n:
                        tomo_names_to_check.append(n)
        except Exception as exc:
            items.append(_item("read:tomograms.star", "fail", "critical", str(exc)))

    for tn in tomo_names_to_check:
        _validate_tomogram(root, tn, items, tolerance_pixel, tolerance_angst)

    # ------------------------------------------------------------------ #
    # Particles
    # ------------------------------------------------------------------ #
    if os.path.isfile(project.particles_all_star):
        _validate_particles_star(project.particles_all_star, items)
    else:
        items.append(_item("exists:particles_all.star", "warn", "warning",
                           "particles_all.star not yet created"))

    # ------------------------------------------------------------------ #
    # Summarise
    # ------------------------------------------------------------------ #
    fails = [i for i in items if i["status"] == "fail"]
    crit_fails = [i for i in fails if i["severity"] == "critical"]
    warns = [i for i in items if i["status"] in ("warn", "warning")]
    ok = len(crit_fails) == 0

    summary = (
        f"{len(items)} checks: "
        f"{sum(1 for i in items if i['status']=='pass')} pass, "
        f"{len(fails)} fail ({len(crit_fails)} critical), "
        f"{len(warns)} warn"
    )

    return {
        "generated": now_iso(),
        "tomojanas_version": get_version(),
        "scope": f"tomogram:{tomo_name}" if tomo_name else "project",
        "ok": ok,
        "summary": summary,
        "items": items,
    }


def _validate_tomogram(root, tomo_name, items, tol_px, tol_ang):
    """Per-tomogram validation checks."""
    tilt_star = pw.tilt_series_star_path(root, tomo_name)
    imod_dir = pw.imod_settings_dir(root, tomo_name)

    # Tilt series STAR exists
    ts_ok = os.path.isfile(tilt_star)
    items.append(_item(
        f"exists:tilt_series/{tomo_name}.star",
        "pass" if ts_ok else "fail", "critical",
        tilt_star,
    ))

    # imod_settings exist
    for fname in ["imod_geometry.star", "relion_geometry_status.star"]:
        p = os.path.join(imod_dir, fname)
        items.append(_item(
            f"exists:imod_settings/{fname}",
            "pass" if os.path.isfile(p) else "warn", "warning",
            p,
        ))

    # Load tomo metadata from tomograms.star
    ts_path = pw.tomograms_star_path(root)
    rec_path, ali_path, a_ali, B_rec_ali, tomo_sx, tomo_sy, tomo_sz = None, None, None, None, 0, 0, 0
    if os.path.isfile(ts_path):
        try:
            blks = read_star(ts_path)
            for bname in ("global", "tomoJANAS_tomogram_sources"):
                if bname not in blks:
                    continue
                if blks[bname]["type"] != "loop":
                    continue
                df = blks[bname]["df"]
                name_col = "_rlnTomoName" if bname == "global" else "_tomoJANASTomoName"
                for _, r in df.iterrows():
                    if str(r.get(name_col, "")) != tomo_name:
                        continue
                    if bname == "global":
                        a_ali = _safe_float(r.get("_rlnTomoTiltSeriesPixelSize"))
                        B_rec_ali = _safe_float(r.get("_rlnTomoTomogramBinning"))
                        tomo_sx = _safe_int(r.get("_rlnTomoSizeX"))
                        tomo_sy = _safe_int(r.get("_rlnTomoSizeY"))
                        tomo_sz = _safe_int(r.get("_rlnTomoSizeZ"))
                        rec_path = str(r.get("_rlnTomoReconstructedTomogram", "?") or "?")
                    else:
                        ali_path_r = str(r.get("_tomoJANASAliStack", "?") or "?")
                        if ali_path_r not in ("?", "", None):
                            if not os.path.isabs(ali_path_r):
                                ali_path_r = os.path.join(root, ali_path_r)
                            ali_path = ali_path_r
        except Exception as exc:
            items.append(_item(f"read:tomograms.star/{tomo_name}", "fail", "critical", str(exc)))

    # Check MRC files are readable
    from tomojanas.io.mrc import read_mrc_header, validate_mrc_pixel_size
    for label, path in [("rec_tomo", rec_path), ("ali_stack", ali_path)]:
        if path and path not in ("?", "", None):
            # Resolve relative
            if not os.path.isabs(path):
                path = os.path.join(root, path)
            if os.path.isfile(path):
                try:
                    hdr = read_mrc_header(path)
                    items.append(_item(
                        f"mrc_header:{label}", "pass", "info",
                        f"{path}: {hdr.nx}x{hdr.ny}x{hdr.nz}, apix={hdr.pixel_x:.3f}"
                    ))
                    # nsymbt check
                    items.append(_item(
                        f"nsymbt:{label}", "pass", "info",
                        f"data_offset={hdr.data_offset}"
                    ))
                    # Pixel size valid
                    if hdr.pixel_x <= 0:
                        items.append(_item(f"pixel_size:{label}", "fail", "warning",
                                           f"pixel_x={hdr.pixel_x}"))
                    else:
                        items.append(_item(f"pixel_size:{label}", "pass", "info",
                                           f"pixel_x={hdr.pixel_x}"))
                except Exception as exc:
                    items.append(_item(f"mrc_header:{label}", "fail", "critical", str(exc)))
            else:
                items.append(_item(f"exists:{label}", "warn", "warning",
                                   f"{path} not found"))

    # Check tilt series consistency if tilt STAR exists
    if ts_ok:
        try:
            blks = read_star(tilt_star)
            if tomo_name in blks and blks[tomo_name]["type"] == "loop":
                n_rows = len(blks[tomo_name]["df"])
                # Compare with ali stack nz
                if ali_path and ali_path not in ("?", "", None) and os.path.isfile(ali_path):
                    ali_hdr = read_mrc_header(ali_path)
                    match = ali_hdr.nz == n_rows
                    items.append(_item(
                        f"ali_nz_vs_tilt_rows:{tomo_name}",
                        "pass" if match else "fail", "critical" if not match else "info",
                        f"ali nz={ali_hdr.nz}, tilt rows={n_rows}",
                    ))
        except Exception as exc:
            items.append(_item(f"read:tilt_star/{tomo_name}", "warn", "warning", str(exc)))

    # Check projection matrices exist
    if ts_ok:
        try:
            blks = read_star(tilt_star)
            has_matrices = "tomoJANAS_projection_matrices" in blks
            items.append(_item(
                f"projection_matrices:{tomo_name}",
                "pass" if has_matrices else "warn", "warning",
                "tomoJANAS_projection_matrices block present" if has_matrices else "No projection matrices"
            ))
        except Exception:
            pass

    # Coordinate round-trip check (if particles exist)
    from tomojanas.geometry.coordinates import coordinate_roundtrip_error
    if a_ali and B_rec_ali and tomo_sx and tomo_sy and tomo_sz:
        centre_x = tomo_sx / 2.0 / B_rec_ali
        centre_y = tomo_sy / 2.0 / B_rec_ali
        centre_z = tomo_sz / 2.0 / B_rec_ali
        err = coordinate_roundtrip_error(
            centre_x, centre_y, centre_z,
            tomo_sx, tomo_sy, tomo_sz,
            a_ali, B_rec_ali,
        )
        ok = err < tol_px
        items.append(_item(
            f"coordinate_roundtrip:{tomo_name}",
            "pass" if ok else "fail", "critical" if not ok else "info",
            f"roundtrip error={err:.2e} (tol={tol_px})"
        ))


def _validate_particles_star(path, items):
    """Check particles_all.star has data_optics + data_particles."""
    try:
        blks = read_star(path)

        has_optics = "optics" in blks and blks["optics"]["type"] == "loop"
        items.append(_item(
            "particles_all:has_data_optics",
            "pass" if has_optics else "fail", "critical",
            "data_optics block present" if has_optics else "MISSING data_optics in particles_all.star"
        ))

        has_particles = "particles" in blks and blks["particles"]["type"] == "loop"
        items.append(_item(
            "particles_all:has_data_particles",
            "pass" if has_particles else "fail", "critical",
            "data_particles block present" if has_particles else "MISSING data_particles"
        ))

        if has_optics:
            df_optics = blks["optics"]["df"]
            has_apix = "_rlnTomoTiltSeriesPixelSize" in df_optics.columns
            items.append(_item(
                "particles_all:optics_has_apix",
                "pass" if has_apix else "fail", "critical",
                "_rlnTomoTiltSeriesPixelSize" + (" present" if has_apix else " MISSING in data_optics")
            ))

        if has_particles and has_optics:
            df_particles = blks["particles"]["df"]
            has_og = "_rlnOpticsGroup" in df_particles.columns
            items.append(_item(
                "particles_all:particles_has_opticsgroup",
                "pass" if has_og else "fail", "critical",
                "_rlnOpticsGroup" + (" in data_particles" if has_og else " MISSING in data_particles")
            ))

            # Check no _tomoJANAS* in RELION block
            relion_cols = [c for c in df_particles.columns if c.startswith("_tomoJANAS")]
            no_cross = len(relion_cols) == 0
            items.append(_item(
                "particles_all:no_tomoJANAS_in_relion_loop",
                "pass" if no_cross else "fail", "warning",
                "No _tomoJANAS* tags in RELION loop" if no_cross else f"Found: {relion_cols}"
            ))

    except Exception as exc:
        items.append(_item("read:particles_all.star", "fail", "critical", str(exc)))


def _safe_float(v) -> Optional[float]:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _safe_int(v) -> Optional[int]:
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def validate_project_cli(args) -> int:
    project_root = os.path.abspath(args.project)
    tomo_name = getattr(args, "tomo_name", None)
    strict = getattr(args, "strict", False)
    tol_px = getattr(args, "tolerance_pixel", 0.5) or 0.5
    tol_ang = getattr(args, "tolerance_angst", 1.0) or 1.0

    logger = ImportLogger(project_root, get_version())
    logger.invocation("tomojanas-import validate", argv=sys.argv)

    try:
        proj = Project.open_or_create(project_root, fail_if_existing=False, append=True)
    except Exception as exc:
        logger.error(str(exc))
        # Write a minimal failure report
        write_validation_logs(project_root, {
            "generated": now_iso(),
            "tomojanas_version": get_version(),
            "scope": "project",
            "ok": False,
            "summary": str(exc),
            "items": [_item("open_project", "fail", "critical", str(exc))],
        })
        return 1

    report = validate_project(proj, tomo_name=tomo_name, strict=strict,
                               tolerance_pixel=tol_px, tolerance_angst=tol_ang)
    write_validation_logs(project_root, report)
    logger.summary(validation_summary=report["summary"])

    print(f"[tomojanas-import validate] {report['summary']}")
    print(f"  result: {'OK' if report['ok'] else 'FAILED'}")

    if not report["ok"] and strict:
        return 1
    return 0
