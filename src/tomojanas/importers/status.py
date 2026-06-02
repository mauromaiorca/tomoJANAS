#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/importers/status.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
``tomojanas-import status`` — scan a project tree and report what particles
are registered vs. what exists on disk, and reconcile after manual edits.

The recommended way to delete a particle is to remove its files by hand from
the OS, e.g.::

    rm tilt_series/<tomo>/individual_particles/P000003.star
    rm tilt_series/<tomo>/individual_particles_recs/P000003_rec.mrc

then run ``tomojanas-import status --project <p> --sync`` to rebuild
``particles_all.star`` from the P*.star files that remain.
"""

from __future__ import annotations

import glob
import os
import re
from typing import Dict, List, Optional

from tomojanas.models.project import Project
from tomojanas.io import project_writer as pw
from tomojanas.io.star_reader import read_star
from tomojanas.io.star_writer import LoopBlock, PairBlock, write_star
from tomojanas.metadata.relion_labels import (
    PARTICLE_OPTICS_COLUMNS, PARTICLES_COLUMNS, PARTICLE_REC_CROP_COLUMNS,
)

_P_STAR_RE = re.compile(r"^(.*)\.star$")


# ------------------------------------------------------------------ #
# helpers
# ------------------------------------------------------------------ #
def _list_tomograms(project: Project) -> List[str]:
    tomos: List[str] = []
    if os.path.isfile(project.tomograms_star):
        try:
            blks = read_star(project.tomograms_star)
            if "global" in blks and blks["global"]["type"] == "loop":
                df = blks["global"]["df"]
                if "_rlnTomoName" in df.columns:
                    tomos = [str(v) for v in df["_rlnTomoName"].tolist()]
        except Exception:
            pass
    # also include any tilt_series/<tomo>/ dirs present on disk
    ts_dir = project.tilt_series_dir
    if os.path.isdir(ts_dir):
        for name in sorted(os.listdir(ts_dir)):
            d = os.path.join(ts_dir, name)
            if os.path.isdir(d) and name not in tomos:
                tomos.append(name)
    return tomos


def _registered_particles(project: Project) -> Dict[str, set]:
    """Map tomo_name -> set of particle names registered in particles_all.star."""
    out: Dict[str, set] = {}
    if not os.path.isfile(project.particles_all_star):
        return out
    try:
        blks = read_star(project.particles_all_star)
        if "particles" in blks and blks["particles"]["type"] == "loop":
            df = blks["particles"]["df"]
            if "_rlnTomoName" in df.columns and "_rlnTomoParticleName" in df.columns:
                for _, r in df.iterrows():
                    out.setdefault(str(r["_rlnTomoName"]), set()).add(
                        str(r["_rlnTomoParticleName"]))
    except Exception:
        pass
    return out


def _disk_names(directory: str, suffix: str) -> set:
    """Return the set of particle names from files in *directory* ending with
    *suffix* (e.g. '.star' -> P000001, '_rec.mrc' -> P000001)."""
    names = set()
    if not os.path.isdir(directory):
        return names
    for path in glob.glob(os.path.join(directory, f"*{suffix}")):
        base = os.path.basename(path)
        if base.endswith(suffix):
            names.add(base[: -len(suffix)])
    return names


def scan_tomogram(project: Project, tomo: str, registry: Dict[str, set]) -> Dict:
    """Scan one tomogram; return a structured report dict."""
    ip = project.individual_particles_dir(tomo)
    recs = os.path.join(project.tomogram_dir(tomo), "individual_particles_recs")
    ali = os.path.join(project.tomogram_dir(tomo), "individual_particles_ali")
    raw = os.path.join(project.tomogram_dir(tomo), "individual_particles_raw")

    registered = registry.get(tomo, set())
    star_files = _disk_names(ip, ".star")
    rec_crops = _disk_names(recs, "_rec.mrc")
    ali_crops = _disk_names(ali, "_ali.mrcs")
    raw_crops = _disk_names(raw, "_raw.mrcs")

    items: List[Dict] = []
    for name in sorted(registered - star_files):
        items.append({"particle": name, "issue": "registered_but_no_star",
                      "detail": "in particles_all.star but P*.star missing on disk"})
    for name in sorted(star_files - registered):
        items.append({"particle": name, "issue": "star_not_registered",
                      "detail": "P*.star on disk but not in particles_all.star (run --sync)"})
    for name in sorted(rec_crops - star_files):
        items.append({"particle": name, "issue": "orphan_rec_crop",
                      "detail": "_rec.mrc on disk but no P*.star"})
    for name in sorted(ali_crops - star_files):
        items.append({"particle": name, "issue": "orphan_ali_crop",
                      "detail": "_ali.mrcs on disk but no P*.star"})
    for name in sorted(raw_crops - star_files):
        items.append({"particle": name, "issue": "orphan_raw_crop",
                      "detail": "_raw.mrcs on disk but no P*.star"})

    return {
        "tomo": tomo,
        "registered": sorted(registered),
        "star_files": sorted(star_files),
        "rec_crops": sorted(rec_crops),
        "ali_crops": sorted(ali_crops),
        "raw_crops": sorted(raw_crops),
        "n_registered": len(registered),
        "n_star": len(star_files),
        "n_rec_crops": len(rec_crops),
        "issues": items,
        "ok": len(items) == 0,
    }


def scan_project(project: Project, tomo_name: Optional[str] = None) -> Dict:
    registry = _registered_particles(project)
    tomos = [tomo_name] if tomo_name else _list_tomograms(project)
    reports = [scan_tomogram(project, t, registry) for t in tomos]
    n_issues = sum(len(r["issues"]) for r in reports)
    return {"tomograms": reports, "n_issues": n_issues, "ok": n_issues == 0}


# ------------------------------------------------------------------ #
# sync: rebuild particles_all.star from on-disk P*.star
# ------------------------------------------------------------------ #
def sync_particles_all(project: Project, tomo_names: List[str]) -> int:
    """Rebuild particles_all.star from the P*.star files present on disk.

    Returns the number of particle rows written.
    """
    optics_by_group: Dict[str, Dict] = {}
    particle_rows: List[Dict] = []

    for tomo in tomo_names:
        ip = project.individual_particles_dir(tomo)
        if not os.path.isdir(ip):
            continue
        for path in sorted(glob.glob(os.path.join(ip, "*.star"))):
            try:
                blks = read_star(path)
            except Exception:
                continue
            if "optics" in blks and blks["optics"]["type"] == "loop":
                odf = blks["optics"]["df"]
                for _, r in odf.iterrows():
                    grp = str(r.get("_rlnOpticsGroup", "1"))
                    optics_by_group.setdefault(grp, {c: r.get(c) for c in odf.columns})
            if "particles" in blks and blks["particles"]["type"] == "loop":
                pdf = blks["particles"]["df"]
                for _, r in pdf.iterrows():
                    particle_rows.append({c: r.get(c) for c in pdf.columns})

    optics_rows = list(optics_by_group.values()) or [{
        "_rlnOpticsGroup": 1, "_rlnOpticsGroupName": "opticsGroup1",
    }]
    blocks = [
        LoopBlock(name="optics", columns=PARTICLE_OPTICS_COLUMNS, rows=optics_rows),
        LoopBlock(name="particles", columns=PARTICLES_COLUMNS, rows=particle_rows),
    ]
    write_star(project.particles_all_star, blocks)
    return len(particle_rows)


# ------------------------------------------------------------------ #
# create rec sub-volumes for already-imported particles
# ------------------------------------------------------------------ #
def _reload_blocks(path):
    """Read a STAR file and return it as writer block objects (round-trip)."""
    blks = read_star(path)
    out = []
    for name, b in blks.items():
        if b["type"] == "loop":
            df = b["df"]
            out.append(LoopBlock(name=name, columns=list(df.columns),
                                 rows=df.values.tolist()))
        else:
            out.append(PairBlock(name=name, pairs=dict(b["pairs"])))
    return out


def _tomo_geometry(project: Project, tomo: str):
    """Return (rec_path, a_ali, B_rec_ali, sx, sy, sz) from tomograms.star, or None."""
    if not os.path.isfile(project.tomograms_star):
        return None
    blks = read_star(project.tomograms_star)
    if "global" not in blks or blks["global"]["type"] != "loop":
        return None
    df = blks["global"]["df"]
    row = df[df["_rlnTomoName"].astype(str) == str(tomo)]
    if row.empty:
        return None
    r = row.iloc[0]

    def _f(col, default=0.0):
        try:
            return float(r.get(col, default))
        except (TypeError, ValueError):
            return default

    rec_path = str(r.get("_rlnTomoReconstructedTomogram", "?"))
    return (
        rec_path,
        _f("_rlnTomoTiltSeriesPixelSize", 1.0),
        _f("_rlnTomoTomogramBinning", 1.0),
        _f("_rlnTomoSizeX"), _f("_rlnTomoSizeY"), _f("_rlnTomoSizeZ"),
    )


def create_rec_crops(project: Project, tomo_names: List[str], args, logger) -> int:
    """Create the 3D rec sub-volume (_rec.mrc) for particles that lack one.

    Coordinates and ROI radius are read from each P*.star; the rec-voxel
    centre is derived from the canonical RELION centered-Angstrom coordinate,
    so it is independent of the original picking convention.
    """
    from tomojanas.io.mrc import read_mrc_data, read_mrc_header
    from tomojanas.geometry.coordinates import relion_centered_angst_to_rec_voxel
    from tomojanas.importers.particle_importer import _write_rec_crop, _resolve_rec_path

    created = 0
    box_override = getattr(args, "crop_storage_box_size", None)
    pad_vox = float(getattr(args, "crop_padding_voxel", 0.0) or 0.0)
    pad_angst = float(getattr(args, "crop_padding_angst", 0.0) or 0.0)
    outside_policy = getattr(args, "crop_outside_policy", "partial") or "partial"
    pad_value = float(getattr(args, "crop_pad_value", 0.0) or 0.0)
    apply_mask = getattr(args, "apply_spherical_mask", False)
    overwrite = getattr(args, "overwrite_crops", False)

    for tomo in tomo_names:
        geom = _tomo_geometry(project, tomo)
        if geom is None:
            logger.warning(f"{tomo}: no tomogram geometry in tomograms.star; skipping crops")
            continue
        rec_path_str, a_ali, B_rec_ali, sx, sy, sz = geom
        rec_path = _resolve_rec_path(rec_path_str, project.root)
        if not rec_path:
            logger.warning(f"{tomo}: rec tomogram not found ({rec_path_str}); skipping crops")
            continue
        try:
            rec_volume, rec_hdr = read_mrc_data(rec_path)
        except Exception as exc:
            logger.warning(f"{tomo}: cannot read rec tomogram: {exc}")
            continue
        a_rec = rec_hdr.pixel_x if rec_hdr.pixel_x > 0 else (a_ali * B_rec_ali)
        crop_dir = os.path.join(project.tomogram_dir(tomo), "individual_particles_recs")

        ip = project.individual_particles_dir(tomo)
        for p_star in sorted(glob.glob(os.path.join(ip, "*.star"))):
            pname = os.path.splitext(os.path.basename(p_star))[0]
            out_mrc = os.path.join(crop_dir, f"{pname}_rec.mrc")
            if os.path.isfile(out_mrc) and not overwrite:
                continue
            try:
                blks = read_star(p_star)
                pdf = blks["particles"]["df"].iloc[0]
                x_ang = float(pdf["_rlnCenteredCoordinateXAngst"])
                y_ang = float(pdf["_rlnCenteredCoordinateYAngst"])
                z_ang = float(pdf["_rlnCenteredCoordinateZAngst"])
                radius_angst = 0.0
                if "tomoJANAS_particle_roi" in blks:
                    rdf = blks["tomoJANAS_particle_roi"]["df"].iloc[0]
                    radius_angst = float(rdf.get("_tomoJANASRoiRadiusAngst", 0.0))
            except Exception as exc:
                logger.warning(f"{pname}: cannot read coordinates/ROI: {exc}")
                continue

            x_rec, y_rec, z_rec = relion_centered_angst_to_rec_voxel(
                x_ang, y_ang, z_ang, sx, sy, sz, a_ali, B_rec_ali, "zero-based")

            crop_row = _write_rec_crop(
                rec_volume, rec_hdr, crop_dir, pname,
                x_rec, y_rec, z_rec, "zero-based",
                radius_angst, a_rec,
                box_override, pad_vox, pad_angst,
                outside_policy, pad_value, apply_mask,
                project.root, rec_path, logger,
            )
            if crop_row is None:
                continue
            # add/refresh the rec-crop block in P*.star
            blocks = [b for b in _reload_blocks(p_star)
                      if getattr(b, "name", None) != "tomoJANAS_particle_rec_crop"]
            blocks.append(LoopBlock(name="tomoJANAS_particle_rec_crop",
                                    columns=PARTICLE_REC_CROP_COLUMNS, rows=[crop_row]))
            write_star(p_star, blocks)
            created += 1

    return created


# ------------------------------------------------------------------ #
# CLI
# ------------------------------------------------------------------ #
def status_cli(args) -> int:
    project_root = os.path.abspath(args.project)
    if not os.path.isdir(project_root):
        print(f"[tomojanas-import status] ERROR: project '{args.project}' not found")
        return 1
    project = Project(root=project_root)

    tomo_name = getattr(args, "tomo_name", None)

    # optionally create the 3D rec sub-volumes for already-imported particles
    if getattr(args, "create_volume", False):
        from tomojanas.io.logs import ImportLogger
        from tomojanas import get_version
        logger = ImportLogger(project_root, version=get_version())
        tomos = [tomo_name] if tomo_name else _list_tomograms(project)
        n = create_rec_crops(project, tomos, args, logger)
        logger.info(f"created {n} rec sub-volume(s)")
        print(f"[create-volume] created {n} rec sub-volume(s) (_rec.mrc)")

    report = scan_project(project, tomo_name=tomo_name)

    # human-readable report
    for r in report["tomograms"]:
        print(f"\n== tomogram: {r['tomo']} ==")
        print(f"   registered in particles_all.star : {r['n_registered']}")
        print(f"   individual P*.star on disk        : {r['n_star']}")
        print(f"   rec crops (_rec.mrc)              : {r['n_rec_crops']}")
        if r["ali_crops"]:
            print(f"   ali crops (_ali.mrcs)             : {len(r['ali_crops'])}")
        if r["raw_crops"]:
            print(f"   raw crops (_raw.mrcs)             : {len(r['raw_crops'])}")
        if r["issues"]:
            print(f"   ISSUES ({len(r['issues'])}):")
            for it in r["issues"]:
                print(f"     - {it['particle']}: {it['issue']} ({it['detail']})")
        else:
            print("   no inconsistencies")

    if getattr(args, "sync", False):
        tomos = [tomo_name] if tomo_name else _list_tomograms(project)
        n = sync_particles_all(project, tomos)
        print(f"\n[sync] rebuilt particles_all.star from on-disk P*.star: {n} particle(s)")
        # log the sync as a command too (handled by cli.main)
        report = scan_project(project, tomo_name=tomo_name)  # re-scan after sync

    print(f"\nTotal inconsistencies: {report['n_issues']}")
    if report["n_issues"] and getattr(args, "strict", False):
        return 1
    return 0
