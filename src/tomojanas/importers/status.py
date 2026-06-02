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
from tomojanas.io.star_writer import LoopBlock, write_star
from tomojanas.metadata.relion_labels import (
    PARTICLE_OPTICS_COLUMNS, PARTICLES_COLUMNS,
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
# CLI
# ------------------------------------------------------------------ #
def status_cli(args) -> int:
    project_root = os.path.abspath(args.project)
    if not os.path.isdir(project_root):
        print(f"[tomojanas-import status] ERROR: project '{args.project}' not found")
        return 1
    project = Project(root=project_root)

    tomo_name = getattr(args, "tomo_name", None)
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
