#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/io/project_writer.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
tomoJANAS project layout: directory creation, canonical paths, and manifest.

A tomoJANAS project is a self-contained directory tree (NOT a relion/ clone):

    <project_root>/
    ├── tomograms.star
    ├── particles_all.star
    ├── optimisation_set.star
    ├── project_manifest.star
    ├── project_manifest.json
    ├── tilt_series/
    │   ├── <tomo>.star
    │   └── <tomo>/
    │       ├── imod_settings/
    │       ├── individual_particles/
    │       ├── individual_particles_recs/
    │       ├── individual_particles_ali/
    │       ├── individual_particles_raw/
    │       └── individual_particles_postProc/
    └── logs/
"""

from __future__ import annotations

import json
import os
from typing import Dict, List, Optional

from .star_writer import PairBlock, LoopBlock, write_star

__all__ = [
    "TOMO_SUBDIRS",
    "ensure_project",
    "ensure_tomogram_dirs",
    "project_paths",
    "tomograms_star_path",
    "particles_all_star_path",
    "optimisation_set_star_path",
    "tilt_series_dir",
    "tilt_series_star_path",
    "tomogram_dir",
    "imod_settings_dir",
    "individual_particles_dir",
    "logs_dir",
    "store_path",
    "write_manifest",
]

TOMO_SUBDIRS = [
    "imod_settings",
    "individual_particles",
    "individual_particles_recs",
    "individual_particles_ali",
    "individual_particles_raw",
    "individual_particles_postProc",
]


# --------------------------------------------------------------------------- #
# canonical paths
# --------------------------------------------------------------------------- #
def logs_dir(root: str) -> str:
    return os.path.join(root, "logs")


def tilt_series_dir(root: str) -> str:
    return os.path.join(root, "tilt_series")


def tomograms_star_path(root: str) -> str:
    return os.path.join(root, "tomograms.star")


def particles_all_star_path(root: str) -> str:
    return os.path.join(root, "particles_all.star")


def optimisation_set_star_path(root: str) -> str:
    return os.path.join(root, "optimisation_set.star")


def tilt_series_star_path(root: str, tomo_name: str) -> str:
    return os.path.join(tilt_series_dir(root), f"{tomo_name}.star")


def tomogram_dir(root: str, tomo_name: str) -> str:
    return os.path.join(tilt_series_dir(root), tomo_name)


def imod_settings_dir(root: str, tomo_name: str) -> str:
    return os.path.join(tomogram_dir(root, tomo_name), "imod_settings")


def individual_particles_dir(root: str, tomo_name: str) -> str:
    return os.path.join(tomogram_dir(root, tomo_name), "individual_particles")


def project_paths(root: str) -> Dict[str, str]:
    return {
        "root": root,
        "tomograms_star": tomograms_star_path(root),
        "particles_all_star": particles_all_star_path(root),
        "optimisation_set_star": optimisation_set_star_path(root),
        "manifest_json": os.path.join(root, "project_manifest.json"),
        "manifest_star": os.path.join(root, "project_manifest.star"),
        "tilt_series_dir": tilt_series_dir(root),
        "logs_dir": logs_dir(root),
    }


def store_path(path: str, root: str, relative: bool = True) -> str:
    """Render ``path`` for storage in a STAR file: relative to the project
    root when ``relative`` (and possible), else absolute."""
    ap = os.path.abspath(path)
    if not relative:
        return ap
    try:
        return os.path.relpath(ap, os.path.abspath(root))
    except ValueError:  # e.g. different drive on Windows
        return ap


# --------------------------------------------------------------------------- #
# directory creation
# --------------------------------------------------------------------------- #
def ensure_project(root: str) -> Dict[str, str]:
    """Create the top-level project skeleton (root, tilt_series/, logs/)."""
    os.makedirs(root, exist_ok=True)
    os.makedirs(tilt_series_dir(root), exist_ok=True)
    os.makedirs(logs_dir(root), exist_ok=True)
    return project_paths(root)


def ensure_tomogram_dirs(
    root: str, tomo_name: str, subdirs: Optional[List[str]] = None
) -> str:
    """Create ``tilt_series/<tomo>/`` and its sub-folders. Returns the tomo dir."""
    base = tomogram_dir(root, tomo_name)
    os.makedirs(base, exist_ok=True)
    for sub in (subdirs if subdirs is not None else TOMO_SUBDIRS):
        os.makedirs(os.path.join(base, sub), exist_ok=True)
    return base


# --------------------------------------------------------------------------- #
# manifest
# --------------------------------------------------------------------------- #
def write_manifest(root: str, manifest: Dict[str, object]) -> Dict[str, str]:
    """Write project_manifest.json (machine-readable provenance) and
    project_manifest.star (human-readable summary)."""
    paths = project_paths(root)
    os.makedirs(root, exist_ok=True)

    with open(paths["manifest_json"], "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, default=str)

    # STAR summary: a pairs block for scalars + a loop block listing tomograms.
    scalars = {
        "_tomoJANASProjectRoot": os.path.abspath(root),
        "_tomoJANASVersion": str(manifest.get("tomojanas_version", "")),
        "_tomoJANASCreated": str(manifest.get("created", "")),
        "_tomoJANASUpdated": str(manifest.get("updated", "")),
    }
    blocks = [PairBlock(name="tomoJANAS_project", pairs=scalars)]

    tomos = manifest.get("tomograms", []) or []
    if tomos:
        cols = ["_tomoJANASTomoName", "_tomoJANASImodDir", "_tomoJANASRecTomogram"]
        rows = [
            [t.get("tomo_name"), t.get("imod_dir"), t.get("rec_tomo")]
            for t in tomos
        ]
        blocks.append(LoopBlock(name="tomoJANAS_project_tomograms", columns=cols, rows=rows))

    write_star(paths["manifest_star"], blocks)
    return {"json": paths["manifest_json"], "star": paths["manifest_star"]}
