#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/importers/cli.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
CLI entry point for ``tomojanas-import``.

Sub-commands:
    tomojanas-import imod       Import an IMOD directory into a tomoJANAS project
    tomojanas-import particles  Import picked coordinates
    tomojanas-import ctf        Import CTF metadata (without modifying pixels)
    tomojanas-import validate   Validate project consistency
"""

from __future__ import annotations

import argparse
import sys

from tomojanas import get_version


# ------------------------------------------------------------------ #
# top-level parser
# ------------------------------------------------------------------ #
def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tomojanas-import",
        description="tomoJANAS import tools for Electron Tomography.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-V", "--version", action="version", version=get_version(),
    )
    sub = parser.add_subparsers(dest="subcommand")

    # ---- imod --------------------------------------------------------- #
    imod_p = sub.add_parser(
        "imod",
        help="import an IMOD directory into a tomoJANAS project",
        description="Import an IMOD tomogram reconstruction directory.",
    )
    _add_project_args(imod_p)
    imod_p.add_argument("--imod-dir", required=True, help="path to the IMOD directory")
    imod_p.add_argument("--tomo-name", required=True, help="tomoJANAS name for this tomogram")
    imod_p.add_argument("--basename", default=None, help="IMOD file basename (default: tomo-name)")
    imod_p.add_argument("--raw-stack", default=None, help="raw tilt-series stack (.mrc)")
    imod_p.add_argument("--mdoc", default=None, help="SerialEM .mdoc file")
    imod_p.add_argument("--ali-stack", default=None, help="aligned stack (.mrc)")
    imod_p.add_argument("--rec-tomo", default=None, help="reconstructed tomogram (.mrc)")
    imod_p.add_argument("--xf", default=None, help="IMOD .xf transform file")
    imod_p.add_argument("--tlt", default=None, help="IMOD .tlt tilt-angle file")
    imod_p.add_argument("--rawtlt", default=None, help="IMOD .rawtlt file")
    imod_p.add_argument("--xtilt", default=None, help="IMOD .xtilt file")
    imod_p.add_argument("--newst-com", default=None, help="newst.com command file")
    imod_p.add_argument("--tilt-com", default=None, help="tilt.com command file")
    # pixel-size / binning
    imod_p.add_argument("--raw-pixel-size", type=float, default=None)
    imod_p.add_argument("--ali-pixel-size", type=float, default=None)
    imod_p.add_argument("--rec-pixel-size", type=float, default=None)
    imod_p.add_argument("--tomo-binning", type=float, default=None)
    imod_p.add_argument("--infer-pixel-size-from-mrc", action="store_true", default=True)
    imod_p.add_argument("--infer-pixel-size-from-mdoc", action="store_true", default=False)
    imod_p.add_argument("--require-ali-rec-same-bin", action="store_true", default=False)
    imod_p.add_argument("--allow-rec-binning", type=float, default=None)
    # CTF
    imod_p.add_argument("--ctf-source", default="auto",
                        choices=["auto", "ctfplotter", "ctffind", "imod", "none"])
    imod_p.add_argument("--ctfplotter-file", default=None)
    imod_p.add_argument("--ctfplotter-log", default=None)
    imod_p.add_argument("--ctfplotter-info", default=None)
    imod_p.add_argument("--ctffind-star", default=None)
    imod_p.add_argument("--do-not-premultiply", action="store_true", default=True)
    imod_p.add_argument("--premultiply-ctf", dest="do_not_premultiply", action="store_false")
    # micrograph reference
    imod_p.add_argument("--micrograph-reference", default="ali", choices=["ali", "raw"])
    # IMOD geometry flags (for RELION mapping port)
    imod_p.add_argument("--flip-yz", action="store_true", default=False)
    imod_p.add_argument("--flip-z", action="store_true", default=False)
    imod_p.add_argument("--flip-angles", action="store_true", default=False)
    imod_p.add_argument("--thickness-override", type=float, default=None)
    imod_p.add_argument("--import-offset-x", type=float, default=0.0)
    imod_p.add_argument("--import-offset-y", type=float, default=0.0)
    imod_p.add_argument("--import-offset-z", type=float, default=0.0)
    # Optional RELION oracle validation
    imod_p.add_argument("--compare-with-relion-import", action="store_true", default=False)
    imod_p.add_argument("--relion-binary", default="relion_tomo_import_tomograms",
                        help="RELION binary for oracle validation")
    imod_p.add_argument("--relion-import-workdir", default=None)
    imod_p.add_argument("--relion-projection-tolerance-pixel", type=float, default=0.5)
    # validation
    _add_validation_args(imod_p)

    # ---- particles ---------------------------------------------------- #
    part_p = sub.add_parser(
        "particles",
        help="import picked particle coordinates",
        description="Import picked coordinates and write particle STAR files.",
    )
    _add_project_args(part_p)
    part_p.add_argument("--tomo-name", required=True)
    part_p.add_argument("--input", default=None, help="coordinates file")
    part_p.add_argument("--input-single-point", default=None,
                        help="single coordinate as X,Y,Z")
    part_p.add_argument("--format", default="csv",
                        choices=["csv", "star", "napari", "imod-model", "imod-point", "mod", "json"])
    part_p.add_argument("--coordinate-system", default="auto",
                        choices=["auto", "rec-voxel", "ali-pixel", "raw-pixel",
                                 "relion-centered-angst", "imod-model"])
    part_p.add_argument("--indexing", default="auto", choices=["auto", "zero-based", "one-based"])
    part_p.add_argument("--axis-order", default="auto",
                        choices=["auto", "xyz", "xzy", "yxz", "yzx", "zxy", "zyx"],
                        help="order of the input coordinate columns. IMOD 3dmod "
                             "readouts of a flipped tomogram are usually 'xzy' "
                             "(Y and Z swapped); napari points are 'zyx'.")
    part_p.add_argument("--roi-shape", default="sphere")
    part_p.add_argument("--projection-roi-shape", default="circle")
    part_p.add_argument("--roi-radius-angst", type=float, default=None)
    part_p.add_argument("--roi-diameter-angst", type=float, default=None)
    part_p.add_argument("--roi-radius-voxel", type=float, default=None)
    part_p.add_argument("--roi-padding-angst", type=float, default=0.0)
    part_p.add_argument("--particle-prefix", default="P")
    part_p.add_argument("--particle-name", default=None)
    part_p.add_argument("--start-particle-id", type=int, default=1)
    part_p.add_argument("--write-individual-particles", action="store_true", default=True)
    part_p.add_argument("--update-particles-all", action="store_true", default=True)
    part_p.add_argument("--write-rec-crops", action="store_true", default=False)
    part_p.add_argument("--write-ali-crops", action="store_true", default=False)
    part_p.add_argument("--write-raw-crops", action="store_true", default=False)
    part_p.add_argument("--crop-storage-box-size", type=int, default=None)
    part_p.add_argument("--crop-padding-voxel", type=float, default=0.0)
    part_p.add_argument("--crop-padding-angst", type=float, default=0.0)
    part_p.add_argument("--crop-outside-policy", default="partial",
                        choices=["error", "pad", "partial", "skip"])
    part_p.add_argument("--crop-pad-value", type=float, default=0.0)
    part_p.add_argument("--apply-spherical-mask", action="store_true", default=False)
    part_p.add_argument("--save-spherical-mask", action="store_true", default=False)
    part_p.add_argument("--save-circular-mask", action="store_true", default=False)
    part_p.add_argument("--coordinate-auto-threshold", type=float, default=0.90)
    part_p.add_argument("--coordinate-auto-margin", type=float, default=0.10)
    part_p.add_argument("--fail-if-coordinate-ambiguous", action="store_true", default=False)
    _add_validation_args(part_p)

    # ---- ctf ---------------------------------------------------------- #
    ctf_p = sub.add_parser(
        "ctf",
        help="import CTF metadata (no pixel modification)",
        description="Import or update CTF metadata without modifying image pixels.",
    )
    _add_project_args(ctf_p)
    ctf_p.add_argument("--tomo-name", required=True)
    ctf_p.add_argument("--ctf-source", default="auto",
                       choices=["auto", "ctfplotter", "ctffind", "imod", "none"])
    ctf_p.add_argument("--ctfplotter-file", default=None)
    ctf_p.add_argument("--ctfplotter-log", default=None)
    ctf_p.add_argument("--ctfplotter-info", default=None)
    ctf_p.add_argument("--ctffind-star", default=None)
    ctf_p.add_argument("--defocus-file", default=None)
    ctf_p.add_argument("--defocus-unit", default="auto",
                       choices=["auto", "angstrom", "nm", "micrometer"])
    ctf_p.add_argument("--phase-shift", type=float, default=None)
    ctf_p.add_argument("--amplitude-contrast", type=float, default=None)
    ctf_p.add_argument("--voltage", type=float, default=None)
    ctf_p.add_argument("--spherical-aberration", type=float, default=None)
    ctf_p.add_argument("--do-not-premultiply", action="store_true", default=True)
    _add_validation_args(ctf_p)

    # ---- validate ----------------------------------------------------- #
    val_p = sub.add_parser(
        "validate",
        help="validate project consistency",
        description="Validate a tomoJANAS project, one tomogram, or one particle.",
    )
    _add_project_args(val_p)
    val_p.add_argument("--tomo-name", default=None,
                       help="validate only this tomogram (omit for whole project)")
    val_p.add_argument("--check-binning", action="store_true", default=True)
    val_p.add_argument("--check-tilt-count", action="store_true", default=True)
    val_p.add_argument("--check-coordinate-roundtrip", action="store_true", default=True)
    val_p.add_argument("--check-projections", action="store_true", default=True)
    val_p.add_argument("--check-raw-links", action="store_true", default=True)
    val_p.add_argument("--check-ctf", action="store_true", default=True)
    val_p.add_argument("--check-particles", action="store_true", default=True)
    val_p.add_argument("--tolerance-pixel", type=float, default=0.5)
    val_p.add_argument("--tolerance-angst", type=float, default=1.0)
    val_p.add_argument("--write-log", action="store_true", default=True)
    _add_validation_args(val_p)

    # ---- status ------------------------------------------------------- #
    status_p = sub.add_parser(
        "status",
        help="scan the project tree and report present/missing particles",
        description=(
            "Scan a project (or one tomogram) and report which particles are "
            "registered, which individual P*.star files and crops exist on "
            "disk, and any inconsistencies. Use --sync to rebuild "
            "particles_all.star from the P*.star files actually present "
            "(e.g. after deleting particle files manually)."
        ),
    )
    status_p.add_argument("--project", required=True, help="path to the tomoJANAS project root")
    status_p.add_argument("--tomo-name", default=None,
                          help="scan only this tomogram (omit for all tomograms)")
    status_p.add_argument("--sync", action="store_true", default=False,
                          help="rebuild particles_all.star from the P*.star files on disk")
    status_p.add_argument("--strict", action="store_true", default=False,
                          help="exit non-zero if any inconsistency is found")
    # create the 3D rec sub-volumes for already-imported particles
    status_p.add_argument("--create-volume", "--create-rec-crops",
                          dest="create_volume", action="store_true", default=False,
                          help="create the 3D rec sub-volume (_rec.mrc) for each "
                               "particle that does not have one yet")
    status_p.add_argument("--overwrite-crops", action="store_true", default=False,
                          help="with --create-volume, regenerate crops even if they exist")
    status_p.add_argument("--crop-storage-box-size", type=int, default=None)
    status_p.add_argument("--crop-padding-voxel", type=float, default=0.0)
    status_p.add_argument("--crop-padding-angst", type=float, default=0.0)
    status_p.add_argument("--crop-outside-policy", default="partial",
                          choices=["error", "pad", "partial", "skip"])
    status_p.add_argument("--crop-pad-value", type=float, default=0.0)
    status_p.add_argument("--apply-spherical-mask", action="store_true", default=False)

    return parser


# ------------------------------------------------------------------ #
# shared argument groups
# ------------------------------------------------------------------ #
def _add_project_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--project", required=True, help="path to the tomoJANAS project root")
    p.add_argument("--create-if-missing", action="store_true", default=False)
    p.add_argument("--append", action="store_true", default=False)
    p.add_argument("--update-existing", action="store_true", default=False)
    p.add_argument("--overwrite", action="store_true", default=False)
    p.add_argument("--fail-if-existing", action="store_true", default=True)
    p.add_argument("--relative-paths", action="store_true", default=True)
    p.add_argument("--absolute-paths", dest="relative_paths", action="store_false")


def _add_validation_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--validate", action="store_true", default=False)
    p.add_argument("--strict", action="store_true", default=False)
    p.add_argument("--check-mrc-headers", action="store_true", default=False)
    p.add_argument("--check-ali-rec-pixel-size", action="store_true", default=False)
    p.add_argument("--check-ali-tilt-count", action="store_true", default=False)
    p.add_argument("--check-rec-dimensions", action="store_true", default=False)


# ------------------------------------------------------------------ #
# sub-command dispatchers (stubs — implementations in later phases)
# ------------------------------------------------------------------ #
def _cmd_imod(args: argparse.Namespace) -> int:
    from tomojanas.importers.imod_importer import import_imod_project
    return import_imod_project(args)


def _cmd_particles(args: argparse.Namespace) -> int:
    from tomojanas.importers.particle_importer import import_particles
    return import_particles(args)


def _cmd_ctf(args: argparse.Namespace) -> int:
    from tomojanas.importers.ctf_importer import import_ctf
    return import_ctf(args)


def _cmd_validate(args: argparse.Namespace) -> int:
    from tomojanas.importers.validators import validate_project_cli
    return validate_project_cli(args)


def _cmd_status(args: argparse.Namespace) -> int:
    from tomojanas.importers.status import status_cli
    return status_cli(args)


_DISPATCH = {
    "imod": _cmd_imod,
    "particles": _cmd_particles,
    "ctf": _cmd_ctf,
    "validate": _cmd_validate,
    "status": _cmd_status,
}


# ------------------------------------------------------------------ #
# entry point
# ------------------------------------------------------------------ #
def main(argv=None) -> None:
    import os
    raw_argv = list(sys.argv[1:]) if argv is None else list(argv)
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.subcommand is None:
        parser.print_help()
        sys.exit(0)
    handler = _DISPATCH.get(args.subcommand)
    if handler is None:
        parser.print_help()
        sys.exit(1)
    rc = handler(args)

    # Append the executed command (with exit status) to the project command log.
    project_root = getattr(args, "project", None)
    if project_root and os.path.isdir(project_root):
        try:
            from tomojanas.io.logs import record_command
            from tomojanas import get_version
            record_command(project_root, raw_argv, rc if rc is not None else 0,
                           version=get_version())
        except Exception:
            pass

    sys.exit(rc if rc is not None else 0)


if __name__ == "__main__":
    main()
