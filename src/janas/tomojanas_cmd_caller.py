#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# File: tomojanas_cmd_caller.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology

"""
Command-line interface for tomoJANAS — sub-volume extraction and
analysis pipeline for Electron Tomography.

Commands:
  import_volume   Import picked coordinates (napari / IMOD / generic CSV),
                  crop sub-volumes from the parent tomogram, and write a
                  STAR file that records origin, coordinates and paths.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

from janas.IO_utils import read_mrc_data, read_mrc_header, write_mrc
from janas.starHandler import writeDataframeToStar_deNovo
from janas.version import get_version


# ------------------------------------------------------------------ #
# coordinate readers
# ------------------------------------------------------------------ #

def _read_coords_napari_csv(path):
    """Read a napari Points-layer CSV (columns: index, axis-0, axis-1, axis-2
    which correspond to Z, Y, X)."""
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    if "axis-0" in df.columns:
        coords = df[["axis-0", "axis-1", "axis-2"]].values.astype(float)
        return coords[:, 2], coords[:, 1], coords[:, 0]  # x, y, z
    elif "z" in [c.lower() for c in df.columns]:
        col_map = {c.lower(): c for c in df.columns}
        z = df[col_map["z"]].values.astype(float)
        y = df[col_map["y"]].values.astype(float)
        x = df[col_map["x"]].values.astype(float)
        return x, y, z
    else:
        raise ValueError(
            f"Cannot interpret napari CSV columns: {list(df.columns)}. "
            "Expected 'axis-0,axis-1,axis-2' or 'x,y,z'."
        )


def _read_coords_imod(path):
    """Read an IMOD model2point text file (space-separated X Y Z, one point
    per line).  Coordinates are in pixel units (as output by
    ``model2point -input model.mod -output points.txt``)."""
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data[:, 0], data[:, 1], data[:, 2]  # x, y, z


def _read_coords_csv(path):
    """Generic CSV/TSV with columns x, y, z (header required)."""
    df = pd.read_csv(path, sep=None, engine="python")
    df.columns = [c.strip().lower() for c in df.columns]
    return (
        df["x"].values.astype(float),
        df["y"].values.astype(float),
        df["z"].values.astype(float),
    )


COORD_READERS = {
    "napari": _read_coords_napari_csv,
    "imod": _read_coords_imod,
    "csv": _read_coords_csv,
}


# ------------------------------------------------------------------ #
# sub-volume extraction
# ------------------------------------------------------------------ #

def _crop_subvolume(volume, cx, cy, cz, box):
    """Extract a cubic sub-volume of side *box* centred on (cx, cy, cz).

    If the box extends outside the tomogram the out-of-bounds region is
    filled with the mean value of the volume (common practice in ET to
    avoid harsh edges).
    """
    half = box // 2
    nz, ny, nx = volume.shape

    z0, z1 = int(round(cz)) - half, int(round(cz)) - half + box
    y0, y1 = int(round(cy)) - half, int(round(cy)) - half + box
    x0, x1 = int(round(cx)) - half, int(round(cx)) - half + box

    pad_val = float(np.mean(volume))
    sub = np.full((box, box, box), pad_val, dtype=np.float32)

    src_z0 = max(z0, 0);       src_z1 = min(z1, nz)
    src_y0 = max(y0, 0);       src_y1 = min(y1, ny)
    src_x0 = max(x0, 0);       src_x1 = min(x1, nx)

    dst_z0 = src_z0 - z0;      dst_z1 = dst_z0 + (src_z1 - src_z0)
    dst_y0 = src_y0 - y0;      dst_y1 = dst_y0 + (src_y1 - src_y0)
    dst_x0 = src_x0 - x0;      dst_x1 = dst_x0 + (src_x1 - src_x0)

    sub[dst_z0:dst_z1, dst_y0:dst_y1, dst_x0:dst_x1] = \
        volume[src_z0:src_z1, src_y0:src_y1, src_x0:src_x1]
    return sub


# ------------------------------------------------------------------ #
# import_volume command
# ------------------------------------------------------------------ #

tomojanas_parser = argparse.ArgumentParser(
    prog="tomojanas",
    usage="%(prog)s [command] [arguments]",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

tomojanas_parser.add_argument(
    "-V", "--version",
    action="version",
    version=get_version(),
)

command = tomojanas_parser.add_subparsers(dest="command")

import_volume_parser = command.add_parser(
    "import_volume",
    description=(
        "Import picked coordinates from napari, IMOD or a generic CSV, "
        "crop sub-volumes from the parent tomogram, and write a STAR "
        "file that records origin, coordinates and output paths."
    ),
    help="crop sub-volumes at picked coordinates and produce a STAR file",
)
import_volume_parser.add_argument(
    "--tomo", required=True, type=str,
    help="path to the input tomogram (MRC)",
)
import_volume_parser.add_argument(
    "--coords", required=True, type=str,
    help="path to the coordinates file (napari CSV, IMOD points, or generic CSV)",
)
import_volume_parser.add_argument(
    "--format", required=False, type=str, default="napari",
    choices=list(COORD_READERS.keys()),
    help="coordinate file format (default: napari)",
)
import_volume_parser.add_argument(
    "--box", required=True, type=int,
    help="box size in voxels (cubic) for the extracted sub-volumes",
)
import_volume_parser.add_argument(
    "--apix", required=False, type=float, default=None,
    help="pixel size in Angstroms (read from MRC header if omitted)",
)
import_volume_parser.add_argument(
    "--o", required=False, type=str, default=None,
    help="output STAR file (default: <tomo_basename>_import.star)",
)
import_volume_parser.add_argument(
    "--outDir", required=False, type=str, default=None,
    help="output directory for sub-volumes (default: ./subvolumes/)",
)


def import_volume(args):
    tomo_path = os.path.abspath(args.tomo)
    if not os.path.isfile(tomo_path):
        print(f'ERROR: tomogram "{args.tomo}" not found')
        sys.exit(1)

    coords_path = os.path.abspath(args.coords)
    if not os.path.isfile(coords_path):
        print(f'ERROR: coordinates file "{args.coords}" not found')
        sys.exit(1)

    reader = COORD_READERS[args.format]
    x_coords, y_coords, z_coords = reader(coords_path)
    n_particles = len(x_coords)
    print(f"Read {n_particles} coordinates from {args.coords} (format={args.format})")

    hdr = read_mrc_header(tomo_path)
    apix = args.apix if args.apix else hdr["pixel_x"]
    print(f"Tomogram: {hdr['nx']}x{hdr['ny']}x{hdr['nz']}, pixel size={apix:.4f} A")

    volume, _ = read_mrc_data(tomo_path)
    box = args.box

    tomo_basename = os.path.splitext(os.path.basename(tomo_path))[0]
    out_dir = args.outDir if args.outDir else os.path.join(".", "subvolumes")
    os.makedirs(out_dir, exist_ok=True)

    out_star = args.o if args.o else f"{tomo_basename}_import.star"

    records = []
    for i in range(n_particles):
        cx, cy, cz = float(x_coords[i]), float(y_coords[i]), float(z_coords[i])
        sub = _crop_subvolume(volume, cx, cy, cz, box)

        sub_name = f"{tomo_basename}_{i:06d}.mrc"
        sub_path = os.path.join(out_dir, sub_name)
        write_mrc(sub_path, sub, apix)

        records.append({
            "_tomoSourceTomogram": tomo_path,
            "_tomoSubVolumePath": os.path.abspath(sub_path),
            "_tomoCoordinateX": f"{cx:.2f}",
            "_tomoCoordinateY": f"{cy:.2f}",
            "_tomoCoordinateZ": f"{cz:.2f}",
            "_tomoBoxSize": str(box),
            "_tomoPixelSize": f"{apix:.4f}",
        })

        if (i + 1) % 100 == 0 or (i + 1) == n_particles:
            print(f"  extracted {i + 1}/{n_particles} sub-volumes")

    df = pd.DataFrame(records)
    writeDataframeToStar_deNovo(df, out_star)
    print(f"STAR file written to {out_star}")
    print(f"Sub-volumes saved in {os.path.abspath(out_dir)}/")


# ------------------------------------------------------------------ #
# main dispatcher
# ------------------------------------------------------------------ #

def main(command_line=None):
    args = tomojanas_parser.parse_args(command_line)
    if args.command == "import_volume":
        import_volume(args)
    else:
        tomojanas_parser.print_help()


if __name__ == "__main__":
    main()
