<h1 align="center">tomoJANAS</h1>
<p align="center"><strong>Sub-volume extraction and analysis pipeline for Electron Tomography</strong></p>

---

tomoJANAS is a command-line toolkit for sub-volume extraction and analysis in Electron Tomography (ET) workflows.

It reads picked coordinates from common software (napari, IMOD, or generic CSV), crops sub-volumes from a parent tomogram, and tracks all metadata in RELION-compatible STAR files.

## Available commands

| Command | Purpose |
|---------|---------|
| `import_volume` | Import picked coordinates, crop sub-volumes from a tomogram, and write a STAR file with source paths, coordinates, box size and pixel size. |

## Installation

Requires Python 3.8+, a C++ compiler, and CMake 3.10+.

```bash
pip install -e .
```

We recommend installing in an isolated environment:

```bash
python3 -m venv ~/.tomojanas_env
source ~/.tomojanas_env/bin/activate
pip install -e .
```

Verify:

```bash
tomojanas --version
```

## Quick start

```bash
# Import sub-volumes from napari picks
tomojanas import_volume \
    --tomo tomogram.mrc \
    --coords picks.csv \
    --format napari \
    --box 64

# Import from IMOD model2point output
tomojanas import_volume \
    --tomo tomogram.mrc \
    --coords points.txt \
    --format imod \
    --box 64

# Import from a generic CSV (columns: x, y, z)
tomojanas import_volume \
    --tomo tomogram.mrc \
    --coords coords.csv \
    --format csv \
    --box 64 \
    --apix 3.42 \
    --outDir ./my_subvolumes \
    --o my_particles.star
```

## Output STAR file

The output STAR file contains one row per extracted sub-volume with the following columns:

| Column | Description |
|--------|-------------|
| `_tomoSourceTomogram` | Absolute path to the source tomogram |
| `_tomoSubVolumePath` | Absolute path to the extracted sub-volume MRC |
| `_tomoCoordinateX` | X coordinate (pixels) in the source tomogram |
| `_tomoCoordinateY` | Y coordinate (pixels) in the source tomogram |
| `_tomoCoordinateZ` | Z coordinate (pixels) in the source tomogram |
| `_tomoBoxSize` | Box size (voxels) used for extraction |
| `_tomoPixelSize` | Pixel size (Angstroms) |

## Contact

For questions or issues: mauro.maiorca@cssb-hamburg.de
