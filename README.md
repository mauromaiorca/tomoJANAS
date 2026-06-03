<h1 align="center">tomoJANAS</h1>
<p align="center"><strong>Import, metadata and analysis framework for Electron Tomography</strong></p>

---

tomoJANAS is a command-line framework for importing IMOD tomogram projects, managing picked particle coordinates, tracking per-particle metadata (spherical/circular ROIs, tilt projections, CTF), and writing RELION-compatible STAR files for downstream sub-tomogram averaging.

tomoJANAS is **not** a RELION clone. It has its own project structure and uses RELION-compatible `_rln*` tags in dedicated STAR blocks, alongside `_tomoJANAS*` tags in separate blocks for ET-specific provenance.

📖 **[Examples of usage](docs/examples_of_usage.md)** — a growing collection of practical recipes (IMOD import, coordinate import, crops, CTF, validation, status/reconcile, command replay).

## Key concepts

- **Biological ROI**: spherical in 3D, circular in 2D projections. Cubes/squares are only MRC storage containers.
- **Coordinates**: the primary representation is RELION centered Angstroms (`_rlnCenteredCoordinate{X,Y,Z}Angst`). Original picked coordinates are preserved in `_tomoJANAS*` tags with explicit indexing and axis-order metadata.
- **CTF**: imported as metadata only; pixel data is never modified during import (no premultiplication by default).
- **Multi-block STAR**: RELION-compatible blocks are kept free of tomoJANAS tags; custom metadata lives in separate `data_tomoJANAS_*` blocks within the same file.

## Commands

### Import commands (`tomojanas-import`)

| Command | Purpose |
|---------|---------|
| `tomojanas-import imod` | Import an IMOD directory (stacks, transforms, tilt angles) into a tomoJANAS project |
| `tomojanas-import particles` | Import picked coordinates (CSV, STAR, napari, IMOD, single point) and compute per-tilt projections |
| `tomojanas-import ctf` | Import CTF metadata from CtfPlotter, CTFFind or IMOD without modifying images |
| `tomojanas-import validate` | Validate project consistency (binning, tilt counts, coordinates, CTF coverage) |
| `tomojanas-import status` | Scan the project tree, report registered vs. on-disk particles, reconcile (`--sync`), and create missing sub-volumes (`--create-volume`) |

**Creating sub-volumes:** add `--write-rec-crops` to `tomojanas-import particles` to extract the 3D sub-volume at import time, or run `tomojanas-import status --create-volume` to backfill sub-volumes for particles imported earlier.

**Coordinate axis order:** the default is `xyz` — the three input numbers map directly to the rec file's X, Y, Z axes (the standard way). Pass `--axis-order xzy` for a flipped IMOD tomogram (thickness on the middle axis) or `--axis-order zyx` for napari single points; `--coord-flip-{x,y,z}` inverts an axis if needed.

Every invocation is appended (with its exit status) to `logs/commands.jsonl` and a replayable `logs/commands.sh`, so a project's import history can be reproduced in a new project.

### Legacy command (`tomojanas`)

| Command | Purpose |
|---------|---------|
| `tomojanas import_volume` | Quick sub-volume extraction from picked coordinates (backward-compatible) |

## Project layout

When you import an IMOD project, tomoJANAS creates a structured directory:

```
my_project/
├── tomograms.star              # RELION-compatible tomogram table
├── particles_all.star          # RELION-compatible global particle table
├── optimisation_set.star       # links tomograms + particles files
├── project_manifest.star       # project provenance (tomoJANAS)
├── project_manifest.json       # machine-readable provenance
├── tilt_series/
│   ├── lam8_ts_006.star        # per-tilt RELION table + tomoJANAS tilt mapping
│   └── lam8_ts_006/
│       ├── imod_settings/      # imported IMOD metadata (.xf, .tlt, CTF, ...)
│       ├── individual_particles/   # one STAR file per particle (P000001.star)
│       ├── individual_particles_recs/  # optional 3D crops
│       ├── individual_particles_ali/   # optional 2D aligned crops
│       └── individual_particles_raw/   # optional 2D raw crops
└── logs/
    ├── tomojanas_import.log
    └── validation_log.{star,json,md}
```

## Quick start

```bash
# 1. Import an IMOD tomogram project
tomojanas-import imod \
    --project my_project --create-if-missing \
    --imod-dir /data/lam8_ts_006 \
    --tomo-name lam8_ts_006 \
    --ali-stack lam8_ts_006_ali.mrc \
    --rec-tomo lam8_ts_006_rec.mrc \
    --xf lam8_ts_006.xf --tlt lam8_ts_006.tlt \
    --validate

# 2. Import a single picked particle
tomojanas-import particles \
    --project my_project \
    --tomo-name lam8_ts_006 \
    --input-single-point 512.3,401.2,188.0 \
    --coordinate-system rec-voxel --indexing zero-based \
    --roi-radius-angst 150 --validate

# 3. Import particles from a CSV
tomojanas-import particles \
    --project my_project \
    --tomo-name lam8_ts_006 \
    --input picks.csv --format csv \
    --roi-radius-angst 150 --validate

# 4. Import CTF metadata
tomojanas-import ctf \
    --project my_project \
    --tomo-name lam8_ts_006 \
    --ctf-source auto --do-not-premultiply --validate

# 5. Validate the whole project
tomojanas-import validate \
    --project my_project --strict
```

### Legacy quick extraction (no project required)

```bash
tomojanas import_volume \
    --tomo tomogram.mrc --coords picks.csv --format napari --box 64
```

## STAR file structure

### tomograms.star

Two blocks: a RELION-compatible `data_global` with `_rln*` tags, and a separate `data_tomoJANAS_tomogram_sources` block with import provenance.

### Individual particle files (P000001.star)

Four blocks per particle:

| Block | Tags | Purpose |
|-------|------|---------|
| `data_particles` | `_rln*` | RELION-compatible coordinates, orientations, optics |
| `data_tomoJANAS_particle_source` | `_tomoJANAS*` | Original picked coordinate, system, indexing, software |
| `data_tomoJANAS_particle_roi` | `_tomoJANAS*` | Spherical 3D ROI, circular 2D projection ROI, storage box |
| `data_tomoJANAS_particle_projections` | `_tomoJANAS*` | Per-tilt projected centre (aligned + raw), visibility, CTF |

## Installation

### Requirements

- Python 3.8+
- A C++ compiler (gcc, clang, or MSVC)
- CMake 3.10+

### First-time install

```bash
# 1. Clone the repository
git clone https://github.com/mauromaiorca/tomoJANAS.git
cd tomoJANAS

# 2. Create and activate a virtual environment (recommended)
python3 -m venv .tomojanas_env
source .tomojanas_env/bin/activate

# 3. Install in editable mode
pip install -e .
```

### Updating to the latest version

```bash
cd tomoJANAS
git pull origin main
pip install -e .
```

### Verify installation

```bash
tomojanas --version
tomojanas-import --help
tomojanas-import imod --help
```

## Running tests

```bash
python tests/test_tomojanas_mrc.py
```

## Contact

For questions or issues: mauro.maiorca@cssb-hamburg.de
