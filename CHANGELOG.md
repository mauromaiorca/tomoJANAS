# Changelog

## [Unreleased]

### Added
- New `tomojanas` Python package (`src/tomojanas/`) for the ET import framework.
- `tomojanas/io/mrc.py` — consolidated MRC I/O with `MRCHeader` dataclass, `1024 + nsymbt` data offsets, memmap, per-slice read/write, spherical/circular crop primitives, and geometry/pixel-size validators. Supports modes 0 (int8), 1 (int16), 2 (float32), 6 (uint16), 12 (float16); complex modes 3/4 raise clear errors.
- `tomojanas/io/star_writer.py` — multi-block STAR writer (loop tables + key/value pair blocks), preserves column order, consistent `?` for missing values.
- `tomojanas/io/star_reader.py` — multi-block STAR reader returning ordered block dicts with DataFrames.
- `tomojanas/io/project_writer.py` — tomoJANAS project directory layout, canonical paths, manifest (JSON + STAR).
- `tomojanas/io/logs.py` — structured import logger and validation report writer (STAR + JSON + Markdown).
- `tomojanas/metadata/relion_labels.py` — canonical RELION and tomoJANAS column orders for all STAR blocks.
- `tomojanas/models/project.py` — project handle with create/append/overwrite/fail-if-existing policy.
- `tests/test_tomojanas_mrc.py` — 13 test cases covering round-trips, nsymbt offsets, mode handling, crops, masks, validators.

### Changed
- `janas/IO_utils.py` reduced to a thin backward-compatibility shim re-exporting from `tomojanas.io.mrc`. Existing code that imports `from janas import IO_utils` continues to work unchanged (dict-based header API preserved).

## [0.1.0] - 2025-06-02

### Added
- `tomojanas import_volume` — import picked coordinates (napari, IMOD, generic CSV), crop sub-volumes from a tomogram, and write a STAR file with metadata.
- Initial fork from JANAS v1.0.0 for Electron Tomography workflows.
