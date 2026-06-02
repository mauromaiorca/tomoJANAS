# Changelog

## [Unreleased]

### Added

**Geometry and metadata modules:**
- `tomojanas/geometry/imod_mapping.py` — faithful Python port of RELION's `ImodImport::import()` algorithm (`imod_import.cpp`). Builds 4×4 `worldToImage` matrices for each non-excluded tilt using the exact RELION matrix composition order. Handles exclusions (EXCLUDELIST/ranges), SHIFT, flip_yz/flip_z/flip_angles, thickness, offsets. Status `relion_imod_algorithm_ported`; requires external oracle for `relion_extract_ready`.
- `tomojanas/geometry/coordinates.py` — RELION centered-Angstrom ↔ reconstruction-voxel conversion with explicit indexing (zero/one-based) and axis-order normalisation (xyz/zyx/yxz).
- `tomojanas/geometry/roi.py` — ROI helpers: sphere/circle scientific ROI, cube/square storage container. `roi_radius_angst_from_args`, `projection_radius_px`, `storage_box_from_radius`, `circle_inside_frame`, `sphere_inside_volume`.
- `tomojanas/metadata/imod_xf.py` — IMOD `.xf` parser; `read_xf`, `apply_xf_to_point`, `invert_xf_for_point`, `xf_matrix`, `xf_inverse_matrix`.
- `tomojanas/metadata/imod_tlt.py` — IMOD `.tlt` parser; `read_tlt`, `read_tlt_angles`.
- `tomojanas/metadata/imod_com.py` — IMOD `.com` file parser; `parse_com_file`, EXCLUDELIST range parsing, SHIFT/XAXISTILT handling.
- `tomojanas/metadata/mdoc.py` — permissive SerialEM `.mdoc` parser; `read_mdoc`.
- `tomojanas/metadata/ctf.py` — CTF auto-detection priority chain, CtfPlotter `.defocus` parser, defocus unit inference and conversion (nm→Å, µm→Å).

**Importers:**
- `tomojanas/importers/imod_importer.py` — `tomojanas-import imod`: reads IMOD directory, builds RELION geometry, writes `tomograms.star` (RELION `data_global` + `data_tomoJANAS_tomogram_sources`), `tilt_series/<tomo>.star` (per-tilt RELION block + `data_tomoJANAS_tilt_mapping` + `data_tomoJANAS_projection_matrices`), `optimisation_set.star`, `imod_settings/*.star`, manifest, and logs.
- `tomojanas/importers/particle_importer.py` — `tomojanas-import particles`: CSV and `--input-single-point`; coordinate conversion to RELION centered Å; per-tilt projection via worldToImage matrices; writes `particles_all.star` (`data_optics` + `data_particles`) and individual `P*.star` files with 5 blocks (optics, particles, source, roi, projections). ROI is sphere/circle; storage is cube/square.
- `tomojanas/importers/ctf_importer.py` — `tomojanas-import ctf`: imports CTF metadata from CtfPlotter/CTFFind files; unit auto-inference; updates per-tilt RELION CTF columns; never modifies image pixels; CTF premultiplication off by default.
- `tomojanas/importers/validators.py` — `tomojanas-import validate`: checks project/tomogram/tilt-series/particle/CTF consistency; writes `validation_log.{star,json,md}`; strict mode returns non-zero on critical failures.

**Labels and status constants:**
- `PARTICLE_OPTICS_COLUMNS`, `RelionGeometryStatus`, `XfDirectionStatus`, `ProjectionStatus` added to `relion_labels.py`. Geometry status vocabulary enforces no fake matrices in RELION columns.
- `PROJECTION_MATRIX_COLUMNS` for scalar 4×4 matrix storage in `data_tomoJANAS_projection_matrices`.

**CLI:**
- `tomojanas-import imod/particles/ctf/validate` commands are now fully functional (no longer stubs).
- Added `--flip-yz`, `--flip-z`, `--flip-angles`, `--thickness-override`, `--import-offset-{x,y,z}`, `--compare-with-relion-import` flags to `tomojanas-import imod`.

**Tests:**
- `tests/test_tomojanas_import.py` — 8 integration tests covering: IMOD parsers (xf/tlt/com/EXCLUDELIST/SHIFT), RELION mapping math (zero-tilt centre, 3-frame identity), full imod import, particle import with single-point and projection tables, validate strict, CTF missing, defocus unit conversion, legacy import_volume. All 21 tests (8+13) pass.

**Infrastructure (previous):**
- `tomojanas/io/mrc.py` — `MRCHeader` dataclass, `1024 + nsymbt` offsets, modes 0/1/2/6/12, crops, masks, validators.
- `tomojanas/io/star_writer.py` / `star_reader.py` — multi-block STAR I/O.
- `tomojanas/io/project_writer.py`, `io/logs.py`, `models/project.py` — project layout and logging.
- `tests/test_tomojanas_mrc.py` — 13 MRC tests.

### Changed
- `janas/IO_utils.py` reduced to a thin backward-compatibility shim re-exporting from `tomojanas.io.mrc`. Existing code that imports `from janas import IO_utils` continues to work unchanged (dict-based header API preserved).

## [0.1.0] - 2025-06-02

### Added
- `tomojanas import_volume` — import picked coordinates (napari, IMOD, generic CSV), crop sub-volumes from a tomogram, and write a STAR file with metadata.
- Initial fork from JANAS v1.0.0 for Electron Tomography workflows.
