[Repository home](../../README.md) · [Docs home](../index.md) · [Installation](../installation.md) · [Quick start](../quick-start.md) · [STAR operations](star-operations.md) · [Troubleshooting](../troubleshooting.md)

---

# CLI commands

JANAS installs seven command-line programs. This page lists each command and its primary subcommands.

## janas

Main entry point for particle scoring and analysis.

```bash
janas <subcommand> [options]
```

| Subcommand | Description |
|------------|-------------|
| `scoreParticles` | Score particles against a reference volume using SCI |
| `eulerHist` | Plot Euler angle distribution histogram |

### Score particles

```bash
janas scoreParticles \
    --i particles.star \
    --map reference.mrc \
    --mask mask.mrc \
    --sigma 1 \
    --mpi 40
```

### Euler histogram

```bash
janas eulerHist --i particles.star
```

## janas_utils

Utilities for masks, image processing, local resolution, and cryoSPARC import.

```bash
janas_utils <subcommand> [options]
```

| Subcommand | Description |
|------------|-------------|
| `randomize_halves` | Randomise half-map assignments in a STAR file |
| `maskedCrop` | Compute automatic crop of an image based on a mask |
| `locresBulk` | Compute local resolution for multiple half-map pairs |
| `locresStats` | Extract statistics from a local resolution map |
| `equalize_images` | Equalise amplitude spectra across multiple maps |
| `csparc2star-stack` | Convert cryoSPARC `.cs` to STAR + assembled stack |
| `backmap_stars` | Restore original `_rlnImageName` in a downstream STAR using the stack-generation STAR (inverse of `create_stack`) |
| `csparc2star` | Convert cryoSPARC `.cs` to STAR (no stack assembly) |
| `update_from_csparc` | Update STAR file metadata from a cryoSPARC `.cs` file |
| `csparc_localnurefinement` | Run cryoSPARC local NU-refinement |
| `clip average` | Average two volumes |

### Examples

```bash
# Randomise half-map assignments
janas_utils randomize_halves --i particles.star --o randomised.star

# Compute automatic crop
janas_utils maskedCrop --map volume.mrc --mask mask.mrc

# Convert cryoSPARC particles
janas_utils csparc2star-stack /path/to/particles.cs output_prefix
```

## janas_session_manager

Create and manage selection and classification sessions.

```bash
janas_session_manager <session_type> [options]
```

### new_select_session

```bash
janas_session_manager new_select_session \
    --name my_selection \
    --particles particles.star \
    --map halfA.mrc \
    --map2 halfB.mrc \
    --mask mask.mrc \
    --sigma 1 \
    --mpi 40 \
    --bootstrap \
    --numRecs 12 \
    --maxSelections 14 \
    --noExternalPrograms --gpu 0 1
```

### classification_session

```bash
janas_session_manager classification_session \
    --name reclassify \
    --particles particles.star \
    --maps class1.mrc class2.mrc class3.mrc \
    --mask mask.mrc \
    --mpi 40 \
    --ctf-mode phaseflip \
    --noExternalPrograms --gpu 0 1
```

`--ctf-mode` chooses how the CTF is applied during particle scoring. Choices: `phaseflip` (default — applies `sign(CTF)`), `modulate` (multiplies by the full CTF), `wiener` (`CTF / (CTF² + 0.1)`). The same option is available on `new_select_session`. The selected mode is forwarded to the `janas scoreParticles` calls in the generated run script.

`--noRecs` (classification_session only): skip per-class reconstruction. The run script performs scoring and assignment, and writes per-class STAR files, but does not call any reconstruction backend. Use this when you prefer to reconstruct each class yourself (RELION, cryoSPARC, `janas_reconstructor`, etc.).

> **About `--noExternalPrograms` and `--gpu`:** these flags affect **only the reconstruction and local resolution steps** — particle scoring always runs on CPU (controlled by `--mpi`). With GPU(s), `--noExternalPrograms --gpu 0 1` (recommended) uses two GPUs, one per half-map; `--gpu 0` uses one. On CPU-only machines, omitting `--noExternalPrograms` lets JANAS call RELION's MPI-based reconstruction and `relion_postprocess`, which is typically faster than JANAS's internal CPU reconstruction — make sure RELION is installed and accessible on your `PATH`.

### random_selection_session

For control experiments with random particle selection:

```bash
janas_session_manager random_selection_session \
    --name random_control \
    --particles particles.star \
    --mask mask.mrc \
    --dir output_dir
```

## janas_reconstructor

3D reconstruction from scored particles. GPU is strongly recommended; CPU-only operation is supported but **very slow** on medium-to-large datasets. On CPU-only machines, prefer RELION (omit `--noExternalPrograms` in session-manager workflows).

Single-volume reconstruction:

```bash
janas_reconstructor \
    --i particles.star \
    --o output.mrc \
    --gpu 0
```

Independent half-map reconstruction (`output_recH1.mrc`, `output_recH2.mrc`):

```bash
janas_reconstructor \
    --i particles.star \
    --o output \
    --gpu 0 \
    --subset 1 2
```

| Option | Description |
|--------|-------------|
| `--gpu` | GPU device indices (e.g., `0` or `0 1`); omit for CPU-only |
| `--gpu-batch` | Batch size for GPU processing (default 20) |
| `--cpu` | Number of CPU workers (for CPU-only mode) |
| `--subset` | Half-map indices to reconstruct (e.g., `1 2` for both half-maps) |
| `--subrec-only` | Particle counts for subset reconstructions |

## janas_optimizer

Overview analysis and optimisation plots.

```bash
janas_optimizer plotOverview --overview overview.txt --plot
```

## janas_app_starProcess

C++ tool for STAR file inspection and manipulation.

```bash
janas_app_starProcess --i input.star [--o output.star] [options]
```

| Option | Description |
|--------|-------------|
| `--info` | Display particle count, labels, and subset distribution |
| `--infoEuler` | Display Euler angle statistics |
| `--csv output.csv` | Export to CSV format |
| `--vem output.vem` | Export to VEM format |
| `--hm h1.mrc h2.mrc [tag]` | Export particles to two half-map stacks |
| `--micrographs [depth]` | Extract unique micrograph list |
| `--checkForSimilarImages` | Find particles with similar coordinates |
| `--invertTagName tag1 tag2` | Swap two column values |
| `--backupImageNameTag [tag]` | Back up `_rlnImageName` to a custom tag |

Run `janas_app_starProcess --h` for the full list.

### Examples

```bash
# Inspect a STAR file
janas_app_starProcess --i particles.star --info

# Export to CSV
janas_app_starProcess --i particles.star --csv particles.csv

# Split into half-map stacks
janas_app_starProcess --i particles.star --o out.star --hm half1.mrc half2.mrc
```

## janas_app_meanMinMax

C++ tool for computing statistics from a local resolution map within a masked region.

```bash
janas_app_meanMinMax locresMap.mrc mask.mrc
```

Returns the mean, minimum, and maximum resolution values within the mask.

---

[Back to documentation index](../index.md)
