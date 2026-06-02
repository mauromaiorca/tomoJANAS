[Repository home](../../README.md) · [Docs home](../index.md) · [Installation](../installation.md) · [Quick start](../quick-start.md) · [CLI reference](cli.md) · [Troubleshooting](../troubleshooting.md)

---

# STAR file operations

JANAS reads and writes RELION 3.1 STAR files. This page covers the tools for inspecting and manipulating STAR file contents.

## Inspecting STAR files

### Basic information

```bash
janas_app_starProcess --i particles.star --info
```

Prints the number of particles, column labels, and subset distribution (half-map assignments).

### Euler angle statistics

```bash
janas_app_starProcess --i particles.star --infoEuler
```

Prints statistics on the `_rlnAngleRot`, `_rlnAngleTilt`, and `_rlnAnglePsi` columns.

### Visual Euler distribution

```bash
janas eulerHist --i particles.star
```

Generates a histogram of the Euler angle distribution.

## Exporting

### CSV

```bash
janas_app_starProcess --i particles.star --csv output.csv
```

### VEM

```bash
janas_app_starProcess --i particles.star --vem output.vem
```

### Half-map stacks

```bash
janas_app_starProcess --i particles.star --o out.star --hm half1.mrc half2.mrc
```

Splits the particle stack into two half-map stacks based on `_rlnRandomSubset` assignments.

## Modifying STAR files

### Randomise half-map assignments

```bash
janas_utils randomize_halves --i particles.star --o randomised.star
```

Randomly reassigns particles to two half-sets for independent reconstruction.

### Swap column values

```bash
janas_app_starProcess --i input.star --o output.star --invertTagName tag1 tag2
```

### Back up image name column

```bash
janas_app_starProcess --i input.star --o output.star --backupImageNameTag custom_tag
```

Copies `_rlnImageName` values to a custom column.

### Find similar particles

```bash
janas_app_starProcess --i particles.star --checkForSimilarImages
```

Identifies particles with similar micrograph coordinates, which may indicate duplicates.

## JANAS-specific columns

JANAS adds custom columns to STAR files during processing:

| Column pattern | Description |
|----------------|-------------|
| `_janas_*_norm*` | SCI score normalised by angular coverage |
| `_janas_csparc_rlnImageName` | Original cryoSPARC image name (provenance) |

For backward compatibility, JANAS also reads columns with the `_emprove_` prefix (the project's former name). New output always uses the `_janas_` prefix.

---

[Back to documentation index](../index.md)
