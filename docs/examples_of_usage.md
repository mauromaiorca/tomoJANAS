# Examples of usage

A growing collection of practical `tomojanas-import` recipes. Replace
`first` / `lam8_ts_006` / paths with your own.

> **Conventions**
> - The biological ROI is a **sphere** in 3D (a **circle** in 2D projections); the box is only a storage container.
> - Particle coordinates are stored as RELION centered Ångströms (`_rlnCenteredCoordinate*Angst`); original picks are kept in `_tomoJANAS*` tags.
> - **External** IMOD source files are stored with **absolute** paths; **internal** project files use relative paths.

---

## 1. Import an IMOD tomogram project

```bash
tomojanas-import imod \
  --project first --create-if-missing \
  --imod-dir /path/to/lam8_ts_006 \
  --tomo-name lam8_ts_006 \
  --rec-tomo lam8_ts_006_rec.mrc \
  --xf lam8_ts_006.xf \
  --tlt lam8_ts_006.tlt \
  --validate
```

Re-import / update an existing tomogram with `--update-existing` (or `--overwrite`).

---

## 2. Import picked coordinates

The key choices are the **coordinate system** (`rec-voxel` if you picked on the
reconstructed tomogram), the **indexing** (`zero-based` for napari, `one-based`
for IMOD models), and the **ROI radius** (`--roi-radius-angst`).

### 2a. Single point (quick test)

```bash
tomojanas-import particles \
  --project first \
  --tomo-name lam8_ts_006 \
  --input-single-point 512,401,188 \
  --coordinate-system rec-voxel \
  --indexing zero-based \
  --axis-order xyz \
  --roi-radius-angst 150 \
  --validate
```

> Re-running with another point **adds** `P000002`, `P000003`, … (particle IDs
> auto-increment from the highest existing one). To overwrite a specific
> particle, pass `--particle-name P000001`.

### 2b. From CSV (columns `x,y,z`)

```bash
tomojanas-import particles \
  --project first \
  --tomo-name lam8_ts_006 \
  --input picks.csv \
  --format csv \
  --coordinate-system rec-voxel \
  --indexing zero-based \
  --roi-radius-angst 150 \
  --validate
```

### 2c. From napari (Points layer CSV: `axis-0/1/2 = z,y,x`, handled automatically)

```bash
tomojanas-import particles \
  --project first \
  --tomo-name lam8_ts_006 \
  --input napari_points.csv \
  --format napari \
  --coordinate-system rec-voxel \
  --indexing zero-based \
  --roi-radius-angst 150 \
  --validate
```

### 2d. From IMOD (`model2point` → text X Y Z, one-based)

```bash
# first: model2point -input picks.mod -output picks.txt
tomojanas-import particles \
  --project first \
  --tomo-name lam8_ts_006 \
  --input picks.txt \
  --format csv \
  --coordinate-system rec-voxel \
  --indexing one-based \
  --axis-order xyz \
  --roi-radius-angst 150 \
  --validate
```

---

## 3. Extract physical sub-volume crops

Add `--write-rec-crops` to write a cubic sub-volume around the spherical ROI
into `individual_particles_recs/P*_rec.mrc`:

```bash
tomojanas-import particles \
  --project first \
  --tomo-name lam8_ts_006 \
  --input-single-point 161,678,170 \
  --coordinate-system rec-voxel \
  --indexing zero-based \
  --roi-radius-angst 150 \
  --write-rec-crops \
  --validate
```

Useful crop options:

| Option | Effect |
|--------|--------|
| `--crop-storage-box-size N` | Force the cubic box size (otherwise derived from the ROI radius) |
| `--apply-spherical-mask` | Zero the voxels outside the spherical ROI (box stays a container) |
| `--crop-outside-policy {error,pad,partial,skip}` | Behaviour when the ROI extends past the tomogram (default `partial`) |
| `--crop-padding-angst` / `--crop-padding-voxel` | Padding added around the ROI |

### 3b. Create sub-volumes for already-imported particles

If you imported particles **without** `--write-rec-crops`, you can backfill the
3D sub-volumes afterwards (no re-import needed) with `status --create-volume`.
It reads the coordinate and ROI radius from each `P*.star`, writes the missing
`individual_particles_recs/P*_rec.mrc`, and adds the crop block to the
`P*.star`:

```bash
tomojanas-import status \
  --project first \
  --tomo-name lam8_ts_006 \
  --create-volume
```

Options mirror the import crop options (`--crop-storage-box-size`,
`--apply-spherical-mask`, `--crop-outside-policy`, `--crop-padding-*`,
`--crop-pad-value`). Use `--overwrite-crops` to regenerate existing crops.

---

## 4. Import CTF metadata (no pixel modification)

```bash
tomojanas-import ctf \
  --project first \
  --tomo-name lam8_ts_006 \
  --ctf-source auto \
  --do-not-premultiply \
  --validate
```

Defocus is converted to Ångström; the detected unit/source decision is recorded
in `imod_settings/ctf_settings.star`. CTF import never alters image pixels.

---

## 5. Validate the project

```bash
tomojanas-import validate \
  --project first \
  --tomo-name lam8_ts_006 \
  --strict
```

Reports are written to `logs/validation_log.{star,json,md}`. `--strict` returns
a non-zero exit code if any critical check fails.

---

## 6. Inspect & reconcile the project (status)

List what is registered vs. what exists on disk, and flag inconsistencies:

```bash
tomojanas-import status --project first --tomo-name lam8_ts_006
```

**Deleting a particle** is best done by hand from the OS, then reconcile:

```bash
# delete the particle files
rm first/tilt_series/lam8_ts_006/individual_particles/P000001.star
rm first/tilt_series/lam8_ts_006/individual_particles_recs/P000001_rec.mrc

# rebuild particles_all.star from the P*.star files that remain
tomojanas-import status --project first --sync
```

---

## 7. Replay commands in a new project (command log)

Every `tomojanas-import` invocation is appended (with exit status) to:

- `logs/commands.jsonl` — machine-readable (timestamp, full argv, exit status, cwd);
- `logs/commands.sh` — a replayable shell script.

To reproduce a project elsewhere, edit the `--project` value in `commands.sh`
and run it:

```bash
sed 's#--project first#--project second#g' first/logs/commands.sh > replay.sh
bash replay.sh
```

---

## Troubleshooting: the sub-volume is wrong / in the wrong place

The crop indexes the reconstructed tomogram as `rec[z, y, x]`, so the
`(x, y, z)` you pass **must match the rec file's `(nx, ny, nz)` axes**. The
most common cause of a wrong crop is a swapped axis order from the picking tool:

| Picking tool | Coordinate order | Indexing | Use |
|--------------|------------------|----------|-----|
| **napari** (Points layer) | `(z, y, x)` = axis-0,1,2 | zero-based | `--format napari` (CSV) **or** `--axis-order zyx` (single point) |
| **IMOD 3dmod**, flipped tomogram (typical) | `(x, z, y)` — Y/Z swapped vs file | one-based | `--axis-order xzy --indexing one-based` |
| **IMOD 3dmod**, non-flipped | `(x, y, z)` | one-based | `--axis-order xyz --indexing one-based` |
| generic CSV `x,y,z` | `(x, y, z)` | as produced | `--axis-order xyz` |

> **IMOD note:** IMOD reconstructions are usually stored "flipped" (the file's
> Y axis is the depth/thickness). 3dmod displays them in the natural
> orientation, so the X/Y/Z you read from the **Zap window** are *display*
> coordinates with **Y and Z swapped** relative to the file. Pasting them
> directly therefore needs `--axis-order xzy`. Also, 3dmod coordinates are
> **1-based**, so use `--indexing one-based`. Confirm with the `[crop]`
> diagnostic that no coordinate is out of range.

All six axis-order permutations are accepted: `xyz, xzy, yxz, yzx, zxy, zyx`.

When you run a crop, tomoJANAS prints a diagnostic line, e.g.:

```
[crop] P000001: rec (nx,ny,nz)=(960,928,300) pixel=10.480 A; picked voxel (x,y,z)=(161,678,170)
```

Check that each picked value falls inside its axis range. If you see a warning
like `picked voxel outside tomogram on Z=678∉[0,300)`, your X and Z are almost
certainly swapped — re-run with the correct `--axis-order` (e.g. `zyx` for
napari). The same diagnostic is written to `logs/tomojanas_import.log`.

> Note: coordinates are always in **voxels (pixels)** of the reconstructed
> tomogram when using `--coordinate-system rec-voxel`. That part is correct;
> it is the *order* of the three numbers that usually needs fixing.

---

## 8. Legacy quick extraction (no project)

The original command still works for ad-hoc sub-volume extraction:

```bash
tomojanas import_volume \
  --tomo tomogram.mrc --coords picks.csv --format napari --box 64
```

In the new architecture the ROI is a sphere/circle and `--box` is only a storage
container size.
