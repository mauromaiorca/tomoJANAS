[Repository home](../../README.md) · [Docs home](../index.md) · [Installation](../installation.md) · [Quick start](../quick-start.md) · [CLI reference](../reference/cli.md) · [Troubleshooting](../troubleshooting.md)

---

# CryoSPARC integration

JANAS operates on RELION-style STAR files. This page covers how to import particle data from cryoSPARC and how to use cryoSPARC's local NU-refinement within a JANAS session.

## Importing particles from cryoSPARC

### Using `csparc2star-stack` (recommended)

From JANAS v0.1.3.2 onwards, a single command converts a cryoSPARC `.cs` file into a RELION STAR file and assembles the particle stack:

```bash
janas_utils csparc2star-stack \
    /path/to/cryosparc/P1/J80/cryosparc_P1_J80_002_particles.cs \
    outStack
```

This produces:

- `outStack.star` — particle metadata with `_rlnImageName` pointing to the new stack. An additional column `_janas_csparc_rlnImageName` records the original cryoSPARC image names for provenance.
- `outStack.mrcs` — assembled 2D particle stack.

The resulting files can be used directly in JANAS or imported into RELION.

!!! tip
    Use the `pwd` command to get the absolute path to your cryoSPARC job directory.

### Legacy method (before v0.1.3.2)

If you prefer scripting, the manual procedure involves:

1. **Convert `.cs` to STAR** using pyem or `janas_utils csparc2star`:

    ```bash
    janas_utils csparc2star particles.cs particles.star
    ```

2. **Create `.mrcs` symlinks** — cryoSPARC uses `.mrc` for stacks; RELION requires `.mrcs`:

    ```bash
    for f in /path/to/extract/*.mrc; do
        ln -s "$f" "$(dirname "$f")/$(basename "$f" .mrc).mrcs"
    done
    ```

3. **Update the STAR file** to reference `.mrcs` filenames:

    ```bash
    sed -i 's/_particles\.mrc /_particles.mrcs /g' particles.star
    ```

4. **Create a single stack** (optional but recommended):

    ```bash
    relion_stack_create --i particles.star --o full_stack --one_by_one
    ```

## Updating metadata from cryoSPARC

To merge updated parameters (angles, shifts, CTF) from a cryoSPARC job into an existing STAR file:

```bash
janas_utils update_from_csparc input.cs input.star output.star
```

This preserves the `data_optics` block and updates only columns present in the STAR header. Supported columns: `_rlnAngleRot`, `_rlnAngleTilt`, `_rlnAnglePsi`, `_rlnOriginXAngst`, `_rlnOriginYAngst`, `_rlnDefocusU`, `_rlnDefocusV`, `_rlnDefocusAngle`, `_rlnPhaseShift`, `_rlnCtfBfactor`, `_rlnOpticsGroup`, `_rlnRandomSubset`, `_rlnClassNumber`.

The `.cs` file and STAR file must contain the same number of particles in the same order. The command will fail if they do not match.

## CryoSPARC local NU-refinement in JANAS sessions

JANAS can dispatch cryoSPARC's local non-uniform refinement between selection iterations.

### Setup

```bash
janas csparc_setup \
    --license-id YOUR_LICENSE_ID \
    --host YOUR_HOST \
    --base-port 39000 \
    --email YOUR_EMAIL

source ~/.janas/cryosparc_env.sh
```

### Running with cryoSPARC refinement

Add `--cs_project` to the session setup, or use `janas_utils csparc_localnurefinement` directly:

```bash
janas_utils csparc_localnurefinement \
    --particle-dir . \
    --project P31 --workspace W1 --lane default \
    --sym C1 \
    --ref "$(pwd)/class_1_rec.mrc" \
    --mask mask.mrc \
    --resplit class_1.star class_1_updatedLNU
```

This adds 15-60 minutes per iteration depending on particle count and GPU hardware. It is skipped automatically when the previous iteration showed no improvement.

---

[Back to documentation index](../index.md)
