[Repository home](../README.md) · [Documentation index](index.md) · [Installation](installation.md) · [CLI reference](reference/cli.md) · [Troubleshooting](troubleshooting.md)

---

# Accessory utilities

JANAS bundles a number of stand-alone utilities for common cryo-EM map and stack operations. They are invoked through `janas_utils <subcommand>` (except `janas_reconstructor`, which is its own command). Use `janas_utils <subcommand> --help` for the full option list.

## sigma_estimate

Estimate a suitable Gaussian sigma for SCI scoring from a pair of half-maps. The value can then be passed to `--sigma` in session-manager and scoring commands.

```bash
janas_utils sigma_estimate halfmap1.mrc halfmap2.mrc [--mask mask.mrc]
```

For the full derivation and the rationale behind the automated estimator, see the dedicated [sigma_estimate page](sigma_estimate.md).

## compare_maps

Compare two 3D maps using cross-correlation and related similarity measures (optionally inside a mask).

```bash
janas_utils compare_maps map1.mrc map2.mrc [--mask mask.mrc]
```

## backmap_stars

Inverse companion of `create_stack`. If a downstream STAR file (refined, classified, etc.) lost the `_janas_source_rlnImageName` provenance column added during stack creation, this utility restores the original source `_rlnImageName` by joining the downstream STAR against the stack-generation STAR. All refined metadata (angles, defocus, class numbers, ...) and row order are preserved.

```bash
janas_utils backmap_stars \
    --processed run_it025_data.star \
    --mapping   EMPIAR_12707_stack.star \
    --output    run_it025_data_backmapped.star
```

By default an audit column `_janas_stack_rlnImageName` is added with the previous stack-based name so the rewrite is reversible. Pass `--stack-reference-tag ""` to skip it. Use `--no-strict` to leave unmapped rows unchanged instead of failing.

## csparc2star-stack

Convert a CryoSPARC `.cs` file into a RELION-style STAR file and assemble a consolidated `.mrcs` particle stack from the original cryoSPARC stack locations.

```bash
janas_utils csparc2star-stack /path/to/particles.cs output_prefix
```

Produces `output_prefix.star` and `output_prefix.mrcs`. See also the [CryoSPARC integration guide](workflows/cryosparc.md).

## clip blur

Gaussian blur a 3D volume with a sigma specified in Ångström.

```bash
janas_utils clip blur in.mrc out.mrc SIGMA_A
```

## clip bfac

B-factor weighting (sharpening) of a 3D volume. Automatic mode estimates the B-factor from the two half-maps; user-driven mode applies a specified B-value.

```bash
# Automatic
janas_utils clip bfac half1.mrc half2.mrc out.mrc

# User-driven
janas_utils clip bfac half1.mrc half2.mrc out.mrc BVALUE
```

## fsc

Compute the Fourier Shell Correlation (FSC) between one or more half-map pairs.

```bash
janas_utils fsc half1.mrc half2.mrc
```

## locres

Compute a local-resolution map from a pair of half-maps. Writes `*_locres.mrc` and auxiliary files.

```bash
janas_utils locres half1.mrc half2.mrc [--mask mask.mrc]
```

For bulk processing of many half-map pairs (used internally by the selection workflow), see `janas_utils locresBulk`.

## project_map

Project a 3D reference map at each particle pose listed in a STAR file and write the resulting 2D reprojections to disk.

```bash
janas_utils project_map --i particles.star --map reference.mrc --o reprojections.mrcs
```

## janas_reconstructor

Internal 3D reconstruction from scored particles. Supports GPU acceleration (PyTorch CUDA, mini-batched) or CPU multiprocessing. Used automatically by the session-manager workflows when `--noExternalPrograms` is set.

Single-volume reconstruction:

```bash
janas_reconstructor \
    --i particles.star \
    --o output.mrc \
    --gpu 0
```

Independent half-map reconstruction (for FSC and local-resolution evaluation):

```bash
janas_reconstructor \
    --i particles.star \
    --o output \
    --gpu 0 \
    --subset 1 2
```

This produces `output_recH1.mrc` and `output_recH2.mrc` from the two particle subsets defined by `_rlnRandomSubset`.

> **GPU vs CPU:** GPU is strongly recommended. Without a GPU (omit `--gpu`), reconstruction falls back to CPU multiprocessing, which is **very slow** on medium-to-large datasets. On CPU-only machines, calling RELION (omit `--noExternalPrograms` in the session-manager workflows) is typically faster than the internal CPU reconstructor.

See the [CLI reference](reference/cli.md#janas_reconstructor) for the full option list.

---

[Back to documentation index](index.md)
