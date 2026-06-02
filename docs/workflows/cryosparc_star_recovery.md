# CryoSPARC particle STAR recovery

CryoSPARC stores particle metadata in `.cs` files, using a binary header format that is different from the text-based STAR format used by RELION and many  tools, including JANAS. A `.cs` file cannot be inspected or edited as directly as a `.star` file, and conversion is often required before the metadata can be used in external workflows.

Several tools can already convert CryoSPARC `.cs` files into STAR files for readability and interoperability. However, conversion alone is not always sufficient. The particle metadata also contains file references, especially `_rlnImageName`, that define how raw 2D particles, particle stacks or micrographs are located on disk. CryoSPARC may rewrite this logic for internal efficiency, for example by creating symbolic links to original files, placing them inside job-specific directories, or adding long numeric prefixes to linked stack names. These changes are valid inside a CryoSPARC project, but they can make exported STAR files difficult to reuse outside that project.

This workflow therefore does more than convert metadata. It helps restore a usable filesystem logic after CryoSPARC export. `csparc2star` can clean or replace CryoSPARC-generated paths during conversion, so that `_rlnImageName` points to files that exist in the user-defined processing environment. When particles have been consolidated into a new stack, as commonly done in JANAS to simplify processing and export, the same conversion step can also prepare the STAR file for provenance recovery. If a JANAS stack-generation STAR file is available, `backmap_stars` can then restore the original source-particle references.

This section covers two related situations.

First, CryoSPARC may export `_rlnImageName` values that contain internal job paths and long numeric prefixes before the stack filename, for example:

```text
000001@J41/extract/008878123945933052272_HTT_46Q_0000_Aug06_00.39.32_aligned_DW_particles.mrc
```

In this case, `csparc2star` can clean the exported image names during conversion.

Second, particles may have been consolidated into a single `.mrcs` stack using JANAS, processed in CryoSPARC, and exported back as metadata. The exported STAR file then points to the consolidated stack rather than to the original source particle stacks. If provenance was recorded during stack creation, `backmap_stars` can restore the original `_rlnImageName` values.

---

## Workflow overview

```text
CryoSPARC .cs
     │
     │  (1) csparc2star
     │      Convert to STAR
     │      Optionally clean or replace CryoSPARC paths
     ▼
STAR file with corrected _rlnImageName
     │
     ├── Use directly in RELION/JANAS
     │
     └── Optional:
          (2) backmap_stars
              Restore original source-particle names
              when the particles came from a JANAS-generated stack
```

In most CryoSPARC export cases, Step 1 is sufficient. Step 2 is needed only when the STAR file points to a consolidated stack and you want to recover the original source-particle references.

---

## Step 1 — Convert the CryoSPARC `.cs` file to STAR

CryoSPARC exports particle metadata in its own `.cs` binary format, which RELION and JANAS cannot read directly. `csparc2star` converts the `.cs` file into a RELION-style `.star` file containing particle image references, orientations, shifts, CTF parameters and class assignments.

Basic conversion:

```bash
janas_utils csparc2star \
  class_J1003_4028particles/J1003_003_particles.cs \
  class_J1003_4028particles/J1003_003_particles.star
```

A typical `_rlnImageName` entry exported from CryoSPARC may look like this:

```text
000001@J41/extract/008878123945933052272_HTT_46Q_0000_Aug06_00.39.32_aligned_DW_particles.mrc
```

This entry contains three components:

```text
000001@J41/extract/008878123945933052272_HTT_46Q_0000_Aug06_00.39.32_aligned_DW_particles.mrc
│      │           │                    │                                      │
│      │           │                    │                                      └── stack extension
│      │           │                    └── actual stack name with _particles suffix
│      │           └── CryoSPARC numeric prefix
│      └── CryoSPARC job path
└── particle index inside the stack
```

The particle index before `@` must be preserved. The options below modify only the path and filename after `@`.

---

## Step 2 — Clean CryoSPARC paths during `csparc2star` conversion

`csparc2star` can clean `_rlnImageName` while writing the STAR file. This avoids a separate STAR editing step in many cases.

### Available options

| Option | Effect |
|------|------|
| `--clean_path` | Removes the CryoSPARC directory from the stack name. |
| `--clean_prefix` | Removes the long numeric CryoSPARC prefix before the filename. |
| `--clean_suffix` | Removes a terminal `_particles` suffix before `.mrc` or `.mrcs`. |
| `--fix_path PATH` | Replaces the CryoSPARC directory with `PATH`. |

These options can be combined. `--fix_path` takes precedence over `--clean_path`, because replacing the path and removing the path are mutually exclusive operations.

---

## Example A — Remove only the CryoSPARC path

Command:

```bash
janas_utils csparc2star \
  particles.cs \
  particles.star \
  --clean_path
```

Input `_rlnImageName`:

```text
000001@J41/extract/008878123945933052272_HTT_46Q_0000_Aug06_00.39.32_aligned_DW_particles.mrc
```

Output `_rlnImageName`:

```text
000001@008878123945933052272_HTT_46Q_0000_Aug06_00.39.32_aligned_DW_particles.mrc
```

Use this when the stack file is in the same directory where the STAR file will be used, but you still want to keep the CryoSPARC-generated filename.

---

## Example B — Remove the CryoSPARC numeric prefix

Command:

```bash
janas_utils csparc2star \
  particles.cs \
  particles.star \
  --clean_prefix
```

Input:

```text
000001@J41/extract/008878123945933052272_HTT_46Q_0000_Aug06_00.39.32_aligned_DW_particles.mrc
```

Output:

```text
000001@J41/extract/HTT_46Q_0000_Aug06_00.39.32_aligned_DW_particles.mrc
```

Use this when the CryoSPARC path is still valid, but the real stack file does not contain the long numeric prefix added by CryoSPARC.

---

## Example C — Remove the `_particles` suffix

Command:

```bash
janas_utils csparc2star \
  particles.cs \
  particles.star \
  --clean_suffix
```

Input:

```text
000001@J41/extract/008878123945933052272_HTT_46Q_0000_Aug06_00.39.32_aligned_DW_particles.mrc
```

Output:

```text
000001@J41/extract/008878123945933052272_HTT_46Q_0000_Aug06_00.39.32_aligned_DW.mrc
```

Use this when CryoSPARC refers to a file ending in `_particles.mrc`, but the real file on disk does not include `_particles`.

---

## Example D — Replace the CryoSPARC path with a usable path

Command:

```bash
janas_utils csparc2star \
  particles.cs \
  particles.star \
  --fix_path data/upload
```

Input:

```text
000001@J41/extract/008878123945933052272_HTT_46Q_0000_Aug06_00.39.32_aligned_DW_particles.mrc
```

Output:

```text
000001@data/upload/008878123945933052272_HTT_46Q_0000_Aug06_00.39.32_aligned_DW_particles.mrc
```

Use this when the stack files are stored in another directory relative to the final STAR file.

---

## Example E — Recommended command for cleaned CryoSPARC exports

If the real stack file is:

```text
data/upload/HTT_46Q_0000_Aug06_00.39.32_aligned_DW.mrc
```

and the CryoSPARC `_rlnImageName` is:

```text
000001@J41/extract/008878123945933052272_HTT_46Q_0000_Aug06_00.39.32_aligned_DW_particles.mrc
```

use:

```bash
janas_utils csparc2star \
  particles.cs \
  particles_clean.star \
  --fix_path data/upload \
  --clean_prefix \
  --clean_suffix
```

The output `_rlnImageName` will be:

```text
000001@data/upload/HTT_46Q_0000_Aug06_00.39.32_aligned_DW.mrc
```

This is the preferred route when the problem is only the CryoSPARC-generated path or filename. It produces a STAR file that can be used directly in RELION or JANAS without running a separate path-renaming command.

---

## Step 3 — alternative: rename the stack path after conversion

As an alternative, you can use `emprove_app_starProcess --stackRename` after `csparc2star` conversion. This is still valid, but it is less flexible than the `csparc2star` options because it replaces the stack name globally after the STAR file has already been written.

Conversion:

```bash
janas_utils csparc2star \
  class_J1003_4028particles/J1003_003_particles.cs \
  class_J1003_4028particles/J1003_003_particles.star
```

Path replacement:

```bash
emprove_app_starProcess \
  --i class_J1003_4028particles/J1003_003_particles.star \
  --stackRename EMPIAR_12707_stack.mrcs \
  --o class_J1003_4028particles/J1003_003_particlesStack.star
```

Before:

```text
000022@J1001/imported/008248594774147691034_EMPIAR_12707_stack.mrcs
```

After:

```text
000022@EMPIAR_12707_stack.mrcs
```

This remains useful when you have an already converted STAR file and only need to replace the whole stack reference. 

---

## Step 4 — Optional STAR recovery with `backmap_stars`

Use `backmap_stars` only when the CryoSPARC STAR file points to a consolidated `.mrcs` stack and you want to restore the original source-particle names.

This applies to workflows where a JANAS-generated stack was created from an original STAR file. In the original STAR file, each particle may have had an `_rlnImageName` such as:

```text
009619@J1149/restack/batch_6_restacked.mrc
```

During stack creation this was replaced by a consolidated-stack reference such as:

```text
000001@EMPIAR_12707_stack.mrcs
```

If provenance tracking was enabled during stack creation, the stack-generation STAR file records the relation between the consolidated-stack particle and the original source particle:

```text
_rlnImageName                    _janas_source_rlnImageName
000001@EMPIAR_12707_stack.mrcs   009619@J1149/restack/batch_6_restacked.mrc
000002@EMPIAR_12707_stack.mrcs   009620@J1149/restack/batch_6_restacked.mrc
000003@EMPIAR_12707_stack.mrcs   009621@J1149/restack/batch_6_restacked.mrc
```

CryoSPARC or other downstream tools may discard `_janas_source_rlnImageName`. The processed STAR file then still points to the consolidated stack, and the original particle names are no longer present. `backmap_stars` restores them using the stack-generation STAR file as a lookup table.

Command:

```bash
janas_utils backmap_stars \
  --mapping ../original_stack/EMPIAR_12707_stack.star \
  --processed class_J1003_4028particles/J1003_003_particlesStack.star \
  --output class_J1003_4028particles/J1003_003_particlesOriginalParticles.star
```

Arguments:

| Argument | Meaning |
|------|------|
| `--mapping` | STAR file produced when the consolidated stack was created. It must contain both `_rlnImageName` and `_janas_source_rlnImageName`. |
| `--processed` | STAR file to fix. Its `_rlnImageName` values point to the consolidated stack. |
| `--output` | Corrected STAR file, with `_rlnImageName` restored to the original source-particle references. |

Example replacement:

```text
000022@EMPIAR_12707_stack.mrcs
```

becomes:

```text
009640@J1149/restack/batch_6_restacked.mrc
```

By default, the command also writes an audit column named `_janas_stack_rlnImageName`. This column records the previous consolidated-stack reference. The output is therefore traceable: `_rlnImageName` points back to the original source particle, while `_janas_stack_rlnImageName` records the stack particle used for the mapping.

---

## Strict mapping

By default, `backmap_stars` runs in strict mode. Every `_rlnImageName` in the processed STAR file must be found in the mapping STAR file. If a particle cannot be mapped, the command stops with an error. This prevents the generation of a partially corrected STAR file.

For debugging, strict mode can be disabled:

```bash
janas_utils backmap_stars \
  --mapping ../original_stack/EMPIAR_12707_stack.star \
  --processed class_J1003_4028particles/J1003_003_particlesStack.star \
  --output class_J1003_4028particles/J1003_003_particlesOriginalParticles.star \
  --no-strict
```

With `--no-strict`, unmapped particles are left unchanged. Use this only to inspect problematic datasets, not for routine processing.

---

## Complete command sequences

### A. Direct CryoSPARC path correction during conversion

Use this when the only problem is the CryoSPARC-generated path or filename:

```bash
janas_utils csparc2star \
  particles.cs \
  particles_clean.star \
  --fix_path data/upload \
  --clean_prefix \
  --clean_suffix
```

Final STAR file:

```text
particles_clean.star
```

This file can be used directly in RELION or JANAS if the stack files exist at the corrected paths.

---

### B. CryoSPARC conversion followed by original-particle recovery

Use this when the particles came from a JANAS-generated consolidated stack and the original `_rlnImageName` values must be restored:

```bash
# 1. Convert CryoSPARC particle metadata to STAR
#    and make the consolidated-stack path match the mapping STAR.
janas_utils csparc2star \
  class_J1003_4028particles/J1003_003_particles.cs \
  class_J1003_4028particles/J1003_003_particlesStack.star \
  --clean_prefix \
  --clean_suffix \
  --fix_path .

# 2. Restore the original source-particle image names.
janas_utils backmap_stars \
  --mapping ../original_stack/EMPIAR_12707_stack.star \
  --processed class_J1003_4028particles/J1003_003_particlesStack.star \
  --output class_J1003_4028particles/J1003_003_particlesOriginalParticles.star
```

The final STAR file to use in subsequent RELION or JANAS steps is:

```text
class_J1003_4028particles/J1003_003_particlesOriginalParticles.star
```

For `backmap_stars` to work, the `_rlnImageName` values in the processed STAR file must match the `_rlnImageName` values in the mapping STAR file. The `csparc2star` cleaning options are therefore useful not only for path correction, but also for STAR recovery: they allow the CryoSPARC-converted STAR file to be normalised to the same stack naming convention used in the mapping STAR.

---

## Summary

| Problem | Recommended command |
|------|------|
| CryoSPARC path should be removed | `csparc2star --clean_path` |
| CryoSPARC numeric prefix should be removed | `csparc2star --clean_prefix` |
| Terminal `_particles` should be removed from the filename | `csparc2star --clean_suffix` |
| CryoSPARC directory should be replaced with another directory | `csparc2star --fix_path PATH` |
| CryoSPARC STAR should be converted and cleaned in one step | `csparc2star --fix_path PATH --clean_prefix --clean_suffix` |
| Existing converted STAR should have its stack name replaced globally | `emprove_app_starProcess --stackRename` |
| Consolidated-stack references should be restored to original source particles | `backmap_stars` |

In short, `csparc2star` fixes the CryoSPARC-generated `_rlnImageName` during conversion. `backmap_stars` restores original source-particle provenance when a JANAS-generated mapping STAR is available.

---

> **Note on command naming.** During the transition from the former name EMPROVE to JANAS, some commands still carry the `emprove_` prefix, for example `emprove_app_starProcess`, as backward-compatible aliases. Equivalent `janas_`-prefixed commands are being introduced. Either form can be used during the transition.
