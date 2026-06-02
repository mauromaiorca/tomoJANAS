[Repository home](../../README.md) · [Docs home](../index.md) · [Installation](../installation.md) · [Quick start](../quick-start.md) · [CLI reference](../reference/cli.md) · [Troubleshooting](../troubleshooting.md)

---

# Iterative particle selection

Particle selection ranks particles by their Structural Cross-correlation Index (SCI) and iteratively determines the subset that maximises mean local resolution within a defined mask.

## How it works

Each iteration performs the following steps:

1. **Randomise half-map assignments** — particles are randomly split into two half-sets for independent reconstruction.
2. **Reconstruct bootstrap half-maps** — 3D reconstruction from the current particle set.
3. **Score particles** — each particle is scored against the reference volume using SCI (3D-to-2D projection at the particle's Euler angles, followed by masked cross-correlation).
4. **Create ranked subsets** — particles are sorted by SCI and grouped into subsets at different cutoff points.
5. **Reconstruct subsets** — each subset is reconstructed independently.
6. **Evaluate local resolution** — local FSC is computed for each subset's half-maps within the mask.
7. **Select the best cutoff** — the subset with the best mean local resolution defines the new particle set.

The loop repeats until no further improvement is observed (controlled by `num_unimproving_iterations_for_early_termination`, default 3).

## Setting up a session

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

> **Note on `--noExternalPrograms` and `--gpu`:** these flags affect **only the reconstruction and local resolution steps**; particle scoring always runs on CPU (controlled by `--mpi`). GPU acceleration applies only to reconstruction. Choose based on your hardware:
>
> - **With GPU(s) (recommended):** `--noExternalPrograms --gpu 0 1` uses two GPUs, one per half-map, for the best throughput. `--gpu 0` uses a single GPU.
> - **No GPU but RELION available:** omit `--noExternalPrograms` — JANAS will call RELION's MPI-based reconstruction and `relion_postprocess`, which on a CPU-only machine is typically faster than JANAS's internal CPU reconstruction. Make sure RELION is installed and accessible on your `PATH`.
> - **No GPU, no RELION:** `--noExternalPrograms` without `--gpu` runs JANAS's internal CPU reconstruction; it works but is slower on large datasets.
>
> You can also choose to use RELION for reconstruction even when a GPU is available — simply omit `--noExternalPrograms`.

This creates:

```
my_selection/
├── my_selection_run.sh
└── session_settings.toml
```

### Key parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--map`, `--map2` | Half-map reconstructions used as initial reference | Required |
| `--mask` | 3D mask defining the region of interest | Required |
| `--sigma` | Gaussian sigma for SCI calculation | 1 |
| `--mpi` | Number of CPU processes for scoring and local resolution | System-dependent |
| `--bootstrap` | Reconstruct initial half-maps before first scoring | Off |
| `--numRecs` | Number of subset reconstructions per iteration | 12 |
| `--maxSelections` | Maximum number of selection iterations | 14 |

You can edit `session_settings.toml` or the run script directly before launching.

## Running

Ensure the raw particle stack is accessible from the working directory (symlink if needed):

```bash
ln -s /path/to/particles.mrcs my_selection/
./my_selection/my_selection_run.sh
```

## Output

| File | Description |
|------|-------------|
| `reference_subset.star` | Selected particle subset |
| `overview.txt` | Per-iteration statistics |
| `_janas_select_iter_*/` | Per-iteration working directories |

## Reviewing results

View the optimisation summary:

```bash
janas_optimizer plotOverview --overview my_selection/overview.txt --plot
```

Check particle count in the final subset:

```bash
janas_app_starProcess --i my_selection/reference_subset.star --info
```

Inspect Euler angle coverage:

```bash
janas eulerHist --i my_selection/reference_subset.star
```

## GPU acceleration

By default, JANAS uses RELION for reconstruction. To use the internal GPU reconstructor instead:

```bash
--noExternalPrograms --gpu 0
```

This removes the dependency on RELION and typically provides 5-20x speedup for reconstruction. See [Computational requirements](../reference/computational-requirements.md) for hardware guidance.

---

[Back to documentation index](../index.md)
