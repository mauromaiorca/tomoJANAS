[Repository home](../../README.md) · [Docs home](../index.md) · [Installation](../installation.md) · [Quick start](../quick-start.md) · [CLI reference](../reference/cli.md) · [Troubleshooting](../troubleshooting.md)

---

# 3D class reassignment

Class reassignment scores particles against multiple reference maps and assigns each particle to the class with the highest Structural Cross-correlation Index (SCI).

## How it works

1. **Equalise amplitudes** — compute the average amplitude spectrum across all K reference maps, then replace each map's amplitudes with the average while preserving phases. This ensures scoring differences reflect structural variation, not amplitude differences.
2. **Score particles against each map** — each particle is scored K times, once per reference map.
3. **Assign particles** — each particle is assigned to the class whose reference map produced the highest SCI score.
4. **Reconstruct per-class volumes** — independent half-map reconstructions for each class.

This is a single-pass workflow (no iterations).

## Setting up a session

```bash
janas_session_manager classification_session \
    --name reclassify \
    --particles particles.star \
    --maps class1.mrc class2.mrc class3.mrc \
    --mask mask.mrc \
    --mpi 40 \
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
reclassify/
├── reclassify_run.sh
├── final_classes/
├── scored_classes/
└── session_classification_settings.txt
```

### Key parameters

| Parameter | Description |
|-----------|-------------|
| `--maps` | Two or more reference maps representing different conformations |
| `--mask` | 3D mask enclosing all relevant densities |
| `--mpi` | Number of CPU processes for scoring |

## Running

```bash
./reclassify/reclassify_run.sh
```

## Output

After completion:

```
reclassify/
├── equalized/
│   ├── class1_equalized.mrc
│   └── class2_equalized.mrc
├── final_classes/
│   ├── class_1.star
│   ├── class_1_recH1.mrc
│   ├── class_1_recH2.mrc
│   ├── class_2.star
│   ├── class_2_recH1.mrc
│   ├── class_2_recH2.mrc
│   └── reclassify_classified.star
└── scored_classes/
    ├── reclassify.csv
    └── *_scoredClass.star
```

| File | Description |
|------|-------------|
| `final_classes/class_N.star` | Particles assigned to class N |
| `final_classes/class_N_recH1.mrc`, `class_N_recH2.mrc` | Half-map reconstructions for class N |
| `final_classes/*_classified.star` | All particles with class assignments |
| `scored_classes/*.csv` | Per-particle SCI scores for all maps |

## Post-classification steps

Average the half-maps for visualisation:

```bash
janas_utils clip average class_1_recH1.mrc class_1_recH2.mrc class_1_rec.mrc
```

Check particle counts:

```bash
janas_app_starProcess --i final_classes/class_1.star --info
```

View Euler angle distributions:

```bash
janas eulerHist --i final_classes/class_1.star
```

## Combining with particle selection

After class reassignment, each class-specific stack may still contain particles that degrade local map quality. Run iterative particle selection on each class separately to further refine the reconstructions. See the [EMPIAR-10308 tutorial](../examples/empiar-10308.md) for a worked example of this combined workflow.

---

[Back to documentation index](../index.md)
