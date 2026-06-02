[Repository home](../../README.md) · [Docs home](../index.md) · [Installation](../installation.md) · [Quick start](../quick-start.md) · [CLI reference](../reference/cli.md) · [Troubleshooting](../troubleshooting.md)

---

# EMPIAR-10308: Serotonin 5-HT1B-Go receptor complex

This tutorial demonstrates both JANAS workflows — 3D class reassignment followed by per-class particle selection — on the serotonin 5-HT1B-Go receptor complex dataset (EMPIAR-10308).

## Data preparation

Download the particle stacks from EMPIAR:

```bash
mkdir -p Particles
lftp ftp://ftp.ebi.ac.uk <<EOF
  mirror --verbose /empiar/world_availability/10308/data/Particles Particles
  quit
EOF
```

Download the reference files from [this link](https://syncandshare.desy.de/index.php/s/7YdCeNegBi5zcW9) and place them in your working directory. You need:

| File | Description | Particles |
|------|-------------|:---------:|
| `reference_Refined_rec.mrc` | Full-dataset reconstruction | 244,565 |
| `J958_003_volume_mapAligned.mrc` | Low-populated conformational state | 14,641 |
| `maskFull_dilatedClose.mrc` | Mask enclosing all densities | - |
| `reference_Refined_notationEMPIAR.star` | Particle metadata | 244,565 |

Create a single particle stack:

```bash
relion_stack_create \
    --i reference_Refined_notationEMPIAR.star \
    --o reference_Refined_stack \
    --one_by_one
```

This produces `reference_Refined_stack.star` and `reference_Refined_stack.mrcs`.

## Step 1: 3D class reassignment

Set up the classification session with two reference maps:

```bash
janas_session_manager classification_session \
    --name demo_JANAS \
    --particles reference_Refined_stack.star \
    --maps reference_Refined_rec.mrc J958_003_volume_mapAligned.mrc \
    --mask maskFull_dilatedClose.mrc \
    --mpi 85 \
    --noExternalPrograms --gpu 0 1
```

> **Note on `--noExternalPrograms` and `--gpu`:** these flags affect **only the reconstruction and local resolution steps**; particle scoring always runs on CPU (controlled by `--mpi`). GPU acceleration applies only to reconstruction. Choose based on your hardware:
>
> - **With GPU(s) (recommended):** `--noExternalPrograms --gpu 0 1` uses two GPUs, one per half-map, for the best throughput. `--gpu 0` uses a single GPU.
> - **No GPU but RELION available:** omit `--noExternalPrograms` — JANAS will call RELION's MPI-based reconstruction and `relion_postprocess`, which on a CPU-only machine is typically faster than JANAS's internal CPU reconstruction. Make sure RELION is installed and accessible on your `PATH`.
> - **No GPU, no RELION:** `--noExternalPrograms` without `--gpu` runs JANAS's internal CPU reconstruction; it works but is slower on large datasets.
>
> You can also choose to use RELION for reconstruction even when a GPU is available — simply omit `--noExternalPrograms`.

Run:

```bash
./demo_JANAS/demo_JANAS_run.sh
```

After completion, the output directory contains:

```
demo_JANAS/
├── equalized/
│   ├── J958_003_volume_map_equalized.mrc
│   └── reference_Refined_rec_equalized.mrc
├── final_classes/
│   ├── class_1.star          (216,884 particles)
│   ├── class_1_recH1.mrc
│   ├── class_1_recH2.mrc
│   ├── class_2.star          (27,681 particles)
│   ├── class_2_recH1.mrc
│   ├── class_2_recH2.mrc
│   └── demo_JANAS_classified.star
└── scored_classes/
    ├── demo_JANAS.csv
    └── *_scoredClass.star
```

Average the half-maps for visualisation:

```bash
cd demo_JANAS/final_classes
janas_utils clip average class_1_recH1.mrc class_1_recH2.mrc class_1_rec.mrc
janas_utils clip average class_2_recH1.mrc class_2_recH2.mrc class_2_rec.mrc
```

Check particle counts and Euler distributions:

```bash
janas_app_starProcess --i class_1.star --info
janas_app_starProcess --i class_2.star --info
janas eulerHist --i class_1.star
janas eulerHist --i class_2.star
```

## Step 2: Per-class particle selection

After class reassignment, refine each class by selecting the particles that contribute most to local map quality.

### Class 1

```bash
cd demo_JANAS/final_classes

janas_session_manager new_select_session \
    --name class1_selection \
    --particles class_1.star \
    --map class_1_recH1.mrc \
    --map2 class_1_recH2.mrc \
    --mask ../../maskFull_dilatedClose.mrc \
    --sigma 1 \
    --mpi 85 --bootstrap \
    --numRecs 12 --maxSelections 14 \
    --noExternalPrograms --gpu 0 1
```

Ensure the raw particle stack is accessible:

```bash
ln -s ../../reference_Refined_stack.mrcs .
```

Run:

```bash
./class1_selection/class1_selection_run.sh
```

The selected particles are in `class1_selection/reference_subset.star`.

View the optimisation overview:

```bash
janas_optimizer plotOverview --overview class1_selection/overview.txt --plot
```

### Class 2

The procedure is identical, substituting class 2 files:

```bash
janas_session_manager new_select_session \
    --name class2_selection \
    --particles class_2.star \
    --map class_2_recH1.mrc \
    --map2 class_2_recH2.mrc \
    --mask ../../maskFull_dilatedClose.mrc \
    --sigma 1 \
    --mpi 85 --bootstrap \
    --numRecs 12 --maxSelections 14 \
    --noExternalPrograms --gpu 0 1

./class2_selection/class2_selection_run.sh
```

## Summary

This tutorial demonstrated:

1. **Class reassignment** with two reference maps, separating 244k particles into two structurally distinct classes (217k and 28k particles).
2. **Per-class particle selection** to further refine each class by removing particles that degrade local map quality.

---

[Back to documentation index](../index.md)
