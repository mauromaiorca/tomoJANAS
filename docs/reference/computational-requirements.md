[Repository home](../../README.md) · [Docs home](../index.md) · [Installation](../installation.md) · [Quick start](../quick-start.md) · [CLI reference](cli.md) · [Troubleshooting](../troubleshooting.md)

---

# Computational requirements

## Hardware requirements

| Component | Minimum | Recommended | Notes |
|-----------|---------|-------------|-------|
| CPU cores | 4 | 32-80 | Scoring and local resolution scale near-linearly with cores |
| RAM | 16 GB | 64-128 GB | Proportional to particle count and box size |
| GPU | Not required | 1x NVIDIA A100 / V100 | 5-20x speedup for reconstruction |
| GPU VRAM | - | 16 GB+ | Required for box sizes above 300 px with CTF correction |
| Storage | 50 GB SSD | 500 GB+ (fast I/O) | Particle stacks dominate disk usage |

GPU is not required. JANAS runs entirely on CPU. GPU is strongly recommended when using `--noExternalPrograms` with datasets above 50k particles.

## Computational steps

### Particle scoring (`janas scoreParticles`)

- **Algorithm:** per-particle 3D-to-2D projection at the particle's Euler angles, followed by masked SCI cross-correlation
- **Complexity:** O(N x B^3 / P)
- **Parallelisation:** CPU multiprocessing, block partitioning across P workers
- **Scaling:** near-linear up to 40-80 cores, limited by I/O beyond that

### 3D reconstruction (`janas_reconstructor`)

- **Algorithm:** direct Fourier insertion with CTF correction, followed by 3D inverse FFT
- **Complexity:** O(N x B^2 log B)
- **Parallelisation:** GPU (PyTorch CUDA, mini-batching) or CPU multiprocessing
- **GPU options:** `--gpu 0` (single) or `--gpu 0 1` (multi-GPU); `--gpu-batch 20` (default)

### Local resolution (`janas_utils locresBulk`)

- **Algorithm:** sliding spherical window with local FSC computation
- **Complexity:** O(S x (V/d)^3 x R^3 log R / P)
- **Parallelisation:** CPU multiprocessing
- **Note:** most repeated step in selection (~15 calls per iteration)

### Amplitude equalisation

- **Complexity:** O(K x B^3 log B)
- **Run once** before classification scoring

## Memory estimates

| Component | B=200 | B=256 | B=320 | B=400 |
|-----------|:-----:|:-----:|:-----:|:-----:|
| Reference volume (float32) | 31 MB | 64 MB | 125 MB | 244 MB |
| GPU accumulator + weight | 92 MB | 192 MB | 375 MB | 732 MB |
| GPU particle batch (batch=20) | 3 MB | 5 MB | 8 MB | 12 MB |
| Two half-maps in memory | 61 MB | 128 MB | 250 MB | 488 MB |

**Particle stack on disk:**

| Particles | B=200 | B=256 | B=320 | B=400 |
|:---------:|:-----:|:-----:|:-----:|:-----:|
| 88k | 13 GB | 22 GB | 34 GB | 54 GB |
| 192k | 29 GB | 48 GB | 75 GB | 117 GB |
| 458k | 70 GB | 114 GB | 179 GB | 280 GB |
| 882k | 135 GB | 220 GB | 345 GB | 539 GB |

## Wall-clock time estimates

Approximate per-iteration times for the selection workflow (B ~256, 40-core CPU, 1x A100 GPU):

| Dataset size | Scoring | Reconstruction | Full iteration |
|:------------:|:-------:|:--------------:|:--------------:|
| ~88k | Minutes | Minutes | < 1 hour |
| ~200k | Tens of minutes | Minutes | 1-2 hours |
| ~450k | ~1 hour | Tens of minutes | 2-4 hours |
| ~730k | ~1.5 hours | Tens of minutes | 3-6 hours |
| ~880k | ~2 hours | Tens of minutes | 4-8 hours |

A full selection session typically converges in 3-6 iterations. Classification is a single pass with K scoring rounds.

## Practical recommendations

### Choosing `--mpi`

Set `--mpi` to the number of physical CPU cores (not hyperthreads). Values of 32-80 are typical for HPC nodes. Beyond ~80 cores, per-particle I/O overhead limits further scaling.

### Choosing `--gpu`

Specify device indices: `--gpu 0` or `--gpu 0 1`. Multi-GPU distributes subset reconstructions across devices. Adjust `--gpu-batch` if encountering out-of-memory errors.

### `--noExternalPrograms` vs RELION mode

| Feature | `--noExternalPrograms` | Without (RELION) |
|---------|----------------------|-----------------|
| Reconstruction | Internal GPU/CPU | RELION (MPI) |
| Local resolution | Internal CPU multiprocessing | `relion_postprocess` |
| Dependencies | Python + PyTorch only | RELION + MPI |
| GPU support | Yes | RELION's own |

`--noExternalPrograms` is self-contained and recommended when a GPU is available.

## Benchmarking

Two benchmark scripts are provided in `benchmarks/`:

```bash
# Selection workflow (single iteration)
./benchmarks/benchmark_single_iteration.sh \
    --particles particles.star \
    --map halfA.mrc --map2 halfB.mrc \
    --mask mask.mrc \
    --mpi 80 --gpu "0 1" --numRecs 15 \
    --output benchmark_results.csv

# Classification workflow
./benchmarks/benchmark_classification.sh \
    --particles particles.star \
    --maps class1.mrc class2.mrc class3.mrc \
    --mask mask.mrc \
    --mpi 85 --gpu "0 1" \
    --output benchmark_classification.csv
```

Both scripts output CSV files and print formatted summary tables.

---

[Back to documentation index](../index.md)
