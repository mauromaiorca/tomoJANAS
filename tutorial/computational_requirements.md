# JANAS v0.1.3 — Computational Requirements and Performance

## 1. Overview

JANAS provides two main computational workflows for cryo-EM particle processing:

1. **Iterative Particle Selection** — scores particles against reference half-maps using the Structural Cross-correlation Index (SCI), then iteratively selects optimal subsets that maximize local resolution within a target mask. Each iteration reconstructs multiple particle subsets and evaluates their local resolution.

2. **3D Class Reassignment** — scores particles against multiple reference maps (amplitude-equalized), assigns each particle to its best-matching class, and reconstructs per-class volumes. This is a single-pass workflow.

Both workflows share common computational primitives — particle scoring, 3D reconstruction, and local resolution estimation — but differ in how these are combined and repeated.

---

## 2. Per-Step Computational Analysis — Selection Workflow

Each selection iteration consists of the following steps, executed sequentially:

### 2.1 Particle Scoring (`janas scoreParticles`)

| Property | Value |
|----------|-------|
| **Algorithm** | Per-particle 3D→2D projection of reference volume at the particle's Euler angles, followed by masked SCI cross-correlation in real space |
| **Complexity** | O(N × B³ / P) |
| **Parallelization** | Python `multiprocessing.Process`, block partitioning of particles across P workers |
| **Hardware** | CPU only (no GPU acceleration) |
| **Key parameters** | `--mpi` (CPU processes), `--rank` (number of Euler views for normalization), `--sigma` (SCI Gaussian sigma) |
| **Memory** | O(B³) shared (reference map + mask) + O(B²) per worker (particle image buffer) |
| **Calls per iteration** | 1 |

**Scaling behaviour:** Near-linear speedup with CPU cores up to approximately 40–80 cores, where I/O overhead from reading MRC slices becomes significant. The per-particle work is dominated by the 3D→2D projection (O(B³)), making scoring the most CPU-intensive step per particle.

**Variables:** N = number of particles, B = box size in pixels, P = number of CPU processes.

### 2.2 3D Reconstruction (`janas_reconstructor`)

| Property | Value |
|----------|-------|
| **Algorithm** | Direct Fourier insertion: per-particle 2D FFT, CTF correction, and accumulation into a 3D Fourier grid at the particle's Euler orientation; followed by 3D inverse FFT |
| **Complexity** | O(N × B² log B) for Fourier insertion + O(B³ log B) for final 3D IFFT |
| **Parallelization** | **GPU path**: PyTorch CUDA with mini-batching (`--gpu-batch`, default 20). **CPU path**: `multiprocessing.Process` with block partitioning |
| **Hardware** | GPU (via `--gpu`) or CPU (via `--cpu`). GPU strongly recommended for datasets > 50k particles |
| **Key parameters** | `--gpu` (device indices), `--gpu-batch` (batch size, memory/speed trade-off), `--cpu` (workers for CPU path), `--subset` (half-map indices), `--subrec-only` (checkpoint particle counts for subset reconstructions) |
| **Memory (GPU)** | Complex64 accumulator: B³ × 8 bytes; Float32 weight grid: B³ × 4 bytes; Batch buffer: batch × B² × 4 bytes |
| **Memory (CPU)** | O(B³) per worker (accumulator) + O(B²) per particle |
| **Calls per iteration** | 1 (bootstrap, 2 half-maps) + ~15 (subset reconstructions via `--subrec-only`) |

**GPU vs CPU:** GPU reconstruction typically provides 5–20× speedup over multi-core CPU for medium-to-large datasets. The GPU path reads particles from disk, transfers to GPU in mini-batches, and maintains the accumulator on-device throughout. Multi-GPU is supported by specifying multiple device indices (e.g., `--gpu 0 1`).

**CTF correction modes:** `phaseflip`, `modulate`, or `wiener`. All add per-particle CTF kernel computation; `wiener` is marginally more expensive due to denominator division. The impact on total runtime is negligible compared to FFT and accumulation costs.

### 2.3 Local Resolution (`janas_utils locresBulk`)

| Property | Value |
|----------|-------|
| **Algorithm** | Sliding spherical window with raised-cosine edge; per-window local FSC computation across frequency shells; interpolation to find FSC threshold crossing |
| **Complexity** | O(S × (V/d)³ × R³ log R / P) |
| **Parallelization** | Python `multiprocessing.Process`, block partitioning of spherical-window centres across P workers |
| **Hardware** | CPU only |
| **Key parameters** | `--cpu` (processes), `--threshold` (FSC threshold, default 0.143), `--gamma` (auto-radius control), `--cycles` (auto-radius iterations) |
| **Memory** | O(V³) shared (two half-map volumes, read-only) + O(R³) per worker (windowed subvolumes) |
| **Calls per iteration** | ~15 (one per half-map pair from subset reconstructions) |

**Variables:** S = number of half-map pairs, V = volume size, d = sampling interval (auto-derived from global resolution), R = window radius in voxels, P = CPU processes.

**Scaling behaviour:** This is the most repeated step in the selection workflow (~15 calls per iteration). The number of regions grows cubically with (V/d), so datasets with high global resolution (small d) produce many more regions. Scales near-linearly with CPU cores.

**Bulk mode:** `locresBulk` processes multiple half-map pairs sequentially. The first pair computes the global FSC and derives auto-parameters (sampling, radius); subsequent pairs reuse these, making them slightly faster.

### 2.4 Auxiliary Steps (Negligible Runtime)

| Step | Complexity | Typical time |
|------|-----------|:------------:|
| `randomize_halves` | O(N) | < 5 s |
| `automaticParticleSubsets` | O(N log N) | < 5 s |
| `locresStats` | O(V³) | < 10 s |
| `getNumParticles` | O(1) | < 1 s |

These steps contribute < 1% of total iteration time and are not computational bottlenecks.

---

## 3. Per-Step Computational Analysis — Classification Workflow

### 3.1 Amplitude Equalization (`janas_utils equalize_images`)

| Property | Value |
|----------|-------|
| **Algorithm** | Compute 3D FFT of each reference map, extract amplitude spectra, average amplitudes across all K maps, replace each map's amplitudes with the average while preserving phases, inverse FFT |
| **Complexity** | O(K × B³ log B) |
| **Parallelization** | Serial (no parallelization) |
| **Hardware** | CPU only |
| **Calls** | 1 (preprocessing, before scoring) |

This step is typically fast (< 1 minute for K ≤ 8 maps at B ≤ 400) and is only run once.

### 3.2 Per-Map Particle Scoring

Identical to Section 2.1, but repeated **K times** (once per reference map). The total scoring time is:

> T_score_total = K × T_score(N)

This is the dominant bottleneck of the classification workflow. Scoring time scales linearly with both K (number of maps) and N (number of particles).

### 3.3 Per-Class Reconstruction

Identical to Section 2.2, called **K times** (once per output class). Since particles are split across classes, each reconstruction processes approximately N/K particles:

> T_recon_total ≈ K × T_recon(N/K) ≈ T_recon(N)

Total reconstruction time is roughly constant regardless of K (same total number of particles reconstructed).

---

## 4. Workflow-Level Time Estimates

### 4.1 Selection — Per Iteration

```
T_iter = T_randomize     (negligible)
       + T_bootstrap      = T_recon(N, 2 half-maps)
       + T_score(N)        = scoring all N particles
       + T_subsets          (negligible)
       + T_subrec           = sum over k of T_recon(N_k, 2 half-maps)
       + T_locresBulk       = S × T_locres_single(V)
       + T_stats            (negligible)
```

Where S ≈ 15 (number of subset reconstructions = `--numRecs`).

**Full session:** T_total = I × T_iter, where I = 3–8 iterations. Early termination occurs after `num_unimproving_iterations_for_early_termination` (default 3) consecutive non-improving iterations. Typical convergence requires 3–6 iterations.

### 4.2 Classification — Single Pass

```
T_classif = T_equalize(K)   (one-time, typically fast)
          + K × T_score(N)    (dominant bottleneck)
          + T_merge           (negligible)
          + K × T_recon(N/K)  (≈ T_recon(N) total)
```

---

## 5. Computational Complexity Summary

| Step | Algorithm | Complexity | Parallelization | Hardware | Selection (per iter.) | Classification (K maps) |
|------|-----------|-----------|-----------------|----------|:---------------------:|:----------------------:|
| Particle scoring | SCI via 3D→2D projection | O(NB³/P) | CPU multiprocessing | CPU | 1 call | K calls |
| 3D reconstruction | Direct Fourier insertion | O(NB² log B) | CUDA mini-batch or CPU multiprocessing | GPU or CPU | 1 + ~15 calls | K calls |
| Local resolution | Windowed local FSC | O(S(V/d)³R³ log R/P) | CPU multiprocessing | CPU | ~15 pairs | — |
| Amplitude equalization | Fourier amplitude averaging | O(KB³ log B) | Serial | CPU | — | 1 call |

N = particles, B = box size (px), P = CPU processes, S = half-map pairs, V = volume size (px), d = sampling interval (px), R = window radius (px), K = reference maps.

---

## 6. Hardware Requirements

| Component | Minimum | Recommended | Rationale |
|-----------|---------|-------------|-----------|
| CPU cores | 4 | 32–80 | Particle scoring and local resolution scale near-linearly with cores |
| RAM | 16 GB | 64–128 GB | Proportional to N × B²; each worker loads particle images |
| GPU | — | 1× NVIDIA A100 / V100 | Accelerates reconstruction 5–20× vs multi-core CPU |
| GPU VRAM | — | ≥ 16 GB | Required for B > 300 px with CTF correction; accumulator alone uses B³ × 12 bytes |
| Storage | 50 GB SSD | ≥ 500 GB (fast I/O) | Particle stacks dominate: N × B² × 4 bytes (e.g., 200k particles at B=256 ≈ 50 GB) |

### When is a GPU required?

- GPU is **not** required — JANAS can run entirely on CPU.
- GPU is **strongly recommended** when using `--noExternalPrograms` with datasets > 50k particles, as reconstruction is called ~16 times per selection iteration.
- Without `--noExternalPrograms`, reconstruction uses RELION (MPI-based) and local resolution uses `relion_postprocess`, which have their own parallelization strategies.

---

## 7. Memory Estimates

| Component | Formula | B=200 | B=256 | B=320 | B=400 |
|-----------|---------|:-----:|:-----:|:-----:|:-----:|
| Reference volume (float32) | B³ × 4 bytes | 31 MB | 64 MB | 125 MB | 244 MB |
| GPU accumulator + weight | B³ × 12 bytes | 92 MB | 192 MB | 375 MB | 732 MB |
| GPU particle batch (batch=20) | 20 × B² × 4 bytes | 3.1 MB | 5.0 MB | 7.8 MB | 12.2 MB |
| LocRes: two half-maps in memory | 2 × B³ × 4 bytes | 61 MB | 128 MB | 250 MB | 488 MB |

**Particle stack on disk:**

| N particles | B=200 | B=256 | B=320 | B=400 |
|:-----------:|:-----:|:-----:|:-----:|:-----:|
| 88k | 13 GB | 22 GB | 34 GB | 54 GB |
| 192k | 29 GB | 48 GB | 75 GB | 117 GB |
| 458k | 70 GB | 114 GB | 179 GB | 280 GB |
| 730k | 111 GB | 182 GB | 285 GB | 446 GB |
| 882k | 135 GB | 220 GB | 345 GB | 539 GB |

---

## 8. Representative Wall-Clock Times

The table below should be populated by running the provided benchmark scripts (`benchmarks/benchmark_single_iteration.sh` and `benchmarks/benchmark_classification.sh`) on the target hardware. Placeholder rows are shown for reference.

### 8.1 Selection Workflow — Per Iteration

| Step | Wall-clock (s) | % of iteration | Hardware |
|------|:--------------:|:--------------:|----------|
| randomize_halves | — | — | CPU |
| reconstruct (bootstrap) | — | — | GPU × 2 |
| scoreParticles | — | — | CPU × 80 |
| reconstruct (15 subsets) | — | — | GPU × 2 |
| locresBulk (15 pairs) | — | — | CPU × 80 |
| locresStats + getNumParticles | — | — | CPU |
| **Total (1 iteration)** | **—** | **100%** | |

### 8.2 Classification Workflow

| Step | Wall-clock (s) | % of total | Hardware |
|------|:--------------:|:----------:|----------|
| equalize_images | — | — | CPU |
| scoreParticles (K maps) | — | — | CPU × 80 |
| scores_to_csv + extract | — | — | CPU |
| reconstruct (K classes) | — | — | GPU × 2 |
| **Total** | **—** | **100%** | |

---

## 9. Practical Recommendations

### Choosing `--mpi` (CPU parallelism)

- Set `--mpi` to the number of physical CPU cores available (not hyperthreads).
- Values of 32–80 are typical for HPC nodes.
- Beyond ~80 cores, per-particle I/O overhead limits further scaling for scoring.
- The same value is used for both `scoreParticles` (via `--mpi`) and `locresBulk` (via `--cpu`).

### Choosing `--gpu` (reconstruction)

- Specify GPU device indices: `--gpu 0` (single GPU) or `--gpu 0 1` (two GPUs).
- Multi-GPU distributes subset reconstructions across devices.
- `--gpu-batch 20` (default) balances memory and throughput. Increase for GPUs with more VRAM; decrease if encountering OOM errors.

### `--noExternalPrograms` vs RELION mode

| Feature | `--noExternalPrograms` | Without (RELION) |
|---------|----------------------|-----------------|
| Reconstruction | Internal GPU/CPU | RELION (MPI) |
| Local resolution | Internal CPU multiprocessing | `relion_postprocess` |
| Dependencies | Python + PyTorch only | RELION + MPI installation |
| GPU support | Yes | RELION's own GPU support |

The `--noExternalPrograms` mode is self-contained and recommended when a GPU is available. It avoids external dependency management and typically provides faster reconstruction via GPU.

### Dataset-size guidelines

| Dataset size (particles) | Scoring (CPU × 80) | Reconstruction (1× GPU) | Full iteration (est.) |
|:------------------------:|:-------------------:|:-----------------------:|:---------------------:|
| ~88k | Minutes | Minutes | < 1 hour |
| ~200k | Tens of minutes | Minutes | 1–2 hours |
| ~450k | ~1 hour | Tens of minutes | 2–4 hours |
| ~730k | ~1.5 hours | Tens of minutes | 3–6 hours |
| ~880k | ~2 hours | Tens of minutes | 4–8 hours |

*Estimates are for B ≈ 256 px on a 40-core CPU node with 1× A100 GPU. Actual times depend on hardware, I/O speed, and dataset-specific parameters. Run `benchmarks/benchmark_single_iteration.sh` for precise values.*

### CryoSPARC Local NU-refinement (optional)

When enabled via `--multipleLocalRefineCS`, JANAS can run CryoSPARC's local non-uniform refinement between iterations. This is an external GPU job dispatched to CryoSPARC's scheduler and typically adds 15–60 minutes per iteration depending on particle count and GPU hardware. It is skipped automatically when the previous iteration showed no improvement.

---

## 10. Benchmarking

Two benchmark scripts are provided in `benchmarks/`:

### `benchmark_single_iteration.sh`

Runs one complete selection iteration with per-step timing:

```bash
./benchmarks/benchmark_single_iteration.sh \
  --particles class_1.star \
  --map class_1_recH1.mrc \
  --map2 class_1_recH2.mrc \
  --mask masks/mask_partB.mrc \
  --assessMask masks/maskFlexibleRegion.mrc \
  --apix 1.350 --sigma 1.0 \
  --mpi 80 --gpu "0 1" --numRecs 15 \
  --output benchmark_results.csv
```

### `benchmark_classification.sh`

Runs the full classification workflow with per-step timing:

```bash
./benchmarks/benchmark_classification.sh \
  --particles HTT_46Q_stack.star \
  --maps class1.mrc class2.mrc class3.mrc class4.mrc \
  --mask reference_mask/maskOpenDilateClosing.mrc \
  --apix 1.350 --sigma 1.0 \
  --mpi 85 --gpu "0 1" \
  --output benchmark_classification_results.csv
```

Both scripts output CSV files and print formatted summary tables to stdout.
