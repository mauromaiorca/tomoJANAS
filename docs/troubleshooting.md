[Repository home](../README.md) · [Docs home](index.md) · [Installation](installation.md) · [Quick start](quick-start.md) · [CLI reference](reference/cli.md)

---

# Troubleshooting

## Installation issues

| Problem | Solution |
|---------|----------|
| `pip: command not found` | Use `pip3` or `python3 -m pip` |
| `externally-managed-environment` on macOS | Install in a venv or use pipx. See [Installation](installation.md). |
| CMake or compiler errors during install | Ensure `cmake` and `g++` (Linux) or `clang++` (macOS) are installed |
| Commands not found after pipx install | Run `pipx ensurepath` and restart your terminal |
| `janas_app_starProcess` not found | The C++ apps were not compiled. Reinstall with `pip install janas` or build from source. |

## Runtime issues

### RELION commands not found

If you see errors about `relion_reconstruct` or `relion_postprocess` not being found, either:

- Install RELION and ensure it is on your `PATH`, or
- Use `--noExternalPrograms` to run JANAS without external dependencies (requires a GPU for reasonable performance on large datasets)

### `ModuleNotFoundError: No module named 'cmake'` (during `pip install .`)

Some conda environments (e.g. colabfold's bundled conda) ship a Python wrapper at `<conda>/bin/cmake` that depends on the `cmake` Python package, but ship without that package installed. The build then fails inside our CMake step.

Fix:

```bash
pip install 'cmake>=3.10'
```

Then retry `pip install .`. From v1.0.7 onwards, `cmake>=3.10` is declared as a build dependency in `pyproject.toml`, so this happens automatically under `pip`'s default build isolation.

### `ModuleNotFoundError: No module named 'torch'` / "PyTorch is required for GPU reconstruction"

PyTorch is not in JANAS's `requirements.txt` by default because the correct wheel depends on your CUDA driver version. Two options:

1. **Recommended** — install the torch wheel matching your CUDA driver (check with `nvidia-smi`):

   ```bash
   python -m pip install torch --index-url https://download.pytorch.org/whl/cu128
   ```

2. **Generic install** — uses the default torch wheel from PyPI:

   ```bash
   pip install 'janas[gpu]'
   ```

See the [Installation guide](installation.md#gpu-support-pytorch) for full details and other CUDA versions.

### Out of GPU memory

Reduce `--gpu-batch` (default 20). For very large box sizes (B > 350), values of 5-10 may be needed.

### Scoring is slow

- Increase `--mpi` to use more CPU cores (up to ~80).
- Check that the particle stack is on fast storage (SSD preferred).
- Scoring is CPU-only; more cores directly reduces wall-clock time.

### Stack file not found during selection

The run script expects the particle stack to be accessible from the session directory. Create a symbolic link:

```bash
ln -s /absolute/path/to/particles.mrcs session_directory/
```

### Particle count mismatch with `update_from_csparc`

The `.cs` file and STAR file must contain the same number of particles in the same order. Check for trailing blank lines in the STAR file's `data_particles` section and remove them.

## Getting help

Open an issue on [GitHub](https://github.com/mauromaiorca/janas/issues) or contact mauro.maiorca@cssb-hamburg.de.

---

[Back to documentation index](index.md)
