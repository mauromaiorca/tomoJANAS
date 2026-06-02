# Install JANAS from source

This guide covers building JANAS from the repository. Use this if you want the latest development version, plan to modify the code, or simply prefer a clean local install (recommended over PyPI when you want full control over the source).

A single `pip install` from the source tree compiles and installs everything: the Python package, the C++ extension (`janas_core`), and the standalone C++ apps (`janas_app_starProcess`, `janas_app_meanMinMax`). No separate CMake step is needed.

## Prerequisites

- Python 3.8+
- C++ compiler (GCC on Linux, Clang on macOS)
- CMake 3.10+
- git

On Ubuntu/WSL2:

```bash
sudo apt update
sudo apt install -y python3 python3-pip python3-venv g++ cmake git
```

On macOS (Homebrew):

```bash
brew install cmake git
xcode-select --install   # provides the C++ compiler
```

## Full installation from scratch

This sequence creates a fresh virtual environment and installs JANAS in one go:

```bash
# 1. Create a working directory
mkdir -p ~/janas_install
cd ~/janas_install

# 2. Create and activate a new virtual environment
python3 -m venv ~/.janas_env
source ~/.janas_env/bin/activate

# 3. Upgrade pip
pip install --upgrade pip

# 4. Clone the repository
git clone https://github.com/mauromaiorca/janas.git
cd janas

# 5. (Optional, GPU only) Install PyTorch with the right CUDA wheel.
#    Pick the line that matches your CUDA driver (see `nvidia-smi`):
# python -m pip install torch --index-url https://download.pytorch.org/whl/cu128
# python -m pip install torch --index-url https://download.pytorch.org/whl/cu121
# python -m pip install torch --index-url https://download.pytorch.org/whl/cpu
# See the "GPU support (PyTorch)" section below for details.

# 6. Install JANAS (Python package + C++ extension + C++ apps)
pip install .

# 7. Verify
janas --version
janas_utils --version
janas_app_starProcess
janas_app_meanMinMax
```

In future sessions, re-activate the environment with:

```bash
source ~/.janas_env/bin/activate
```

To deactivate:

```bash
deactivate
```

## Updating an existing clone

```bash
source ~/.janas_env/bin/activate
cd ~/janas_install/janas
git pull origin main
pip install --upgrade .
janas --version
```

## GPU support (PyTorch)

GPU-accelerated 3D reconstruction (used when `--noExternalPrograms --gpu N` is set) requires **PyTorch**. PyTorch is intentionally **not** in JANAS's `requirements.txt` by default, because the right wheel depends on your CUDA driver.

### Recommended: install the torch wheel matching your CUDA driver

Check your CUDA driver with `nvidia-smi`, then pick the matching wheel from the [official PyTorch selector](https://pytorch.org/get-started/locally/). Examples:

```bash
# CUDA 12.8
python -m pip install torch --index-url https://download.pytorch.org/whl/cu128

# CUDA 12.1
python -m pip install torch --index-url https://download.pytorch.org/whl/cu121

# CUDA 11.8
python -m pip install torch --index-url https://download.pytorch.org/whl/cu118

# CPU-only (no GPU; reconstruction will work but be very slow)
python -m pip install torch --index-url https://download.pytorch.org/whl/cpu
```

### Alternative: install via the `gpu` extra (from the source checkout)

```bash
pip install '.[gpu]'
```

Uses the default torch wheel from PyPI. Fine on many systems but may not match unusual CUDA driver versions — in that case use the explicit `--index-url` command above instead.

If you only use RELION for reconstruction (omit `--noExternalPrograms`), PyTorch is **not** needed at all.

## Editable (development) install

For active development — changes to Python code take effect immediately without re-installing:

```bash
source ~/.janas_env/bin/activate
cd ~/janas_install/janas
pip install -e .
```

## Rebuild after C++ changes

If you modify files in `src/src_cpp/`, rebuild the C++ extension with:

```bash
python setup.py build_ext --inplace
```

## Uninstall

```bash
pip uninstall janas
```

To remove the virtual environment entirely:

```bash
deactivate
rm -rf ~/.janas_env
```

## Optional: shell activation shortcut

To avoid typing `source ~/.janas_env/bin/activate` every time:

```bash
cd ~/janas_install/janas
chmod +x ./src/janas/setup_environment_shell.sh
./src/janas/setup_environment_shell.sh install
source ~/.bashrc   # or source ~/.zshrc
```

Then from any directory, just run:

```bash
janas_activate_environment
```

## What gets installed

| Command | Description |
|---|---|
| `janas` | Main CLI: particle scoring, selection, classification |
| `janas_utils` | Utilities: masks, image crops, FSC, local resolution, half-map randomisation |
| `janas_optimizer` | Optimisation and overview analysis |
| `janas_session_manager` | Create and manage selection/classification sessions |
| `janas_reconstructor` | 3D reconstruction from scored particles |
| `janas_app_starProcess` | STAR file manipulation (C++) |
| `janas_app_meanMinMax` | Local resolution statistics (C++) |

## C++ apps usage reference

### janas_app_starProcess

Manipulates STAR files: inspect metadata, extract subsets, export formats, compare particles.

```bash
janas_app_starProcess --i input.star --o output.star [options]
```

| Option | Description |
|---|---|
| `--info` | Display particle count, labels, and subset distribution |
| `--infoEuler` | Display Euler angle statistics |
| `--hm h1.mrc h2.mrc [tag]` | Export particles to two half-map stacks |
| `--csv output.csv` | Export to CSV format |
| `--vem output.vem` | Export to VEM format |
| `--micrographs [depth]` | Extract unique micrograph list |
| `--checkForSimilarImages` | Find particles with similar coordinates |
| `--invertTagName tag1 tag2` | Swap two column values |
| `--backupImageNameTag [tag]` | Back up `_rlnImageName` to a custom tag |

Run `janas_app_starProcess --h` for the full list.

Examples:

```bash
# Inspect a STAR file
janas_app_starProcess --i particles.star --info

# Export to CSV
janas_app_starProcess --i particles.star --csv particles.csv

# Split into half-map stacks
janas_app_starProcess --i particles.star --o out.star --hm half1.mrc half2.mrc
```

### janas_app_meanMinMax

Computes statistics (mean, min, max) from a local resolution map within a masked region.

```bash
janas_app_meanMinMax locresMap.mrc mask.mrc
```

## Manual C++ compilation (optional)

If you only need the C++ apps without the Python package:

```bash
git clone https://github.com/mauromaiorca/janas.git
cd janas
mkdir build && cd build
cmake ..
make
make install   # installs to ~/.local/bin
```
