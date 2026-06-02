# Contributing

Contributions to tomoJANAS are welcome. This guide covers the basics of setting up a development environment and submitting changes.

## Development setup

```bash
git clone <repo-url>
cd tomoJANAS
python3 -m venv .tomojanas_env
source .tomojanas_env/bin/activate
pip install -e .
```

Editable mode (`-e`) means changes to Python code take effect immediately. If you modify C++ code in `src/src_cpp/`, rebuild with:

```bash
python setup.py build_ext --inplace
```

## Project structure

The codebase has two Python packages under `src/`:

- **`janas/`** — legacy SPA code inherited from the JANAS fork (C extension, reconstruction, scoring).
- **`tomojanas/`** — new ET import framework (MRC I/O, multi-block STAR, project layout, importers, geometry).

New ET-related code should go into `tomojanas/`. The `janas` package is kept for backward compatibility (the `tomojanas import_volume` legacy command still uses it).

## Running tests

Tests are self-contained scripts (no pytest runner required):

```bash
python tests/test_tomojanas_mrc.py
```

Each test file can also be discovered by pytest if you prefer:

```bash
pip install pytest
pytest tests/
```

## Submitting changes

1. Fork the repository and create a branch for your changes.
2. Make your changes and test them locally.
3. Open a pull request against `main` with a clear description of what changed and why.

## Reporting issues

Open an issue with:

- What you were trying to do
- The command you ran
- The error message or unexpected behaviour
- Your OS, Python version, and tomoJANAS version (`tomojanas --version`)
