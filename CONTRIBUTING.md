# Contributing

Contributions to JANAS are welcome. This guide covers the basics of setting up a development environment and submitting changes.

## Development setup

```bash
git clone https://github.com/mauromaiorca/janas.git
cd janas
python3 -m venv .janas_env
source .janas_env/bin/activate
pip install -e .
```

Editable mode (`-e`) means changes to Python code take effect immediately. If you modify C++ code in `src/src_cpp/`, rebuild with:

```bash
python setup.py build_ext --inplace
```

## Submitting changes

1. Fork the repository and create a branch for your changes.
2. Make your changes and test them locally.
3. Open a pull request against `main` with a clear description of what changed and why.

## Reporting issues

Open an issue on [GitHub](https://github.com/mauromaiorca/janas/issues) with:

- What you were trying to do
- The command you ran
- The error message or unexpected behaviour
- Your OS, Python version, and JANAS version (`janas --version`)
