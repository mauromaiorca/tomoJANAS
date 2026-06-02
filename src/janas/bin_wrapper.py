"""Thin wrappers that exec the compiled C++ binaries bundled in janas/bin/."""

import os
import sys
from pathlib import Path


def _bin_dir():
    return Path(__file__).resolve().parent / "bin"


def _exec_binary(name):
    binary = _bin_dir() / name
    if not binary.exists():
        print(f"Error: {name} not found at {binary}", file=sys.stderr)
        print("You may need to reinstall janas: pip install --force-reinstall janas", file=sys.stderr)
        sys.exit(1)
    os.execv(str(binary), [str(binary)] + sys.argv[1:])


def meanMinMax_main():
    _exec_binary("janas_app_meanMinMax")


def starProcess_main():
    _exec_binary("janas_app_starProcess")
