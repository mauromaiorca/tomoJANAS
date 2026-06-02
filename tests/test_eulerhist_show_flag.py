"""
Regression tests for ``janas eulerHist --show`` parsing.

The previous declaration used ``type=bool`` in argparse, which silently
turned ``--show False`` back into ``True`` because ``bool("False")`` is
truthy. The fix is a small ``_str2bool`` helper used as ``type=`` of
the argument. These tests confirm both the helper and the argparse
behaviour end-to-end.

Run with::

    python tests/test_eulerhist_show_flag.py
"""

from __future__ import annotations

import os
import sys
import traceback
import types
from typing import Callable, List, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.normpath(os.path.join(HERE, "..", "src"))
if os.path.isdir(SRC) and SRC not in sys.path:
    sys.path.insert(0, SRC)

# Stub the C++ extension so we can import the CLI module without building.
if "janas.janas_core" not in sys.modules:
    sys.modules["janas.janas_core"] = types.ModuleType("janas.janas_core")

from janas import janas_cmd_caller as C  # noqa: E402


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def test_str2bool_true_variants() -> None:
    for s in ("true", "True", "TRUE", "t", "T", "yes", "y", "1", "on", True):
        assert C._str2bool(s) is True, s


def test_str2bool_false_variants() -> None:
    for s in ("false", "False", "FALSE", "f", "no", "n", "0", "off", False):
        assert C._str2bool(s) is False, s


def test_str2bool_invalid_raises() -> None:
    import argparse
    for s in ("maybe", "", "2", "bla"):
        try:
            C._str2bool(s)
        except argparse.ArgumentTypeError:
            continue
        raise AssertionError(f"expected ArgumentTypeError for {s!r}")


# ---------------------------------------------------------------------------
# argparse end-to-end
# ---------------------------------------------------------------------------


def _parse_show(value: str) -> bool:
    """Run the real top-level parser with a synthetic eulerHist call and
    return the resolved ``args.show`` value."""
    args = C.janas_parser.parse_args(
        ["eulerHist", "--i", "/tmp/_nonexistent.star",
         "--show", value]
    )
    return args.show


def test_argparse_show_false_actually_parses_false() -> None:
    assert _parse_show("False") is False
    assert _parse_show("false") is False
    assert _parse_show("0") is False
    assert _parse_show("no") is False
    assert _parse_show("off") is False


def test_argparse_show_true_parses_true() -> None:
    assert _parse_show("True") is True
    assert _parse_show("true") is True
    assert _parse_show("1") is True
    assert _parse_show("yes") is True
    assert _parse_show("on") is True


def test_argparse_show_default_is_true() -> None:
    args = C.janas_parser.parse_args(
        ["eulerHist", "--i", "/tmp/_nonexistent.star"]
    )
    assert args.show is True


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------


TESTS: List[Tuple[str, Callable[[], None]]] = [
    ("_str2bool: 'true'/'yes'/'1'/'on' all map to True",
        test_str2bool_true_variants),
    ("_str2bool: 'false'/'no'/'0'/'off' all map to False",
        test_str2bool_false_variants),
    ("_str2bool: invalid strings raise ArgumentTypeError",
        test_str2bool_invalid_raises),
    ("argparse: --show False / false / 0 / no / off all parse to False",
        test_argparse_show_false_actually_parses_false),
    ("argparse: --show True / true / 1 / yes / on all parse to True",
        test_argparse_show_true_parses_true),
    ("argparse: --show default is True",
        test_argparse_show_default_is_true),
]


def main() -> int:
    failed: List[str] = []
    for name, fn in TESTS:
        try:
            fn()
            print(f"[ OK ] {name}")
        except Exception as exc:  # noqa: BLE001
            failed.append(name)
            print(f"[FAIL] {name}: {exc}")
            traceback.print_exc(limit=3)
    print()
    if failed:
        print(f"{len(failed)}/{len(TESTS)} test(s) failed: {failed}")
        return 1
    print(f"All {len(TESTS)} tests passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
