"""
Unit tests for :func:`janas.utils._clean_csparc_blob_path`.

Run from the repository root with::

    python tests/test_clean_csparc_blob_path.py

Exits with code 0 if all tests pass, 1 otherwise. Self-contained: stubs
``janas.janas_core`` so it does not require the C++ extension to be built.
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

if "janas.janas_core" not in sys.modules:
    sys.modules["janas.janas_core"] = types.ModuleType("janas.janas_core")

from janas import utils  # noqa: E402

CLEAN = utils._clean_csparc_blob_path


# Reusable inputs
ORIGINAL = (
    "J41/extract/008878123945933052272_HTT_46Q_0000_Aug06_00.39.32_aligned_DW_particles.mrc"
)
FNAME = "008878123945933052272_HTT_46Q_0000_Aug06_00.39.32_aligned_DW_particles.mrc"
NO_PREFIX = "J41/extract/HTT_46Q_0000_Aug06_00.39.32_aligned_DW_particles.mrc"
NO_SUFFIX = "J41/extract/008878123945933052272_HTT_46Q_0000_Aug06_00.39.32_aligned_DW.mrc"


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_default_passthrough() -> None:
    assert CLEAN(ORIGINAL) == ORIGINAL


def test_clean_path_keeps_only_filename() -> None:
    assert CLEAN(ORIGINAL, clean_path=True) == FNAME


def test_clean_prefix_strips_only_numeric_prefix_of_filename() -> None:
    assert CLEAN(ORIGINAL, clean_prefix=True) == NO_PREFIX


def test_clean_suffix_strips_terminal_particles() -> None:
    assert CLEAN(ORIGINAL, clean_suffix=True) == NO_SUFFIX


def test_fix_path_replaces_directory() -> None:
    assert (
        CLEAN(ORIGINAL, fix_path="data/upload")
        == "data/upload/" + FNAME
    )


def test_fix_path_with_clean_prefix_and_clean_suffix() -> None:
    expected = "data/upload/HTT_46Q_0000_Aug06_00.39.32_aligned_DW.mrc"
    got = CLEAN(
        ORIGINAL,
        fix_path="data/upload",
        clean_prefix=True,
        clean_suffix=True,
    )
    assert got == expected, (got, expected)


def test_fix_path_trailing_slash_no_double_slash() -> None:
    out = CLEAN(ORIGINAL, fix_path="data/upload/")
    assert "//" not in out, out
    assert out == "data/upload/" + FNAME


def test_fix_path_takes_precedence_over_clean_path() -> None:
    # When both are passed, fix_path wins (directory is replaced, not dropped).
    assert (
        CLEAN(ORIGINAL, fix_path="data/upload", clean_path=True)
        == "data/upload/" + FNAME
    )


def test_clean_prefix_does_not_strip_non_numeric_prefix() -> None:
    # Non-numeric prefix must be preserved.
    p = "J41/extract/HTT_46Q_0000_aligned_DW_particles.mrc"
    assert CLEAN(p, clean_prefix=True) == p


def test_clean_prefix_does_not_touch_directory_digits() -> None:
    # Numeric leading digits in the directory part must not be affected;
    # only the filename component is considered.
    p = "008878_J41/extract/HTT_46Q_aligned_DW_particles.mrc"
    assert (
        CLEAN(p, clean_prefix=True)
        == "008878_J41/extract/HTT_46Q_aligned_DW_particles.mrc"
    )


def test_clean_suffix_only_strips_terminal_particles_before_extension() -> None:
    # "_particles" appearing earlier in the name must stay.
    p = "J41/extract/HTT_particles_aligned.mrc"
    assert CLEAN(p, clean_suffix=True) == p
    # And a name that does NOT end in "_particles" is untouched.
    p2 = "J41/extract/HTT_aligned.mrc"
    assert CLEAN(p2, clean_suffix=True) == p2


def test_clean_suffix_preserves_extension() -> None:
    p = "J41/extract/X_particles.mrcs"
    assert CLEAN(p, clean_suffix=True) == "J41/extract/X.mrcs"


def test_leading_gt_stripped() -> None:
    # raw_paths.lstrip(">") behaviour mirrored here.
    assert CLEAN(">" + ORIGINAL) == ORIGINAL


def test_backslashes_converted() -> None:
    win = ORIGINAL.replace("/", "\\")
    assert CLEAN(win) == ORIGINAL


def test_filename_only_input() -> None:
    # No '/' in input: behaves correctly under each flag.
    assert CLEAN(FNAME) == FNAME
    assert CLEAN(FNAME, clean_path=True) == FNAME
    assert CLEAN(FNAME, clean_prefix=True) == FNAME.replace(
        "008878123945933052272_", ""
    )
    assert CLEAN(FNAME, fix_path="data/upload") == "data/upload/" + FNAME


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

TESTS: List[Tuple[str, Callable[[], None]]] = [
    ("default passthrough", test_default_passthrough),
    ("clean_path keeps only filename", test_clean_path_keeps_only_filename),
    ("clean_prefix strips only numeric prefix of filename",
        test_clean_prefix_strips_only_numeric_prefix_of_filename),
    ("clean_suffix strips terminal _particles", test_clean_suffix_strips_terminal_particles),
    ("fix_path replaces directory", test_fix_path_replaces_directory),
    ("fix_path + clean_prefix + clean_suffix",
        test_fix_path_with_clean_prefix_and_clean_suffix),
    ("fix_path trailing slash no double slash",
        test_fix_path_trailing_slash_no_double_slash),
    ("fix_path takes precedence over clean_path",
        test_fix_path_takes_precedence_over_clean_path),
    ("clean_prefix does NOT strip non-numeric prefix",
        test_clean_prefix_does_not_strip_non_numeric_prefix),
    ("clean_prefix does not touch directory digits",
        test_clean_prefix_does_not_touch_directory_digits),
    ("clean_suffix only strips terminal _particles before extension",
        test_clean_suffix_only_strips_terminal_particles_before_extension),
    ("clean_suffix preserves extension (.mrcs)", test_clean_suffix_preserves_extension),
    ("leading '>' stripped", test_leading_gt_stripped),
    ("backslashes converted to forward slashes", test_backslashes_converted),
    ("filename-only input handled", test_filename_only_input),
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
