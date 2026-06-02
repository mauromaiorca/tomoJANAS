"""
Unit tests for the new ``path_mode`` resolution helpers in
:mod:`janas.utils`:

- :func:`_candidate_stack_paths`
- :func:`_resolve_existing_stack_path`

Run from the repository root with::

    python tests/test_create_stack_path_mode.py

Stubs ``janas.janas_core`` so it does not require the C++ extension.
"""

from __future__ import annotations

import os
import sys
import tempfile
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

CAND = utils._candidate_stack_paths
RESOLVE = utils._resolve_existing_stack_path


# ---------------------------------------------------------------------------
# _candidate_stack_paths
# ---------------------------------------------------------------------------


def test_absolute_path_taken_as_is_in_every_mode() -> None:
    abs_path = os.path.abspath("/data/foo/bar.mrcs")
    for mode in ("auto", "root", "as_is", "star_dir"):
        assert CAND(abs_path, "/proj", "/star", mode) == [abs_path], mode


def test_root_mode_joins_to_project_root() -> None:
    out = CAND("J41/extract/x.mrc", project_root="/p", star_dir="/s", path_mode="root")
    assert out == [os.path.normpath("/p/J41/extract/x.mrc")]


def test_root_mode_without_root_falls_back_to_as_written() -> None:
    out = CAND("J41/extract/x.mrc", project_root=None, star_dir="/s", path_mode="root")
    assert out == [os.path.normpath("J41/extract/x.mrc")]


def test_as_is_mode_does_not_prepend_root_or_star_dir() -> None:
    out = CAND("J41/x.mrc", project_root="/p", star_dir="/s", path_mode="as_is")
    assert out == [os.path.normpath("J41/x.mrc")]


def test_star_dir_mode_joins_to_star_dir() -> None:
    out = CAND("rel/x.mrc", project_root="/p", star_dir="/s", path_mode="star_dir")
    assert out == [os.path.normpath("/s/rel/x.mrc")]


def test_star_dir_mode_without_star_dir_falls_back() -> None:
    out = CAND("rel/x.mrc", project_root="/p", star_dir=None, path_mode="star_dir")
    assert out == [os.path.normpath("rel/x.mrc")]


def test_auto_mode_priority_root_then_as_written_then_star_dir() -> None:
    out = CAND("rel/x.mrc", project_root="/p", star_dir="/s", path_mode="auto")
    assert out == [
        os.path.normpath("/p/rel/x.mrc"),
        os.path.normpath("rel/x.mrc"),
        os.path.normpath("/s/rel/x.mrc"),
    ]


def test_auto_mode_no_root_no_star_dir() -> None:
    out = CAND("rel/x.mrc", project_root=None, star_dir=None, path_mode="auto")
    assert out == [os.path.normpath("rel/x.mrc")]


def test_auto_mode_deduplicates_when_root_or_star_dir_match_cwd_pathing() -> None:
    # If project_root == star_dir, we should not emit the same candidate twice.
    out = CAND("rel/x.mrc", project_root="/same", star_dir="/same", path_mode="auto")
    abs1 = os.path.normpath("/same/rel/x.mrc")
    assert out == [abs1, os.path.normpath("rel/x.mrc")] or out == [
        abs1,
        os.path.normpath("rel/x.mrc"),
    ]
    assert len(out) == len(set(out))


def test_unknown_path_mode_raises() -> None:
    try:
        CAND("rel/x.mrc", "/p", "/s", "bogus")
    except ValueError as e:
        assert "Unknown path_mode" in str(e), str(e)
        return
    raise AssertionError("unknown path_mode should raise ValueError")


# ---------------------------------------------------------------------------
# _resolve_existing_stack_path
# ---------------------------------------------------------------------------


def test_resolve_decodes_index_and_path() -> None:
    idx, p = RESOLVE("000005@rel/x.mrc", project_root=None, star_dir=None, path_mode="as_is")
    assert idx == 4  # idx0 = image_no - 1
    assert p == os.path.normpath("rel/x.mrc")


def test_resolve_invalid_image_name_raises() -> None:
    try:
        RESOLVE("no_at_sign", None, None, "auto")
    except ValueError as e:
        assert "no '@'" in str(e), str(e)
        return
    raise AssertionError("missing @ should raise")


def test_resolve_invalid_index_raises() -> None:
    try:
        RESOLVE("abc@x.mrc", None, None, "auto")
    except ValueError as e:
        assert "Invalid image index" in str(e), str(e)
        return
    raise AssertionError("non-integer index should raise")


def _touch(p: str) -> None:
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "wb") as f:
        f.write(b"")


def test_resolve_auto_picks_first_existing() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        # Layout:
        #   tmp/proj/J41/x.mrc     <- should be picked (root candidate)
        #   tmp/star/J41/x.mrc     <- also exists (star_dir candidate)
        proj = os.path.join(tmp, "proj")
        star = os.path.join(tmp, "star")
        in_proj = os.path.join(proj, "J41", "x.mrc")
        in_star = os.path.join(star, "J41", "x.mrc")
        _touch(in_proj)
        _touch(in_star)

        idx, resolved = RESOLVE(
            "000001@J41/x.mrc",
            project_root=proj,
            star_dir=star,
            path_mode="auto",
        )
        assert idx == 0
        assert resolved == os.path.normpath(in_proj), resolved


def test_resolve_auto_falls_back_to_star_dir_when_root_missing() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        proj = os.path.join(tmp, "proj")           # exists, but the file does NOT
        star = os.path.join(tmp, "star")
        in_star = os.path.join(star, "J41", "x.mrc")
        os.makedirs(proj, exist_ok=True)
        _touch(in_star)

        _, resolved = RESOLVE(
            "000001@J41/x.mrc",
            project_root=proj,
            star_dir=star,
            path_mode="auto",
        )
        assert resolved == os.path.normpath(in_star), resolved


def test_resolve_auto_returns_first_candidate_when_none_exist() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        _, resolved = RESOLVE(
            "000001@J41/x.mrc",
            project_root=os.path.join(tmp, "proj"),
            star_dir=os.path.join(tmp, "star"),
            path_mode="auto",
        )
        # First candidate is the root-joined one.
        assert resolved == os.path.normpath(
            os.path.join(tmp, "proj", "J41", "x.mrc")
        ), resolved


def test_resolve_root_mode_returns_root_joined_path_without_filesystem_check() -> None:
    _, resolved = RESOLVE(
        "000001@J41/x.mrc",
        project_root="/no/such/proj",
        star_dir="/no/such/star",
        path_mode="root",
    )
    assert resolved == os.path.normpath("/no/such/proj/J41/x.mrc")


def test_resolve_as_is_mode_returns_as_written() -> None:
    _, resolved = RESOLVE(
        "000001@J41/x.mrc",
        project_root="/proj",
        star_dir="/star",
        path_mode="as_is",
    )
    assert resolved == os.path.normpath("J41/x.mrc")


def test_resolve_star_dir_mode_returns_star_dir_joined() -> None:
    _, resolved = RESOLVE(
        "000001@J41/x.mrc",
        project_root="/proj",
        star_dir="/star",
        path_mode="star_dir",
    )
    assert resolved == os.path.normpath("/star/J41/x.mrc")


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

TESTS: List[Tuple[str, Callable[[], None]]] = [
    ("absolute path taken as is in every mode", test_absolute_path_taken_as_is_in_every_mode),
    ("root mode joins to project root", test_root_mode_joins_to_project_root),
    ("root mode without root falls back to as-written", test_root_mode_without_root_falls_back_to_as_written),
    ("as_is mode does not prepend root or star_dir", test_as_is_mode_does_not_prepend_root_or_star_dir),
    ("star_dir mode joins to star_dir", test_star_dir_mode_joins_to_star_dir),
    ("star_dir mode without star_dir falls back", test_star_dir_mode_without_star_dir_falls_back),
    ("auto priority: root, then as-written, then star_dir", test_auto_mode_priority_root_then_as_written_then_star_dir),
    ("auto mode no root no star_dir", test_auto_mode_no_root_no_star_dir),
    ("auto dedup when root == star_dir", test_auto_mode_deduplicates_when_root_or_star_dir_match_cwd_pathing),
    ("unknown path_mode raises ValueError", test_unknown_path_mode_raises),
    ("resolve decodes index and path", test_resolve_decodes_index_and_path),
    ("resolve invalid image name raises", test_resolve_invalid_image_name_raises),
    ("resolve invalid index raises", test_resolve_invalid_index_raises),
    ("auto picks first existing candidate", test_resolve_auto_picks_first_existing),
    ("auto falls back to star_dir when root file missing", test_resolve_auto_falls_back_to_star_dir_when_root_missing),
    ("auto returns first candidate when none exist", test_resolve_auto_returns_first_candidate_when_none_exist),
    ("root mode does not probe filesystem", test_resolve_root_mode_returns_root_joined_path_without_filesystem_check),
    ("as_is mode returns as written", test_resolve_as_is_mode_returns_as_written),
    ("star_dir mode returns star_dir joined", test_resolve_star_dir_mode_returns_star_dir_joined),
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
