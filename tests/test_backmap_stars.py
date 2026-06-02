"""
Self-contained tests for :func:`janas.utils.backmap_stars`.

Run from the repository root with::

    python tests/test_backmap_stars.py

Exits with code 0 if all tests pass, 1 otherwise. Avoids any test framework so
it can be executed in any environment where ``janas`` is importable.
"""

from __future__ import annotations

import os
import sys
import tempfile
import traceback
import types
from textwrap import dedent
from typing import Callable, List, Tuple


# Make sure we import the in-tree janas package, not a possibly older one
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.normpath(os.path.join(HERE, "..", "src"))
if os.path.isdir(SRC) and SRC not in sys.path:
    sys.path.insert(0, SRC)

# Stub the compiled C++ extension so we can import janas.utils without first
# building janas_core. backmap_stars itself does not touch any janas_core
# symbol — it only depends on starHandler. This makes the test runnable on a
# clean clone without `pip install .`.
if "janas.janas_core" not in sys.modules:
    sys.modules["janas.janas_core"] = types.ModuleType("janas.janas_core")

from janas import utils  # noqa: E402
from janas import starHandler  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

OPTICS_BLOCK = dedent(
    """\
    data_optics

    loop_
    _rlnOpticsGroup #1
    _rlnOpticsGroupName #2
    _rlnAmplitudeContrast #3
    _rlnSphericalAberration #4
    _rlnVoltage #5
    _rlnImagePixelSize #6
    _rlnImageSize #7
    _rlnImageDimensionality #8
    1 opticsGroup1 0.100000 2.700000 300.000000 0.840000 256 2
    """
)


def _write_mapping_star(path: str, rows: List[Tuple[str, str]]) -> None:
    """Write a minimal RELION 3.1 STAR with _rlnImageName + _janas_source_rlnImageName."""
    lines: List[str] = [OPTICS_BLOCK, "data_particles", "", "loop_"]
    lines += ["_rlnImageName #1", "_janas_source_rlnImageName #2"]
    for img, src in rows:
        lines.append(f"{img} {src}")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _write_processed_star(path: str, rows: List[Tuple[str, str, str]]) -> None:
    """
    Write a minimal RELION 3.1 STAR with three columns:
    _rlnImageName, _rlnAngleRot, _rlnClassNumber.
    """
    lines: List[str] = [OPTICS_BLOCK, "data_particles", "", "loop_"]
    lines += [
        "_rlnImageName #1",
        "_rlnAngleRot #2",
        "_rlnClassNumber #3",
    ]
    for img, ang, cls in rows:
        lines.append(f"{img} {ang} {cls}")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _read_particles(path: str):
    """Read the 'particles' block back into a DataFrame for verification."""
    return utils._read_particle_section_df(path, "particles")  # noqa: SLF001


# ---------------------------------------------------------------------------
# Test cases
# ---------------------------------------------------------------------------


def test_basic_mapping_ignores_row_order(tmp: str) -> None:
    mapping = os.path.join(tmp, "mapping.star")
    processed = os.path.join(tmp, "processed.star")
    out = os.path.join(tmp, "out.star")

    _write_mapping_star(
        mapping,
        [
            ("000001@stack.mrcs", "009619@orig.mrc"),
            ("000002@stack.mrcs", "009620@orig.mrc"),
            ("000003@stack.mrcs", "009621@orig.mrc"),
        ],
    )
    _write_processed_star(
        processed,
        [
            ("000003@stack.mrcs", "10.0", "2"),
            ("000001@stack.mrcs", "20.0", "1"),
        ],
    )

    report = utils.backmap_stars(
        processed_star=processed,
        mapping_star=mapping,
        output_star=out,
        stack_reference_tag=None,  # keep output minimal for this test
    )

    df = _read_particles(out)
    assert list(df["_rlnImageName"]) == [
        "009621@orig.mrc",
        "009619@orig.mrc",
    ], df["_rlnImageName"].tolist()
    assert report["n_processed"] == 2
    assert report["n_mapped"] == 2
    assert report["n_missing"] == 0
    assert report["n_mapping_rows"] == 3


def test_preserves_metadata_columns(tmp: str) -> None:
    mapping = os.path.join(tmp, "mapping.star")
    processed = os.path.join(tmp, "processed.star")
    out = os.path.join(tmp, "out.star")

    _write_mapping_star(
        mapping,
        [
            ("000001@stack.mrcs", "001@orig.mrc"),
            ("000002@stack.mrcs", "002@orig.mrc"),
        ],
    )
    _write_processed_star(
        processed,
        [
            ("000001@stack.mrcs", "12.5", "1"),
            ("000002@stack.mrcs", "78.25", "3"),
        ],
    )

    utils.backmap_stars(processed, mapping, out, stack_reference_tag=None)

    df = _read_particles(out)
    assert list(df["_rlnAngleRot"]) == ["12.5", "78.25"]
    assert list(df["_rlnClassNumber"]) == ["1", "3"]


def test_adds_stack_reference_audit_column(tmp: str) -> None:
    mapping = os.path.join(tmp, "mapping.star")
    processed = os.path.join(tmp, "processed.star")
    out = os.path.join(tmp, "out.star")

    _write_mapping_star(
        mapping,
        [
            ("000001@stack.mrcs", "001@orig.mrc"),
            ("000002@stack.mrcs", "002@orig.mrc"),
        ],
    )
    _write_processed_star(
        processed,
        [
            ("000002@stack.mrcs", "0", "1"),
            ("000001@stack.mrcs", "0", "1"),
        ],
    )

    utils.backmap_stars(processed, mapping, out)  # default audit column

    df = _read_particles(out)
    assert "_janas_stack_rlnImageName" in df.columns
    assert list(df["_janas_stack_rlnImageName"]) == [
        "000002@stack.mrcs",
        "000001@stack.mrcs",
    ]
    assert list(df["_rlnImageName"]) == ["002@orig.mrc", "001@orig.mrc"]


def test_strict_missing_raises(tmp: str) -> None:
    mapping = os.path.join(tmp, "mapping.star")
    processed = os.path.join(tmp, "processed.star")
    out = os.path.join(tmp, "out.star")

    _write_mapping_star(
        mapping,
        [("000001@stack.mrcs", "001@orig.mrc")],
    )
    _write_processed_star(
        processed,
        [
            ("000001@stack.mrcs", "0", "1"),
            ("999999@stack.mrcs", "0", "1"),  # absent
        ],
    )
    try:
        utils.backmap_stars(processed, mapping, out, strict=True)
    except ValueError as e:
        assert "999999@stack.mrcs" in str(e), str(e)
        assert not os.path.exists(out), "output must not be written when strict and missing"
        return
    raise AssertionError("strict=True must raise on missing keys")


def test_non_strict_leaves_missing_unchanged(tmp: str) -> None:
    mapping = os.path.join(tmp, "mapping.star")
    processed = os.path.join(tmp, "processed.star")
    out = os.path.join(tmp, "out.star")

    _write_mapping_star(
        mapping,
        [("000001@stack.mrcs", "001@orig.mrc")],
    )
    _write_processed_star(
        processed,
        [
            ("000001@stack.mrcs", "0", "1"),
            ("999999@stack.mrcs", "0", "1"),
        ],
    )

    report = utils.backmap_stars(processed, mapping, out, strict=False)

    df = _read_particles(out)
    # mapped row -> original, missing row -> kept unchanged
    assert list(df["_rlnImageName"]) == ["001@orig.mrc", "999999@stack.mrcs"]
    assert report["n_missing"] == 1
    assert report["missing_examples"] == ["999999@stack.mrcs"]
    # audit column carries the previous (stack) values in both cases
    assert list(df["_janas_stack_rlnImageName"]) == [
        "000001@stack.mrcs",
        "999999@stack.mrcs",
    ]


def test_missing_source_tag_in_mapping_raises(tmp: str) -> None:
    mapping = os.path.join(tmp, "mapping.star")
    processed = os.path.join(tmp, "processed.star")
    out = os.path.join(tmp, "out.star")

    # mapping STAR with only _rlnImageName -> no provenance column
    with open(mapping, "w") as f:
        f.write(OPTICS_BLOCK)
        f.write(
            "\n".join(
                [
                    "data_particles",
                    "",
                    "loop_",
                    "_rlnImageName #1",
                    "000001@stack.mrcs",
                    "",
                ]
            )
        )
    _write_processed_star(
        processed,
        [("000001@stack.mrcs", "0", "1")],
    )

    try:
        utils.backmap_stars(processed, mapping, out)
    except ValueError as e:
        msg = str(e)
        assert "_janas_source_rlnImageName" in msg, msg
        return
    raise AssertionError("missing source_tag in mapping must raise ValueError")


def test_conflicting_duplicate_mapping_raises(tmp: str) -> None:
    mapping = os.path.join(tmp, "mapping.star")
    processed = os.path.join(tmp, "processed.star")
    out = os.path.join(tmp, "out.star")

    _write_mapping_star(
        mapping,
        [
            ("000001@stack.mrcs", "001@orig.mrc"),
            ("000001@stack.mrcs", "999@other.mrc"),  # SAME key, DIFFERENT source
        ],
    )
    _write_processed_star(processed, [("000001@stack.mrcs", "0", "1")])

    try:
        utils.backmap_stars(processed, mapping, out)
    except ValueError as e:
        assert "more than one source name" in str(e), str(e)
        return
    raise AssertionError("conflicting duplicate mapping must raise ValueError")


def test_consistent_duplicate_mapping_is_tolerated(tmp: str) -> None:
    """Duplicates that all resolve to the same source are accepted and reported."""
    mapping = os.path.join(tmp, "mapping.star")
    processed = os.path.join(tmp, "processed.star")
    out = os.path.join(tmp, "out.star")

    _write_mapping_star(
        mapping,
        [
            ("000001@stack.mrcs", "001@orig.mrc"),
            ("000001@stack.mrcs", "001@orig.mrc"),  # same key, same value
            ("000002@stack.mrcs", "002@orig.mrc"),
        ],
    )
    _write_processed_star(processed, [("000001@stack.mrcs", "0", "1")])

    report = utils.backmap_stars(processed, mapping, out)
    df = _read_particles(out)
    assert list(df["_rlnImageName"]) == ["001@orig.mrc"]
    assert "000001@stack.mrcs" in report["duplicate_keys"]


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

TESTS: List[Tuple[str, Callable[[str], None]]] = [
    ("basic mapping ignores row order", test_basic_mapping_ignores_row_order),
    ("preserves metadata columns", test_preserves_metadata_columns),
    ("adds stack reference audit column", test_adds_stack_reference_audit_column),
    ("strict missing raises", test_strict_missing_raises),
    ("non-strict leaves missing unchanged", test_non_strict_leaves_missing_unchanged),
    ("missing source_tag in mapping raises", test_missing_source_tag_in_mapping_raises),
    ("conflicting duplicate mapping raises", test_conflicting_duplicate_mapping_raises),
    ("consistent duplicate mapping is tolerated", test_consistent_duplicate_mapping_is_tolerated),
]


def main() -> int:
    failed: List[str] = []
    for name, fn in TESTS:
        with tempfile.TemporaryDirectory() as tmp:
            try:
                fn(tmp)
                print(f"[ OK ] {name}")
            except Exception as exc:  # noqa: BLE001
                failed.append(name)
                print(f"[FAIL] {name}: {exc}")
                traceback.print_exc(limit=4)
    print()
    if failed:
        print(f"{len(failed)}/{len(TESTS)} test(s) failed: {failed}")
        return 1
    print(f"All {len(TESTS)} tests passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
