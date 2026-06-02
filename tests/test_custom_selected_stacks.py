"""
Tests for the custom_selected_stacks feature.

Two surfaces are exercised:

1. ``janas selectBestRanked --exact`` — the new flag must cause
   ``extractBest`` (and the CLI wrapper) to refuse silently truncating
   to fewer than the requested ``N`` particles.

2. The generated ``extract_custom_selected_stack.sh`` script — the
   session manager must:
     a. embed the iteration's ``ite`` and ``tag`` as literals,
     b. promote the script to the session root when an iteration
        becomes the new ``_janas_target_selection``.

The bash side is only checked at the string-template level; we do not
spawn bash here.

Run with::

    python tests/test_custom_selected_stacks.py
"""

from __future__ import annotations

import io
import os
import sys
import traceback
import types
from typing import Callable, List, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.normpath(os.path.join(HERE, "..", "src"))
if os.path.isdir(SRC) and SRC not in sys.path:
    sys.path.insert(0, SRC)

# Stub the C++ extension so we can import without building.
if "janas.janas_core" not in sys.modules:
    sys.modules["janas.janas_core"] = types.ModuleType("janas.janas_core")

from janas import starHandler  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_synthetic_star(path: str, n_rows: int, rank_tag: str) -> None:
    """Write a minimal RELION-style STAR file with N rows and a ranking col."""
    header = [
        "",
        "data_particles",
        "",
        "loop_",
        "_rlnImageName #1",
        f"{rank_tag} #2",
    ]
    body = [f"img_{i:05d}.mrcs  {n_rows - i}" for i in range(n_rows)]
    with open(path, "w") as f:
        f.write("\n".join(header + body) + "\n")


# ---------------------------------------------------------------------------
# --exact tests
# ---------------------------------------------------------------------------


def test_extractBest_exact_ok(tmpdir: str) -> None:
    rank_tag = "_janas_SCI__0.99_scored_selection_1_norm7846"
    src = os.path.join(tmpdir, "src.star")
    dst = os.path.join(tmpdir, "dst.star")
    _write_synthetic_star(src, n_rows=10, rank_tag=rank_tag)
    # exact=True with N==available: must succeed.
    starHandler.extractBest(src, dst, 10, rank_tag, exact=True)
    assert os.path.isfile(dst), "extractBest should write the output STAR"


def test_extractBest_exact_raises(tmpdir: str) -> None:
    rank_tag = "_janas_SCI__0.99_scored_selection_1_norm7846"
    src = os.path.join(tmpdir, "src.star")
    dst = os.path.join(tmpdir, "dst.star")
    _write_synthetic_star(src, n_rows=5, rank_tag=rank_tag)
    try:
        starHandler.extractBest(src, dst, 100, rank_tag, exact=True)
    except ValueError as exc:
        assert "100" in str(exc) and "5" in str(exc), (
            f"ValueError should mention requested and available counts; got: {exc}"
        )
        return
    raise AssertionError("extractBest(exact=True) should raise when N > available")


def test_extractBest_non_exact_truncates(tmpdir: str) -> None:
    """Without --exact, the historical truncation behaviour must be preserved."""
    rank_tag = "_janas_SCI__0.99_scored_selection_1_norm7846"
    src = os.path.join(tmpdir, "src.star")
    dst = os.path.join(tmpdir, "dst.star")
    _write_synthetic_star(src, n_rows=5, rank_tag=rank_tag)
    # exact=False (default): must succeed silently even though N > available.
    starHandler.extractBest(src, dst, 100, rank_tag, exact=False)
    assert os.path.isfile(dst)


def test_selectBestRanked_cli_has_exact_flag() -> None:
    from janas import janas_cmd_caller as C
    args = C.janas_parser.parse_args(
        ["selectBestRanked", "--i", "x.star", "--num", "10", "--exact"]
    )
    assert args.exact is True, "--exact must parse to True"

    args2 = C.janas_parser.parse_args(
        ["selectBestRanked", "--i", "x.star", "--num", "10"]
    )
    assert args2.exact is False, "--exact must default to False"


# ---------------------------------------------------------------------------
# Generated-script template tests
# ---------------------------------------------------------------------------


def _session_manager_source() -> str:
    """Return the session manager source as a string (no execution)."""
    path = os.path.join(SRC, "janas", "janas_cmd_session_manager.py")
    with open(path) as f:
        return f.read()


def test_session_manager_emits_extract_script_block() -> None:
    src = _session_manager_source()
    # The cat-heredoc that materialises the script must be present.
    assert "extract_custom_selected_stack.sh" in src
    assert "<<EXTRACT_EOF" in src
    assert "janas selectBestRanked --i " in src
    assert "--exact" in src


def test_session_manager_promotes_on_target_match() -> None:
    src = _session_manager_source()
    # The promotion hook compares current_target_ID to ${ite} and copies.
    assert 'getTarget --overviewFile "${workingDir}/overview.txt" --current_selection_ID' in src
    assert 'cp -f "${extractScript}" "${workingDir}/extract_custom_selected_stack.sh"' in src


def test_session_manager_eulerhist_uses_fontscale_2() -> None:
    src = _session_manager_source()
    # Inside the extractor heredoc, the eulerHist invocation must request
    # --fontScale 2.0 and --show False (the headless-safe combination).
    assert "--fontScale 2.0" in src
    assert "--show False" in src


def test_session_manager_n_is_mandatory() -> None:
    src = _session_manager_source()
    # The script must error out when no positional N is provided.
    assert 'if [ \\$# -lt 1 ]; then' in src


def test_session_manager_refreshes_progress_after_extract() -> None:
    """Each invocation of the extractor must refresh the companion HTML
    (custom_selected_stacks.html + progress.html) by calling
    ``janas_optimizer progress --quiet`` at the end."""
    src = _session_manager_source()
    assert 'janas_optimizer progress --session "\\${SESSION_DIR}" --quiet' in src


# ---------------------------------------------------------------------------
# Doc + README integration
# ---------------------------------------------------------------------------


def test_readme_links_custom_selected_stacks_between_rows() -> None:
    repo = os.path.normpath(os.path.join(HERE, ".."))
    with open(os.path.join(repo, "README.md")) as f:
        readme = f.read()
    assert "docs/custom_selected_stacks.md" in readme
    # Must sit between iterative selection and class reassignment.
    i_iter = readme.find("docs/ITERATIVE_SELECTION.md")
    i_custom = readme.find("docs/custom_selected_stacks.md")
    i_class = readme.find("docs/CLASS_REASSIGNMENT.md")
    assert -1 < i_iter < i_custom < i_class, (
        "custom_selected_stacks row must sit between iterative selection "
        f"and class reassignment (got iter={i_iter}, custom={i_custom}, class={i_class})"
    )


def test_doc_mentions_repicking_and_workflow() -> None:
    repo = os.path.normpath(os.path.join(HERE, ".."))
    with open(os.path.join(repo, "docs", "custom_selected_stacks.md")) as f:
        doc = f.read()
    # Reference to JANAS-based repicking section of the article.
    assert "JANAS-based repicking" in doc
    # Workflow: selection first, then extract best, then run the script.
    assert "Iterative particle selection" in doc or "particle selection" in doc.lower()
    assert "./extract_custom_selected_stack.sh" in doc
    assert "subset_" in doc and "ite" in doc


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------


def _with_tmp(fn):
    def wrapper():
        import tempfile, shutil
        td = tempfile.mkdtemp(prefix="janas_test_custom_")
        try:
            fn(td)
        finally:
            shutil.rmtree(td, ignore_errors=True)
    return wrapper


TESTS: List[Tuple[str, Callable[[], None]]] = [
    ("extractBest(exact=True): success when N == available",
        _with_tmp(test_extractBest_exact_ok)),
    ("extractBest(exact=True): raises ValueError when N > available",
        _with_tmp(test_extractBest_exact_raises)),
    ("extractBest(exact=False): preserves historical truncation",
        _with_tmp(test_extractBest_non_exact_truncates)),
    ("selectBestRanked CLI exposes --exact flag",
        test_selectBestRanked_cli_has_exact_flag),
    ("session manager emits extract_custom_selected_stack.sh block",
        test_session_manager_emits_extract_script_block),
    ("session manager promotes script to session root on target match",
        test_session_manager_promotes_on_target_match),
    ("generated eulerHist call uses --fontScale 2.0 --show False",
        test_session_manager_eulerhist_uses_fontscale_2),
    ("generated script makes N mandatory",
        test_session_manager_n_is_mandatory),
    ("generated script refreshes companion HTML via janas_optimizer progress",
        test_session_manager_refreshes_progress_after_extract),
    ("README links custom_selected_stacks between iterative selection and 3D class reassignment",
        test_readme_links_custom_selected_stacks_between_rows),
    ("doc references JANAS-based repicking, workflow and script usage",
        test_doc_mentions_repicking_and_workflow),
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
