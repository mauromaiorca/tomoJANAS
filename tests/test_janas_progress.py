"""
Tests for janas.janas_progress: progress.html generator.

Run::

    python tests/test_janas_progress.py
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
import traceback
import types
from pathlib import Path
from typing import Callable, List, Tuple


HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.normpath(os.path.join(HERE, "..", "src"))
if os.path.isdir(SRC) and SRC not in sys.path:
    sys.path.insert(0, SRC)

if "janas.janas_core" not in sys.modules:
    sys.modules["janas.janas_core"] = types.ModuleType("janas.janas_core")

from janas import janas_progress as P  # noqa: E402


# ---------------------------------------------------------------------------
# Tiny session-dir factory
# ---------------------------------------------------------------------------


def _make_status(text_kwargs):
    lines = [
        "========================================",
        " JANAS Session Status",
        "========================================",
    ]
    for k, v in text_kwargs.items():
        lines.append(f" {k.capitalize()} : {v}")
    lines.append("========================================")
    return "\n".join(lines) + "\n"


def _make_session(tmp: Path, kind: str = "selection",
                  events=None, timings=None, status=None,
                  overview=None) -> Path:
    sd = tmp / "session_x"
    sd.mkdir(parents=True, exist_ok=True)
    if kind == "selection":
        (sd / "session_settings.toml").write_text("# settings\n", encoding="utf-8")
    elif kind == "classification":
        (sd / "session_classification_settings.txt").write_text(
            "# settings\n", encoding="utf-8"
        )
    if overview is not None:
        (sd / "overview.txt").write_text(overview, encoding="utf-8")
    rt = sd / "runtime"
    rt.mkdir(parents=True, exist_ok=True)
    if status is not None:
        (rt / "status.txt").write_text(_make_status(status), encoding="utf-8")
    if events is not None:
        with (rt / "events.ndjson").open("w", encoding="utf-8") as f:
            for ev in events:
                f.write(json.dumps(ev) + "\n")
    if timings is not None:
        with (rt / "step_timings.csv").open("w", encoding="utf-8") as f:
            f.write("iteration,step,t_start,t_end,elapsed_s,rc\n")
            for r in timings:
                f.write(",".join(r) + "\n")
    return sd


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_step_to_phase_known_steps() -> None:
    assert P._step_to_phase("01_randomize_halves") == 1
    assert P._step_to_phase("02_bootstrap_reconstruct") == 1
    assert P._step_to_phase("03_score_particles") == 2
    assert P._step_to_phase("04_subsets") == 3
    assert P._step_to_phase("05_subset_reconstructions") == 3
    assert P._step_to_phase("06_locres") == 4
    assert P._step_to_phase("06_locres_bulk") == 4
    assert P._step_to_phase("07_locres_stats") == 4
    assert P._step_to_phase("08_get_num_particles") == 4


def test_step_to_phase_unknown() -> None:
    assert P._step_to_phase("99_unknown") is None
    assert P._step_to_phase("") is None


def test_pick_stage_unstarted() -> None:
    stage = P._pick_stage(events=[], status={}, session_kind="selection")
    assert stage["state"] == "unstarted"
    assert stage["image"] == "selection_unstarted.png"


def test_pick_stage_session_started_no_step_yet_uses_step1_image() -> None:
    # Only session_start has been emitted (script has been launched but the
    # first step_start has not been written yet). The dashboard must show
    # the preprocessing image, not the unstarted one.
    events = [
        {"event": "session_start", "iteration": "0", "status": "started"},
    ]
    stage = P._pick_stage(events, status={}, session_kind="selection")
    assert stage["state"] == "running"
    assert stage["image"] == "selection_step1.png"


def test_pick_stage_status_only_no_events_uses_step1_image() -> None:
    # status.txt has been written (init_runtime_logging called) but the
    # events file has not been flushed yet from the caller's point of
    # view. Treat as preprocessing, not unstarted.
    status = {"iteration": "0", "step": "init", "status": "starting"}
    stage = P._pick_stage(events=[], status=status, session_kind="selection")
    assert stage["state"] == "running"
    assert stage["image"] == "selection_step1.png"


def test_pick_stage_running_maps_step_to_image() -> None:
    events = [
        {"event": "session_start", "iteration": "0", "status": "started"},
        {"event": "step_start", "iteration": "1", "step": "03_score_particles",
         "status": "running", "t_start": "2026-05-28T02:00:00Z"},
    ]
    stage = P._pick_stage(events, status={}, session_kind="selection")
    assert stage["state"] == "running"
    assert stage["image"] == "selection_step2.png"
    assert stage["current_step"] == "03_score_particles"
    assert stage["current_iter"] == "1"


def test_pick_stage_between_steps_keeps_last_phase_image() -> None:
    # The dashboard is regenerated AFTER step_end is written to
    # events.ndjson (see _runtime_logging_shell_block.run_step), so the
    # most common state when the HTML is built is "last event is step_end".
    # The image must reflect the phase of the just-completed step, not
    # silently fall back to the unstarted picture.
    events_step_end_03 = [
        {"event": "session_start", "iteration": "0", "status": "started"},
        {"event": "step_start", "iteration": "2",
         "step": "03_score_particles", "status": "running",
         "t_start": "2026-05-28T02:00:00Z"},
        {"event": "step_end", "iteration": "2",
         "step": "03_score_particles", "status": "success", "rc": "0",
         "elapsed_s": "120", "t_start": "2026-05-28T02:00:00Z",
         "t_end": "2026-05-28T02:02:00Z"},
    ]
    stage = P._pick_stage(events_step_end_03, status={}, session_kind="selection")
    assert stage["state"] == "running"
    assert stage["image"] == "selection_step2.png", stage
    assert "Step done" in stage["label"]

    # Same check for a step in phase 4 (locres family)
    events_step_end_07 = [
        {"event": "session_start", "iteration": "0", "status": "started"},
        {"event": "step_start", "iteration": "3",
         "step": "07_locres_stats", "status": "running",
         "t_start": "2026-05-28T02:10:00Z"},
        {"event": "step_end", "iteration": "3",
         "step": "07_locres_stats", "status": "success", "rc": "0",
         "elapsed_s": "12", "t_start": "2026-05-28T02:10:00Z",
         "t_end": "2026-05-28T02:10:12Z"},
    ]
    stage = P._pick_stage(events_step_end_07, status={}, session_kind="selection")
    assert stage["image"] == "selection_step4.png", stage


def test_pick_stage_finished() -> None:
    events = [
        {"event": "session_end", "status": "finished",
         "t_end": "2026-05-28T03:00:00Z"},
    ]
    stage = P._pick_stage(events, status={}, session_kind="selection")
    assert stage["state"] == "finished"
    assert stage["image"] == "selection_finished.png"


def test_pick_stage_aborted_uses_finished_image_with_aborted_state() -> None:
    events = [
        {"event": "session_end", "status": "aborted",
         "t_end": "2026-05-28T03:00:00Z"},
    ]
    stage = P._pick_stage(events, status={}, session_kind="selection")
    assert stage["state"] == "aborted"
    # Same picture but the badge state differentiates aborted from finished
    assert stage["image"] == "selection_finished.png"


def test_pick_stage_classification_has_no_image() -> None:
    events = [
        {"event": "step_start", "iteration": "1", "step": "03_score_particles",
         "status": "running", "t_start": "2026-05-28T02:00:00Z"},
    ]
    stage = P._pick_stage(events, status={}, session_kind="classification")
    assert stage["image"] == ""
    assert stage["state"] == "running"


def test_parse_status_txt() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        sd = _make_session(
            Path(tmp), status={"Iteration": "3", "Step": "03_score_particles",
                               "Status": "running"},
        )
        st = P._parse_status_txt(sd / "runtime" / "status.txt")
        assert st["iteration"] == "3"
        assert st["step"] == "03_score_particles"
        assert st["status"] == "running"


def test_iter_events_skips_malformed_lines() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        sd = Path(tmp) / "s"
        rt = sd / "runtime"
        rt.mkdir(parents=True, exist_ok=True)
        (rt / "events.ndjson").write_text(
            '{"event": "session_start"}\n'
            'this is not json\n'
            '{"event": "step_start", "step": "01_randomize_halves"}\n'
            '{"event": "step_end", "incomplete...\n',
            encoding="utf-8",
        )
        evs = P._iter_events(rt / "events.ndjson")
        assert len(evs) == 2
        assert evs[0]["event"] == "session_start"
        assert evs[1]["step"] == "01_randomize_halves"


def test_write_progress_html_creates_file_and_copies_images() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        sd = _make_session(
            Path(tmp),
            kind="selection",
            status={"Iteration": "2", "Step": "03_score_particles",
                    "Status": "running", "Host": "test-host"},
            events=[
                {"event": "session_start", "iteration": "0",
                 "status": "started", "hostname": "test-host"},
                {"event": "step_start", "iteration": "2",
                 "step": "03_score_particles", "status": "running",
                 "t_start": "2026-05-28T02:00:00Z",
                 "hostname": "test-host", "slurm_job_id": "12345"},
            ],
            timings=[
                ("1", "01_randomize_halves",
                 "2026-05-28T01:00:00Z", "2026-05-28T01:00:05Z", "5", "0"),
                ("1", "02_bootstrap_reconstruct",
                 "2026-05-28T01:00:05Z", "2026-05-28T01:01:00Z", "55", "0"),
            ],
            overview="iter sigma nParticles meanResolution\n1 1.00 100000 3.5\n",
        )
        out = P.write_progress_html(sd, refresh_seconds=10, max_events=50)
        assert out.exists()
        text = out.read_text(encoding="utf-8")
        # Sanity checks on the rendered HTML
        assert '<meta http-equiv="refresh" content="10">' in text
        assert "JANAS — session_x" in text or "session_x" in text
        assert "03_score_particles" in text
        assert "test-host" in text
        # SLURM and CUDA fields were removed from the page (JANAS does not
        # currently integrate with the cluster scheduler).
        assert ">SLURM job<" not in text
        assert ">SLURM nodes<" not in text
        assert ">SLURM CPUs<" not in text
        assert ">CUDA devices<" not in text
        assert "Scoring particles" in text
        assert "selection_step2.png" in text
        assert "01_randomize_halves" in text  # in timings table
        # Stage images were copied to runtime/imgs/
        imgs_dir = sd / "runtime" / "imgs"
        if P._PKG_IMAGES_DIR.is_dir():
            # In the dev tree the images ARE present; verify the copy
            assert (imgs_dir / "selection_step2.png").exists()


def test_write_progress_html_classification_session_text_only() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        sd = _make_session(
            Path(tmp),
            kind="classification",
            events=[
                {"event": "session_start", "iteration": "0",
                 "status": "started"},
            ],
        )
        out = P.write_progress_html(sd)
        text = out.read_text(encoding="utf-8")
        # 'Type' label is in the right card meta line; the value sits in
        # a <strong> tag so the substring is split across the tag.
        assert "Type:" in text
        assert ">classification<" in text
        # No image filename should appear in the page for classification
        assert "selection_step" not in text
        assert "selection_unstarted.png" not in text
        assert "selection_finished.png" not in text


def test_timing_table_uses_pass_fail_and_iteration_banding() -> None:
    """The Step timings table must:
      - use 'Return code' as the column header,
      - render rc=0 as 'PASS (rc=0)' (green) and non-zero as 'FAIL (rc=N)'
        (red),
      - band rows by iteration parity using iter-odd / iter-even.
    """
    with tempfile.TemporaryDirectory() as tmp:
        sd = _make_session(
            Path(tmp),
            timings=[
                ("1", "01_randomize_halves",
                 "2026-05-28T01:00:00Z", "2026-05-28T01:00:05Z", "5", "0"),
                ("1", "03_score_particles",
                 "2026-05-28T01:00:05Z", "2026-05-28T01:01:00Z", "55", "0"),
                ("2", "03_score_particles",
                 "2026-05-28T01:01:00Z", "2026-05-28T01:02:00Z", "60", "1"),
            ],
        )
        out = P.write_progress_html(sd)
        text = out.read_text(encoding="utf-8")

        # Column header
        assert "Return code" in text
        assert ">rc<" not in text  # the old short header is gone

        # Elapsed time header is centred, value format is human-readable
        assert 'class="elapsed-col">Elapsed time<' in text
        assert ">Elapsed (s)<" not in text   # the old short header is gone
        # Each elapsed cell is left-aligned and renders the days/hours/
        # mins/secs decomposition (the second row had elapsed_s=55).
        assert ">0</span> days," in text
        assert ">55</span> secs" in text

        # PASS / FAIL rendering
        assert "PASS (rc=0)" in text
        assert "FAIL (rc=1)" in text

        # Iteration banding classes are emitted on the <tr>
        assert "iter-odd" in text   # iteration 1 -> odd
        assert "iter-even" in text  # iteration 2 -> even

        # And the CSS rules backing them are present
        assert "tr.iter-odd  td" in text or "tr.iter-odd" in text
        assert "tr.iter-even td" in text or "tr.iter-even" in text


def test_overview_data_extraction_and_rendering() -> None:
    """End-to-end: overview.txt is parsed, iterations bar and particle
    counts appear in the HTML, the target iteration is highlighted."""
    with tempfile.TemporaryDirectory() as tmp:
        sd = Path(tmp) / "session_x"
        sd.mkdir(parents=True, exist_ok=True)
        (sd / "session_settings.toml").write_text("# settings\n", encoding="utf-8")
        (sd / "overview.txt").write_text(
            '[[_janas_target_selection]]\n'
            'reference_starFile = "selected.star"\n'
            'reference_num_particles = 45678\n'
            'selection_number = 2\n'
            '\n'
            '[[_janas_selection_0]]\n'
            'reference_starFile = "all.star"\n'
            'reference_num_particles = 123456\n'
            'selection_number = 0\n'
            '\n'
            '[[_janas_selection_1]]\n'
            'reference_starFile = "iter1.star"\n'
            'reference_num_particles = 90000\n'
            'selection_number = 1\n'
            '\n'
            '[[_janas_selection_2]]\n'
            'reference_starFile = "iter2.star"\n'
            'reference_num_particles = 45678\n'
            'selection_number = 2\n'
            '\n'
            '[[_janas_selection_3]]\n'
            'reference_starFile = "iter3.star"\n'
            'reference_num_particles = 60000\n'
            'selection_number = 3\n',
            encoding="utf-8",
        )
        # Minimal runtime/ so the rest of the renderer is happy
        (sd / "runtime").mkdir(parents=True, exist_ok=True)

        # Parser — iteration 0 (the reference) is now included in the bar
        data = P._read_overview_data(sd / "overview.txt")
        assert data["iterations"] == [0, 1, 2, 3], data
        assert data["target_iter"] == 2, data
        assert data["full_dataset_np"] == 123456, data
        assert data["target_np"] == 45678, data

        # Renderer
        out = P.write_progress_html(sd)
        text = out.read_text(encoding="utf-8")

        # Iterations bar: three numbers; '2' is the current/green one
        assert "Iterations:" in text
        assert ">1<" in text and ">2<" in text and ">3<" in text
        assert 'class="iter-num current">2<' in text
        # '1' and '3' are NOT marked as current
        assert 'class="iter-num">1<' in text
        assert 'class="iter-num">3<' in text

        # Particle counts with thousand separators
        assert "Full dataset:" in text
        assert "123,456" in text
        assert "Selection:" in text
        assert "45,678" in text


def test_overview_missing_skips_iterations_and_counts() -> None:
    """No overview.txt → neither the iterations bar nor the particle
    counts should appear in the HTML."""
    with tempfile.TemporaryDirectory() as tmp:
        sd = Path(tmp) / "session_no_overview"
        sd.mkdir(parents=True, exist_ok=True)
        (sd / "session_settings.toml").write_text("# settings\n", encoding="utf-8")
        (sd / "runtime").mkdir(parents=True, exist_ok=True)

        out = P.write_progress_html(sd)
        text = out.read_text(encoding="utf-8")
        assert "Iterations:" not in text
        assert "Full dataset:" not in text
        assert "Selection:" not in text


def test_overview_data_no_target_skips_highlight() -> None:
    """If `[[_janas_target_selection]]` is missing, the iterations bar
    must still render (with no green highlight) and particle counts
    fall back to the full dataset only."""
    with tempfile.TemporaryDirectory() as tmp:
        sd = Path(tmp) / "session_no_target"
        sd.mkdir(parents=True, exist_ok=True)
        (sd / "session_settings.toml").write_text("# settings\n", encoding="utf-8")
        (sd / "overview.txt").write_text(
            '[[_janas_selection_0]]\n'
            'reference_num_particles = 100000\n'
            '\n'
            '[[_janas_selection_1]]\n'
            'reference_num_particles = 80000\n',
            encoding="utf-8",
        )
        (sd / "runtime").mkdir(parents=True, exist_ok=True)

        data = P._read_overview_data(sd / "overview.txt")
        assert data["iterations"] == [0, 1]
        assert data["target_iter"] is None
        assert data["full_dataset_np"] == 100000
        assert data["target_np"] is None

        text = P.write_progress_html(sd).read_text(encoding="utf-8")
        assert "Iterations:" in text
        # No target → iteration 0 is highlighted as the current reference
        assert 'class="iter-num current">0<' in text
        assert "Full dataset:" in text
        assert "Selection:" not in text   # no target, no selection count


def test_default_refresh_is_fifteen_seconds() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        sd = _make_session(Path(tmp))
        text = P.write_progress_html(sd).read_text(encoding="utf-8")
        assert '<meta http-equiv="refresh" content="15">' in text
        assert "refreshed every 15s" in text


def test_finished_session_disables_meta_refresh() -> None:
    """Once the session is finished, the page should not auto-refresh."""
    with tempfile.TemporaryDirectory() as tmp:
        sd = _make_session(
            Path(tmp),
            events=[
                {"event": "session_end", "status": "finished",
                 "t_end": "2026-05-28T03:00:00Z"},
            ],
        )
        text = P.write_progress_html(sd).read_text(encoding="utf-8")
        assert "http-equiv=\"refresh\"" not in text
        assert "session finished, auto-refresh disabled" in text


def test_aborted_session_keeps_meta_refresh() -> None:
    """Aborted sessions may be restarted in place, so we keep the
    auto-refresh so the dashboard picks up the new session_start
    event when the user retries."""
    with tempfile.TemporaryDirectory() as tmp:
        sd = _make_session(
            Path(tmp),
            events=[
                {"event": "session_end", "status": "aborted",
                 "t_end": "2026-05-28T03:00:00Z"},
            ],
        )
        text = P.write_progress_html(sd).read_text(encoding="utf-8")
        assert '<meta http-equiv="refresh" content="15">' in text


def test_star_file_link_strips_session_dir_prefix() -> None:
    """Paths in overview.txt typically start with the session-dir name
    (because overview.txt is written from one level above). The link
    in progress.html must therefore have that prefix stripped so the
    href is relative to progress.html itself."""
    with tempfile.TemporaryDirectory() as tmp:
        sd = Path(tmp) / "janas_selection_example"
        sd.mkdir(parents=True)
        (sd / "session_settings.toml").write_text(
            "# settings\n", encoding="utf-8"
        )
        (sd / "overview.txt").write_text(
            '[[_janas_target_selection]]\n'
            'reference_starFile = "janas_selection_example/'
            '_janas_SCI__0.99_scored_selection_1/'
            'norm__janas_SCI__0.99_scored_selection_1_best8498.star"\n'
            'reference_num_particles = 8498\n'
            'selection_number = 1\n'
            '\n'
            '[[_janas_selection_0]]\n'
            'reference_num_particles = 9951\n'
            '\n'
            '[[_janas_selection_1]]\n'
            'reference_num_particles = 8498\n',
            encoding="utf-8",
        )
        (sd / "runtime").mkdir(parents=True)

        # Parser exposes the raw value
        data = P._read_overview_data(sd / "overview.txt")
        assert data["target_starfile"].endswith("best8498.star")

        # _starfile_relative_to_session strips the session-dir name
        rel = P._starfile_relative_to_session(
            data["target_starfile"], sd
        )
        assert rel == (
            "_janas_SCI__0.99_scored_selection_1/"
            "norm__janas_SCI__0.99_scored_selection_1_best8498.star"
        )

        # HTML carries an <a href=...> with the stripped path and the
        # basename as link text
        text = P.write_progress_html(sd).read_text(encoding="utf-8")
        assert (
            'href="_janas_SCI__0.99_scored_selection_1/'
            'norm__janas_SCI__0.99_scored_selection_1_best8498.star"'
            in text
        )
        assert ">norm__janas_SCI__0.99_scored_selection_1_best8498.star<" in text


def test_star_file_link_absolute_path_passes_through() -> None:
    rel = P._starfile_relative_to_session(
        "/abs/path/to/file.star", Path("/tmp/janas_selection_foo")
    )
    assert rel == "/abs/path/to/file.star"


def test_settings_html_generated_once_and_linked_from_progress() -> None:
    """First call to write_progress_html generates settings.html (from
    session_settings.toml) and the progress page links to it. A second
    call is idempotent — the existing settings.html is not overwritten."""
    with tempfile.TemporaryDirectory() as tmp:
        sd = Path(tmp) / "janas_selection_demo"
        sd.mkdir(parents=True)
        (sd / "session_settings.toml").write_text(
            'session_name = "demo"\n'
            'particles = "x.star"\n'
            'mpi = "80"\n'
            'gpu = "0 1"\n'
            'autoSigma = "True"\n'
            'sigma = "1.0"\n'
            'noExternalPrograms = "True"\n',
            encoding="utf-8",
        )
        (sd / "runtime").mkdir()

        # First call: generates settings.html and links to it
        P.write_progress_html(sd)
        assert (sd / "settings.html").exists()
        progress = (sd / "progress.html").read_text(encoding="utf-8")
        assert 'href="settings.html"' in progress
        # Raw TOML link should NOT be present once the HTML view exists
        assert 'href="session_settings.toml"' not in progress

        settings = (sd / "settings.html").read_text(encoding="utf-8")
        # Page must include the key/value table and the raw TOML block
        assert "session_name" in settings
        assert "mpi" in settings
        assert '<pre>' in settings
        # The raw-TOML block is HTML-escaped, so the literal quotes show
        # up as &quot;.
        assert 'sigma = &quot;1.0&quot;' in settings

        # Second call is idempotent: file is not overwritten
        first_mtime = (sd / "settings.html").stat().st_mtime
        # Sleep enough to detect mtime changes on coarse-resolution FS
        import time as _time
        _time.sleep(0.02)
        P.write_progress_html(sd)
        assert (sd / "settings.html").stat().st_mtime == first_mtime


def test_settings_link_falls_back_to_raw_toml_when_html_absent() -> None:
    """If somehow settings.html is missing but the .toml is there, the
    progress header still has a clickable link, pointing to the raw
    TOML."""
    with tempfile.TemporaryDirectory() as tmp:
        sd = Path(tmp) / "janas_selection_demo"
        sd.mkdir(parents=True)
        (sd / "session_settings.toml").write_text(
            'session_name = "demo"\n', encoding="utf-8"
        )
        (sd / "runtime").mkdir()

        # Bypass the helper that would generate settings.html — render the
        # HTML directly
        text = P._render_html(sd)
        assert 'href="session_settings.toml"' in text
        assert 'href="settings.html"' not in text


def test_settings_html_renders_booleans_with_classes() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        sd = Path(tmp) / "session_x"
        sd.mkdir(parents=True)
        # Real TOML booleans (not strings)
        (sd / "session_settings.toml").write_text(
            'autoSigma = true\n'
            'maskingCrop = false\n',
            encoding="utf-8",
        )
        out = P.write_settings_html(sd)
        text = out.read_text(encoding="utf-8")
        assert 'class="bool-true">true' in text
        assert 'class="bool-false">false' in text


def test_settings_html_skipped_when_no_toml() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        sd = Path(tmp) / "session_no_toml"
        sd.mkdir(parents=True)
        out = P.write_settings_html(sd)
        assert out is None
        assert not (sd / "settings.html").exists()


def test_progress_layout_session_info_in_header() -> None:
    """The Type/Host/Generated trio sits in the header (not in a
    separate card); the right card now hosts the Euler histograms."""
    with tempfile.TemporaryDirectory() as tmp:
        sd = _make_session(
            Path(tmp),
            events=[
                {"event": "session_start", "iteration": "0",
                 "status": "started", "hostname": "PC-587054"},
            ],
        )
        text = P.write_progress_html(sd).read_text(encoding="utf-8")
        assert 'class="grid grid-2"' in text
        # Type/Host/Generated now in the header
        assert "Type:" in text
        assert "PC-587054" in text
        assert "Generated:" in text
        # The old card headings are gone
        assert ">Runtime<" not in text
        assert ">Session info<" not in text
        # The new right card is present
        assert ">Euler angle distribution<" in text


def test_iterations_bar_shows_zero_highlighted_when_no_target() -> None:
    """With only _janas_selection_0 present (no real iterations yet), the
    bar shows just '0' in green."""
    with tempfile.TemporaryDirectory() as tmp:
        sd = Path(tmp) / "janas_selection_demo"
        sd.mkdir()
        (sd / "session_settings.toml").write_text("# settings\n", encoding="utf-8")
        (sd / "overview.txt").write_text(
            '[[_janas_selection_0]]\n'
            'reference_starFile = "demo/input.star"\n'
            'reference_num_particles = 9951\n',
            encoding="utf-8",
        )
        (sd / "runtime").mkdir()
        text = P.write_progress_html(sd).read_text(encoding="utf-8")
        assert "Iterations:" in text
        # Only the 0 chip is rendered, and it is highlighted
        assert 'class="iter-num current">0<' in text
        # No further iteration chips should appear
        for n in (1, 2, 3):
            assert f'>{n}<' not in text.split("Iterations:", 1)[1].split("</div>", 1)[0]


def test_iterations_bar_includes_zero_with_target_at_two() -> None:
    """0 1 2: iteration 2 is the target so it is the green one; 0 and 1
    stay neutral."""
    with tempfile.TemporaryDirectory() as tmp:
        sd = Path(tmp) / "janas_selection_demo"
        sd.mkdir()
        (sd / "session_settings.toml").write_text("# settings\n", encoding="utf-8")
        (sd / "overview.txt").write_text(
            '[[_janas_target_selection]]\n'
            'reference_starFile = "demo/iter2.star"\n'
            'reference_num_particles = 8498\n'
            'selection_number = 2\n'
            '\n'
            '[[_janas_selection_0]]\nreference_num_particles = 9951\n\n'
            '[[_janas_selection_1]]\nreference_num_particles = 9000\n\n'
            '[[_janas_selection_2]]\nreference_num_particles = 8498\n',
            encoding="utf-8",
        )
        (sd / "runtime").mkdir()
        text = P.write_progress_html(sd).read_text(encoding="utf-8")
        # 2 is the highlighted iteration; 0 and 1 are neutral
        assert 'class="iter-num">0<' in text
        assert 'class="iter-num">1<' in text
        assert 'class="iter-num current">2<' in text


def test_eulerhist_subprocess_is_invoked_for_input_and_target() -> None:
    """Calling write_progress_html should attempt to generate both euler
    histograms by spawning 'janas eulerHist --show False'. We patch
    subprocess.run to capture the calls."""
    import subprocess as _sp
    from unittest.mock import patch

    with tempfile.TemporaryDirectory() as tmp:
        sd = Path(tmp) / "janas_selection_demo"
        sd.mkdir()
        # An input star and a (different) target star, both present on disk
        input_star = sd / "input.star"
        input_star.write_text("# fake\n", encoding="utf-8")
        target_dir = sd / "_janas_SCI"
        target_dir.mkdir()
        target_star = target_dir / "best.star"
        target_star.write_text("# fake\n", encoding="utf-8")

        (sd / "session_settings.toml").write_text("# settings\n", encoding="utf-8")
        (sd / "overview.txt").write_text(
            '[[_janas_target_selection]]\n'
            f'reference_starFile = "janas_selection_demo/_janas_SCI/best.star"\n'
            'reference_num_particles = 8498\n'
            'selection_number = 1\n'
            '\n'
            '[[_janas_selection_0]]\n'
            'reference_starFile = "janas_selection_demo/input.star"\n'
            'reference_num_particles = 9951\n'
            '\n'
            '[[_janas_selection_1]]\n'
            'reference_starFile = "janas_selection_demo/_janas_SCI/best.star"\n'
            'reference_num_particles = 8498\n',
            encoding="utf-8",
        )
        (sd / "runtime").mkdir()

        # Patch subprocess.run AND make sure the resulting PNG exists so
        # the renderer believes the histogram was produced.
        def _fake_run(cmd, **kwargs):
            # The command structure is: ["janas", "eulerHist", "--i", ...,
            # "--outImage", out_png, "--show", "False"]
            try:
                out_idx = cmd.index("--outImage") + 1
                Path(cmd[out_idx]).write_text("fake png\n", encoding="utf-8")
            except (ValueError, IndexError, OSError):
                pass
            class _R: returncode = 0
            return _R()

        with patch("subprocess.run", side_effect=_fake_run) as run_mock:
            text = P.write_progress_html(sd).read_text(encoding="utf-8")

        # Two invocations: one for the input star, one for the target star.
        # Path.resolve() may rewrite /tmp -> /private/tmp on macOS, so we
        # compare on the resolved form on both sides.
        called_cmds = [call.args[0] for call in run_mock.call_args_list]
        input_resolved = str(input_star.resolve())
        target_resolved = str(target_star.resolve())
        assert any("eulerHist" in c and input_resolved in c
                   for c in called_cmds), called_cmds
        assert any("eulerHist" in c and target_resolved in c
                   for c in called_cmds), called_cmds
        # Both PNGs are referenced in the page
        assert "eulerhist_input.png" in text
        assert "eulerhist_target.png" in text


def test_eulerhist_skipped_when_png_newer_than_star() -> None:
    """Idempotency: if the PNG exists and its mtime is at least as new as
    the star, _ensure_eulerhist does not spawn the subprocess again."""
    from unittest.mock import patch
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        star = d / "x.star"
        star.write_text("hi", encoding="utf-8")
        png = d / "out.png"
        png.write_text("existing", encoding="utf-8")
        # Force png mtime ahead of the star
        import os as _os, time as _time
        _time.sleep(0.02)
        _os.utime(png, None)
        with patch("subprocess.run") as run_mock:
            out = P._ensure_eulerhist(star, png)
        assert out == png
        assert run_mock.call_count == 0


def test_resolve_star_path_uses_session_parent_for_relative_paths() -> None:
    """Regression: the run script's cwd is session_dir.parent, so a
    relative star path (typical of session_settings.toml.particles) must
    be resolved against session_dir.parent, not session_dir."""
    with tempfile.TemporaryDirectory() as tmp:
        parent = Path(tmp)
        sd = parent / "janas_selection_example"
        sd.mkdir()
        # The actual input star sits next to the session, in the parent.
        star = parent / "reference_subset.star"
        star.write_text("# fake\n", encoding="utf-8")

        # session_settings.toml-style relative path: no session-dir prefix
        out = P._resolve_star_path("reference_subset.star", sd)
        assert out is not None
        assert out.exists()
        assert out.resolve() == star.resolve()


def test_resolve_star_path_overview_style_with_session_prefix() -> None:
    """overview.txt records paths like 'janas_selection_example/foo.star'
    — must still resolve correctly via session_dir.parent."""
    with tempfile.TemporaryDirectory() as tmp:
        parent = Path(tmp)
        sd = parent / "janas_selection_example"
        (sd / "_janas_SCI").mkdir(parents=True)
        star = sd / "_janas_SCI" / "best.star"
        star.write_text("# fake\n", encoding="utf-8")

        out = P._resolve_star_path(
            "janas_selection_example/_janas_SCI/best.star", sd
        )
        assert out is not None
        assert out.exists()
        assert out.resolve() == star.resolve()


def test_eulerhist_input_star_resolved_from_session_settings() -> None:
    """End-to-end: with no overview.txt yet (early in the run), the input
    star comes from session_settings.toml.particles and must be resolved
    against the session parent so 'janas eulerHist' actually finds it."""
    from unittest.mock import patch

    with tempfile.TemporaryDirectory() as tmp:
        parent = Path(tmp)
        sd = parent / "janas_selection_example"
        sd.mkdir()
        # Input star next to the session, as in the real workflow
        input_star = parent / "reference_subset.star"
        input_star.write_text("# fake\n", encoding="utf-8")

        (sd / "session_settings.toml").write_text(
            'particles = "reference_subset.star"\n',
            encoding="utf-8",
        )
        # No overview.txt yet — runtime/ may not even exist
        (sd / "runtime").mkdir()

        def _fake_run(cmd, **kwargs):
            try:
                out_idx = cmd.index("--outImage") + 1
                Path(cmd[out_idx]).write_text("png\n", encoding="utf-8")
            except (ValueError, IndexError, OSError):
                pass
            class _R: returncode = 0
            return _R()

        with patch("subprocess.run", side_effect=_fake_run) as run_mock:
            text = P.write_progress_html(sd).read_text(encoding="utf-8")

        called_cmds = [call.args[0] for call in run_mock.call_args_list]
        # janas eulerHist must have been called with the ABSOLUTE path
        # to the star in the session parent, not session_dir.
        assert any(
            "eulerHist" in c and str(input_star.resolve()) in c
            for c in called_cmds
        ), called_cmds
        # And the page renders the input PNG
        assert "eulerhist_input.png" in text
        assert "Not available yet" not in text or "Current selection" in text


def test_eulerhist_subprocess_forwards_fontScale_when_not_default() -> None:
    """When ``font_scale != 1.0`` the subprocess command must contain
    '--fontScale <value>'."""
    from unittest.mock import patch
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        star = d / "x.star"; star.write_text("hi", encoding="utf-8")
        png = d / "out.png"
        with patch("subprocess.run") as run_mock:
            P._ensure_eulerhist(star, png, font_scale=2.0)
        cmd = run_mock.call_args.args[0]
        assert "--fontScale" in cmd
        idx = cmd.index("--fontScale")
        assert cmd[idx + 1] == "2.0"


def test_eulerhist_subprocess_omits_fontScale_when_default() -> None:
    from unittest.mock import patch
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        star = d / "x.star"; star.write_text("hi", encoding="utf-8")
        png = d / "out.png"
        with patch("subprocess.run") as run_mock:
            P._ensure_eulerhist(star, png)   # default 1.0
        cmd = run_mock.call_args.args[0]
        assert "--fontScale" not in cmd


def test_card_skips_current_selection_image_at_iteration_zero() -> None:
    """Iteration 0: no [[_janas_target_selection]] in overview.txt yet —
    the bottom slot should render a placeholder, not a duplicate of the
    input histogram."""
    from unittest.mock import patch
    with tempfile.TemporaryDirectory() as tmp:
        parent = Path(tmp)
        sd = parent / "janas_selection_demo"
        sd.mkdir()
        star = parent / "input.star"; star.write_text("# fake\n", encoding="utf-8")
        (sd / "session_settings.toml").write_text(
            'particles = "input.star"\n', encoding="utf-8")
        (sd / "runtime").mkdir()

        def _fake_run(cmd, **kwargs):
            try:
                out_idx = cmd.index("--outImage") + 1
                Path(cmd[out_idx]).write_text("png\n", encoding="utf-8")
            except (ValueError, IndexError, OSError):
                pass
            class _R: returncode = 0
            return _R()

        with patch("subprocess.run", side_effect=_fake_run) as run_mock:
            text = P.write_progress_html(sd).read_text(encoding="utf-8")

        # Only ONE eulerHist invocation: the input star
        ext_calls = [
            call for call in run_mock.call_args_list
            if "eulerHist" in call.args[0]
        ]
        assert len(ext_calls) == 1, ext_calls

        # The bottom slot still renders a heading, but with a placeholder
        # text and no <img>
        assert "Waiting for the first selection iteration" in text
        assert "eulerhist_target.png" not in text
        # Top slot is normal
        assert "eulerhist_input.png" in text
        # At iteration 0 the subhead is the prefix-only "Current selection"
        # (no '% of full dataset' suffix yet).
        assert ">Current selection<" in text
        assert "of full dataset" not in text.split("Current selection")[1].split("</div>")[0]


def test_particle_counts_show_percentage_of_full_dataset() -> None:
    """Selection: NNN particles (XX.X% of full dataset)"""
    with tempfile.TemporaryDirectory() as tmp:
        sd = Path(tmp) / "janas_selection_demo"
        sd.mkdir()
        (sd / "session_settings.toml").write_text("# settings\n", encoding="utf-8")
        (sd / "overview.txt").write_text(
            '[[_janas_target_selection]]\n'
            'reference_starFile = "demo/iter1.star"\n'
            'reference_num_particles = 5000\n'
            'selection_number = 1\n'
            '\n'
            '[[_janas_selection_0]]\nreference_num_particles = 20000\n\n'
            '[[_janas_selection_1]]\nreference_num_particles = 5000\n',
            encoding="utf-8",
        )
        (sd / "runtime").mkdir()
        text = P.write_progress_html(sd).read_text(encoding="utf-8")
        # 5000 / 20000 = 25.0%
        # Particle-counts line carries the percentage:
        assert "(25.0% of full dataset)" in text
        # ...and so does the bottom-slot subhead in the eulerhist card
        # ('Current selection — (25.0% of full dataset)').
        assert (
            "Current selection — (25.0% of full dataset)" in text
            or "Current selection &#x2014; (25.0% of full dataset)" in text
        )
        # The old generic 'Euler angle distribution' subhead is gone
        # from the bottom slot (still present on the top slot).
        assert "Current selection — Euler angle distribution" not in text


def test_particle_counts_no_percentage_when_full_dataset_zero() -> None:
    """Avoid ZeroDivisionError when the full-dataset count is missing/0."""
    with tempfile.TemporaryDirectory() as tmp:
        sd = Path(tmp) / "janas_selection_demo"
        sd.mkdir()
        (sd / "session_settings.toml").write_text("# settings\n", encoding="utf-8")
        (sd / "overview.txt").write_text(
            '[[_janas_target_selection]]\n'
            'reference_num_particles = 5000\n'
            'selection_number = 1\n',
            encoding="utf-8",
        )
        (sd / "runtime").mkdir()
        text = P.write_progress_html(sd).read_text(encoding="utf-8")
        # We still show the absolute count
        assert "5,000" in text
        # But no percentage line
        assert "of full dataset" not in text


def test_eulerhist_returns_none_when_star_missing() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        png = d / "out.png"
        out = P._ensure_eulerhist(d / "no_such.star", png)
        assert out is None
        assert not png.exists()


def test_format_elapsed_time_decomposition() -> None:
    # Pure unit decomposition: total = D*86400 + H*3600 + M*60 + S
    # 2 mins 32 secs
    assert ">2</span> mins, <span class=\"elapsed-num\">32</span> secs" in (
        P._format_elapsed_time(152)
    )
    # 1 hour exactly
    out = P._format_elapsed_time(3600)
    assert ">1</span> hours" in out
    assert ">0</span> mins" in out
    assert ">0</span> secs" in out
    # 1 day, 2 hours, 3 mins, 4 secs = 86400 + 7200 + 180 + 4 = 93784
    out = P._format_elapsed_time(93784)
    assert ">1</span> days" in out
    assert ">2</span> hours" in out
    assert ">3</span> mins" in out
    assert ">4</span> secs" in out
    # Strings from CSV are accepted
    out = P._format_elapsed_time("90")
    assert ">1</span> mins, <span class=\"elapsed-num\">30</span> secs" in out
    # Negative and non-numeric are handled gracefully
    out = P._format_elapsed_time(-5)
    assert ">0</span> secs" in out
    assert P._format_elapsed_time("") == '<span class="meta">--</span>'
    assert P._format_elapsed_time(None) == '<span class="meta">--</span>'


def test_atomic_write() -> None:
    """Re-writing must not leave a stray progress.html.tmp."""
    with tempfile.TemporaryDirectory() as tmp:
        sd = _make_session(Path(tmp))
        P.write_progress_html(sd)
        P.write_progress_html(sd)  # second call overwrites
        assert (sd / "progress.html").exists()
        assert not (sd / "progress.html.tmp").exists()


def test_html_is_escaped() -> None:
    """An injected '<' in events must not break the page."""
    with tempfile.TemporaryDirectory() as tmp:
        sd = _make_session(
            Path(tmp),
            events=[
                {"event": "step_start", "iteration": "1",
                 "step": "<script>alert('x')</script>",
                 "status": "running"},
            ],
        )
        out = P.write_progress_html(sd)
        text = out.read_text(encoding="utf-8")
        # The malicious step name must NOT appear verbatim as a tag
        assert "<script>alert('x')</script>" not in text
        # But its HTML-escaped form should
        assert "&lt;script&gt;" in text


# ---------------------------------------------------------------------------
# custom_selected_stacks/custom_selected_stacks.html
# ---------------------------------------------------------------------------


def _make_subset_dir(session_dir: Path, n: int, ite: int,
                     write_star: bool = True, write_png: bool = True) -> Path:
    folder = session_dir / "custom_selected_stacks" / f"subset_{n}_ite{ite}"
    folder.mkdir(parents=True, exist_ok=True)
    if write_star:
        (folder / f"subset_{n}_ite{ite}.star").write_text(
            "data_\nloop_\n_rlnImageName #1\nfoo.mrcs\n", encoding="utf-8"
        )
    if write_png:
        # Minimal valid 1x1 PNG (just for existence checks, not a real image).
        (folder / f"subset_{n}_ite{ite}_eulerhist.png").write_bytes(b"\x89PNG\r\n\x1a\n")
    return folder


def test_custom_stacks_index_html_created_empty() -> None:
    """With no subsets on disk, the index page must still be written and
    show the placeholder row plus the link back to progress.html."""
    with tempfile.TemporaryDirectory() as tmp:
        sd = _make_session(Path(tmp))
        out = P.write_custom_selected_stacks_html(sd)
        assert out is not None and out.exists()
        text = out.read_text(encoding="utf-8")
        # Empty-table placeholder row.
        assert '<tr class="empty-row">' in text
        assert "extract_custom_selected_stack.sh" in text
        # The folder is created even before any subset exists.
        assert (sd / "custom_selected_stacks").is_dir()
        # Header link back to progress.
        assert 'href="../progress.html"' in text


def test_custom_stacks_index_html_lists_subsets() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        sd = _make_session(Path(tmp))
        # Two subsets, different iterations and counts.
        _make_subset_dir(sd, n=3000, ite=1)
        _make_subset_dir(sd, n=5000, ite=2)
        # A stray non-matching dir must be silently ignored.
        (sd / "custom_selected_stacks" / "garbage").mkdir()
        out = P.write_custom_selected_stacks_html(sd)
        assert out is not None
        text = out.read_text(encoding="utf-8")
        # Both subsets must appear; the placeholder row must not.
        assert "subset_3000_ite1" in text
        assert "subset_5000_ite2" in text
        # The placeholder ROW must be absent (the CSS class itself is OK).
        assert '<tr class="empty-row">' not in text
        # Numbers formatted with thousands separator.
        assert "3,000" in text and "5,000" in text
        # The selected-subset thumbnail references the per-subdir PNG.
        assert 'src="subset_3000_ite1/subset_3000_ite1_eulerhist.png"' in text
        # The STAR link uses a relative path the browser can follow.
        assert 'href="subset_3000_ite1/subset_3000_ite1.star"' in text
        # Stray dir not listed.
        assert "garbage" not in text


def test_progress_html_links_to_custom_stacks() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        sd = _make_session(Path(tmp))
        out = P.write_progress_html(sd)
        text = out.read_text(encoding="utf-8")
        assert 'href="custom_selected_stacks/custom_selected_stacks.html"' in text


def test_custom_stacks_index_links_back_to_progress() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        sd = _make_session(Path(tmp))
        # Generate progress first so the index has the back-link to follow.
        P.write_progress_html(sd)
        index = sd / "custom_selected_stacks" / "custom_selected_stacks.html"
        assert index.exists()
        text = index.read_text(encoding="utf-8")
        assert 'href="../progress.html"' in text
        # Documentation pointer must also be present.
        assert "custom_selected_stacks.md" in text


def test_write_progress_html_creates_custom_stacks_index() -> None:
    """write_progress_html must always materialise the companion folder
    and custom_selected_stacks.html — even before any subset has been extracted."""
    with tempfile.TemporaryDirectory() as tmp:
        sd = _make_session(Path(tmp))
        P.write_progress_html(sd)
        assert (sd / "custom_selected_stacks").is_dir()
        assert (sd / "custom_selected_stacks" / "custom_selected_stacks.html").exists()


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

TESTS: List[Tuple[str, Callable[[], None]]] = [
    ("step → phase: known steps map correctly", test_step_to_phase_known_steps),
    ("step → phase: unknown steps return None", test_step_to_phase_unknown),
    ("pick stage: unstarted by default", test_pick_stage_unstarted),
    ("pick stage: session_start only -> step1 (preprocessing)",
        test_pick_stage_session_started_no_step_yet_uses_step1_image),
    ("pick stage: status only (no events) -> step1 (preprocessing)",
        test_pick_stage_status_only_no_events_uses_step1_image),
    ("pick stage: running maps step to image", test_pick_stage_running_maps_step_to_image),
    ("pick stage: between steps keeps last phase image",
        test_pick_stage_between_steps_keeps_last_phase_image),
    ("pick stage: finished", test_pick_stage_finished),
    ("pick stage: aborted uses finished image with aborted badge",
        test_pick_stage_aborted_uses_finished_image_with_aborted_state),
    ("pick stage: classification has no image", test_pick_stage_classification_has_no_image),
    ("parse status.txt", test_parse_status_txt),
    ("iter_events skips malformed JSON lines", test_iter_events_skips_malformed_lines),
    ("write progress.html + copies stage images",
        test_write_progress_html_creates_file_and_copies_images),
    ("classification session: text-only HTML",
        test_write_progress_html_classification_session_text_only),
    ("timing table: 'Return code' header, PASS/FAIL cells, iteration banding",
        test_timing_table_uses_pass_fail_and_iteration_banding),
    ("overview.txt: iterations bar + particle counts + target highlight",
        test_overview_data_extraction_and_rendering),
    ("overview.txt missing: bar and counts skipped",
        test_overview_missing_skips_iterations_and_counts),
    ("overview.txt without [[_janas_target_selection]]: bar without highlight",
        test_overview_data_no_target_skips_highlight),
    ("default refresh is 15 seconds",
        test_default_refresh_is_fifteen_seconds),
    ("finished session disables meta refresh",
        test_finished_session_disables_meta_refresh),
    ("aborted session keeps meta refresh",
        test_aborted_session_keeps_meta_refresh),
    ("STAR link: strip session-dir prefix",
        test_star_file_link_strips_session_dir_prefix),
    ("STAR link: absolute path is passed through",
        test_star_file_link_absolute_path_passes_through),
    ("settings.html: generated once + linked from progress header",
        test_settings_html_generated_once_and_linked_from_progress),
    ("settings link falls back to raw TOML when HTML absent",
        test_settings_link_falls_back_to_raw_toml_when_html_absent),
    ("settings.html: booleans rendered with bool-true/bool-false classes",
        test_settings_html_renders_booleans_with_classes),
    ("settings.html: skipped when no session_settings.toml",
        test_settings_html_skipped_when_no_toml),
    ("progress layout: Type/Host/Generated moved to header",
        test_progress_layout_session_info_in_header),
    ("iterations bar: '0' alone, highlighted, when no target yet",
        test_iterations_bar_shows_zero_highlighted_when_no_target),
    ("iterations bar: '0 1 2' with target=2 highlighted",
        test_iterations_bar_includes_zero_with_target_at_two),
    ("eulerhist: subprocess invoked for input + target stars",
        test_eulerhist_subprocess_is_invoked_for_input_and_target),
    ("eulerhist: skipped when PNG newer than star (idempotent)",
        test_eulerhist_skipped_when_png_newer_than_star),
    ("eulerhist: returns None when star file missing",
        test_eulerhist_returns_none_when_star_missing),
    ("eulerhist subprocess: --fontScale forwarded when != 1.0",
        test_eulerhist_subprocess_forwards_fontScale_when_not_default),
    ("eulerhist subprocess: --fontScale omitted when default",
        test_eulerhist_subprocess_omits_fontScale_when_default),
    ("card: at iteration 0 the bottom slot is a placeholder, not an image",
        test_card_skips_current_selection_image_at_iteration_zero),
    ("particle counts: '(XX.X% of full dataset)' is appended",
        test_particle_counts_show_percentage_of_full_dataset),
    ("particle counts: skip percentage when full dataset is 0/missing",
        test_particle_counts_no_percentage_when_full_dataset_zero),
    ("resolve_star_path: relative path uses session parent (regression)",
        test_resolve_star_path_uses_session_parent_for_relative_paths),
    ("resolve_star_path: overview-style prefix still resolves",
        test_resolve_star_path_overview_style_with_session_prefix),
    ("eulerhist: input star from session_settings.particles, no overview",
        test_eulerhist_input_star_resolved_from_session_settings),
    ("elapsed time: 'days, hours, mins, secs' decomposition + edge cases",
        test_format_elapsed_time_decomposition),
    ("atomic write leaves no .tmp file", test_atomic_write),
    ("HTML escapes injected step names", test_html_is_escaped),
    ("custom_selected_stacks: index.html created with empty table by default",
        test_custom_stacks_index_html_created_empty),
    ("custom_selected_stacks: index.html lists discovered subset_*_ite* subdirs",
        test_custom_stacks_index_html_lists_subsets),
    ("custom_selected_stacks: progress.html links to index.html",
        test_progress_html_links_to_custom_stacks),
    ("custom_selected_stacks: index.html links back to progress.html",
        test_custom_stacks_index_links_back_to_progress),
    ("custom_selected_stacks: write_progress_html ensures folder + index.html",
        test_write_progress_html_creates_custom_stacks_index),
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
            traceback.print_exc(limit=4)
    print()
    if failed:
        print(f"{len(failed)}/{len(TESTS)} test(s) failed: {failed}")
        return 1
    print(f"All {len(TESTS)} tests passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
