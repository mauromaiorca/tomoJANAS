# Changelog

## 2.1.5

- Fix `janas eulerHist --show False` so it actually skips
  `plt.show()`. The argument was declared as `type=bool` in argparse,
  which is broken: `bool("False")` returns `True` (any non-empty
  string is truthy). The new `_str2bool` helper accepts
  `true/false/yes/no/1/0/on/off` (case-insensitive) and raises
  `ArgumentTypeError` for anything else.
- Make `plotRoundEulerHist` headless-safe. When `toShow=False` the
  figure is built via `matplotlib.figure.Figure` directly (no GUI
  backend init), so `janas eulerHist --outImage out.png --show False`
  now works on headless nodes and on SSH sessions with broken X11
  forwarding (same pattern as `janas_optimizer.predict_min_particles`).
- Cast `--maxValue` from CLI string to float before passing to
  `vmax`.
- New test file `tests/test_eulerhist_show_flag.py` with 6 cases.

## 2.1.4

- Restore the two-column layout on `progress.html`. The Current stage
  card (image + iterations bar + step info + particle counts) now
  sits in the left column at a `1.4 fr` width, and the stage image
  itself is capped at 360 px so it stays compact even on wide
  monitors.
- The right column carries a new "Session info" card with:
  - the `Type`, `Host`, `Generated` meta line that briefly lived in
    the page header in 2.1.3;
  - the full contents of `[[_janas_target_selection]]` from
    `overview.txt`, rendered as a sorted key/value table inside a
    scrollable container so a tall selection record cannot push the
    rest of the card off-screen.
- The page header is trimmed back to just the H1 and the
  `Session directory: ... Settings: ...` line.
- The right card is intentionally left "thin-but-stable" for now —
  more useful widgets (parameter highlights, mini-plots) can be
  added in subsequent versions without re-shuffling the layout.
- Three new tests: two-column layout sanity, target-selection block
  rendered inside a `scroll-area`, placeholder when the block is
  missing.

## 2.1.3

- `progress.html` no longer auto-refreshes once the session has been
  marked `finished`: the `<meta http-equiv="refresh">` tag is dropped
  and the footer reads "session finished, auto-refresh disabled".
  Aborted sessions still refresh so the dashboard picks up a restart
  in place.
- Default refresh cadence changed from 10 s to 15 s
  (`DEFAULT_REFRESH_SECONDS = 15` in code, `--refresh` default in CLI,
  documentation kept in sync). The explicit `--refresh N` override is
  unchanged.
- Render a clickable link to the selected STAR file in the particle
  counts block. Source: `reference_starFile` in
  `[[_janas_target_selection]]` (with the fall-back to the matching
  selection block when the target does not carry it). The path stored
  in `overview.txt` typically starts with the session-directory name
  (because `overview.txt` is written from one level above); the new
  helper `_starfile_relative_to_session` strips that prefix so the
  href is relative to `progress.html` and works equally well under
  `file://` and under `python -m http.server`. Absolute paths and
  paths that do not start with the session name are passed through
  unchanged.
- Remove the Runtime card from `progress.html` and the SLURM
  (job/nodes/CPUs) and `CUDA_VISIBLE_DEVICES` fields it carried:
  JANAS does not currently integrate with the cluster scheduler, so
  on most installs those rows were empty. The hostname is now
  surfaced inline in the page header as `Host: <name>`. The values
  are still captured in `runtime/events.ndjson` for downstream
  tooling.
- New `settings.html` companion page next to `progress.html`. It is
  rendered once per session (idempotent) from
  `session_settings.toml`. The page presents the parameters as a
  sorted key/value table — booleans coloured green/red, lists shown
  comma-separated — plus the raw TOML in a `<pre>` block for
  copy-paste. `progress.html` links to it in the header as
  `Settings: session_settings`. When `session_settings.toml` exists
  but `settings.html` has not been generated yet, the header link
  falls back to the raw `.toml` so the user always has something
  clickable.
- Step timings table: the column header `Elapsed (s)` is renamed to
  `Elapsed time` (centered), and each cell now reads
  `X days, Y hours, Z mins, W secs` (left-aligned) instead of the
  raw integer second count. The numeric portion of each unit is
  wrapped in a fixed-width inline-block of tabular numerals so
  successive rows visually line up.
- Eleven new tests in total covering: the new 15 s default, the
  meta-refresh drop on finished, the meta-refresh kept on aborted,
  the prefix-stripped relative href, the absolute-path fall-through,
  the settings.html generation + linking, the raw-TOML fallback link,
  boolean rendering classes, the no-TOML skip, the host-in-header
  move, and the elapsed-time decomposition.

## 2.1.2

- Show selection-iteration progress and dataset sizes in `progress.html`,
  driven directly from `overview.txt`:
  - **Iterations bar** under the current-stage image: one numbered chip
    per `[[_janas_selection_N]]` block with `N > 0`. The iteration
    referenced by `selection_number` in `[[_janas_target_selection]]`
    is highlighted in green; all others are neutral.
  - **Particle counts** below the step-info line:
    `Full dataset: NNN,NNN particles · Selection: NNN,NNN particles`,
    sourced from `_janas_selection_0.reference_num_particles` and
    `_janas_target_selection.reference_num_particles` (with a
    fall-back to the selection block pointed to by `selection_number`
    when the target block does not carry the count itself).
- Both sections degrade gracefully: if `overview.txt` is absent, or
  contains no selections, or has no target, the relevant block is
  omitted instead of rendered empty.
- Three new tests covering the full extraction + rendering, the
  missing-overview case and the missing-target case.

## 2.1.1

- Show `selection_step1.png` (preprocessing) as soon as the run script
  has started instead of keeping `selection_unstarted.png` until the
  first scientific step emits `step_start`. `selection_unstarted.png`
  is now reserved for sessions that have not been launched at all
  (no `status.txt`, no events). This covers two early states:
  - only `runtime/status.txt` exists (init_runtime_logging has just
    written its initial status but the events file has not been
    flushed yet from the caller's point of view), and
  - `events.ndjson` contains only `session_start` and no `step_start`
    has been emitted yet.
- **Patch (no version bump):** also keep the just-completed phase
  image between steps. The `run_step` shell helper regenerates
  `progress.html` right after `step_end` is written, so the most
  common state when the HTML is built is "last event is `step_end`".
  Previously the `step_end` branch in `_pick_stage` updated the
  label to "Step done (rc=…), awaiting next" but left the image
  untouched, so it silently reverted to the unstarted picture
  between every step. Now the image stays on the phase of the
  just-completed step until the next `step_start` arrives.
- **Patch (no version bump):** Step timings table UI polish:
  - column header renamed from "rc" to "Return code";
  - each cell rendered as `PASS (rc=0)` (green) or `FAIL (rc=N)`
    (red) instead of just the raw integer;
  - alternating row backgrounds keyed off iteration parity
    (`iter-odd` light-grey, `iter-even` white) so it is easy to
    visually separate one iteration from the next.
- Four new tests covering the early-startup transitions, the
  between-steps regression, and the timing-table UI guarantees.

## 2.1.0

- New subcommand `janas_optimizer progress` that writes a single
  `progress.html` in the session directory summarising
  `runtime/events.ndjson`, `runtime/status.txt`,
  `runtime/step_timings.csv` and `overview.txt`. Includes a visual
  "current stage" card (one of 6 bundled PNGs for the selection
  workflow), per-step timings, recent events, resource info
  (host / SLURM job / CUDA devices) and an `<meta http-equiv="refresh">`
  for auto-reload. Classification sessions render a text-only variant
  (no image card) so they work without per-step illustrations.
- The 6 stage images ship inside the package
  (`src/janas/images/selection_*.png`) and are copied once into
  `<session>/runtime/imgs/`.
- The runtime-logging shell helpers in the generated run script now
  call `janas_optimizer progress --quiet` after every state change
  (`init_runtime_logging`, end of every `run_step`,
  `finish_runtime_logging`, and the abort EXIT trap), so the dashboard
  stays fresh during the session without any manual polling. Each call
  is `|| true`-guarded and silenced, so dashboard errors can never
  abort the scientific run.

## 2.0.0

- **Headless-safe plotting in `janas_optimizer`.** `predict_min_particles` and
  `plot_standard_whiskers_with_side_table` now build the figure via
  `matplotlib.figure.Figure` when no interactive display is requested
  (`--plot` not set, or `show=False` with an output file). `Figure()` is a
  pure in-memory object and never initialises a GUI backend, so the steps
  that save PNG/PDF complete cleanly on headless nodes and on SSH sessions
  with broken X11 forwarding (`_tkinter.TclError: couldn't connect to
  display ...`). Interactive paths (`--plot`) still go through
  `pyplot.subplots()` / `pyplot.figure()` so `plt.show()` can manage the
  window through the configured backend.
- Replace the runtime backend-swap fallback added in 1.0.x
  (`_safe_make_figure`) with `_make_figure_and_axes`, a small helper that
  routes batch callers to `Figure()` and interactive callers to
  `pyplot.subplots()`.
- Switch `plt.savefig` to `fig.savefig` in `predict_min_particles` so the
  PNG output works on a bare Figure.
- The headless guard in `janas/__init__.py` (`MPLBACKEND=Agg` when neither
  `MPLBACKEND` nor `DISPLAY` is set) is retained as defence-in-depth for
  other modules that still call `pyplot` directly.

## 1.0.8

- Add `utils.backmap_stars()` — the inverse companion of `create_stack_from_star`. Restores the original source `_rlnImageName` in a downstream STAR file by joining against the stack-generation STAR (which carries `_janas_source_rlnImageName`), preserving all refined metadata. Optionally writes `_janas_stack_rlnImageName` as an audit column so the rewrite is reversible.
- Expose it as `janas_utils backmap_stars --processed P.star --mapping S.star --output O.star` (`--no-strict`, `--stack-reference-tag ""`, `--section-name` available).
- Add `tests/test_backmap_stars.py` covering: row-order-independent mapping, metadata preservation, audit column, strict vs non-strict missing keys, missing provenance column in the mapping STAR, conflicting duplicate keys, tolerated consistent duplicate keys.

## 1.0.7

- Declare `cmake>=3.10` as a build-system requirement in `pyproject.toml`. Some Python environments (notably colabfold's bundled conda) ship a Python wrapper at `<env>/bin/cmake` that depends on the `cmake` Python package without installing it, breaking the build with `ModuleNotFoundError: No module named 'cmake'`. With this change, pip's default build isolation will install a working CMake before invoking our build.

## 1.0.6

- Add a `gpu` extra so `pip install 'janas[gpu]'` pulls a generic PyTorch wheel. Users who need a specific CUDA build should still install torch manually via `--index-url`.
- `janas_reconstructor --gpu` now fails fast with a clear error message and install hints (including the `--index-url cu128 / cu121 / cu118` examples and `pip install 'janas[gpu]'`) when PyTorch is not available, instead of crashing inside torch-using code paths.

## 1.0.5

- Revert the default CTF application mode for particle scoring from `modulate` back to `phaseflip`. This restores the manuscript-described behaviour for `janas scoreParticles`, `janas_session_manager new_select_session`, and `janas_session_manager classification_session`. `modulate` and `wiener` remain available via `--ctf-mode`.

## 1.0.4

- Add `--noRecs` option to `janas_session_manager classification_session`. When set, the generated run script skips all per-class reconstructions and only performs scoring and class assignment. Per-class STAR files are still written, so the user can reconstruct each class independently with their preferred software (RELION, cryoSPARC, `janas_reconstructor`, etc.).

## 1.0.3

- Change the default CTF application mode for particle scoring from `phaseflip` to `modulate` (multiply by the full CTF). This applies to `janas scoreParticles`, `janas_session_manager classification_session`, and `janas_session_manager new_select_session`.
- Fix `new_select_session --ctf-mode`: align choices with `janas scoreParticles` (`modulate`, `phaseflip`, `wiener`) — previous choices (`none`, `image`, `phaseflip`, `ref`) did not match the scoring backend.
- The selected CTF mode is now actually propagated to all `janas scoreParticles` invocations in the generated selection run script.

## 1.0.2

- Add `--ctf-mode` option to `janas_session_manager classification_session`. Choices: `modulate`, `phaseflip` (default), `wiener`. The flag is propagated to the `janas scoreParticles` calls in the generated run script, so the chosen CTF handling is applied during class scoring.

## 1.0.1

- Wrap `os.chmod` calls in `try / except PermissionError` for Windows compatibility (some filesystems do not permit `chmod`).

## 1.0.0

- First public release on PyPI.
- Integrated C++ app compilation into `pip install` (no separate CMake step).
- Added `janas_utils csparc2star-stack` for direct cryoSPARC particle import.
- Added `janas_utils update_from_csparc` for updating STAR metadata from cryoSPARC jobs.
- Added GPU-accelerated 3D reconstruction (`janas_reconstructor`).
- Added `--noExternalPrograms` mode for self-contained operation without RELION.
- Renamed project from emprove to JANAS. Backward compatibility with `_emprove_` STAR tags is preserved.

## 0.1.3.x

- Development releases under the name emprove.
- Iterative particle selection and 3D class reassignment workflows.
- CryoSPARC integration for local NU-refinement.
- Local resolution estimation (`locresBulk`).
