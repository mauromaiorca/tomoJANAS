# --------------------------------------------------------------------------
# Headless matplotlib defaults
# --------------------------------------------------------------------------
# JANAS scripts (especially those produced by janas_session_manager) run in
# batch mode and save plots to PNG without ever opening a GUI window. The
# default matplotlib backend on most Linux installs is TkAgg, which tries to
# connect to an X server and crashes on headless / broken-X-forwarding nodes
# (`_tkinter.TclError: couldn't connect to display ...`).
#
# Policy:
#   - Honour an explicit MPLBACKEND from the environment if the user set one.
#   - Otherwise, if DISPLAY is unset (true headless), force Agg.
#   - Leave DISPLAY-is-set cases untouched at import time. If the configured
#     interactive backend then fails at runtime (e.g. broken SSH X11
#     forwarding), the call sites fall back to Agg gracefully (see
#     janas_cmd_optimizer._safe_make_figure).
import os as _os
if "MPLBACKEND" not in _os.environ and not _os.environ.get("DISPLAY"):
    _os.environ["MPLBACKEND"] = "Agg"

from .starHandler import (
    header_columns,
    infoStarFile,
    dataOptics,
    merge_star_section,
    extract_particles_from_label_from_sections,
    read_star_sections,
    replace_star_columns_from_sections,
    delete_star_columns_from_sections,
    read_star_columns_from_sections,
    readColumns,
    readStar,
    removeColumns,
    removeColumnsTagsStartingWith,
    addColumns,
    writeDataframeToStar,
    extractBest,
    extractWorst,
    extractRandom,
    extractCategory,
    mergeRefinements,
)
from .starDisplay import resolutionPlot, plotEulerHist
from .assessParticles import ParticleVsReprojectionScores
from .utils import get_MRC_map_pixel_spacing

