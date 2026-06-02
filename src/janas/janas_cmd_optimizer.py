#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# File: janas_cmd_optimizer.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology


"""
Module: janas_optimizer.py

Provides a suite of routines for optimizing particle subset size based on
local-resolution evaluations, summarizing iterative runs, and managing
JANAS optimization workflows. Key capabilities include:

- Spline‐based prediction of the minimal number of particles needed
  to reach a target local resolution.
- Generation of summary figures across multiple iterations.
- Analysis of log and overview files in TOML and CSV formats.
- Automatic sampling of new particle subset sizes guided by previous
  local‐resolution results.
- Command‐line wrappers for each of the above functions.
"""

# Standard library
import argparse
import csv
import os.path
import shutil
import time
import re


# Third-party
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.interpolate import UnivariateSpline
import toml


def _make_figure_and_axes(showPlot: bool, **fig_kwargs):
    """
    Return ``(fig, ax)`` without touching the matplotlib GUI backend when
    the caller does not intend to display the figure interactively.

    The batch path (``showPlot=False``) returns a bare
    :class:`matplotlib.figure.Figure`, which is purely in-memory and never
    triggers backend initialisation — so a broken/absent display server
    (``_tkinter.TclError: couldn't connect to display ...``) cannot crash
    a save-to-PNG-only step.

    The interactive path (``showPlot=True``) goes through
    :func:`matplotlib.pyplot.subplots`, which produces a *managed* figure
    that :func:`matplotlib.pyplot.show` can later display via the
    configured interactive backend (any backend error there is honestly
    raised, since the user explicitly asked for a window).

    Both paths produce an object with the same ``ax`` interface, so the
    plotting code that follows is identical.
    """
    if showPlot:
        return plt.subplots(**fig_kwargs)
    from matplotlib.figure import Figure  # noqa: WPS433 — deliberate lazy import
    fig = Figure(**fig_kwargs)
    ax = fig.subplots()
    return fig, ax

# Local
import janas.janas_core as janas_core
from janas import starHandler
from janas.version import get_version

janas_parser = argparse.ArgumentParser(
    prog="janas_optimizer",
    usage="%(prog)s [command] [arguments]",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
janas_parser.add_argument(
    "-V", "--version",
    action="version",
    version=get_version(),
    help="show program’s version number and exit"
)
command = janas_parser.add_subparsers(dest="command")

#################################
#################################
## accessory functions
def _to_float(x, default=None):
    try:
        return float(x)
    except Exception:
        return default

def _to_int(x, default=None):
    try:
        return int(x)
    except Exception:
        # Also handle strings that include commas (e.g., "194,527")
        try:
            return int(str(x).replace(",", ""))
        except Exception:
            return default

def _to_bool(x, default=False):
    if isinstance(x, bool):
        return x
    s = str(x).strip().lower()
    if s in ("true", "1", "yes", "y"):
        return True
    if s in ("false", "0", "no", "n"):
        return False
    return default

#################################
#################################
## predict_min_particles
def predict_min_particles(
    file_path="",
    outputImageFile="",
    outputSplineFile="",
    showPlot=True,
    predicted_particles=None,
    ax=None,
    resolutionBestTarget_string="mean",
):
    sns.set_style("whitegrid")
    if ax is None:
        fig, ax = _make_figure_and_axes(showPlot=showPlot)
    else:
        fig = ax.get_figure()

    colors = sns.color_palette("colorblind", 10)
    data_color = colors[0]
    outlier_color = colors[2]
    spline_color = colors[3]

    if not file_path:
        return None

    df = pd.read_csv(file_path)

    if "numParticles" not in df.columns:
        raise KeyError(f"CSV is missing required column 'numParticles': {file_path}")

    # Decide which y-column to use, with robust mean<->median fallback
    requested = str(resolutionBestTarget_string or "mean").strip()
    if requested not in df.columns:
        if requested == "mean" and "median" in df.columns:
            requested = "median"
        elif requested == "median" and "mean" in df.columns:
            requested = "mean"
        else:
            raise KeyError(
                f"Requested column '{resolutionBestTarget_string}' not in CSV. "
                f"Available columns: {list(df.columns)}"
            )

    # x/y as numeric, drop non-finite rows
    x = pd.to_numeric(df["numParticles"], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df[requested], errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x = x[keep]
    y = y[keep]

    if x.size == 0:
        print("No valid numeric rows for spline fitting.")
        return None

    # Outlier detection (your original approach)
    Q1 = np.percentile(y, 10)
    Q3 = np.percentile(y, 70)
    IQR = Q3 - Q1
    outlier_threshold_upper = Q3 + 1.5 * IQR
    outlier_threshold_lower = Q1 - 1.5 * IQR

    fit_mask = (y >= outlier_threshold_lower) & (y <= outlier_threshold_upper)
    outlier_mask = ~fit_mask

    x_fit = x[fit_mask]
    y_fit = y[fit_mask]
    x_outlier = x[outlier_mask]
    y_outlier = y[outlier_mask]

    if len(x_fit) <= 1:
        print("Insufficient unique data points for spline fitting.")
        return None

    # Fit spline (weights favour lower y-values)
    k = max(1, min(4, len(x_fit) - 1))
    weights = 1.0 / np.maximum(y_fit, 1e-12)
    s = UnivariateSpline(x_fit, y_fit, w=weights, k=k, s=len(x_fit))
    x_smooth = np.linspace(float(np.min(x_fit)), float(np.max(x_fit)), 1000)
    y_smooth = s(x_smooth)

    # Whiskers: pick a centre column that exists (and matches requested when requested is mean/median)
    centre_col = "mean" if "mean" in df.columns else ("median" if "median" in df.columns else requested)
    if requested in ("mean", "median") and requested in df.columns:
        centre_col = requested

    # Whisker columns: require these to draw error bars.
    whisker_cols = ["max", "highQuartile", centre_col, "lowQuartile", "min"]
    missing = [c for c in whisker_cols if c not in df.columns]
    if missing:
        # If the only missing is centre_col, try the other central column
        if centre_col in missing and centre_col in ("mean", "median"):
            alt = "median" if centre_col == "mean" else "mean"
            if alt in df.columns:
                centre_col = alt
                whisker_cols = ["max", "highQuartile", centre_col, "lowQuartile", "min"]
                missing = [c for c in whisker_cols if c not in df.columns]

        if missing:
            raise KeyError(
                f"CSV is missing required whisker columns {missing}. "
                f"Available columns: {list(df.columns)}"
            )

    # Build whisker_data aligned with the *original* df rows, then apply keep mask
    whisker_data_all = df[whisker_cols].apply(pd.to_numeric, errors="coerce").to_numpy().T
    whisker_data_all = whisker_data_all[:, keep]  # now aligned with x/y

    # IMPORTANT: in your locresStats CSV, the header order is:
    # numParticles, max, highQuartile, centre, lowQuartile, min
    # but numerically "max" is typically the *best/lowest Å* and "min" the *worst/highest Å* (because of legacy naming).
    # So the lower errorbar should go from centre down to "max", and upper from centre up to "min".
    low_whisker_fit = whisker_data_all[0, fit_mask]   # "max" column
    high_whisker_fit = whisker_data_all[4, fit_mask]  # "min" column

    lower_err_fit = np.maximum(y_fit - low_whisker_fit, 0.0)
    upper_err_fit = np.maximum(high_whisker_fit - y_fit, 0.0)
    yerr_fit = np.vstack((lower_err_fit, upper_err_fit))

    ax.errorbar(
        x_fit,
        y_fit,
        yerr=yerr_fit,
        fmt="o",
        color=data_color,
        label="Local Resolution Variation",
        capsize=5,
        capthick=1.5,
    )

    # Centre markers (mean/median)
    centre_label = "Median" if centre_col == "median" else ("Mean" if centre_col == "mean" else centre_col)
    centre_vals_fit = whisker_data_all[2, fit_mask]
    for i, (xi, w_centre) in enumerate(zip(x_fit, centre_vals_fit)):
        ax.plot(
            xi,
            w_centre,
            "o",
            color=outlier_color,
            label=centre_label if i == 0 else "",
        )

    # Outliers (optional)
    if len(x_outlier) > 0:
        low_whisker_out = whisker_data_all[0, outlier_mask]
        high_whisker_out = whisker_data_all[4, outlier_mask]
        lower_err_out = np.maximum(y_outlier - low_whisker_out, 0.0)
        upper_err_out = np.maximum(high_whisker_out - y_outlier, 0.0)
        yerr_out = np.vstack((lower_err_out, upper_err_out))

        ax.errorbar(
            x_outlier,
            y_outlier,
            yerr=yerr_out,
            fmt="o",
            color=outlier_color,
            label="Outliers",
            capsize=5,
            capthick=1.5,
        )

    ax.plot(x_smooth, y_smooth, label="Fitted Spline", color=spline_color, linewidth=2.5)

    idx = int(np.argmin(y_smooth))
    predicted_particle_number = float(x_smooth[idx])

    # Best tested point (after outlier filtering) in the requested statistic
    best_observed_particle_number = float(x_fit[np.argmin(y_fit)])

    # If the spline minimum is effectively at the boundary, fall back to the best observed point.
    # This preserves the intended smoothing while preventing monotonic/oversmoothed splines from
    # incorrectly selecting the largest (or smallest) tested subset.
    x_min = float(np.min(x_fit))
    x_max = float(np.max(x_fit))
    boundary_tol = 0.01 * (x_max - x_min)  # 1% of the explored range

    if (
        predicted_particle_number <= x_min + boundary_tol
        or predicted_particle_number >= x_max - boundary_tol
    ):
        predicted_particle_number = best_observed_particle_number

    ax.axvline(
        x=predicted_particle_number,
        color=colors[1],
        linestyle="--",
        linewidth=1.5,
        label="Estimation",
    )

    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks([predicted_particle_number])
    ax2.set_xticklabels([f"{int(predicted_particle_number)}"])

    if predicted_particles:
        for particle in predicted_particles:
            ax.axvline(x=particle, color=colors[4], linestyle=":", linewidth=0.8)

        ax.set_xlim(
            min(float(np.min(x)), float(np.min(predicted_particles))),
            max(float(np.max(x)), float(np.max(predicted_particles))),
        )

    leg = ax.legend(loc="best", frameon=True)
    leg.get_frame().set_alpha(1.0)

    ax.set_title(f"{requested} Local Resolution per Number of Particles")
    ax.set_xlabel("Number of Particles")
    ax.set_ylabel(f"{requested} Local Resolution Estimation Value")

    if outputImageFile:
        directory = os.path.dirname(outputImageFile)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)
        # Use fig.savefig (not plt.savefig) so the batch path works on a
        # bare Figure that never registered with pyplot.
        fig.savefig(outputImageFile, dpi=300, format="png", bbox_inches="tight")

    if showPlot:
        plt.show()
        plt.close(fig)
    # When not showing, the bare Figure is discarded on function exit; no
    # explicit plt.close() is needed because the figure was never managed.

    if outputSplineFile:
        spline_data = pd.DataFrame({"numParticles": x_smooth, "estimatedResolution": y_smooth})
        output_directory = os.path.dirname(outputSplineFile)
        if output_directory and not os.path.exists(output_directory):
            os.makedirs(output_directory)
        spline_data.to_csv(outputSplineFile, index=False)

    return int(predicted_particle_number)

def _resolve_metric_column(df: pd.DataFrame, requested: str) -> str:
    """
    Return a column name that exists in df and matches the requested metric.
    Handles 'mean'/'median' mismatch and common variants.
    """
    req = (requested or "").strip().lower()

    # Common normalisations (extend if your CSV uses other names)
    alias = {
        "meanresolution": "mean",
        "medianresolution": "median",
        "avg": "mean",
        "average": "mean",
    }
    req = alias.get(req, req)

    cols = {c.lower(): c for c in df.columns}  # map lowercase->original

    # 1) Exact match
    if req in cols:
        return cols[req]

    # 2) Specific fallback between mean and median (your failure mode)
    if req == "mean" and "median" in cols:
        return cols["median"]
    if req == "median" and "mean" in cols:
        return cols["mean"]

    # 3) Known JANAS columns as last resort (best effort)
    for candidate in ("mean", "median", "min", "max", "lowquartile", "highquartile"):
        if candidate in cols:
            return cols[candidate]

    raise KeyError(
        f"None of the expected metric columns are present. Requested='{requested}', columns={list(df.columns)}"
    )


def create_summary_figure(basename, rows=3, columns=3):
    import matplotlib.gridspec as gridspec

    print("basename=", os.path.basename(basename))

    def extract_number(dir_name, prefix=os.path.basename(basename)):
        # Get the part of the string after the known prefix
        if prefix in dir_name:
            part_after_prefix = dir_name.split(prefix)[
                -1
            ]  # This gets the part of the string after the prefix
            # Now, extract consecutive digits from the beginning of this substring
            number_str = ""
            for char in part_after_prefix:
                if char.isdigit():
                    number_str += char
                else:
                    break  # Break at the first non-digit character

            if number_str:
                return int(number_str)

        return 0  # Return 0 if no number is found or handle it differently

    base_folder = "."
    if os.path.sep in basename:
        base_folder = os.path.dirname(basename)
    all_directories = [
        d
        for d in os.listdir(base_folder)
        if os.path.isdir(os.path.join(base_folder, d)) and "_janas_SCI__" in d
    ]

    # Sort directories based on the embedded number
    all_directories.sort(key=extract_number)

    # Break the selected directories into chunks for multiple figures
    figure_count = 0
    for i in range(0, len(all_directories), rows * columns):
        figure_count += 1

        selected_dirs = all_directories[i : i + rows * columns]

        # Create a figure and GridSpec object
        fig = plt.figure(figsize=(16, 5 * rows))
        gs = gridspec.GridSpec(rows, columns)  # Grid based on input rows and columns

        for j, dir_name in enumerate(selected_dirs):
            file_path = os.path.join(
                base_folder, dir_name, "bestRanked_locres_values.csv"
            )

            ax = fig.add_subplot(gs[j])

            # Call your plotting function using the given axes
            predict_min_particles(file_path=file_path, showPlot=False, ax=ax)
            iteration_number = extract_number(dir_name)
            ax.set_title(
                f"Iteration {iteration_number}"
            )  # Use the extracted number as title for each subplot

        plt.tight_layout()  # Adjust layout

        # Save the figure with an incremented name for multiple figures
        fig.savefig(
            os.path.join(base_folder, f"summary_{figure_count}.png"),
            dpi=300,
            format="png",
            bbox_inches="tight",
        )
        fig.savefig(
            os.path.join(base_folder, f"summary_{figure_count}.pdf"),
            format="pdf",
            bbox_inches="tight",
        )
        plt.show()


#################################
## log analyzer
janas_logAnalyzer = command.add_parser(
    "logAnalyzer", description="compute logAnalyzer", help="logAnalyzer"
)
janas_logAnalyzer.add_argument(
    "--prefix",
    required=True,
    type=str,
    help="folders prefix, e.g. _janas_SCI__1.00_scored_selection_",
)
janas_logAnalyzer.add_argument(
    "--numRows",
    required=False,
    type=int,
    default=3,
    help="now of rows in reported figure (default=3)",
)
janas_logAnalyzer.add_argument(
    "--numColumns",
    required=False,
    type=int,
    default=3,
    help="now of columns in reported figure (default=3)",
)


def logAnalyzer(args):
    create_summary_figure(args.prefix, args.numRows, args.numColumns)


#################################
## Analyse reconstruction script
janas_getNumParticles = command.add_parser(
    "getNumParticles",
    description="compute the optimal Target Number Of Particles, it requires the locres summary file",
    help="compute the optimal Target Number Of Particles, it requires the locres summary file",
)
janas_getNumParticles.add_argument(
    "--locres", required=True, type=str, help="file with locres evaluation"
)
janas_getNumParticles.add_argument(
    "--plot", action="store_true", help="Display the plot"
)
janas_getNumParticles.add_argument(
    "--mean_res", action="store_true", help="predict the best mean local resolution"
)
janas_getNumParticles.add_argument(
    "--plotOnFile", required=False, default="", type=str, help="Save the plot on file"
)
janas_getNumParticles.add_argument(
    "--saveSplineOnCsv",
    required=False,
    default="",
    type=str,
    help="Save the plot on file",
)
janas_getNumParticles.add_argument(
    "--save",
    required=False,
    default="",
    type=str,
    help="Save the best reconstruction particle on file",
)
janas_getNumParticles.add_argument(
    "--median_res",
    action="store_true",
    help="Use the 'median' column from locresStats CSV (instead of 'mean').",
)
janas_getNumParticles.add_argument(
    "--resolutionBestTarget",
    required=False,
    type=lambda v: v.lower()
    if v.lower()
    in [
        "meanresolution",
        "highresolution",
        "highresolutionquartile",
        "lowresolution",
        "lowresolutionquartile",
    ]
    else argparse.ArgumentTypeError(
        f"Invalid choice: {v}. Choose from ['meanResolution', 'highResolution', 'highresolutionquartile', 'lowResolution', 'lowresolutionquartile']"
    ),
    default="meanResolution",
    help="Choose the best resolution target from: [meanResolution, highResolution, highresolutionquartile, lowResolution, lowresolutionquartile]. Case insensitive.",
)

def getNumParticles(args):
    if not os.path.isfile(args.locres):
        print('ERROR: file "', args.locres, '" not existing')
        exit()

    # Default mapping from resolutionBestTarget to CSV column
    resolutionBestTarget_string = {
        "meanresolution": "mean",
        "highresolution": "max",
        "highresolutionquartile": "highQuartile",
        "lowresolution": "min",
        "lowresolutionquartile": "lowQuartile",
    }.get((args.resolutionBestTarget or "meanresolution").lower(), "mean")

    # Explicit overrides via flags (these are clearer than resolutionBestTarget)
    if getattr(args, "median_res", False):
        resolutionBestTarget_string = "median"
    elif getattr(args, "mean_res", False):
        resolutionBestTarget_string = "mean"

    # If the requested column is missing, auto-fallback mean<->median
    try:
        cols = pd.read_csv(args.locres, nrows=1).columns
        if resolutionBestTarget_string not in cols:
            if resolutionBestTarget_string == "mean" and "median" in cols:
                resolutionBestTarget_string = "median"
            elif resolutionBestTarget_string == "median" and "mean" in cols:
                resolutionBestTarget_string = "mean"
    except Exception:
        pass

    # Use the predict_min_particles function with the chosen column
    result = predict_min_particles(
        args.locres,
        outputImageFile=args.plotOnFile,
        outputSplineFile=args.saveSplineOnCsv,
        showPlot=args.plot,
        resolutionBestTarget_string=resolutionBestTarget_string,
    )

    # Override with best-observed tested subset ONLY when the user
    # explicitly requested it via --mean_res or --median_res.
    #
    # Previously this branch also checked `args.resolutionBestTarget`,
    # which always has a default value ("meanResolution") and is therefore
    # always truthy.  That caused the spline result to be unconditionally
    # replaced with the raw idxmin() lookup, making the spline prediction
    # unreachable.  The --resolutionBestTarget flag now controls which
    # metric column is used for BOTH the spline fitting and the fallback,
    # but does not by itself force the best-observed override.
    if args.mean_res or getattr(args, "median_res", False):
        data = pd.read_csv(args.locres)
        if resolutionBestTarget_string not in data.columns:
            print(
                f"ERROR: column '{resolutionBestTarget_string}' not found in {args.locres}. "
                f"Available columns: {list(data.columns)}"
            )
            exit(1)

        idx_min_value = data[resolutionBestTarget_string].idxmin()
        result = data.loc[idx_min_value, "numParticles"]

    print(result)

    if args.save:
        directory = os.path.dirname(args.save)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)
        with open(args.save, "w") as f:
            f.write(str(int(result)))


#################################
## Analyse reconstruction script


janas_plotOverview = command.add_parser(
    "plotOverview", description="plot the overview", help="plot the overview"
)
janas_plotOverview.add_argument(
    "--overview", required=True, type=str, help="overviewFile"
)
janas_plotOverview.add_argument(
    "--plot", action="store_true", help="Display the figure on screen"
)
janas_plotOverview.add_argument(
    "--o", required=False, default="", type=str, help="output image prefix (PNG + PDF)"
)

def _parse_locres_stats(stats_str):
    """
    Parse reference_locres_stats from overview.
    Expected order (per your note):
      num_particles, best, best_quartile, mean, worst_quartile, worst
    Returns dict with keys:
      np, best, q1, mean, q3, worst
    """
    if stats_str is None:
        return None
    parts = [p.strip() for p in str(stats_str).split(",")]
    if len(parts) < 6:
        return None
    npart = _to_int(parts[0])
    best  = _to_float(parts[1])
    q1    = _to_float(parts[2])
    mean  = _to_float(parts[3])
    q3    = _to_float(parts[4])
    worst = _to_float(parts[5])
    if None in (npart, best, q1, mean, q3, worst):
        return None
    return {"np": npart, "best": best, "q1": q1, "mean": mean, "q3": q3, "worst": worst}

def _load_overview_tables(overview_path):
    """
    Return:
      target: dict or None (from [[_janas_target_selection]])
      selections: list of (idx:int, row:dict) for each [[_janas_selection_<idx>]] block
    """
    data = toml.load(overview_path)

    target_list = data.get("_janas_target_selection", [])
    target = target_list[0] if isinstance(target_list, list) and target_list else None

    selections = []
    for key, val in data.items():
        if not key.startswith("_janas_selection_"):
            continue
        m = re.match(r"_janas_selection_(\d+)$", key)
        if not m:
            continue
        idx = int(m.group(1))
        if isinstance(val, list) and val:
            row = val[0]
        elif isinstance(val, dict):
            row = val
        else:
            continue
        selections.append((idx, row))
    selections.sort(key=lambda x: x[0])
    return target, selections


def print_and_save_overview_summary(overview_path):
    """
    Prints a compact summary:
      - Start (_janas_selection_0): star file, num particles, target metric
      - Final target ([_janas_target_selection]): star file, num particles,
        last_consecutive_non_improving_selections, percentage_particles_retained,
        target metric and delta vs start.
    """
    if not os.path.isfile(overview_path):
        print(f'ERROR: overview file "{overview_path}" not existing')
        return None

    target, selections = _load_overview_tables(overview_path)
    if not selections:
        print("No [[_janas_selection_*]] sections found.")
        return None

    # Start = selection_0
    start_row = dict(selections[0][1])  # assumes index 0 exists
    start_star = start_row.get("reference_starFile")
    start_np   = _to_int(start_row.get("reference_num_particles"))
    start_rt   = _to_float(start_row.get("reference_locres_ResolutionTarget"))

    # Target
    if not target:
        print("Missing [[_janas_target_selection]] in overview.")
        return None
    tgt_star = target.get("reference_starFile")
    tgt_np   = _to_int(target.get("reference_num_particles"))
    tgt_rt   = _to_float(target.get("reference_locres_ResolutionTarget"))
    tgt_last_unimpr = _to_int(target.get("last_consecutive_non_improving_selections"), 0)
    tgt_pct_retained = _to_float(target.get("percentage_particles_retained"))

    # Δ (negative means improvement)
    delta_rt = None
    if start_rt is not None and tgt_rt is not None:
        delta_rt = tgt_rt - start_rt

    # Print summary for copy/paste in terminal
    print("\n=== JANAS Overview Summary ===\n")
    print("Start (selection_0):")
    print(f"  reference_starFile: {start_star}")
    print(f"  reference_num_particles: {start_np}")
    print(f"  reference_locres_ResolutionTarget: {start_rt}")

    print("\nFinal target ([_janas_target_selection]):")
    print(f"  reference_starFile: {tgt_star}")
    print(f"  reference_num_particles: {tgt_np}")
    print(f"  last_consecutive_non_improving_selections: {tgt_last_unimpr}")
    print(f"  percentage_particles_retained: {tgt_pct_retained}")
    print(f"  reference_locres_ResolutionTarget: {tgt_rt}")
    if delta_rt is not None:
        sign = "+" if delta_rt >= 0 else ""
        print(f"  Δ ResolutionTarget vs start: {sign}{delta_rt:.5f}")

    # Return a dict to reuse in the side table
    return {
        "start_star": start_star,
        "start_np": start_np,
        "start_rt": start_rt,
        "tgt_star": tgt_star,
        "tgt_np": tgt_np,
        "tgt_rt": tgt_rt,
        "tgt_last_unimpr": tgt_last_unimpr,
        "tgt_pct_retained": tgt_pct_retained,
        "delta_rt": delta_rt,
    }

def plot_standard_whiskers_with_side_table(overview_path, out_image_path="", show=False):
    """
    Left: standard box/whisker plot per selection (Q1–Q3 box, whiskers at best/worst).
          Mean is drawn as the horizontal line inside the box (substituting 'median').
          A dotted horizontal line marks the best mean achieved across selections.
    Right: summary table (start vs final target).
    Also prints the same summary to stdout.
    """
    if not os.path.isfile(overview_path):
        print(f'ERROR: overview file "{overview_path}" not existing')
        return

    target, selections = _load_overview_tables(overview_path)
    if not selections:
        print("No [[_janas_selection_*]] sections found.")
        return

    # Collect stats
    rows = []
    for idx, row in selections:
        stats = _parse_locres_stats(row.get("reference_locres_stats"))
        if stats is None:
            continue
        rows.append({
            "iteration": idx,
            "best": stats["best"],
            "q1": stats["q1"],
            "mean": stats["mean"],
            "q3": stats["q3"],
            "worst": stats["worst"],
            "np": stats["np"],
        })

    if not rows:
        print("No usable reference_locres_stats found.")
        return

    df = pd.DataFrame(rows).sort_values("iteration")

    # Prepare summary
    summary = print_and_save_overview_summary(overview_path)
    if summary is None:
        return

    # --- Figure ---
    # Build the figure without touching the matplotlib GUI backend when the
    # caller will not display it (typical: --o set, --plot not set).
    # When showing is requested, go through pyplot so plt.show() can find
    # the managed figure.
    import matplotlib.gridspec as gridspec
    _do_show = bool(show) or not bool(out_image_path)
    if _do_show:
        fig = plt.figure(figsize=(16, 6), constrained_layout=True)
    else:
        from matplotlib.figure import Figure  # noqa: WPS433
        fig = Figure(figsize=(16, 6), constrained_layout=True)
    gs = gridspec.GridSpec(1, 2, width_ratios=[2.2, 1], figure=fig)

    # Left panel: box/whisker
    ax = fig.add_subplot(gs[0])

    # Build data for ax.bxp
    box_data = []
    for _, r in df.iterrows():
        box_data.append({
            'label': str(int(r['iteration'])),
            'whislo': r['best'],    # min
            'q1':    r['q1'],
            'med':   r['mean'],     # put mean where median normally goes
            'q3':    r['q3'],
            'whishi':r['worst'],    # max
            'fliers': []
        })
    ax.bxp(box_data, showfliers=False)

    # Horizontal line for best mean achieved
    best_mean = df['mean'].min()
    ax.axhline(best_mean, linestyle=':', color='red', linewidth=1.5,
               label=f'Best mean = {best_mean:.2f} Å')
    ax.legend(loc='best')

    ax.set_xlabel("Selection iteration")
    ax.set_ylabel("Local resolution (Å)")
    ax.set_title("Local-resolution box/whisker per selection")

    # Right: summary table text
    ax2 = fig.add_subplot(gs[1])
    ax2.axis('off')

    lines = []
    # Bold headings using fontweight
    lines.append(("JANAS Overview Summary", "bold"))
    lines.append(("", None))

    lines.append(("Starting dataset", "bold"))
    lines.append((f"  reference_starFile:\n    {summary['start_star']}", None))
    lines.append(("", None))    
    lines.append((f"  reference_num_particles: {summary['start_np']}", None))
    lines.append((f"  reference_locres_ResolutionTarget: {summary['start_rt']}", None))
    lines.append(("", None))

    lines.append(("Selected subset of the data:", "bold"))
    lines.append((f"  reference_starFile:\n    {summary['tgt_star']}", None))
    lines.append(("", None))
    lines.append((f"  reference_num_particles: {summary['tgt_np']}", None))
    lines.append((f"  percentage_particles_retained: {summary['tgt_pct_retained']}", None))
    lines.append((f"  reference_locres_ResolutionTarget: {summary['tgt_rt']}", None))
    if summary['delta_rt'] is not None:
        sign = "+" if summary['delta_rt'] >= 0 else ""
        lines.append((f"  Δ ResolutionTarget vs start: {sign}{summary['delta_rt']:.5f}", None))

    # Render line by line, with bold where needed
    y = 1.0
    dy = 0.06
    for text, weight in lines:
        if text == "" and weight is None:
            y -= dy / 2
            continue
        ax2.text(0, y, text, va='top', ha='left', fontsize=10,
                 family='monospace',
                 fontweight=weight if weight else 'normal')
        y -= dy

    # Save
    if out_image_path:
        png_path = out_image_path if out_image_path.lower().endswith(".png") else out_image_path + ".png"
        fig.savefig(png_path, dpi=300, bbox_inches="tight")
        pdf_path = png_path[:-4] + ".pdf"
        fig.savefig(pdf_path, bbox_inches="tight")
        print(f"Saved: {png_path}")
        print(f"Saved: {pdf_path}")

    if _do_show:
        plt.show()
        plt.close(fig)
    # else: the figure is a bare Figure() and goes out of scope cleanly



def plotOverview(args):
    """
    Produce a single figure with a whisker plot (left) and a summary table (right),
    and print the summary to stdout. Use --o to save (PNG+PDF), --plot to display.
    """
    overview_path = args.overview
    plot_standard_whiskers_with_side_table(
        overview_path,
        out_image_path=args.o,
        show=args.plot or not bool(args.o)
    )
    



#################################
#################################

from scipy.stats import skewnorm as _skewnorm


def skewed_gaussian(mean, skew_factor, size, seed=0, sampling_density_factor=0.15):
    """
    Sample from an Azzalini skew-normal distribution.

    Parameters
    ----------
    mean : float
        Location parameter (centre of the distribution).
    skew_factor : float
        Shape parameter α of the Azzalini skew-normal.  Negative values
        produce a left-skewed (negatively skewed) distribution, biasing
        samples toward values smaller than *mean*.
    size : int
        Number of samples to draw.
    seed : int
        Random seed for reproducibility (0 = fixed seed).
    sampling_density_factor : float
        The scale parameter is ``mean * sampling_density_factor``.  Larger
        values produce a wider distribution.

    Returns
    -------
    np.ndarray of shape (size,)
    """
    scale = mean * sampling_density_factor
    rng = np.random.default_rng(seed)
    values = _skewnorm.rvs(
        a=skew_factor, loc=mean, scale=scale, size=size, random_state=rng,
    )
    return values


def generate_particle_samples(
    predicted_particle_number,
    n_samples,
    min_val,
    max_val,
    sampling_density_factor=0.25,
    extraSamples_num=1,
    extraSamples_randomness=0.1,
    extraSamples_audacity=0.5,   # DEPRECATED — retained for one release cycle
    seed=0,
    force_include=None,
    include_max=True,
    hard_max=None,
    skew_alpha=-4.0,
):
    """
    Generate a set of candidate particle-subset sizes.

    The candidates are drawn from a negatively skewed Azzalini skew-normal
    distribution centred on *predicted_particle_number*, clipped to
    [min_val, max_val], and thinned **from the centre outward** so that
    the neighbourhood of the current optimum is always well-sampled.

    Parameters
    ----------
    predicted_particle_number : int or float
        Centre (location) of the sampling distribution — the current
        best estimate of the optimum (spline-predicted or best observed).
    n_samples : int
        Desired number of candidate subset sizes.
    min_val, max_val : float
        Feasible particle-count interval.
    sampling_density_factor : float
        Scale of the skew-normal as a fraction of the centre.
    extraSamples_num : int
        Number of exploratory samples placed well below the centre.
    extraSamples_randomness : float
        Width of the random jitter applied to the extra samples.
    extraSamples_audacity : float
        .. deprecated::
            This parameter has no effect.  It is retained in the
            function signature and CLI for one release cycle to avoid
            breaking existing scripts.  It will be removed in a
            future version.
    seed : int
        Random seed.
    force_include : list[int] or None
        Values that *must* appear in the output (e.g. best previous
        subset).  **Forced points bypass the minimum-spacing rule.**
        They are merged into the candidate list after thinning and
        are never removed by the spacing filter.  This guarantees
        that the best previous subset is always present, even if it
        would otherwise violate the nominal gap.
    include_max : bool
        If True, *max_val* is appended to the output (also bypasses
        spacing).
    hard_max : int or None
        If set, all candidates above this value are dropped after
        merging.
    skew_alpha : float
        Shape parameter for the Azzalini skew-normal (negative = left skew).

    Returns
    -------
    list[int]  – sorted candidate particle counts.

    Notes
    -----
    **Spacing semantics.**  The minimum-gap rule is applied only
    during the centre-outward thinning of the oversampled skew-normal
    draw (step 3).  Forced points (``force_include``, exploratory
    extra samples, optional ``max_val``) are added *after* thinning
    and are unconditionally retained.  This means the final output
    may contain pairs of adjacent values closer than ``min_gap`` when
    one of them is a forced point.  This is by design: the previous
    best subset and the upper-bound anchor must never be discarded
    by the spacing filter.

    **Output size.**  The returned list is not guaranteed to contain
    exactly ``n_samples`` elements.  Thinning stops at ``n_samples``
    accepted candidates, but forced points are merged afterwards
    (potentially increasing the count), and ``hard_max`` filtering
    may remove entries (potentially decreasing it).  Callers should
    not assume a fixed output length.
    """

    # ---- Deprecation notice for extraSamples_audacity ----
    # FutureWarning (not DeprecationWarning) is used because Python's
    # default filters show FutureWarning to end users, whereas
    # DeprecationWarning is suppressed unless the caller opts in.
    import warnings
    if extraSamples_audacity != 0.5:
        warnings.warn(
            "extraSamples_audacity has no effect and will be removed in a "
            "future release.  You can safely drop this parameter.",
            FutureWarning,
            stacklevel=2,
        )

    centre = float(predicted_particle_number)
    global_max_val = float(max_val)

    # ---- 1. Draw over-sampled raw candidates from skew-normal ----
    n_oversample = max(n_samples * 10, 500)
    raw = skewed_gaussian(
        centre,
        skew_factor=skew_alpha,
        size=n_oversample,
        seed=seed,
        sampling_density_factor=sampling_density_factor,
    )
    raw = np.clip(raw, max(0, min_val), global_max_val)

    # ---- 2. Collect forced / special points ----
    forced_set = set()
    if force_include:
        for v in force_include:
            if v is not None:
                forced_set.add(int(v))

    # Exploratory "extra" samples well below centre
    rng = np.random.default_rng(seed)
    for i in range(extraSamples_num):
        lower_bound = centre / (2 + i + 1)
        upper_bound = centre / (2 + i - 1) if (2 + i - 1) > 0 else centre
        mid = (upper_bound + lower_bound) / 2.0
        half_w = (upper_bound - lower_bound) * extraSamples_randomness / 2.0
        extra = rng.uniform(mid - half_w, mid + half_w)
        forced_set.add(int(np.clip(extra, min_val, global_max_val)))

    if include_max:
        forced_set.add(int(global_max_val))

    # ---- 3. Centre-outward greedy thinning ----
    # Minimum gap: prevent redundant evaluations while allowing density
    # near the centre.  We use a fraction of the feasible range scaled
    # by the requested number of samples.
    feasible_range = global_max_val - min_val
    min_gap = feasible_range / max(n_samples * 2, 4)

    # Unique integer candidates, sorted by proximity to the centre
    candidates = sorted(set(int(round(v)) for v in raw), key=lambda v: abs(v - centre))

    accepted = []
    for c in candidates:
        if len(accepted) >= n_samples:
            break
        if c < min_val or c > global_max_val:
            continue
        if all(abs(c - a) >= min_gap for a in accepted):
            accepted.append(c)

    # ---- 4. Merge forced points, enforce hard_max, sort ----
    merged = set(accepted) | forced_set
    merged = sorted(merged)

    if hard_max is not None:
        try:
            hm = int(hard_max)
            merged = [v for v in merged if v <= hm]
        except Exception:
            pass

    return merged


def automaticParticleSubsetsCore(
    locresResultsCsvFile,
    maxNumberOfParticles,
    number_of_sampling,
    randomSeed=True,
    showPlot=False,
    outputImageFile="",
    sampling_density_factor=0.25,
    extraSamples_num=1,
    extraSamples_randomness=0.1,
    extraSamples_audacity=0.5,
    resolutionBestTarget_string="mean",
    aggressive=False,
):
    """
    Core routine for generating candidate particle-subset sizes.

    The sampling centre is the current best estimate of the optimum,
    determined by the following precedence:

      1. If ``predict_min_particles()`` returns a non-None value, that
         value is used as the centre.  Because ``predict_min_particles()``
         already applies a boundary-fallback internally (replacing
         boundary-clipped spline minima with the best observed point),
         a non-None return is always a usable estimate — either a
         genuine interior spline optimum or the best observed count.
      2. If ``predict_min_particles()`` returns None (spline fitting
         failed entirely), the centre is set to ``best_observed_n``.
      3. On the first iteration (no previous data), the centre defaults
         to ``maxNumberOfParticles``.

    When ``aggressive`` is True (the typical iterative path set by the
    session manager) the upper bound is capped at the best previous
    count and the maximum tested value is *not* retained, so exploration
    is focused below the current optimum.

    The ``aggressive`` flag controls only two things:

      1. Whether the candidate list is capped at *best_observed_n*
         (``hard_max``).
      2. Whether the largest tested subset is appended (``include_max``).

    It does **not** change the sampling centre.
    """

    if locresResultsCsvFile and locresResultsCsvFile != "":
        df = pd.read_csv(locresResultsCsvFile)

        # Column to optimise in the locres CSV
        target_col = _resolve_metric_column(df, resolutionBestTarget_string)

        # ---- Numeric coercion and NaN/non-finite filtering ----
        # Mirrors the robustness of predict_min_particles(): coerce to
        # numeric and drop rows where either column is non-finite, so
        # that malformed CSV values cannot propagate into the sampling
        # logic.
        df_clean = df.copy()
        df_clean["numParticles"] = pd.to_numeric(
            df_clean["numParticles"], errors="coerce"
        )
        df_clean[target_col] = pd.to_numeric(
            df_clean[target_col], errors="coerce"
        )
        finite_mask = (
            np.isfinite(df_clean["numParticles"].to_numpy(dtype=float))
            & np.isfinite(df_clean[target_col].to_numpy(dtype=float))
        )
        df_clean = df_clean[finite_mask]

        if df_clean.empty:
            # All rows were invalid — fall back to first-iteration defaults
            print(
                "WARNING: locres CSV has no valid numeric rows; "
                "falling back to first-iteration defaults."
            )
            predicted_particle_number = maxNumberOfParticles
            best_observed_n = predicted_particle_number
            min_val = maxNumberOfParticles / 3.0
            max_val = float(maxNumberOfParticles)
            predicted_opt_n = None
        else:
            best_observed_n = int(
                df_clean.loc[df_clean[target_col].idxmin(), "numParticles"]
            )

            # Spline-based prediction
            predicted_opt_n = predict_min_particles(
                locresResultsCsvFile,
                showPlot=showPlot,
                outputImageFile=outputImageFile,
                resolutionBestTarget_string=resolutionBestTarget_string,
            )
            print("predicted_particle_number = ", predicted_opt_n)

            # ---- Sampling centre: optimum-estimate precedence ----
            #
            # The methods text specifies that the location parameter is
            # the best previous particle count, and "when spline
            # interpolation is employed" it corresponds to the spline-
            # predicted optimum.  predict_min_particles() already falls
            # back to the best observed point when the spline minimum
            # is at a boundary, so a non-None return that differs from
            # best_observed_n indicates a genuine interior spline optimum.
            #
            # Precedence:
            #   1. valid spline prediction (interior)  → predicted_opt_n
            #   2. spline unavailable or at boundary   → best_observed_n
            #   3. no previous data                    → maxNumberOfParticles
            #
            if predicted_opt_n is not None:
                centre_n = int(predicted_opt_n)
            else:
                centre_n = best_observed_n

            # Sampling range
            min_val = float(df_clean["numParticles"].min()) * 2.0 / 3.0
            max_val = float(df_clean["numParticles"].max())

            predicted_particle_number = int(centre_n)

    else:
        # First iteration: no previous data
        predicted_particle_number = maxNumberOfParticles
        best_observed_n = predicted_particle_number
        min_val = maxNumberOfParticles * 1.0 / 3.0
        max_val = maxNumberOfParticles

    if randomSeed:
        seed = int(time.time())
    else:
        seed = 0

    predicted_particles = generate_particle_samples(
        predicted_particle_number,
        number_of_sampling,
        min_val,
        max_val,
        sampling_density_factor=sampling_density_factor,
        extraSamples_num=extraSamples_num,
        extraSamples_randomness=extraSamples_randomness,
        extraSamples_audacity=extraSamples_audacity,
        seed=seed,
        force_include=[best_observed_n],
        include_max=(not aggressive),
        hard_max=(best_observed_n if aggressive else None),
    )
    return predicted_particles


#################################
## automaticParticleSubsets
janas_automaticParticleSubsets = command.add_parser(
    "automaticParticleSubsets",
    description="compute automaticParticleSubsets",
    help="compute automaticParticleSubsets",
)
janas_automaticParticleSubsets.add_argument(
    "--starFile", required=True, type=str, help="file with the input star file"
)
janas_automaticParticleSubsets.add_argument(
    "--locres",
    required=False,
    type=str,
    default="",
    help="file with the previous locres evaluation file is",
)
janas_automaticParticleSubsets.add_argument(
    "--save",
    required=False,
    type=str,
    default="",
    help="file where to save the particle to reconstruct and compute locres, stored as comma separated file (csv)",
)
janas_automaticParticleSubsets.add_argument(
    "--plot", action="store_true", help="Display the plot"
)
janas_automaticParticleSubsets.add_argument(
    "--plotPrediction",
    action="store_true",
    help="Display the predictions for next plot",
)
janas_automaticParticleSubsets.add_argument(
    "--plotOnFile",
    required=False,
    default="",
    type=str,
    help="Save the plot on an image file",
)
janas_automaticParticleSubsets.add_argument(
    "--numSamples", required=False, type=int, default=10, help="number of samples"
)
janas_automaticParticleSubsets.add_argument(
    "--samplingDensityFactor",
    required=False,
    type=float,
    default=0.15,
    help="sampling density factor: default=0.15, higher values produce a wider distribution",
)
janas_automaticParticleSubsets.add_argument(
    "--extraSamples_num",
    required=False,
    type=int,
    default=1,
    help="number of extra samples for more correct optimization (to reduce the local minima trap)",
)
janas_automaticParticleSubsets.add_argument(
    "--extraSamples_randomness",
    required=False,
    type=float,
    default=0.1,
    help="number of extra samples randomness, between 0 and 1: default 0.1, max=1",
)
janas_automaticParticleSubsets.add_argument(
    "--extraSamples_audacity",
    required=False,
    type=float,
    default=0.5,
    help="(Deprecated, no effect.) Will be removed in a future release.",
)
janas_automaticParticleSubsets.add_argument(
    "--aggressive",
    required=False,
    action="store_true",
    help="If set, cap candidate subset sizes at the best previous particle count and exclude the upper bound, focusing exploration below the current optimum. Does not change the sampling centre.",
)
janas_automaticParticleSubsets.add_argument(
    "--resolutionBestTarget",
    required=False,
    type=lambda v: v.lower()
    if v.lower()
    in [
        "meanresolution",
        "highresolution",
        "highresolutionquartile",
        "lowresolution",
        "lowresolutionquartile",
    ]
    else argparse.ArgumentTypeError(
        f"Invalid choice: {v}. Choose from ['meanResolution', 'highResolution', 'highresolutionquartile', 'lowResolution',  'lowresolutionquartile']"
    ),
    default="meanResolution",
    help="Choose the best resolution target from: [meanResolution, highResolution, highresolutionquartile, lowResolution, lowresolutionquartile]. Case insensitive.",
)


# janas_automaticParticleSubsets.add_argument("--pureRandom", action="store_true", help="pure random number, not reproducible")
def automaticParticleSubsets(args):
    if not os.path.isfile(args.starFile):
        print('ERROR: file "', args.starFile, '" not existing')
        exit()

    # ---- Argument validation ----
    if args.numSamples < 1:
        print("ERROR: --numSamples must be a positive integer (got %d)" % args.numSamples)
        exit(1)
    if args.samplingDensityFactor <= 0:
        print("ERROR: --samplingDensityFactor must be positive (got %g)" % args.samplingDensityFactor)
        exit(1)
    if not (0.0 <= args.extraSamples_randomness <= 1.0):
        print("ERROR: --extraSamples_randomness must be between 0 and 1 (got %g)" % args.extraSamples_randomness)
        exit(1)
    if args.extraSamples_num < 0:
        print("ERROR: --extraSamples_num must be non-negative (got %d)" % args.extraSamples_num)
        exit(1)

    from janas import starHandler

    starFile = starHandler.readStar(args.starFile)

    num_non_null_items = int(starFile.count().iloc[0])

    resolutionBestTarget_string = {
        "meanresolution": "mean",
        "highresolution": "max",
        "highresolutionquartile": "highQuartile",
        "lowresolution": "min",
        "lowresolutionquartile": "lowQuartile",
    }.get(args.resolutionBestTarget.lower(), "meanresolution")
    if getattr(args, "median_res", False):
        resolutionBestTarget_string = "median"

    if os.path.isfile(args.locres):
        expectedEstimatedParticlesNumber = predict_min_particles(
            args.locres,
            resolutionBestTarget_string=resolutionBestTarget_string,
            showPlot=False,
        )
        print("expected number of particles=", expectedEstimatedParticlesNumber)
    else:
        args.locres = ""
        print(
            "WARNING: you might want to specify a valid locres file, however it is ok for the first iteration\n"
        )
        expectedEstimatedParticlesNumber = num_non_null_items

    # print ("plot on file=",args.plotOnFile)
    result = automaticParticleSubsetsCore(
        args.locres,
        num_non_null_items,
        args.numSamples,
        showPlot=args.plot,
        outputImageFile=args.plotOnFile,
        sampling_density_factor=args.samplingDensityFactor,
        extraSamples_num=args.extraSamples_num,
        extraSamples_randomness=args.extraSamples_randomness,
        extraSamples_audacity=args.extraSamples_audacity,
        resolutionBestTarget_string=resolutionBestTarget_string,
        aggressive=getattr(args, "aggressive", False),
    )

    #    predict_min_particles(args.locres, showPlot=args.plot, predicted_particles=result, outputImageFile=args.plotOnFile)
    #    print("Particles for performing selection",result)
    if not args.save == "":
        directory = os.path.dirname(args.save)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)
        with open(args.save, "w") as f:
            f.write(",".join(map(str, result)))


#####################################################
######## OVERVIEW FILE
#####################################################
def overview_read_item_at_iteration(overview_file_path, iteration_number):
    print("Load the TOML file")
    try:
        data = toml.load(overview_file_path)
    except Exception as e:
        print(f"Error reading TOML file: {e}")
        return None

    print("Iterate through each section to find the matching iteration")
    for section, value in data.items():
        # Only attempt membership if this is a mapping (i.e. table) rather than a primitive
        if isinstance(value, dict):
            if value.get("iteration") == iteration_number:
                # Return the target_num_particles field, or None if it does not exist
                return value.get("target_num_particles")
        # If value is not a dict, skip it

    print(f"Iteration {iteration_number} not found in the file.")
    return None



def remove_duplicates_toml(file_path):
    # Read the TOML file as plain text
    with open(file_path, "r") as file:
        lines = file.readlines()

    # Process each line to remove duplicates, keeping only the last occurrence
    new_lines = []
    seen_keys = set()
    for line in reversed(lines):
        if "=" in line:
            key = line.split("=", 1)[0].strip()
            if key not in seen_keys:
                seen_keys.add(key)
                new_lines.append(line)
        else:
            new_lines.append(line)
            seen_keys.clear()  # Reset for a new table

    # Reverse the lines back to original order
    new_lines = new_lines[::-1]

    # Write the processed lines back to the file
    with open(file_path, "w") as file:
        file.writelines(new_lines)


janas_generate_overview = command.add_parser(
    "generate_overview", description="generate overview", help="generate overview"
)
janas_generate_overview.add_argument(
    "--directory",
    required=True,
    type=str,
    help="directory with the overview to be generated",
)
janas_generate_overview.add_argument(
    "--verbose", action="store_true", help="verbose overview display"
)


def extractSigmaFromTag(tag):
    """
    Given a directory name (tag) of the form "_janas_<score>__<sigma>_…" (or legacy "_emprove_"),
    return the "<sigma>" portion as a string. If no pattern is found,
    return "0".
    """
    match = re.search(r'(_janas_|_emprove_)[^_]+__([0-9.]+)_', tag)
    return match.group(2) if match else "0"

def generate_overview(args):
    prefix = "_janas_"

    def _lookup_fullmask_stats(fullmask_csv_path, num_particles, resolution_key):
        """Look up the row matching *num_particles* in the full-mask CSV.

        Returns (stats_string, resolution_value) or (None, None) on failure.
        Matching is by particle count, not row order, because the adaptive
        and full-mask CSVs may not share row ordering.
        """
        try:
            if not os.path.isfile(fullmask_csv_path):
                return None, None
            with open(fullmask_csv_path, mode="r") as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    if str(row.get("numParticles", "")).strip() == str(num_particles).strip():
                        row_str = ",".join(map(str, row.values()))
                        if resolution_key in row:
                            res_val = row[resolution_key]
                        elif resolution_key == "mean" and "median" in row:
                            res_val = row["median"]
                        elif resolution_key == "median" and "mean" in row:
                            res_val = row["mean"]
                        else:
                            return None, None
                        return row_str, res_val
        except Exception as exc:
            import warnings
            warnings.warn(
                f"[generate_overview] Could not read full-mask CSV "
                f"'{fullmask_csv_path}': {exc}"
            )
        return None, None

    entries = os.listdir(args.directory)

    # find candidate iteration dirs (accept both _janas_ and legacy _emprove_ prefixes)
    janas_dirs = [
        entry for entry in entries
        if (entry.startswith(prefix) or entry.startswith("_emprove_"))
           and "_" in entry
           and entry.split("_")[-1].isdigit()
    ]

    # drop any that don’t have the locres file we need
    required = "bestRanked_locres_values.csv"
    janas_dirs = [
        d for d in janas_dirs
        if os.path.isfile(os.path.join(args.directory, d, required))
    ]

    # nothing valid? bail out
    if not janas_dirs:
        print("ERROR_EMPTY")
        return
     
    # sort by the trailing number
    janas_dirs.sort(key=lambda x: int(x.split("_")[-1]))

    settings_path = os.path.join(args.directory, "session_settings.toml")
    with open(settings_path, "r") as file:
        dataSettings = toml.load(file)

    resolution_mapping = {
        "meanresolution": "mean",
        "highresolution": "max",
        "highresolutionquartile": "highQuartile",
        "lowresolution": "min",
        "lowresolutionquartile": "lowQuartile",
    }
    resolutionBestTarget_string = resolution_mapping.get(
        dataSettings.get("resolutionBestTarget", "meanResolution").lower(),
        "mean"
    )
    assessmentMethod = str(dataSettings.get("assessmentMethod", "mean")).strip().lower()
    resolutionTarget_string = assessmentMethod
    if assessmentMethod == "median" and resolutionBestTarget_string == "mean":
        resolutionBestTarget_string = "median"

    comparison_mode = str(dataSettings.get("comparison_mode", "auto")).strip().lower()

    adaptive_mask_enabled = str(
        dataSettings.get("adaptiveMask", "False")
    ).strip().lower() in ("true", "1", "yes")

    if comparison_mode == "auto":
        if adaptive_mask_enabled:
            comparison_mode = "target_stability"
        else:
            comparison_mode = "resolution"

    # Helper to read the chosen statistic from a CSV row, with mean<->median fallback
    def _get_stat_from_row(row: dict, key: str) -> str:
        if key in row:
            return row[key]
        if key == "mean" and "median" in row:
            return row["median"]
        if key == "median" and "mean" in row:
            return row["mean"]
        raise KeyError(f"'{key}'")


    output_log_file_string = ""
    lastSelectionsNoImprovements = 0
    baseline_num_particles = None
    first_processed = True

    for index, ii in enumerate(janas_dirs):
        if args.verbose:
            print("VERBOSE: checking directory", os.path.join(args.directory, ii))

        sigma_value = extractSigmaFromTag(ii)

        bestRanked_path = os.path.join(args.directory, ii, "bestRanked_locres_values.csv")
        targetNum_path = os.path.join(args.directory, ii, "target_num_of_particles.csv")
        if not os.path.exists(bestRanked_path) or not os.path.exists(targetNum_path):
            print(
                "WARNING: missing relevant files in the directory",
                os.path.join(args.directory, ii),
                ", ignoring it"
            )
            continue

        # Read and merge rows from both CSVs, keyed uniquely by numParticles and the chosen resolution field
        with open(bestRanked_path, mode="r") as file:
            reader = csv.DictReader(file)
            data1 = [row for row in reader]

        with open(targetNum_path, mode="r") as file:
            reader = csv.DictReader(file)
            data2 = [row for row in reader]

        combined_data = {
            f"{row['numParticles']}_{_get_stat_from_row(row, resolutionBestTarget_string)}": row
            for row in (data1 + data2)
        }.values()
        combined_data = sorted(combined_data, key=lambda x: int(x["numParticles"]))

        # Write the merged data back into each iteration folder
        fullPrediction_locres_file = os.path.join(
            args.directory, ii, "fullPrediction_bestRanked_locres_values.csv"
        )
        headers = combined_data[0].keys() if combined_data else []
        with open(fullPrediction_locres_file, mode="w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=headers)
            writer.writeheader()
            writer.writerows(combined_data)

        # Identify the row with the lowest chosen‐resolution value
        lowest_ResolutionTarget = float("inf")
        lowest_ResolutionTarget_row = None
        for row in combined_data:
            value = float(_get_stat_from_row(row, resolutionBestTarget_string))
            if value < lowest_ResolutionTarget:
                lowest_ResolutionTarget = value
                lowest_ResolutionTarget_row = row

        # Identify the row with the largest numParticles
        largest_num_particles = 0
        largest_num_particles_row = None
        for row in combined_data:
            np_val = int(row["numParticles"])
            if np_val > largest_num_particles:
                largest_num_particles = np_val
                largest_num_particles_row = row

        resolution_of_lowest = _get_stat_from_row(lowest_ResolutionTarget_row, resolutionBestTarget_string)        
        numParticles_of_lowest = lowest_ResolutionTarget_row["numParticles"]

        if first_processed:
            resolution_of_largest = _get_stat_from_row(largest_num_particles_row, resolutionBestTarget_string)

        lowest_row_str = ",".join(map(str, lowest_ResolutionTarget_row.values()))
        largest_row_str = ",".join(map(str, largest_num_particles_row.values()))

        if args.verbose:
            print("   lowest_ResolutionTarget_row:", lowest_row_str)
            print("   resolution_of_lowest_ResolutionTarget:", resolution_of_lowest)
            print("   largest_num_particles_row:", largest_row_str)
            if first_processed:
                print("   resolution_of_largest_ResolutionTarget:", resolution_of_largest)
            print("   numParticles_of_lowest_ResolutionTarget:", numParticles_of_lowest)

        toCheckDir = os.path.join(args.directory, ii)

        session_assess_mask = str(dataSettings.get("assessMask", "") or "").strip()
        session_main_mask = str(dataSettings.get("mask", "") or "").strip()
        adaptive_mask_path = os.path.join(args.directory, ii, "adaptive_mask.mrc")
        adaptive_mask_exists = os.path.isfile(adaptive_mask_path)
        bestRanked_fullMask_path = os.path.join(args.directory, ii, "bestRanked_locres_values_fullMask.csv")
        if adaptive_mask_exists:
            assessment_mask_file = adaptive_mask_path
            assessment_mask_mode = "adaptive_mask"
        elif session_assess_mask:
            assessment_mask_file = session_assess_mask
            assessment_mask_mode = "assessMask"
        else:
            assessment_mask_file = session_main_mask
            assessment_mask_mode = "mask"

        # Build this iteration's selection block, including the extracted sigma value
        selection_block = f'tag = "{ii}"\n'
        selection_block += f'SCI_sigma = "{sigma_value}"\n'
        selection_block += f'working_directory = "{toCheckDir}"\n'
        selection_block += f"selection_number = {int(ii.split('_')[-1])}\n"
        selection_block += f"reference_num_particles = {numParticles_of_lowest}\n"
        selection_block += (
            f'reference_starFile = "{os.path.join(args.directory, ii, "norm_" + ii + "_best" + str(numParticles_of_lowest) + ".star")}"\n'
        )
        selection_block += (
            f'reference_mapA = "{os.path.join(args.directory, ii, "norm_" + ii + "_best" + str(numParticles_of_lowest) + "_recH1.mrc")}"\n'
        )
        selection_block += (
            f'reference_mapB = "{os.path.join(args.directory, ii, "norm_" + ii + "_best" + str(numParticles_of_lowest) + "_recH2.mrc")}"\n'
        )
        selection_block += f'reference_locres_stats = "{lowest_row_str}"\n'
        selection_block += f"reference_locres_ResolutionTarget = {resolution_of_lowest}\n"
        selection_block += f"ResolutionTarget = '{resolutionTarget_string}'\n"
        selection_block += (
            f'computed_locres_file = "{os.path.join(args.directory, ii, "fullPrediction_bestRanked_locres_values.csv")}"\n'
        )
        selection_block += f'assessment_mask_file = "{assessment_mask_file}"\n'
        selection_block += f'assessment_mask_mode = "{assessment_mask_mode}"\n'
        if os.path.isfile(bestRanked_fullMask_path):
            selection_block += f'computed_locres_file_fullMask = "{bestRanked_fullMask_path}"\n'
            fm_stats, fm_res = _lookup_fullmask_stats(
                bestRanked_fullMask_path,
                numParticles_of_lowest,
                resolutionBestTarget_string,
            )
            if fm_stats is not None and fm_res is not None:
                selection_block += f'reference_locres_stats_fullMask = "{fm_stats}"\n'
                selection_block += f"reference_locres_ResolutionTarget_fullMask = {fm_res}\n"
            else:
                import warnings
                warnings.warn(
                    f"[generate_overview] Could not resolve full-mask stats for "
                    f"numParticles={numParticles_of_lowest} in {bestRanked_fullMask_path}"
                )
        selection_block += f"selected_target_num_particles = {numParticles_of_lowest}\n"
        selection_block += f'comparison_mode = "{comparison_mode}"\n'

        # If this is the first directory, write the initial “target” block:
        if first_processed:
            baseline_num_particles = largest_num_particles
            output_log_file_string += "\n[[_janas_selection_0]]\n"
            output_log_file_string += f'tag = "{ii}"\n'
            output_log_file_string += f'SCI_sigma = "{sigma_value}"\n'
            output_log_file_string += f'working_directory = "{toCheckDir}"\n'
            output_log_file_string += f"reference_num_particles = {largest_num_particles}\n"
            output_log_file_string += (
                f'reference_starFile = "{os.path.join(args.directory, ii, "norm_" + ii + "_best" + str(largest_num_particles) + ".star")}"\n'
            )
            output_log_file_string += (
                f'reference_mapA = "{os.path.join(args.directory, ii, "norm_" + ii + "_best" + str(largest_num_particles) + "_recH1.mrc")}"\n'
            )
            output_log_file_string += (
                f'reference_mapB = "{os.path.join(args.directory, ii, "norm_" + ii + "_best" + str(largest_num_particles) + "_recH2.mrc")}"\n'
            )
            output_log_file_string += f'reference_locres_stats = "{largest_row_str}"\n'
            output_log_file_string += (
                f"reference_locres_ResolutionTarget = {resolution_of_largest}\n"
            )
            output_log_file_string += f"ResolutionTarget = '{resolutionTarget_string}'\n"
            output_log_file_string += (
                f'computed_locres_file = "{os.path.join(args.directory, ii, "bestRanked_locres_values.csv")}"\n'
            )
            output_log_file_string += f'assessment_mask_file = "{assessment_mask_file}"\n'
            output_log_file_string += f'assessment_mask_mode = "{assessment_mask_mode}"\n'
            if os.path.isfile(bestRanked_fullMask_path):
                output_log_file_string += f'computed_locres_file_fullMask = "{bestRanked_fullMask_path}"\n'
                fm_stats_0, fm_res_0 = _lookup_fullmask_stats(
                    bestRanked_fullMask_path,
                    largest_num_particles,
                    resolutionBestTarget_string,
                )
                if fm_stats_0 is not None and fm_res_0 is not None:
                    output_log_file_string += f'reference_locres_stats_fullMask = "{fm_stats_0}"\n'
                    output_log_file_string += f"reference_locres_ResolutionTarget_fullMask = {fm_res_0}\n"
                else:
                    import warnings
                    warnings.warn(
                        f"[generate_overview] Could not resolve full-mask stats for "
                        f"numParticles={largest_num_particles} in {bestRanked_fullMask_path}"
                    )
            output_log_file_string += f"selected_target_num_particles = {largest_num_particles}\n"
            output_log_file_string += f'comparison_mode = "{comparison_mode}"\n'
            output_log_file_string += "improve_previouses_locres_ResolutionTarget = true\n"
            output_log_file_string += (
                f'session_settings = "{os.path.join(args.directory, "session_settings.toml")}"\n'
            )

            # Record these initial selections so they can be copied later if they remain best
            target_reference_starFile = os.path.join(
                args.directory, ii, f"norm_{ii}_best{largest_num_particles}.star"
            )
            target_reference_mapA = os.path.join(
                args.directory, ii, f"norm_{ii}_best{largest_num_particles}_recH1.mrc"
            )
            target_reference_mapB = os.path.join(
                args.directory, ii, f"norm_{ii}_best{largest_num_particles}_recH2.mrc"
            )
            best_locres = float(resolution_of_largest)
            best_selectionBlock = selection_block
            best_selectionBlock_size = largest_num_particles
            previous_best_num_particles = largest_num_particles
            first_processed = False

        # Decide whether this iteration improved upon the previous best
        improve_previous = "false"

        if comparison_mode == "target_stability":
            # Stability signal (internal): the winning subset size changed.
            # Drives sigma decay via lastSelectionsNoImprovements. This is
            # NOT what gets written to improve_previouses_locres_ResolutionTarget.
            if int(numParticles_of_lowest) != int(previous_best_num_particles):
                previous_best_num_particles = int(numParticles_of_lowest)
                lastSelectionsNoImprovements = 0
            else:
                lastSelectionsNoImprovements += 1

            # Target selection: pick the iteration with the best resolution.
            # The reported improve_previous flag reflects resolution improvement,
            # which is the meaningful signal to the user even in target_stability
            # mode — the stability signal is tracked separately above.
            if float(resolution_of_lowest) < best_locres:
                improve_previous = "true"
                best_locres = float(resolution_of_lowest)
                best_selectionBlock = selection_block
                best_selectionBlock_size = numParticles_of_lowest
                target_reference_starFile = os.path.join(
                    args.directory, ii, f"norm_{ii}_best{numParticles_of_lowest}.star"
                )
                target_reference_mapA = os.path.join(
                    args.directory, ii, f"norm_{ii}_best{numParticles_of_lowest}_recH1.mrc"
                )
                target_reference_mapB = os.path.join(
                    args.directory, ii, f"norm_{ii}_best{numParticles_of_lowest}_recH2.mrc"
                )
        else:
            # Fixed mask: improvement = better resolution
            if float(resolution_of_lowest) < best_locres:
                improve_previous = "true"
                best_locres = float(resolution_of_lowest)
                best_selectionBlock = selection_block
                best_selectionBlock_size = numParticles_of_lowest
                target_reference_starFile = os.path.join(
                    args.directory, ii, f"norm_{ii}_best{numParticles_of_lowest}.star"
                )
                target_reference_mapA = os.path.join(
                    args.directory, ii, f"norm_{ii}_best{numParticles_of_lowest}_recH1.mrc"
                )
                target_reference_mapB = os.path.join(
                    args.directory, ii, f"norm_{ii}_best{numParticles_of_lowest}_recH2.mrc"
                )
                lastSelectionsNoImprovements = 0
            else:
                lastSelectionsNoImprovements += 1

        output_log_file_string += f"\n[[_janas_selection_{int(ii.split('_')[-1])}]]\n"
        output_log_file_string += selection_block
        output_log_file_string += f"improve_previouses_locres_ResolutionTarget = {improve_previous}\n\n"

    # Copy the final chosen reference files into the top‐level directory
    shutil.copy(
        target_reference_starFile,
        os.path.join(args.directory, "reference_subset.star")
    )
    shutil.copy(
        target_reference_mapA,
        os.path.join(args.directory, "reference_subset_mapA.mrc")
    )
    shutil.copy(
        target_reference_mapB,
        os.path.join(args.directory, "reference_subset_mapB.mrc")
    )
    starHandler.removeColumnsTagsStartingWith(
        os.path.join(args.directory, "reference_subset.star"),
        os.path.join(args.directory, "reference_subset_clean.star"),
        "_janas",
    )

    # Average the two half‐maps to produce reference_subset_map.mrc
    sizeMap = janas_core.sizeMRC(target_reference_mapA)
    map1 = np.array(janas_core.ReadMRC(target_reference_mapA))
    map2 = np.array(janas_core.ReadMRC(target_reference_mapB))
    map12 = (map1 + map2) / 2
    spacingMRC = round(janas_core.spacingMRC(target_reference_mapA), 4)
    janas_core.WriteMRC(
        map12.flatten().tolist(),
        os.path.join(args.directory, "reference_subset_map.mrc"),
        sizeMap[2],
        sizeMap[1],
        sizeMap[0],
        spacingMRC,
    )
    janas_core.replaceMrcHeader(
        target_reference_mapA,
        os.path.join(args.directory, "reference_subset_map.mrc")
    )

    # Compute percentage of particles retained
    percentage_particles_retained = 0
    denom = baseline_num_particles if baseline_num_particles else largest_num_particles  # <-- CHANGE THIS
    if denom and denom > 0:
        percentage_particles_retained = (float(best_selectionBlock_size) / float(denom)) * 100.0

#    if largest_num_particles > 0:
#        percentage_particles_retained = (
#            float(best_selectionBlock_size) / float(largest_num_particles)
#        ) * 100

    # Finally, write the overview.txt with the accumulated blocks and the best‐selection block
    with open(os.path.join(args.directory, "overview.txt"), "w") as toml_file:
        toml_file.write("# Selections overview\n\n")
        toml_file.write("[[_janas_target_selection]]\n")
        toml_file.write(best_selectionBlock)
        toml_file.write(
            f"last_consecutive_non_improving_selections = {lastSelectionsNoImprovements}\n"
        )
        toml_file.write(
            f"percentage_particles_retained = {percentage_particles_retained:.2f}\n"
        )
        toml_file.write(output_log_file_string)

    remove_duplicates_toml(os.path.join(args.directory, "overview.txt"))
    return janas_dirs



janas_getTarget = command.add_parser(
    "getTarget",
    description="get the target from a specific overview file (overview.txt)",
    help="get the target from a specific overview file (overview.txt)",
)
janas_getSettings = command.add_parser(
    "getSettings",
    description="Retrieve a single setting from a session_settings.toml",
    help="getSettings"
)
janas_getSettings.add_argument(
    "--settingsFile",
    required=True,
    type=str,
    help="Path to the session_settings.toml file"
)
janas_getSettings.add_argument(
    "--setting",
    required=True,
    type=str,
    help="Name of the setting to retrieve"
)
def getSettings(args):
    """
    Load the TOML file at args.settingsFile, look up args.setting,
    and print its value (or None if the key is not present).
    """
    try:
        data = toml.load(args.settingsFile)
    except Exception as e:
        print(f"Error reading settings file: {e}")
        return
    value = data.get(args.setting)
    print(value if value is not None else None)

janas_getTarget.add_argument(
    "--overviewFile",
    required=True,
    type=str,
    help="toml overview file to acquire target information",
)
janas_getTarget.add_argument(
    "--particles",
    action="store_true",
    help="gets the target star file with relevant particles",
)
janas_getTarget.add_argument(
    "--map1", action="store_true", help="gets the target map1 file"
)
janas_getTarget.add_argument(
    "--map2", action="store_true", help="gets the target map2 file"
)
janas_getTarget.add_argument(
    "--stats",
    action="store_true",
    help="num particles, min, quartile, mean, quartile, max",
)
janas_getTarget.add_argument(
    "--unimproved",
    action="store_true",
    help="tell if the result has improved, if so it gets true, if not gets false",
)
janas_getTarget.add_argument("--sigma", action="store_true", help="get the target sigma")
janas_getTarget.add_argument("--current_sigma", action="store_true", help="get the sigma of the current selection id")
janas_getTarget.add_argument("--locres_file", action="store_true", help="get the target locres file")
janas_getTarget.add_argument("--current_selection_ID", action="store_true", help="get the current selection ID")
janas_getTarget.add_argument(
    "--reference_num_particles",
    action="store_true",
    help="print the reference_num_particles from the target selection"
)
janas_getTarget.add_argument(
    "--assessment_mask_file",
    action="store_true",
    help="print the assessment_mask_file path from the target selection"
)
def getTarget(args):
    """Retrieve target information from an overview file.

    This helper is robust to "minimal" or partially-written overview files.
    When the overview does not yet contain the reference_starFile entry, it
    falls back to a session_settings.toml located in the same directory as
    the overview file and returns the 'particles' field.
    """
    overview_path = args.overviewFile
    overview_dir = os.path.dirname(os.path.abspath(overview_path))

    # Try to read the TOML overview, but do not fail hard if it is empty
    # or only contains section headers.
    try:
        data = toml.load(overview_path)
    except Exception:
        data = {}

    # Safely extract the target selection, supporting both the "array of tables"
    # TOML representation and a single table.
    target_selection = {}
    if isinstance(data, dict):
        section = data.get("_janas_target_selection")
        if isinstance(section, list) and section:
            target_selection = section[0] or {}
        elif isinstance(section, dict):
            target_selection = section or {}

    # ------------------------------------------------------------------
    # Particles: support minimal overviews by falling back to the
    # session_settings.toml file when reference_starFile is absent.
    # ------------------------------------------------------------------
    if args.particles:
        target_starfile = target_selection.get("reference_starFile")

        # If the overview is still minimal (e.g. only headers) the entry
        # will be missing or empty. In that case use the particles field
        # from the session_settings.toml that was used to start the run.
        if not target_starfile:
            settings_path = os.path.join(overview_dir, "session_settings.toml")
            if os.path.exists(settings_path):
                try:
                    settings = toml.load(settings_path)
                    target_starfile = settings.get("particles")
                except Exception:
                    # Silence TOML errors here; we will just return None.
                    pass

        print(target_starfile if target_starfile is not None else None)

    # ------------------------------------------------------------------
    # Maps and stats: keep the original behaviour but guard against
    # missing keys so that an incomplete overview does not crash.
    # ------------------------------------------------------------------
    if args.map1:
        reference_mapA = target_selection.get("reference_mapA")
        print(reference_mapA if reference_mapA is not None else None)

    if args.map2:
        reference_mapB = target_selection.get("reference_mapB")
        print(reference_mapB if reference_mapB is not None else None)

    if args.stats:
        reference_stats = target_selection.get("reference_locres_stats")
        print(reference_stats if reference_stats is not None else None)

    if args.unimproved:
        last_consecutive_non_improving_selections = target_selection.get(
            "last_consecutive_non_improving_selections"
        )
        print(
            last_consecutive_non_improving_selections
            if last_consecutive_non_improving_selections is not None
            else None
        )

    if args.sigma:
        sigma_str = target_selection.get("SCI_sigma")
        if sigma_str is not None:
            try:
                print(f"{float(sigma_str):.2f}")
            except ValueError:
                # If it cannot be cast to float, just echo the raw value.
                print(sigma_str)
        else:
            print(None)

    # ------------------------------------------------------------------
    # current_selection_ID and current_sigma work directly from the
    # text of the overview file and do not require target_selection.
    # ------------------------------------------------------------------
    if args.current_selection_ID:
        max_id = None
        pattern = re.compile(r'^\[\[_janas_selection_(\d+)\]\]')
        with open(args.overviewFile, "r") as f:
            for line in f:
                m = pattern.match(line.strip())
                if m:
                    idx = int(m.group(1))
                    if max_id is None or idx > max_id:
                        max_id = idx
        print(max_id if max_id is not None else None)
        return

    if args.current_sigma:
        max_id = None
        pattern = re.compile(r'^\[\[_janas_selection_(\d+)\]\]')
        with open(args.overviewFile, "r") as f:
            for line in f:
                m = pattern.match(line.strip())
                if m:
                    idx = int(m.group(1))
                    if max_id is None or idx > max_id:
                        max_id = idx

        if max_id is None:
            print(None)
            return

        try:
            data = toml.load(args.overviewFile)
        except Exception:
            print(None)
            return

        section = f"_janas_selection_{max_id}"
        section_data = data.get(section)

        if isinstance(section_data, list) and section_data:
            entry = section_data[0]
        elif isinstance(section_data, dict):
            entry = section_data
        else:
            entry = {}

        sigma_str = entry.get("SCI_sigma")
        if sigma_str is not None:
            try:
                print(f"{float(sigma_str):.2f}")
            except ValueError:
                print(sigma_str)
        else:
            print(None)
        return

    if args.locres_file:
        # Try to read the explicit field if it exists (new overviews)
        computed_locres_file = target_selection.get("computed_locres_file")

        # Fallback for older / minimal overviews that do not have this field
        if not computed_locres_file:
            working_dir = target_selection.get("working_directory")
            if working_dir:
                # Prefer the "new" merged file if present, otherwise the original one
                candidate_full = os.path.join(
                    working_dir, "fullPrediction_bestRanked_locres_values.csv"
                )
                candidate_best = os.path.join(
                    working_dir, "bestRanked_locres_values.csv"
                )

                if os.path.exists(candidate_full):
                    computed_locres_file = candidate_full
                elif os.path.exists(candidate_best):
                    computed_locres_file = candidate_best

        # Final guard: do not raise, just print None if we still have nothing
        print(computed_locres_file if computed_locres_file else None)

    if args.reference_num_particles:
        val = target_selection.get("reference_num_particles")
        print(val if val is not None else "")

    if args.assessment_mask_file:
        val = target_selection.get("assessment_mask_file")
        print(val if val is not None else "")


# -------------------------------------------------------------------
# settingBasedOptimization: compute next sigma according to session settings
janas_settingBasedOptimization = command.add_parser(
    "settingBasedOptimization",
    description="Adjust sigma using session settings",
    help="settingBasedOptimization"
)
# add its own sub-commands
setting_subparsers = janas_settingBasedOptimization.add_subparsers(
    dest="setting_action"
)

# -------------------------------------------------------------------
# 'next_sigma' sub-command
next_sigma_parser = setting_subparsers.add_parser(
    "next_sigma",
    description="Compute the next sigma based on improvement history",
    help="next_sigma"
)
next_sigma_parser.add_argument(
    "--settingsFile",
    required=True,
    type=str,
    help="Path to session_settings.toml"
)
next_sigma_parser.add_argument(
    "--lastNonImproving",
    required=True,
    type=int,
    help="Number of consecutive non-improving selections"
)
next_sigma_parser.add_argument(
    "--overviewFile",
    required=True,
    type=str,
    help="Path to overview.txt to extract the current SCI_sigma (overrides sigma in settings)"
)


def settingBasedOptimization_next_sigma(args):
    """
    Load the TOML at args.settingsFile, then load args.overviewFile to
    extract the current selection’s SCI_sigma. Apply sigma_decreasing_step
    according to the improvement flags, clamp to minimum, and print to two
    decimal places.
    """
    # 1. Load session settings
    try:
        cfg = toml.load(args.settingsFile)
    except Exception as e:
        print(f"Error reading settings file: {e}")
        return

    # 2. Load overview and extract current SCI_sigma
    max_id = None
    pattern = re.compile(r'^\[\[_janas_selection_(\d+)\]\]')
    with open(args.overviewFile, 'r') as f:
        for line in f:
            m = pattern.match(line.strip())
            if m:
                idx = int(m.group(1))
                if max_id is None or idx > max_id:
                    max_id = idx

    if max_id is None:
        print("No selections found in overview file")
        return

    ov = toml.load(args.overviewFile)
    section = f"_janas_selection_{max_id}"
    sec_data = ov.get(section, [])
    if isinstance(sec_data, list) and sec_data:
        entry = sec_data[0]
    elif isinstance(sec_data, dict):
        entry = sec_data
    else:
        entry = {}

    try:
        current_sigma = float(entry.get("SCI_sigma"))
    except (TypeError, ValueError):
        print("Invalid or missing SCI_sigma in overview file")
        return

    # 3. Parse step and minimum
    try:
        step    = float(cfg.get("sigma_decreasing_step", 0.0))
        minimum = float(cfg.get("minimum_sigma_allowed", 0.0))
    except ValueError:
        print("Configuration error: non-numeric sigma_decreasing_step or minimum_sigma_allowed")
        return

    # 4. Normalise boolean flags
    def to_bool(v):
        if isinstance(v, bool):
            return v
        return str(v).lower() in ("true", "1", "yes")

    decay_with_improve    = to_bool(cfg.get("sigma_decrease_with_improvement", False))
    decay_without_improve = to_bool(cfg.get("sigma_decrease_without_improvement", False))

    non_improved = args.lastNonImproving > 0

    # 5. Compute new sigma
    if non_improved and decay_without_improve:
        new_sigma = current_sigma - step
    elif not non_improved and decay_with_improve:
        new_sigma = current_sigma - step
    else:
        new_sigma = current_sigma

    # 6. Clamp and print
    if new_sigma < minimum:
        new_sigma = minimum

    print(f"{new_sigma:.2f}")

# -------------------------------------------------------------------
# 'check_early_termination' sub-command
check_early_termination_parser = setting_subparsers.add_parser(
    "check_early_termination",
    description="Return true if unimproving iterations ≥ threshold",
    help="check_early_termination"
)
check_early_termination_parser.add_argument(
    "--settingsFile",
    required=True,
    type=str,
    help="Path to session_settings.toml"
)
check_early_termination_parser.add_argument(
    "--current_unimproving_iterations",
    required=True,
    type=int,
    help="Number of consecutive non-improving iterations"
)

def settingBasedOptimization_check_early_termination(args):
    """
    Load num_unimproving_iterations_for_early_termination from settingsFile,
    compare to args.current_unimproving_iterations, and print 'true'
    or 'false'.
    """
    try:
        cfg = toml.load(args.settingsFile)
    except Exception as e:
        print(f"Error reading settings file: {e}")
        return

    threshold = cfg.get("num_unimproving_iterations_for_early_termination")
    try:
        threshold = int(threshold)
    except Exception:
        print("Configuration error: num_unimproving_iterations_for_early_termination must be an integer")
        return

    result = args.current_unimproving_iterations >= threshold
    print("true" if result else "false")

# -------------------------------------------------------------------
# 'subset_particle_check' sub-command
subset_particle_check_parser = setting_subparsers.add_parser(
    "subset_particle_check",
    description="Return true if unimproving iterations > threshold for selecting best subset",
    help="subset_particle_check"
)
subset_particle_check_parser.add_argument(
    "--settingsFile",
    required=True,
    type=str,
    help="Path to session_settings.toml"
)
subset_particle_check_parser.add_argument(
    "--current_unimproving_iterations",
    required=True,
    type=int,
    help="Number of consecutive non-improving iterations"
)

def settingBasedOptimization_subset_particle_check(args):
    """
    Load num_unimproving_iterations_before_selecting_from_best_subset from settingsFile,
    compare to args.current_unimproving_iterations, and print 'true' or 'false'.
    """
    try:
        cfg = toml.load(args.settingsFile)
    except Exception as e:
        print(f"Error reading settings file: {e}")
        return

    key = "num_unimproving_iterations_before_selecting_from_best_subset"
    threshold = cfg.get(key)
    try:
        threshold = int(threshold)
    except Exception:
        print(f"Configuration error: {key} must be an integer")
        return

    result = args.current_unimproving_iterations > threshold
    print("true" if result else "false")


#################################
## progress (HTML dashboard)
janas_progress_parser = command.add_parser(
    "progress",
    description=(
        "Write a single progress.html in the session directory summarising "
        "events.ndjson, status.txt, step_timings.csv and overview.txt."
    ),
    help="generate progress.html for the JANAS session",
)
janas_progress_parser.add_argument(
    "--session", default=None,
    help="JANAS session directory (the one containing overview.txt and "
         "runtime/). If omitted, defaults to the current directory or, "
         "with --overview, the parent of the overview file."
)
janas_progress_parser.add_argument(
    "--overview", default=None,
    help="Path to overview.txt; the session directory is inferred from it. "
         "Mutually exclusive with --session."
)
janas_progress_parser.add_argument(
    "--refresh", type=int, default=15,
    help="Seconds for the HTML <meta http-equiv='refresh'> tag. 0 disables "
         "auto-refresh. The dashboard also drops the meta tag automatically "
         "once the session is marked finished, regardless of this value "
         "(default: 15)."
)
janas_progress_parser.add_argument(
    "--max-events", dest="max_events", type=int, default=100,
    help="Maximum number of events rendered in the 'Recent events' section "
         "(default: 100)."
)
janas_progress_parser.add_argument(
    "--quiet", action="store_true",
    help="Do not print the output file path on success."
)


def main(command_line=None):
    args = janas_parser.parse_args(command_line)
    if args.command == "getNumParticles":
        getNumParticles(args)
    elif args.command == "automaticParticleSubsets":
        automaticParticleSubsets(args)
    elif args.command == "logAnalyzer":
        logAnalyzer(args)
    elif args.command == "plotOverview":
        plotOverview(args)
    elif args.command == "generate_overview":
        generate_overview(args)
    elif args.command == "getTarget":
        getTarget(args)
    elif args.command == "getSettings":
        getSettings(args)
    elif args.command == "progress":
        from janas import janas_progress
        return janas_progress.cmd_progress(args)
    elif args.command == "settingBasedOptimization":
        if args.setting_action == "next_sigma":
            settingBasedOptimization_next_sigma(args)
        elif args.setting_action == "check_early_termination":
            settingBasedOptimization_check_early_termination(args)
        elif args.setting_action == "subset_particle_check":
            settingBasedOptimization_subset_particle_check(args)
        else:
            janas_parser.print_help()
    else:
        janas_parser.print_help()


if __name__ == "__main__":
    main()

