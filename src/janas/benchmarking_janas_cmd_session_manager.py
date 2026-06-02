
# -*- coding: utf-8 -*-

# File: janas_cmd_session_manager.py 
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology


"""
Module: janas_session_manager.py

Manages JANAS selection sessions and automates reconstruction workflows:
- Merge and sort CSV results
- Generate shell scripts for particle reconstruction and local‐resolution evaluation
- Initialize new iterative selection sessions based on user parameters
- Support random and classification‐based sessions
"""

# Standard library
import argparse
import os
import stat

# Third-party
import numpy as np
import pandas as pd
import scipy.stats as stats
import toml

# Local
from janas import starHandler
from janas import utils
from janas.version import get_version

class JANASHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawTextHelpFormatter,
):
    pass

janas_parser = argparse.ArgumentParser(
    prog="janas_session_manager",
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


def merge_and_sort_csv(file1, file2, output_file):
    # Read the CSV files
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    # Merge the two dataframes
    merged_df = pd.concat([df1, df2])

    # Sort by 'numParticles'
    sorted_df = merged_df.sort_values(by="numParticles")

    # Save the sorted dataframe to a new CSV file
    sorted_df.to_csv(output_file, index=False)

    return sorted_df


##########################################
##########################################
##### janas_produce_reconstructions_script
##########################################
janas_produce_reconstructions_script = command.add_parser(
    "produce_reconstructions_script",
    description="produce_reconstructions_script",
    help="Reconstruction Evaluation",
)
janas_produce_reconstructions_script.add_argument(
    "--i", required=True, type=str, help="input star file with particles to score"
)
janas_produce_reconstructions_script.add_argument(
    "--outDir", required=False, default="./", type=str, help="output Dir"
)
janas_produce_reconstructions_script.add_argument(
    "--tagRank", required=False, type=str, help="tag for the ranked particles"
)
janas_produce_reconstructions_script.add_argument(
    "--mask", required=True, type=str, help="mask to use for evaluation"
)
janas_produce_reconstructions_script.add_argument(
    "--manualParticleSubsets",
    required=False,
    type=str,
    default=None,
    help="Comma separated list of manual particles subsets",
)
janas_produce_reconstructions_script.add_argument(
    "--scriptName",
    required=False,
    default="script_reconstructions.sh",
    type=str,
    help="script name as output",
)
janas_produce_reconstructions_script.add_argument(
    "--resultFilename",
    required=False,
    default="bestRanked_locres_values.csv",
    type=str,
    help="filename with results",
)
janas_produce_reconstructions_script.add_argument(
    "--mode",
    required=False,
    default="bestRanked",
    type=str,
    help="mode for selecting particles choose: bestRanked(default), or random",
)
# janas_produce_reconstructions_script.add_argument("--mask_XYZ_boxsize", required=False, default='FALSE', type=str, help="X,Y,Z,boxsize coordinates for masked locres")
janas_produce_reconstructions_script.add_argument(
    "--masked_crop", action="store_true", help="check if automatic masked crop"
)


def produce_reconstructions_script(args):
    if not os.path.isfile(args.i):
        print('ERROR: file "', args.i, '" not existing')
        exit()
    elif not os.path.exists(args.outDir):
        try:
            os.makedirs(args.outDir)
        except OSError as e:
            print(f"Error: {e.strerror}")

    # print("num_non_null_items=",num_non_null_items)

    if args.manualParticleSubsets is not None:
        print("Manual Particle Subset Selection")
        numParticlesList = args.manualParticleSubsets.split(",")

    #    elif args.automaticParticleSubsets:
    # print ("Automatic Particle Subset Selection,")
    # starFile=starHandler.readStar(args.i)
    # num_non_null_items = int(starFile.count()[0])
    # numParticleSubsetsSelected=int(args.automaticParticleSubsets[0])
    # expectedNumberOfPartilces=int(args.automaticParticleSubsets[1])
    # standardDeviation=int(args.automaticParticleSubsets[2])
    # print("numParticleSubsetsSelected=", args.automaticParticleSubsets[0])
    # print("expectedNumberOfPartilces=", args.automaticParticleSubsets[1])
    # print("num_non_null_items=", num_non_null_items)
    # numParticlesList = assessParticles.automaticParticleSubsetSelection(numParticleSubsetsSelected, expectedNumberOfPartilces, num_non_null_items, standardDeviation)

    print(numParticlesList)

    file_selection_suffix = "best"
    if args.mode == "random":
        file_selection_suffix = "random"



    #######################
    #####RECONSTRUCTION
    reconstruction_command = """#!/bin/bash
\n\n
##############################
#######  RECONSTRUCTIONS
rec_subset() {
        fileIn=$1
        fileOut_basename=$2
        subset=$3
        if [ -f ${fileOut_basename}_recH${subset}.mrc ]; then
            echo "DOING NOTHING: Reconstructed file ${fileOut_basename}_recH${subset}.mrc exists"
        else
            echo "DOING Reconstruction for file ${fileOut_basename}_recH${subset}.mrc"
            #mpirun --np 28 relion_reconstruct --i ${fileIn} --o  ${fileOut_basename}_recH${subset}.mrc --subset ${subset} --ctf
            relion_reconstruct --i ${fileIn} --o  ${fileOut_basename}_recH${subset}.mrc --subset ${subset} --ctf &
            sleep 40
        fi
}
    \n\n"""
    numParticlesListStr = ",".join(map(str, numParticlesList))
    reconstruction_command += "numParticlesCsv=" + numParticlesListStr + "\n"
    reconstruction_command += """\n
for numParticles in $(echo $numParticlesCsv | sed "s/,/ /g")
do\n"""
    outFile = os.path.join(
        args.outDir,
        "norm_"
        + str(os.path.split(args.outDir)[-1])
        + "_"
        + file_selection_suffix
        + "${numParticles}",
    )
    if args.mode == "bestRanked":
        reconstruction_command += (
            "    janas selectBestRanked --i "
            + args.i
            + " --o "
            + outFile
            + ".star --num  ${numParticles} \n"
        )
    elif args.mode == "random":
        reconstruction_command += (
            "    janas selectRandom --i "
            + args.i
            + " --o "
            + outFile
            + ".star --num  ${numParticles} \n"
        )
    reconstruction_command += (
        "    rec_subset  " + outFile + ".star  " + outFile + " 1  \n"
    )
    reconstruction_command += (
        "    rec_subset  " + outFile + ".star  " + outFile + " 2  \n"
    )
    reconstruction_command += """\n
    echo ${scorelabelToNormalize}
done
    \n\n"""
    reconstruction_command += "wait\n"

    reconstruction_file_path = os.path.join(args.outDir, args.scriptName)
    with open(reconstruction_file_path, "w") as f:
        f.write(reconstruction_command)
    try:
        os.chmod(
            reconstruction_file_path,
            os.stat(reconstruction_file_path).st_mode | stat.S_IXUSR,
        )
    except PermissionError:
        pass

    ########################
    #####LOCRES COMMAND
    locres_command = """
##############################
#######  LOCRES
locres() {
        file_basename=$1        
        if [ -f ${file_basename}_locres.mrc ]; then
            echo "DOING NOTHING: locres file ${file_basename}_locres.mrc exists"
        else
            echo "DOING locres for file ${file_basename}_locres.mrc"
            mpirun --np 28 relion_postprocess_mpi --i ${file_basename}_recH1.mrc --i2 ${file_basename}_recH2.mrc --o  ${file_basename} --locres #--locres_thresholdFSC 0.5
            rm ${file_basename}_locres_fscs.star
            rm ${file_basename}_locres_filtered.mrc
        fi
}

crop_image() {
  local imageIn="$1"
  local imageOut="$2"
  trimvol -x 100,300 -y 100,300 -z 100,300  "$imageIn" "$imageOut"
}

    \n\n"""
    numParticlesListStr = ",".join(map(str, numParticlesList))
    locres_command += "numParticlesCsv=" + numParticlesListStr + "\n"
    locres_command += """
for numParticles in $(echo $numParticlesCsv | sed "s/,/ /g")
do\n"""

    outFile = os.path.join(
        args.outDir,
        "norm_"
        + str(os.path.split(args.outDir)[-1])
        + "_"
        + file_selection_suffix
        + "${numParticles}",
    )
    if args.masked_crop:
        locres_command += (
            "    janas_utils maskedCrop --mask   "
            + args.mask
            + "  --padding 8 --i  "
            + outFile
            + "_recH1.mrc --o "
            + outFile
            + "_crop_recH1.mrc\n"
        )
        locres_command += (
            "    janas_utils maskedCrop --mask   "
            + args.mask
            + "  --padding 8 --i  "
            + outFile
            + "_recH2.mrc --o "
            + outFile
            + "_crop_recH2.mrc\n"
        )
        locres_command += "    locres  " + outFile + "_crop   \n"
    else:
        locres_command += "    locres  " + outFile + "   \n"
    locres_command += """
    echo ${scorelabelToNormalize}
done
    \n\n"""
    reconstruction_file_path = os.path.join(args.outDir, args.scriptName)
    with open(reconstruction_file_path, "a") as f:
        f.write(locres_command)

    ########################
    #####ASSESS LOCRES
    assess_command = """
############################
#######  ASSESS locres
"""
    maskLocresFilename = args.mask
    print("mask location=", args.mask)
    locresSuffixfix = "_locres"
    if args.masked_crop:
        locresSuffixfix = "_locres"
        print("DOING CROP!")
        maskLocresFilename = os.path.join(args.outDir, "mask_crop.mrc")
        assess_command += (
            "janas_utils maskedCrop --mask "
            + args.mask
            + " --padding 8 --i "
            + args.mask
            + "  --o "
            + maskLocresFilename
            + "\n"
        )
        locresSuffixfix = "_crop_locres"
    else:
        print("Not doing DOING CROP!")

    numParticlesListStr = ",".join(map(str, numParticlesList))
    result_filename = os.path.join(args.outDir, args.resultFilename)
    assess_command += "numParticlesCsv=" + numParticlesListStr + "\n"
    assess_command += (
        "echo numParticles,max,highQuartile,mean,lowQuartile,min >"
        + result_filename
        + "\n"
    )
    assess_command += """
for numParticles in $(echo $numParticlesCsv | sed "s/,/ /g")
do\n"""
    assess_command += 'printf "%s," ${numParticles} >> ' + result_filename + "\n"
    outFile = os.path.join(
        args.outDir,
        "norm_"
        + str(os.path.split(args.outDir)[-1])
        + "_"
        + file_selection_suffix
        + "${numParticles}",
    )
    assess_command += (
        "janas_app_meanMinMax  "
        + outFile
        + locresSuffixfix
        + ".mrc  "
        + maskLocresFilename
        + " >> "
        + result_filename
        + "\n"
    )
    assess_command += """
    echo ${scorelabelToNormalize}
done
    \n\n"""
    reconstruction_file_path = os.path.join(args.outDir, args.scriptName)
    with open(reconstruction_file_path, "a") as f:
        f.write(assess_command)


def find_janas_ranking_tag(elements):
    """
    Search for a column name that starts with '_janas_' or '_emprove_' (backward compat) and ends with '_norm' followed by a number.
    """
    import re

    for element in elements:
        if re.match(r"(_janas_|_emprove_).*_norm\d+$", element):
            return element
    return None


#######################################################
## GENERATE PARTICLE SELECTION run script
def generate_run_script(fileSettings):
    # Reuse the unified generator and only switch the backend
    # for reconstructions and local-resolution estimation to RELION.
    return generate_run_script_no_external(
        fileSettings,
        gpu_list=None,
        use_external_relion=True,
    )

def _assert_replace_once(script: str, old: str, new: str, label: str) -> str:
    """Replace exactly one occurrence of *old* in *script*; raise RuntimeError otherwise.

    All script-injection points use this helper so that a stale or mis-typed anchor
    fails loudly at script-generation time rather than producing a silently broken
    run script.
    """
    count = script.count(old)
    if count != 1:
        raise RuntimeError(
            f"[janas-log] Script injection '{label}': "
            f"expected exactly 1 occurrence of anchor string, found {count}.\n"
            f"Anchor: {old!r}"
        )
    return script.replace(old, new, 1)


def _runtime_logging_shell_block() -> str:
    """Return the shell runtime-logging helper block to inject into generated run scripts.

    Produces five bash functions:
      emit_runtime_event   - appends one NDJSON record to workingDir/runtime/events.ndjson
      write_runtime_status - overwrites workingDir/runtime/status.txt (human-readable)
      init_runtime_logging - creates runtime/ dir, initialises CSV header, emits session_start
      run_step             - wraps a scientific command with per-step log + events + CSV row
      finish_runtime_logging - emits session_end event and writes final status
    """
    return r"""
# ==============================
# Runtime logging helpers
# ==============================
# Artifacts:   workingDir/runtime/status.txt        (overwritten each step; watch -n 2 cat)
#              workingDir/runtime/events.ndjson      (append-only NDJSON event stream)
#              workingDir/runtime/step_timings.csv   (one row per completed step)
# Per-step:    workingDir/<tag>/steps/<NN>_<step>.log

emit_runtime_event() {
  # Caller sets _EV_* env vars; this function appends one JSON line to events.ndjson.
  # Uses python3 heredoc to avoid shell-JSON escaping hazards.
  # try/except inside Python ensures filesystem errors never abort the scientific run.
  _EV_WORKDIR="${workingDir}" python3 - <<'PY'
import json, os
from datetime import datetime, timezone
e = os.environ
rec = {k: v for k, v in {
    "event":          e.get("_EV_TYPE"),
    "timestamp":      datetime.now(timezone.utc).isoformat(),
    "iteration":      e.get("_EV_ITE"),
    "step":           e.get("_EV_STEP"),
    "status":         e.get("_EV_STATUS"),
    "rc":             e.get("_EV_RC"),
    "t_start":        e.get("_EV_T_START"),
    "t_end":          e.get("_EV_T_END"),
    "elapsed_s":      e.get("_EV_ELAPSED"),
    "log_path":       e.get("_EV_LOG"),
    "cmd":            e.get("_EV_CMD"),
    "tag":            e.get("_EV_TAG"),
    "hostname":       os.uname().nodename,
    "slurm_job_id":   e.get("SLURM_JOB_ID") or None,
    "slurm_nodelist": e.get("SLURM_NODELIST") or None,
    "slurm_cpus":     e.get("SLURM_CPUS_ON_NODE") or None,
    "cuda_devices":   e.get("CUDA_VISIBLE_DEVICES") or None,
    "working_dir":    e.get("_EV_WORKDIR"),
}.items() if v is not None}
wd = e.get("_EV_WORKDIR", ".")
try:
    with open(wd + "/runtime/events.ndjson", "a") as f:
        f.write(json.dumps(rec) + "\n")
except Exception as ex:
    import sys; print("[janas-log] Warning: " + str(ex), file=sys.stderr)
PY
}

write_runtime_status() {
  local _wrs_ite="$1" _wrs_step="$2" _wrs_status="$3"
  local _wrs_log="$4" _wrs_tstart="$5" _wrs_tend="${6:-running...}"
  {
    echo "========================================"
    echo " JANAS Session Status"
    echo "========================================"
    echo " Session   : ${workingDir}"
    echo " Iteration : ${_wrs_ite}"
    echo " Step      : ${_wrs_step}"
    echo " Status    : ${_wrs_status}"
    echo " Log       : ${_wrs_log}"
    echo " Started   : ${_wrs_tstart}"
    echo " Finished  : ${_wrs_tend}"
    echo " Host      : $(hostname 2>/dev/null || echo unknown)"
    [ -n "${SLURM_JOB_ID:-}" ] && echo " SLURM Job : ${SLURM_JOB_ID}"
    [ -n "${SLURM_NODELIST:-}" ] && echo " Nodes     : ${SLURM_NODELIST}"
    echo " Updated   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "========================================"
  } > "${workingDir}/runtime/status.txt" 2>/dev/null || true
}

_JANAS_SESSION_DONE=0

_janas_exit_trap() {
  # Fires on all exits: set -e abort, normal exit, or signal.
  # If finish_runtime_logging was not called (i.e. the session aborted mid-run),
  # emit a session_end event with status "aborted" so the event log is never truncated.
  if [ "${_JANAS_SESSION_DONE:-0}" -eq 0 ]; then
    local _et; _et=$(date -u +%Y-%m-%dT%H:%M:%SZ)
    _EV_TYPE="session_end" _EV_ITE="--" _EV_STEP="session_end" \
      _EV_STATUS="aborted" _EV_T_END="${_et}" \
      emit_runtime_event
    write_runtime_status "--" "session_end" "aborted" "" "" "${_et}"
    echo "[janas-log] Session aborted: ${workingDir}/runtime/" >&2
  fi
}

init_runtime_logging() {
  local _rt="${workingDir}/runtime"
  mkdir -p "${_rt}"
  # Only write CSV header if file does not exist (supports mid-run resume)
  [ -f "${_rt}/step_timings.csv" ] || \
    printf 'iteration,step,t_start,t_end,elapsed_s,rc\n' > "${_rt}/step_timings.csv"
  [ -f "${_rt}/events.ndjson" ] || touch "${_rt}/events.ndjson"
  # Install EXIT trap BEFORE emitting session_start so every abort is covered.
  _JANAS_SESSION_DONE=0
  trap '_janas_exit_trap' EXIT
  local _t; _t=$(date -u +%Y-%m-%dT%H:%M:%SZ)
  write_runtime_status "--" "init" "starting" "" "${_t}" ""
  _EV_TYPE="session_start" _EV_ITE="0" _EV_STEP="session_start" \
    _EV_STATUS="started" _EV_T_START="${_t}" \
    emit_runtime_event
  echo "[janas-log] Runtime logging started: ${_rt}/"
}

run_step() {
  # Usage: run_step <ite> <step_name> <log_dir> <command...>
  # Wraps one scientific command:
  #   - Creates per-step log file at log_dir/step_name.log
  #   - Emits step_start / step_end NDJSON events
  #   - Updates status.txt on start and completion
  #   - Appends a timing row to step_timings.csv
  #   - Propagates the command's exit code to the caller (set -e will abort on failure)
  #
  # Shell-mode assumption: this helper is designed for scripts running under
  #   set -eo pipefail
  # It temporarily does  set +e  to capture PIPESTATUS[0] from the tee pipeline,
  # then immediately restores  set -e  and returns the captured rc.  When the
  # caller is running with  set -e , a non-zero return will abort the script,
  # which is the desired behaviour.  Do NOT use this function in a script that
  # runs without  set -e , as the  set -e  restoration would change the caller's
  # error-handling semantics.
  #
  # PIPESTATUS note: uses set +e / set -e sandwich around the tee pipeline.
  # This is the only correct pattern under set -eo pipefail to capture PIPESTATUS[0].
  local _rs_ite="$1" _rs_step="$2" _rs_log_dir="$3"
  shift 3
  # $@ is now the command + all its arguments

  mkdir -p "${_rs_log_dir}"
  local _rs_log="${_rs_log_dir}/${_rs_step}.log"

  local _rs_cmd_display
  printf -v _rs_cmd_display '%q ' "$@"

  local _rs_t_start _rs_ep_start _rs_t_end _rs_ep_end _rs_elapsed _rs_rc=0
  _rs_t_start=$(date -u +%Y-%m-%dT%H:%M:%SZ)
  _rs_ep_start=$(date +%s)

  _EV_TYPE="step_start" _EV_ITE="${_rs_ite}" _EV_STEP="${_rs_step}" \
    _EV_STATUS="running" _EV_LOG="${_rs_log}" \
    _EV_CMD="${_rs_cmd_display}" _EV_T_START="${_rs_t_start}" \
    _EV_TAG="${tag:-}" \
    emit_runtime_event
  write_runtime_status "${_rs_ite}" "${_rs_step}" "running" \
    "${_rs_log}" "${_rs_t_start}" ""

  printf '\n==> [%s] %s\n\n' "$(date +%T)" "${_rs_step}" | tee -a "${_rs_log}" || true

  # Temporarily disable set -e to capture PIPESTATUS correctly from the pipe.
  # set -e is restored immediately; return ${_rs_rc} then propagates failure.
  set +e
  "$@" 2>&1 | tee -a "${_rs_log}"
  _rs_rc=${PIPESTATUS[0]}
  set -e

  _rs_ep_end=$(date +%s)
  _rs_t_end=$(date -u +%Y-%m-%dT%H:%M:%SZ)
  _rs_elapsed=$(( _rs_ep_end - _rs_ep_start ))

  local _rs_status="success"
  [ "${_rs_rc}" -ne 0 ] && _rs_status="failed"

  _EV_TYPE="step_end" _EV_ITE="${_rs_ite}" _EV_STEP="${_rs_step}" \
    _EV_STATUS="${_rs_status}" _EV_RC="${_rs_rc}" \
    _EV_LOG="${_rs_log}" _EV_CMD="${_rs_cmd_display}" \
    _EV_T_START="${_rs_t_start}" _EV_T_END="${_rs_t_end}" \
    _EV_ELAPSED="${_rs_elapsed}" _EV_TAG="${tag:-}" \
    emit_runtime_event
  write_runtime_status "${_rs_ite}" "${_rs_step}" "${_rs_status}" \
    "${_rs_log}" "${_rs_t_start}" "${_rs_t_end}"

  printf '%s,%s,%s,%s,%d,%d\n' \
    "${_rs_ite}" "${_rs_step}" "${_rs_t_start}" "${_rs_t_end}" \
    "${_rs_elapsed}" "${_rs_rc}" \
    >> "${workingDir}/runtime/step_timings.csv" 2>/dev/null || true

  printf '\n==> [%s] %s done  (%ds, rc=%d)\n\n' \
    "$(date +%T)" "${_rs_step}" "${_rs_elapsed}" "${_rs_rc}" \
    | tee -a "${_rs_log}" || true

  return "${_rs_rc}"
}

finish_runtime_logging() {
  # Mark done FIRST so the EXIT trap (installed by init_runtime_logging) does not
  # emit a duplicate session_end event when the script exits normally after this call.
  _JANAS_SESSION_DONE=1
  local _ft; _ft=$(date -u +%Y-%m-%dT%H:%M:%SZ)
  _EV_TYPE="session_end" _EV_ITE="--" _EV_STEP="session_end" \
    _EV_STATUS="finished" _EV_T_END="${_ft}" \
    emit_runtime_event
  write_runtime_status "--" "session_end" "finished" "" "" "${_ft}"
  echo "[janas-log] Session finished: ${workingDir}/runtime/"
}

"""


#######################################################
#######################################################
#####      No external code
def generate_run_script_no_external(fileSettings, gpu_list=None, use_external_relion=False):
    with open(fileSettings, "r") as f:
        data = toml.load(f)

    def as_str(v, default=""):
        return str(v) if v is not None else str(default)

    # Build GPU argument only from CLI list, ignore TOML
    if gpu_list:
        gpu_arg = "--gpu " + " ".join(str(g) for g in gpu_list)
    else:
        gpu_arg = ""

    external_flag = "True" if use_external_relion else "False"

    # ------------------------------
    # Extract knobs
    # ------------------------------
    autoSigma_flag = as_str(data.get("autoSigma", "False"))
    bootstrap_flag = as_str(data.get("bootstrap", "False"))
    maxSel         = as_str(data.get("maxSelections", 8))
    numViews       = as_str(data.get("numViews", 350))
    numMpi         = as_str(data.get("mpi", 5))
    numRecs        = as_str(data.get("numRecs", 10))
    score_method   = as_str(data.get("score", "SCIFIR")).upper()
    subtractionFlag= as_str(data.get("subtractionFlag", "False"))
    aggressive_flag = as_str(data.get("aggressive_mode", "False"))



    # optimiser parameters (mostly passed by default)
    samplingDensityFactor = as_str(data.get("samplingDensityFactor", 0.5))
    extraSamples_num      = as_str(data.get("extraSamples_num", 5))
    extraSamples_random   = as_str(data.get("extraSamples_randomness", 0.05))
    extraSamples_audacity = as_str(data.get("extraSamples_audacity", 0.8))

    # ------------------------------
    # Build bash script (single large string)
    # ------------------------------
    run_script_cmd = """#!/bin/bash
set -eo pipefail

# Locale for numeric parsing
export LC_ALL=C
export LC_NUMERIC=C

readonly True=1
readonly False=0

# ------------------------------
# Utilities
# ------------------------------
str2bool() {
  case "${1,,}" in
    true|1)  echo "$True" ;;
    false|0) echo "$False" ;;
    *)       echo "$False" ;;
  esac
}

ensure_directory() {
  local dir="${1:?ensure_directory: missing path}"
  mkdir -p -- "$dir"
}

# ------------------------------
# Reconstruction / locres backends
# ------------------------------
rec_subsets() {
  local starFileIn="$1"
  local fileOut_basename="$2"

  rm -f "${fileOut_basename}_recH1.mrc" "${fileOut_basename}_recH2.mrc"

  local all_done=false
  while [ "$all_done" = false ]; do
    all_done=true
    for subset in 1 2; do
      local fileOut="${fileOut_basename}_recH${subset}.mrc"
      if [ -f "$fileOut" ]; then
        echo "DOING NOTHING: Reconstructed file $fileOut exists"
      else
        all_done=false
        echo "DOING Reconstruction for file $fileOut"
        relion_reconstruct --i "$starFileIn" --o "$fileOut" --subset "$subset" --ctf &
      fi
    done
    wait
  done
}

locres_one_relion() {
  local file_basename="$1"

  if [ -f "${file_basename}_locres.mrc" ]; then
    echo "DOING NOTHING: locres file ${file_basename}_locres.mrc exists"
  else
    echo "DOING locres for file ${file_basename}_locres.mrc"
    mpirun --np "${NUM_MPI}" relion_postprocess_mpi \
      --i "${file_basename}_recH1.mrc" \
      --i2 "${file_basename}_recH2.mrc" \
      --o "${file_basename}" \
      --locres
    rm -f "${file_basename}_locres_fscs.star" "${file_basename}_locres_filtered.mrc"
  fi
}

reconstruct_and_locres_external() {
  local scored_star="$1"
  local out_base="$2"
  local listParticles="$3"

  for numParticles in $(echo "$listParticles" | sed 's/,/ /g'); do
    local subset_base="${out_base}${numParticles}"
    janas selectBestRanked --i "${scored_star}" --o "${subset_base}.star" --num "${numParticles}"
    rec_subsets "${subset_base}.star" "${subset_base}"
    locres_one_relion "${subset_base}"
  done
}

# Normalise to two decimals; accept only numeric; fallback to 1.00
format2d() {
  local v="${1:-}"
  v="${v/,/.}"
  if [[ "$v" =~ ^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$ ]]; then
    printf '%.2f' "$v"
  else
    printf '%.2f' 1.00
  fi
}

# ------------------------------
# Auto-sigma helpers
# ------------------------------
is_nonempty_file() {
  local f="$1"
  [[ -n "$f" && "$f" != "None" && -f "$f" ]]
}

compute_autosigma() {
  local h1="$1"
  local h2="$2"
  local m="$3"

  if ! is_nonempty_file "$h1" || ! is_nonempty_file "$h2"; then
    echo ""
    return 0
  fi

  local cmd=(janas_utils sigma_estimate "$h1" "$h2")
  if [[ -n "$m" && "$m" != "None" && "${m,,}" != "none" ]]; then
    cmd+=("$m")
  fi

  local out
  # sigma_estimate can print verbose logs; extract the last line that is purely numeric
  if ! out="$("${cmd[@]}" 2>/dev/null | awk '
      /^[[:space:]]*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?[[:space:]]*$/ {v=$0}
      END{
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", v);
        print v
      }'
    )"; then
    echo ""
    return 0
  fi

  out="${out/,/.}"
  if [[ "$out" =~ ^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$ ]]; then
    printf '%s' "$out"
  else
    echo ""
  fi
}


# Refresh overview-derived globals used by the loop
update_global_variables() {
	local workingDir=$1
	DIRS=$(janas_optimizer generate_overview --directory "${workingDir}")

	if [[ "$DIRS" == "ERROR_EMPTY" ]]; then
		echo "[[_janas_target_selection]]" > "${workingDir}/overview.txt"
		echo "[[_janas_selection_0]]" >> "${workingDir}/overview.txt"
		notImprovingIterations=0
		next_sigma=${sigma}
	else
		previousLocresFile=$(janas_optimizer getTarget --overviewFile "${workingDir}/overview.txt" --locres_file)
		sigma=$(janas_optimizer getTarget --overviewFile "${workingDir}/overview.txt" --sigma)
		notImprovingIterations=$(janas_optimizer getTarget --overviewFile "${workingDir}/overview.txt" --unimproved)
		next_sigma=$(janas_optimizer settingBasedOptimization next_sigma --overviewFile "${workingDir}/overview.txt" --settingsFile "${workingDir}/session_settings.toml" --lastNonImproving "${notImprovingIterations}")
		targetMap1=$(janas_optimizer getTarget --overviewFile "${workingDir}/overview.txt" --map1)
		targetMap2=$(janas_optimizer getTarget --overviewFile "${workingDir}/overview.txt" --map2)
        
	fi
}


#main iteration
process_iteration() {
    local ite=$1
    local workingDir=$2
    local maskedCropOption=$3

    local tag="_janas_SCI__${sigma}_scored_selection_${ite}"
    ensure_directory ${workingDir}
    ensure_directory "${workingDir}/${tag}"
    echo "#############" > ${workingDir}/${tag}/log.txt
    echo "ite=${ite}" >> ${workingDir}/${tag}/log.txt
    echo "directory=${tag}" >> ${workingDir}/${tag}/log.txt
    signalSubtractionFlag=""" + subtractionFlag + """
    DO_WORK_WITH_SIGNAL_SUBTRACTION=$signalSubtractionFlag  # or $False
    local fileWithNormalization=${workingDir}/${tag}/scored_selection_${ite}.star #for regular operations
"""
    

    if data.get("singleMapReferenceFlag", "True") == "True":
        run_script_cmd += (
            "    local _map_for_score=\"${singleReferenceFile}\"\n"
            "    if awk 'BEGIN{exit !(s>0)}' s=\"${preGaussianBlur}\"; then\n"
            "        if [ -n \"${_map_for_score}\" ] && [ \"${_map_for_score}\" != \"None\" ] && [ -f \"${_map_for_score}\" ]; then\n"
            "            local _blur1=\"${workingDir}/${tag}/target_map1_preblur.mrc\"\n"
            "            janas_utils clip blur \"${_map_for_score}\" \"${_blur1}\" \"${preGaussianBlur}\"\n"
            "            _map_for_score=\"${_blur1}\"\n"
            "        fi\n"
            "    fi\n"
            "    run_step \"${ite}\" \"03_score_particles\" \"${workingDir}/${tag}/steps\" \\\n"
            "        janas scoreParticles --i $targetStar --map \"${_map_for_score}\"  --mask $targetMask --apix ${apix} --sigma ${sigma} --selectionName ${ite} --o ${fileWithNormalization} --mpi "
            + numMpi
            + " --rank "
            + numViews
            + "\n"
        )
    else:
        if data.get("bootstrap", "True") == "True":
            run_script_cmd += "\n    #bootstrap, randomized particles\n"
            run_script_cmd += "    local particlesStar=\"${targetStar}\"\n"
            run_script_cmd += "    if [ -f \"${workingDir}/overview.txt\" ]; then\n"
            run_script_cmd += "        particlesStar=$(janas_optimizer getTarget --overviewFile \"${workingDir}/overview.txt\" --particles)\n"
            run_script_cmd += "    fi\n"
            run_script_cmd += (
                '    run_step "${ite}" "01_randomize_halves" "${workingDir}/${tag}/steps" \\\n'
                '        janas_utils randomize_halves --i "${particlesStar}" --o "${workingDir}/target_particles_randomizedHalf.star"\n'
            )

            run_script_cmd += (
                '    if [ "$(str2bool "$USE_EXTERNAL_TOOLS")" -eq "$True" ]; then\n'
                '        run_step "${ite}" "02_bootstrap_reconstruct" "${workingDir}/${tag}/steps" \\\n'
                '            rec_subsets "${workingDir}/target_particles_randomizedHalf.star" "${workingDir}/target_particles_randomizedHalf"\n'
                '    else\n'
                '        run_step "${ite}" "02_bootstrap_reconstruct" "${workingDir}/${tag}/steps" \\\n'
                '            janas_reconstructor "${workingDir}/target_particles_randomizedHalf.star" --subset 1 2 --ctf ${GPU_ARG_OPT} --out-basename "${workingDir}/target_particles_randomizedHalf"\n'
                '    fi\n'
            )
            run_script_cmd += (
                "    local _map1_for_score=\"${workingDir}/target_particles_randomizedHalf_recH1.mrc\"\n"
                "    local _map2_for_score=\"${workingDir}/target_particles_randomizedHalf_recH2.mrc\"\n"
                "    if awk 'BEGIN{exit !(s>0)}' s=\"${preGaussianBlur}\"; then\n"
                "        if [ -f \"${_map1_for_score}\" ]; then\n"
                "            local _b1=\"${workingDir}/${tag}/target_map1_preblur.mrc\"\n"
                "            janas_utils clip blur \"${_map1_for_score}\" \"${_b1}\" \"${preGaussianBlur}\"\n"
                "            _map1_for_score=\"${_b1}\"\n"
                "        fi\n"
                "        if [ -f \"${_map2_for_score}\" ]; then\n"
                "            local _b2=\"${workingDir}/${tag}/target_map2_preblur.mrc\"\n"
                "            janas_utils clip blur \"${_map2_for_score}\" \"${_b2}\" \"${preGaussianBlur}\"\n"
                "            _map2_for_score=\"${_b2}\"\n"
                "        fi\n"
                "    fi\n"
                "    run_step \"${ite}\" \"03_score_particles\" \"${workingDir}/${tag}/steps\" \\\n"
                "        janas scoreParticles --i ${workingDir}/target_particles_randomizedHalf.star --map \"${_map1_for_score}\"  --map2 \"${_map2_for_score}\"  --mask $targetMask --apix ${apix} --sigma ${sigma} --selectionName ${ite} --o ${fileWithNormalization} --mpi "
                + numMpi
                + " --rank "
                + numViews
                + "\n"
            )

        else:
            run_script_cmd += (
                "    local _map1_for_score=\"${targetMap1}\"\n"
                "    local _map2_for_score=\"${targetMap2}\"\n"
                "    if awk 'BEGIN{exit !(s>0)}' s=\"${preGaussianBlur}\"; then\n"
                "        if [ -n \"${_map1_for_score}\" ] && [ \"${_map1_for_score}\" != \"None\" ] && [ -f \"${_map1_for_score}\" ]; then\n"
                "            local _b1=\"${workingDir}/${tag}/target_map1_preblur.mrc\"\n"
                "            janas_utils clip blur \"${_map1_for_score}\" \"${_b1}\" \"${preGaussianBlur}\"\n"
                "            _map1_for_score=\"${_b1}\"\n"
                "        fi\n"
                "        if [ -n \"${_map2_for_score}\" ] && [ \"${_map2_for_score}\" != \"None\" ] && [ -f \"${_map2_for_score}\" ]; then\n"
                "            local _b2=\"${workingDir}/${tag}/target_map2_preblur.mrc\"\n"
                "            janas_utils clip blur \"${_map2_for_score}\" \"${_b2}\" \"${preGaussianBlur}\"\n"
                "            _map2_for_score=\"${_b2}\"\n"
                "        fi\n"
                "    fi\n"
                "    run_step \"${ite}\" \"03_score_particles\" \"${workingDir}/${tag}/steps\" \\\n"
                "        janas scoreParticles --i $targetStar --map \"${_map1_for_score}\"  --map2 \"${_map2_for_score}\"  --mask $targetMask --apix ${apix} --sigma ${sigma} --selectionName ${ite} --o ${fileWithNormalization} --mpi "
                + numMpi
                + " --rank "
                + numViews
                + "\n"
            )

    run_script_cmd += "    \n"

    run_script_cmd += (
        """    if [ $DO_WORK_WITH_SIGNAL_SUBTRACTION -eq $True ]; then
    	echo "Working with signal subtraction"
    	janas_app_starProcess --i ${fileWithNormalization}  --invertTagName _janas_backup_rlnImageName _rlnImageName --o ${workingDir}/${tag}/tmp_inverted.star
    	local fileWithNormalization=${workingDir}/${tag}/tmp_inverted.star #for signal subtraction, replace the file with subtraction with the whole particles
    else
    	echo "Working without signal subtraction (regular)"
    fi


    #produce reconstruction script for selected subset of particles
    run_step "${ite}" "04_particle_subsets" "${workingDir}/${tag}/steps" \
        janas_optimizer automaticParticleSubsets """ + (" --aggressive " if aggressive_flag else " ") + """ --starFile ${fileWithNormalization} --locres ${previousLocresFile} --numSamples """ + numRecs
        + """ --samplingDensityFactor """ + samplingDensityFactor
        + """  --extraSamples_num  """ + extraSamples_num
        + """  --extraSamples_randomness   """ + extraSamples_random
        + """  --extraSamples_audacity   """ + extraSamples_audacity
        + """    --save ${workingDir}/${tag}/reconstructions_list.csv

    listParticles=$(<${workingDir}/${tag}/reconstructions_list.csv)
    echo "listParticles=$listParticles" >>${workingDir}/${tag}/log.txt

    if [ "$(str2bool "$USE_EXTERNAL_TOOLS")" -eq "$True" ]; then
        run_step "${ite}" "05_subset_reconstruct" "${workingDir}/${tag}/steps" \
            reconstruct_and_locres_external "${fileWithNormalization}" "${workingDir}/${tag}/norm_${tag}_best" "$listParticles"
    else
        run_step "${ite}" "05_subset_reconstruct" "${workingDir}/${tag}/steps" \
            janas_reconstructor ${fileWithNormalization} --subset 1 2 --subrec-only $listParticles --ctf --out-basename ${workingDir}/${tag}/norm_${tag}_best --sort ${tag}_norm""" + numViews + """ ${GPU_ARG_OPT}
        run_step "${ite}" "06_locres_bulk" "${workingDir}/${tag}/steps" \
            janas_utils locresBulk --threshold 0.143 --gamma 1.5 --cycles 15 --cpu """ + numMpi + """ --mask ${targetMask} ${workingDir}/${tag}/norm_${tag}_best $listParticles
    fi

    mask_for_locresStats="$targetMask"
    if [ -n "${assessMask}" ] && [ "${assessMask}" != "None" ]; then
        mask_for_locresStats="${assessMask}"
    fi
    run_step "${ite}" "07_locres_stats" "${workingDir}/${tag}/steps" \
        janas_utils locresStats --locres-files "${workingDir}/${tag}/norm_${tag}_best*locres.mrc" --mask ${mask_for_locresStats}  --assessmentMethod ${assessmentMethod} --out-stats ${workingDir}/${tag}/bestRanked_locres_values.csv --out-local-best-particles ${workingDir}/${tag}/partialLocres.mrc

    #save the target
    run_step "${ite}" "08_get_num_particles" "${workingDir}/${tag}/steps" \
        janas_optimizer getNumParticles --locres ${workingDir}/${tag}/bestRanked_locres_values.csv --save ${workingDir}/${tag}/target_num_of_particles.csv   ${NUMPART_METHOD_FLAG} --saveSplineOnCsv ${workingDir}/${tag}/spline_locres_predictions.csv

    targetNumOfParticles=$(<${workingDir}/${tag}/target_num_of_particles.csv)
    previousLocresFile=${workingDir}/${tag}/bestRanked_locres_values.csv

    janas_optimizer generate_overview --directory ${workingDir}
"""
    )


    run_script_cmd += "    targetMap1=$(janas_optimizer getTarget --overviewFile ${workingDir}/overview.txt --map1)\n"
    run_script_cmd += "    targetMap2=$(janas_optimizer getTarget --overviewFile ${workingDir}/overview.txt --map2)\n"
    if aggressive_flag == "True":
        run_script_cmd += "    targetStar=$(janas_optimizer getTarget --overviewFile ${workingDir}/overview.txt --particles)\n"
        
    if data.get("singleMapReferenceFlag", "True") == "True":
        run_script_cmd += (
            "    local singleReferenceFileBasename=${workingDir}/${tag}/norm_${tag}\n"
        )
        if data.get("postprocess", "avg") == "autobfac":
            run_script_cmd += "    relion_postprocess --i $targetMap1 --i2 $targetMap2  --auto_bfac --o ${singleReferenceFileBasename}_autobfac \n"
            run_script_cmd += (
                "    singleReferenceFile=${singleReferenceFileBasename}_autobfac.mrc \n"
            )
        else:
            run_script_cmd += "    janas_utils clip average $targetMap1  $targetMap2 ${singleReferenceFileBasename}.mrc \n"
            run_script_cmd += (
                "    singleReferenceFile=${singleReferenceFileBasename}.mrc \n"
            )

    run_script_cmd += """
    

    local unimprovingIterations=$(janas_optimizer getTarget --overviewFile ${workingDir}/overview.txt --unimproved)
    useSubsetParticleCheckStr=$( janas_optimizer settingBasedOptimization subset_particle_check  --settingsFile "${workingDir}/session_settings.toml" --current_unimproving_iterations "${notImprovingIterations}")
    useSubsetParticleCheck=$(str2bool "$useSubsetParticleCheckStr")
    if [ "$useSubsetParticleCheck" -eq "$True" ]; then
  	targetStar=$(janas_optimizer getTarget  --overviewFile "${workingDir}/overview.txt"  --particles )
    fi


"""

    run_script_cmd += "}\n\n\n"
    # Static config values (one per line)
    run_script_cmd += f'sigma="{as_str(data.get("sigma", "1.00"))}"\n'
    run_script_cmd += f'preGaussianBlur="{as_str(data.get("preGaussianBlur", "0.0"))}"\n'
    run_script_cmd += f'targetStar="{as_str(data.get("particles", "None"))}"\n'
    run_script_cmd += f'targetMap1="{as_str(data.get("map", "None"))}"\n'
    run_script_cmd += f'targetMap2="{as_str(data.get("map2", "None"))}"\n'
    run_script_cmd += f'aggressive_mode="{as_str(data.get("aggressive_mode", "False"))}"\n'
    if as_str(data.get("singleMapReferenceFlag", "False")).lower() == "true":
        run_script_cmd += f'singleReferenceFile="{as_str(data.get("map", "None"))}"\n'
    else:
        run_script_cmd += 'singleReferenceFile=""\n'
    run_script_cmd += f'autoSigma="{autoSigma_flag}"\n'
    run_script_cmd += f'autoSigmaInitialHalfMaps="{as_str(data.get("autoSigmaInitialHalfMaps", ""))}"\n'
    run_script_cmd += f'autoSigmaMask="{as_str(data.get("autoSigmaMask", ""))}"\n'
    run_script_cmd += f'bootstrap="{bootstrap_flag}"\n'
    run_script_cmd += f'targetMask="{as_str(data.get("mask", "None"))}"\n'
    run_script_cmd += f'assessMask="{as_str(data.get("assessMask", ""))}"\n'
    run_script_cmd += 'assessmentMethod="' + str(data.get("assessmentMethod", "mean")).lower() + '"\n'
    # NUMPART_METHOD_FLAG depends on assessmentMethod and must be set after it.
    run_script_cmd += 'NUMPART_METHOD_FLAG="--mean_res"\n'
    run_script_cmd += 'if [ "${assessmentMethod}" = "median" ]; then\n'
    run_script_cmd += '    NUMPART_METHOD_FLAG="--median_res"\n'
    run_script_cmd += 'fi\n'
    run_script_cmd += f'AssessCtfMode="{as_str(data.get("ctf_mode", "image"))}"\n'
    run_script_cmd += f'do_subtraction="{as_str(data.get("do_subtraction", "False"))}"\n'
    run_script_cmd += f'subtraction_mask="{as_str(data.get("subtraction_mask", ""))}"\n'
    run_script_cmd += f'workingDir="{as_str(data.get("session_name"))}"\n'
    run_script_cmd += f'SCORE_METHOD="{score_method}"\n'
    run_script_cmd += f'NUM_VIEWS="{numViews}"\n'
    run_script_cmd += f'NUM_MPI="{numMpi}"\n'
    run_script_cmd += f'NUM_RECS="{numRecs}"\n'
    run_script_cmd += f'GPU_ARG_OPT="{gpu_arg}"\n'
    run_script_cmd += f'USE_EXTERNAL_TOOLS="{external_flag}"\n'
    run_script_cmd += f'apix="{as_str(data.get("apix", "None"))}"\n'
    run_script_cmd += 'previousLocresFile="None"\n'
    run_script_cmd += "\n"

    # Main loop
    run_script_cmd += f"""
# ==============================
# Main loop
# ==============================
ensure_directory "${{workingDir}}"

update_global_variables $workingDir

# Auto-sigma initialisation (optional):
# Prefer --autoSigmaInitialHalfMaps if provided; otherwise use session --map/--map2.
if [ "$(str2bool "$autoSigma")" -eq "$True" ]; then
  est_h1=""
  est_h2=""
  if [[ -n "$autoSigmaInitialHalfMaps" ]]; then
    read -r est_h1 est_h2 <<< "$autoSigmaInitialHalfMaps"
  fi

  if ! is_nonempty_file "$est_h1" || ! is_nonempty_file "$est_h2"; then
    est_h1="$targetMap1"
    est_h2="$targetMap2"
  fi

  auto_sig="$(compute_autosigma "$est_h1" "$est_h2" "$autoSigmaMask")"
  if [[ -n "$auto_sig" ]]; then
    sigma="$auto_sig"
    next_sigma="$auto_sig"
  fi
fi


sigma="$(format2d "${{sigma}}")"
if [ -z "${{next_sigma:-}}" ] || [ "${{next_sigma}}" = "None" ]; then
  next_sigma="${{sigma}}"
fi
next_sigma="$(format2d "${{next_sigma}}")"

current_selection_ID=$(janas_optimizer getTarget --overviewFile "${{workingDir}}/overview.txt" --current_selection_ID)
next_selection_ID=$(( current_selection_ID + 1 ))

for (( ite=next_selection_ID; ite<={maxSel}; ite++ )); do
  sigma="$(format2d "${{next_sigma}}")"
  process_iteration "$ite" "$workingDir" ""
  update_global_variables "$workingDir"

  # Auto-sigma update from current best half-maps (only if a valid pair exists)
  if [ "$(str2bool "$autoSigma")" -eq "$True" ]; then
    if is_nonempty_file "$targetMap1" && is_nonempty_file "$targetMap2"; then
      auto_sig="$(compute_autosigma "$targetMap1" "$targetMap2" "$autoSigmaMask")"
      if [[ -n "$auto_sig" ]]; then
        next_sigma="$auto_sig"
      fi
    fi
  fi
  
  next_sigma="$(format2d "${{next_sigma}}")"
  tmp_unimprovingIterations=$(janas_optimizer getTarget --overviewFile "${{workingDir}}/overview.txt" --unimproved)
  earlyTerminationStr=$(janas_optimizer settingBasedOptimization check_early_termination --settingsFile "${{workingDir}}/session_settings.toml" --current_unimproving_iterations "${{tmp_unimprovingIterations}}")
  earlyTermination=$(str2bool "${{earlyTerminationStr}}")
  if [ "${{earlyTermination}}" -eq "${{True}}" ]; then
    echo "Early termination after iteration $ite"
    break
  fi
done
"""

    # ------------------------------
    # Write script to disk
    # ------------------------------
    runScriptName = os.path.join(
        as_str(data.get("dir", "./")),
        as_str(data.get("session_name")),
        as_str(data.get("session_name")) + "_run.sh",
    )

    # Inject numeric knobs while leaving bash ${...} untouched
    run_script_cmd = (
        run_script_cmd
        .replace('{samplingDensityFactor}', str(samplingDensityFactor))
        .replace('{extraSamples_num}',      str(extraSamples_num))
        .replace('{extraSamples_random}',   str(extraSamples_random))
        .replace('{extraSamples_audacity}', str(extraSamples_audacity))
    )

    # ------------------------------------------
    # Inject runtime logging helpers (before update_global_variables)
    # ------------------------------------------
    run_script_cmd = _assert_replace_once(
        run_script_cmd,
        "# Refresh overview-derived globals used by the loop\n",
        "# ------------------------------\n"
        "# Runtime logging helpers\n"
        "# ------------------------------\n"
        + _runtime_logging_shell_block()
        + "# Refresh overview-derived globals used by the loop\n",
        "runtime-helpers-inject",
    )

    # Add init_runtime_logging call right after ensure_directory in the main loop
    run_script_cmd = _assert_replace_once(
        run_script_cmd,
        'ensure_directory "${workingDir}"\n\nupdate_global_variables $workingDir\n',
        'ensure_directory "${workingDir}"\ninit_runtime_logging\n\nupdate_global_variables $workingDir\n',
        "init-runtime-logging-inject",
    )

    # Add finish_runtime_logging after the main iteration loop
    run_script_cmd = _assert_replace_once(
        run_script_cmd,
        '    echo "Early termination after iteration $ite"\n    break\n  fi\ndone\n',
        '    echo "Early termination after iteration $ite"\n    break\n  fi\ndone\nfinish_runtime_logging\n',
        "finish-runtime-logging-inject",
    )

    # ------------------------------------------
    # Optional: inject CryoSPARC Local NU-refinement into the run script
    # ------------------------------------------
    run_cs = as_str(data.get("run_csparc_localnu", "false")).lower() == "true"
    if run_cs:
        cs_project = as_str(data.get("cs_project", ""))
        cs_workspace = as_str(data.get("cs_workspace", ""))
        ln_mask = as_str(data.get("localNURefinementMask", ""))
        out_basename = as_str(data.get("cs_localnu_out_basename", "reference_LNU"))

        # 1) Source CryoSPARC environment
        run_script_cmd = _assert_replace_once(
            run_script_cmd,
            "export LC_NUMERIC=C\n\n",
            "export LC_NUMERIC=C\n\nsource ~/.janas/cryosparc_env.sh\n\n",
            "cs-source-env-inject",
        )

        # 2) Insert CryoSPARC settings block after readonly constants
        cryo_settings = f'''
# ==============================
# Optional: CryoSPARC Local NU-refinement
# ==============================
RUN_CSPARC_LOCALNU="True"

CSPARC_PROJECT="{cs_project}"
CSPARC_WORKSPACE="{cs_workspace}"
CSPARC_LANE="default"
CSPARC_SYM="C1"

# Optional args
CSPARC_PARTICLE_DIR="."
CSPARC_RESPLIT="True"
CSPARC_MIN_ANGULAR_STEP="0.2"
CSPARC_USER=""
CSPARC_LOGLEVEL="info"
CSPARC_PRECOMPUTED=""

# Provided by session manager (from session_settings.toml):
localNURefinementMask="{ln_mask}"
CSPARC_OUT_BASENAME="{out_basename}"

# Optimisation:
CSPARC_SKIP_IF_NO_IMPROVEMENT="True"
'''

        run_script_cmd = _assert_replace_once(
            run_script_cmd,
            "readonly False=0\n\n",
            "readonly False=0\n" + cryo_settings + "\n",
            "cs-settings-inject",
        )

        # 2b) Harden update_global_variables(): never clobber targets with "None"
        # Anchors use the EXACT indentation (tabs) from the triple-quoted bash function body.
        run_script_cmd = _assert_replace_once(
            run_script_cmd,
            '\t\ttargetMap1=$(janas_optimizer getTarget --overviewFile "${workingDir}/overview.txt" --map1)\n',
            '\t\tlocal __prev_targetMap1="${targetMap1}"\n'
            '\t\ttargetMap1=$(janas_optimizer getTarget --overviewFile "${workingDir}/overview.txt" --map1)\n'
            '\t\tif [ -z "${targetMap1}" ] || [ "${targetMap1}" = "None" ]; then targetMap1="${__prev_targetMap1}"; fi\n',
            "cs-harden-targetmap1",
        )
        run_script_cmd = _assert_replace_once(
            run_script_cmd,
            '\t\ttargetMap2=$(janas_optimizer getTarget --overviewFile "${workingDir}/overview.txt" --map2)\n',
            '\t\tlocal __prev_targetMap2="${targetMap2}"\n'
            '\t\ttargetMap2=$(janas_optimizer getTarget --overviewFile "${workingDir}/overview.txt" --map2)\n'
            '\t\tif [ -z "${targetMap2}" ] || [ "${targetMap2}" = "None" ]; then targetMap2="${__prev_targetMap2}"; fi\n',
            "cs-harden-targetmap2",
        )
        # targetStar is updated inside process_iteration() when useSubsetParticleCheck fires.
        # The anchor uses "  \t" (2 spaces + 1 tab) to match the literal indentation, plus the
        # double-spaces and trailing space that are present in the raw script string.
        run_script_cmd = _assert_replace_once(
            run_script_cmd,
            '  \ttargetStar=$(janas_optimizer getTarget  --overviewFile "${workingDir}/overview.txt"  --particles )\n',
            '  \tlocal __prev_targetStar="${targetStar}"\n'
            '  \ttargetStar=$(janas_optimizer getTarget  --overviewFile "${workingDir}/overview.txt"  --particles )\n'
            '  \tif [ -z "${targetStar}" ] || [ "${targetStar}" = "None" ]; then targetStar="${__prev_targetStar}"; fi\n',
            "cs-harden-targetstar",
        )



        # 3) Inject function definition before process_iteration()
        cryo_func = r'''
run_csparc_localnu_if_enabled() {
  local ite="${1:?missing ite}"
  local tagDir="${2:?missing tagDir}"

  local do_localnu
  do_localnu=$(str2bool "${RUN_CSPARC_LOCALNU}")
  if [ "${do_localnu}" -ne "${True}" ]; then
    return 0
  fi


  # Optional optimisation: skip local refinement when the previous selection did not improve.
  # This avoids re-running local refinement on the same particle subset.
  local session_dir
  session_dir="$(dirname -- "${tagDir}")"
  local overview_file="${session_dir}/overview.txt"
  if [ "${ite}" -gt 1 ] && [ -f "${overview_file}" ] && [ "$(str2bool "${CSPARC_SKIP_IF_NO_IMPROVEMENT}")" -eq "${True}" ]; then
    local unimproved
    unimproved="$(janas_optimizer getTarget --overviewFile "${overview_file}" --unimproved 2>/dev/null || true)"
    if [ -n "${unimproved}" ] && [ "${unimproved}" -gt 0 ] 2>/dev/null; then
      {
        echo "Skipping CryoSPARC local NU-refinement (ite=${ite}): no improvement in previous selection (unimproved=${unimproved})"
      } >> "${tagDir}/log.txt"
      local _cs_skip_ts; _cs_skip_ts=$(date -u +%Y-%m-%dT%H:%M:%SZ)
      _EV_TYPE="step_skipped" _EV_ITE="${ite}" _EV_STEP="00_csparc_lnu" \
        _EV_STATUS="skipped" _EV_T_START="${_cs_skip_ts}" \
        _EV_TAG="${tagDir##*/}" \
        emit_runtime_event
      write_runtime_status "${ite}" "00_csparc_lnu" "skipped" \
        "" "${_cs_skip_ts}" "${_cs_skip_ts}"
      printf '%s,%s,%s,%s,%d,%d\n' \
        "${ite}" "00_csparc_lnu" "${_cs_skip_ts}" "${_cs_skip_ts}" \
        0 0 \
        >> "${workingDir}/runtime/step_timings.csv" 2>/dev/null || true
      return 0
    fi
  fi

  # Runtime timing for this CryoSPARC LNU step.
  # Note: elapsed time measures shell-side wrapper (submission + synchronous wait),
  # not the actual cluster job compute time.
  local _cs_t_start _cs_ep_start
  _cs_t_start=$(date -u +%Y-%m-%dT%H:%M:%SZ)
  _cs_ep_start=$(date +%s)
  local _cs_log_dir="${tagDir}/steps"
  mkdir -p "${_cs_log_dir}"
  local _cs_step_log="${_cs_log_dir}/00_csparc_lnu.log"
  _EV_TYPE="step_start" _EV_ITE="${ite}" _EV_STEP="00_csparc_lnu" \
    _EV_STATUS="running" _EV_LOG="${_cs_step_log}" \
    _EV_T_START="${_cs_t_start}" _EV_TAG="${tagDir##*/}" \
    emit_runtime_event
  write_runtime_status "${ite}" "00_csparc_lnu" "running" \
    "${_cs_step_log}" "${_cs_t_start}" ""

  # Store outputs inside the iteration directory
  local out_base="${tagDir}/${CSPARC_OUT_BASENAME}"

  # Use current targetMap1 as reference volume.
  # If targetMap2 exists, average the two half-maps into a temporary reference map.
  local ref_map="${targetMap1}"
  if [ -n "${targetMap2:-}" ] && [ "${targetMap2}" != "None" ] && [ -f "${targetMap1}" ] && [ -f "${targetMap2}" ]; then
    local ref_avg="${tagDir}/${CSPARC_OUT_BASENAME}_refavg.mrc"
    janas_utils clip average "${targetMap1}" "${targetMap2}" "${ref_avg}"
    ref_map="${ref_avg}"
  fi





  # Use mask provided by the session settings
  local mask_map="${localNURefinementMask}"

  # Work on a private copy: csparc_localnurefinement sanitises STAR in-place.
  local star_in="${tagDir}/csparc_input_${ite}.star"
  cp -f -- "${targetStar}" "${star_in}"


  # Basic validation.
  # step_start was already emitted above, so every early-exit path MUST emit a
  # matching step_end to keep the event log balanced.  Use a single consolidated
  # exit point (_cs_val_rc pattern) to avoid duplicating the event-emission code.
  local _cs_val_rc=0
  if [ ! -f "${targetStar}" ]; then
    echo "ERROR: CryoSPARC localNU: targetStar not found: ${targetStar}" >&2
    _cs_val_rc=1
  elif [ ! -f "${ref_map}" ]; then
    echo "ERROR: CryoSPARC localNU: reference map not found: ${ref_map}" >&2
    _cs_val_rc=1
  elif [ ! -f "${mask_map}" ]; then
    echo "ERROR: CryoSPARC localNU: mask not found: ${mask_map}" >&2
    _cs_val_rc=1
  elif [ ! -f "${star_in}" ]; then
    echo "ERROR: CryoSPARC localNU: failed to create STAR copy: ${star_in}" >&2
    _cs_val_rc=1
  fi
  if [ "${_cs_val_rc}" -ne 0 ]; then
    local _cs_t_end_v; _cs_t_end_v=$(date -u +%Y-%m-%dT%H:%M:%SZ)
    local _cs_ep_end_v _cs_elapsed_v
    _cs_ep_end_v=$(date +%s)
    _cs_elapsed_v=$(( _cs_ep_end_v - _cs_ep_start ))
    _EV_TYPE="step_end" _EV_ITE="${ite}" _EV_STEP="00_csparc_lnu" \
      _EV_STATUS="failed" _EV_RC="${_cs_val_rc}" \
      _EV_T_START="${_cs_t_start}" _EV_T_END="${_cs_t_end_v}" \
      _EV_ELAPSED="${_cs_elapsed_v}" _EV_TAG="${tagDir##*/}" \
      emit_runtime_event
    write_runtime_status "${ite}" "00_csparc_lnu" "failed" \
      "${_cs_step_log}" "${_cs_t_start}" "${_cs_t_end_v}"
    printf '%s,%s,%s,%s,%d,%d\n' \
      "${ite}" "00_csparc_lnu" "${_cs_t_start}" "${_cs_t_end_v}" \
      "${_cs_elapsed_v}" "${_cs_val_rc}" \
      >> "${workingDir}/runtime/step_timings.csv" 2>/dev/null || true
    return "${_cs_val_rc}"
  fi

  # Optional args
  local opt_particle_dir=()
  [ -n "${CSPARC_PARTICLE_DIR:-}" ] && opt_particle_dir=(--particle-dir "${CSPARC_PARTICLE_DIR}")

  local opt_user=()
  [ -n "${CSPARC_USER:-}" ] && opt_user=(--user "${CSPARC_USER}")

  local opt_precomputed=()
  [ -n "${CSPARC_PRECOMPUTED:-}" ] && opt_precomputed=(--precomputed "${CSPARC_PRECOMPUTED}")

  local opt_resplit=()
  if [ "$(str2bool "${CSPARC_RESPLIT}")" -eq "${True}" ]; then
    opt_resplit=(--resplit)
  fi

  {
    echo "Running CryoSPARC local NU-refinement (ite=${ite})"
    echo "  input_star = ${targetStar}"
    echo "  input_star_copy = ${star_in}"
    echo "  ref_map    = ${ref_map}"
    echo "  mask       = ${mask_map}"
    echo "  out_base   = ${out_base}"
  } >> "${tagDir}/log.txt"

  local cs_log="${tagDir}/csparc_localnu_${ite}.log"
  rm -f "${cs_log}" 2>/dev/null || true
  echo "Running CryoSPARC local NU-refinement (ite=${ite})"
  janas_utils csparc_localnurefinement \
    "${opt_particle_dir[@]}" \
    --project "${CSPARC_PROJECT}" \
    --workspace "${CSPARC_WORKSPACE}" \
    --lane "${CSPARC_LANE}" \
    --sym "${CSPARC_SYM}" \
    --ref "${ref_map}" \
    --mask "${mask_map}" \
    "${opt_resplit[@]}" \
    --min-angular-step "${CSPARC_MIN_ANGULAR_STEP}" \
    "${opt_user[@]}" \
    --loglevel "${CSPARC_LOGLEVEL}" \
    "${opt_precomputed[@]}" \
    "${star_in}" \
    "${out_base}" 2>&1 | tee -a "${cs_log}" | tee -a "${tagDir}/log.txt"

  # PIPESTATUS[0] is the exit code of janas_utils csparc_localnurefinement.
  # Must be captured IMMEDIATELY after the pipeline — any intervening command resets PIPESTATUS.
  local rc=${PIPESTATUS[0]}

  # Copy CryoSPARC log into the steps per-step log
  cat "${cs_log}" >> "${_cs_step_log}" 2>/dev/null || true
  if [ "${rc}" -ne 0 ]; then
    local _cs_ep_end _cs_elapsed
    _cs_ep_end=$(date +%s)
    _cs_elapsed=$(( _cs_ep_end - _cs_ep_start ))
    local _cs_t_end; _cs_t_end=$(date -u +%Y-%m-%dT%H:%M:%SZ)
    _EV_TYPE="step_end" _EV_ITE="${ite}" _EV_STEP="00_csparc_lnu" \
      _EV_STATUS="failed" _EV_RC="${rc}" \
      _EV_T_START="${_cs_t_start}" _EV_T_END="${_cs_t_end}" \
      _EV_ELAPSED="${_cs_elapsed}" _EV_TAG="${tagDir##*/}" \
      emit_runtime_event
    write_runtime_status "${ite}" "00_csparc_lnu" "failed" \
      "${_cs_step_log}" "${_cs_t_start}" "${_cs_t_end}"
    printf '%s,%s,%s,%s,%d,%d\n' \
      "${ite}" "00_csparc_lnu" "${_cs_t_start}" "${_cs_t_end}" \
      "${_cs_elapsed}" "${rc}" \
      >> "${workingDir}/runtime/step_timings.csv" 2>/dev/null || true
    return "${rc}"
  fi

  # Save the CryoSPARC refinement job ID for tracking
  local job_id=""
  if [ -n "${CSPARC_PRECOMPUTED:-}" ]; then
    job_id="${CSPARC_PRECOMPUTED}"
  else
    job_id="$(grep -E 'Queuing job J[0-9]+ \(new_local_refine\)' "${cs_log}" | tail -n 1 | sed -E 's/.*Queuing job (J[0-9]+).*/\1/')"
    if [ -z "${job_id}" ]; then
      job_id="$(grep -Eo 'Queuing job J[0-9]+' "${cs_log}" | tail -n 1 | sed -E 's/Queuing job //')"
    fi
  fi
  if [ -n "${job_id}" ]; then
    printf '%s\n' "${job_id}" > "${tagDir}/cs_jobID.txt"
  fi


  # Switch targets to CryoSPARC outputs
  targetStar="${out_base}.star"
  targetMap1="${out_base}_recH1.mrc"
  targetMap2="${out_base}_recH2.mrc"

  {
    echo "CryoSPARC localNU outputs set as new targets:"
    echo "  targetStar = ${targetStar}"
    echo "  targetMap1 = ${targetMap1}"
    echo "  targetMap2 = ${targetMap2}"
  } >> "${tagDir}/log.txt"

  # Emit step_end event for 00_csparc_lnu (success path)
  local _cs_ep_end _cs_elapsed
  _cs_ep_end=$(date +%s)
  _cs_elapsed=$(( _cs_ep_end - _cs_ep_start ))
  local _cs_t_end; _cs_t_end=$(date -u +%Y-%m-%dT%H:%M:%SZ)
  _EV_TYPE="step_end" _EV_ITE="${ite}" _EV_STEP="00_csparc_lnu" \
    _EV_STATUS="success" _EV_RC="0" \
    _EV_T_START="${_cs_t_start}" _EV_T_END="${_cs_t_end}" \
    _EV_ELAPSED="${_cs_elapsed}" _EV_TAG="${tagDir##*/}" \
    emit_runtime_event
  write_runtime_status "${ite}" "00_csparc_lnu" "success" \
    "${_cs_step_log}" "${_cs_t_start}" "${_cs_t_end}"
  printf '%s,%s,%s,%s,%d,%d\n' \
    "${ite}" "00_csparc_lnu" "${_cs_t_start}" "${_cs_t_end}" \
    "${_cs_elapsed}" 0 \
    >> "${workingDir}/runtime/step_timings.csv" 2>/dev/null || true
}
'''

        run_script_cmd = _assert_replace_once(
            run_script_cmd,
            "process_iteration() {",
            cryo_func + "\n\nprocess_iteration() {",
            "cs-cryo-func-inject",
        )

        # 4) Inject call at the beginning of process_iteration()
        run_script_cmd = _assert_replace_once(
            run_script_cmd,
            'echo "directory=${tag}" >> ${workingDir}/${tag}/log.txt\n',
            'echo "directory=${tag}" >> ${workingDir}/${tag}/log.txt\n'
            '\n'
            '    # Aggressive mode: force targets to the current best selection from overview.txt\n'
            '    # Do NOT overwrite initial targets with "None" (overview may not contain a valid target yet).\n'
            '    if [ -f "${workingDir}/overview.txt" ] && [ "${aggressive_mode}" = "True" ]; then\n'
            '        local bestStar bestMap1 bestMap2\n'
            '        bestStar=$(janas_optimizer getTarget --overviewFile "${workingDir}/overview.txt" --particles)\n'
            '        bestMap1=$(janas_optimizer getTarget --overviewFile "${workingDir}/overview.txt" --map1)\n'
            '        bestMap2=$(janas_optimizer getTarget --overviewFile "${workingDir}/overview.txt" --map2)\n'
            '        if [ -n "${bestStar}" ] && [ "${bestStar}" != "None" ]; then targetStar="${bestStar}"; fi\n'
            '        if [ -n "${bestMap1}" ] && [ "${bestMap1}" != "None" ]; then targetMap1="${bestMap1}"; fi\n'
            '        if [ -n "${bestMap2}" ] && [ "${bestMap2}" != "None" ]; then targetMap2="${bestMap2}"; fi\n'
            '    fi\n'
            '\n'
            '    # Optional: run CryoSPARC local NU-refinement and store outputs in this iteration directory\n'
            '    run_csparc_localnu_if_enabled "${ite}" "${workingDir}/${tag}"\n',
            "cs-call-inject",
        )



    os.makedirs(os.path.dirname(runScriptName), exist_ok=True)
    with open(runScriptName, "w") as f:
        f.write(run_script_cmd)
    try:
        os.chmod(runScriptName, 0o755)
    except PermissionError:
        pass









#######################################################
## CREATE NEW SESSION, what the user has to call first
#######################################################
#######################################################
## CREATE NEW SESSION, what the user has to call first
#######################################################
janas_new_select_session = command.add_parser(
    "new_select_session",
    description="Create a new JANAS selection session.",
    help="Create a new selection session with experiment parameters",
    formatter_class=JANASHelpFormatter,
)

# -------------------------
# Core inputs
# -------------------------
core = janas_new_select_session.add_argument_group("Core inputs")
core.add_argument(
    "--name",
    required=True,
    type=str,
    help="Name of the new session. Creates a directory and a TOML settings file inside it.",
)
core.add_argument(
    "--particles",
    required=True,
    type=str,
    help="STAR file with the list of particles.",
)
core.add_argument(
    "--map",
    required=True,
    type=str,
    help="First reference, first half-map, or full map.",
)
core.add_argument(
    "--map2",
    required=False,
    type=str,
    help="Second half-map.",
)
core.add_argument(
    "--angpix",
    required=False,
    type=str,
    help="Pixel spacing in Angstrom. If omitted, it is inferred from --map.",
)

# -------------------------
# Masks
# -------------------------
masks = janas_new_select_session.add_argument_group("Masks")
masks.add_argument(
    "--mask",
    required=True,
    type=str,
    help="Main mask for the target region.",
)
masks.add_argument(
    "--subtractionMask",
    required=False,
    type=str,
    default=None,
    help=(
        "Optional 3D MRC mask defining the region to subtract during scoring "
        "(map-based signal subtraction)."
    ),
)
masks.add_argument(
    "--assessMask",
    required=False,
    type=str,
    default=None,
    help=(
        "Optional MRC mask used only for local-resolution assessment. "
        "If omitted, the main --mask is used."
    ),
)

# -------------------------
# Scoring and assessment
# -------------------------
scoring = janas_new_select_session.add_argument_group("Scoring and assessment")
scoring.add_argument(
    "--sigma",
    required=False,
    type=float,
    default=1.0,
    help="Sigma used for SCI scoring.",
)
scoring.add_argument(
    "--preGaussianBlur",
    required=False,
    type=float,
    default=0.0,
    help=(
        "Optional Gaussian blur applied to input map(s) before scoring "
        "(sigma in Angstrom). Use 0 to disable."
    ),
)
scoring.add_argument(
    "--sigma_decreasing_step",
    required=False,
    type=float,
    default=0.02,
    help="Decrease of sigma at each iteration.",
)
scoring.add_argument(
    "--CC",
    action="store_true",
    help="Use CC scoring instead of SCI scoring.",
)
scoring.add_argument(
    "--postprocessing",
    required=False,
    default="avg",
    type=str,
    help="Postprocessing mode: avg or autobfac.",
)
scoring.add_argument(
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
        f"Invalid choice: {v}. Choose from ['meanResolution', 'highResolution', "
        f"'highresolutionquartile', 'lowResolution', 'lowresolutionquartile']"
    ),
    default="meanResolution",
    help=(
        "Resolution target to optimise: meanResolution, highResolution, "
        "highresolutionquartile, lowResolution, or lowresolutionquartile."
    ),
)
scoring.add_argument(
    "--assessmentMethod",
    required=False,
    default="mean",
    choices=["mean", "median"],
    help="Statistic used for local-resolution assessment.",
)
scoring.add_argument(
    "--ctf-mode",
    required=False,
    default="phaseflip",
    choices=["none", "image", "phaseflip", "ref"],
    help=(
        "CTF handling for scoring: image, ref, phaseflip, or none."
    ),
)
scoring.add_argument(
    "--maskingCrop",
    action="store_true",
    help="Crop to the mask to accelerate local-resolution computation.",
)

# -------------------------
# Auto-sigma
# -------------------------
autosigma = janas_new_select_session.add_argument_group("Auto-sigma")
autosigma.add_argument(
    "--autoSigma",
    action="store_true",
    help=(
        "Estimate sigma automatically from two half-maps. "
        "Requires --map and --map2 with different files, or --bootstrap."
    ),
)
autosigma.add_argument(
    "--autoSigmaInitialHalfMaps",
    required=False,
    nargs=2,
    metavar=("HALFM1", "HALFM2"),
    type=str,
    help=(
        "Optional half-map pair used only for the initial auto-sigma estimate. "
        "If omitted, --map and --map2 are used."
    ),
)
autosigma.add_argument(
    "--autoSigmaMask",
    required=False,
    type=str,
    default="",
    help="Optional mask used only for auto-sigma estimation. Use 'none' for no mask.",
)

# -------------------------
# Iteration and optimisation
# -------------------------
optim = janas_new_select_session.add_argument_group("Iteration and optimisation")
optim.add_argument(
    "--bootstrap",
    action="store_true",
    help="Randomise half-maps at each iteration.",
)
optim.add_argument(
    "--aggressive",
    action="store_true",
    help="Aggressively update target particles from the current overview selection at each iteration.",
)
optim.add_argument(
    "--maxSelections",
    required=False,
    default=8,
    type=int,
    help="Maximum number of selections.",
)
optim.add_argument(
    "--numRecs",
    required=False,
    default=10,
    type=int,
    help="Number of sampling reconstructions per selection.",
)
optim.add_argument(
    "--samplingDensityFactor",
    required=False,
    type=float,
    default=0.5,
    help="Sampling density factor for reconstruction sampling.",
)
optim.add_argument(
    "--extraSamples_num",
    required=False,
    type=int,
    default=5,
    help="Number of extra reconstruction samples used to reduce local minima traps.",
)
optim.add_argument(
    "--extraSamples_randomness",
    required=False,
    type=float,
    default=0.05,
    help="Randomness of extra reconstruction samples, between 0 and 1.",
)
optim.add_argument(
    "--extraSamples_audacity",
    required=False,
    type=float,
    default=0.8,
    help="Audacity of extra reconstruction samples, between 0 and 1.",
)
optim.add_argument(
    "--randomSeed",
    action="store_true",
    help="Use deterministic random initialisation for reproducible sampling.",
)
optim.add_argument(
    "--numViews",
    required=False,
    default=350,
    type=str,
    help="Number of Euler views for ranking.",
)

# -------------------------
# Runtime and backend
# -------------------------
runtime = janas_new_select_session.add_argument_group("Runtime and backend")
runtime.add_argument(
    "--mpi",
    required=False,
    default=5,
    type=str,
    help="Number of MPI processes.",
)
runtime.add_argument(
    "--particleSubtraction",
    action="store_true",
    help=(
        "Input stack comes from particle subtraction. "
        "The STAR file should contain a backup field linking to unsubtracted particles."
    ),
)
runtime.add_argument(
    "--noExternalPrograms",
    action="store_true",
    help=(
        "Generate a self-contained run script that avoids external software packages "
        "(for example RELION for local-resolution estimation and reconstructions)."
    ),
)
runtime.add_argument(
    "--gpu",
    nargs="*",
    type=int,
    default=[],
    help="GPU indices to use, for example: --gpu 0 1 2. Leave empty for CPU-only.",
)



# -------------------------
# CryoSPARC integration
# -------------------------
cs = janas_new_select_session.add_argument_group("CryoSPARC integration")
cs.add_argument(
    "--cs_project",
    required=False,
    type=str,
    help="CryoSPARC project UID associated with this selection session, for example P31.",
)
cs.add_argument(
    "--cs_workspace",
    required=False,
    type=str,
    help="CryoSPARC workspace UID associated with this selection session, for example W1.",
)
cs.add_argument(
    "--multipleLocalRefineCS",
    nargs=4,
    metavar=("PROJECT", "WORKSPACE", "MASK_MRC", "OUT_BASENAME"),
    required=False,
    type=str,
    help=(
        "Enable iterative CryoSPARC Local NU-refinement in the generated *_run.sh. "
        "Arguments: <PROJECT> <WORKSPACE> <MASK_MRC> <OUT_BASENAME>."
    ),
)


# janas_new_select_session.add_argument("--comparisonType", required=False, default="halfmaps", type=str, help="type of selection, select from [halfmaps,singlemap,]")
# janas_new_select_session.add_argument("--typeSession", required=False, default="halfmaps", type=str, help="type of selection, select from [halfmaps,singlemap,]")


def new_select_session(args):

    haveReferenceHalfMaps = False
    fullMapReference = "None"
    if os.path.isfile(args.map):
        if args.map2:
            if os.path.isfile(args.map2):
                haveReferenceHalfMaps = True
    if haveReferenceHalfMaps:
        if not os.path.isfile(args.map2):
            print("ERROR: second half-map defined byt does not exist, exiting")
            exit(0)

    # Optional: user-provided half-map pair used only for the INITIAL auto-sigma estimate
    autoSigmaInitialHalfMaps = getattr(args, "autoSigmaInitialHalfMaps", None)
    if autoSigmaInitialHalfMaps is not None:
        if len(autoSigmaInitialHalfMaps) != 2:
            print("ERROR: --autoSigmaInitialHalfMaps requires exactly two files")
            exit(1)
        hm1, hm2 = autoSigmaInitialHalfMaps
        if not os.path.isfile(hm1) or not os.path.isfile(hm2):
            print("ERROR: --autoSigmaInitialHalfMaps files must exist")
            exit(1)
        if os.path.abspath(hm1) == os.path.abspath(hm2):
            print("ERROR: --autoSigmaInitialHalfMaps must be two different files")
            exit(1)
    # Optional: mask used ONLY for auto-sigma estimation (default: no mask)
    autoSigmaMask = getattr(args, "autoSigmaMask", "")
    if autoSigmaMask is None:
        autoSigmaMask = ""
    autoSigmaMask = str(autoSigmaMask).strip()
    if autoSigmaMask.lower() == "none":
        autoSigmaMask = ""
    if autoSigmaMask != "" and not os.path.isfile(autoSigmaMask):
        print("ERROR: --autoSigmaMask file must exist, or use 'none' / omit for no mask")
        exit(1)
    if getattr(args, "autoSigma", False):
        safe_autosigma = False
        if args.bootstrap:
            safe_autosigma = True
        else:
            # standard case: session provides a true half-map pair
            if args.map and args.map2 and os.path.isfile(args.map) and os.path.isfile(args.map2):
                if os.path.abspath(args.map) != os.path.abspath(args.map2):
                    safe_autosigma = True

            # allow initial sigma estimate from explicit initial halfmaps (even if session maps are not a pair)
            if not safe_autosigma and autoSigmaInitialHalfMaps is not None:
                safe_autosigma = True

        if not safe_autosigma:
            print(
                "ERROR: --autoSigma requires two independent half-maps (--map and --map2, different files), "
                "--bootstrap, or --autoSigmaInitialHalfMaps."
            )
            exit(1)



    haveReferenceFullMap = False
    if os.path.isfile(args.map) and not args.map2:
        haveReferenceFullMap = True
        fullMapReference = args.map

    if not haveReferenceFullMap and not haveReferenceHalfMaps:
        print("ERROR: you have define valid reference maps")
        exit(0)
    singleMapReferenceFlag = haveReferenceFullMap

    # Create the directory if it doesn't exist
    if not os.path.exists(args.name):
        os.makedirs(args.name)

    # Path for the settings file
    settings_file_path = os.path.join(args.name, "session_settings.toml")
    if not os.path.isfile(args.map):
        print("map file ", args.map, "not valid")

    if not os.path.isfile(args.mask):
        print("mask file ", args.mask, "not valid")

    # optional map-based subtraction mask
    do_subtraction_flag = "False"
    subtraction_mask_path = ""
    if getattr(args, "subtractionMask", None):
        if not os.path.isfile(args.subtractionMask):
            print("subtraction mask file ", args.subtractionMask, "not valid")
            exit(1)
        do_subtraction_flag = "True"
        subtraction_mask_path = args.subtractionMask

    assess_mask_path = ""
    if getattr(args, "assessMask", None):
        if not os.path.isfile(args.assessMask):
            print("assess mask file ", args.assessMask, "not valid")
            exit(1)
        assess_mask_path = args.assessMask

    if not args.angpix:
        raw_spacing = utils.get_MRC_map_pixel_spacing(args.map)
        try:
            pixel_size = raw_spacing[0]
        except (TypeError, IndexError):
            pixel_size = raw_spacing
        args.angpix = f"{pixel_size:.3f}"

    if args.particleSubtraction:
        particleSubtractionFlag = "True"
    else:
        particleSubtractionFlag = "False"

    if args.bootstrap:
        bootstrap = "True"
    else:
        bootstrap = "False"

    if args.postprocessing not in ["avg", "autobfac"]:
        args.postprocessing = "avg"

    # Create the settings file if it doesn't exist
    # if not os.path.isfile(settings_file_path):
    with open(settings_file_path, "w") as file:
        # Write a comment and the session_name variable
        # file.write("# name of the session\n")
        # toml.dump({'session_name': args.name}, file)
        file.write("###########################n\n")
        file.write("\n# name of the session\n")
        file.write(f'session_name = "{args.name}"\n')
        file.write("\n# particles stack\n")
        file.write(f'particles = "{args.particles}"\n')
        file.write("\n# Pixel size in angstrom for the map\n")
        file.write(f'apix = "{args.angpix}"\n')
        file.write("\n# Path of the mask file\n")
        file.write(f'mask= "{args.mask}"\n')

        file.write("\n# Optional map-based signal subtraction\n")
        file.write(f'do_subtraction = "{do_subtraction_flag}"\n')
        file.write(f'subtraction_mask = "{subtraction_mask_path}"\n')
        file.write("\n# Optional mask used only for locresStats (assessment)\n")
        file.write(f'assessMask = "{assess_mask_path}"\n')
        file.write("\n# locres assessment method for locresStats (mean or median)\n")
        file.write(f'assessmentMethod = "{str(args.assessmentMethod).lower()}"\n')

        file.write(
            "\n# maps. If only one map it is given, janas uses only this map, but this will break the assumption of independency for the two half maps, angular re-assignment might be required.\n"
        )
        file.write(f'map = "{args.map}"\n')
        file.write(f'map2 = "{args.map2}"\n')
        file.write(f'bootstrap = "{bootstrap}"\n')
        file.write(f'resolutionBestTarget = "{args.resolutionBestTarget}"\n')
        file.write(f'fullMap = "{fullMapReference}"\n')
        file.write(f'postprocess = "{args.postprocessing}"\n')
        file.write(f'singleMapReferenceFlag= "{singleMapReferenceFlag}"\n')

        file.write("\n# particles selection options\n")
        file.write('num_unimproving_iterations_before_selecting_from_best_subset = "1"\n')
        file.write('num_unimproving_iterations_for_early_termination = "3"\n')
        file.write(f'maxSelections = "{args.maxSelections}"\n')
        
        file.write("\n# sigma used for the SCI score\n")
        file.write(f'autoSigma = "{str(bool(args.autoSigma)).replace("True","True").replace("False","False")}"\n')
        file.write(f'autoSigmaInitialHalfMaps = "{(" ".join(autoSigmaInitialHalfMaps) if autoSigmaInitialHalfMaps else "")}"\n')
        file.write(f'autoSigmaMask = "{autoSigmaMask}"\n')
        file.write(f'sigma = "{args.sigma}"\n')
        file.write("\n# Optional pre-scoring Gaussian blur on input map(s) (sigma in Angstrom; 0 disables)\n")
        file.write(f'preGaussianBlur = "{args.preGaussianBlur}"\n')
        file.write(f"minimum_sigma_allowed = 0.6\n")
        file.write(f'sigma_decreasing_step = "{args.sigma_decreasing_step}"\n')        
        sigma_decay_with_improvement = "True" if bool(getattr(args, "aggressive", False)) else "False"
        file.write(f'sigma_decrease_with_improvement = "{sigma_decay_with_improvement}"\n')
        file.write('sigma_decrease_without_improvement = "True"\n')
        file.write(
            "\n# maskingCrop, if true crop the file according to mask for computing locres, in the sake of speed\n"
        )
        file.write(f'maskingCrop = "{args.maskingCrop}"\n')
        file.write('\n# aggressive mode: if True, always reload targetStar from overview\n')
        file.write(
            f'aggressive_mode = "{str(bool(getattr(args, "aggressive", False))).replace("True","True").replace("False","False")}"\n'
        )


        file.write(
            "\n# Optimization: numRecs, Number of reconstructions for each selections. More reconstruction, more precise is the selection\n"
        )
        file.write(f'numRecs = "{args.numRecs}"\n')
        file.write(f'samplingDensityFactor = "{args.samplingDensityFactor}"\n')
        file.write(f'extraSamples_num = "{args.extraSamples_num}"\n')
        file.write(f'extraSamples_randomness = "{args.extraSamples_randomness}"\n')
        file.write(f'extraSamples_audacity = "{args.extraSamples_audacity}"\n')
        file.write("\n# num of mpi processes\n")
        file.write(f'mpi = "{args.mpi}"\n')
        file.write("\n# num of eulerian views\n")
        file.write(f'numViews = "{args.numViews}"\n')
        file.write("\n# working with particle subtraction:\n")
        file.write(f'subtractionFlag= "{particleSubtractionFlag}"\n')


        gpu_str = " ".join(map(str, args.gpu)) if args.gpu else ""
        file.write("\n# GPU list (enables no-external path if non-empty):\n")
        file.write(f'gpu= "{gpu_str}"\n')
        file.write("\n# No-external flag:\n")
        file.write(f'no_external= "{str(bool(args.noExternalPrograms)).lower()}"\n')
        file.write("\n# ctf mode:\n")
        file.write(f'ctf_mode= "{str(args.ctf_mode).lower()}"\n')
        # --- CryoSPARC integration (optional) ---
        cs_project = args.cs_project if getattr(args, "cs_project", None) else ""
        cs_workspace = args.cs_workspace if getattr(args, "cs_workspace", None) else ""

        # --- CryoSPARC Local NU-refinement (optional, generated into *_run.sh) ---
        run_csparc_localnu = False
        localNURefinementMask = ""
        cs_localnu_out_basename = "reference_LNU"
        if getattr(args, "multipleLocalRefineCS", None):
            cs_project, cs_workspace, localNURefinementMask, cs_localnu_out_basename = args.multipleLocalRefineCS
            run_csparc_localnu = True

        file.write("\n# CryoSPARC integration (optional):\n")
        file.write(f'cs_project = "{cs_project}"\n')
        file.write(f'cs_workspace = "{cs_workspace}"\n')

        file.write("\n# CryoSPARC Local NU-refinement (optional):\n")
        file.write(f'run_csparc_localnu = "{str(run_csparc_localnu).lower()}"\n')
        file.write(f'localNURefinementMask = "{localNURefinementMask}"\n')
        file.write(f'cs_localnu_out_basename = "{cs_localnu_out_basename}"\n')


    if args.noExternalPrograms:
        generate_run_script_no_external(settings_file_path, gpu_list=args.gpu)
    else:
        generate_run_script(settings_file_path)





#######################################################
#######################################################
## GENERATE run script for random selection section
def generate_random_run_script(fileSettings):
    with open(fileSettings, "r") as file:
        data = toml.load(fileSettings)
    run_script_cmd = "#!/bin/bash\n\n"
    run_script_cmd += """
ensure_directory() {
    local dir="$1"
    [ ! -d "$dir" ] && mkdir -p "$dir"
}

process_random() {
    local name=$1
    local workingDir=$2
    local particles=$3
    local listSplits=$4
    local maskedCropOption=$5

    local tag=${name}
    ensure_directory ${workingDir}
    ensure_directory "${workingDir}/${tag}"
    janas_session_manager produce_reconstructions_script --mode random --i ${particles} --mask ${targetMask}  --outDir ${workingDir}/${tag} --manualParticleSubsets $listSplits --scriptName script_reconstructions.sh ${maskedCropOption} --resultFilename randomRanked_locres_values.csv
    ./${workingDir}/${tag}/script_reconstructions.sh
"""

    run_script_cmd += "}\n\n\n"
    run_script_cmd += 'targetStar="' + data.get("particles", "None") + '"\n'
    run_script_cmd += 'targetMask="' + data.get("mask", "None") + '"\n'
    run_script_cmd += 'particles="' + data.get("particles", "None") + '"\n'
    run_script_cmd += 'listSplits="' + data.get("splits", "") + '"\n'
    run_script_cmd += (
        'workingDir="'
        + os.path.join(
            data.get("dir", "./"), data.get("session_name", "randomSelection")
        )
        + '"\n'
    )
    run_script_cmd += 'masked_crop_otpion=" --masked_crop "\n'
    run_script_cmd += 'process_random "$ite" "$workingDir" "$particles" "$listSplits" "$masked_crop_otpion" \n'
    run_script_cmd += "./${workingDir}/script_reconstructions.sh"
    run_script_cmd += "\n\n\n"

    runScriptName = os.path.join(
        data.get("dir", "./"),
        data.get("session_name"),
        str(data.get("session_name")) + "_run.sh",
    )
    with open(runScriptName, "w") as f:
        f.write(run_script_cmd)
    try:
        os.chmod(runScriptName, 0o755)
    except PermissionError:
        pass


#######################################################
## CREATE random selection section
#######################################################
janas_random_selection_session = command.add_parser(
    "random_selection_session",
    description="random_selection_session",
    help="Create a new session with random selection",
)
janas_random_selection_session.add_argument(
    "--name",
    required=False,
    default="random_selection",
    type=str,
    help="name of the random_selection section, it creates a dir with that name, in not given it will assign the default name 'random_selection'",
)
janas_random_selection_session.add_argument(
    "--dir",
    required=False,
    default="./",
    type=str,
    help="working directory, if not defined will go with the current directory",
)
janas_random_selection_session.add_argument(
    "--particles",
    required=True,
    type=str,
    help="filename with the particles to be created",
)
janas_random_selection_session.add_argument(
    "--mask", required=True, type=str, help="filename with the mask"
)
janas_random_selection_session.add_argument(
    "--subsets",
    required=False,
    type=int,
    default=10,
    help="number of subsets to reconstruct (default=10)",
)


def random_selection_session(args):
    # Create the directory if it doesn't exist
    if not os.path.exists(args.name):
        os.makedirs(args.name)

    # Path for the settings file
    if not os.path.isfile(args.particles):
        print("particles star file ", args.particles, "not valid")

    if not os.path.isfile(args.mask):
        print("mask file ", args.mask, "not valid")

    params = starHandler.readStar(args.particles)
    # print('size=',params.shape[0])
    # print('splits=',args.subsets)
    interval_values = np.linspace(0, params.shape[0], args.subsets)
    interval_values_integer = map(int, interval_values[1:])
    interval_values_string = ",".join(map(str, interval_values_integer))
    # print('values=',interval_values_string)
    settings_file_path = os.path.join(args.name, "session_random_settings.txt")

    with open(settings_file_path, "w") as file:
        # Write a comment and the session_name variable
        # file.write("# name of the session\n")
        # toml.dump({'session_name': args.name}, file)
        file.write("###########################\n")
        file.write("\n# name of the session\n")
        file.write(f'session_name = "{args.name}"\n')
        file.write("\n# particles stack\n")
        file.write(f'particles = "{args.particles}"\n')
        file.write("\n# directory\n")
        file.write(f'dir = "{args.dir}"\n')
        file.write("\n# Path of the mask file\n")
        file.write(f'mask= "{args.mask}"\n')
        file.write("\n# number of particles to reconstruct:\n")
        file.write(f'splits= "{interval_values_string}"\n')

    generate_random_run_script(str(settings_file_path))


#######################################################
#######################################################
## GENERATE run script for random selection section
def generate_classification_run_script(fileSettings):
    with open(fileSettings, "r") as file:
        data = toml.load(fileSettings)

    # NEW: build GPU and no-external flags for the bash script
    gpu_cli = str(data.get("gpu", "")).strip()
    if gpu_cli:
        gpu_opt = "--gpu " + gpu_cli
    else:
        gpu_opt = ""
    no_external_flag = str(data.get("no_external", "false")).strip().lower() in (
        "true",
        "1",
        "yes",
    )

    run_script_cmd = "#!/bin/bash\n\n"
    run_script_cmd += f'GPU_ARG_OPT="{gpu_opt}"\n'
    run_script_cmd += f'NO_EXTERNAL={"1" if no_external_flag else "0"}\n'
    run_script_cmd += """
ensure_directory() {
    local dir="$1"
    [ ! -d "$dir" ] && mkdir -p "$dir"
}

rec_subset() {
        fileIn=$1
        fileOut_basename=$2
        subset=$3

        if [ "$NO_EXTERNAL" -eq 1 ]; then
            # Internal reconstructor path: build both half-maps in one call
            if [ -f "${fileOut_basename}_recH1.mrc" ] && [ -f "${fileOut_basename}_recH2.mrc" ]; then
                echo "DOING NOTHING: Internal reconstructed files ${fileOut_basename}_recH1.mrc / _recH2.mrc exist"
            else
                echo "DOING Internal Reconstruction for files ${fileOut_basename}_recH1.mrc / _recH2.mrc"
                janas_reconstructor "${fileIn}" --subset 1 2 --ctf wiener --out-basename "${fileOut_basename}" ${GPU_ARG_OPT}
            fi
        else
            # Original external RELION path
            if [ -f "${fileOut_basename}_recH${subset}.mrc" ]; then
                echo "DOING NOTHING: Reconstructed file ${fileOut_basename}_recH${subset}.mrc exists"
            else
                echo "DOING Reconstruction for file ${fileOut_basename}_recH${subset}.mrc"
                relion_reconstruct --i "${fileIn}" --o "${fileOut_basename}_recH${subset}.mrc" --subset "${subset}" --ctf &
                sleep 40
            fi
        fi
}

process_classification() {
    local name=$1
    local workingDir=$2
    local scoringDir=$3
    local classesDir=$4
    local particles=$5
    local angpix=$6
    local sigma=$7
    local mask=$8
    local listEqualizedMaps=$9
    local listClassNames=${10}
    local numMPI=${11}

    ensure_directory ${workingDir}

    
    # Convert space-separated strings to arrays
    IFS=' ' read -r -a mapsArray <<< "$listEqualizedMaps"
    IFS=' ' read -r -a classNamesArray <<< "$listClassNames"

    results_classification_array=""
    for i in "${!mapsArray[@]}"; do
        local map=${mapsArray[$i]}
        local className=${classNamesArray[$i]}
        echo "Processing map: $map with class name: $className"
        janas scoreParticles --i ${particles} --mask ${mask} --map ${map} --apix ${angpix} --sigma  ${sigma} --o "${scoringDir}/${className}".star --mpi "$numMPI"
        results_classification_array="${results_classification_array} ${scoringDir}/${className}.star"
    done

    janas_utils scores_to_csv --i $results_classification_array --o ${classesDir}/${name}_classified.star --csv ${scoringDir}/${name}.csv
    
    for i in "${!mapsArray[@]}"; do
        class_number=$((i + 1))
        janas_utils extract_particles_from_label_value --i  ${classesDir}/${name}_classified.star --o ${classesDir}/class_${class_number}.star --label _rlnClassNumber --value ${class_number}

        if [ "$NO_EXTERNAL" -eq 1 ]; then
            # Internal reconstructor builds both half-maps in one call
            rec_subset  ${classesDir}/class_${class_number}.star ${classesDir}/class_${class_number} 1
        else
            # RELION path: one call per half-map
            rec_subset  ${classesDir}/class_${class_number}.star ${classesDir}/class_${class_number} 1
            rec_subset  ${classesDir}/class_${class_number}.star ${classesDir}/class_${class_number} 2
        fi
    done

    
"""
    # janas_session_manager produce_reconstructions_script --mode random --i ${particles} --mask ${targetMask}  --outDir ${workingDir}/${tag} --manualParticleSubsets $listSplits --scriptName script_reconstructions.sh ${maskedCropOption} --resultFilename randomRanked_locres_values.csv
    # ./${workingDir}/${tag}/script_reconstructions.sh

    run_script_cmd += "}\n\n\n"
    run_script_cmd += 'name="' + data.get("session_name", "classification") + '"\n'
    run_script_cmd += 'targetStar="' + data.get("particles", "None") + '"\n'
    run_script_cmd += 'targetMask="' + data.get("mask", "None") + '"\n'
    run_script_cmd += 'particles="' + data.get("particles", "None") + '"\n'
    run_script_cmd += 'listEqualizedMaps="' + data.get("equalized_maps", "") + '"\n'
    run_script_cmd += 'listClassNames="' + data.get("scored_classes", "") + '"\n'
    run_script_cmd += 'workingDir="' + data.get("dir", "./") + '"\n'
    run_script_cmd += (
        'scoring_dir="'
        + os.path.join(
            data.get("dir", "./"), str(data.get("scoring_dir", "scoring_dir"))
        )
        + '"\n'
    )
    run_script_cmd += (
        'classes_dir="'
        + os.path.join(
            data.get("dir", "./"), str(data.get("classes_dir", "classes_dir"))
        )
        + '"\n'
    )
    run_script_cmd += 'angpix="' + data.get("angpix", "1.0") + '"\n'
    run_script_cmd += 'sigma="' + data.get("sigma", "1.0") + '"\n'
    run_script_cmd += 'numMPI="' + data.get("numMPI", "8") + '"\n'
    run_script_cmd += (
        "janas_utils equalize_images  --i "
        + data.get("input_maps", "None")
        + " --o_suffix _"
        + data.get("equalized_suffix", "equalized")
        + " --dir "
        + os.path.join(data.get("dir", "./"), data.get("equalized_suffix", "equalized"))
        + "\n"
    )
    run_script_cmd += 'process_classification "$name" "$workingDir" "$scoring_dir" "$classes_dir" "$particles" "$angpix" "$sigma" "$targetMask" "$listEqualizedMaps" "$listClassNames" "$numMPI"\n'
    # run_script_cmd += './${workingDir}/script_reconstructions.sh'
    run_script_cmd += "\n\n\n"

    runScriptName = os.path.join(
        data.get("dir", "./"), str(data.get("session_name")) + "_run.sh"
    )
    with open(runScriptName, "w") as f:
        f.write(run_script_cmd)
    try:
        os.chmod(runScriptName, 0o755)
    except PermissionError:
        pass


#######################################################
## CREATE classification selection section
#######################################################
janas_classification_session = command.add_parser(
    "classification_session",
    description="classification_session",
    help="Create a new session with classification",
)
janas_classification_session.add_argument(
    "--name",
    required=False,
    default="classification",
    type=str,
    help="name of the classification_selection section, it creates a dir with that name, in not given it will assign the default name 'random_selection'",
)
janas_classification_session.add_argument(
    "--dir",
    required=False,
    default="./",
    type=str,
    help="working directory, if not defined will go with the current directory",
)
janas_classification_session.add_argument(
    "--particles",
    required=True,
    type=str,
    help="filename with the particles to be created",
)
janas_classification_session.add_argument(
    "--mask", required=True, type=str, help="filename with the mask"
)
janas_classification_session.add_argument(
    "--maps", required=True, nargs="+", type=str, help="files with the input MRC maps"
)
janas_classification_session.add_argument(
    "--scoring_dir",
    required=False,
    default="scored_classes",
    type=str,
    help="directory with the scores to be created",
)
janas_classification_session.add_argument(
    "--classes_dir",
    required=False,
    default="final_classes",
    type=str,
    help="directory with the final classes to be created",
)
janas_classification_session.add_argument(
    "--o",
    required=False,
    type=str,
    default="classification.star",
    help="filename with the updated classes",
)
janas_classification_session.add_argument(
    "--angpix", required=False, type=str, help="pixel spacing"
)
janas_classification_session.add_argument(
    "--sigma", required=False, type=float, default="1", help="sigma for the SCI SCORE"
)
janas_classification_session.add_argument(
    "--mpi", required=False, type=int, default="8", help="number of MPI values"
)

# NEW: mirror new_select_session flags
janas_classification_session.add_argument(
    "--noExternalPrograms",
    action="store_true",
    help="if set, use internal reconstructor (janas_reconstructor) instead of external packages (e.g. relion) for class reconstructions",
)
janas_classification_session.add_argument(
    "--gpu",
    nargs="*",
    type=int,
    default=[],
    help="GPU indices to use for janas_reconstructor (e.g. --gpu 0 1). Leave empty for CPU-only.",
)
janas_classification_session.add_argument(
    "--cs_project",
    required=False,
    type=str,
    help="CryoSPARC project UID to associate with this classification session (e.g. P31).",
)

janas_classification_session.add_argument(
    "--cs_workspace",
    required=False,
    type=str,
    help="CryoSPARC workspace UID to associate with this classification session (e.g. W1).",
)

def classification_session(args):

    # Path for the settings file
    if not os.path.isfile(args.particles):
        print("particles star file ", args.particles, "not valid")
        return

    if not os.path.isfile(args.mask):
        print("mask file ", args.mask, "not valid")
        return

    for map_file in args.maps:
        if not os.path.isfile(map_file):
            print("Map file", map_file, "not valid")
            return

    # Create the directory if it doesn't exist
    if not os.path.exists(args.dir):
        os.makedirs(args.dir)
    workingDir = os.path.join(args.dir, args.name)
    if not os.path.exists(workingDir):
        os.makedirs(workingDir)
    scores_dir = os.path.join(workingDir, args.scoring_dir)
    if not os.path.exists(scores_dir):
        os.makedirs(scores_dir)
    classes_dir = os.path.join(workingDir, args.classes_dir)
    if not os.path.exists(classes_dir):
        os.makedirs(classes_dir)
    if not args.angpix:
        raw_spacing = utils.get_MRC_map_pixel_spacing(args.maps[0])
        try:
            pixel_size = raw_spacing[0]
        except (TypeError, IndexError):
            pixel_size = raw_spacing
        args.angpix = f"{pixel_size:.3f}"

    params = starHandler.readStar(args.particles)
    # print('size=',params.shape[0])
    # print('splits=',args.subsets)
    settings_file_path = os.path.join(args.name, "session_classification_settings.txt")
    maps_string = " ".join(args.maps)
    equalizedDir = os.path.join(workingDir, "equalized")
    modified_maps = [
        f"{equalizedDir}/{os.path.basename(file).split('.')[0]}_{'equalized'}.mrc"
        for file in args.maps
    ]
    modified_maps_string = " ".join(modified_maps)
    scored_classes_names = [
        f"{os.path.basename(file).split('.')[0]}_{'scoredClass'}" for file in args.maps
    ]
    scored_classes_names_string = " ".join(scored_classes_names)

    with open(settings_file_path, "w") as file:
        # Write a comment and the session_name variable
        # file.write("# name of the session\n")
        # toml.dump({'session_name': args.name}, file)
        file.write("###########################\n")
        file.write("\n# name of the session\n")
        file.write(f'session_name = "{args.name}"\n')
        file.write("\n# particles stack\n")
        file.write(f'particles = "{args.particles}"\n')
        file.write("\n# working directory\n")
        file.write(f'dir = "{workingDir}"\n')
        file.write("\n# Path of the mask file\n")
        file.write(f'mask= "{args.mask}"\n')
        file.write("\n# maps to analyze:\n")
        file.write(f'input_maps= "{maps_string}"\n')
        file.write(f'equalized_maps= "{modified_maps_string}"\n')
        file.write(f'scored_classes= "{scored_classes_names_string}"\n')
        file.write(f'equalized_suffix= "equalized"\n')
        file.write("\n# scoring dir:\n")
        file.write(f'scoring_dir= "{args.scoring_dir}"\n')
        file.write(f'classes_dir= "{args.classes_dir}"\n')
        file.write("\n# output file with updated classification:\n")
        file.write(f'output_classificated_file= "{args.o}"\n')
        file.write("\n# angpix, pixel spacing:\n")
        file.write(f'angpix= "{args.angpix}"\n')
        file.write("\n# sigma for the SCI score:\n")
        file.write(f'sigma= "{args.sigma}"\n')
        file.write("\n# number of MPI processing:\n")
        file.write(f'numMPI= "{args.mpi}"\n')
        gpu_str = " ".join(map(str, args.gpu)) if args.gpu else ""
        file.write("\n# GPU list for internal reconstructor:\n")
        file.write(f'gpu= "{gpu_str}"\n')
        file.write("\n# No-external flag for using internal reconstructor:\n")
        file.write(f'no_external= "{str(bool(args.noExternalPrograms)).lower()}"\n')
        # --- CryoSPARC integration (optional) ---
        cs_project = args.cs_project if getattr(args, "cs_project", None) else ""
        cs_workspace = args.cs_workspace if getattr(args, "cs_workspace", None) else ""

        file.write("\n# CryoSPARC integration (optional):\n")
        file.write(f'cs_project = "{cs_project}"\n')
        file.write(f'cs_workspace = "{cs_workspace}"\n')

    generate_classification_run_script(str(settings_file_path))


def main(command_line=None):
    args = janas_parser.parse_args(command_line)
    if args.command == "produce_reconstructions_script":
        produce_reconstructions_script(args)
    elif args.command == "random_selection_session":
        random_selection_session(args)
    elif args.command == "classification_session":
        classification_session(args)
    elif args.command == "new_select_session":
        new_select_session(args)
    else:
        janas_parser.print_help()


if __name__ == "__main__":
    main()
