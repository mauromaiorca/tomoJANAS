#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/io/logs.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
tomoJANAS logging and validation-log writers.

* Every command appends a structured entry to ``logs/tomojanas_import.log``.
* The validate command (and any command run with ``--validate``) writes a
  machine + human readable validation report to ``logs/validation_log.{star,json,md}``.

Timestamps use the local wall clock (``datetime.now``); this is ordinary
runtime code (not a workflow script), so that is fine.
"""

from __future__ import annotations

import json
import os
import shlex
from datetime import datetime
from typing import Any, Dict, List, Optional, Sequence

from . import project_writer
from .star_writer import LoopBlock, PairBlock, write_star

__all__ = [
    "now_iso",
    "ImportLogger",
    "write_validation_logs",
    "record_command",
]


def record_command(
    project_root: str,
    argv: Sequence[str],
    exit_status: int,
    *,
    prog: str = "tomojanas-import",
    version: str = "",
    cwd: Optional[str] = None,
) -> None:
    """Append an executed command to the project command log.

    Writes two files under ``<root>/logs``:
      * ``commands.jsonl`` — one JSON object per line (machine-readable, with
        timestamp, full argv, exit status, cwd, version);
      * ``commands.sh`` — a replayable shell script (one command per entry,
        prefixed by a comment with timestamp + exit status). Replay in a new
        project by editing the ``--project`` value and running ``bash``.

    Failures here never raise (logging must not break the command).
    """
    try:
        ld = project_writer.logs_dir(project_root)
        os.makedirs(ld, exist_ok=True)
        ts = now_iso()
        argv = list(argv)
        cmd_str = prog + " " + " ".join(shlex.quote(a) for a in argv)

        jsonl = os.path.join(ld, "commands.jsonl")
        with open(jsonl, "a", encoding="utf-8") as f:
            f.write(json.dumps({
                "timestamp": ts,
                "tomojanas_version": version,
                "prog": prog,
                "argv": argv,
                "command": cmd_str,
                "exit_status": int(exit_status),
                "cwd": cwd or os.getcwd(),
            }) + "\n")

        sh = os.path.join(ld, "commands.sh")
        new_file = not os.path.isfile(sh)
        with open(sh, "a", encoding="utf-8") as f:
            if new_file:
                f.write("#!/usr/bin/env bash\n")
                f.write("# tomoJANAS command replay log.\n")
                f.write("# Re-run these commands (edit --project for a new project).\n")
                f.write("# NOTE: only the command lines are replayed; comments record exit status.\n\n")
            f.write(f"# [{ts}] exit={int(exit_status)}\n")
            f.write(cmd_str + "\n\n")
    except Exception:
        pass

_LOG_NAME = "tomojanas_import.log"


def now_iso() -> str:
    return datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


class ImportLogger:
    """Append-only logger writing to ``<root>/logs/tomojanas_import.log``.

    Collects warnings/errors in memory too, so a command can summarise them
    at the end and feed them into the manifest / validation report.
    """

    def __init__(self, root: str, version: str = ""):
        self.root = root
        self.version = version
        self.warnings: List[str] = []
        self.errors: List[str] = []
        os.makedirs(project_writer.logs_dir(root), exist_ok=True)
        self.path = os.path.join(project_writer.logs_dir(root), _LOG_NAME)

    # -- low level ---------------------------------------------------------- #
    def _append(self, text: str) -> None:
        with open(self.path, "a", encoding="utf-8") as f:
            f.write(text)
            if not text.endswith("\n"):
                f.write("\n")

    def line(self, level: str, message: str) -> None:
        self._append(f"[{now_iso()}] {level.upper():7s} {message}")

    def info(self, message: str) -> None:
        self.line("INFO", message)

    def warning(self, message: str) -> None:
        self.warnings.append(message)
        self.line("WARNING", message)

    def error(self, message: str) -> None:
        self.errors.append(message)
        self.line("ERROR", message)

    # -- structured invocation header -------------------------------------- #
    def invocation(
        self,
        command: str,
        argv: Optional[List[str]] = None,
        resolved_inputs: Optional[Dict[str, Any]] = None,
        generated_files: Optional[List[str]] = None,
    ) -> None:
        lines = [
            "=" * 72,
            f"[{now_iso()}] tomoJANAS {self.version} :: command='{command}'",
        ]
        if argv is not None:
            lines.append("  argv: " + " ".join(str(a) for a in argv))
        if resolved_inputs:
            lines.append("  resolved inputs:")
            for k, v in resolved_inputs.items():
                lines.append(f"    {k} = {v}")
        if generated_files:
            lines.append("  generated files:")
            for g in generated_files:
                lines.append(f"    {g}")
        self._append("\n".join(lines))

    def summary(self, validation_summary: Optional[str] = None) -> None:
        self._append(
            f"  summary: {len(self.warnings)} warning(s), "
            f"{len(self.errors)} error(s)"
        )
        if validation_summary:
            self._append(f"  validation: {validation_summary}")


# --------------------------------------------------------------------------- #
# validation report writers
# --------------------------------------------------------------------------- #
def _md_table(items: List[Dict[str, Any]]) -> str:
    if not items:
        return "_No checks recorded._\n"
    cols = ["check", "status", "severity", "detail"]
    out = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for it in items:
        out.append(
            "| "
            + " | ".join(
                str(it.get(c, "")).replace("|", "\\|").replace("\n", " ")
                for c in cols
            )
            + " |"
        )
    return "\n".join(out) + "\n"


def write_validation_logs(root: str, report: Dict[str, Any]) -> Dict[str, str]:
    """Write ``validation_log.{star,json,md}`` under ``<root>/logs``.

    ``report`` shape::

        {
          "generated": iso8601,
          "tomojanas_version": "...",
          "scope": "project|tomogram|particle",
          "ok": bool,
          "summary": "...",
          "items": [ {check, status(pass|fail|warn|skip), severity, detail}, ... ],
        }
    """
    ld = project_writer.logs_dir(root)
    os.makedirs(ld, exist_ok=True)
    json_path = os.path.join(ld, "validation_log.json")
    star_path = os.path.join(ld, "validation_log.star")
    md_path = os.path.join(ld, "validation_log.md")

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, default=str)

    items = report.get("items", []) or []
    scalars = {
        "_tomoJANASValidationGenerated": report.get("generated", now_iso()),
        "_tomoJANASValidationScope": report.get("scope", "project"),
        "_tomoJANASValidationOk": bool(report.get("ok", False)),
        "_tomoJANASValidationSummary": report.get("summary", ""),
    }
    blocks = [PairBlock(name="tomoJANAS_validation", pairs=scalars)]
    if items:
        cols = [
            "_tomoJANASCheckName",
            "_tomoJANASCheckStatus",
            "_tomoJANASCheckSeverity",
            "_tomoJANASCheckDetail",
        ]
        rows = [
            [it.get("check"), it.get("status"), it.get("severity", "info"),
             it.get("detail", "")]
            for it in items
        ]
        blocks.append(LoopBlock(name="tomoJANAS_validation_items", columns=cols, rows=rows))
    write_star(star_path, blocks)

    md = [
        f"# tomoJANAS validation report",
        "",
        f"- generated: {report.get('generated', now_iso())}",
        f"- scope: {report.get('scope', 'project')}",
        f"- result: {'OK' if report.get('ok') else 'FAILED'}",
        f"- summary: {report.get('summary', '')}",
        "",
        _md_table(items),
    ]
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("\n".join(md))

    return {"json": json_path, "star": star_path, "md": md_path}
