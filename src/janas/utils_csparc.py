# File: utils_csparc.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology

"""
CryoSPARC integration utilities for JANAS, via cryosparc-tools.

Features:
- open a CryoSPARC session using cryosparc.tools.CryoSPARC
- submit import / ab-initio / non-uniform / local refinement jobs
- fetch refined particle metadata and update a Relion STAR using
  utils.update_star_from_csparc

This module is optional. If cryosparc-tools is not installed or a CryoSPARC
instance is not available, callers will get a RuntimeError only when
they invoke the high-level functions here; import of the module itself
remains safe.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional
import time
import datetime
import numpy as np
import shutil

from janas import utils
# Optional import: do not explode if cryosparc-tools is missing.
try:
    from cryosparc.tools import CryoSPARC
except ImportError:
    CryoSPARC = None  # type: ignore[assignment]


def sanitise_star_optics_for_cryosparc(star_path: str) -> None:
    """
    In-place STAR sanity check for CryoSPARC import.

    Some RELION 3.1+ STAR files contain duplicated microscope parameters in both
    data_optics and data_particles (e.g. _rlnVoltage, _rlnImagePixelSize, etc.).
    CryoSPARC's STAR reader is fragile in these cases.

    Strategy:
      - Parse the data_optics and data_particles loop headers as raw text.
      - Find STAR tags that appear in both sections.
      - In data_optics, drop all duplicated tags *except*:
            _rlnOpticsGroup
            _rlnImageSize
            _rlnImageDimensionality
      - Keep all non-duplicated optics tags.
      - Rebuild the optics header and records with renumbered column indices.
      - The file is modified in place only if a change is actually needed.

    If the file has no data_optics / data_particles sections, or is already
    compatible, the function is a no-op.
    """
    star_file = Path(star_path)
    if not star_file.is_file():
        # Let the caller fail later in the usual way.
        return

    text = star_file.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines(keepends=True)
    n = len(lines)

    def find_section(section_name: str) -> int:
        target = f"data_{section_name}"
        for i, l in enumerate(lines):
            if l.strip() == target:
                return i
        return -1

    optics_idx = find_section("optics")
    particles_idx = find_section("particles")
    if optics_idx < 0 or particles_idx < 0:
        # Not a RELION 3.1+ particles STAR or no optics section: nothing to do.
        return

    def parse_loop_header(start_idx: int):
        """
        From a 'data_xxx' line index, locate the following:

          data_xxx
          ...
          loop_
          _rlnSomething #1
          _rlnSomethingElse #2
          ...

        Returns (loop_idx, header_start_idx, header_end_idx, tag_list)
        or None if no loop_ header is found.
        """
        loop_idx = -1
        for i in range(start_idx + 1, n):
            stripped = lines[i].strip()
            if stripped == "loop_":
                loop_idx = i
                break
            if stripped.startswith("data_") and i != start_idx:
                break

        if loop_idx < 0:
            return None

        header_start = loop_idx + 1
        header_end = header_start
        tags: list[str] = []
        while header_end < n:
            stripped = lines[header_end].lstrip()
            if stripped.startswith("_"):
                tag = stripped.split()[0]
                tags.append(tag)
                header_end += 1
            else:
                break

        if not tags:
            return None

        return loop_idx, header_start, header_end, tags

    optics_parsed = parse_loop_header(optics_idx)
    particles_parsed = parse_loop_header(particles_idx)
    if optics_parsed is None or particles_parsed is None:
        # No loop_ in optics or particles: keep file untouched.
        return

    loop_idx_o, header_start_o, header_end_o, optics_tags = optics_parsed
    _, _, _, particles_tags = particles_parsed

    ALWAYS_KEEP = {"_rlnOpticsGroup", "_rlnImageSize", "_rlnImageDimensionality"}
    duplicates = set(optics_tags) & set(particles_tags)

    # Decide which optics columns to keep
    keep_indices: list[int] = []
    for idx, tag in enumerate(optics_tags):
        if tag in ALWAYS_KEEP:
            keep_indices.append(idx)
        elif tag in duplicates:
            # Drop duplicated microscope parameters that also appear per-particle
            continue
        else:
            # Non-duplicated optics-only metadata: keep it
            keep_indices.append(idx)

    # Nothing to change
    if len(keep_indices) == len(optics_tags):
        return

    # Rebuild optics header
    new_header_lines: list[str] = []
    for new_col, old_idx in enumerate(keep_indices, start=1):
        old_line = lines[header_start_o + old_idx]
        stripped = old_line.lstrip()
        tag = stripped.split()[0]
        leading = old_line[: old_line.find(tag)]
        new_header_lines.append(f"{leading}{tag} #{new_col}\n")

    # Data block for optics: from header_end_o up to (but not including) next data_*
    data_start = header_end_o
    data_end = data_start
    while data_end < n and not lines[data_end].strip().startswith("data_"):
        data_end += 1

    new_data_lines: list[str] = []
    for i in range(data_start, data_end):
        line = lines[i]
        stripped = line.strip()
        if not stripped or line.lstrip().startswith("#"):
            # Preserve blank/comment lines as they are
            new_data_lines.append(line)
            continue

        tokens = line.split()
        if len(tokens) < len(optics_tags):
            # Unexpected row layout, do not touch this line
            new_data_lines.append(line)
            continue

        new_tokens = [tokens[j] for j in keep_indices if j < len(tokens)]
        leading = line[: len(line) - len(line.lstrip())]
        new_data_lines.append(leading + " ".join(new_tokens) + "\n")

    # Stitch new STAR
    new_lines: list[str] = []
    new_lines.extend(lines[:header_start_o])
    new_lines.extend(new_header_lines)
    new_lines.extend(new_data_lines)
    new_lines.extend(lines[data_end:])

    new_text = "".join(new_lines)
    if new_text != text:
        star_file.write_text(new_text, encoding="utf-8")

class CryoSPARCToolsSession:
    """
    Thin wrapper around cryosparc.tools.CryoSPARC for a given
    (project_uid, workspace_uid, lane, user_email) context.

    Only relies on the HTTP API and environment variables like
    CRYOSPARC_LICENSE_ID, CRYOSPARC_MASTER_HOSTNAME, CRYOSPARC_BASE_PORT.
    """

    def __init__(
        self,
        user_email: str,
        project_uid: str,
        workspace_uid: str,
        lane: str,
        license_id: Optional[str] = None,
        host: Optional[str] = None,
        base_port: Optional[int] = None,
        password: Optional[str] = None,
    ) -> None:
        if CryoSPARC is None:
            raise RuntimeError(
                "CryoSPARC integration is not available: the 'cryosparc-tools' "
                "Python package is not installed. Install it with "
                "'pip install cryosparc-tools' or avoid using the "
                "'csparc_nurefinement' command."
            )

        # Defaults: let CryoSPARC pick up values from environment if not given.
        license_id = license_id or os.getenv("CRYOSPARC_LICENSE_ID", "")
        host = host or os.getenv("CRYOSPARC_MASTER_HOSTNAME", "localhost")
        base_port = base_port or int(os.getenv("CRYOSPARC_BASE_PORT", "39000"))
        # For password: either pass explicit, or rely on CRYOSPARC_PASSWORD
        password = password or os.getenv("CRYOSPARC_PASSWORD", "test")

        # Construct session. If authentication or connection fails, this will raise.
        self._cs = CryoSPARC(
            license=license_id,
            host=host,
            base_port=base_port,
            email=user_email,
            password=password,
        )

        # Check that the requested lane exists, to fail with a clear message.
        lanes_info = self._cs.cli.get_scheduler_lanes()

        # cryosparc-tools può restituire:
        #  - una lista di dict, ciascuno con 'name'
        #  - oppure un dict con chiave 'lanes'
        if isinstance(lanes_info, dict):
            lanes = lanes_info.get("lanes", [])
        elif isinstance(lanes_info, list):
            lanes = lanes_info
        else:
            lanes = []

        lane_names = [
            ln.get("name")
            for ln in lanes
            if isinstance(ln, dict) and "name" in ln
        ]

        if lane not in lane_names:
            raise RuntimeError(
                f"Requested lane '{lane}' does not exist in CryoSPARC. "
                f"Available lanes: {', '.join(lane_names) if lane_names else '(none)'}"
            )

        self.project_uid = project_uid
        self.workspace_uid = workspace_uid
        self.lane = lane

    # ----- job helpers -----
    def get_project_dir(self) -> Path:
        """
        Return the absolute CryoSPARC project directory for this session.
        """
        proj_dir = self._cs.cli.get_project_dir_abs(self.project_uid)
        return Path(proj_dir)

    def _queue_and_wait(self, job, label: Optional[str] = None):
        """
        Queue a job to the configured lane and wait until completion.
        Prints simple status updates with elapsed time.
        """
        name = label or job.type
        print(f"  [csparc] Queuing job {job.uid} ({name}) on lane '{self.lane}'", flush=True)

        job.queue(self.lane)

        start = time.time()
        last_status = None
        while True:
            job.refresh()  # ricarica status dal server
            status = job.status
            elapsed = time.time() - start
            if status != last_status:
                ts = datetime.datetime.now().strftime("%H:%M:%S")
                print(
                    f"[csparc] {ts} | {job.uid} ({name}) status: {status} "
                    f"(elapsed {elapsed:5.1f}s)"
                )
                last_status = status

            if status in ("failed", "killed"):
                try:
                    streamlog = self._cs.cli.get_job_streamlog(self.project_uid, job.uid)
                    print(f"[csparc] --- last log lines for {job.uid} ---")
                    last = [e.get("text", "").rstrip() for e in streamlog if "text" in e]
                    for line in last[-20:]:
                        if line:
                            print(f"[csparc] {line}")
                except Exception as e:
                    print(f"[csparc] Could not fetch streamlog for {job.uid}: {e}")
                raise RuntimeError(
                    f"CryoSPARC job {job.uid} ended with status={status}"
                )
            if status == "completed":
                break

            time.sleep(15)

        return job

    def import_particles_from_star(self, star_path: str, particle_dir: str):
        """
        Import a Relion STAR into CryoSPARC as an Import Particles job.
        Returns the Job instance.
        """
        star_path = str(Path(star_path).absolute())
        particle_dir = str(Path(particle_dir).absolute())

        job = self._cs.create_job(
            self.project_uid,
            self.workspace_uid,
            "import_particles",
            params={
                "particle_meta_path": star_path,
                "particle_blob_path": particle_dir,
            },
        )
        return self._queue_and_wait(job, label="import_particles")

    def abinit(self, particles_job):
        """
        Ab-initio reconstruction from imported particles.
        Returns (Job, output_group_name).
        """
        job = self._cs.create_job(
            self.project_uid,
            self.workspace_uid,
            "homo_abinit",
            connections={"particles": (particles_job.uid, "imported_particles")},
        )
        job = self._queue_and_wait(job, label="homo_abinit")
        return job, "volume_class_0"

    def import_reference_volume(self, ref_map: str):
        """
        Import a reference volume.
        Returns (Job, output_group_name).
        """
        ref_map = str(Path(ref_map).absolute())
        job = self._cs.create_job(
            self.project_uid,
            self.workspace_uid,
            "import_volumes",
            params={"volume_blob_path": ref_map},
        )
        job = self._queue_and_wait(job, label="import_volumes(ref)")
        return job, "imported_volume_1"

    def import_mask(self, mask_map: str):
        mask_map = str(Path(mask_map).absolute())

        # TODO: sostituire 'PARAM_KEY' e VALUE con quelli che leggi da doc["params"]
        params = {
            "volume_blob_path": mask_map,
            "volume_out_name": "mask",
        }

        job = self._cs.create_job(
            self.project_uid,
            self.workspace_uid,
            "import_volumes",
            params=params,
        )
        job = self._queue_and_wait(job, label="import_volumes(mask)")

        return job, "imported_mask_1"


    def local_refine_from_imports(
        self,
        particles_job,
        volume_job,
        volume_output: str,
        mask_job=None,
        mask_output: Optional[str] = None,
        symmetry: str = "C1",
        min_angular_step: float = 0.2,
        resplit: bool = False,
    ):
        """
        Run a 'new_local_refine' job directly from:

        - particles_job: import_particles job (with alignments3D)
        - volume_job / volume_output: imported reference volume
        - mask_job / mask_output: imported static mask (optional)

        This corresponds to the GUI 'Local refinement' job where the inputs
        are particles + volume + mask.
        """

        # Build CryoSPARC-style connection dict
        connections = {
            "particles": (particles_job.uid, "imported_particles"),
            "volume": (volume_job.uid, volume_output),
        }
        if mask_job is not None and mask_output is not None:
            connections["mask"] = (mask_job.uid, mask_output)

        # Parameters corresponding to the job.json you showed
        params = {
            "refine_symmetry": symmetry,
            # This is implicit for this job type, but we can be explicit:
            "do_local_refine": True,
            "min_angular_step": float(min_angular_step),
            # Useful in practice – same meaning as in the GUI
            "use_alignment_prior": False,
            "refine_gs_resplit": bool(resplit),
        }

        job = self._cs.create_job(
            self.project_uid,
            self.workspace_uid,
            "new_local_refine",
            params=params,
            connections=connections,
        )
        return self._queue_and_wait(job, label="new_local_refine")


    def nonuniform_refine(
        self,
        particles_job,
        volume_job,
        volume_output: str,
        mask_job=None,
        mask_output: Optional[str] = None,
        symmetry: str = "C1",
        ini_high: Optional[float] = None,
        resplit: bool = False,
    ):
        """
        Run non-uniform refinement connecting a particles job and a volume job.
        Se mask_job/mask_output sono forniti e l'output del job di mask è di tipo
        'mask', la mask viene collegata allo slot 'mask' del NU.
        """
        params = {
            "refine_symmetry": symmetry,
            "refine_gs_resplit": bool(resplit),
        }
        if ini_high is not None:
            params["refine_res_init"] = float(ini_high)

        connections = {
            "particles": (particles_job.uid, "imported_particles"),
            "volume": (volume_job.uid, volume_output),
        }

        # Qui usiamo ESATTAMENTE mask_job e mask_output se ci sono
        if mask_job is not None and mask_output is not None:
            connections["mask"] = (mask_job.uid, mask_output)

        job = self._cs.create_job(
            self.project_uid,
            self.workspace_uid,
            "nonuniform_refine_new",
            connections=connections,
            params=params,
        )
        return self._queue_and_wait(job, label="nonuniform_refine_new")

    def local_refine(
        self,
        refine_job,
        symmetry: str = "C1",
        min_angular_step: float = 0.01,
    ):
        """
        Run Local Refinement on top of an existing refinement job.
        Returns the Job instance.
        """
        job = self._cs.create_job(
            self.project_uid,
            self.workspace_uid,
            "new_local_refine",
            connections={
                "particles": (refine_job.uid, "particles"),
                "volume": (refine_job.uid, "volume"),
                "mask": (refine_job.uid, "mask"),
            },
            params={
                "refine_symmetry": symmetry,
                "use_alignment_prior": True,
                "min_angular_step": float(min_angular_step),
            },
        )
        return self._queue_and_wait(job, label="new_local_refine")


import shutil

def export_refine_maps_from_jobdir(
    job_dir: Path,
    refine_job_uid: str,
    output_basename: str,
):
    """
    Dato job_dir (es. .../P31/J17) ed un basename (es. 'class_1_updatedNU'),
    copia le mappe finali di nonuniform_refine_new in:

        <basename>_rec.mrc
        <basename>_recH1.mrc
        <basename>_recH2.mrc
    """
    base = Path(output_basename)
    out_dir = base.parent if base.parent != Path("") else Path(".")
    stem = base.name

    # full map
    full_candidates = sorted(job_dir.glob(f"{refine_job_uid}_*_volume_map.mrc"))
    # half-maps
    halfA_candidates = sorted(job_dir.glob(f"{refine_job_uid}_*_volume_map_half_A.mrc"))
    halfB_candidates = sorted(job_dir.glob(f"{refine_job_uid}_*_volume_map_half_B.mrc"))

    def last_or_none(cands):
        return cands[-1] if cands else None

    full_map = last_or_none(full_candidates)
    halfA = last_or_none(halfA_candidates)
    halfB = last_or_none(halfB_candidates)

    def safe_copy(src: Optional[Path], dest: Path):
        if src is None:
            return
        try:
            shutil.copy2(src, dest)
        except Exception as e:
            print(f"[csparc_nurefinement] WARNING: failed to copy {src} -> {dest}: {e}")

    safe_copy(full_map, out_dir / f"{stem}_rec.mrc")
    safe_copy(halfA,   out_dir / f"{stem}_recH1.mrc")
    safe_copy(halfB,   out_dir / f"{stem}_recH2.mrc")



def _finalise_refinement_job(
    session: CryoSPARCToolsSession,
    refine_job_uid: str,
    star_in: str,
    star_out: str,
    loglevel: str,
    output_basename: Optional[str] = None,
) -> str:
    """
    Common post-processing for a CryoSPARC refinement job:
    - locate job dir
    - find refined particles .cs
    - update Relion STAR
    - export maps (if requested)
    Returns the path to the .cs file.
    """
    proj_dir = session.get_project_dir()
    job_dir = proj_dir / refine_job_uid

    if not job_dir.is_dir():
        raise FileNotFoundError(
            f"CryoSPARC job directory does not exist: {job_dir}"
        )

    # Find refined particle .cs (e.g. J89_002_particles.cs)
    particle_cs_candidates = sorted(
        p for p in job_dir.glob(f"{refine_job_uid}_*_particles.cs")
        if "passthrough" not in p.name
    )
    if not particle_cs_candidates:
        particle_cs_candidates = sorted(
            p for p in job_dir.glob("*_particles.cs")
            if "passthrough" not in p.name
        )

    if not particle_cs_candidates:
        raise FileNotFoundError(
            f"No refined particle .cs files found in {job_dir} for job {refine_job_uid}"
        )

    csfile = particle_cs_candidates[-1]

    # Update STAR from .cs
    utils.update_star_from_csparc(
        csfile=str(csfile),
        starfile_in=star_in,
        starfile_out=star_out,
        loglevel=loglevel,
    )

    # Export maps if requested
    if output_basename is not None:
        export_refine_maps_from_jobdir(
            job_dir=job_dir,
            refine_job_uid=refine_job_uid,
            output_basename=output_basename,
        )

    return str(csfile)


def run_csparc_nu_refinement_and_update_star(
    star_in: str,
    star_out: str,
    particle_dir: str,
    user_email: str,
    project_uid: str,
    workspace_uid: str,
    lane: str,
    symmetry: str = "C1",
    ref_map: Optional[str] = None,
    ini_high: Optional[float] = None,
    resplit: bool = False,
    do_local: bool = False,
    min_angular_step: float = 0.01,
    loglevel: str = "WARNING",
    license_id: Optional[str] = None,
    host: Optional[str] = None,
    base_port: Optional[int] = None,
    password: Optional[str] = None,
    mask_map: Optional[str] = None,
    output_basename: Optional[str] = None,
    precomputed_job_uid: Optional[str] = None,
) -> str:
    """
    High-level pipeline:

      Relion STAR  -> Import Particles
                    -> [ Ab-initio | Import Volumes ]
                    -> Non-uniform Refinement
                    -> [ Local Refinement ]
                    -> extract refined particle metadata
                    -> update Relion STAR in-place using update_star_from_csparc.

    Returns the path of the .cs file that was used to update the STAR.
    """
    sanitise_star_optics_for_cryosparc(star_in)

    # 1) open CryoSPARC session
    session = CryoSPARCToolsSession(
        user_email=user_email,
        project_uid=project_uid,
        workspace_uid=workspace_uid,
        lane=lane,
        license_id=license_id,
        host=host,
        base_port=base_port,
        password=password,
    )

    # ------------------------------------------------------------------
    # FAST PATH: reuse an existing NU/local refinement job
    # ------------------------------------------------------------------
    if precomputed_job_uid is not None:
        refine_job_uid = precomputed_job_uid.strip()
        proj_dir = session.get_project_dir()
        job_dir = proj_dir / refine_job_uid

        print(
            f"[csparc_nurefinement] Using precomputed job {refine_job_uid} "
            f"in project dir {proj_dir}, job dir {job_dir}"
        )

        return _finalise_refinement_job(
            session=session,
            refine_job_uid=refine_job_uid,
            star_in=star_in,
            star_out=star_out,
            loglevel=loglevel,
            output_basename=output_basename,
        )

    # 2) Import particles from STAR
    import_job = session.import_particles_from_star(
        star_path=star_in,
        particle_dir=particle_dir,
    )

    # 3) expose CryoSPARC project dir + symlink
    proj_dir = session.get_project_dir()
    print(f"[csparc_nurefinement] CryoSPARC project directory: {proj_dir}")
    link = Path.cwd() / f"{project_uid}_cryosparc"
    if not link.exists():
        try:
            os.symlink(proj_dir, link)
            print(f"[csparc_nurefinement] Created symlink: {link} -> {proj_dir}")
        except OSError as e:
            print(f"[csparc_nurefinement] Could not create symlink {link}: {e}")

    # 4) reference volume: either import or ab-initio
    if ref_map is not None:
        volume_job, volume_output = session.import_reference_volume(ref_map)
    else:
        volume_job, volume_output = session.abinit(import_job)

    # 4b) mask, if provided
    mask_job = None
    mask_output = None
    if mask_map is not None:
        mask_job, mask_output = session.import_mask(mask_map)

    # 5) non-uniform refinement
    refine_job = session.nonuniform_refine(
        particles_job=import_job,
        volume_job=volume_job,
        volume_output=volume_output,
        mask_job=mask_job,
        mask_output=mask_output,
        symmetry=symmetry,
        ini_high=ini_high,
        resplit=resplit,
    )

    # 5b) optional local refinement *on top* of NU
    if do_local:
        refine_job = session.local_refine(
            refine_job=refine_job,
            symmetry=symmetry,
            min_angular_step=min_angular_step,
        )

    # 6–8) common post-processing
    return _finalise_refinement_job(
        session=session,
        refine_job_uid=refine_job.uid,
        star_in=star_in,
        star_out=star_out,
        loglevel=loglevel,
        output_basename=output_basename,
    )




def run_csparc_local_refinement_and_update_star(
    star_in: str,
    star_out: str,
    particle_dir: str,
    user_email: str,
    project_uid: str,
    workspace_uid: str,
    lane: str,
    symmetry: str = "C1",
    ref_map: Optional[str] = None,
    resplit: bool = False,
    min_angular_step: float = 0.2,
    loglevel: str = "WARNING",
    license_id: Optional[str] = None,
    host: Optional[str] = None,
    base_port: Optional[int] = None,
    password: Optional[str] = None,
    mask_map: Optional[str] = None,
    output_basename: Optional[str] = None,
    precomputed_job_uid: Optional[str] = None,
) -> str:
    """
    Local refinement pipeline:

      Relion STAR  -> Import Particles
                    -> Import Reference Volume
                    -> Import Mask (static)
                    -> new_local_refine
                    -> extract refined particle metadata
                    -> update Relion STAR via update_star_from_csparc

    Returns the path of the .cs file used to update the STAR.
    """

    if ref_map is None and precomputed_job_uid is None:
        # Only require ref_map when we actually need to run a new refinement
        raise ValueError(
            "run_csparc_local_refinement_and_update_star requires a reference map "
            "(ref_map must not be None) for local refinement, unless a precomputed "
            "job UID is provided."
        )

    sanitise_star_optics_for_cryosparc(star_in)

    # 1) open CryoSPARC session
    session = CryoSPARCToolsSession(
        user_email=user_email,
        project_uid=project_uid,
        workspace_uid=workspace_uid,
        lane=lane,
        license_id=license_id,
        host=host,
        base_port=base_port,
        password=password,
    )

    # ------------------------------------------------------------------
    # FAST PATH: reuse an existing refinement job (no new computation)
    # ------------------------------------------------------------------
    if precomputed_job_uid is not None:
        refine_job_uid = precomputed_job_uid.strip()
        proj_dir = session.get_project_dir()
        job_dir = proj_dir / refine_job_uid

        print(
            f"[csparc_localnurefinement] Using precomputed job {refine_job_uid} "
            f"in project dir {proj_dir}, job dir {job_dir}"
        )

        return _finalise_refinement_job(
            session=session,
            refine_job_uid=refine_job_uid,
            star_in=star_in,
            star_out=star_out,
            loglevel=loglevel,
            output_basename=output_basename,
        )

    # 2) import particles from STAR
    import_job = session.import_particles_from_star(
        star_path=star_in,
        particle_dir=particle_dir,
    )

    # 3) expose CryoSPARC project dir + symlink
    proj_dir = session.get_project_dir()
    print(f"[csparc_localnurefinement] CryoSPARC project directory: {proj_dir}")
    link = Path.cwd() / f"{project_uid}_cryosparc"
    if not link.exists():
        try:
            os.symlink(proj_dir, link)
            print(f"[csparc_localnurefinement] Created symlink: {link} -> {proj_dir}")
        except OSError as e:
            print(f"[csparc_localnurefinement] Could not create symlink {link}: {e}")

    # 4) import reference volume
    volume_job, volume_output = session.import_reference_volume(ref_map)

    # 5) import mask if provided
    mask_job = None
    mask_output = None
    if mask_map is not None:
        mask_job, mask_output = session.import_mask(mask_map)

    # 6) run local refinement directly from imported particles + volume + mask
    refine_job = session.local_refine_from_imports(
        particles_job=import_job,
        volume_job=volume_job,
        volume_output=volume_output,
        mask_job=mask_job,
        mask_output=mask_output,
        symmetry=symmetry,
        min_angular_step=min_angular_step,
        resplit=resplit,
    )

    # 7–9) common post-processing
    return _finalise_refinement_job(
        session=session,
        refine_job_uid=refine_job.uid,
        star_in=star_in,
        star_out=star_out,
        loglevel=loglevel,
        output_basename=output_basename,
    )


