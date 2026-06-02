#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/models/project.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
Lightweight project handle that wraps a project root and exposes the
canonical tomoJANAS paths plus existence/create policy helpers.
"""

from __future__ import annotations

import os
from dataclasses import dataclass

from ..io import project_writer


class ProjectExistsError(RuntimeError):
    pass


class ProjectMissingError(RuntimeError):
    pass


@dataclass
class Project:
    root: str

    # -- paths -------------------------------------------------------------- #
    @property
    def tomograms_star(self) -> str:
        return project_writer.tomograms_star_path(self.root)

    @property
    def particles_all_star(self) -> str:
        return project_writer.particles_all_star_path(self.root)

    @property
    def optimisation_set_star(self) -> str:
        return project_writer.optimisation_set_star_path(self.root)

    @property
    def logs_dir(self) -> str:
        return project_writer.logs_dir(self.root)

    @property
    def tilt_series_dir(self) -> str:
        return project_writer.tilt_series_dir(self.root)

    def tilt_series_star(self, tomo_name: str) -> str:
        return project_writer.tilt_series_star_path(self.root, tomo_name)

    def tomogram_dir(self, tomo_name: str) -> str:
        return project_writer.tomogram_dir(self.root, tomo_name)

    def imod_settings_dir(self, tomo_name: str) -> str:
        return project_writer.imod_settings_dir(self.root, tomo_name)

    def individual_particles_dir(self, tomo_name: str) -> str:
        return project_writer.individual_particles_dir(self.root, tomo_name)

    # -- state -------------------------------------------------------------- #
    @property
    def exists(self) -> bool:
        return os.path.isdir(self.root) and os.path.isfile(self.tomograms_star)

    def store(self, path: str, relative: bool = True) -> str:
        return project_writer.store_path(path, self.root, relative=relative)

    # -- lifecycle ---------------------------------------------------------- #
    @classmethod
    def open_or_create(
        cls,
        root: str,
        *,
        create_if_missing: bool = False,
        overwrite: bool = False,
        fail_if_existing: bool = True,
        append: bool = False,
        update_existing: bool = False,
    ) -> "Project":
        """Resolve project create/append/overwrite policy and ensure the
        skeleton directories exist.

        Default policy (matching the CLI): ``fail_if_existing`` for a project
        that already has tomograms.star, unless ``append``/``update_existing``/
        ``overwrite``/``create_if_missing`` say otherwise.
        """
        proj = cls(root=os.path.abspath(root))
        already = proj.exists
        root_dir_present = os.path.isdir(proj.root)

        if already:
            if overwrite:
                # Caller is responsible for clearing stale files; we just allow it.
                pass
            elif append or update_existing:
                pass
            elif fail_if_existing:
                raise ProjectExistsError(
                    f"project already exists at '{root}'. Use --append, "
                    f"--update-existing or --overwrite (or --fail-if-existing to abort)."
                )
        else:
            if not root_dir_present and not create_if_missing:
                raise ProjectMissingError(
                    f"project '{root}' does not exist. Pass --create-if-missing "
                    f"to create a new project."
                )

        project_writer.ensure_project(proj.root)
        return proj
