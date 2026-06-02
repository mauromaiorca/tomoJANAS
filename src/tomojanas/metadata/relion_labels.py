#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: tomojanas/metadata/relion_labels.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
"""
Canonical STAR column orders for tomoJANAS output.

RELION-compatible blocks use ``_rln*`` tags ONLY (so they can be read by
RELION). tomoJANAS-specific blocks use ``_tomoJANAS*`` tags and are kept in
SEPARATE blocks — never mixed into a RELION loop.
"""

from __future__ import annotations

# --------------------------------------------------------------------------- #
# tomograms.star :: data_global (RELION-compatible)
# --------------------------------------------------------------------------- #
TOMOGRAMS_GLOBAL_COLUMNS = [
    "_rlnTomoName",
    "_rlnTomoTiltSeriesStarFile",
    "_rlnTomoSizeX",
    "_rlnTomoSizeY",
    "_rlnTomoSizeZ",
    "_rlnTomoHand",
    "_rlnMicrographOriginalPixelSize",
    "_rlnTomoTiltSeriesPixelSize",
    "_rlnTomoTomogramBinning",
    "_rlnTomoReconstructedTomogram",
    "_rlnTomoDenoisedTomogram",
    "_rlnTomoFrameCount",
    "_rlnEtomoDirectiveFile",
    "_rlnOpticsGroupName",
    "_rlnVoltage",
    "_rlnSphericalAberration",
    "_rlnAmplitudeContrast",
]

# Minimum required tags for a tomogram to be considered valid.
TOMOGRAMS_REQUIRED = [
    "_rlnTomoName",
    "_rlnTomoTiltSeriesStarFile",
    "_rlnTomoSizeX",
    "_rlnTomoSizeY",
    "_rlnTomoSizeZ",
    "_rlnTomoTiltSeriesPixelSize",
    "_rlnTomoTomogramBinning",
    "_rlnTomoReconstructedTomogram",
]

# tomograms.star :: data_tomoJANAS_tomogram_sources (tomoJANAS-specific)
TOMOGRAMS_SOURCE_COLUMNS = [
    "_tomoJANASTomoName",
    "_tomoJANASImodDir",
    "_tomoJANASRawStack",
    "_tomoJANASMdocFile",
    "_tomoJANASAliStack",
    "_tomoJANASRecTomogram",
    "_tomoJANASXfFile",
    "_tomoJANASTltFile",
    "_tomoJANASNewstCom",
    "_tomoJANASTiltCom",
    "_tomoJANASMicrographReference",
    "_tomoJANASCtfSource",
    "_tomoJANASCtfPremultiplied",
    "_tomoJANASImportOffsetX",
    "_tomoJANASImportOffsetY",
    "_tomoJANASImportOffsetZ",
    "_tomoJANASRelionGeometryStatus",
]

# --------------------------------------------------------------------------- #
# tilt_series/<tomo>.star :: data_<tomo> (RELION-compatible, one row per tilt)
# --------------------------------------------------------------------------- #
TILT_SERIES_COLUMNS = [
    "_rlnMicrographName",
    "_rlnMicrographMovieName",
    "_rlnMicrographMetadata",
    "_rlnTomoTiltMovieFrameCount",
    "_rlnTomoTiltMovieIndex",
    "_rlnTomoNominalStageTiltAngle",
    "_rlnTomoNominalTiltAxisAngle",
    "_rlnTomoNominalDefocus",
    "_rlnMicrographPreExposure",
    "_rlnDefocusU",
    "_rlnDefocusV",
    "_rlnDefocusAngle",
    "_rlnCtfAstigmatism",
    "_rlnCtfFigureOfMerit",
    "_rlnCtfMaxResolution",
    "_rlnCtfScalefactor",
    "_rlnPhaseShift",
    "_rlnTomoXTilt",
    "_rlnTomoYTilt",
    "_rlnTomoZRot",
    "_rlnTomoXShiftAngst",
    "_rlnTomoYShiftAngst",
    "_rlnTomoProjX",
    "_rlnTomoProjY",
    "_rlnTomoProjZ",
    "_rlnTomoProjW",
]

# tilt_series/<tomo>.star :: data_tomoJANAS_tilt_mapping (tomoJANAS-specific)
TILT_MAPPING_COLUMNS = [
    "_tomoJANASTiltIndex",
    "_tomoJANASAliStack",
    "_tomoJANASAliSlice",
    "_tomoJANASRawStack",
    "_tomoJANASRawSlice",
    "_tomoJANASXfA11",
    "_tomoJANASXfA12",
    "_tomoJANASXfA21",
    "_tomoJANASXfA22",
    "_tomoJANASXfDX",
    "_tomoJANASXfDY",
    "_tomoJANASXfDirection",
    "_tomoJANASXfDirectionStatus",
]

# tilt_series/<tomo>.star :: data_tomoJANAS_projection_matrices (4x4 world->image,
# stored as scalar columns when the STAR writer cannot emit vector-valued fields)
PROJECTION_MATRIX_COLUMNS = (
    ["_tomoJANASTiltIndex"]
    + [f"_tomoJANASProj{r}{c}" for r in range(4) for c in range(4)]
)

# --------------------------------------------------------------------------- #
# CTF per-tilt RELION-compatible tags
# --------------------------------------------------------------------------- #
CTF_COLUMNS = [
    "_rlnDefocusU",
    "_rlnDefocusV",
    "_rlnDefocusAngle",
    "_rlnCtfAstigmatism",
    "_rlnCtfFigureOfMerit",
    "_rlnCtfMaxResolution",
    "_rlnCtfScalefactor",
    "_rlnPhaseShift",
    "_rlnCtfPowerSpectrum",
    "_rlnCtfImage",
]

# --------------------------------------------------------------------------- #
# optimisation_set.star :: data_optimisation_set
# --------------------------------------------------------------------------- #
OPTIMISATION_SET_COLUMNS = [
    "_rlnTomoTomogramsFile",
    "_rlnTomoParticlesFile",
]

# --------------------------------------------------------------------------- #
# particles_all.star and individual P*.star :: data_optics (RELION block)
# RELION tomography particle STAR files require BOTH data_optics and
# data_particles. data_optics must carry at least the pixel size + CTF basics.
# --------------------------------------------------------------------------- #
PARTICLE_OPTICS_COLUMNS = [
    "_rlnOpticsGroup",
    "_rlnOpticsGroupName",
    "_rlnTomoTiltSeriesPixelSize",
    "_rlnVoltage",
    "_rlnSphericalAberration",
    "_rlnAmplitudeContrast",
]

# --------------------------------------------------------------------------- #
# particles_all.star and individual P*.star :: data_particles (RELION block)
# --------------------------------------------------------------------------- #
PARTICLES_COLUMNS = [
    "_rlnTomoName",
    "_rlnTomoParticleName",
    "_rlnTomoParticleId",
    "_rlnCenteredCoordinateXAngst",
    "_rlnCenteredCoordinateYAngst",
    "_rlnCenteredCoordinateZAngst",
    "_rlnOriginXAngst",
    "_rlnOriginYAngst",
    "_rlnOriginZAngst",
    "_rlnTomoSubtomogramRot",
    "_rlnTomoSubtomogramTilt",
    "_rlnTomoSubtomogramPsi",
    "_rlnAngleRot",
    "_rlnAngleTilt",
    "_rlnAnglePsi",
    "_rlnOpticsGroup",
]

# individual P*.star :: data_tomoJANAS_particle_source
PARTICLE_SOURCE_COLUMNS = [
    "_tomoJANASParticleName",
    "_tomoJANASPickedVolume",
    "_tomoJANASPickedCoordinateX",
    "_tomoJANASPickedCoordinateY",
    "_tomoJANASPickedCoordinateZ",
    "_tomoJANASPickedCoordinateSystem",
    "_tomoJANASPickedIndexing",
    "_tomoJANASPickedAxisOrder",
    "_tomoJANASPickedSoftware",
]

# individual P*.star :: data_tomoJANAS_particle_roi
PARTICLE_ROI_COLUMNS = [
    "_tomoJANASParticleName",
    "_tomoJANASRoiShape3D",
    "_tomoJANASRoiShape2D",
    "_tomoJANASRoiRadiusAngst",
    "_tomoJANASRoiDiameterAngst",
    "_tomoJANASRoiPaddingAngst",
    "_tomoJANASStorageBoxShape3D",
    "_tomoJANASStorageBoxSizeVoxel",
    "_tomoJANASProjectionRadiusPixel",
    "_tomoJANASProjectionStorageBoxSizePixel",
]

# individual P*.star :: data_tomoJANAS_particle_projections
PARTICLE_PROJECTION_COLUMNS = [
    "_tomoJANASTiltIndex",
    "_tomoJANASStageTiltAngle",
    "_tomoJANASAlignedStack",
    "_tomoJANASAlignedSlice",
    "_tomoJANASAlignedCenterX",
    "_tomoJANASAlignedCenterY",
    "_tomoJANASAlignedRadiusPixel",
    "_tomoJANASAlignedCircleInsideFrame",
    "_tomoJANASRawStack",
    "_tomoJANASRawSlice",
    "_tomoJANASRawCenterX",
    "_tomoJANASRawCenterY",
    "_tomoJANASRawRadiusPixel",
    "_tomoJANASRawCircleInsideFrame",
    "_tomoJANASVisibleInTilt",
    "_tomoJANASProjectionStatus",
    "_tomoJANASDefocusU",
    "_tomoJANASDefocusV",
    "_tomoJANASDefocusAngle",
]

# Optional rec-crop columns (added when --write-rec-crops)
PARTICLE_REC_CROP_COLUMNS = [
    "_tomoJANASParticleRecPath",
    "_tomoJANASParticleRecSourceTomogram",
    "_tomoJANASParticleRecCropOriginX",
    "_tomoJANASParticleRecCropOriginY",
    "_tomoJANASParticleRecCropOriginZ",
    "_tomoJANASParticleRecBoxSizeVoxel",
    "_tomoJANASParticleRecRadiusVoxel",
    "_tomoJANASParticleRecRadiusAngst",
    "_tomoJANASParticleRecPixelSize",
    "_tomoJANASParticleRecStorageShape",
    "_tomoJANASParticleRecMaskShape",
]


# --------------------------------------------------------------------------- #
# Status vocabularies
# --------------------------------------------------------------------------- #
class RelionGeometryStatus:
    """Allowed values for _tomoJANASRelionGeometryStatus.

    The project is only marked RELION-extraction-ready when faithful 4x4
    world->image projection matrices are present in the RELION columns AND
    validated. Matrices built by the ported IMOD algorithm but not externally
    validated are flagged ``relion_imod_algorithm_ported``.
    """
    EXTRACT_READY = "relion_extract_ready"
    ALGORITHM_PORTED = "relion_imod_algorithm_ported"
    MATRICES_MISSING = "projection_matrices_missing"
    MATRICES_UNVALIDATED = "projection_matrices_unvalidated"
    MATRICES_APPROXIMATE = "projection_matrices_approximate"
    VALIDATION_FAILED = "relion_validation_failed"
    TOMOJANAS_VALID_ONLY = "tomojanas_valid_only"


class XfDirectionStatus:
    """Allowed values for _tomoJANASXfDirectionStatus."""
    EXPLICIT = "explicit"
    INFERRED_FROM_NEWST = "inferred_from_newst"
    ASSUMED_DEFAULT = "assumed_default"
    AMBIGUOUS = "ambiguous"


class ProjectionStatus:
    """Allowed values for _tomoJANASProjectionStatus (per particle/tilt)."""
    OK = "ok"
    OUTSIDE_FRAME = "outside_frame"
    MISSING = "projection_unavailable"
    APPROXIMATE = "approximate"
