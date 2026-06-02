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
]

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
