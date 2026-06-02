#
# File: assessParticles.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology
#

"""
Module: assessParticles.py
Implements particle assessment routines for the JANAS pipeline, including:
- Gaussian-based sampling of candidate particle subset sizes
- Optional CTF modulation of 2D particle images
- Block-wise per-particle scoring against reference volumes or half-maps
- Multiprocessing orchestration and integration of scoring results into STAR metadata
- Creation of differential stacks for diagnostic visualization

Aligns with:
  - “Particle subset selection and optimization” (Methods) for automaticParticleSubsetSelection
  - “Structural Cross Correlation Index (SCI)” (Methods) for per-particle scoring functions 
"""

# Standard library
import math
import operator
import timeit
import multiprocessing
from os import PathLike




# Third-party
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
from typing import Optional, Union

# Local
import janas.janas_core as janas_core
from janas import starHandler
from janas import janas_mapProcess

# eventually not used (implemented in C++):
# from skimage.metrics import structural_similarity as ssim
# from skimage.metrics import peak_signal_noise_ratio


def automaticParticleSubsetSelection(
    numParticleSubsetsSelected,
    expectedNumberOfPartilces,
    totalNumberOfPartilces,
    standardDeviation,
    seed=0,
):
    """
    Generate candidate particle subset sizes sampled from a skewed Gaussian distribution.

    Implements the adaptive sampling strategy described in the JANAS pipeline
    (Methods: Particle subset selection and optimization), exploring map quality
    across different particle counts to identify the optimal reconstruction threshold.

    Parameters
    ----------
    numParticleSubsetsSelected : int
        Number of subset sizes to sample.
    expectedNumberOfPartilces : int
        Mean of the Gaussian distribution (initial threshold).
    totalNumberOfPartilces : int
        Total particles in the input stack.
    standardDeviation : int
        Standard deviation controlling sampling spread.
    seed : int, optional
        Random seed for reproducibility (default: 0).

    Returns
    -------
    numpy.ndarray
        Sorted array of candidate subset sizes, always including the full stack size.
    """

    np.random.seed(seed)  # Set the random seed

    n_samples = numParticleSubsetsSelected
    total_samples = totalNumberOfPartilces
    mean = expectedNumberOfPartilces
    std_dev = standardDeviation  # You can adjust this value
    lower_bound = 20
    min_separation = total_samples / 5000  # Define minimum separation

    # Initialize an empty set to store the unique points
    unique_points = set()

    # Keep generating points until we have the desired number
    while len(unique_points) < n_samples:
        # Create a uniform distribution of points between 0 and 1
        uniform_point = np.random.uniform(0, 1)

        # Use the percent point function (inverse of CDF) to map the uniform point to a Gaussian distribution
        gaussian_point = stats.norm.ppf(uniform_point, mean, std_dev)

        # Clip the value to lie between 1 and the total number of samples
        gaussian_point = np.clip(gaussian_point, lower_bound, total_samples)

        # Round the value to the nearest integer
        gaussian_point = round(gaussian_point)

        # Check if the point is not too close to any existing point
        if all(math.fabs(gaussian_point - x) >= min_separation for x in unique_points):
            unique_points.add(gaussian_point)

    gaussian_points = np.array(list(unique_points))
    gaussian_points = np.sort(gaussian_points)
    if gaussian_points[-1] != totalNumberOfPartilces:
        gaussian_points = np.append(gaussian_points, totalNumberOfPartilces)
    return gaussian_points


def transformCtfImage(
    I,
    nx,
    ny,
    angpix,
    Voltage,
    DefocusU,
    DefocusV,
    DefocusAngle,
    SphericalAberration,
    CtfBfactor,
    PhaseShift,
    AmplitudeContrast,
    DetectorPixelSize,
    ctfMode="phaseflip",
    wiener_tau=0.1,
):
    """
    Apply CTF-based filtering to a single 2D particle image in Fourier space.

    Computes the microscope Contrast Transfer Function (CTF) using provided
    optical parameters and applies one of several modes in Fourier space:

      - 'modulate'  : multiply by the full CTF (current behaviour).
      - 'phaseflip': apply only the sign of the CTF (phase-flip).
      - 'wiener'   : apply a Wiener-like filter H / (H^2 + tau) with tau=wiener_tau.

    Parameters
    ----------
    I : list or array-like
        Flat list of pixel intensities for the particle image.
    nx, ny : int
        Image dimensions.
    angpix : float
        Pixel size (Å).
    Voltage : float
        Microscope accelerating voltage (kV).
    DefocusU, DefocusV : float
        Defocus values along U and V axes (Å).
    DefocusAngle : float
        Defocus astigmatism angle (degrees).
    SphericalAberration : float
        Spherical aberration coefficient (mm).
    CtfBfactor : float
        B-factor used in CTF.
    PhaseShift : float
        Phase plate shift (degrees; converted to radians internally).
    AmplitudeContrast : float
        Amplitude contrast ratio.
    DetectorPixelSize : float
        Detector pixel size (Å).
    ctfMode : {'modulate', 'phaseflip', 'wiener'}
        CTF application mode.
    wiener_tau : float
        Tau parameter for the Wiener filter, used when ctfMode='wiener'.

    Returns
    -------
    list
        Flat list of pixel values after applying the selected CTF mode.
        Returns zeros on size mismatch.
    """
    # PhaseShift is in degrees here; C++ converts to radians internally.
    ctfI = janas_core.CtfCenteredImage(
        nx,
        ny,
        angpix,
        SphericalAberration,
        Voltage,
        DefocusAngle,
        DefocusU,
        DefocusV,
        AmplitudeContrast,
        CtfBfactor,
        PhaseShift,
    )
    if (not len(ctfI) == nx * ny) or (not len(I) == nx * ny):
        return np.zeros(nx * ny).tolist()

    mapI = np.array(I).reshape((ny, nx))
    ctfI = np.array(ctfI).reshape((ny, nx))

    mapFFT = np.fft.fftshift(np.fft.fft2(mapI))

    mode = str(ctfMode).lower()
    if mode == "modulate":
        filt = ctfI
    elif mode == "phaseflip":
        ctf_sign = np.sign(ctfI)
        ctf_sign[ctf_sign == 0] = 1.0
        filt = ctf_sign
    elif mode == "wiener":
        # Simple Wiener-like filter: H / (H^2 + tau)
        H = ctfI
        denom = H * H + float(wiener_tau)
        filt = H / denom
    else:
        # Fallback to phase-flip if an unexpected mode is passed.
        ctf_sign = np.sign(ctfI)
        ctf_sign[ctf_sign == 0] = 1.0
        filt = ctf_sign

    mapIFFT = np.real(
        np.fft.ifft2(np.fft.ifftshift(mapFFT * filt))
    ).flatten().tolist()
    return mapIFFT



def scoreBlockParticles_original(
    particleIdxStart,
    particleIdxEnd,
    listScoresTags,
    mapI,
    maskI,
    MaskSubtractionAndI,   # NEW: 3D subtraction mask (or None / empty)
    sizeMap,
    angpix,
    pandasStarHeader,
    subsetCtfParameters,
    procnum,
    return_dict,
    ctfMode="phaseflip",
):
    """
    Compute per-particle similarity scores against a single reference volume, with
    optional RELION-like signal subtraction.

    If MaskSubtractionAndI is provided, for each particle:
      1. Project a 2D subtraction mask MSubI from MaskSubtractionAndI.
      2. Project the corresponding 2D signal RSubI from the 3D map.
      3. Compute a 2D subtraction mask in image space:
           mask_subtract_2d = MSubI * (1 - MI)
         so that only signal present in MaskSubtractionAndI but *outside* maskI is removed.
      4. If CTF parameters are available:
           - Build the CTF image once via janas_core.CtfCenteredImage.
           - CTF-modulate the subtraction projection (modulate mode).
           - Subtract this from the raw particle.
           - Apply the requested ctfMode (phaseflip, modulate, wiener) to the
             subtracted particle for scoring.
         If no CTF is provided:
           - Subtract RSubI * mask_subtract_2d directly in real space.

    Downstream MaskedImageComparison calls are unchanged, but see a modified particle I.
    """
    if procnum == 0:
        print("scoreBlockParticles_original, processes:", end=" ")
    print(procnum, end=" ", flush=True)
    is_last_proc = int(particleIdxEnd) == len(pandasStarHeader["_rlnImageName"])
    if is_last_proc:
        print()
        print("scoreBlockParticles")

    tmpScore = np.zeros(
        [int(particleIdxEnd) - int(particleIdxStart), len(listScoresTags) + 1]
    )

    counter = -1
    for ii in range(int(particleIdxStart), int(particleIdxEnd)):
        counter += 1
        if procnum == 0:
            percentage = 100 * ii / float(particleIdxEnd)
            print(
                "percentage= ",
                "{:10.4f}".format(percentage),
                "%  ",
                end="\r",
                flush=True,
            )

        tmpLine = pandasStarHeader["_rlnImageName"][ii]
        atPosition = tmpLine.find("@")
        imageNo = int(tmpLine[:atPosition])
        stackName = tmpLine[atPosition + 1 :]

        phi = pandasStarHeader.at[ii, "_rlnAngleRot"]
        theta = pandasStarHeader.at[ii, "_rlnAngleTilt"]
        psi = pandasStarHeader.at[ii, "_rlnAnglePsi"]
        tx = pandasStarHeader.at[ii, "_rlnOriginX"]
        ty = pandasStarHeader.at[ii, "_rlnOriginY"]

        # 1) Read particle
        I = janas_core.ReadMrcSlice(stackName, imageNo - 1)

        # 2) Project main scoring mask to 2D
        MI = janas_core.projectMask(
            maskI,
            sizeMap[0],
            sizeMap[1],
            sizeMap[2],
            phi,
            theta,
            psi,
            tx,
            ty,
            0,
            0.5,
        )

        # 3) Optional subtraction mask and corresponding reprojection
        MSubI = None
        RSubI = None
        use_subtraction = (
            MaskSubtractionAndI is not None and len(MaskSubtractionAndI) > 0
        )
        if use_subtraction:
            MSubI = janas_core.projectMask(
                MaskSubtractionAndI,
                sizeMap[0],
                sizeMap[1],
                sizeMap[2],
                phi,
                theta,
                psi,
                tx,
                ty,
                0,
                0.5,
            )
            RSubI = janas_core.projectMap_with2DMask(
                mapI,
                sizeMap[0],
                sizeMap[1],
                sizeMap[2],
                phi,
                theta,
                psi,
                tx,
                ty,
                0,
                MSubI,
                10,
            )

        # 4) Reprojection for scoring
        RI = janas_core.projectMap_with2DMask(
            mapI,
            sizeMap[0],
            sizeMap[1],
            sizeMap[2],
            phi,
            theta,
            psi,
            tx,
            ty,
            0,
            MI,
            10,
        )

        nx, ny = sizeMap[0], sizeMap[1]

        # 5) Apply subtraction (with or without CTF)
        if not subsetCtfParameters.empty:
            Voltage = subsetCtfParameters.iloc[counter]["_rlnVoltage"]
            DefocusU = subsetCtfParameters.iloc[counter]["_rlnDefocusU"]
            DefocusV = subsetCtfParameters.iloc[counter]["_rlnDefocusV"]
            DefocusAngle = subsetCtfParameters.iloc[counter]["_rlnDefocusAngle"]
            SphericalAberration = subsetCtfParameters.iloc[counter][
                "_rlnSphericalAberration"
            ]
            CtfBfactor = subsetCtfParameters.iloc[counter]["_rlnCtfBfactor"]
            PhaseShift = subsetCtfParameters.iloc[counter]["_rlnPhaseShift"]
            AmplitudeContrast = subsetCtfParameters.iloc[counter][
                "_rlnAmplitudeContrast"
            ]
            DetectorPixelSize = subsetCtfParameters.iloc[counter][
                "_rlnDetectorPixelSize"
            ]

            if use_subtraction:
                I2d = np.asarray(I, dtype=np.float32).reshape(ny, nx)
                MI_2d = np.asarray(MI, dtype=np.float32).reshape(ny, nx)
                MSub_2d = np.asarray(MSubI, dtype=np.float32).reshape(ny, nx)
                RSub_2d = np.asarray(RSubI, dtype=np.float32).reshape(ny, nx)

                # only subtract where in MaskSubtractionAndI but NOT in maskI
                mask_subtract_2d = MSub_2d * (1.0 - MI_2d)

                # Build CTF image ONCE
                ctfI = janas_core.CtfCenteredImage(
                    nx,
                    ny,
                    angpix,
                    SphericalAberration,
                    Voltage,
                    DefocusAngle,
                    DefocusU,
                    DefocusV,
                    AmplitudeContrast,
                    CtfBfactor,
                    PhaseShift,
                )
                if len(ctfI) != nx * ny:
                    I = np.zeros(nx * ny, dtype=np.float32).tolist()
                else:
                    ctf2d = np.asarray(ctfI, dtype=np.float32).reshape(ny, nx)

                    # 5a) RELION-like subtraction: CTF-modulate the projection to subtract
                    P2d = RSub_2d * mask_subtract_2d
                    P_fft = np.fft.fftshift(np.fft.fft2(P2d))
                    P_filt = P_fft * ctf2d  # "modulate" mode
                    P_ctf_2d = np.real(
                        np.fft.ifft2(np.fft.ifftshift(P_filt))
                    )

                    # Subtract from raw particle
                    I2d_sub = I2d - P_ctf_2d

                    # 5b) Apply requested ctfMode for scoring to the subtracted particle
                    I_fft = np.fft.fftshift(np.fft.fft2(I2d_sub))
                    mode = str(ctfMode).lower()
                    if mode == "modulate":
                        filt_I = ctf2d
                    elif mode == "phaseflip":
                        ctf_sign = np.sign(ctf2d)
                        ctf_sign[ctf_sign == 0] = 1.0
                        filt_I = ctf_sign
                    elif mode == "wiener":
                        H = ctf2d
                        denom = H * H + 0.1
                        filt_I = H / denom
                    else:
                        ctf_sign = np.sign(ctf2d)
                        ctf_sign[ctf_sign == 0] = 1.0
                        filt_I = ctf_sign

                    I2d_final = np.real(
                        np.fft.ifft2(np.fft.ifftshift(I_fft * filt_I))
                    )
                    I = I2d_final.flatten().tolist()
            else:
                # Original behaviour: just CTF-transform the particle for scoring
                I = transformCtfImage(
                    I,
                    nx,
                    ny,
                    angpix,
                    Voltage,
                    DefocusU,
                    DefocusV,
                    DefocusAngle,
                    SphericalAberration,
                    CtfBfactor,
                    PhaseShift,
                    AmplitudeContrast,
                    DetectorPixelSize,
                    ctfMode=ctfMode,
                )
        else:
            # No CTF → pure real-space subtraction if requested
            if use_subtraction:
                I2d = np.asarray(I, dtype=np.float32).reshape(ny, nx)
                MI_2d = np.asarray(MI, dtype=np.float32).reshape(ny, nx)
                MSub_2d = np.asarray(MSubI, dtype=np.float32).reshape(ny, nx)
                RSub_2d = np.asarray(RSubI, dtype=np.float32).reshape(ny, nx)
                mask_subtract_2d = MSub_2d * (1.0 - MI_2d)
                I2d_sub = I2d - RSub_2d * mask_subtract_2d
                I = I2d_sub.flatten().tolist()
            # else: leave I as read

        # 6) scoring
        tmpScore[counter][0] = ii
        for kk in range(0, len(listScoresTags)):
            comparisonMethod = (
                listScoresTags[kk].split("_")[2]
                if len(listScoresTags[kk].split("_")) > 2
                else "CC"
            )
            preprocessingMethod = (
                listScoresTags[kk].split("_")[3]
                if len(listScoresTags[kk].split("_")) > 3
                else "unprocessed"
            )
            sigmaBlur = (
                str(listScoresTags[kk].split("_")[4])
                if len(listScoresTags[kk].split("_")) > 4
                else "1"
            )

            if comparisonMethod in ("CC", "MI"):
                tmpScore[counter][kk + 1] = janas_core.MaskedImageComparison(
                    I,
                    RI,
                    MI,
                    sizeMap[0],
                    sizeMap[1],
                    1,
                    comparisonMethod,
                    preprocessingMethod,
                    sigmaBlur,
                )
            elif comparisonMethod == "SSIM":
                scoreSSIM = janas_core.MaskedImageComparison(
                    RI,
                    I,
                    MI,
                    sizeMap[0],
                    sizeMap[1],
                    1,
                    comparisonMethod,
                    preprocessingMethod,
                    sigmaBlur,
                )
                tmpScore[counter][kk + 1] = scoreSSIM
            elif comparisonMethod == "PSNR":
                scorePSNR = janas_core.MaskedImageComparison(
                    RI,
                    I,
                    MI,
                    sizeMap[0],
                    sizeMap[1],
                    1,
                    comparisonMethod,
                    preprocessingMethod,
                    sigmaBlur,
                )
                tmpScore[counter][kk + 1] = scorePSNR
            elif comparisonMethod == "SCI":
                RI_2d = np.reshape(RI, [sizeMap[1], sizeMap[0]])
                I_fft = np.fft.fftn(np.reshape(I, [sizeMap[1], sizeMap[0]]))
                RI_fft = np.fft.fftn(RI_2d)
                I_abs_fft = np.abs(I_fft) + 1e-7
                RI_abs_fft = np.abs(RI_fft) + 1e-7
                ampAvg = 0.5 * (I_abs_fft + RI_abs_fft)
                I_out = (
                    np.fft.ifftn(I_fft * (ampAvg / I_abs_fft)).real.flatten().tolist()
                )
                RI_out = (
                    np.fft.ifftn(RI_fft * (ampAvg / RI_abs_fft)).real.flatten().tolist()
                )
                scoreSCI = janas_core.MaskedImageComparison(
                    RI_out,
                    I_out,
                    MI,
                    sizeMap[0],
                    sizeMap[1],
                    1,
                    comparisonMethod,
                    preprocessingMethod,
                    sigmaBlur,
                )
                tmpScore[counter][kk + 1] = scoreSCI

    return_dict[procnum] = tmpScore
    return tmpScore


def scoreBlockParticles_halfMaps(
    particleIdxStart,
    particleIdxEnd,
    listScoresTags,
    mapI_A,
    mapI_B,
    maskI,
    MaskSubtractionAndI,   # NEW
    sizeMap,
    angpix,
    pandasStarHeader,
    subsetCtfParameters,
    procnum,
    return_dict,
    ctfMode="phaseflip",
):
    """
    As scoreBlockParticles_original but choosing mapI_A or mapI_B per particle
    according to _rlnRandomSubset, and optionally subtracting a map region
    defined by MaskSubtractionAndI in a RELION-like way with CTF.
    """
    if procnum == 0:
        print("scoreBlockParticles_halfMaps, processes:", end=" ")
    print(procnum, end=" ", flush=True)

    is_last_proc = int(particleIdxEnd) == len(pandasStarHeader["_rlnImageName"])
    if is_last_proc:
        print()
        print("scoreBlockParticles_halfMaps")

    tmpScore = np.zeros(
        [int(particleIdxEnd) - int(particleIdxStart), len(listScoresTags) + 1]
    )

    counter = -1
    nx, ny = sizeMap[0], sizeMap[1]

    for ii in range(int(particleIdxStart), int(particleIdxEnd)):
        counter += 1
        if procnum == 0:
            percentage = 100 * ii / float(particleIdxEnd)
            print(
                "percentage= ",
                "{:10.4f}".format(percentage),
                "%  ",
                end="\r",
                flush=True,
            )

        tmpLine = pandasStarHeader["_rlnImageName"][ii]
        atPosition = tmpLine.find("@")
        imageNo = int(tmpLine[:atPosition])
        stackName = tmpLine[atPosition + 1 :]

        phi = pandasStarHeader.at[ii, "_rlnAngleRot"]
        theta = pandasStarHeader.at[ii, "_rlnAngleTilt"]
        psi = pandasStarHeader.at[ii, "_rlnAnglePsi"]
        tx = pandasStarHeader.at[ii, "_rlnOriginX"]
        ty = pandasStarHeader.at[ii, "_rlnOriginY"]
        hm = int(pandasStarHeader.at[ii, "_rlnRandomSubset"])

        I = janas_core.ReadMrcSlice(stackName, imageNo - 1)

        MI = janas_core.projectMask(
            maskI,
            sizeMap[0],
            sizeMap[1],
            sizeMap[2],
            phi,
            theta,
            psi,
            tx,
            ty,
            0,
            0.5,
        )

        use_subtraction = (
            MaskSubtractionAndI is not None and len(MaskSubtractionAndI) > 0
        )
        MSubI = None
        RSubI = None
        if use_subtraction:
            MSubI = janas_core.projectMask(
                MaskSubtractionAndI,
                sizeMap[0],
                sizeMap[1],
                sizeMap[2],
                phi,
                theta,
                psi,
                tx,
                ty,
                0,
                0.5,
            )

        # Choose half-map
        if hm == 1:
            RI = janas_core.projectMap_with2DMask(
                mapI_A,
                sizeMap[0],
                sizeMap[1],
                sizeMap[2],
                phi,
                theta,
                psi,
                tx,
                ty,
                0,
                MI,
                10,
            )
            if use_subtraction:
                RSubI = janas_core.projectMap_with2DMask(
                    mapI_A,
                    sizeMap[0],
                    sizeMap[1],
                    sizeMap[2],
                    phi,
                    theta,
                    psi,
                    tx,
                    ty,
                    0,
                    MSubI,
                    10,
                )
        else:
            RI = janas_core.projectMap_with2DMask(
                mapI_B,
                sizeMap[0],
                sizeMap[1],
                sizeMap[2],
                phi,
                theta,
                psi,
                tx,
                ty,
                0,
                MI,
                10,
            )
            if use_subtraction:
                RSubI = janas_core.projectMap_with2DMask(
                    mapI_B,
                    sizeMap[0],
                    sizeMap[1],
                    sizeMap[2],
                    phi,
                    theta,
                    psi,
                    tx,
                    ty,
                    0,
                    MSubI,
                    10,
                )

        # CTF + subtraction
        if not subsetCtfParameters.empty:
            Voltage = subsetCtfParameters.iloc[counter]["_rlnVoltage"]
            DefocusU = subsetCtfParameters.iloc[counter]["_rlnDefocusU"]
            DefocusV = subsetCtfParameters.iloc[counter]["_rlnDefocusV"]
            DefocusAngle = subsetCtfParameters.iloc[counter]["_rlnDefocusAngle"]
            SphericalAberration = subsetCtfParameters.iloc[counter][
                "_rlnSphericalAberration"
            ]
            CtfBfactor = subsetCtfParameters.iloc[counter]["_rlnCtfBfactor"]
            PhaseShift = subsetCtfParameters.iloc[counter]["_rlnPhaseShift"]
            AmplitudeContrast = subsetCtfParameters.iloc[counter][
                "_rlnAmplitudeContrast"
            ]
            DetectorPixelSize = subsetCtfParameters.iloc[counter][
                "_rlnDetectorPixelSize"
            ]

            if use_subtraction:
                I2d = np.asarray(I, dtype=np.float32).reshape(ny, nx)
                MI_2d = np.asarray(MI, dtype=np.float32).reshape(ny, nx)
                MSub_2d = np.asarray(MSubI, dtype=np.float32).reshape(ny, nx)
                RSub_2d = np.asarray(RSubI, dtype=np.float32).reshape(ny, nx)
                mask_subtract_2d = MSub_2d * (1.0 - MI_2d)

                ctfI = janas_core.CtfCenteredImage(
                    nx,
                    ny,
                    angpix,
                    SphericalAberration,
                    Voltage,
                    DefocusAngle,
                    DefocusU,
                    DefocusV,
                    AmplitudeContrast,
                    CtfBfactor,
                    PhaseShift,
                )
                if len(ctfI) != nx * ny:
                    I = np.zeros(nx * ny, dtype=np.float32).tolist()
                else:
                    ctf2d = np.asarray(ctfI, dtype=np.float32).reshape(ny, nx)

                    P2d = RSub_2d * mask_subtract_2d
                    P_fft = np.fft.fftshift(np.fft.fft2(P2d))
                    P_filt = P_fft * ctf2d
                    P_ctf_2d = np.real(
                        np.fft.ifft2(np.fft.ifftshift(P_filt))
                    )

                    I2d_sub = I2d - P_ctf_2d

                    I_fft = np.fft.fftshift(np.fft.fft2(I2d_sub))
                    mode = str(ctfMode).lower()
                    if mode == "modulate":
                        filt_I = ctf2d
                    elif mode == "phaseflip":
                        ctf_sign = np.sign(ctf2d)
                        ctf_sign[ctf_sign == 0] = 1.0
                        filt_I = ctf_sign
                    elif mode == "wiener":
                        H = ctf2d
                        denom = H * H + 0.1
                        filt_I = H / denom
                    else:
                        ctf_sign = np.sign(ctf2d)
                        ctf_sign[ctf_sign == 0] = 1.0
                        filt_I = ctf_sign
                    I2d_final = np.real(
                        np.fft.ifft2(np.fft.ifftshift(I_fft * filt_I))
                    )
                    I = I2d_final.flatten().tolist()
            else:
                I = transformCtfImage(
                    I,
                    nx,
                    ny,
                    angpix,
                    Voltage,
                    DefocusU,
                    DefocusV,
                    DefocusAngle,
                    SphericalAberration,
                    CtfBfactor,
                    PhaseShift,
                    AmplitudeContrast,
                    DetectorPixelSize,
                    ctfMode=ctfMode,
                )
        else:
            if use_subtraction:
                I2d = np.asarray(I, dtype=np.float32).reshape(ny, nx)
                MI_2d = np.asarray(MI, dtype=np.float32).reshape(ny, nx)
                MSub_2d = np.asarray(MSubI, dtype=np.float32).reshape(ny, nx)
                RSub_2d = np.asarray(RSubI, dtype=np.float32).reshape(ny, nx)
                mask_subtract_2d = MSub_2d * (1.0 - MI_2d)
                I2d_sub = I2d - RSub_2d * mask_subtract_2d
                I = I2d_sub.flatten().tolist()

        tmpScore[counter][0] = ii
        for kk in range(0, len(listScoresTags)):
            comparisonMethod = (
                listScoresTags[kk].split("_")[2]
                if len(listScoresTags[kk].split("_")) > 2
                else "CC"
            )
            preprocessingMethod = (
                listScoresTags[kk].split("_")[3]
                if len(listScoresTags[kk].split("_")) > 3
                else "unprocessed"
            )
            sigmaBlur = (
                str(listScoresTags[kk].split("_")[4])
                if len(listScoresTags[kk].split("_")) > 4
                else "1"
            )

            if comparisonMethod in ("CC", "MI"):
                tmpScore[counter][kk + 1] = janas_core.MaskedImageComparison(
                    I,
                    RI,
                    MI,
                    sizeMap[0],
                    sizeMap[1],
                    1,
                    comparisonMethod,
                    preprocessingMethod,
                    sigmaBlur,
                )
            elif comparisonMethod == "SSIM":
                scoreSSIM = janas_core.MaskedImageComparison(
                    RI,
                    I,
                    MI,
                    sizeMap[0],
                    sizeMap[1],
                    1,
                    comparisonMethod,
                    preprocessingMethod,
                    sigmaBlur,
                )
                tmpScore[counter][kk + 1] = scoreSSIM
            elif comparisonMethod == "PSNR":
                scorePSNR = janas_core.MaskedImageComparison(
                    RI,
                    I,
                    MI,
                    sizeMap[0],
                    sizeMap[1],
                    1,
                    comparisonMethod,
                    preprocessingMethod,
                    sigmaBlur,
                )
                tmpScore[counter][kk + 1] = scorePSNR
            elif comparisonMethod == "SCI":
                RI_2d = np.reshape(RI, [sizeMap[1], sizeMap[0]])
                I_fft = np.fft.fftn(np.reshape(I, [sizeMap[1], sizeMap[0]]))
                RI_fft = np.fft.fftn(RI_2d)
                I_abs_fft = np.abs(I_fft) + 1e-7
                RI_abs_fft = np.abs(RI_fft) + 1e-7
                ampAvg = 0.5 * (I_abs_fft + RI_abs_fft)
                I_out = (
                    np.fft.ifftn(I_fft * (ampAvg / I_abs_fft)).real.flatten().tolist()
                )
                RI_out = (
                    np.fft.ifftn(RI_fft * (ampAvg / RI_abs_fft)).real.flatten().tolist()
                )
                scoreSCI = janas_core.MaskedImageComparison(
                    RI_out,
                    I_out,
                    MI,
                    sizeMap[0],
                    sizeMap[1],
                    1,
                    comparisonMethod,
                    preprocessingMethod,
                    sigmaBlur,
                )
                tmpScore[counter][kk + 1] = scoreSCI

    return_dict[procnum] = tmpScore
    return tmpScore


def ParticleVsReprojectionScores(
    particlesStarFile: PathLike,
    scoredParticlesStarFile: PathLike,
    referenceMap: PathLike,
    referenceMask: PathLike,
    angpix,
    listScoresTags="",
    numProcesses=1,
    numViews=[50, 200, 400],
    doCTF=False,
    ctfMode="phaseflip",
    referenceSubtractionMask: Union[PathLike, None] = None,
):
    """
    Orchestrate multi-process per-particle scoring against a single reference map.

    Parses the input STAR file (Relion v3.1 vs older), optionally merges CTF
    parameters, divides particles across processes, runs scoreBlockParticles_original,
    aggregates and sorts scores, and writes new score columns into the output STAR.

    Parameters
    ----------
    particlesStarFile : PathLike
        Input STAR filename.
    scoredParticlesStarFile : PathLike
        Output STAR filename with added score columns.
    referenceMap : PathLike
        Reference MRC volume for reprojections.
    referenceMask : PathLike
        MRC volume mask.
    angpix : float
        Pixel size (Å).
    listScoresTags : list of str or ""
        Score identifiers to compute; defaults to ['_janas_SCI__1'].
    numProcesses : int
        Number of parallel worker processes.
    numViews : list of int
        Placeholder for future orientation coverage logic.
    doCTF : bool
        Whether to apply per-particle CTF correction.
    """
    if listScoresTags == "":
        listScoresTags = ["_janas_SCI__1"]
    print("listScoresTags=", listScoresTags)

    mapI = janas_core.ReadMRC(referenceMap)
    sizeMap = janas_core.sizeMRC(referenceMap)
    maskI = janas_core.ReadMRC(referenceMask)
    #    sizeMask=janas_core.sizeMRC(referenceMask)

    MaskSubtractionAndI = None
    if referenceSubtractionMask is not None:
        MaskSubtractionAndI = janas_core.ReadMRC(referenceSubtractionMask)

    version = starHandler.infoStarFile(particlesStarFile)[2]
    if version == "relion_v31":
        coordinatesFULL = starHandler.readColumns(
            particlesStarFile,
            [
                "_rlnRandomSubset",
                "_rlnImageName",
                "_rlnAngleRot",
                "_rlnAngleTilt",
                "_rlnAnglePsi",
                "_rlnOriginXAngst",
                "_rlnOriginYAngst",
                "_rlnOpticsGroup",
            ],
        )
        idx = [x for x in range(0, len(coordinatesFULL))]
        coordinatesFULL["idx"] = idx
        coordinatesDataOptics = starHandler.dataOptics(particlesStarFile)[
            ["_rlnOpticsGroup", "_rlnImagePixelSize"]
        ]
        coordinates = pd.merge(
            coordinatesFULL, coordinatesDataOptics, on=["_rlnOpticsGroup"]
        ).sort_values(["idx"])
        coordinates["_rlnOriginX"] = (
            coordinates["_rlnOriginXAngst"] / coordinates["_rlnImagePixelSize"]
        )
        coordinates["_rlnOriginY"] = (
            coordinates["_rlnOriginYAngst"] / coordinates["_rlnImagePixelSize"]
        )
        coordinates = coordinates.drop(["_rlnOpticsGroup"], axis=1).reindex()
        coordinates = coordinates.set_index("idx")
        coordinates = coordinates.drop(["_rlnOriginXAngst", "_rlnOriginYAngst"], axis=1)
    else:
        coordinates = starHandler.readColumns(
            particlesStarFile,
            [
                "_rlnRandomSubset",
                "_rlnImageName",
                "_rlnAngleRot",
                "_rlnAngleTilt",
                "_rlnAnglePsi",
                "_rlnOriginX",
                "_rlnOriginY",
            ],
        )

    ctfParameters = pd.DataFrame([])
    if doCTF:
        columns = starHandler.header_columns(particlesStarFile)
        print(len(coordinates))
        if not "_rlnPhaseShift" in columns:
            PhaseShift = pd.DataFrame(np.zeros(len(coordinates)))
        else:
            PhaseShift = starHandler.readColumns(particlesStarFile, ["_rlnPhaseShift"])

        # print (columns)

        print("doing CTF...")
        if version == "relion_v31":
            print("READ parameters")
            parametersFULL = starHandler.readColumns(
                particlesStarFile,
                [
                    "_rlnImageName",
                    "_rlnDefocusU",
                    "_rlnDefocusV",
                    "_rlnDefocusAngle",
                    "_rlnOpticsGroup",
                    "_rlnCtfBfactor",
                ],
            )
            print("GOT the parameters")
            idx = [x for x in range(0, len(parametersFULL))]
            parametersFULL["idx"] = idx
            parametersDataOptics = starHandler.dataOptics(particlesStarFile)[
                [
                    "_rlnImagePixelSize",
                    "_rlnVoltage",
                    "_rlnAmplitudeContrast",
                    "_rlnSphericalAberration",
                    "_rlnOpticsGroup",
                ]
            ]
            ctfParameters = pd.merge(
                parametersFULL, parametersDataOptics, on=["_rlnOpticsGroup"]
            ).sort_values(["idx"])
            ctfParameters = ctfParameters.drop(["_rlnOpticsGroup"], axis=1).reindex()
            ctfParameters = ctfParameters.set_index("idx")
            ctfParameters.rename(
                columns={"_rlnImagePixelSize": "_rlnDetectorPixelSize"}, inplace=True
            )
        else:
            (ctfParameters) = starHandler.readColumns(
                particlesStarFile,
                [
                    "_rlnImageName",
                    "_rlnDefocusU",
                    "_rlnDefocusV",
                    "_rlnDefocusAngle",
                    "_rlnDetectorPixelSize",
                    "_rlnVoltage",
                    "_rlnAmplitudeContrast",
                    "_rlnSphericalAberration",
                    "_rlnCtfBfactor",
                ],
            )
        ctfParameters["_rlnPhaseShift"] = PhaseShift
    numParticles = len(coordinates["_rlnImageName"])
    print("num particles=", numParticles, "  num Processes=", numProcesses, ", loading process:")
    start_time = timeit.default_timer()


    if numProcesses > multiprocessing.cpu_count():
        numProcesses = multiprocessing.cpu_count()
    elif numProcesses < 1:
        numProcesses = 1
    if numProcesses > numParticles:
        numProcesses = numParticles

    # if numProcesses == 1:
    # manager = multiprocessing.Manager()
    # return_dict = manager.dict()
    # OK   subsetCtfParameters=ctfParameters[0:numParticles]
    # OK   #print (subsetCtfParameters)
    # OK    scoresOut = scoreBlockParticles(0,numParticles,listScoresTags, mapI, maskI, sizeMap, angpix, coordinates, subsetCtfParameters, 0,0)
    # OK   df_full_scores = pd.DataFrame(scoresOut)
    # OK    df_full_scores = df_full_scores.iloc[: , 1:]
    # OK    df_full_scores.columns=listScoresTags
    # OK    starHandler.removeColumnsTagsStartingWith(particlesStarFile, scoredParticlesStarFile, "_emprove_")
    # OK    starHandler.addDataframeColumns(scoredParticlesStarFile, scoredParticlesStarFile, listScoresTags, df_full_scores)
    # OK    return

    blockSize = np.floor(numParticles / numProcesses)
    idxMatrix = np.zeros([numProcesses, 2])
    for mm in range(0, numProcesses):
        idxMatrix[mm][0] = mm * blockSize
        idxMatrix[mm][1] = (mm + 1) * (blockSize)
        if mm == numProcesses - 1:
            idxMatrix[mm][1] = int(numParticles)

    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    jobs = []
    for ii in range(0, len(idxMatrix)):
        # print ('**************\n***********\ndebug=',idxMatrix[ii][0]),'   ',int(idxMatrix[ii][1])
        subsetCtfParameters = ctfParameters[
            int(idxMatrix[ii][0]) : int(idxMatrix[ii][1])
        ]
        p = multiprocessing.Process(
            target=scoreBlockParticles_original,
            args=(
                idxMatrix[ii][0],
                idxMatrix[ii][1],
                listScoresTags,
                mapI,
                maskI,
                MaskSubtractionAndI,   # NEW
                sizeMap,
                angpix,
                coordinates,
                subsetCtfParameters,
                ii,
                return_dict,
                ctfMode,
            ),
        )
        jobs.append(p)
        p.start()

    for proc in jobs:
        proc.join()
    print(f"\nTotal workers: {len(jobs)}\n")

    scoresOut = []
    for ii in range(0, len(idxMatrix)):
        scoresOut = scoresOut + return_dict.values()[ii].tolist()
    scoresOut = np.array(scoresOut).reshape([numParticles, len(listScoresTags) + 1])
    listScoreOut = (
        np.array(sorted(scoresOut, key=operator.itemgetter(0), reverse=False))
        .reshape([numParticles, len(listScoresTags) + 1])[:, 1:]
        .tolist()
    )
    df_full_scores = pd.DataFrame(data=listScoreOut, columns=listScoresTags)
    headers_columns = starHandler.header_columns(particlesStarFile)
    columnsToRemove = [item for item in headers_columns if item in listScoresTags]
    starHandler.removeColumns(
        particlesStarFile, scoredParticlesStarFile, columnsToRemove
    )
    starHandler.addDataframeColumns(
        scoredParticlesStarFile, scoredParticlesStarFile, listScoresTags, df_full_scores
    )
    elapsed = timeit.default_timer() - start_time
    hours = int(elapsed // 3600)
    minutes = int((elapsed % 3600) // 60)
    seconds = elapsed % 60
    if elapsed > 0:
        speed = numParticles / elapsed
        print(f"[ScoreParticles] Elapsed time  "
              f"{hours:02d}h:{minutes:02d}m:{seconds:05.2f}s  "
              f"({speed:.2f} particles/s)")
    else:
        print("[ScoreParticles] Elapsed time: < 1e-6 s")

def ParticleVsReprojectionScores_HalfMaps(
    particlesStarFile: PathLike,
    scoredParticlesStarFile: PathLike,
    referenceMap1: PathLike,
    referenceMap2: PathLike,
    referenceMask: PathLike,
    angpix,
    listScoresTags="",
    numProcesses=1,
    numViews=[50, 200, 400],
    doCTF=False,
    ctfMode="phaseflip",
    referenceSubtractionMask: Union[PathLike, None] = None,
):
    """
    Orchestrate multi-process per-particle scoring against two gold-standard half-maps.

    Similar to ParticleVsReprojectionScores but splits scoring between two
    half-maps to maintain gold-standard validation.

    Parameters
    ----------
    particlesStarFile : PathLike
        Input STAR filename.
    scoredParticlesStarFile : PathLike
        Output STAR filename with added score columns.
    referenceMap1, referenceMap2 : PathLike
        Filenames for the two half-map volumes.
    referenceMask : PathLike
        MRC volume mask.
    angpix : float
        Pixel size (Å).
    listScoresTags : list of str or ""
        Score identifiers; defaults to ['_janas_SCI__1'].
    numProcesses : int
        Number of parallel workers.
    numViews : list of int
        Placeholder parameter.
    doCTF : bool
        Whether to apply CTF correction.
    """
    if listScoresTags == "":
        listScoresTags = ["_janas_SCI__1"]
    print("listScoresTags=", listScoresTags)

    mapI_A = janas_core.ReadMRC(referenceMap1)
    mapI_B = janas_core.ReadMRC(referenceMap2)
    sizeMap = janas_core.sizeMRC(referenceMap1)
    maskI = janas_core.ReadMRC(referenceMask)
    #    sizeMask=janas_core.sizeMRC(referenceMask)

    # --- Amplitude normalisation of half-maps
    amplitudeNormalisation=False
    if amplitudeNormalisation:
        try:
            volA_norm, volB_norm, Bnorm = janas_mapProcess.normalize_amplitudes(
                referenceMap1,
                referenceMap2,
            )
            # Replace flat map lists with normalised volumes for all subsequent projections
            mapI_A = volA_norm.astype(np.float32, copy=False).ravel().tolist()
            mapI_B = volB_norm.astype(np.float32, copy=False).ravel().tolist()
        except Exception as e:
            print(
                "[ParticleVsReprojectionScores_HalfMaps] WARNING: "
                f"amplitude normalisation failed: {e}"
            )


    MaskSubtractionAndI = None
    if referenceSubtractionMask is not None:
        MaskSubtractionAndI = janas_core.ReadMRC(referenceSubtractionMask)

    version = starHandler.infoStarFile(particlesStarFile)[2]
    if version == "relion_v31":
        coordinatesFULL = starHandler.readColumns(
            particlesStarFile,
            [
                "_rlnRandomSubset",
                "_rlnImageName",
                "_rlnAngleRot",
                "_rlnAngleTilt",
                "_rlnAnglePsi",
                "_rlnOriginXAngst",
                "_rlnOriginYAngst",
                "_rlnOpticsGroup",
            ],
        )
        idx = [x for x in range(0, len(coordinatesFULL))]
        coordinatesFULL["idx"] = idx
        coordinatesDataOptics = starHandler.dataOptics(particlesStarFile)[
            ["_rlnOpticsGroup", "_rlnImagePixelSize"]
        ]
        coordinates = pd.merge(
            coordinatesFULL, coordinatesDataOptics, on=["_rlnOpticsGroup"]
        ).sort_values(["idx"])
        coordinates["_rlnOriginX"] = (
            coordinates["_rlnOriginXAngst"] / coordinates["_rlnImagePixelSize"]
        )
        coordinates["_rlnOriginY"] = (
            coordinates["_rlnOriginYAngst"] / coordinates["_rlnImagePixelSize"]
        )
        coordinates = coordinates.drop(["_rlnOpticsGroup"], axis=1).reindex()
        coordinates = coordinates.set_index("idx")
        coordinates = coordinates.drop(["_rlnOriginXAngst", "_rlnOriginYAngst"], axis=1)
    else:
        coordinates = starHandler.readColumns(
            particlesStarFile,
            [
                "_rlnRandomSubset",
                "_rlnImageName",
                "_rlnAngleRot",
                "_rlnAngleTilt",
                "_rlnAnglePsi",
                "_rlnOriginX",
                "_rlnOriginY",
            ],
        )

    ctfParameters = pd.DataFrame([])
    if doCTF:
        columns = starHandler.header_columns(particlesStarFile)
        print(len(coordinates))
        if not "_rlnPhaseShift" in columns:
            PhaseShift = pd.DataFrame(np.zeros(len(coordinates)))
        else:
            PhaseShift = starHandler.readColumns(particlesStarFile, ["_rlnPhaseShift"])

        # print (columns)

        print("doing CTF...")
        if version == "relion_v31":
            print("READ parameters")
            parametersFULL = starHandler.readColumns(
                particlesStarFile,
                [
                    "_rlnImageName",
                    "_rlnDefocusU",
                    "_rlnDefocusV",
                    "_rlnDefocusAngle",
                    "_rlnOpticsGroup",
                    "_rlnCtfBfactor",
                ],
            )
            print("GOT the parameters")
            idx = [x for x in range(0, len(parametersFULL))]
            parametersFULL["idx"] = idx
            parametersDataOptics = starHandler.dataOptics(particlesStarFile)[
                [
                    "_rlnImagePixelSize",
                    "_rlnVoltage",
                    "_rlnAmplitudeContrast",
                    "_rlnSphericalAberration",
                    "_rlnOpticsGroup",
                ]
            ]
            ctfParameters = pd.merge(
                parametersFULL, parametersDataOptics, on=["_rlnOpticsGroup"]
            ).sort_values(["idx"])
            ctfParameters = ctfParameters.drop(["_rlnOpticsGroup"], axis=1).reindex()
            ctfParameters = ctfParameters.set_index("idx")
            ctfParameters.rename(
                columns={"_rlnImagePixelSize": "_rlnDetectorPixelSize"}, inplace=True
            )
        else:
            (ctfParameters) = starHandler.readColumns(
                particlesStarFile,
                [
                    "_rlnImageName",
                    "_rlnDefocusU",
                    "_rlnDefocusV",
                    "_rlnDefocusAngle",
                    "_rlnDetectorPixelSize",
                    "_rlnVoltage",
                    "_rlnAmplitudeContrast",
                    "_rlnSphericalAberration",
                    "_rlnCtfBfactor",
                ],
            )
        ctfParameters["_rlnPhaseShift"] = PhaseShift
    numParticles = len(coordinates["_rlnImageName"])
    print("num particles=", numParticles, "  num Processes=", numProcesses, ", loading process:")
    start_time = timeit.default_timer()


    if numProcesses > multiprocessing.cpu_count():
        numProcesses = multiprocessing.cpu_count()
    elif numProcesses < 1:
        numProcesses = 1
    if numProcesses > numParticles:
        numProcesses = numParticles

    # if numProcesses == 1:
    # manager = multiprocessing.Manager()
    # return_dict = manager.dict()
    # OK   subsetCtfParameters=ctfParameters[0:numParticles]
    # OK   #print (subsetCtfParameters)
    # OK    scoresOut = scoreBlockParticles(0,numParticles,listScoresTags, mapI, maskI, sizeMap, angpix, coordinates, subsetCtfParameters, 0,0)
    # OK   df_full_scores = pd.DataFrame(scoresOut)
    # OK    df_full_scores = df_full_scores.iloc[: , 1:]
    # OK    df_full_scores.columns=listScoresTags
    # OK    starHandler.removeColumnsTagsStartingWith(particlesStarFile, scoredParticlesStarFile, "_emprove_")
    # OK    starHandler.addDataframeColumns(scoredParticlesStarFile, scoredParticlesStarFile, listScoresTags, df_full_scores)
    # OK    return

    blockSize = np.floor(numParticles / numProcesses)
    idxMatrix = np.zeros([numProcesses, 2])
    for mm in range(0, numProcesses):
        idxMatrix[mm][0] = mm * blockSize
        idxMatrix[mm][1] = (mm + 1) * (blockSize)
        if mm == numProcesses - 1:
            idxMatrix[mm][1] = int(numParticles)

    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    jobs = []
    for ii in range(0, len(idxMatrix)):
        # print ('**************\n***********\ndebug=',idxMatrix[ii][0]),'   ',int(idxMatrix[ii][1])
        subsetCtfParameters = ctfParameters[
            int(idxMatrix[ii][0]) : int(idxMatrix[ii][1])
        ]
        p = multiprocessing.Process(
            target=scoreBlockParticles_halfMaps,
            args=(
                idxMatrix[ii][0],
                idxMatrix[ii][1],
                listScoresTags,
                mapI_A,
                mapI_B,
                maskI,
                MaskSubtractionAndI,  # NEW
                sizeMap,
                angpix,
                coordinates,
                subsetCtfParameters,
                ii,
                return_dict,
                ctfMode,
            ),
        )
        jobs.append(p)
        p.start()

    for proc in jobs:
        proc.join()
    print(f"\nTotal workers: {len(jobs)}\n")


    scoresOut = []
    for ii in range(0, len(idxMatrix)):
        scoresOut = scoresOut + return_dict.values()[ii].tolist()
    scoresOut = np.array(scoresOut).reshape([numParticles, len(listScoresTags) + 1])
    listScoreOut = (
        np.array(sorted(scoresOut, key=operator.itemgetter(0), reverse=False))
        .reshape([numParticles, len(listScoresTags) + 1])[:, 1:]
        .tolist()
    )
    df_full_scores = pd.DataFrame(data=listScoreOut, columns=listScoresTags)
    headers_columns = starHandler.header_columns(particlesStarFile)
    columnsToRemove = [item for item in headers_columns if item in listScoresTags]
    starHandler.removeColumns(
        particlesStarFile, scoredParticlesStarFile, columnsToRemove
    )
    starHandler.addDataframeColumns(
        scoredParticlesStarFile, scoredParticlesStarFile, listScoresTags, df_full_scores
    )
    elapsed = timeit.default_timer() - start_time
    hours = int(elapsed // 3600)
    minutes = int((elapsed % 3600) // 60)
    seconds = elapsed % 60
    if elapsed > 0:
        speed = numParticles / elapsed
        print(f"[ScoreParticles] Elapsed time  "
              f"{hours:02d}h:{minutes:02d}m:{seconds:05.2f}s  "
              f"({speed:.2f} particles/s)")
    else:
        print("[ScoreParticles] Elapsed time: < 1e-6 s")



def createDiffStack(
    particlesStarFile: PathLike,
    outputBasename,
    referenceMap: PathLike,
    referenceMask: PathLike,
    useCTF=False,
):
    """
    Build a differential MRC stack by weighting each particle image by its local SDIM score.

    For each particle:
      1. Reads the image slice and corresponding reprojection/mask.
      2. Normalizes both image and reprojection.
      3. Computes a per-pixel SDIM score via janas_core.SDIM.
      4. Multiplies the raw image by mask and score, zeros negatives.
      5. Writes the processed slice into a new .mrcs stack.
      6. Updates the STAR file to reference the new differential stack.

    Parameters
    ----------
    particlesStarFile : PathLike
        Input STAR filename.
    outputBasename : str
        Base name for generated .mrcs and .star files.
    referenceMap : PathLike
        Reference MRC volume.
    referenceMask : PathLike
        MRC volume mask.
    useCTF : bool, unused
        Placeholder for future CTF logic.
    """
    mapI = janas_core.ReadMRC(referenceMap)
    sizeMap = janas_core.sizeMRC(referenceMap)
    maskI = janas_core.ReadMRC(referenceMask)
    #    sizeMask=janas_core.sizeMRC(referenceMask)
    imageNameTag = "_rlnImageName"
    imageNames = starHandler.readColumns(particlesStarFile, [imageNameTag])
    tmpLine = imageNames[imageNameTag][0]
    stackName = tmpLine[tmpLine.find("@") + 1 :]
    sizeI = janas_core.sizeMRC(stackName)
    janas_core.WriteEmptyMRC(
        outputBasename + ".mrcs", sizeI[0], sizeI[1], len(imageNames[imageNameTag])
    )

    version = starHandler.infoStarFile(particlesStarFile)[2]
    if version == "relion_v31":
        coordinatesFULL = starHandler.readColumns(
            particlesStarFile,
            [
                "_rlnRandomSubset",
                "_rlnImageName",
                "_rlnAngleRot",
                "_rlnAngleTilt",
                "_rlnAnglePsi",
                "_rlnOriginXAngst",
                "_rlnOriginYAngst",
                "_rlnOpticsGroup",
            ],
        )
        idx = [x for x in range(0, len(coordinatesFULL))]
        coordinatesFULL["idx"] = idx
        coordinatesDataOptics = starHandler.dataOptics(particlesStarFile)[
            ["_rlnOpticsGroup", "_rlnImagePixelSize"]
        ]
        coordinates = pd.merge(
            coordinatesFULL, coordinatesDataOptics, on=["_rlnOpticsGroup"]
        ).sort_values(["idx"])
        coordinates["_rlnOriginX"] = (
            coordinates["_rlnOriginXAngst"] / coordinates["_rlnImagePixelSize"]
        )
        coordinates["_rlnOriginY"] = (
            coordinates["_rlnOriginYAngst"] / coordinates["_rlnImagePixelSize"]
        )
        coordinates = coordinates.drop(["_rlnOpticsGroup"], axis=1).reindex()
        coordinates = coordinates.set_index("idx")
        coordinates = coordinates.drop(["_rlnOriginXAngst", "_rlnOriginYAngst"], axis=1)
    else:
        coordinates = starHandler.readColumns(
            particlesStarFile,
            [
                "_rlnRandomSubset",
                "_rlnImageName",
                "_rlnAngleRot",
                "_rlnAngleTilt",
                "_rlnAnglePsi",
                "_rlnOriginX",
                "_rlnOriginY",
            ],
        )

    nx = sizeI[0]
    ny = sizeI[1]
    outImageNames = []
    numIterations = len(imageNames[imageNameTag])
    for ii in range(0, len(coordinates["_rlnImageName"])):
        print("iteration ", ii, " out of ", numIterations, end="\r")

        # print('iteration ',ii,' out of ', numIterations)
        tmpLine = coordinates["_rlnImageName"][ii]
        atPosition = tmpLine.find("@")
        imageNo = int(tmpLine[:atPosition])
        stackName = tmpLine[atPosition + 1 :]
        newFileNames = tmpLine[:atPosition] + "@"

        phi = coordinates.at[ii, "_rlnAngleRot"]
        theta = coordinates.at[ii, "_rlnAngleTilt"]
        psi = coordinates.at[ii, "_rlnAnglePsi"]
        tx = coordinates.at[ii, "_rlnOriginX"]
        ty = coordinates.at[ii, "_rlnOriginY"]
        outImageNames.append(str(tmpLine[:atPosition] + "@" + outputBasename + ".mrcs"))
        I = janas_core.ReadMrcSlice(stackName, imageNo - 1)
        MI = janas_core.projectMask(
            maskI, sizeMap[0], sizeMap[1], sizeMap[2], phi, theta, psi, tx, ty, 0, 0.5
        )

        #RI = janas_core.projectMap( mapI, sizeMap[0], sizeMap[1], sizeMap[2], phi, theta, psi, tx, ty, 0, maskI, 0.5)
        RI = janas_core.projectMap_with2DMask( mapI, sizeMap[0], sizeMap[1], sizeMap[2], phi, theta, psi, tx, ty, 0, MI, 10)


        I_out = np.array(I)
        meanValue1 = np.mean(I_out)
        stdValue1 = np.std(I_out)
        if stdValue1 > 0:
            I_out = ((I_out - meanValue1) / (stdValue1)).tolist()

        RI_out = np.array(RI)
        meanValue2 = np.mean(RI_out)
        stdValue2 = np.std(RI_out)
        if stdValue2 > 0:
            RI_out = ((RI_out - meanValue2) / (stdValue2)).tolist()

        # janas_core.WriteMRC( I_out, "I1_preproc1.mrc",  shapeImg1[0], shapeImg1[1],  1, 1)
        # janas_core.WriteMRC( RI_out, "I2_preproc1.mrc", shapeImg1[0], shapeImg1[1],  1, 1)

        # shapeImg1 = shapeImg1 + (1,1,)
        scoreImage = janas_core.SDIM(RI_out, I_out, nx, ny, 1, str(1), MI)

        Iout = np.array(I) * np.array(MI) * np.array(scoreImage)
        Iout[Iout < 0.0] = 0.0
        Iout = Iout.tolist()
        janas_core.ReplaceMrcSlice(
            Iout, outputBasename + ".mrcs", sizeI[0], sizeI[1], ii
        )

    inputStarFile = starHandler.readStar(particlesStarFile)
    inputStarFile[imageNameTag] = outImageNames
    starHandler.writeDataframeToStar(
        particlesStarFile, outputBasename + ".star", inputStarFile
    )
    # Report elapsed time and speed
    elapsed = timeit.default_timer() - start_time
