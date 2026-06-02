# File: utils.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology

"""
Module: utils.py

Provides utilities for MRC map handling, including pixel spacing extraction,
CTF stack creation, frequency separation, map differencing, and reprojection
correction. Integrates with janas_core and starHandler for reading/writing
MRC and STAR files.
"""

# Standard library
import os
from os import PathLike, path

# Third-party
import numpy as np
import pandas as pd
import json
import logging
import re
from typing import List, Optional, Tuple, Dict, Any
from scipy.spatial.transform import Rotation
from numpy.fft import fftn, ifftn, fftfreq
import struct


# Local
import janas.janas_core as janas_core
from janas import starHandler



# ---------- Minimal MRC header/data helpers (little-endian CCP4/MRC) ----------

def _read_mrc_header(fp) -> Dict[str, Any]:
    """Read 1024-byte CCP4/MRC header at current file position (assumes LE)."""
    h = fp.read(1024)
    if len(h) < 1024:
        raise ValueError("File too short for an MRC header.")
    nx, ny, nz, mode = struct.unpack("<4i", h[0:16])
    nxstart, nystart, nzstart = struct.unpack("<3i", h[16:28])
    mx, my, mz = struct.unpack("<3i", h[28:40])
    xlen, ylen, zlen = struct.unpack("<3f", h[40:52])
    alpha, beta, gamma = struct.unpack("<3f", h[52:64])
    mapc, mapr, maps = struct.unpack("<3i", h[64:76])
    amin, amax, amean = struct.unpack("<3f", h[76:88])
    ispg, nsymbt = struct.unpack("<2i", h[88:96])

    # ORIGIN floats: standard offsets 196–208 (words 49–51)
    try:
        ox, oy, oz = struct.unpack("<3f", h[196:208])
    except struct.error:
        ox = oy = oz = 0.0

    # Stamp + machine: bytes 208–216
    stamp = h[208:212]
    mach  = h[212:216]

    # guard against bogus or missing unit cell entries
    if mx <= 0: mx = nx
    if my <= 0: my = ny
    if mz <= 0: mz = nz

    pix_x = (xlen / mx) if mx else 1.0
    pix_y = (ylen / my) if my else pix_x
    pix_z = (zlen / mz) if mz else pix_x

    return {
        "nx": nx, "ny": ny, "nz": nz, "mode": mode,
        "nxstart": nxstart, "nystart": nystart, "nzstart": nzstart,
        "mx": mx, "my": my, "mz": mz,
        "xlen": float(xlen), "ylen": float(ylen), "zlen": float(zlen),
        "alpha": float(alpha), "beta": float(beta), "gamma": float(gamma),
        "mapc": mapc, "mapr": mapr, "maps": maps,
        "amin": float(amin), "amax": float(amax), "amean": float(amean),
        "ispg": ispg, "nsymbt": nsymbt,
        "origin_x": float(ox), "origin_y": float(oy), "origin_z": float(oz),
        "pixel_x": float(pix_x), "pixel_y": float(pix_y), "pixel_z": float(pix_z),
        "stamp": stamp, "mach": mach,
        "_raw": h  # retain original if ever needed
    }

def _dtype_from_mode(mode: int):
    # Standard modes: 0=int8, 1=int16, 2=float32, 3=complex int16, 4=complex float32, 6=uint16
    if mode == 2: return np.float32
    if mode == 1: return np.int16
    if mode == 0: return np.int8
    if mode == 6: return np.uint16
    # Non-standard but seen in practice
    if mode == 12: return np.float16  # some stacks use 12 to mean float16
    # Complex types not supported by this copier:
    if mode in (3, 4):
        raise ValueError(f"Unsupported complex MRC mode {mode}.")
    raise ValueError(f"Unsupported MRC mode {mode}.")

def _write_mrc_header(fp, nx: int, ny: int, nz: int, apix: float,
                      vmin: float, vmax: float, vmean: float, vrms: float,
                      mode: int = 2,
                      origin_angs: Tuple[float, float, float] = (0.0, 0.0, 0.0)):
    """Write a 1024-byte header (LE)."""
    header = bytearray(1024)
    def pw(word_off: int, fmt: str, *vals):
        struct.pack_into("<" + fmt, header, word_off * 4, *vals)
    mx, my, mz = nx, ny, nz
    # words 0–3
    pw(0, "i", nx); pw(1, "i", ny); pw(2, "i", nz); pw(3, "i", mode)
    # starts
    pw(4, "i", 0); pw(5, "i", 0); pw(6, "i", 0)
    # mx,my,mz
    pw(7, "i", mx); pw(8, "i", my); pw(9, "i", mz)
    # cell dims
    pw(10, "f", mx * apix); pw(11, "f", my * apix); pw(12, "f", mz * apix)
    # angles
    pw(13, "f", 90.0); pw(14, "f", 90.0); pw(15, "f", 90.0)
    # axis order
    pw(16, "i", 1); pw(17, "i", 2); pw(18, "i", 3)
    # stats
    pw(19, "f", vmin); pw(20, "f", vmax); pw(21, "f", vmean)
    # spacegroup & nsymbt
    pw(22, "i", 0); pw(23, "i", 0)
    # origin (words 49–51)
    ox, oy, oz = origin_angs
    pw(49, "f", float(ox)); pw(50, "f", float(oy)); pw(51, "f", float(oz))
    # 'MAP ' + machine stamp 'DA\x00\x00' (little-endian)
    header[208:212] = b"MAP "
    header[212:216] = b"DA\x00\x00"
    # rms and labels count
    pw(54, "f", vrms); pw(55, "i", 0)
    # write out
    fp.seek(0)
    fp.write(header)

def _slice_offset_bytes(hdr: Dict[str, Any], slice_index0: int) -> int:
    """Byte offset of a given z-slice plane (0-based) in a stack file."""
    nx, ny = hdr["nx"], hdr["ny"]
    nsymbt  = hdr["nsymbt"]
    elsz    = np.dtype(_dtype_from_mode(hdr["mode"])).itemsize
    return 1024 + nsymbt + slice_index0 * nx * ny * elsz


# get_pixel_spacing
def get_MRC_map_pixel_spacing(filename):
    """
    Parse an MRC header to extract pixel dimensions and compute angstrom-per-voxel.

    Args:
        filename: Path to the .mrc file.
    Returns:
        (apix_x, apix_y, apix_z): Pixel spacing along each axis.
    """
    import struct

    with open(filename, "rb") as f:
        # Read the number of columns, rows, and sections (bytes 0-11)
        nx = struct.unpack("i", f.read(4))[0]
        ny = struct.unpack("i", f.read(4))[0]
        nz = struct.unpack("i", f.read(4))[0]

        # Seek to the cell dimensions (bytes 40-51)
        f.seek(40)
        x_dim = struct.unpack("f", f.read(4))[0]
        y_dim = struct.unpack("f", f.read(4))[0]
        z_dim = struct.unpack("f", f.read(4))[0]

    # Calculate the actual pixel spacing by dividing dimensions by the number
    # of voxels
    apix_x = x_dim / nx
    apix_y = y_dim / ny
    apix_z = z_dim / nz
    return apix_x, apix_y, apix_z


#############################
# FSC functions
def compute_fsc_3d(vol1: np.ndarray, vol2: np.ndarray, angpix: float, n_bins: Optional[int] = None, eps: float = 1e-12):
    """
    Compute Fourier Shell Correlation between two 3D volumes.

    vol1, vol2: arrays shaped (Z, Y, X)
    angpix: pixel size in Å/pixel
    n_bins: number of concentric shells (default: half of smallest dimension)
    returns: freqs (Å^-1), fsc (dimensionless)
    """
    if vol1.shape != vol2.shape:
        raise ValueError(f"Shape mismatch: {vol1.shape} vs {vol2.shape}")

    # Forward FFTs
    F1 = fftn(vol1)
    F2 = fftn(vol2)

    nz, ny, nx = vol1.shape
    # Spatial frequency axes in Å^-1 (fftfreq uses sample spacing 'd' in Å)
    fx = fftfreq(nx, d=angpix)
    fy = fftfreq(ny, d=angpix)
    fz = fftfreq(nz, d=angpix)
    # Build frequency radius grid matching (Z, Y, X)
    FZ, FY, FX = np.meshgrid(fz, fy, fx, indexing='ij')
    radii = np.sqrt(FX*FX + FY*FY + FZ*FZ).ravel()

    S1 = F1.ravel()
    S2 = F2.ravel()

    if n_bins is None:
        n_bins = int(min(vol1.shape) // 2)
        n_bins = max(n_bins, 1)

    max_freq = radii.max()
    bins = np.linspace(0.0, max_freq, n_bins + 1, dtype=np.float64)
    idx = np.digitize(radii, bins) - 1

    fsc = np.zeros(n_bins, dtype=np.float64)
    freqs = 0.5 * (bins[:-1] + bins[1:])

    for i in range(n_bins):
        m = (idx == i)
        if not np.any(m):
            fsc[i] = np.nan
            continue
        num = np.sum(S1[m] * np.conj(S2[m]))
        den = np.sqrt(np.sum(np.abs(S1[m])**2) * np.sum(np.abs(S2[m])**2)) + eps
        fsc[i] = np.real(num / den)

    return freqs, fsc

def compute_fsc_2d(img1: np.ndarray, img2: np.ndarray, angpix: float,
                   n_bins: Optional[int] = None, eps: float = 1e-12):
    """
    Compute Fourier Ring Correlation (2D analogue of FSC) between two 2D images.

    img1, img2: arrays shaped (Y, X)
    angpix: pixel size in Å/pixel
    n_bins: number of concentric rings (default: half of smallest image dimension)
    returns: freqs (Å^-1), frc (dimensionless)
    """
    if img1.shape != img2.shape:
        raise ValueError(f"Shape mismatch: {img1.shape} vs {img2.shape}")
    if img1.ndim != 2 or img2.ndim != 2:
        raise ValueError(f"compute_fsc_2d expects 2D arrays; got {img1.ndim}D and {img2.ndim}D")

    # Forward FFTs
    F1 = fftn(img1)
    F2 = fftn(img2)

    ny, nx = img1.shape
    # Spatial frequency axes in Å^-1
    fx = fftfreq(nx, d=angpix)
    fy = fftfreq(ny, d=angpix)

    # Build frequency radius grid matching (Y, X)
    FY, FX = np.meshgrid(fy, fx, indexing='ij')
    radii = np.sqrt(FX*FX + FY*FY).ravel()

    S1 = F1.ravel()
    S2 = F2.ravel()

    if n_bins is None:
        n_bins = int(min(img1.shape) // 2)
        n_bins = max(n_bins, 1)

    max_freq = radii.max()
    bins = np.linspace(0.0, max_freq, n_bins + 1, dtype=np.float64)
    idx = np.digitize(radii, bins) - 1

    frc = np.zeros(n_bins, dtype=np.float64)
    freqs = 0.5 * (bins[:-1] + bins[1:])

    for i in range(n_bins):
        m = (idx == i)
        if not np.any(m):
            frc[i] = np.nan
            continue
        num = np.sum(S1[m] * np.conj(S2[m]))
        den = np.sqrt(np.sum(np.abs(S1[m])**2) * np.sum(np.abs(S2[m])**2)) + eps
        frc[i] = np.real(num / den)

    return freqs, frc

def find_resolution_at_threshold(freqs: np.ndarray, fsc_vals: np.ndarray, threshold: float):
    """
    Linear interpolation of the first downward crossing of FSC relative to 'threshold'.
    Returns (f_cross [Å^-1], resolution [Å]); returns (None, None) if no crossing.
    """
    # Ensure finite values for robust crossing detection
    fsc = np.array(fsc_vals, dtype=float)
    finite = np.isfinite(fsc)
    if not np.any(finite):
        return None, None
    fsc[~finite] = -1.0

    above = fsc >= threshold
    crossings = np.where(above[:-1] & (~above[1:]))[0]
    if crossings.size == 0:
        return None, None

    i = int(crossings[0])
    f1, f2 = float(freqs[i]), float(freqs[i+1])
    y1, y2 = float(fsc[i]), float(fsc[i+1])
    if y1 == y2:
        return None, None
    frac = (y1 - threshold) / (y1 - y2)
    f_cross = f1 + frac * (f2 - f1)
    if f_cross <= 0:
        return None, None
    return f_cross, 1.0 / f_cross

#############################
# create ctfStack
def ctfStack(
        particlesStarFile: PathLike,
        outputStackBasename,
        saveOriginal=False):
    imageNameTag = "_rlnImageName"
    imageNames = starHandler.readColumns(particlesStarFile, [imageNameTag])

    #    if '_rlnPhaseShift' in imageNames:
    #    	print ('we have _rlnPhaseShift')
    #    else:
    #    	print ('we DO NOT have _rlnPhaseShift')

    # CTFParameters ctf_parameters(rlnSphAberration[ii], rlnVoltage[ii],
    # rlnDefocusAngle[ii], rlnDefocusU[ii], rlnDefocusV[ii],
    # rlnAmplitudeContrast[ii], 0, 0);
    version = starHandler.infoStarFile(particlesStarFile)[2]
    if version == "relion_v31":

        if "_rlnPhaseShift" in imageNames:
            parametersFULL = starHandler.readColumns(
                particlesStarFile,
                [
                    "_rlnImageName",
                    "_rlnDefocusU",
                    "_rlnDefocusV",
                    "_rlnDefocusAngle",
                    "_rlnOpticsGroup",
                    "_rlnCtfBfactor",
                    "_rlnPhaseShift",
                ],
            )
        else:
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
        ctfParameters = ctfParameters.drop(
            ["_rlnOpticsGroup"], axis=1).reindex()
        ctfParameters = ctfParameters.set_index("idx")
        ctfParameters.rename(
            columns={
                "_rlnImagePixelSize": "_rlnDetectorPixelSize"},
            inplace=True)
    else:
        if "_rlnPhaseShift" in imageNames:
            ctfParameters = starHandler.readColumns(
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
                    "_rlnPhaseShift",
                ],
            )
        else:
            ctfParameters = starHandler.readColumns(
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

    if "_rlnPhaseShift" not in imageNames:
        ctfParameters["_rlnPhaseShift"] = np.zeros(len(ctfParameters))

    tmpLine = imageNames[imageNameTag][0]
    stackName = tmpLine[tmpLine.find("@") + 1:]
    sizeI = janas_core.sizeMRC(stackName)
    janas_core.WriteEmptyMRC(
        outputStackBasename + ".mrcs", sizeI[0], sizeI[1], len(imageNames[imageNameTag])
    )
    nx = sizeI[0]
    ny = sizeI[1]
    outImageNames = []
    numIterations = len(imageNames[imageNameTag])
    for ii in range(0, len(imageNames[imageNameTag])):
        print("iteration ", ii, " out of ", numIterations, end="\r")
        tmpLine = imageNames[imageNameTag][ii]
        atPosition = tmpLine.find("@")
        imageNo = int(tmpLine[:atPosition])
        stackName = tmpLine[atPosition + 1:]
        # print (ii,' >> ',imageNo,' >> ',imageNames[imageNameTag][ii])
        angpix = ctfParameters.at[ii, "_rlnDetectorPixelSize"]
        SphericalAberration = ctfParameters.at[ii, "_rlnSphericalAberration"]
        voltage = ctfParameters.at[ii, "_rlnVoltage"]
        DefocusAngle = ctfParameters.at[ii, "_rlnDefocusAngle"]
        DefocusU = ctfParameters.at[ii, "_rlnDefocusU"]
        DefocusV = ctfParameters.at[ii, "_rlnDefocusV"]
        AmplitudeContrast = ctfParameters.at[ii, "_rlnAmplitudeContrast"]
        Bfac = ctfParameters.at[ii, "_rlnCtfBfactor"]
        phase_shift = ctfParameters.at[ii, "_rlnPhaseShift"]
        ctfI = janas_core.CtfCenteredImage(
            nx,
            ny,
            angpix,
            SphericalAberration,
            voltage,
            DefocusAngle,
            DefocusU,
            DefocusV,
            AmplitudeContrast,
            Bfac,
            phase_shift,
        )
        ctfI = np.array(ctfI).reshape((sizeI[1], sizeI[0]))
        mapI = np.array(janas_core.ReadMrcSlice(stackName, imageNo - 1))
        mapI = mapI.reshape((sizeI[1], sizeI[0]))
        mapFFT = np.fft.fftshift(np.fft.fft2(mapI))
        mapIFFT = np.real(
            np.fft.ifft2(
                np.fft.ifftshift(
                    mapFFT *
                    ctfI))).flatten()
        janas_core.ReplaceMrcSlice(
            mapIFFT.tolist(),
            outputStackBasename +
            ".mrcs",
            sizeI[0],
            sizeI[1],
            ii)
        outImageNames.append(
            str(str(ii + 1).zfill(7) + "@" + outputStackBasename + ".mrcs")
        )

    #    for ii in range(0,len(outImageNames)):
    #        print (outImageNames[ii])

    originalStarDataframe = starHandler.readStar(particlesStarFile)
    if saveOriginal:
        print("keeping the original")
        # originalStarDataframe.rename(columns = {'_rlnImageName':'_original_rlnImageName'}, inplace = True)
        originalStarDataframe["_original_rlnImageName"] = imageNames
    originalStarDataframe["_rlnImageName"] = outImageNames
    # print (originalStarDataframe)
    starHandler.writeDataframeToStar(
        particlesStarFile, outputStackBasename + ".star", originalStarDataframe
    )


def separate_frequencies_by_two_thresholds(
        I, apix, low_threshold, high_threshold):
    """
    Separate frequencies in a Fourier-transformed image into three components:
    frequencies lower than low_threshold, between low_threshold and high_threshold,
    and higher than high_threshold, assuming I is already correctly shaped.

    Parameters:
    - I: 2D numpy array, the input image, already reshaped.
    - apix: Pixel size in Angstroms.
    - low_threshold, high_threshold: The low and high threshold values in Angstroms.

    Returns:
    - Three numpy arrays corresponding to the separated spatial components.
    """

    # Fourier transform
    I_fft = np.fft.fftn(I)
    I_fft_shifted = np.fft.fftshift(I_fft)

    # Frequency grid
    rows, cols = I.shape
    kx = np.fft.fftfreq(cols, d=apix).reshape(-1, 1)
    ky = np.fft.fftfreq(rows, d=apix).reshape(1, -1)
    kx, ky = np.fft.fftshift(kx), np.fft.fftshift(ky)
    frequency_grid = np.sqrt(kx ** 2 + ky ** 2)

    # Masks for frequency separation
    low_freq_mask = frequency_grid < (1 / low_threshold)
    mid_freq_mask = (frequency_grid >= (1 / low_threshold)) & (
        frequency_grid < (1 / high_threshold)
    )
    high_freq_mask = frequency_grid >= (1 / high_threshold)

    # Apply masks
    low_freq_component = np.fft.ifftshift(I_fft_shifted * low_freq_mask)
    mid_freq_component = np.fft.ifftshift(I_fft_shifted * mid_freq_mask)
    high_freq_component = np.fft.ifftshift(I_fft_shifted * high_freq_mask)

    # Inverse Fourier transform to spatial domain
    low_freq_image = np.fft.ifftn(low_freq_component).real
    mid_freq_image = np.fft.ifftn(mid_freq_component).real
    high_freq_image = np.fft.ifftn(high_freq_component).real

    return low_freq_image, mid_freq_image, high_freq_image
    # return low_freq_image


def separate_frequencies_by_threshold(I, apix, threshold):
    """
    Separate frequencies in a Fourier-transformed image into two components:
    frequencies lower than the threshold, and frequencies higher than the threshold,
    assuming I is already correctly shaped.

    Parameters:
    - I: 2D numpy array, the input image.
    - apix: Pixel size in Angstroms.
    - threshold: The threshold value in Angstroms.

    Returns:
    - Two numpy arrays corresponding to the separated spatial components: the low-frequency
      component and the high-frequency component.
    """

    # Fourier transform
    I_fft = np.fft.fftn(I)
    I_fft_shifted = np.fft.fftshift(I_fft)

    # Frequency grid
    rows, cols = I.shape
    kx = np.fft.fftfreq(cols, d=apix).reshape(-1, 1)
    ky = np.fft.fftfreq(rows, d=apix).reshape(1, -1)
    kx, ky = np.fft.fftshift(kx), np.fft.fftshift(ky)
    frequency_grid = np.sqrt(kx ** 2 + ky ** 2)

    # Masks for frequency separation
    low_freq_mask = frequency_grid < (1 / threshold)
    # mid_freq_mask = (frequency_grid >= (1 / threshold)) & (frequency_grid < (1 / high_threshold))
    high_freq_mask = frequency_grid >= (1 / threshold)

    # Apply masks
    low_freq_component = np.fft.ifftshift(I_fft_shifted * low_freq_mask)
    # mid_freq_component = np.fft.ifftshift(I_fft_shifted * mid_freq_mask)
    high_freq_component = np.fft.ifftshift(I_fft_shifted * high_freq_mask)

    # Inverse Fourier transform to spatial domain
    low_freq_image = np.fft.ifftn(low_freq_component).real
    # mid_freq_image = np.fft.ifftn(mid_freq_component).real
    high_freq_image = np.fft.ifftn(high_freq_component).real

    return low_freq_image, high_freq_image
    # return low_freq_image


#############################
# Read MRC
def readMRC(mrcFile: PathLike):
    print("read MRC file ", mrcFile)
    if not path.exists(mrcFile):
        print("ERROR: file ", mrcFile, " does not exists")
        return None
    if (
        not mrcFile.endswith(".mrc")
        or not mrcFile.endswith(".mrcs")
        or not mrcFile.endswith(".st")
        or not not mrcFile.endswith(".rawst")
    ):
        print(
            "ERROR: file ",
            mrcFile,
            " does not have a recognized extension. If you sure it is a mrc file, just rename it as .mrc",
        )
        return None

    # if (particlesStarFile.endswith(suffix))


#############################
# maps Difference
def mapsDifference(mrcFileList, mapout):
    for ii in range(0, len(mrcFileList)):
        mrcImage = mrcFileList[ii]
        if not path.exists(mrcImage):
            print("ERROR: file ", mrcImage, " does not exists")
            return None
        # if not  mrcImage.endswith('.mrc') or not mrcImage.endswith('.mrcs') or not mrcImage.endswith('.st') or not not mrcImage.endswith('.rawst'):
        #    print ('ERROR: file ',mrcImage,' does not have a recognized extension. If you sure it is a mrc file, just rename it as .mrc')
        #    return None
    AvgMap = np.array(janas_core.ReadMRC(mrcFileList[0]))
    sizeMap = janas_core.sizeMRC(mrcFileList[0])
    for ii in range(1, len(mrcFileList)):
        AvgMap = AvgMap + np.array(janas_core.ReadMRC(mrcFileList[ii]))
    AvgMap = AvgMap / len(mrcFileList)
    janas_core.WriteMRC(
        AvgMap.tolist(), "tmpAvg.mrc", sizeMap[0], sizeMap[1], sizeMap[2], 1
    )
    janas_core.replaceMrcHeader(mrcFileList[0], "tmpAvg.mrc")

    # amplitude comp
    AvgMap = np.reshape(AvgMap, sizeMap)
    # I_fft=np.fft.fftn(np.reshape(AvgMap,sizeMap))
    # I_abs_fft=np.abs(I_fft)+0.0000001
    diffMap = np.zeros(sizeMap)
    for ii in range(0, len(mrcFileList)):
        mapTmp = np.reshape(
            np.array(
                janas_core.ReadMRC(
                    mrcFileList[ii])),
            sizeMap)
        # tmpI_fft=np.fft.fftn(np.reshape(mapTmp,sizeMap))
        # tmpI=np.fft.ifftn(I_fft*(I_abs_fft/abs(tmpI_fft))).real
        # diffMap=diffMap+np.square(tmpI-AvgMap)
        diffMap = diffMap + np.square(mapTmp - AvgMap)

    # diffMap=diffMap/(len(mrcFileList)+np.square(AvgMap))

    # I_out=np.fft.ifftn(I_fft*(RI_abs_fft/I_abs_fft)).real.flatten().tolist()
    janas_core.WriteMRC(
        diffMap.flatten().tolist(),
        mapout,
        sizeMap[0],
        sizeMap[1],
        sizeMap[2],
        1)
    janas_core.replaceMrcHeader(mrcFileList[0], mapout)


def normalize_preserving_sign(arr):
    # Create an empty array of the same shape as arr to hold the normalized
    # values
    normalized = np.zeros_like(arr, dtype=float)

    # Process negative values (if any)
    neg_indices = arr < 0
    if np.any(neg_indices):
        negative_values = arr[neg_indices]
        normalized[neg_indices] = (negative_values - np.min(negative_values)) / (
            0 - np.min(negative_values)
        ) - 1

    # Process positive values (if any)
    pos_indices = arr > 0
    if np.any(pos_indices):
        positive_values = arr[pos_indices]
        normalized[pos_indices] = (
            positive_values - 0) / (np.max(positive_values) - 0)

    return normalized

#############################
# map_histogram_utils
#############################
#############################
# map_histogram_stats (no plotting)
#############################
def map_histogram_stats(map_path, mask_path=None, bins=256, clip=None):
    """
    Compute histogram, summary statistics, and local extrema for an MRC volume
    (optionally masked). Returns a dict with everything needed for plotting.

    Parameters
    ----------
    map_path : str
        Path to input .mrc volume.
    mask_path : Optional[str]
        Path to a mask .mrc with same shape; voxels >0.1 are selected. If None, use all voxels.
    bins : int
        Number of histogram bins.
    clip : Optional[Tuple[float, float]]
        Optional percentile range (low, high) in [0, 100] to clip data before histogramming.

    Returns
    -------
    out : dict
        {
          "region_desc": str,
          "data": np.ndarray,           # selected voxel values (after mask/clip), for boxplot/rug
          "N": int,                     # number of voxels used
          "edges": np.ndarray,          # bin edges
          "centers": np.ndarray,        # bin centers
          "counts": np.ndarray,         # counts per bin (integers)
          "y_plot": np.ndarray,         # percentage per bin (counts/N*100)
          "dmin": float, "dmax": float, # plotting x-limits
          "peak_idx": List[int],        # local maxima bin indices (non-zero neighbours)
          "trough_idx": List[int],      # local minima bin indices (non-zero neighbours)
          "mean": float, "std": float,  # summary
          "p1": float, "p25": float, "p50": float, "p75": float, "p99": float
        }
    """
    import numpy as np

    # Load map
    size_map = janas_core.sizeMRC(map_path)  # (nz, ny, nx)
    vol = np.array(janas_core.ReadMRC(map_path)).reshape(size_map).astype(np.float32, copy=False)

    # Optional mask
    if mask_path is not None:
        size_mask = janas_core.sizeMRC(mask_path)
        if tuple(size_mask) != tuple(size_map):
            raise ValueError(f"Mask size {size_mask} does not match map size {size_map}")
        m = np.array(janas_core.ReadMRC(mask_path)).reshape(size_map)
        sel = m > 0.1
        data = vol[sel]
        region_desc = f"Masked voxels (>0.1) [{np.count_nonzero(sel)} / {vol.size}]"
    else:
        data = vol.ravel()
        region_desc = f"All voxels [{vol.size}]"

    # Drop NaNs / infs
    data = data[np.isfinite(data)]
    if data.size == 0:
        raise ValueError("No finite voxels to analyse after masking / filtering.")

    # Optional percentile clipping
    if clip is not None:
        low, high = float(clip[0]), float(clip[1])
        if not (0.0 <= low < high <= 100.0):
            raise ValueError("clip must be two percent values in [0,100] with low < high")
        lo = np.percentile(data, low)
        hi = np.percentile(data, high)
        data = data[(data >= lo) & (data <= hi)]

    # Histogram (true counts)
    nbins = int(bins)
    dmin, dmax = float(np.min(data)), float(np.max(data))
    if dmin == dmax:
        dmin, dmax = dmin - 0.5, dmax + 0.5

    edges = np.linspace(dmin, dmax, nbins + 1, dtype=np.float64)
    counts, edges = np.histogram(data, bins=edges, density=False)
    centers = 0.5 * (edges[:-1] + edges[1:])
    N = data.size

    # Percentages for plotting
    y_plot = (counts.astype(np.float64) / N) * 100.0

    # Local maxima & minima among non-zero neighbours
    nz_idx = np.flatnonzero(counts > 0)
    peak_idx, trough_idx = [], []
    for k in range(1, len(nz_idx) - 1):
        i_prev, i, i_next = nz_idx[k - 1], nz_idx[k], nz_idx[k + 1]
        if counts[i] > counts[i_prev] and counts[i] > counts[i_next]:
            peak_idx.append(i)
        if counts[i] < counts[i_prev] and counts[i] < counts[i_next]:
            trough_idx.append(i)

    # Summary stats
    mean = float(np.mean(data))
    std = float(np.std(data, ddof=0))
    p1, p25, p50, p75, p99 = np.percentile(data, [1, 25, 50, 75, 99])

    # Console report (counts reported as integers)
    print(f"[map_histogram] {region_desc}")
    print(f"[map_histogram] N={N}  mean={mean:.6g}  std={std:.6g}  "
          f"p1={p1:.6g}  p25={p25:.6g}  median={p50:.6g}  p75={p75:.6g}  p99={p99:.6g}")
    if peak_idx:
        print("[map_histogram] Local maxima (counts):")
        for idx in peak_idx:
            print(f"  - intensity={centers[idx]:.6g}, count={int(counts[idx])}")
    if trough_idx:
        print("[map_histogram] Local minima (counts):")
        for idx in trough_idx:
            print(f"  - intensity={centers[idx]:.6g}, count={int(counts[idx])}")

    return {
        "region_desc": region_desc,
        "data": data,
        "N": N,
        "edges": edges,
        "centers": centers,
        "counts": counts,
        "y_plot": y_plot,
        "dmin": dmin, "dmax": dmax,
        "peak_idx": peak_idx, "trough_idx": trough_idx,
        "mean": mean, "std": std,
        "p1": p1, "p25": p25, "p50": p50, "p75": p75, "p99": p99
    }


#############################
# map_histogram_plot (plotting only)
#############################
def map_histogram_plot(args):
    import os.path
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    stats = map_histogram_stats(
        map_path=args.map,
        mask_path=args.mask,
        bins=int(args.bins),
        clip=(tuple(args.clip) if args.clip is not None else None)
    )

    centers   = stats["centers"]
    counts    = stats["counts"]
    y_plot    = stats["y_plot"]
    data      = stats["data"]
    N         = stats["N"]
    dmin      = stats["dmin"]
    dmax      = stats["dmax"]
    peak_idx  = stats["peak_idx"]
    trough_idx= stats["trough_idx"]
    region_desc = stats["region_desc"]

    # Figure
    fig = plt.figure(figsize=(10, 10.0))
    gs = fig.add_gridspec(nrows=4, ncols=1,
                          height_ratios=[3.0, 1.2, 1.6, 1.6],
                          hspace=0.55)

    # Panel 1: histogram + dotted line
    ax = fig.add_subplot(gs[0, 0])
    markerline, stemlines, baseline = ax.stem(centers, y_plot)
    markerline.set_visible(False)
    baseline.set_visible(False)
    try:
        stemlines.set_linewidth(0.8); stemlines.set_alpha(0.85)
    except Exception:
        for st in stemlines:
            st.set_linewidth(0.8); st.set_alpha(0.85)
    nonzero = counts > 0
    ax.plot(centers[nonzero], y_plot[nonzero], linestyle=":", linewidth=1.8, color="C1")
    ax.set_xlim(dmin, dmax)
    ax.set_ylim(bottom=0.0)
    ax.set_ylabel("% of voxels")
    ax.set_title("Histogram bins (stems) with dotted connection")
    ax.grid(True, alpha=0.3)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=100.0, decimals=0))

    # Panel 2: boxplot + rug
    bx = fig.add_subplot(gs[1, 0], sharex=ax)
    bx.boxplot(data, vert=False, whis=(5, 95), showfliers=False)
    bx.set_title("Box-and-whisker (5–95%)", pad=6)
    bx.set_yticks([])
    bx.grid(True, axis="x", alpha=0.3)
    max_rug = 2000
    if N > max_rug:
        rng = np.random.default_rng(0)
        rug_vals = rng.choice(data, size=max_rug, replace=False)
    else:
        rug_vals = data
    bx.plot(rug_vals, np.zeros_like(rug_vals), "|", alpha=0.25)

    # Compact, non-scientific intensity formatter for labels in Panels 3 & 4
    def _fmt_compact_intensity(x):
        axv = abs(x)
        if axv >= 1000:
            return str(int(round(x)))
        if axv >= 100:
            return f"{x:.0f}"
        if axv >= 10:
            s = f"{x:.1f}"
        else:
            s = f"{x:.2f}"
        return s.rstrip("0").rstrip(".")

    # Panel 3: local maxima
    cx = fig.add_subplot(gs[2, 0], sharex=ax)
    cx.set_title("Local maxima")
    cx.set_ylabel("% of voxels")
    cx.grid(True, axis="y", alpha=0.25)
    cx.set_ylim(0.0, max(y_plot) * 1.2 if peak_idx else 1.0)
    cx.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=100.0, decimals=0))

    # draw peaks + horizontal intensity labels (no tick changes)
    for i in peak_idx:
        x0 = centers[i]
        y0 = y_plot[i]
        cx.axvline(x0, linestyle="--", linewidth=1.0)
        cx.plot([x0], [y0], marker="o", markersize=4)
        ymin, ymax = cx.get_ylim()
        dy = 0.04 * (ymax - ymin)
        y_txt = min(y0 + dy, ymax * 0.98)
        cx.text(x0, y_txt, _fmt_compact_intensity(x0),
                rotation=0, va="bottom", ha="center",
                fontsize=8, clip_on=True)

    # Panel 4: local minima
    dx = fig.add_subplot(gs[3, 0], sharex=ax)
    dx.set_title("Local minima")
    dx.set_xlabel("Intensity")
    dx.set_ylabel("% of voxels")
    dx.grid(True, axis="y", alpha=0.25)
    dx.set_ylim(0.0, max(y_plot) * 1.2 if trough_idx else 1.0)
    dx.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=100.0, decimals=0))

    # draw minima + horizontal intensity labels (no tick changes)
    for i in trough_idx:
        x0 = centers[i]
        y0 = y_plot[i]
        dx.axvline(x0, linestyle="--", linewidth=1.0)
        dx.plot([x0], [y0], marker="o", markersize=4, color="C2")
        ymin, ymax = dx.get_ylim()
        dy = 0.04 * (ymax - ymin)
        y_txt = min(y0 + dy, ymax * 0.98)
        dx.text(x0, y_txt, _fmt_compact_intensity(x0),
                rotation=0, va="bottom", ha="center",
                fontsize=8, color="C2", clip_on=True)

    fig.suptitle(f"{os.path.basename(args.map)}  |  {region_desc}", y=0.995, fontsize=10)
    plt.tight_layout()
    plt.show()


#############################
# Backward-compatible alias
#############################
def map_histogram_utils(args):
    return map_histogram_plot(args)




#############################
# create CTF corrected stack
def create_CTF_reprojectionCorrected_stack(
    particlesStarFile,
    halfMap1,
    halfMap2,
    maskIn,
    outputStackBasename,
    saveOriginal=False,
):
    imageNameTag = "_rlnImageName"
    suffix_particles = ""
    doCTF = True

    # print("project_volume")
    sizeMap = janas_core.sizeMRC(halfMap1)
    # print("sizeMap=",sizeMap)

    inputMapH1_original = janas_core.ReadMRC(halfMap1)
    inputMapH1 = np.array(inputMapH1_original)

    inputMapH2_original = janas_core.ReadMRC(halfMap2)
    inputMapH1 = np.array(inputMapH2_original)

    inputMask_original = janas_core.ReadMRC(maskIn)

    version = starHandler.infoStarFile(particlesStarFile)[2]
    image_names = starHandler.read_star_columns_from_sections(
        particlesStarFile, suffix_particles, "_rlnImageName"
    )
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
        angpix = coordinates["_rlnImagePixelSize"]
        coordinates = coordinates.drop(["_rlnOpticsGroup"], axis=1).reindex()
        coordinates = coordinates.set_index("idx")
        coordinates = coordinates.drop(
            ["_rlnOriginXAngst", "_rlnOriginYAngst"], axis=1)
    else:
        angpix = ctfParameters.at[0, "_rlnDetectorPixelSize"]
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
        if "_rlnPhaseShift" not in columns:
            PhaseShift = pd.DataFrame(np.zeros(len(coordinates)))
        else:
            PhaseShift = starHandler.readColumns(
                particlesStarFile, ["_rlnPhaseShift"])

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
            ctfParameters = ctfParameters.drop(
                ["_rlnOpticsGroup"], axis=1).reindex()
            ctfParameters = ctfParameters.set_index("idx")
            ctfParameters.rename(
                columns={
                    "_rlnImagePixelSize": "_rlnDetectorPixelSize"},
                inplace=True)
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
    print("num particles=", numParticles)

    imageNames = starHandler.readColumns(particlesStarFile, [imageNameTag])
    tmpLine = imageNames[imageNameTag][0]
    stackName = tmpLine[tmpLine.find("@") + 1:]
    sizeI = janas_core.sizeMRC(stackName)
    janas_core.WriteEmptyMRC(
        outputStackBasename + ".mrcs", sizeI[0], sizeI[1], len(imageNames[imageNameTag])
    )
    nx = sizeI[0]
    ny = sizeI[1]

    for ii in range(0, numParticles):
        # print('iteration ',ii,' out of ', numParticles)
        print("iteration ", ii, " out of ", numParticles, end="\r")
        tmpLine = imageNames[imageNameTag][ii]
        atPosition = tmpLine.find("@")
        imageNo = int(tmpLine[:atPosition])
        stackName = tmpLine[atPosition + 1:]
        angpix = ctfParameters.at[ii, "_rlnDetectorPixelSize"]
        SphericalAberration = ctfParameters.at[ii, "_rlnSphericalAberration"]
        voltage = ctfParameters.at[ii, "_rlnVoltage"]
        DefocusAngle = ctfParameters.at[ii, "_rlnDefocusAngle"]
        DefocusU = ctfParameters.at[ii, "_rlnDefocusU"]
        DefocusV = ctfParameters.at[ii, "_rlnDefocusV"]
        AmplitudeContrast = ctfParameters.at[ii, "_rlnAmplitudeContrast"]
        Bfac = ctfParameters.at[ii, "_rlnCtfBfactor"]
        phase_shift = ctfParameters.at[ii, "_rlnPhaseShift"]
        # print("nx, ny=",nx,"  , ",ny)
        ctfI = janas_core.CtfCenteredImage(
            nx,
            ny,
            angpix,
            SphericalAberration,
            voltage,
            DefocusAngle,
            DefocusU,
            DefocusV,
            AmplitudeContrast,
            Bfac,
            phase_shift,
        )
        ctfI = np.array(ctfI).reshape((sizeI[1], sizeI[0]))

        # janas_core.WriteMRC(ctfI.real.flatten().tolist(), 'ctf.mrc' ,sizeMap[1],sizeMap[0],1,angpix)

        # OK block for CTF correction
        mapI = np.array(janas_core.ReadMrcSlice(stackName, imageNo - 1))
        mapI = mapI.reshape((sizeI[1], sizeI[0]))
        I_fft = np.fft.fftn(mapI)
        I_abs_fft = np.abs(I_fft) + 0.0000001

        phi = coordinates.at[ii, "_rlnAngleRot"]
        theta = coordinates.at[ii, "_rlnAngleTilt"]
        psi = coordinates.at[ii, "_rlnAnglePsi"]
        tx = coordinates.at[ii, "_rlnOriginX"]
        ty = coordinates.at[ii, "_rlnOriginY"]
        hm = coordinates.at[ii, "_rlnRandomSubset"]
        if hm == 1:
            RI = janas_core.projectMap(
                inputMapH1_original,
                sizeMap[0],
                sizeMap[1],
                sizeMap[2],
                phi,
                theta,
                psi,
                tx,
                ty,
                0,
            )
        else:
            RI = janas_core.projectMap(
                inputMapH2_original,
                sizeMap[0],
                sizeMap[1],
                sizeMap[2],
                phi,
                theta,
                psi,
                tx,
                ty,
                0,
            )
        RM = janas_core.projectMap(
            inputMask_original,
            sizeMap[0],
            sizeMap[1],
            sizeMap[2],
            phi,
            theta,
            psi,
            tx,
            ty,
            0,
        )
        RI = np.reshape(RI, [sizeMap[1], sizeMap[0]])
        RM = np.reshape(RM, [sizeMap[1], sizeMap[0]])
        minVal = RM.min()
        maxVal = RM.max()
        RM = (RM - minVal) / (maxVal - minVal)
        RI_subtracted = RM * RI
        I_subtracted = RM * mapI

        I_full_fft = np.fft.fftn(mapI)
        I_subtracted_fft = np.fft.fftn(I_subtracted)
        RI_subtracted_fft = np.fft.fftn(RI_subtracted)
        diffI = I_subtracted_fft - (
            RI_subtracted_fft
            * np.abs(I_subtracted)
            / (np.abs(RI_subtracted) + 0.0000001)
        )

        diffI = diffI * np.abs(I_full_fft) / (np.abs(diffI) + 0.0000001)
        # diffI_real_space=np.fft.ifftn(diffI)

        imageToSave = ctfI * np.abs(np.fft.fftshift(diffI))
        imageToSave = normalize_preserving_sign(imageToSave)
        imageToSave = np.fft.ifftn(
            np.fft.ifftshift(imageToSave * np.fft.fftshift(I_full_fft))
        )

        janas_core.ReplaceMrcSlice(
            imageToSave.real.flatten().tolist(),
            outputStackBasename + ".mrcs",
            sizeI[0],
            sizeI[1],
            ii,
        )
        # janas_core.ReplaceMrcSlice(mapI.tolist(),outputStackBasename+'.mrcs',sizeI[0],sizeI[1],ii)

    originalStarDataframe = starHandler.readStar(particlesStarFile)
    if saveOriginal:
        print("keeping the original")
        originalStarDataframe["_original_rlnImageName"] = imageNames

    outImageNames = []
    for ii in range(0, numParticles):
        outImageNames.append(
            str(str(ii + 1).zfill(7) + "@" + outputStackBasename + ".mrcs")
        )

    originalStarDataframe["_rlnImageName"] = outImageNames
    # print (originalStarDataframe)
    starHandler.writeDataframeToStar(
        particlesStarFile, outputStackBasename + ".star", originalStarDataframe
    )


# ===== Shared helpers for stack creation =====
def infer_cs_project_root_from_path(cs_path: str) -> str:
    import re
    d = path.dirname(path.abspath(cs_path))
    while True:
        base = path.basename(d)
        if not re.fullmatch(r"J\d+", base):
            return d
        d = path.dirname(d)

def resolve_stack_path_from_image_name(image_name: str, project_root: Optional[str]) -> Tuple[int, str]:
    at = image_name.find("@")
    if at < 0:
        raise ValueError(f"Invalid _rlnImageName entry (no '@'): {image_name}")
    idx_str = image_name[:at].strip()
    stack_rel = image_name[at + 1:].strip()
    try:
        image_no = int(idx_str)
    except Exception:
        raise ValueError(f"Invalid image index in _rlnImageName: {image_name}")
    if path.isabs(stack_rel):
        stack_path = stack_rel
    else:
        stack_path = path.normpath(path.join(project_root, stack_rel)) if project_root else path.normpath(stack_rel)
    return image_no - 1, stack_path


# Supported path-resolution modes for create_stack_from_star. See the
# `--path_mode` CLI flag for the user-facing description.
_VALID_PATH_MODES = ("auto", "root", "as_is", "star_dir")


def _candidate_stack_paths(
    stack_rel: str,
    project_root: Optional[str],
    star_dir: Optional[str],
    path_mode: str,
) -> List[str]:
    """
    Return the candidate absolute paths to try for ``stack_rel`` under the
    requested ``path_mode``, in priority order.

    Absolute ``stack_rel`` values are always honoured verbatim regardless of
    mode. Relative values are joined to one or more bases depending on the
    mode:

    - ``"root"``    — join to ``project_root`` (or fall back to ``stack_rel``
      as-is if no root was given).
    - ``"as_is"``   — use ``stack_rel`` literally (CWD-relative when not
      absolute).
    - ``"star_dir"`` — join to ``star_dir`` (or fall back to ``stack_rel``).
    - ``"auto"``    — try ``project_root``, then the path as-written, then
      ``star_dir``. The caller picks the first existing one.

    Strict modes return a single candidate; ``"auto"`` returns up to three.
    """
    if path_mode not in _VALID_PATH_MODES:
        raise ValueError(
            f"Unknown path_mode {path_mode!r}. "
            f"Expected one of: {', '.join(_VALID_PATH_MODES)}."
        )

    # Absolute paths are taken at face value in every mode.
    if path.isabs(stack_rel):
        return [stack_rel]

    if path_mode == "root":
        if project_root:
            return [path.normpath(path.join(project_root, stack_rel))]
        return [path.normpath(stack_rel)]

    if path_mode == "as_is":
        return [path.normpath(stack_rel)]

    if path_mode == "star_dir":
        if star_dir:
            return [path.normpath(path.join(star_dir, stack_rel))]
        return [path.normpath(stack_rel)]

    # path_mode == "auto"
    candidates: List[str] = []
    if project_root:
        candidates.append(path.normpath(path.join(project_root, stack_rel)))
    candidates.append(path.normpath(stack_rel))
    if star_dir:
        candidates.append(path.normpath(path.join(star_dir, stack_rel)))
    # de-duplicate while preserving order
    seen = set()
    uniq: List[str] = []
    for c in candidates:
        if c not in seen:
            seen.add(c)
            uniq.append(c)
    return uniq


def _resolve_existing_stack_path(
    image_name: str,
    project_root: Optional[str],
    star_dir: Optional[str],
    path_mode: str,
) -> Tuple[int, str]:
    """
    Decode an ``_rlnImageName`` entry and return ``(slice_idx0, stack_path)``.

    In ``"auto"`` mode this probes the filesystem and returns the first
    candidate that exists; if none exists, the first candidate is returned
    so that the caller can produce a clear error message listing every
    attempted location.

    In strict modes (``"root"``, ``"as_is"``, ``"star_dir"``) the single
    mode-dictated path is returned without filesystem checks; the caller
    decides whether to error on a missing file.
    """
    at = image_name.find("@")
    if at < 0:
        raise ValueError(f"Invalid _rlnImageName entry (no '@'): {image_name}")
    idx_str = image_name[:at].strip()
    stack_rel = image_name[at + 1:].strip()
    try:
        image_no = int(idx_str)
    except Exception:
        raise ValueError(f"Invalid image index in _rlnImageName: {image_name}")

    candidates = _candidate_stack_paths(stack_rel, project_root, star_dir, path_mode)
    if path_mode == "auto":
        for c in candidates:
            if path.exists(c):
                return image_no - 1, c
        # Nothing exists; return the first candidate so the caller can emit
        # an informative error listing all attempted locations.
        return image_no - 1, candidates[0]
    return image_no - 1, candidates[0]



def _auto_particles_section_name(star_path: str) -> str:
    version = starHandler.infoStarFile(star_path)[2]
    return "" if version == "relion_v30" else "particles"





def _open_src_mm(src_stack: str, hdr: Dict[str, Any]) -> np.memmap:
    """Memmap the image payload as (nz, ny, nx)."""
    nx, ny, nz = int(hdr["nx"]), int(hdr["ny"]), int(hdr["nz"])
    dt = _dtype_from_mode(int(hdr["mode"]))
    off = 1024 + int(hdr["nsymbt"])
    return np.memmap(src_stack, mode="r", dtype=dt, offset=off, shape=(nz, ny, nx), order="C")

def _probe_first_stack(
    image_names: List[str],
    project_root: Optional[str],
    star_dir: Optional[str] = None,
    path_mode: str = "auto",
) -> Tuple[int, int, float]:
    idx0, first_stack = _resolve_existing_stack_path(
        image_names[0], project_root, star_dir, path_mode
    )
    if not path.exists(first_stack):
        attempted = _candidate_stack_paths(
            image_names[0].split("@", 1)[1].strip()
            if "@" in image_names[0] else image_names[0],
            project_root, star_dir, path_mode,
        )
        raise FileNotFoundError(
            f"Cannot locate source stack for particle '{image_names[0]}'. "
            f"Path mode: {path_mode}. Tried: {attempted}."
        )
    with open(first_stack, "rb") as f0:
        hdr0 = _read_mrc_header(f0)
    nx, ny = int(hdr0["nx"]), int(hdr0["ny"])
    apix   = float(hdr0["pixel_x"])
    return nx, ny, apix

def create_stack_from_star(
    star_in: str,
    out_root: str,
    project_root: Optional[str] = None,
    provenance_tag: Optional[str] = None,
    zfill_width: int = 6,
    override_section_name: Optional[str] = None,
    chunk_size: int = 128,  # number of slices to move per batch
    path_mode: str = "auto",
) -> Tuple[str, str]:
    """
    Build a consolidated .mrcs with minimal open files:
      - group by source stack
      - open one memmap per group, copy slices (in chunks), close it
      - keep a single output memmap open

    Path resolution for the source stacks referenced by ``_rlnImageName`` is
    controlled by ``path_mode``:

    - ``"auto"`` (default): try ``project_root``, then the path as written
      (CWD-relative), then the directory of ``star_in``. The first existing
      candidate is used.
    - ``"root"``: resolve relative paths against ``project_root`` only
      (matches the historical CryoSPARC ``--root`` behaviour).
    - ``"as_is"``: use the path as written, never prepending ``project_root``.
      Relative paths are taken to be CWD-relative.
    - ``"star_dir"``: resolve relative paths against the directory containing
      ``star_in``.

    The generated output STAR always points to the new consolidated stack;
    ``path_mode`` controls only how the *input* stacks are located on disk.
    """
    image_tag = "_rlnImageName"
    names_df = starHandler.readColumns(star_in, [image_tag])
    image_names = list(names_df[image_tag])
    if not image_names:
        raise ValueError("No _rlnImageName entries found in STAR.")

    star_dir = path.dirname(path.abspath(star_in)) or None

    nx, ny, apix = _probe_first_stack(image_names, project_root, star_dir, path_mode)
    N = len(image_names)
    out_stack = f"{out_root}.mrcs"
    out_star  = f"{out_root}.star"

    # Pre-create output file (mode=2 float32) and map as (N, ny, nx)
    with open(out_stack, "wb") as f:
        _write_mrc_header(f, nx=nx, ny=ny, nz=N, apix=apix,
                          vmin=0.0, vmax=0.0, vmean=0.0, vrms=0.0, mode=2)
        f.seek(1024 + (N * ny * nx * 4) - 1)
        f.write(b"\0")
    out_mm = np.memmap(out_stack, mode="r+", dtype=np.float32,
                       offset=1024, shape=(N, ny, nx), order="C")

    # Build work lists grouped by source stack → [(out_idx, slice_idx0), ...]
    by_stack: Dict[str, List[Tuple[int, int]]] = {}
    for ii, img_name in enumerate(image_names):
        slice_idx0, src_stack = _resolve_existing_stack_path(
            img_name, project_root, star_dir, path_mode
        )
        if not path.exists(src_stack):
            stack_rel = img_name.split("@", 1)[1].strip() if "@" in img_name else img_name
            attempted = _candidate_stack_paths(
                stack_rel, project_root, star_dir, path_mode
            )
            raise FileNotFoundError(
                f"Missing source stack for particle '{img_name}'. "
                f"Path mode: {path_mode}. Tried: {attempted}."
            )
        by_stack.setdefault(src_stack, []).append((ii, slice_idx0))

    # Running stats (vectorised)
    total_count = 0
    sum_vals    = 0.0
    sum_squares = 0.0
    vmin        = np.inf
    vmax        = -np.inf

    # Process each source stack independently: open → copy → close
    processed = 0
    for src_stack, pairs in by_stack.items():
        # Header + geometry/type validation
        with open(src_stack, "rb") as fs:
            hdr = _read_mrc_header(fs)
        if int(hdr["nx"]) != nx or int(hdr["ny"]) != ny:
            raise ValueError(
                f"Geometry mismatch: first stack is {nx}x{ny}, but {src_stack} is {hdr['nx']}x{hdr['ny']}"
            )

        mm = _open_src_mm(src_stack, hdr)  # open only for this group

        # Sort by slice index for better locality (optional but cheap)
        pairs.sort(key=lambda p: p[1])
        out_idx_all = np.fromiter((p[0] for p in pairs), dtype=np.int64)
        slice_idx_all = np.fromiter((p[1] for p in pairs), dtype=np.int64)

        # Copy in chunks to cap peak RAM
        for start in range(0, len(pairs), chunk_size):
            stop = min(start + chunk_size, len(pairs))
            out_idx_chunk   = out_idx_all[start:stop]
            slice_idx_chunk = slice_idx_all[start:stop]

            # Read a (k, ny, nx) view from the source, cast once
            batch = np.asarray(mm[slice_idx_chunk, :, :], dtype=np.float32, order="C")

            # Write directly to the correct output planes (random access)
            out_mm[out_idx_chunk, :, :] = batch

            # Stats for this chunk
            # Use nan-aware reductions just in case
            vmin = min(vmin, float(np.nanmin(batch)))
            vmax = max(vmax, float(np.nanmax(batch)))
            sum_vals    += float(np.nansum(batch))
            sum_squares += float(np.nansum(batch * batch))
            total_count += int(batch.size)

            processed += int(len(out_idx_chunk))
            if (processed % 1000) == 0 or processed == N:
                print(f"Copied {processed}/{N} particles", end="\r")

            # free batch early
            del batch

        # Close this source memmap immediately
        del mm

    # Flush and close output map
    out_mm.flush()
    del out_mm

    # Final stats
    if total_count > 0:
        vmean = sum_vals / total_count
        var = (sum_squares / total_count) - (vmean * vmean)
        var = max(var, 0.0)
        vrms = float(np.sqrt(var))
    else:
        vmean = 0.0
        vrms = 0.0

    # Rewrite header with final stats
    with open(out_stack, "r+b") as f:
        _write_mrc_header(f, nx=nx, ny=ny, nz=N, apix=apix,
                          vmin=float(vmin), vmax=float(vmax),
                          vmean=float(vmean), vrms=float(vrms), mode=2)

    # STAR rewrite
    star_df = starHandler.readStar(star_in)
    if len(star_df) != N:
        raise ValueError(f"STAR row count ({len(star_df)}) does not match image list length ({N}).")

    if provenance_tag:
        star_df[provenance_tag] = image_names
    star_df["_rlnImageName"] = [f"{str(i+1).zfill(zfill_width)}@{out_stack}" for i in range(N)]
    if provenance_tag and provenance_tag in star_df.columns:
        ordered_cols = [c for c in star_df.columns if c != provenance_tag] + [provenance_tag]
        star_df = star_df[ordered_cols]

    section_name = override_section_name if override_section_name is not None else _auto_particles_section_name(star_in)
    starHandler.update_star_columns_from_sections(
        filenameIn=star_in,
        filenameOut=out_star,
        section_name=section_name,
        df=star_df
    )

    print(f"Wrote stack: {out_stack}")
    print(f"Wrote STAR:  {out_star}")
    return out_stack, out_star


# ———————————————————————————————————————————————————————————————
# backmap_stars — inverse companion of create_stack_from_star
# ———————————————————————————————————————————————————————————————


def _read_particle_section_df(star_path: str, section_name: str) -> pd.DataFrame:
    """
    Read the full particle section of ``star_path`` into a pandas DataFrame,
    preserving every column and the on-disk row order.

    This is a thin reader on top of ``starHandler.read_star_sections`` /
    ``process_section`` that does **not** mutate or reformat values: every cell
    is returned as a string. Used internally by :func:`backmap_stars` to read
    the processed STAR file when we need to inspect (and not just rewrite)
    the particle table.
    """
    sections = starHandler.read_star_sections(star_path)
    if section_name not in sections:
        available = ", ".join(repr(k) for k in sections.keys())
        raise ValueError(
            f"STAR file '{star_path}' has no data block named 'data_{section_name}'. "
            f"Available blocks: {available or '(none)'}."
        )

    section_lines = sections[section_name]
    header_lines = [ln for ln in section_lines if ln.startswith("_")]
    data_lines = [
        ln
        for ln in section_lines
        if not ln.startswith("_") and not ln.startswith("loop_")
    ]
    headers = [ln.split()[0] for ln in header_lines]
    rows = [ln.split() for ln in data_lines]
    return pd.DataFrame(rows, columns=headers)


def backmap_stars(
    processed_star: str,
    mapping_star: str,
    output_star: str,
    image_tag: str = "_rlnImageName",
    source_tag: str = "_janas_source_rlnImageName",
    stack_reference_tag: Optional[str] = "_janas_stack_rlnImageName",
    section_name: Optional[str] = None,
    strict: bool = True,
) -> Dict[str, Any]:
    """
    Restore original source particle image names in a downstream STAR file.

    :func:`create_stack_from_star` can rewrite ``_rlnImageName`` so that particles
    point to a consolidated ``.mrcs`` stack (for example
    ``000022@EMPIAR_12707_stack.mrcs``). When that helper is called with a
    ``provenance_tag``, the original source reference (for example
    ``009640@J1149/restack/batch_6_restacked.mrc``) is preserved in a parallel
    column, usually ``_janas_source_rlnImageName``.

    Downstream processing (RELION refinement, cryoSPARC round-trips, manual
    selection scripts, ...) sometimes preserves the rewritten
    ``_rlnImageName`` column but drops the provenance column. The resulting
    STAR file then no longer carries the link back to the original micrograph-
    level particles.

    ``backmap_stars`` restores that link by joining the downstream STAR file
    against the stack-generation STAR file:

        ``processed_star[_rlnImageName]``
            → looked up in ``mapping_star[_rlnImageName]``
            → replaced with ``mapping_star[_janas_source_rlnImageName]``

    The processed STAR supplies the particle rows and all refined metadata
    (angles, origins, defocus, optics group, random subset, class assignments,
    JANAS scores, ...) that must be preserved exactly. The mapping STAR
    supplies only the relation between consolidated stack names and original
    source names.

    The function uses the stack ``_rlnImageName`` as the lookup key — never
    row number, angles, defocus, origins, or class number — so it is safe to
    apply to a re-ordered subset of the original particles.

    Parameters
    ----------
    processed_star : str
        Downstream STAR file to fix. Its ``_rlnImageName`` currently points
        to the consolidated stack. All other columns and row order are
        preserved in the output.
    mapping_star : str
        STAR file produced by :func:`create_stack_from_star` with
        ``provenance_tag=_janas_source_rlnImageName`` (or whatever is passed
        as ``source_tag``). Must contain both ``image_tag`` and ``source_tag``
        in its particle section.
    output_star : str
        Path of the new STAR file to write. Existing files at this path are
        overwritten.
    image_tag : str, default ``"_rlnImageName"``
        Column that holds the consolidated-stack reference (both in
        ``processed_star`` and ``mapping_star``).
    source_tag : str, default ``"_janas_source_rlnImageName"``
        Column in ``mapping_star`` that holds the original source reference.
    stack_reference_tag : str or None, default ``"_janas_stack_rlnImageName"``
        If not ``None``, an audit column with this name is added (or replaced)
        in the output. Each row carries the previous (stack-based) value of
        ``image_tag`` before the replacement, making the operation reversible.
        Pass ``None`` to skip the audit column.
    section_name : str or None, default ``None``
        Name of the particle data block (without the ``data_`` prefix). If
        ``None``, it is inferred per file via
        :func:`_auto_particles_section_name` (``"particles"`` for
        RELION 3.1, ``""`` for older STAR files).
    strict : bool, default ``True``
        If ``True``, raise :class:`ValueError` when at least one value in
        ``processed_star[image_tag]`` is missing from the mapping. If
        ``False``, leave such values unchanged, write the output, and report
        them in the returned dictionary.

    Returns
    -------
    dict
        Report dictionary with the following keys:

        - ``processed_star``, ``mapping_star``, ``output_star`` — input/output
          file paths echoed back.
        - ``section_name`` — particle section name used for the output write.
        - ``image_tag``, ``source_tag``, ``stack_reference_tag`` — column
          names actually used.
        - ``n_processed`` — number of rows in the processed STAR.
        - ``n_mapping_rows`` — number of rows in the mapping STAR.
        - ``n_mapped`` — number of processed rows successfully remapped.
        - ``n_missing`` — number of processed rows with no entry in the
          mapping table (always 0 when ``strict=True`` succeeds).
        - ``missing_examples`` — up to 10 example missing keys (only
          populated when ``strict=False``).
        - ``duplicate_keys`` — list of consolidated keys that appear more
          than once in the mapping STAR but always map to the same source
          (kept for audit; conflicting duplicates raise instead).

    Raises
    ------
    ValueError
        - ``mapping_star`` does not contain ``image_tag`` or ``source_tag``.
        - either STAR file has no usable particle section.
        - either STAR file has zero rows.
        - the mapping table has duplicated keys that resolve to *different*
          source names.
        - some ``processed_star[image_tag]`` values are missing from the
          mapping when ``strict=True``.

    Examples
    --------
    >>> from janas import utils
    >>> report = utils.backmap_stars(
    ...     processed_star="run_it025_data.star",
    ...     mapping_star="EMPIAR_12707_stack.star",
    ...     output_star="run_it025_data_backmapped.star",
    ... )
    >>> print(report["n_mapped"], "/", report["n_processed"])
    """
    # ---------- 1. Resolve section names per file ----------
    mapping_section = (
        section_name
        if section_name is not None
        else _auto_particles_section_name(mapping_star)
    )
    processed_section = (
        section_name
        if section_name is not None
        else _auto_particles_section_name(processed_star)
    )

    # ---------- 2. Read mapping STAR (only the two columns we need) ----------
    map_df = starHandler.read_star_columns_from_sections(
        mapping_star, mapping_section, [image_tag, source_tag]
    )
    if map_df is None or len(map_df) == 0:
        raise ValueError(
            f"Mapping STAR '{mapping_star}' has no particle rows in section "
            f"'data_{mapping_section}' (or the section was not found)."
        )
    for required in (image_tag, source_tag):
        if required not in map_df.columns:
            raise ValueError(
                f"Mapping STAR '{mapping_star}' is missing column '{required}' "
                f"in section 'data_{mapping_section}'. "
                f"Found columns: {list(map_df.columns)}."
            )

    n_mapping_rows = int(len(map_df))

    # ---------- 3. Build the lookup, detecting duplicates ----------
    lookup: Dict[str, str] = {}
    duplicate_keys: List[str] = []
    conflicting: Dict[str, set] = {}
    for key, value in zip(map_df[image_tag], map_df[source_tag]):
        if key in lookup:
            if lookup[key] != value:
                conflicting.setdefault(key, set()).add(lookup[key])
                conflicting[key].add(value)
            elif key not in duplicate_keys:
                duplicate_keys.append(key)
        else:
            lookup[key] = value

    if conflicting:
        examples = list(conflicting.items())[:5]
        bullets = "\n  - ".join(
            f"{k!r} → {sorted(v)}" for k, v in examples
        )
        raise ValueError(
            f"Mapping STAR '{mapping_star}' has {len(conflicting)} stack "
            f"key(s) that resolve to more than one source name. "
            f"Examples:\n  - {bullets}"
        )

    # ---------- 4. Read the processed STAR particle section in full ----------
    processed_df = _read_particle_section_df(processed_star, processed_section)
    if len(processed_df) == 0:
        raise ValueError(
            f"Processed STAR '{processed_star}' has no particle rows in "
            f"section 'data_{processed_section}'."
        )
    if image_tag not in processed_df.columns:
        raise ValueError(
            f"Processed STAR '{processed_star}' is missing column "
            f"'{image_tag}' in section 'data_{processed_section}'. "
            f"Found columns: {list(processed_df.columns)}."
        )

    n_processed = int(len(processed_df))

    # ---------- 5. Apply the mapping ----------
    original_image_values = processed_df[image_tag].tolist()
    missing_keys: List[str] = []
    new_image_values: List[str] = []
    for value in original_image_values:
        if value in lookup:
            new_image_values.append(lookup[value])
        else:
            missing_keys.append(value)
            new_image_values.append(value)  # leave unchanged (used only if not strict)

    n_missing = len(missing_keys)
    n_mapped = n_processed - n_missing

    if n_missing > 0 and strict:
        sample = missing_keys[:10]
        bullets = "\n  - ".join(repr(k) for k in sample)
        raise ValueError(
            f"{n_missing}/{n_processed} particles in '{processed_star}' have "
            f"an '{image_tag}' value that is absent from the mapping STAR "
            f"'{mapping_star}' (column '{image_tag}'). "
            f"Example missing keys:\n  - {bullets}\n"
            f"Pass strict=False to leave unmapped rows unchanged and report "
            f"them in the result dictionary."
        )

    # ---------- 6. Prepare the update DataFrame ----------
    # process_section overwrites/adds the columns we pass via positional
    # alignment, so we hand it only the columns we actually want to change.
    update_df = pd.DataFrame({image_tag: new_image_values})
    if stack_reference_tag:
        update_df[stack_reference_tag] = original_image_values

    # ---------- 7. Write the output STAR ----------
    starHandler.update_star_columns_from_sections(
        filenameIn=processed_star,
        filenameOut=output_star,
        section_name=processed_section,
        df=update_df,
    )

    return {
        "processed_star": processed_star,
        "mapping_star": mapping_star,
        "output_star": output_star,
        "section_name": processed_section,
        "image_tag": image_tag,
        "source_tag": source_tag,
        "stack_reference_tag": stack_reference_tag,
        "n_processed": n_processed,
        "n_mapping_rows": n_mapping_rows,
        "n_mapped": n_mapped,
        "n_missing": n_missing,
        "missing_examples": missing_keys[:10],
        "duplicate_keys": duplicate_keys,
    }


# ———————————————————————————————————————————————————————————————
# ———————————————————————————————————————————————————————————————
# Helper for rewriting CryoSPARC blob paths in _rlnImageName
# ———————————————————————————————————————————————————————————————

# Matches a numeric CryoSPARC prefix at the start of a filename: digits
# followed by a single underscore. Only used against the filename component
# of a blob path, never against the directory part.
_CSPARC_NUMERIC_PREFIX_RE = re.compile(r"^\d+_")

# Matches a terminal "_particles" immediately before the file extension, e.g.
# "..._particles.mrc" -> "....mrc" (stem). Multi-component extensions are
# preserved because we operate on os.path.splitext (which only splits the
# last dot, as RELION/cryoSPARC files use ".mrc" / ".mrcs").
_CSPARC_PARTICLES_SUFFIX_RE = re.compile(r"_particles$")


def _clean_csparc_blob_path(
    blob_path: str,
    clean_path: bool = False,
    clean_prefix: bool = False,
    clean_suffix: bool = False,
    fix_path: Optional[str] = None,
) -> str:
    """
    Rewrite a CryoSPARC ``blob/path`` value before it is used as a RELION
    ``_rlnImageName``.

    The operations are applied in a fixed order and only affect the
    filename / directory components of ``blob_path``:

    1. backslashes are converted to forward slashes;
    2. leading ``>`` characters (CryoSPARC stream marker) are stripped;
    3. the path is split into directory and filename on the last ``/``;
    4. if ``clean_prefix``, the leading ``<digits>_`` prefix is stripped from
       the *filename*;
    5. if ``clean_suffix``, a terminal ``_particles`` is stripped from the
       filename *stem* (the extension is preserved);
    6. directory selection:

       - if ``fix_path`` is given, the original directory is replaced by it
         (it takes precedence over ``clean_path``); a trailing ``/`` in
         ``fix_path`` is normalised so the output does not contain ``//``;
       - else if ``clean_path``, the directory is dropped entirely and only
         the filename is returned;
       - otherwise the original directory is preserved.

    When no option is supplied the path is returned essentially unchanged
    (other than the backslash / leading-``>`` cleanups, which mirror what
    :func:`csparc2star` already does on ``raw_paths``).
    """
    path_str = str(blob_path).replace("\\", "/").lstrip(">")

    # Split on the last '/' to separate directory and filename.
    if "/" in path_str:
        directory, filename = path_str.rsplit("/", 1)
    else:
        directory, filename = "", path_str

    if clean_prefix:
        filename = _CSPARC_NUMERIC_PREFIX_RE.sub("", filename)

    if clean_suffix:
        stem, ext = os.path.splitext(filename)
        new_stem = _CSPARC_PARTICLES_SUFFIX_RE.sub("", stem)
        filename = new_stem + ext

    if fix_path is not None:
        # fix_path takes precedence over clean_path.
        new_dir = str(fix_path).replace("\\", "/").rstrip("/")
        return f"{new_dir}/{filename}" if new_dir else filename
    if clean_path:
        return filename
    if directory:
        return f"{directory}/{filename}"
    return filename


def csparc2star(infile: str,
                outfile: str,
                transform: Optional[str] = None,
                loglevel: str = "WARNING",
                clean_path: bool = False,
                clean_prefix: bool = False,
                clean_suffix: bool = False,
                fix_path: Optional[str] = None,
                missing_pose_to_zero: bool = False) -> None:
    """
    Convert a CryoSPARC .cs to a Relion .star with both data_optics
    and data_particles sections, using alignments3D/pose & shift.

    The optional ``clean_path``, ``clean_prefix``, ``clean_suffix`` and
    ``fix_path`` arguments rewrite the ``blob/path`` portion of the
    generated ``_rlnImageName`` values (see :func:`_clean_csparc_blob_path`
    for the exact semantics). All other STAR columns are unaffected.

    The ``missing_pose_to_zero`` flag controls behaviour when the input
    ``.cs`` file lacks 3D alignment metadata. CryoSPARC extraction,
    picking, passthrough or coordinate-only jobs typically produce files
    without an ``alignments3D/pose`` (and often without
    ``alignments3D/shift``) field. By default this is a hard error, so the
    caller is not silently given an unaligned STAR file. When
    ``missing_pose_to_zero=True``, any missing pose is replaced by zero
    Euler angles and any missing shift by zero origin shifts, producing
    an explicitly unaligned STAR file. The resulting
    ``_rlnOriginXAngst`` / ``_rlnOriginYAngst`` are refinement shifts
    (not particle extraction coordinates) — use this option only when an
    unaligned STAR is intended, or when alignment will be supplied by a
    subsequent processing step.
    """
    log = logging.getLogger("csparc2star")

    # Attach a handler only once (avoid duplicate lines)
    if not log.handlers:
        handler = logging.StreamHandler()
        log.addHandler(handler)

    # Normalise loglevel (accept "info", "INFO", etc.)
    if isinstance(loglevel, str):
        loglevel = loglevel.upper()
    log.setLevel(loglevel)

    if not infile.endswith(".cs"):
        raise ValueError("csparc2star only supports .cs files")

    # 1) Load the CryoSPARC .cs array
    cs = np.load(infile, max_header_size=100000)
    N = cs.shape[0]

    # 2) Build an optics table per micrograph
    raw_paths = pd.Series(cs["blob/path"].astype(str), name="_raw_path")
    # remove stray '>' at start of any path
    raw_paths = raw_paths.str.lstrip(">")
    # collapse to parent micrograph folder
    micro = raw_paths.str.split(pat="/extract", n=1).str[0]
    micro.name = "_micrograph"

    raw_optics = pd.DataFrame({
        "_micrograph":             micro,
        "_rlnVoltage":             cs["ctf/accel_kv"],
        "_rlnSphericalAberration": cs["ctf/cs_mm"],
        "_rlnAmplitudeContrast":   cs["ctf/amp_contrast"],
        "_rlnImageSize":           cs["blob/shape"][:, -1].astype(int),
        "_rlnImagePixelSize":      cs["blob/psize_A"],
        "_rlnImageDimensionality": np.full(N, 2, int),
    })

    # deduplicate to one row per micrograph
    optics = (
        raw_optics
        .drop_duplicates(subset="_micrograph", keep="first")
        .reset_index(drop=True)
    )
    optics["_rlnOpticsGroup"] = np.arange(1, len(optics) + 1)
    micro_to_group = dict(zip(optics["_micrograph"], optics["_rlnOpticsGroup"]))
    optics = optics.drop(columns="_micrograph")

    # reorder to put OpticsGroup fourth
    optics = optics[[
        "_rlnVoltage",
        "_rlnSphericalAberration",
        "_rlnAmplitudeContrast",
        "_rlnOpticsGroup",
        "_rlnImageSize",
        "_rlnImagePixelSize",
        "_rlnImageDimensionality",
    ]]

    # 3) Helpers for pose → rotation matrix → ZYZ Euler
    def _expmap_to_R(v: np.ndarray) -> np.ndarray:
        θ = np.linalg.norm(v)
        if θ < 1e-8:
            return np.eye(3)
        k = v / θ
        K = np.array([[   0, -k[2],  k[1]],
                      [ k[2],    0, -k[0]],
                      [-k[1], k[0],     0]])
        return np.eye(3) + np.sin(θ)*K + (1-np.cos(θ))*(K @ K)

    def _R_to_zyz(Rm: np.ndarray):
        β = np.arccos(np.clip(Rm[2,2], -1, 1))
        α = np.arctan2(Rm[1,2], Rm[0,2])
        γ = np.arctan2(Rm[2,1], -Rm[2,0])
        return np.degrees([α, β, γ])

    # 4) Build rotation matrices from CryoSPARC pose
    #
    # Extraction / picking / passthrough .cs files commonly do not contain
    # alignments3D/pose. Without --missing_pose_to_zero this is a hard error,
    # so the user is not silently given an unaligned STAR.
    cs_fields = set(cs.dtype.names)
    pose_present = "alignments3D/pose" in cs_fields
    if pose_present:
        pose_arr = np.stack(cs["alignments3D/pose"])
        Kdim = pose_arr.shape[1]
        if Kdim == 3:
            Rmats = [_expmap_to_R(v) for v in pose_arr]
        elif Kdim == 4:
            q = pose_arr[:, [1, 2, 3, 0]]  # reorder to x,y,z,w
            Rmats = Rotation.from_quat(q).as_matrix()
        elif Kdim == 9:
            Rmats = pose_arr.reshape(-1, 3, 3)
        else:
            raise ValueError(f"Unrecognized pose length {Kdim}")

        # 5) Convert each rotation to φ, θ, ψ
        eulers = np.vstack([_R_to_zyz(Rm) for Rm in Rmats])
        phi, theta, psi = eulers[:, 0], eulers[:, 1], eulers[:, 2]
    elif missing_pose_to_zero:
        log.warning(
            "alignments3D/pose missing in '%s'; writing zero Euler angles "
            "(_rlnAngleRot=_rlnAngleTilt=_rlnAnglePsi=0) for %d particles "
            "because --missing_pose_to_zero was set.",
            infile, N,
        )
        phi = np.zeros(N, dtype=float)
        theta = np.zeros(N, dtype=float)
        psi = np.zeros(N, dtype=float)
    else:
        raise ValueError(
            f"'{infile}' has no 'alignments3D/pose' field. "
            "This is normal for CryoSPARC extraction, picking, passthrough or "
            "coordinate-only jobs, which carry particle blob/location metadata "
            "but no 3D refinement alignment. "
            "csparc2star normally writes the RELION orientation columns "
            "_rlnAngleRot, _rlnAngleTilt, _rlnAnglePsi, _rlnOriginXAngst and "
            "_rlnOriginYAngst from alignments3D/{pose,shift}. "
            "Rerun with '--missing_pose_to_zero' if you intentionally want a "
            "STAR with zero angles and zero origins (the result is an "
            "unaligned STAR file — only do this when that is intended or "
            "when alignments will be supplied by a later processing step). "
            "_rlnOriginXAngst / _rlnOriginYAngst are refinement shifts, "
            "not particle extraction coordinates."
        )

    # 6) Compute origins in Å from shift × pixel size
    #
    # Extraction-style .cs files typically also lack alignments3D/shift. With
    # --missing_pose_to_zero we substitute zero origins and skip the pixel-
    # size lookup entirely (it is only needed to convert shifts to Å). If
    # alignments3D/psize_A is missing but a shift IS present (uncommon), fall
    # back first to blob/psize_A and then to location/micrograph_psize_A.
    shift_present = "alignments3D/shift" in cs_fields
    if shift_present:
        shifts = np.stack(cs["alignments3D/shift"])
        if "alignments3D/psize_A" in cs_fields:
            psize3 = cs["alignments3D/psize_A"]
        elif "blob/psize_A" in cs_fields:
            log.warning(
                "alignments3D/psize_A missing in '%s'; using blob/psize_A to "
                "convert shifts to Å.",
                infile,
            )
            psize3 = cs["blob/psize_A"]
        elif "location/micrograph_psize_A" in cs_fields:
            log.warning(
                "alignments3D/psize_A missing in '%s'; using "
                "location/micrograph_psize_A to convert shifts to Å.",
                infile,
            )
            psize3 = cs["location/micrograph_psize_A"]
        else:
            raise ValueError(
                f"'{infile}' has alignments3D/shift but no pixel-size field "
                "(tried alignments3D/psize_A, blob/psize_A, "
                "location/micrograph_psize_A). Shifts cannot be converted "
                "to Å for _rlnOriginXAngst / _rlnOriginYAngst."
            )
        orig_x = shifts[:, 0] * psize3
        orig_y = shifts[:, 1] * psize3
    elif missing_pose_to_zero:
        log.warning(
            "alignments3D/shift missing in '%s'; writing zero origins "
            "(_rlnOriginXAngst=_rlnOriginYAngst=0) for %d particles "
            "because --missing_pose_to_zero was set.",
            infile, N,
        )
        orig_x = np.zeros(N, dtype=float)
        orig_y = np.zeros(N, dtype=float)
    else:
        raise ValueError(
            f"'{infile}' has no 'alignments3D/shift' field. "
            "Pass --missing_pose_to_zero to write zero origins "
            "(_rlnOriginXAngst=_rlnOriginYAngst=0)."
        )

    # 7) Construct Relion‐style image names
    cleaned_paths = raw_paths.map(
        lambda p: _clean_csparc_blob_path(
            p,
            clean_path=clean_path,
            clean_prefix=clean_prefix,
            clean_suffix=clean_suffix,
            fix_path=fix_path,
        )
    )
    idx_series = pd.Series(cs["blob/idx"].astype(int) + 1)
    image_names = (
        idx_series.astype(str)
                  .str.zfill(6)
                  .str.cat(cleaned_paths, sep="@")
    ).to_numpy()

    # 8) Build the particle table with optics columns
    particles = pd.DataFrame({
        "_rlnImageName":    image_names,
        "_rlnAngleRot":     phi,
        "_rlnAngleTilt":    theta,
        "_rlnAnglePsi":     psi,
        "_rlnOriginXAngst": orig_x,
        "_rlnOriginYAngst": orig_y,
        "_rlnDefocusU":     cs["ctf/df1_A"],
        "_rlnDefocusV":     cs["ctf/df2_A"],
        "_rlnDefocusAngle": np.degrees(cs["ctf/df_angle_rad"]),
        "_rlnPhaseShift":   np.degrees(cs["ctf/phase_shift_rad"]),
        "_rlnCtfBfactor":   cs["ctf/bfactor"],
        **{col: raw_optics[col].values for col in raw_optics.columns},
        "_rlnRandomSubset": (
            cs["alignments3D/split"].astype(int) + 1
            if "alignments3D/split" in cs.dtype.names
            else np.ones(N, int)
        ),
        "_rlnClassNumber": (
            cs["alignments3D/class"].astype(int) + 1
            if "alignments3D/class" in cs.dtype.names
            else np.ones(N, int)
        ),
    })

    # 9) Assign optics group to each particle
    particles["_rlnOpticsGroup"] = micro.map(micro_to_group).astype(int).values

    # 10) Retain only the 14 fields Relion expects
    particle_fields = [
        "_rlnImageName", "_rlnAngleRot", "_rlnAngleTilt", "_rlnAnglePsi",
        "_rlnOriginXAngst", "_rlnOriginYAngst",
        "_rlnDefocusU", "_rlnDefocusV", "_rlnDefocusAngle",
        "_rlnPhaseShift", "_rlnCtfBfactor",
        "_rlnOpticsGroup", "_rlnRandomSubset", "_rlnClassNumber",
    ]
    particles = particles[particle_fields]

    # 11) Format particle rows and compute widths
    formatted_rows = []
    for _, row in particles.iterrows():
        vals = []
        for c in particle_fields:
            v = row[c]
            if c == "_rlnImageName":
                s = str(v)
            elif isinstance(v, float):
                s = f"{v:.6f}"
            else:
                s = str(int(v))
            vals.append(s)
        formatted_rows.append(vals)

    col_widths = {
        c: max(len(c), max(len(r[i]) for r in formatted_rows))
        for i, c in enumerate(particle_fields)
    }

    # 12) Write out the STAR file
    with open(outfile, "w") as f:
        # data_optics
        f.write("data_optics\n\nloop_\n")
        for i, c in enumerate(optics.columns, start=1):
            f.write(f"{c} #{i}\n")
        for _, row in optics.iterrows():
            vals = []
            for c in optics.columns:
                v = row[c]
                if c in ("_rlnOpticsGroup", "_rlnImageSize", "_rlnImageDimensionality"):
                    vals.append(str(int(v)))
                else:
                    vals.append(f"{float(v):.6f}")
            f.write(" ".join(vals) + "\n")

        # data_particles
        f.write("\n\ndata_particles\n\nloop_\n")
        for idx, c in enumerate(particle_fields, start=1):
            f.write(f"{c} #{idx}\n")
        for row_vals in formatted_rows:
            pieces = []
            for val, c in zip(row_vals, particle_fields):
                if c == "_rlnImageName":
                    pieces.append(val.ljust(col_widths[c]))
                else:
                    pieces.append(val.rjust(col_widths[c]))
            f.write(" ".join(pieces) + "\n")

    log.info(f"Wrote STAR with {len(optics)} optics groups and {len(particles)} particles")


# ———————————————————————————————————————————————————————————————
def update_star_from_csparc(csfile: str,
                            starfile_in: str,
                            starfile_out: str,
                            loglevel: str = "WARNING") -> None:
    """
    Update a Relion .star (>=3.1 layout with data_optics / data_particles) using a CryoSPARC .cs.
    The two inputs must address the same set of particles (same length). If not, raise an error.

    Columns updated per particle (only if present in the STAR header):
        _rlnAngleRot        #2
        _rlnAngleTilt       #3
        _rlnAnglePsi        #4
        _rlnOriginXAngst    #5
        _rlnOriginYAngst    #6
        _rlnDefocusU        #7
        _rlnDefocusV        #8
        _rlnDefocusAngle    #9
        _rlnPhaseShift      #10
        _rlnCtfBfactor      #11
        _rlnOpticsGroup     #12
        _rlnRandomSubset    #13
        _rlnClassNumber     #14
    """
    import logging
    from scipy.spatial.transform import Rotation

    log = logging.getLogger("update_star_from_csparc")
    if not log.handlers:
        handler = logging.StreamHandler()
        log.addHandler(handler)

    if isinstance(loglevel, str):
        loglevel = loglevel.upper()
    log.setLevel(loglevel)

    if not csfile.endswith(".cs"):
        raise ValueError("First input must be a CryoSPARC .cs file")
    if not starfile_in.endswith(".star"):
        raise ValueError("Second input must be a Relion .star file")

    # 1) Load CryoSPARC .cs and derive quantities identically to csparc2star
    cs = np.load(csfile, max_header_size=100000)
    N = cs.shape[0]

    # micrograph mapping for optics group (same rule as csparc2star)
    raw_paths = pd.Series(cs["blob/path"].astype(str), name="_raw_path").str.lstrip(">")
    micro = raw_paths.str.split(pat="/extract", n=1).str[0]
    micro.name = "_micrograph"
    raw_optics = pd.DataFrame({"_micrograph": micro})
    optics_unique = raw_optics.drop_duplicates(subset="_micrograph", keep="first").reset_index(drop=True)
    optics_unique["_rlnOpticsGroup"] = np.arange(1, len(optics_unique) + 1)
    micro_to_group = dict(zip(optics_unique["_micrograph"], optics_unique["_rlnOpticsGroup"]))
    optics_group_per_particle = micro.map(micro_to_group).astype(int).values

    # pose → rotation matrices
    pose_arr = np.stack(cs["alignments3D/pose"])
    Kdim = pose_arr.shape[1]

    def _expmap_to_R(v: np.ndarray) -> np.ndarray:
        θ = np.linalg.norm(v)
        if θ < 1e-8:
            return np.eye(3)
        k = v / θ
        K = np.array([[   0, -k[2],  k[1]],
                      [ k[2],    0, -k[0]],
                      [-k[1], k[0],     0]])
        return np.eye(3) + np.sin(θ)*K + (1-np.cos(θ))*(K @ K)

    def _R_to_zyz(Rm: np.ndarray):
        β = np.arccos(np.clip(Rm[2,2], -1, 1))
        α = np.arctan2(Rm[1,2], Rm[0,2])
        γ = np.arctan2(Rm[2,1], -Rm[2,0])
        return np.degrees([α, β, γ])

    if Kdim == 3:
        Rmats = [_expmap_to_R(v) for v in pose_arr]
    elif Kdim == 4:
        q = pose_arr[:, [1,2,3,0]]  # x,y,z,w
        Rmats = Rotation.from_quat(q).as_matrix()
    elif Kdim == 9:
        Rmats = pose_arr.reshape(-1, 3, 3)
    else:
        raise ValueError(f"Unrecognised pose length {Kdim}")

    eulers = np.vstack([_R_to_zyz(Rm) for Rm in Rmats])
    phi, theta, psi = eulers[:,0], eulers[:,1], eulers[:,2]

    shifts = np.stack(cs["alignments3D/shift"])
    psize3 = cs["alignments3D/psize_A"]
    orig_x = shifts[:,0] * psize3
    orig_y = shifts[:,1] * psize3

    # derived arrays (same tags used by csparc2star)
    new_values = {
        "_rlnAngleRot":       phi,
        "_rlnAngleTilt":      theta,
        "_rlnAnglePsi":       psi,
        "_rlnOriginXAngst":   orig_x,
        "_rlnOriginYAngst":   orig_y,
        "_rlnDefocusU":       cs["ctf/df1_A"],
        "_rlnDefocusV":       cs["ctf/df2_A"],
        "_rlnDefocusAngle":   np.degrees(cs["ctf/df_angle_rad"]),
        "_rlnPhaseShift":     np.degrees(cs["ctf/phase_shift_rad"]),
        "_rlnCtfBfactor":     cs["ctf/bfactor"],
        "_rlnRandomSubset": (
            cs["alignments3D/split"].astype(int) + 1
            if "alignments3D/split" in cs.dtype.names
            else np.ones(N, int)
        ),
        "_rlnClassNumber": (
            cs["alignments3D/class"].astype(int) + 1
            if "alignments3D/class" in cs.dtype.names
            else np.ones(N, int)
        ),
        "_rlnOpticsGroup":    optics_group_per_particle,  # updated only if present in STAR
    }

    # 2) Read STAR, split into optics block (verbatim) and particles block
    with open(starfile_in, "r") as f:
        lines = f.readlines()

    # Track optics and particle blocks
    optics_block = []
    particles_block = []

    section = None
    for ln in lines:
        s = ln.strip()
        if s.startswith("data_optics"):
            section = "optics"
        elif s.startswith("data_particles"):
            section = "particles"

        if section == "optics":
            optics_block.append(ln)
        elif section == "particles":
            particles_block.append(ln)

    if not particles_block:
        raise ValueError("Could not find data_particles section in STAR file.")

    # 3) Parse particle header → column index map; locate first data row
    header_map = {}      # tag -> zero-based column index
    particle_header_rows = []
    first_data_idx = None

    for i, ln in enumerate(particles_block):
        st = ln.strip()
        if st.startswith("_"):
            parts = st.split()
            if len(parts) >= 2 and parts[1].startswith("#"):
                tag = parts[0]
                idx = int(parts[1][1:]) - 1
                header_map[tag] = idx
                particle_header_rows.append(ln)
        elif (not st) or st == "loop_" or st.startswith("data_particles"):
            continue
        else:
            # first non-header line inside particles block is data start
            first_data_idx = i
            break

    if first_data_idx is None:
        raise ValueError("Particles header found but no data rows detected in STAR file.")

    # 4) Load particle rows (robust to blank lines and full-line comments)
    #    - ignore empty/whitespace-only lines
    #    - ignore lines starting with '#' or ';'
    #    - ignore lines that don't have at least as many columns as in the header
    #      (protects against stray lines)
    min_cols = max(header_map.values()) + 1 if header_map else 0

    particle_table = []
    for ln in particles_block[first_data_idx:]:
        st = ln.strip()
        if not st:
            continue                       # skip blank line
        if st.startswith("#") or st.startswith(";"):
            continue                       # skip full-line comment
        if st.startswith("data_"):
            break                          # safety: stop at any next data_ section
        if st.startswith("loop_") or st.startswith("_"):
            # safety: if a second header appears by mistake, stop before duplicating
            break
        row = st.split()
        if len(row) < min_cols:
            # Pad missing columns instead of dropping; keeps row count stable.
            row = row + ([""] * (min_cols - len(row)))
        particle_table.append(row)

    # 5) Apply updates only for columns present in STAR
    #    Strategy:
    #      - Prefer matching by _rlnImageName (robust to row order and count differences).
    #      - Update only particles found in BOTH CS and STAR.
    #      - Leave unmatched STAR particles unchanged; log a warning.

    img_col = header_map.get("_rlnImageName", None)
    if img_col is None:
        # Fallback: if we cannot match by name, require strict row identity.
        if len(particle_table) != N:
            raise ValueError(
                f"Mismatch in particle number: CS has {N}, STAR has {len(particle_table)} "
                "and STAR has no _rlnImageName column to match on."
            )

        updated = 0
        for i in range(N):
            for tag, arr in new_values.items():
                if tag not in header_map:
                    continue
                j = header_map[tag]
                val = arr[i]
                if isinstance(val, (np.floating, float)):
                    particle_table[i][j] = f"{float(val):.6f}"
                elif isinstance(val, (np.integer, int)):
                    particle_table[i][j] = str(int(val))
                else:
                    try:
                        fv = float(val)
                        particle_table[i][j] = f"{fv:.6f}"
                    except Exception:
                        particle_table[i][j] = str(val)
            updated += 1

        log.info(f"Updated {updated}/{N} particles by row index (no _rlnImageName in STAR).")

    else:
        def _normalise_stack_basename(p: str) -> str:
            """
            Convert a stack path to a stable basename for matching.

            CryoSPARC often rewrites stacks as:
              J40/imported/015980123305114813081_HTT_46Q_stack.mrcs

            while RELION STAR may use:
              HTT_46Q_stack.mrcs

            We normalise by:
              - basename()
              - strip a long leading digit prefix followed by '_' (>=8 digits)
            """
            base = os.path.basename(p.strip())
            base = re.sub(r"^\d{8,}_", "", base)
            return base

        def _key_from_rln_image_name(s: str):
            """
            Return (particle_index_1based, normalised_stack_basename) from _rlnImageName.
            """
            st = str(s).strip()
            if "@" not in st:
                return None
            idx_str, stack = st.split("@", 1)
            try:
                idx_1based = int(idx_str)
            except Exception:
                return None
            return (idx_1based, _normalise_stack_basename(stack))

        # Build CS lookup: (idx_1based, stack_basename) -> cs_row_index
        cs_keys = {}
        for i in range(N):
            idx_1based = int(cs["blob/idx"][i]) + 1
            stack_path = str(raw_paths.iat[i])
            key = (idx_1based, _normalise_stack_basename(stack_path))
            # keep first occurrence if duplicates exist
            if key not in cs_keys:
                cs_keys[key] = i

        updated = 0
        missing_in_cs = 0

        for star_i in range(len(particle_table)):
            key = _key_from_rln_image_name(particle_table[star_i][img_col])
            if key is None:
                continue
            cs_i = cs_keys.get(key, None)
            if cs_i is None:
                missing_in_cs += 1
                continue

            # Update STAR row star_i using CS row cs_i
            for tag, arr in new_values.items():
                if tag not in header_map:
                    continue
                j = header_map[tag]
                val = arr[cs_i]
                if isinstance(val, (np.floating, float)):
                    particle_table[star_i][j] = f"{float(val):.6f}"
                elif isinstance(val, (np.integer, int)):
                    particle_table[star_i][j] = str(int(val))
                else:
                    try:
                        fv = float(val)
                        particle_table[star_i][j] = f"{fv:.6f}"
                    except Exception:
                        particle_table[star_i][j] = str(val)

            updated += 1

        if updated == 0:
            raise ValueError(
                "No particles could be matched between CS and STAR via _rlnImageName "
                "after normalising stack basenames. Check that CS was generated from "
                "the same particle stack and that _rlnImageName is consistent."
            )

        if updated != len(particle_table) or updated != N:
            log.warning(
                f"Updated {updated} STAR particles using CS matches; "
                f"CS has {N} particles, STAR has {len(particle_table)} particles; "
                f"{missing_in_cs} STAR particles had no match in CS and were left unchanged."
            )


    # 6) Compute column widths (max of header name and any data entry)
    #    This retains neat alignment similar to csparc2star writer.
    # Determine number of columns in the particle table from the longest row
    num_cols = max(len(r) for r in particle_table)
    col_widths = [0] * num_cols

    # include header tags in width
    for tag, idx in header_map.items():
        if 0 <= idx < num_cols:
            col_widths[idx] = max(col_widths[idx], len(tag))

    # include data
    for row in particle_table:
        for j, val in enumerate(row):
            if j < num_cols:
                col_widths[j] = max(col_widths[j], len(val))

    # 7) Reconstruct output .star: optics block verbatim; particles: header+formatted rows
    with open(starfile_out, "w") as f:
        # optics as-is
        for ln in optics_block:
            f.write(ln)

        # ensure a blank line between sections (only if optics existed)
        if optics_block and not optics_block[-1].endswith("\n"):
            f.write("\n")

        # particles header (verbatim up to first_data_idx)
        for ln in particles_block[:first_data_idx]:
            f.write(ln)

        # formatted particle rows: left-justify ImageName if present, right-justify others
        # detect image-name column, if any
        img_col = header_map.get("_rlnImageName", None)

        for row in particle_table:
            pieces = []
            for j in range(num_cols):
                val = row[j] if j < len(row) else ""
                width = col_widths[j]
                if img_col is not None and j == img_col:
                    pieces.append(val.ljust(width))
                else:
                    pieces.append(val.rjust(width))
            f.write(" ".join(pieces) + "\n")

    log.info(f"Wrote updated STAR with {N} particles → {starfile_out}")


