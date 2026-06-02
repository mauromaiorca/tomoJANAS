#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# File: janas_cmd_caller.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology


"""
Module: janas_cmd_caller.py

Defines the command-line interface for the JANAS toolkit, providing:
- Particle scoring and ranking commands
- Selection and classification utilities
- Reconstruction script generation
- STAR file manipulation (extract, duplicate removal, parameter editing)
- Diagnostic plotting of Euler angle distributions

Commands map directly to core JANAS methods:
  - assessParticles.ParticleVsReprojectionScores
  - janas_core.EqualizedParticlesRank
  - starHandler utilities
"""

# Standard library
import argparse
import os.path
import stat
from os import PathLike

import json
import sys
import getpass

# Third-party
import numpy as np
import scipy.stats as stats

# Local
import janas.janas_core as janas_core
from janas import assessParticles
from janas import starHandler
from janas.version import get_version

# from janas import utils


janas_parser = argparse.ArgumentParser(
    prog="janas",
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


##########################################
##########################################
##### janas_rankParticles
##########################################
janas_rankParticles = command.add_parser(
    "rankParticles",
    description="Normalize Score in star file tag",
    help="normalize score",
)
janas_rankParticles.add_argument(
    "--i", required=True, type=str, help="input star file to normalize"
)
janas_rankParticles.add_argument(
    "--o", required=False, type=str, help="output star file"
)
janas_rankParticles.add_argument(
    "--tag", required=True, type=str, help="tag score to normalize"
)
janas_rankParticles.add_argument(
    "--tagOut", required=True, type=str, help="output tag with normalized score"
)
janas_rankParticles.add_argument(
    "--avgViews",
    required=False,
    default="",
    type=str,
    help="mrc stack with averaged views",
)
janas_rankParticles.add_argument(
    "--views",
    required=False,
    default=50,
    type=int,
    help="number of views partitioned for normalization",
)


def rankParticles(args):
    views = args.views
    # tagPrefix=args.prefix
    if not os.path.isfile(args.i):
        print('ERROR: file "', args.i, '" not existing')
    if args.o is None:
        args.o = args.i
    elif not os.path.exists(os.path.dirname(args.o)):
        try:
            os.makedirs(os.path.dirname(args.o))
        except OSError as e:
            print(f"Error: {e.strerror}")

    columns = starHandler.header_columns(args.i)
    # columns=starHandler.readColumns(args.i, [args.tagOut], sortColumnsNameOrder=True)
    tagOut_view = args.tagOut
    if args.tagOut in columns:
        starHandler.removeColumnsTagsStartingWith(args.i, args.o, args.tagOut)
    if tagOut_view in columns:
        starHandler.removeColumnsTagsStartingWith(args.i, args.o, tagOut_view)

    # print ('args.tag=',args.tag)
    coordinates = starHandler.readColumns(
        args.i, ["_rlnAngleRot", "_rlnAngleTilt", "_rlnRandomSubset", args.tag]
    )
    phiListParticle = coordinates["_rlnAngleRot"].tolist()
    thetaListParticle = coordinates["_rlnAngleTilt"].tolist()
    scores = coordinates[args.tag].tolist()
    randomSubset = coordinates["_rlnRandomSubset"].tolist()
    rankedParticles = janas_core.EqualizedParticlesRank(
        phiListParticle, thetaListParticle, scores, randomSubset, int(views)
    )
    starHandler.addColumns(args.i, args.o, [args.tagOut], [rankedParticles])
    viewsParticles = janas_core.GetEulerClassGroup(
        phiListParticle, thetaListParticle, int(views)
    )
    starHandler.addColumns(args.i, args.o, [tagOut_view], [viewsParticles])


##########################################
##########################################
##### janas_scoreParticles
##########################################
janas_scoreParticles = command.add_parser(
    "scoreParticles",
    description="score Particles",
    help="score Particles described in star file",
)
janas_scoreParticles.add_argument(
    "--i", required=True, type=str, help="input star file with particles to score"
)
janas_scoreParticles.add_argument(
    "--mask", required=True, type=str, help="mrc volumetric file with mask"
)
janas_scoreParticles.add_argument(
    "--maskToSubtract",
    required=False,
    default=None,
    type=str,
    help=(
        "optional 3D MRC mask defining the region whose signal should be "
        "subtracted in the particle images before scoring "
        "(projection-weighted, RELION-like)."
    ),
)
janas_scoreParticles.add_argument(
    "--map",
    required=True,
    type=str,
    help="mrc volumetric file with reference map (or first half reference map if map2 is given)",
)
janas_scoreParticles.add_argument(
    "--map2",
    required=False,
    type=str,
    help="mrc volumetric file with second half map as reference",
)
janas_scoreParticles.add_argument("--apix", required=True, type=float, help="angpix")
janas_scoreParticles.add_argument(
    "--sigma", required=False, type=float, default=1, help="sigma blurring in pixels"
)
janas_scoreParticles.add_argument(
    "--selectionName",
    required=False,
    type=str,
    default="0",
    help="Name for the selection",
)
janas_scoreParticles.add_argument(
    "--o", required=False, default=None, type=str, help="output file with scores"
)
janas_scoreParticles.add_argument(
    "--mpi",
    required=False,
    default=4,
    type=int,
    help="number of mpi parallel processes",
)
janas_scoreParticles.add_argument(
    "--rank",
    required=False,
    default="350",
    type=int,
    help="rank particles using a given number of Euler viewss",
)
janas_scoreParticles.add_argument(
    "--ctf-mode",
    required=False,
    type=str,
    choices=["modulate", "phaseflip", "wiener"],
    default="phaseflip",
    help=(
        "CTF application mode for particle scoring: "
        "'modulate' (multiply by full CTF), "
        "'phaseflip' (sign of CTF, default), or "
        "'wiener' (CTF / (CTF^2 + 0.1)). "
        "Default: phaseflip."
    ),
)
janas_scoreParticles.add_argument(
    "--scoreMethod",
    required=False,
    type=str,
    default="SCI",
    choices=["SCI", "CC", "MI", "SSIM", "PSNR"],
    help="scoring method to use (default: SCI)",
)


def scoreParticles(args):
    inputFile = args.i
    outputFile = args.o
    templateMask = args.mask
    templateMap = args.map
    templateMap2 = args.map2
    subtractionMask = args.maskToSubtract
    rankViews = args.rank
    # tagPrefix=args.prefix
    if not os.path.isfile(inputFile):
        print('ERROR: file "', inputFile, '" not existing')
        exit()
    if not os.path.isfile(templateMask):
        print('ERROR: mask file "', templateMask, '" not existing')
        exit()
    if subtractionMask is not None and len(subtractionMask) > 0:
        if not os.path.isfile(subtractionMask):
            print('ERROR: maskToSubtract file "', subtractionMask, '" not existing')
            exit()
    else:
        subtractionMask = None
    if not os.path.isfile(templateMap):
        print('ERROR: map file "', templateMap, '" not existing')
        exit()
    if outputFile == None:
        outputFile = inputFile

    score_method = getattr(args, "scoreMethod", "SCI")
    if score_method is None or str(score_method).strip() == "":
        score_method = "SCI"
    score_method = str(score_method).upper()

    selection_name = args.selectionName
    if selection_name is None or str(selection_name).strip() == "":
        selection_name = "0"
    else:
        selection_name = str(selection_name).strip()
        if selection_name.startswith("selection_"):
            selection_name = selection_name[len("selection_"):]

    tag = (
        f"_janas_{score_method}__"
        + "{:.2f}".format(args.sigma)
        + "_scored_selection_"
        + selection_name
    )
    listScoresTags = [tag]
    print(listScoresTags)



    if args.map2 is None:
        assessParticles.ParticleVsReprojectionScores(
            inputFile,
            outputFile,
            templateMap,
            templateMask,
            angpix=args.apix,
            numProcesses=args.mpi,
            listScoresTags=listScoresTags,
            doCTF=True,
            ctfMode=args.ctf_mode,
            referenceSubtractionMask=subtractionMask,
        )
    else:
        if not os.path.isfile(templateMap2):
            print('ERROR: map file "', templateMap2, '" not existing')
            exit()
        assessParticles.ParticleVsReprojectionScores_HalfMaps(
            inputFile,
            outputFile,
            templateMap,
            templateMap2,
            templateMask,
            angpix=args.apix,
            numProcesses=args.mpi,
            listScoresTags=listScoresTags,
            doCTF=True,
            ctfMode=args.ctf_mode,
            referenceSubtractionMask=subtractionMask,
        )

    # ranking
    rankingTag = tag + "_norm" + str(rankViews)
    coordinates = starHandler.readColumns(
        outputFile, ["_rlnAngleRot", "_rlnAngleTilt", "_rlnRandomSubset", tag]
    )
    phiListParticle = coordinates["_rlnAngleRot"].tolist()
    thetaListParticle = coordinates["_rlnAngleTilt"].tolist()
    scores = coordinates[tag].tolist()
    randomSubset = coordinates["_rlnRandomSubset"].tolist()
    rankedParticles = janas_core.EqualizedParticlesRank(
        phiListParticle, thetaListParticle, scores, randomSubset, int(rankViews)
    )
    columns = starHandler.header_columns(outputFile)
    if rankingTag in columns:
        starHandler.removeColumnsTagsStartingWith(outputFile, outputFile, rankingTag)
    starHandler.addColumns(outputFile, outputFile, [rankingTag], [rankedParticles])


##########################################
##########################################
##### scoreClassify
##########################################
janas_scoreClassify = command.add_parser(
    "scoreClassify",
    description="score classify",
    help="classify particles based on the measured score",
)
janas_scoreClassify.add_argument(
    "--i", required=True, type=str, help="input star file with scored particles"
)
janas_scoreClassify.add_argument(
    "--classLabelOut",
    required=False,
    default="_rlnClassNumber",
    type=str,
    help="label with the class to be assigned, default=_rlnClassNumber",
)
janas_scoreClassify.add_argument(
    "--classScoredLabels",
    required=False,
    nargs="+",
    type=str,
    help="labels with scores for each reference class",
)
janas_scoreClassify.add_argument(
    "--o",
    required=False,
    default=None,
    type=str,
    help="output star file with updated classification",
)


def scoreClassify(args):
    inputFile = args.i
    outputFile = args.o
    if not os.path.isfile(inputFile):
        print('ERROR: file "', inputFile, '" not existing')
        exit()
    if outputFile == None:
        outputFile = inputFile
    print("classLabelOut=", args.classLabelOut)
    print("classScoredLabels=", args.classScoredLabels)
    class_columns = args.classScoredLabels
    columns = starHandler.header_columns(args.i)
    if class_columns == None:
        print(
            "WARNING: --classScoredLabels not properly given by the user, using the array labels ending with _class and a classID"
        )
        class_columns = [
            col
            for col in columns
            if col.split("_")[-1].startswith("class")
            and col.split("_")[-1][5:].isdigit()
        ]
    if class_columns == None:
        print(
            "ERROR: no suitable classes given by the user, please check the star file"
        )
        exit()
    # print ("class_columns=",class_columns)
    import pandas as pd

    dataDF = starHandler.readColumns(args.i, class_columns)
    # print ("dataDF=",dataDF)
    max_column = dataDF.idxmax(axis=1)
    class_numbers = max_column.str.extract("class(\d+)$")[0]
    # if all the values are not positive, assign class -1, as the program can't decide
    all_non_positive = (dataDF <= 0).all(axis=1)
    class_numbers[all_non_positive] = "-1"
    resultDF = pd.DataFrame({args.classLabelOut: class_numbers})
    # print(resultDF)

    columnsToRemove = [item for item in columns if item in [args.classLabelOut]]
    if not columnsToRemove == []:
        starHandler.removeColumns(inputFile, outputFile, columnsToRemove)
    starHandler.addDataframeColumns(outputFile, outputFile, columnsToRemove, resultDF)


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
    "--tagRank", required=True, type=str, help="tag for the ranked particles"
)
janas_produce_reconstructions_script.add_argument(
    "--mask", required=True, type=str, help="mask to use for evaluation"
)
# janas_produce_reconstructions_script.add_argument("--numReconstructions", required=False, type=int, default="5",  help="Number of Reconstructions")
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

    if args.manualParticleSubsets is not None:
        print("Manual Particle Subset Selection")
        numParticlesList = args.manualParticleSubsets.split(",")

    print(numParticlesList)

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
        "norm_" + str(os.path.split(args.outDir)[-1]) + "_best${numParticles}",
    )
    reconstruction_command += (
        "    janas selectBestRanked --i "
        + args.i
        + " --o "
        + outFile
        + ".star --num  ${numParticles} --tag "
        + args.tagRank
        + "\n"
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

    #######################
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
            mpirun --np 28 relion_postprocess_mpi --i ${file_basename}_recH1.mrc --i2 ${file_basename}_recH2.mrc --o  ${file_basename} --locres --locres_thresholdFSC 0.5
            rm ${file_basename}_locres_fscs.star
            rm ${file_basename}_locres_filtered.mrc
        fi
}

#IMOD example for cropping, we are using janas_utils maskedCrop insted
#crop_image_IMOD() {
#  local imageIn="$1"
#  local imageOut="$2"
#  trimvol -x 100,300 -y 100,300 -z 100,300  "$imageIn" "$imageOut"
#}

    \n\n"""
    numParticlesListStr = ",".join(map(str, numParticlesList))
    locres_command += "numParticlesCsv=" + numParticlesListStr + "\n"
    locres_command += """
for numParticles in $(echo $numParticlesCsv | sed "s/,/ /g")
do\n"""

    outFile = os.path.join(
        args.outDir,
        "norm_" + str(os.path.split(args.outDir)[-1]) + "_best${numParticles}",
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
    reconstruction_file_path = os.path.join(args.outDir, "script_reconstructions.sh")
    with open(reconstruction_file_path, "a") as f:
        f.write(locres_command)

    #######################
    #####ASSESS LOCRES
    assess_command = """
##############################
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
    result_filename = os.path.join(args.outDir, "bestRanked_locres_values.csv")
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
        "norm_" + str(os.path.split(args.outDir)[-1]) + "_best${numParticles}",
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
    reconstruction_file_path = os.path.join(args.outDir, "script_reconstructions.sh")
    with open(reconstruction_file_path, "a") as f:
        f.write(assess_command)

    #######################
    #####SELECT OPTIMAL SUBSET
    select_command = """
##############################
#######  SELECT OPTIMAL SUBSET
"""
    numParticlesListStr = ",".join(map(str, numParticlesList))
    result_filename = os.path.join(args.outDir, "bestRanked_locres_values.csv")


#    with open(reconstruction_file_path, 'a') as f:
#        f.write(assess_command)


def find_janas_ranking_tag(elements):
    """
    Search for a column name that starts with '_janas_' or '_emprove_' (backward compat) and ends with '_norm' followed by a number.
    """
    import re

    for element in elements:
        if re.match(r"(_janas_|_emprove_).*_norm\d+$", element):
            return element
    return None


janas_selectWorstRanked = command.add_parser(
    "selectWorstRanked",
    description="selectWorstRanked",
    help="select worst scoring particles",
)
janas_selectWorstRanked.add_argument(
    "--i", required=True, type=str, help="input file"
)
janas_selectWorstRanked.add_argument(
    "--num", required=True, type=str, help="number of ranked particles to select"
)
janas_selectWorstRanked.add_argument(
    "--tag", required=False, type=str, help="tag to select worst"
)
janas_selectWorstRanked.add_argument(
    "--o", required=False, default=None, type=str, help="output file"
)


def selectWorstRanked(args):
    inputFile = args.i
    tag = str(args.tag)
    outputFile = args.o
    if outputFile == None:
        outputFile = inputFile
    if not os.path.isfile(inputFile):
        print('ERROR: file "', inputFile, '" not existing')
        exit(0)
    columns = starHandler.header_columns(inputFile)
    if not args.tag:
        args.tag = args.tag if args.tag else find_janas_ranking_tag(columns)
    if args.tag not in columns:
        print("target tag [", args.tag, "] not in ", inputFile)
        exit(0)
    print("target tag=", args.tag)
    starHandler.extractWorst(inputFile, outputFile, int(args.num), args.tag)


janas_selectRandom = command.add_parser(
    "selectRandom", description="selectRandom", help="select random particles"
)
janas_selectRandom.add_argument("--i", required=True, type=str, help="input file")
janas_selectRandom.add_argument(
    "--num", required=True, type=str, help="number of particles to select"
)
janas_selectRandom.add_argument(
    "--o", required=False, default=None, type=str, help="output file"
)


def selectRandom(args):
    inputFile = args.i
    outputFile = args.o
    if outputFile == None:
        outputFile = inputFile
    if not os.path.isfile(inputFile):
        print('ERROR: file "', inputFile, '" not existing')
        exit(0)
    starHandler.extractRandom(inputFile, outputFile, int(args.num))


janas_selectBestRanked = command.add_parser(
    "selectBestRanked", description="selectBestRanked", help="select best particles"
)
janas_selectBestRanked.add_argument("--i", required=True, type=str, help="input file")
janas_selectBestRanked.add_argument(
    "--num", required=True, type=str, help="number of ranked particles to select"
)
janas_selectBestRanked.add_argument(
    "--tag", required=False, type=str, help="tag to select best"
)
janas_selectBestRanked.add_argument(
    "--o", required=False, default=None, type=str, help="output file"
)
janas_selectBestRanked.add_argument(
    "--exact",
    required=False,
    action="store_true",
    help="Fail with a non-zero exit code if the input STAR file contains fewer particles than --num.",
)


def selectBestRanked(args):
    inputFile = args.i
    outputFile = args.o
    if outputFile == None:
        outputFile = inputFile
    if not os.path.isfile(inputFile):
        print('ERROR: file "', inputFile, '" not existing')
        exit(0)
    columns = starHandler.header_columns(inputFile)
    if not args.tag:
        args.tag = args.tag if args.tag else find_janas_ranking_tag(columns)
    if args.tag not in columns:
        print("target tag [", args.tag, "] not in ", inputFile)
        exit(0)
    print("target tag=", args.tag)
    try:
        starHandler.extractBest(
            inputFile, outputFile, int(args.num), args.tag, exact=bool(getattr(args, "exact", False))
        )
    except ValueError as exc:
        print("ERROR:", exc)
        sys.exit(2)


janas_removeDuplicates = command.add_parser(
    "removeDuplicates",
    description="removeDuplicates",
    help="remove duplicate particles from star file",
)
janas_removeDuplicates.add_argument("--i", required=True, type=str, help="input file")
janas_removeDuplicates.add_argument(
    "--o", required=False, default=None, type=str, help="output file"
)


def removeDuplicates(args):
    inputFile = args.i
    outputFile = args.o
    if outputFile == None:
        outputFile = inputFile
    if not os.path.isfile(inputFile):
        print('ERROR: file "', inputFile, '" not existing')
        exit(0)
    starHandler.removeStarDuplicates(inputFile, outputFile)


janas_assignClassName = command.add_parser(
    "assignClassName", description="assignClassName", help="assign class name"
)
janas_assignClassName.add_argument("--i", required=True, type=str, help="input file")
janas_assignClassName.add_argument(
    "--className", required=True, type=str, help="Class Name (usually a number)"
)
janas_assignClassName.add_argument(
    "--o", required=False, default=None, type=str, help="output file"
)


def assignClassName(args):
    inputFile = args.i
    outputFile = args.o
    if outputFile == None:
        outputFile = inputFile
    if not os.path.isfile(inputFile):
        print('ERROR: file "', inputFile, '" not existing')
        exit(0)
    starHandler.AssignStarClassName(inputFile, outputFile, args.className)


janas_changeParamValue = command.add_parser(
    "changeParamValue",
    description="changeParamValue",
    help="change value of a certain column (parameter) in a star file",
)
janas_changeParamValue.add_argument("--i", required=True, type=str, help="input file")
janas_changeParamValue.add_argument(
    "--columnName",
    required=True,
    type=str,
    help="Name of The column(parameter) to be changed",
)
janas_changeParamValue.add_argument(
    "--newValue", required=True, type=str, help="new value to be inserted"
)
janas_changeParamValue.add_argument(
    "--o", required=False, default=None, type=str, help="output file"
)


def changeParamValue(args):
    inputFile = args.i
    outputFile = args.o
    if outputFile == None:
        outputFile = inputFile
    if not os.path.isfile(inputFile):
        print('ERROR: file "', inputFile, '" not existing')
        exit(0)
    starHandler.changeStarFileParamValue(
        inputFile, outputFile, args.columnName, args.newValue
    )


janas_removeParamValue = command.add_parser(
    "removeParamValue",
    description="removeParamValue",
    help="remove particles in a star file where a parameter is of a certain value",
)
janas_removeParamValue.add_argument("--i", required=True, type=str, help="input file")
janas_removeParamValue.add_argument(
    "--columnName",
    required=True,
    type=str,
    help="Name of The column(parameter) to be removed",
)
janas_removeParamValue.add_argument(
    "--value", required=True, type=str, help="value to be removed"
)
janas_removeParamValue.add_argument(
    "--o", required=False, default=None, type=str, help="output file"
)


def removeParamValue(args):
    inputFile = args.i
    outputFile = args.o
    if outputFile == None:
        outputFile = inputFile
    if not os.path.isfile(inputFile):
        print('ERROR: file "', inputFile, '" not existing')
        exit(0)
    starHandler.deleteStarFileParamValue(
        inputFile, outputFile, args.columnName, args.value
    )


janas_plotParamValue = command.add_parser(
    "plotParamValue",
    description="plotParamValue",
    help="plot values in a star file at a certain column",
)
janas_plotParamValue.add_argument("--i", required=True, type=str, help="input file")
janas_plotParamValue.add_argument(
    "--columnName",
    required=True,
    type=str,
    help="Name of The column(parameter) to be changed",
)


def plotParamValue(args):
    inputFile = args.i
    if not os.path.isfile(inputFile):
        print('ERROR: file "', inputFile, '" not existing')
        exit(0)
    starHandler.plotStarFileParamValueInteractive(inputFile, args.columnName)


###################
##  DISPLAY


def _str2bool(value):
    """
    Robust string-to-bool for argparse arguments.

    ``argparse(type=bool)`` is broken: ``bool("False")`` returns ``True``
    because any non-empty string is truthy. As a result something like
    ``--show False`` was silently turning the flag back on. This helper
    accepts the usual textual variants and rejects everything else with
    a clear argparse error.
    """
    if isinstance(value, bool):
        return value
    s = str(value).strip().lower()
    if s in ("true", "t", "yes", "y", "1", "on"):
        return True
    if s in ("false", "f", "no", "n", "0", "off"):
        return False
    import argparse
    raise argparse.ArgumentTypeError(
        f"Expected a boolean value (true/false), got: {value!r}"
    )


def plotRoundEulerHist(
    Phi, Theta, titlePlot, maxValue, numBins,
    outImage: PathLike = None, toShow=True, fontScale=1.0,
):
    """Render the (phi, theta) histogram on a Mollweide projection.

    When ``toShow`` is False, the figure is built via
    :class:`matplotlib.figure.Figure` directly, which does not initialise
    any matplotlib GUI backend. This makes ``janas eulerHist --outImage
    out.png --show False`` work cleanly on headless nodes and on SSH
    sessions with broken X11 forwarding (same pattern used in
    janas_optimizer's predict_min_particles).
    """
    # Phi=Phi.astype(np.float)
    # Theta=Theta.astype(np.float)
    eulers1 = (Phi.to_numpy()).reshape((len(Phi),)) * np.pi / 180.0
    eulers2 = (Theta.to_numpy()).reshape((len(Theta),)) * np.pi / 180.0
    eulers1 = np.mod(eulers1, 2 * np.pi) - np.pi
    eulers2 = np.mod(eulers2, np.pi) - (np.pi / 2.0)

    H, xedges, yedges = np.histogram2d(
        eulers1,
        eulers2,
        numBins,
        range=[[-3.1415926535, 3.1415926535], [-1.570796, 1.570796]],
    )
    H = H.T

    if toShow:
        import matplotlib.pyplot as plt  # noqa: WPS433
        fig = plt.figure(figsize=(10, 5))
    else:
        from matplotlib.figure import Figure  # noqa: WPS433
        fig = Figure(figsize=(10, 5))

    ax = fig.add_subplot(
        111, title="pcolormesh: actual edges", aspect="equal", projection="mollweide"
    )
    # Font sizes are derived from matplotlib defaults and scaled by
    # ``fontScale`` so the dashboard can request a 2x-zoom when it
    # renders the histogram in a small card (where the default font
    # would be unreadable).
    try:
        fs = float(fontScale) if fontScale else 1.0
    except (TypeError, ValueError):
        fs = 1.0
    if fs <= 0:
        fs = 1.0
    title_fs = 14.0 * fs
    label_fs = 12.0 * fs
    tick_fs = 10.0 * fs

    ax.tick_params(
        axis="x",
        direction="out",
        length=6,
        width=4,
        colors="w",
        grid_color="w",
        grid_alpha=0.5,
        label1On=False,
    )
    ax.tick_params(axis="y", labelsize=tick_fs)
    X, Y = np.meshgrid(xedges, yedges)
    # Cast maxValue to float here so a CLI-supplied string works as vmax.
    try:
        vmax = float(maxValue) if maxValue is not None else None
    except (TypeError, ValueError):
        vmax = None
    pcm = ax.pcolormesh(X, Y, H, cmap="RdBu_r", vmin=-1, vmax=vmax)

    # colormaps https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
    cbar = fig.colorbar(pcm, ax=ax, extend="both")
    cbar.ax.tick_params(labelsize=tick_fs)
    ax.grid(color="w", linestyle=":", linewidth=1)
    ax.set_title(titlePlot, pad=20, fontweight="bold", fontsize=title_fs)
    # Use the axes-level setters (not pyplot's current-axes shortcuts) so
    # the labels are still applied when fig is a bare Figure().
    ax.set_xlabel(r"Rot Angles ($\phi$)", fontweight="bold", fontsize=label_fs)
    ax.set_ylabel(r"Tilt Angles ($\theta$)", fontweight="bold", fontsize=label_fs)

    if outImage:
        fig.savefig(outImage)

    if toShow:
        plt.show()
        plt.close(fig)
    # When not showing, the bare Figure is discarded on function return;
    # no explicit plt.close() is needed because it was never managed.


janas_eulerHist = command.add_parser(
    "eulerHist",
    description="produce euler histogram for the data",
    help="produce euler histogram for the data",
)
janas_eulerHist.add_argument("--i", required=True, type=str, help="input star file")
janas_eulerHist.add_argument(
    "--title",
    required=False,
    default="phi/theta histogram",
    type=str,
    help="Title of the plot",
)
janas_eulerHist.add_argument(
    "--maxValue",
    required=False,
    default=None,
    type=str,
    help="maximum value to display",
)
janas_eulerHist.add_argument(
    "--numBins",
    required=False,
    default=40,
    type=float,
    help="Number of bins for the histogram",
)
janas_eulerHist.add_argument(
    "--outImage",
    required=False,
    default=None,
    type=str,
    help="output file where storing image",
)
janas_eulerHist.add_argument(
    "--show",
    required=False,
    default=True,
    type=_str2bool,
    metavar="{true,false}",
    help="Whether or not to display the image on screen. Accepts true/false "
         "(also yes/no, 1/0, on/off; case-insensitive). Default: true. "
         "Set to false when only saving via --outImage on a headless node.",
)
janas_eulerHist.add_argument(
    "--fontScale",
    required=False,
    default=1.0,
    type=float,
    help="Multiplier applied to the title, axis-label and tick-label font "
         "sizes. Useful when saving a small PNG (e.g. for the progress "
         "dashboard) where the default text would be unreadable. "
         "Default: 1.0; pass 2.0 to double every font size.",
)


def eulerHist(args):
    if not os.path.isfile(args.i):
        print('ERROR: file "', args.i, '" not existing')
        exit(0)
    Phi = starHandler.readColumns(args.i, ["_rlnAngleRot"])
    Theta = starHandler.readColumns(args.i, ["_rlnAngleTilt"])
    plotRoundEulerHist(
        Phi, Theta, args.title, args.maxValue, args.numBins,
        args.outImage, args.show, args.fontScale,
    )


janas_deleteTags = command.add_parser(
    "deleteTags",
    description="Cleans StarTags",
    help="delete tags with a certain prefix",
)
janas_deleteTags.add_argument(
    "--i", required=True, type=str, help="input file for deleting tags"
)
janas_deleteTags.add_argument(
    "--o",
    required=False,
    default=None,
    type=str,
    help="output file with deleted tags file",
)
janas_deleteTags.add_argument(
    "--prefix",
    required=False,
    default="_scorem_",
    type=str,
    help="prefix for the tag to remove",
)


def deleteTags(args):
    inputFile = args.i
    outputFile = args.o
    tagPrefix = args.prefix
    if outputFile == None:
        outputFile = inputFile
    if not os.path.isfile(inputFile):
        print('ERROR: file "', inputFile, '" not existing')
        exit(0)
    starHandler.removeColumnsTagsStartingWith(inputFile, outputFile, tagPrefix)

##########################################
##########################################
##### janas_csparc_setup
##########################################
janas_csparc_setup = command.add_parser(
    "csparc_setup",
    description=(
        "Configure CryoSPARC integration for JANAS by generating a shell "
        "environment file with the required CRYOSPARC_* variables."
    ),
    help="configure CryoSPARC integration (license, host, ports, user, password)",
)

janas_csparc_setup.add_argument(
    "--license-id",
    required=False,
    type=str,
    help="CryoSPARC license ID (CRYOSPARC_LICENSE_ID).",
)
janas_csparc_setup.add_argument(
    "--host",
    required=False,
    type=str,
    help="CryoSPARC master hostname (CRYOSPARC_MASTER_HOSTNAME).",
)
janas_csparc_setup.add_argument(
    "--base-port",
    required=False,
    type=int,
    help="CryoSPARC base port (CRYOSPARC_BASE_PORT, default 39000).",
)
janas_csparc_setup.add_argument(
    "--email",
    required=False,
    type=str,
    help="CryoSPARC user e-mail (login e-mail).",
)
janas_csparc_setup.add_argument(
    "--password",
    required=False,
    type=str,
    help=(
        "CryoSPARC user password. If omitted, you will be prompted. "
        "NOTE: stored in plain text in the generated env file."
    ),
)
janas_csparc_setup.add_argument(
    "--config-dir",
    required=False,
    type=str,
    default=None,
    help=(
        "Directory where cryosparc_env.sh and cryosparc_config.json will "
        "be written. Default: ~/.janas"
    ),
)
def csparc_setup(args):
    """
    Configure CryoSPARC integration for JANAS.

    This will:
    - verify that the current Python environment is reasonably recent,
    - verify that 'cryosparc-tools' can connect to the CryoSPARC instance,
    - write a shell file (cryosparc_env.sh) with CRYOSPARC_* exports,
    - optionally write a JSON file with the same configuration.

    The idea is that users can then do:
        source ~/.janas/cryosparc_env.sh
        janas_utils csparc_nurefinement ...
    """
    # 1. Verifica versione di Python (non blocca, ma avvisa)
    if sys.version_info < (3, 10):
        print(
            "[csparc_setup] WARNING: Python >= 3.10 is recommended for JANAS "
            "and CryoSPARC integration (current: "
            f"{sys.version_info.major}.{sys.version_info.minor})."
        )

    # 2. Verifica che cryosparc-tools sia installato
    try:
        from cryosparc.tools import CryoSPARC
    except ImportError:
        print(
            "[csparc_setup] ERROR: 'cryosparc-tools' is not installed in this "
            "Python environment.\n"
            "Install it with:\n\n"
            "    pip install 'cryosparc-tools>=4.3,<5.0'\n"
        )
        sys.exit(1)

    # 3. Raccogli parametri, usando argomenti o variabili d'ambiente
    license_id = args.license_id or os.environ.get("CRYOSPARC_LICENSE_ID")
    if not license_id:
        license_id = input("CryoSPARC license ID (CRYOSPARC_LICENSE_ID): ").strip()
    if not license_id:
        print("[csparc_setup] ERROR: license ID is required.")
        sys.exit(1)

    host = args.host or os.environ.get("CRYOSPARC_MASTER_HOSTNAME", "localhost")
    base_port = args.base_port or int(os.environ.get("CRYOSPARC_BASE_PORT", "39000"))

    email = args.email or os.environ.get("CRYOSPARC_EMAIL")
    if not email:
        email = input("CryoSPARC user e-mail (login e-mail): ").strip()
    if not email:
        print("[csparc_setup] ERROR: CryoSPARC user e-mail is required.")
        sys.exit(1)

    password = args.password or os.environ.get("CRYOSPARC_PASSWORD")
    if not password:
        # Input senza eco a schermo
        password = getpass.getpass(
            "CryoSPARC user password (will be stored in plain text): "
        )
    if not password:
        print("[csparc_setup] ERROR: CryoSPARC password is required.")
        sys.exit(1)

    # 4. Testa la connessione a CryoSPARC
    print(
        f"[csparc_setup] Testing connection to CryoSPARC at "
        f"{host}:{base_port} as {email}..."
    )
    try:
        cs = CryoSPARC(
            license=license_id,
            host=host,
            base_port=base_port,
            email=email,
            password=password,
        )
        lanes = cs.cli.get_scheduler_lanes()
    except Exception as e:
        print(
            "[csparc_setup] ERROR: could not connect to CryoSPARC with the "
            "provided credentials and host.\n"
            f"Details: {e}"
        )
        sys.exit(1)

    num_lanes = len(lanes) if isinstance(lanes, (list, tuple)) else "unknown"
    print(f"[csparc_setup] Connection OK. CryoSPARC version: {num_lanes}")

    # 5. Determina directory di configurazione
    config_dir = (
        args.config_dir
        if args.config_dir is not None
        else os.path.join(os.path.expanduser("~"), ".janas")
    )
    os.makedirs(config_dir, exist_ok=True)

    # 6. Scrive il file di ambiente shell
    env_path = os.path.join(config_dir, "cryosparc_env.sh")
    with open(env_path, "w", encoding="utf-8") as f:
        f.write("#!/usr/bin/env bash\n")
        f.write("# Auto-generated by 'janas csparc_setup'\n")
        f.write(f'export CRYOSPARC_LICENSE_ID="{license_id}"\n')
        f.write(f'export CRYOSPARC_MASTER_HOSTNAME="{host}"\n')
        f.write(f"export CRYOSPARC_BASE_PORT={base_port}\n")
        f.write(f'export CRYOSPARC_EMAIL="{email}"\n')
        f.write(f'export CRYOSPARC_PASSWORD="{password}"\n')

    # Rende lo script eseguibile per l'utente
    try:
        os.chmod(env_path, os.stat(env_path).st_mode | stat.S_IXUSR)
    except PermissionError:
        pass

    # 7. Scrive anche un JSON di configurazione (opzionale, utile in futuro)
    cfg_path = os.path.join(config_dir, "cryosparc_config.json")
    cfg = {
        "license_id": license_id,
        "host": host,
        "base_port": base_port,
        "email": email,
        # per sicurezza si potrebbe omettere la password dal JSON,
        # dato che è già nello shell env file.
    }
    with open(cfg_path, "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=2)

    print()
    print("[csparc_setup] Configuration files written to:")
    print(f"    {env_path}")
    print(f"    {cfg_path}")
    print()
    print("To enable CryoSPARC integration in a shell, run:")
    print(f"    source {env_path}")
    print("before running JANAS commands that use CryoSPARC, e.g.:")
    print("    janas_utils csparc_nurefinement ...")
    print()
    print("[csparc_setup] Done.")


def main(command_line=None):
    args = janas_parser.parse_args(command_line)
    if args.command == "scoreParticles":
        scoreParticles(args)
    elif args.command == "rankParticles":
        rankParticles(args)
    elif args.command == "scoreClassify":
        scoreClassify(args)
    elif args.command == "produce_reconstructions_script":
        produce_reconstructions_script(args)
    elif args.command == "selectRandom":
        selectRandom(args)
    elif args.command == "selectWorstRanked":
        selectWorstRanked(args)
    elif args.command == "selectBestRanked":
        selectBestRanked(args)
    elif args.command == "removeDuplicates":
        removeDuplicates(args)
    elif args.command == "assignClassName":
        assignClassName(args)
    elif args.command == "changeParamValue":
        changeParamValue(args)
    elif args.command == "removeParamValue":
        removeParamValue(args)
    elif args.command == "plotParamValue":
        plotParamValue(args)
    elif args.command == "eulerHist":
        eulerHist(args)
    elif args.command == "deleteTags":
        deleteTags(args)
    elif args.command == "csparc_setup":
        csparc_setup(args)
    else:
        janas_parser.print_help()


if __name__ == "__main__":
    main()
