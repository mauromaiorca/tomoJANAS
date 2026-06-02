# File: star_handler.py
# (C) 2025 Mauro Maiorca - Leibniz Institute of Virology


"""
star_handler.py

Utilities for reading, modifying, and writing Relion STAR files.
Provides functions to:
  - Inspect STAR file structure (infoStarFile, dataOptics)
  - Read and update specific sections (particles, optics)
  - Extract, add, remove, and replace STAR columns
  - Randomize half-sets and merge refinements

Author: janas
"""


from os import PathLike
import pandas as pd
import numpy as np

# from sklearn.tree import DecisionTreeClassifier
# from sklearn import metrics

# from matplotlib import pyplot
# from sklearn.tree import plot_tre
# from matplotlib.pyplot import figure
# import matplotlib.pyplot as plt


####################
# ACCESSORY FUNCTIONS
# START

# get starfile structure info
def infoStarFile(filename: PathLike):
    MAX_STAR_HEADER_SIZE = 500
    startLabels = 0
    startRawInfo = 0
    header = ""
    with open(filename) as f:
        # content = f.readlines()
        header = f.readlines()[:MAX_STAR_HEADER_SIZE]
    counter_loops = 0
    counter_data_optics = 0
    # search for the last loop_
    for ii in range(0, len(header)):
        tmpStr = header[ii].strip().replace(" ", "").replace("\t", "")
        if len(tmpStr) > 0:
            # print (tmpStr)
            if tmpStr.startswith("data_optics"):
                counter_data_optics += 1
            if tmpStr.startswith("loop_"):
                startLabels = ii + 1
                counter_loops += 1
            if tmpStr.startswith("_"):
                startRawInfo = ii + 2
    version = "relion_v30"
    if counter_loops > 1 and counter_data_optics > 0:
        version = "relion_v31"
    return startRawInfo, startLabels, version


# get dataOptics info
def dataOptics(filename: PathLike):
    version = infoStarFile(filename)[2]
    if version == "relion_v30":
        return

    MAX_STAR_HEADER_SIZE = 500
    with open(filename) as f:
        header = f.readlines()[:MAX_STAR_HEADER_SIZE]

    foundDataOptics = False
    foundDataOpticsLoop = False
    dataOpticsTags = []
    startRawIdx = 0
    endRawIdx = 0

    for ii in range(0, len(header)):
        tmpStr = header[ii].strip().replace(" ", "").replace("\t", "")
        if len(tmpStr) > 0:
            if tmpStr.startswith("data_optics"):
                foundDataOptics = True
            if (not foundDataOpticsLoop) and foundDataOptics and tmpStr.startswith("loop_"):
                foundDataOpticsLoop = True
            if foundDataOpticsLoop and foundDataOptics and tmpStr.startswith("_") and startRawIdx == 0:
                dataOpticsTags.append(header[ii].split(" ")[0].strip())
            if (
                foundDataOpticsLoop
                and foundDataOptics
                and not tmpStr.startswith("_")
                and startRawIdx == 0
                and not tmpStr.startswith("loop_")
            ):
                startRawIdx = ii

        if startRawIdx > 0 and endRawIdx == 0 and len(tmpStr) == 0:
            endRawIdx = ii

    rows_to_keep = [x for x in range(startRawIdx, endRawIdx)]
    optics_df = pd.read_csv(
        filename,
        skiprows=lambda x: x not in rows_to_keep,
        names=dataOpticsTags,
        skipinitialspace=True,
        sep=r"\s+",
    )

    # ---- Enrich optics_df with per-particle optics fields if missing ----
    required = [
        "_rlnImagePixelSize",
        "_rlnVoltage",
        "_rlnAmplitudeContrast",
        "_rlnSphericalAberration",
    ]
    missing = [c for c in required if c not in optics_df.columns]

    if missing:
        # Try to recover from data_particles
        desired = ["_rlnOpticsGroup"] + missing
        try:
            part = read_star_columns_from_sections(
                str(filename), section_name="particles", desired_columns=desired
            )
        except Exception:
            # If we cannot read particles section, return original optics_df
            return optics_df

        if part is None or part.empty or "_rlnOpticsGroup" not in part.columns:
            return optics_df

        # Convert numerics where possible
        # --- normalise merge key dtype: _rlnOpticsGroup must match on both sides ---
        # particles section is read as strings → coerce to numeric
        part["_rlnOpticsGroup"] = pd.to_numeric(part["_rlnOpticsGroup"], errors="coerce")
        optics_df["_rlnOpticsGroup"] = pd.to_numeric(optics_df["_rlnOpticsGroup"], errors="coerce")

        # If either side failed coercion, drop those rows (cannot be merged reliably)
        part = part.dropna(subset=["_rlnOpticsGroup"]).copy()
        optics_df = optics_df.dropna(subset=["_rlnOpticsGroup"]).copy()

        # Use a stable integer-like key (Relion optics groups are integers even if written as 1.000000)
        part["_rlnOpticsGroup"] = part["_rlnOpticsGroup"].round().astype("Int64")
        optics_df["_rlnOpticsGroup"] = optics_df["_rlnOpticsGroup"].round().astype("Int64")

        # Convert missing optics fields to numeric where possible
        for c in missing:
            if c in part.columns:
                part[c] = pd.to_numeric(part[c], errors="coerce")

        # One row per optics group (take first non-null)
        part_grp = (
            part.sort_values("_rlnOpticsGroup")
                .groupby("_rlnOpticsGroup", as_index=False)
                .first()
        )

        # Merge (prefer values already present in optics_df)
        optics_df = pd.merge(optics_df, part_grp, on="_rlnOpticsGroup", how="left")

    return optics_df



def read_star_columns_from_sections(
    filename: str, section_name: str, desired_columns: list
) -> pd.DataFrame:
    if not isinstance(
        desired_columns, list
    ):  # if desired_columns is a single string, convert it to a list
        desired_columns = [desired_columns]

    with open(filename, "r") as f:
        lines = f.readlines()

    in_section = False
    headers = []
    column_positions = []
    data = []

    for line in lines:
        stripped = line.strip()

        if stripped.startswith("data_"):
            in_section = stripped == f"data_{section_name}"
            continue

        if in_section:
            if stripped.startswith("loop_"):
                continue

            if stripped.startswith("_"):
                column_name = stripped.split()[0]
                headers.append(column_name)
                if column_name in desired_columns:
                    column_positions.append(headers.index(column_name))
            elif (
                stripped and column_positions
            ):  # if column_positions is not empty and line contains data
                row_data = stripped.split()
                desired_data = [row_data[i] for i in column_positions]
                data.append(desired_data)

    if not data:
        return None

    return pd.DataFrame(data, columns=[headers[i] for i in column_positions])


def read_star_sections(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    sections = {}
    current_section = None
    current_data = []

    for line in lines:
        stripped_line = line.strip()

        # If line is empty, skip processing this iteration
        if not stripped_line or stripped_line.startswith("#"):
            continue

        if stripped_line.startswith("data_"):
            if current_section:
                sections[current_section] = current_data
                current_data = []
            current_section = stripped_line.split("data_")[1]
        else:
            current_data.append(stripped_line)

    if current_section:
        sections[current_section] = current_data

    return sections


def process_section(section_lines, df):
    # Check if the section has "loop_"
    has_loop = any(line.startswith("loop_") for line in section_lines)

    # Extract header lines and data lines
    if has_loop:
        header_lines = [line for line in section_lines if line.startswith("_")]
        data_lines = [
            line
            for line in section_lines
            if not line.startswith("_") and not line.startswith("loop_")
        ]
    else:
        header_lines = [
            line.split()[0] for line in section_lines if line.startswith("_")
        ]
        data_lines = [
            line
            for line in section_lines
            if not line.startswith("_") and not line.startswith("loop_")
        ]

    # Extract header names without the hash numbers
    header_names = [line.split()[0] for line in header_lines]

    # Create a new dataframe from the section data
    section_df = pd.DataFrame(
        [line.split() for line in data_lines], columns=header_names
    )

    # Update the section dataframe with the passed dataframe
    for col in df.columns:
        if col in section_df.columns:
            section_df[col] = df[col].values
        else:
            section_df[col] = df[col].values

    # Generate the updated section lines
    updated_header_lines = [f"{col} #{i+1}" for i, col in enumerate(section_df.columns)]
    #    updated_data_lines = [" ".join(map(str, row)) for _, row in section_df.iterrows()]

    # Calculate maximum widths for each column
    max_widths = section_df.astype(str).apply(lambda col: col.map(len)).max()

    # Format data lines based on column widths
    updated_data_lines = [
        " ".join(str(row[col]).rjust(max_widths[col]) for col in section_df.columns)
        for _, row in section_df.iterrows()
    ]

    # Return the updated section
    if has_loop:
        return ["loop_"] + updated_header_lines + updated_data_lines
    else:
        return updated_header_lines + updated_data_lines


def update_star_columns_from_sections(filenameIn, filenameOut, section_name, df):
    # Load the STAR file into sections
    sections = read_star_sections(filenameIn)

    # Process the section named 'section_name' using the function we defined earlier
    updated_section = process_section(sections[section_name], df)

    # Replace the old section content with the updated one
    sections[section_name] = updated_section

    # Now, reconstruct the entire content
    new_content = []
    for key, value in sections.items():
        new_content.append("\ndata_" + key + "\n\n")
        for line in value:
            new_content.append(line + "\n")

    with open(filenameOut, "w") as f:
        f.writelines(new_content)


def randomize_halves(filenameIn, filenameOut):
    version = infoStarFile(filenameIn)[2]
    if version == "relion_v31":
        stringBlock = "particles"
    else:
        stringBlock = "images"
    outValues = read_star_columns_from_sections(
        filenameIn, stringBlock, ["_rlnRandomSubset"]
    )
    for col in ["_rlnRandomSubset"]:
        outValues[col] = np.random.choice([1, 2], size=len(outValues))
    update_star_columns_from_sections(filenameIn, filenameOut, stringBlock, outValues)


def delete_star_columns_from_sections(
    filenameIn, filenameOut, section_name, column_prefix
):
    # Load the STAR file into sections
    sections = read_star_sections(filenameIn)

    # Extract the section lines for the specified section name
    section_lines = sections[section_name]

    # Check if the section has "loop_"
    has_loop = any(line.startswith("loop_") for line in section_lines)

    # Extract header lines and data lines
    if has_loop:
        header_lines = [line for line in section_lines if line.startswith("_")]
        data_lines = [
            line
            for line in section_lines
            if not line.startswith("_") and not line.startswith("loop_")
        ]
    else:
        header_lines = [
            line.split()[0] for line in section_lines if line.startswith("_")
        ]
        data_lines = [
            line
            for line in section_lines
            if not line.startswith("_") and not line.startswith("loop_")
        ]

    # Extract header names without the hash numbers
    header_names = [line.split()[0] for line in header_lines]

    # Create a dataframe from the section data
    section_df = pd.DataFrame(
        [line.split() for line in data_lines], columns=header_names
    )

    # Drop columns with the given prefix
    columns_to_drop = [
        col for col in section_df.columns if col.startswith(column_prefix)
    ]
    section_df.drop(columns=columns_to_drop, inplace=True)

    # Generate the updated section lines
    updated_header_lines = [f"{col} #{i+1}" for i, col in enumerate(section_df.columns)]
    updated_data_lines = [" ".join(map(str, row)) for _, row in section_df.iterrows()]

    # Replace the old section content with the updated one
    sections[section_name] = (
        (["loop_"] + updated_header_lines + updated_data_lines)
        if has_loop
        else (updated_header_lines + updated_data_lines)
    )

    # Now, reconstruct the entire content
    new_content = []
    for key, value in sections.items():
        new_content.append("\ndata_" + key + "\n\n")
        for line in value:
            new_content.append(line + "\n")

    with open(filenameOut, "w") as f:
        f.writelines(new_content)


def replace_star_columns_from_sections(
    filenameIn, filenameOut, section_name, column_prefix, df
):
    # Load the STAR file into sections
    sections = read_star_sections(filenameIn)

    # Extract the section lines for the specified section name
    section_lines = sections[section_name]

    # Check if the section has "loop_"
    has_loop = any(line.startswith("loop_") for line in section_lines)

    # Extract header lines and data lines
    if has_loop:
        header_lines = [line for line in section_lines if line.startswith("_")]
        data_lines = [
            line
            for line in section_lines
            if not line.startswith("_") and not line.startswith("loop_")
        ]
    else:
        header_lines = [
            line.split()[0] for line in section_lines if line.startswith("_")
        ]
        data_lines = [
            line
            for line in section_lines
            if not line.startswith("_") and not line.startswith("loop_")
        ]

    # Extract header names without the hash numbers
    header_names = [line.split()[0] for line in header_lines]

    # Create a dataframe from the section data
    section_df = pd.DataFrame(
        [line.split() for line in data_lines], columns=header_names
    )

    # Replace values in columns with the given prefix
    for col in section_df.columns:
        if col.startswith(column_prefix):
            # Ensure the column exists in the provided DataFrame 'df'
            if col in df.columns:
                section_df[col] = df[col].values
            else:
                raise ValueError(f"Column '{col}' not found in the provided DataFrame")

    # Generate the updated section lines
    updated_header_lines = [f"{col} #{i+1}" for i, col in enumerate(section_df.columns)]
    updated_data_lines = [" ".join(map(str, row)) for _, row in section_df.iterrows()]

    # Replace the old section content with the updated one
    sections[section_name] = (
        (["loop_"] + updated_header_lines + updated_data_lines)
        if has_loop
        else (updated_header_lines + updated_data_lines)
    )

    # Now, reconstruct the entire content
    new_content = []
    for key, value in sections.items():
        new_content.append("\ndata_" + key + "\n\n")
        for line in value:
            new_content.append(line + "\n")

    with open(filenameOut, "w") as f:
        f.writelines(new_content)


def extract_particles_from_label_from_sections(
    filenameIn, filenameOut, section_name, column_prefix, column_prefix_value
):
    # Load the STAR file into sections
    sections = read_star_sections(filenameIn)

    # Extract the section lines for the specified section name
    section_lines = sections[section_name]

    # Check if the section has "loop_"
    has_loop = any(line.startswith("loop_") for line in section_lines)

    # Extract header lines and data lines
    if has_loop:
        header_lines = [line for line in section_lines if line.startswith("_")]
        data_lines = [
            line
            for line in section_lines
            if not line.startswith("_") and not line.startswith("loop_")
        ]
    else:
        header_lines = [
            line.split()[0] for line in section_lines if line.startswith("_")
        ]
        data_lines = [
            line
            for line in section_lines
            if not line.startswith("_") and not line.startswith("loop_")
        ]

    # Extract header names without the hash numbers
    header_names = [line.split()[0] for line in header_lines]

    # Create a dataframe from the section data
    section_df = pd.DataFrame(
        [line.split() for line in data_lines], columns=header_names
    )

    # Filter rows where the column with the prefix has a specific value
    if column_prefix in section_df.columns:
        section_df = section_df[section_df[column_prefix] == column_prefix_value]
    else:
        raise ValueError(f"Column '{column_prefix}' not found in the DataFrame")

    # Generate the updated section lines
    updated_header_lines = [f"{col} #{i+1}" for i, col in enumerate(section_df.columns)]
    updated_data_lines = [" ".join(map(str, row)) for _, row in section_df.iterrows()]

    # Replace the old section content with the updated one
    sections[section_name] = (
        (["loop_"] + updated_header_lines + updated_data_lines)
        if has_loop
        else (updated_header_lines + updated_data_lines)
    )

    # Now, reconstruct the entire content
    new_content = []
    for key, value in sections.items():
        new_content.append("\ndata_" + key + "\n\n")
        for line in value:
            new_content.append(line + "\n")

    with open(filenameOut, "w") as f:
        f.writelines(new_content)


def merge_star_section(
    inStarFile,
    nameMainSection="particles",
    nameIndexedSection="optics",
    indexingTag="_rlnOpticsGroup",
):
    # Read the STAR file sections
    sections = read_star_sections(inStarFile)

    # Extract data from main and indexed sections
    main_section_data = sections.get(nameMainSection, [])
    indexed_section_data = sections.get(nameIndexedSection, [])

    # Convert sections data to dataframes
    main_df_headers = [
        line.split()[0] for line in main_section_data if line.startswith("_")
    ]
    indexed_df_headers = [
        line.split()[0] for line in indexed_section_data if line.startswith("_")
    ]

    main_df_data = [
        line.split()
        for line in main_section_data
        if not line.startswith("_") and not line.startswith("loop_")
    ]
    indexed_df_data = [
        line.split()
        for line in indexed_section_data
        if not line.startswith("_") and not line.startswith("loop_")
    ]

    main_df = pd.DataFrame(main_df_data, columns=main_df_headers)
    indexed_df = pd.DataFrame(indexed_df_data, columns=indexed_df_headers)

    # Perform merge
    merged_df = pd.merge(main_df, indexed_df, on=[indexingTag])

    # Drop '_rlnOpticsGroup'
    merged_df = merged_df.drop([indexingTag], axis=1)

    return merged_df


# END
# ACCESSORY FUNCTIONS
#######################


###########################
# get columns in header file
def header_columns(filename: PathLike):
    """get columns in star file

    Arguments:
        Inputs:
            filename: PathLike
                star file name

        Outputs:
            outVect: float, 1D array
                Unit cell
 """
    startRawInfo, startLabels, version = infoStarFile(filename)
    with open(filename) as f:
        header = f.readlines()[: startRawInfo - 1]
    outVect = []
    for ii in range(startLabels, len(header)):
        tmpStr = header[ii].lstrip().split()
        if len(tmpStr[0]) > 0:
            outVect.append(tmpStr[0])
    return outVect


###########################
# read one or multiple columns in star file
def readColumns(filename: PathLike, columnsToRead, sortColumnsNameOrder=False):
    # print ("readColumns ", filename)
    header_list = header_columns(filename)
    startRawInfo = infoStarFile(filename)[0]
    df = pd.read_csv(
        filename,
        skiprows=startRawInfo - 1,
        names=header_list,
        usecols=columnsToRead,
        skipinitialspace=True,
        sep="\s+",
    )
    if sortColumnsNameOrder:  # output of the same order of columnsToRead
        outDf = df[columnsToRead]
    else:
        outDf = df
    return outDf


####################################
# read all columns of the star file
def readStar(filename: PathLike):
    # print ("readColumns ", filename)
    header_list = header_columns(filename)
    startRawInfo = infoStarFile(filename)[0]
    df = pd.read_csv(
        filename,
        skiprows=startRawInfo - 1,
        names=header_list,
        skipinitialspace=True,
        sep="\s+",
    )
    return df


def extractColumns(filenameIn: PathLike, filenameOut: PathLike, listColumnsnameToRead):
    # print ("readColumns ", filenameIn)
    startRawInfo, startLabels, version = infoStarFile(filenameIn)[0]
    pre_header = ""
    if startLabels > 0:
        with open(filenameIn) as f:
            pre_header = f.readlines()[:startLabels]
    pre_header = "".join(pre_header)
    columns = readColumns(filenameIn, listColumnsnameToRead, sortColumnsNameOrder=True)
    # print (columns.head)
    # print(columns.keys())
    # print(df.columns.tolist())
    for ii in range(0, len(listColumnsnameToRead)):
        pre_header += str(listColumnsnameToRead[ii]) + " #" + str(ii + 1) + "\n"
    with open(filenameOut, "w") as fw:
        fw.write(pre_header)
        columns.to_csv(fw, header=False, sep=" ", index=False)
        fw.close()


###############################
# remove selected column (drop function)
def removeColumns(
    filenameIn: PathLike, filenameOut: PathLike, listColumnsToRemoveNames
):
    # print ("removeColumns")
    startRawInfo, startLabels, version = infoStarFile(filenameIn)
    columns = header_columns(filenameIn)
    # get pre-header information
    pre_header = ""
    if startLabels > 0:
        with open(filenameIn) as f:
            pre_header = f.readlines()[:startLabels]
    pre_header = "".join(pre_header)
    df = pd.read_csv(
        filenameIn,
        skiprows=startRawInfo - 1,
        names=columns,
        skipinitialspace=True,
        sep="\s+",
    )
    columnsDeletedNames = []
    mask = np.isin(columns, listColumnsToRemoveNames, invert=True)
    for ii in range(0, len(columns)):
        if mask[ii]:
            columnsDeletedNames.append(columns[ii])
    for ii in range(0, len(columnsDeletedNames)):
        pre_header += str(columnsDeletedNames[ii]) + " #" + str(ii + 1) + "\n"
    with open(filenameOut, "w") as fw:
        fw.write(pre_header)
        X = df.drop(listColumnsToRemoveNames, axis="columns", inplace=False)
        X.to_csv(fw, header=False, sep=" ", index=False)
        fw.close()


###############################################################
# remove Columns with a Tags Starting With "StartingWithTags"
def removeColumnsTagsStartingWith(
    filenameIn: PathLike, filenameOut: PathLike, StartingWithTags
):
    # print ("removeColumnsTagsStartingWith")
    # startRawInfo,startLabels=infoStarFile(filenameIn)
    columns = header_columns(filenameIn)
    columnsToRemove = []
    for iiObj in columns:
        if iiObj.startswith(StartingWithTags):
            columnsToRemove.append(iiObj)
    removeColumns(filenameIn, filenameOut, columnsToRemove)


#########################
# add columns
def addDataframeColumns(
    filenameIn: PathLike, filenameOut: PathLike, columnsToAddName, columnToAddContent
):
    # print ("addColumns")
    startRawInfo, startLabels, version = infoStarFile(filenameIn)
    columns = header_columns(filenameIn)
    # get pre-header information
    pre_header = ""
    if startLabels > 0:
        with open(filenameIn) as f:
            pre_header = f.readlines()[:startLabels]
    pre_header = "".join(pre_header)
    df = pd.read_csv(
        filenameIn,
        skiprows=startRawInfo - 1,
        names=columns,
        skipinitialspace=True,
        sep="\s+",
    )
    columnsAdded = columns
    for ii in range(0, len(columnsAdded)):
        pre_header += str(columnsAdded[ii]) + " #" + str(ii + 1) + "\n"
    counter = 0
    for kk in range(0, len(columnsToAddName)):
        if str(columnsToAddName[kk]) not in columns:
            pre_header += (
                str(columnsToAddName[kk])
                + " #"
                + str(len(columnsAdded) + 1 + counter)
                + "\n"
            )
            counter += 1

    # df=df.drop(columnsToAddName, axis='columns', inplace=True)
    df[columnsToAddName] = columnToAddContent.values
    # fullColumns=np.append(columns, columnToAddContent, axis=1)
    with open(filenameOut, "w") as fw:
        fw.write(pre_header)
        # df=df.drop(columns, axis='columns', inplace=False)
        df.to_csv(fw, header=False, sep=" ", index=False)
        fw.close()


#########################
# add columns
def addColumns(
    filenameIn: PathLike, filenameOut: PathLike, columnToAddName, columnToAddContent
):
    # print ("addColumns")
    startRawInfo, startLabels, version = infoStarFile(filenameIn)
    columns = header_columns(filenameIn)
    # get pre-header information
    pre_header = ""
    if startLabels > 0:
        with open(filenameIn) as f:
            pre_header = f.readlines()[:startLabels]
    pre_header = "".join(pre_header)
    df = pd.read_csv(
        filenameIn,
        skiprows=startRawInfo - 1,
        names=columns,
        skipinitialspace=True,
        sep="\s+",
    )
    columnsAdded = columns
    for ii in range(0, len(columnsAdded)):
        pre_header += str(columnsAdded[ii]) + " #" + str(ii + 1) + "\n"
    for kk in range(0, len(columnToAddName)):
        pre_header += (
            str(columnToAddName[kk])
            + " #"
            + str(len(columnsAdded) + len(columnToAddName) + kk)
            + "\n"
        )
        df[columnToAddName[kk]] = columnToAddContent[kk]
    with open(filenameOut, "w") as fw:
        fw.write(pre_header)
        df.to_csv(fw, header=False, sep=" ", index=False)
        fw.close()


#################################################
# writeDataframe To Star using a template starfile
def writeDataframeToStar(filenameIn: PathLike, filenameOut: PathLike, dataframe):
    startRawInfo, startLabels, version = infoStarFile(filenameIn)
    # print (df)
    header = ""
    if startLabels > 0:
        with open(filenameIn) as f:
            header = f.readlines()[: startRawInfo - 1]
    header = "".join(header)
    with open(filenameOut, "w") as fw:
        fw.write(header)
        dataframe.to_csv(fw, header=False, sep=" ", index=False)
        fw.close()


#################################################
# writeDataframe To Star without a template starfile
# uses information stored in dataframe
def writeDataframeToStar_deNovo(dataframe, filenameOut: PathLike):
    header = "data_images\n\nloop_ \n"
    for ii in range(0, len(dataframe.columns)):
        header += str(dataframe.columns[ii]) + " #" + str(ii + 1) + "\n"
    with open(filenameOut, "w") as fw:
        fw.write(header)
        dataframe.to_csv(fw, header=False, sep=" ", index=False)
        fw.close()


def extractWorst(filenameIn: PathLike, filenameOut: PathLike, numItems, tagToSelect):
    # print ("selectWorst ",numItems)
    header_list = header_columns(filenameIn)
    startRawInfo = infoStarFile(filenameIn)[0]
    df = pd.read_csv(
        filenameIn,
        skiprows=startRawInfo - 1,
        names=header_list,
        skipinitialspace=True,
        sep="\s+",
    )

    # print("Data before sorting:")
    # print(df.head())

    df = df.sort_values([tagToSelect], ascending=False, kind="quicksort").head(numItems)

    # print("Data after sorting:")
    # print(df.head())  # Print the first few rows after sorting for inspection

    df = df.sort_index()
    writeDataframeToStar(filenameIn, filenameOut, df)


# extract best/worst
def extractBest(filenameIn: PathLike, filenameOut: PathLike, numItems, tagToSelect, exact: bool = False):
    # print ("selectBest ",numItems)
    header_list = header_columns(filenameIn)
    startRawInfo = infoStarFile(filenameIn)[0]
    df = pd.read_csv(
        filenameIn,
        skiprows=startRawInfo - 1,
        names=header_list,
        skipinitialspace=True,
        sep="\s+",
    )
    available = len(df)
    if exact and int(numItems) > available:
        raise ValueError(
            "selectBestRanked --exact: requested {} particles but only {} are available in {}".format(
                int(numItems), available, filenameIn
            )
        )
    df = df.sort_values([tagToSelect], ascending=True, kind="quicksort").head(numItems)
    df = df.sort_index()
    writeDataframeToStar(filenameIn, filenameOut, df)


# extract Random
def extractRandom(filenameIn: PathLike, filenameOut: PathLike, numItems):
    # print ("selectRandom ",numItems)
    header_list = header_columns(filenameIn)
    startRawInfo = infoStarFile(filenameIn)[0]
    df = pd.read_csv(
        filenameIn,
        skiprows=startRawInfo - 1,
        names=header_list,
        skipinitialspace=True,
        sep="\s+",
    )
    # df = df.sort_values([tagToSelect],ascending=True, kind='quicksort').head(numItems)
    df = df.iloc[np.random.permutation(len(df))].head(numItems)
    df = df.sort_index()
    writeDataframeToStar(filenameIn, filenameOut, df)


def removeStarDuplicates(filenameIn: PathLike, filenameOut: PathLike):
    header_list = header_columns(filenameIn)
    # print(header_list)
    startRawInfo = infoStarFile(filenameIn)[0]
    df = pd.read_csv(
        filenameIn,
        skiprows=startRawInfo - 1,
        names=header_list,
        skipinitialspace=True,
        sep="\s+",
    )
    df = df.drop_duplicates(subset="_rlnImageName", keep="first")
    writeDataframeToStar(filenameIn, filenameOut, df)


def AssignStarClassName(filenameIn: PathLike, filenameOut: PathLike, className):
    header_list = header_columns(filenameIn)
    startRawInfo = infoStarFile(filenameIn)[0]
    df = pd.read_csv(
        filenameIn,
        skiprows=startRawInfo - 1,
        names=header_list,
        skipinitialspace=True,
        sep="\s+",
    )
    df["_rlnClassNumber"] = className
    writeDataframeToStar(filenameIn, filenameOut, df)


def changeStarFileParamValue(
    filenameIn: PathLike, filenameOut: PathLike, columnName, newParamValue
):
    header_list = header_columns(filenameIn)
    startRawInfo = infoStarFile(filenameIn)[0]
    df = pd.read_csv(
        filenameIn,
        skiprows=startRawInfo - 1,
        names=header_list,
        skipinitialspace=True,
        sep="\s+",
    )
    df[columnName] = newParamValue
    writeDataframeToStar(filenameIn, filenameOut, df)


def deleteStarFileParamValue0(
    filenameIn: PathLike, filenameOut: PathLike, columnName, valueToDelete
):
    header_list = header_columns(filenameIn)
    startRawInfo = infoStarFile(filenameIn)[0]
    df = pd.read_csv(
        filenameIn,
        skiprows=startRawInfo - 1,
        names=header_list,
        skipinitialspace=True,
        sep="\s+",
    )
    df = df[df[columnName] != valueToDelete]
    writeDataframeToStar(filenameIn, filenameOut, df)


def deleteStarFileParamValue(
    filenameIn: PathLike, filenameOut: PathLike, columnName, valueToDelete
):
    header_list = header_columns(filenameIn)
    startRawInfo = infoStarFile(filenameIn)[0]
    df = pd.read_csv(
        filenameIn,
        skiprows=startRawInfo - 1,
        names=header_list,
        skipinitialspace=True,
        sep="\s+",
    )
    if isinstance(valueToDelete, str):
        try:
            valueToDelete = int(valueToDelete)
        except ValueError:
            print("Error: valueToDelete cannot be converted to an integer.")
            return
    print("value to remove=", valueToDelete)
    df_filtered = df[df[columnName] == valueToDelete]
    # Debug: Check the effect of filtering
    print("Rows before filtering:", len(df))
    print("Rows after filtering:", len(df_filtered))
    writeDataframeToStar(filenameIn, filenameOut, df_filtered)


def plotStarFileParamValue(filenameIn: PathLike, columnName):
    import matplotlib.pyplot as plt

    header_list = header_columns(filenameIn)
    startRawInfo = infoStarFile(filenameIn)[0]
    df = pd.read_csv(
        filenameIn,
        skiprows=startRawInfo - 1,
        names=header_list,
        skipinitialspace=True,
        sep="\s+",
    )
    sorted_values = df[columnName].sort_values()
    plt.scatter(range(len(sorted_values)), sorted_values, linestyle="-")
    plt.yscale("log")
    plt.xlabel(columnName)
    plt.ylabel("Values")
    plt.title(f"Plot of {columnName}")
    plt.grid(True)
    plt.show()


def plotStarFileParamValueInteractive(filenameIn: PathLike, columnName):
    import matplotlib.pyplot as plt

    header_list = header_columns(filenameIn)
    startRawInfo = infoStarFile(filenameIn)[0]
    df = pd.read_csv(
        filenameIn,
        skiprows=startRawInfo - 1,
        names=header_list,
        skipinitialspace=True,
        sep="\s+",
    )

    sorted_values = df[columnName].sort_values().reset_index(drop=True)

    fig, ax = plt.subplots()
    (line,) = ax.plot(
        sorted_values, linestyle="-", marker="o", pickradius=5
    )  # `pickradius` determines the sensitivity of the click

    def onpick(event):
        # When the line is clicked, print the coordinates
        ind = event.ind[0]
        x, y = event.artist.get_data()
        print(f"Value at index {ind}: ({x[ind]}, {y[ind]})")

    # Connect the click handler to the figure
    fig.canvas.mpl_connect("pick_event", onpick)

    plt.xlabel(columnName)
    plt.ylabel("Values")
    plt.title(f"Plot of {columnName}")
    plt.grid(True)
    plt.show()


def extractCategory(
    filenameIn: PathLike, filenameOut: PathLike, categoryName, categoryValue
):
    # print ("extractCategory")
    header_list = header_columns(filenameIn)
    startRawInfo = infoStarFile(filenameIn)[0]
    df = pd.read_csv(
        filenameIn,
        skiprows=startRawInfo - 1,
        names=header_list,
        skipinitialspace=True,
        sep="\s+",
    )
    df = df.loc[(df[categoryName].astype(str) == str(categoryValue))]
    # df = df.sort_values(['_LRA_CC_unprocessed_simple'],ascending=False, kind='quicksort').head(numItems)
    # df = df.sort_index()
    # _LRA_particle_category
    writeDataframeToStar(filenameIn, filenameOut, df)
    # removeColumns(filenameOut, filenameOut, [categoryName])


def extractImageNameInfo_starFile(InputStarFile: PathLike):
    imageNameTag = "_rlnImageName"
    referenceColumns = readColumns(InputStarFile, [imageNameTag])
    listOfDatasets = []
    listOccurrencesNames = []
    listImagesIdx = []
    listOccurrencesNamesIdx = []
    for ii in range(0, len(referenceColumns[imageNameTag])):
        tmpLine = referenceColumns[imageNameTag][ii]
        atPosition = tmpLine.find("@")
        imageNo = int(tmpLine[:atPosition])
        stackName = tmpLine[atPosition + 1 :]
        listOccurrencesNames.append(stackName)
        listImagesIdx.append(imageNo)
        if not (stackName in listOfDatasets):
            listOfDatasets.append(stackName)
        listOccurrencesNamesIdx.append(listOfDatasets.index(stackName))
    return listImagesIdx, listOccurrencesNamesIdx, listOfDatasets


################################################
#  SIMPLEST WAY TO MERGE MULTIPLE REFINEMENTS
#  mergeRefinementsWithHighestScore
def mergeRefinements(
    referenceRefinementStarFile: PathLike,
    fileNameOut: PathLike,
    listRefinementFiles,
    primaryScore="_janas_SCI__1",
    refinedTags=[
        "_rlnOriginXAngst",
        "_rlnOriginYAngst",
        "_rlnAngleRot",
        "_rlnAngleTilt",
        "_rlnAnglePsi",
    ],
):
    df = readColumns(referenceRefinementStarFile, [primaryScore])
    (
        listImagesIdx,
        listOccurrencesNamesIdx,
        listOfDatasets,
    ) = extractImageNameInfo_starFile(referenceRefinementStarFile)

    listDataFilenameIdx = [0] * len(listImagesIdx)
    listDataIdx = [-1] * len(listImagesIdx)
    for kk in range(0, len(listRefinementFiles)):
        pos1 = [-1] * len(listImagesIdx)
        pos1_scores = [-10] * len(listImagesIdx)
        (
            listImagesIdx1,
            listOccurrencesNamesIdx1,
            listOfDatasets1,
        ) = extractImageNameInfo_starFile(listRefinementFiles[kk])
        referenceScores_tmp = readColumns(listRefinementFiles[kk], [primaryScore])
        # print ( kk, ' ->',referenceScore_tmp )
        for ii in range(0, len(listImagesIdx)):
            for jj in range(0, len(listImagesIdx1)):
                if (
                    listImagesIdx1[jj] == listImagesIdx[ii]
                    and listOccurrencesNamesIdx1[jj] == listOccurrencesNamesIdx[ii]
                ):
                    pos1[ii] = jj
                    pos1_scores[ii] = referenceScores_tmp.iloc[jj][primaryScore]
                    # print ('ii=',ii,'   jj=',jj, '  => ', referenceScore_tmp.iloc[jj]['_LRA_CC_unprocessed_simple'])
                    break
        df["idx_ParameterLine_" + str(kk)] = pos1
        df["scores_file_" + str(kk)] = pos1_scores

    for ii in range(0, len(listImagesIdx)):
        # print (ii, '  ==> ',df.iloc[ii])
        listTargetParameterFile = [referenceRefinementStarFile]
        listTargetParameterLine = [ii]
        listScores = [df.iloc[ii][primaryScore]]
        # OK   meanTargetEulerId=[df.iloc[ii]['_LRA_CC_unprocessed_simple']]
        # OK   varianceTargetEulerId=[df.iloc[ii]['_LRA_CC_unprocessed_simple']]
        # print ('values =')
        for kk in range(0, len(listRefinementFiles)):
            tmp_idx = df.iloc[ii]["idx_ParameterLine_" + str(kk)]
            if tmp_idx >= 0:  # other files have other values
                listScores.append(df.iloc[ii]["scores_file_" + str(kk)])
                listTargetParameterFile.append("scores_file_" + str(kk))
                listTargetParameterLine.append(tmp_idx)
        sortedIdx = np.array(listScores).argsort()[::-1]
        # sorted_scores = np.array(listScores)[sortedIdx]
        sorted_filenames = np.array(listTargetParameterFile)[sortedIdx]
        sorted_paramsIdx = np.array(listTargetParameterLine)[sortedIdx]

        #####
        # criterium for selecting the particle:
        #   thresholdNumParticles=10
        listDataFilenameIdx[ii] = sorted_filenames[0]
        listDataIdx[ii] = sorted_paramsIdx[0]
        # targetBestParticle=0
        # targetZscore=0
        # for ii in range(0, len(sorted_filenames)):

        ######
        # zscoretest
    #   print (ii, '  ==>  eulerID=',eulerId[ii],'    numParticlesPerEuler=',len(sorted_filenames),'    ParticlesEulerId=',countEulerID[eulerIdx],'   density%=',100.0*countEulerID[eulerIdx]/len(listImagesIdx),'   MeanScoreEulerID=',meanScoreEulerID[eulerIdx],'   stdScoreEulerID=',stdScoreEulerID[eulerIdx] )
    #   for tt in range(0,len(sorted_filenames)):
    #    zscore=0
    #    if stdScoreEulerID[eulerIdx]>0:
    #      zscore= (sorted_scores[tt]-meanScoreEulerID[eulerIdx])/stdScoreEulerID[eulerIdx]
    #    print ('    score=',  sorted_scores[tt],'   zscore=',zscore)

    ######

    # print (ii, '  ==> ',listScores, '   list_files=', listTargetParameterFile,'   list_scores=', listScores,'   targetFile=', listDataFilenameIdx[ii],'   targetIdx=', listDataIdx[ii] )

    referenceTags = readColumns(referenceRefinementStarFile, refinedTags)
    referenceTagsOut = referenceTags.copy()
    referenceTags["listDataIdx"] = pd.DataFrame(
        listDataIdx, columns=["listDataIdx"], dtype="int"
    )
    referenceTags["listDataFilenameIdx"] = pd.DataFrame(
        listDataFilenameIdx, columns=["listDataFilenameIdx"]
    )

    # save the file with all the amendments
    for kk in range(0, len(listRefinementFiles)):
        tmpTags = readColumns(listRefinementFiles[kk], refinedTags)
        # listDataFilenameIdx
        # tmp0=(referenceTags.loc[ (referenceTags['listDataFilenameIdx'] == 'scores_file_'+str(kk)) ]['listDataIdx'])
        idx = referenceTags.loc[
            (referenceTags["listDataFilenameIdx"] == "scores_file_" + str(kk))
        ]["listDataIdx"]
        tmpTags = tmpTags.iloc[idx.to_numpy()]
        tmpTags.index = idx.index
        referenceTagsOut.loc[tmpTags.index, :] = tmpTags[:]
        # print ('\n+++++++\n idx=', idx, '\n' ,tmpTags.head)
    del referenceTags
    # writeDataframeToStar_deNovo(referenceTagsOut, 'test.star')
    # print (referenceTagsOut.head)
    fullStar = readStar(referenceRefinementStarFile)
    for obj in refinedTags:
        fullStar[obj] = referenceTagsOut[obj]
    writeDataframeToStar(referenceRefinementStarFile, fileNameOut, fullStar)


################################################
################################################
#  MULTIVARIATE ANALYSIS SVD
#  mergeRefinements
def mergeRefinements_SVD(
    referenceRefinementStarFile: PathLike,
    fileNameOut: PathLike,
    listRefinementFiles,
    referenceScoresColumnsName,
    primaryScore="_LRA_CC_unprocessed_simple",
    minNumOfParticlesPerBin=3,
    numBins=10,
    refinedTags=[
        "_rlnOriginXAngst",
        "_rlnOriginYAngst",
        "_rlnAngleRot",
        "_rlnAngleTilt",
        "_rlnAnglePsi",
    ],
):
    # print ("refinementsMerge")

    df = readColumns(referenceRefinementStarFile, referenceScoresColumnsName)
    (
        listImagesIdx,
        listOccurrencesNamesIdx,
        listOfDatasets,
    ) = extractImageNameInfo_starFile(referenceRefinementStarFile)

    # get the euler group
    eulerId = [-1] * len(listImagesIdx)
    meanScoreEulerID = [0] * (numBins * numBins)
    stdScoreEulerID = [0] * (numBins * numBins)
    countEulerID = [0] * (numBins * numBins)
    Phi = readColumns(referenceRefinementStarFile, ["_rlnAngleRot"])
    Theta = readColumns(referenceRefinementStarFile, ["_rlnAngleTilt"])
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
    # H = H.T
    # print(H[:])
    for ii in range(0, len(listImagesIdx)):
        digX = np.digitize(eulers1[ii], xedges) - 1
        digY = np.digitize(eulers2[ii], yedges) - 1
        eulerId[ii] = digX + numBins * digY
    df["euler_ID"] = eulerId
    # get statistics for each group
    classifierPerGroup = []  # *(numBins*numBins)
    for ii in range(0, numBins * numBins):
        counts = (df["euler_ID"] == ii).sum()
        countEulerID[ii] = counts
        # print ('\n*****************************\n full=', df.loc[df['euler_ID'] == ii][referenceScoresColumnsName])
        # NOTE: trick to select the best particles from view, and store the scores into a numpy vector
        ID_scores = (
            df.loc[df["euler_ID"] == ii][referenceScoresColumnsName]
            .sort_values(primaryScore, ascending=False)
            .head(minNumOfParticlesPerBin)
            .to_numpy()
        )
        # print ('partial=', ID_scores)
        # print ('training data shape=',np.shape(ID_scores))
        #   print ('numParameters=',len(referenceScoresColumnsName))
        if counts < minNumOfParticlesPerBin:
            trainedParameters = [[1] * len(referenceScoresColumnsName)]
            U1, D1, V1 = np.linalg.svd(
                np.transpose(trainedParameters, (1, 0)), full_matrices=False
            )
            # print ('U1=',U1)
            classifierPerGroup.append(zip(U1))
        else:
            U1, D1, V1 = np.linalg.svd(
                np.transpose(ID_scores, (1, 0)), full_matrices=False
            )
            classifierPerGroup.append(zip(U1))
    print("\n\n******************\nclassifierPerGroup=", len(classifierPerGroup))

    # mean and std for each view
    # if counts > 0:
    # ID_scores=df.loc[df['euler_ID'] == ii]['_LRA_CC_unprocessed_simple'].to_numpy()
    # ID_scores=np.sort( ID_scores)[::-1]
    # meanScoreEulerID[ii]=np.average(ID_scores)
    # stdScoreEulerID[ii]=np.std(ID_scores)
    # ID_scores=np.sort( ID_scores)[::-1]
    # print ('\n\n---------\n',ii,' => ', counts, '  avg=', meanScoreEulerID[ii], '  std=', stdScoreEulerID[ii] )

    # train classifier for each group

    # plt.imshow(H)
    # plt.show()

    # inspect the files for refinement
    # and compute for each particle of each file the reference score
    listDataFilenameIdx = [0] * len(listImagesIdx)
    listDataIdx = [-1] * len(listImagesIdx)

    for kk in range(0, len(listRefinementFiles)):
        pos1 = [-1] * len(listImagesIdx)
        pos1_scores = [-10] * len(listImagesIdx)
        (
            listImagesIdx1,
            listOccurrencesNamesIdx1,
            listOfDatasets1,
        ) = extractImageNameInfo_starFile(listRefinementFiles[kk])
        referenceScores_tmp = readColumns(
            listRefinementFiles[kk], referenceScoresColumnsName
        )
        # print ( kk, ' ->',referenceScore_tmp )
        print(" ite ", kk, " of ", len(listRefinementFiles))
        for ii in range(0, len(listImagesIdx)):
            for jj in range(0, len(listImagesIdx1)):
                if (
                    listImagesIdx1[jj] == listImagesIdx[ii]
                    and listOccurrencesNamesIdx1[jj] == listOccurrencesNamesIdx[ii]
                ):
                    pos1[ii] = jj
                    y = (
                        (referenceScores_tmp.iloc[jj][referenceScoresColumnsName])
                        .to_numpy()
                        .tolist()
                    )
                    # print ('ii=',ii,'   multivariateScores=',y,'    data shape=',np.shape(y))
                    # unpack U1
                    unzippedMatrix = []
                    tmpU1 = list(classifierPerGroup[eulerId[ii]])
                    # print ('ii=',ii,'     eulerID=',eulerId[ii],'     U1=',tmpU1)
                    for mm in range(0, len(tmpU1)):
                        unzippedMatrix.append(np.array(tmpU1[mm][0]).tolist())
                    U1 = np.array(unzippedMatrix)
                    # print ('     unzipped  U1=',U1)
                    # print ('     type=',type(unzippedMatrix))
                    y1 = np.matmul(U1, np.matmul(np.transpose(U1), y))
                    y1_distance = np.linalg.norm(y - y1)
                    # print ('     distance=',y1_distance)
                    pos1_scores[ii] = y1_distance
                    classifierPerGroup[eulerId[ii]] = zip(U1)

                    # pos1_scores[ii]=0
                    # y=testParameters[ii]
                    # y1=np.matmul(U1,np.matmul(np.transpose(U1),y))
                    # y1_distance=np.linalg.norm(y-y1)
                    # pos1_scores[ii]=
                    # print ('ii=',ii,'   jj=',jj, '  => ', referenceScore_tmp.iloc[jj]['_LRA_CC_unprocessed_simple'])
                    break
        df["idx_ParameterLine_" + str(kk)] = pos1
        df["scores_file_" + str(kk)] = pos1_scores
    # print (df['idx_file_0'])
    # select best particles out of options

    for ii in range(0, len(listImagesIdx)):
        # print (ii, '  ==> ',df.iloc[ii])
        listTargetParameterFile = [referenceRefinementStarFile]
        listTargetParameterLine = [ii]
        # compute the score for the original files
        unzippedMatrix = []
        tmpU1 = list(classifierPerGroup[eulerId[ii]])
        for mm in range(0, len(tmpU1)):
            unzippedMatrix.append(np.array(tmpU1[mm][0]).tolist())
        U1 = np.array(unzippedMatrix)
        y1 = np.matmul(U1, np.matmul(np.transpose(U1), y))
        y1_distance = np.linalg.norm(y - y1)
        listScores = [y1_distance]
        classifierPerGroup[eulerId[ii]] = zip(U1)

        # OK   listScores=[df.iloc[ii]['_LRA_CC_unprocessed_simple']]
        # OK   listScores=[df.iloc[ii][referenceScoresColumnsName]]
        # OK   meanTargetEulerId=[df.iloc[ii]['_LRA_CC_unprocessed_simple']]
        # OK   varianceTargetEulerId=[df.iloc[ii]['_LRA_CC_unprocessed_simple']]
        # print ('values =')
        for kk in range(0, len(listRefinementFiles)):
            tmp_idx = df.iloc[ii]["idx_ParameterLine_" + str(kk)]
            if tmp_idx >= 0:  # other files have other values
                listScores.append(df.iloc[ii]["scores_file_" + str(kk)])
                listTargetParameterFile.append("scores_file_" + str(kk))
                listTargetParameterLine.append(tmp_idx)

        ##########################################################################
        # NOTE: do not need to sort the score, but to compute the L2 distance using statistic trick
        sortedIdx = np.array(listScores).argsort()[::-1]
        sorted_scores = np.array(listScores)[sortedIdx]
        sorted_filenames = np.array(listTargetParameterFile)[sortedIdx]
        sorted_paramsIdx = np.array(listTargetParameterLine)[sortedIdx]

        #####
        # criterium for selecting the particle:
        #   thresholdNumParticles=10
        eulerIdx = eulerId[ii]
        listDataFilenameIdx[ii] = sorted_filenames[0]
        listDataIdx[ii] = sorted_paramsIdx[0]
        # targetBestParticle=0
        # targetZscore=0
        # for ii in range(0, len(sorted_filenames)):

        ######
        # zscoretest
    #   print (ii, '  ==>  eulerID=',eulerId[ii],'    numParticlesPerEuler=',len(sorted_filenames),'    ParticlesEulerId=',countEulerID[eulerIdx],'   density%=',100.0*countEulerID[eulerIdx]/len(listImagesIdx),'   MeanScoreEulerID=',meanScoreEulerID[eulerIdx],'   stdScoreEulerID=',stdScoreEulerID[eulerIdx] )
    #   for tt in range(0,len(sorted_filenames)):
    #    zscore=0
    #    if stdScoreEulerID[eulerIdx]>0:
    #      zscore= (sorted_scores[tt]-meanScoreEulerID[eulerIdx])/stdScoreEulerID[eulerIdx]
    #    print ('    score=',  sorted_scores[tt],'   zscore=',zscore)

    ######

    # print (ii, '  ==> ',listScores, '   list_files=', listTargetParameterFile,'   list_scores=', listScores,'   targetFile=', listDataFilenameIdx[ii],'   targetIdx=', listDataIdx[ii] )

    referenceTags = readColumns(referenceRefinementStarFile, refinedTags)
    referenceTagsOut = referenceTags.copy()
    referenceTags["listDataIdx"] = pd.DataFrame(
        listDataIdx, columns=["listDataIdx"], dtype="int"
    )
    referenceTags["listDataFilenameIdx"] = pd.DataFrame(
        listDataFilenameIdx, columns=["listDataFilenameIdx"]
    )

    # save the file with all the amendments
    for kk in range(0, len(listRefinementFiles)):
        tmpTags = readColumns(listRefinementFiles[kk], refinedTags)
        # listDataFilenameIdx
        # tmp0=(referenceTags.loc[ (referenceTags['listDataFilenameIdx'] == 'scores_file_'+str(kk)) ]['listDataIdx'])
        idx = referenceTags.loc[
            (referenceTags["listDataFilenameIdx"] == "scores_file_" + str(kk))
        ]["listDataIdx"]
        tmpTags = tmpTags.iloc[idx.to_numpy()]
        tmpTags.index = idx.index
        referenceTagsOut.loc[tmpTags.index, :] = tmpTags[:]
        # print ('\n+++++++\n idx=', idx, '\n' ,tmpTags.head)
    del referenceTags
    # writeDataframeToStar_deNovo(referenceTagsOut, 'test.star')
    # print (referenceTagsOut.head)
    fullStar = readStar(referenceRefinementStarFile)
    for obj in refinedTags:
        fullStar[obj] = referenceTagsOut[obj]
    writeDataframeToStar(referenceRefinementStarFile, fileNameOut, fullStar)
