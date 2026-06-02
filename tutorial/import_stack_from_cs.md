# Import stack from cryoSPARC



JANAS operates on particle stacks referenced in RELION-style STAR files. 

To begin, locate the cryoSPARC job directory that produced the particle metadata you want to export. This will be a path such as `CS-directory/P1/J80`, where `P1` is the project folder and `J80` is a job (for example, a Non-uniform Refinement). Inside that job directory you will find a file named like `cryosparc_P1_J80_002_particles.cs`, which contains the particle metadata and references to the corresponding stacks.

```
janas_utils csparc2star-stack /your/absolute/Path/Directory/cryosparc/P1/J80/cryosparc_P1_J80_002_particles.cs outStack
```

This will generate:

- `outStack.star` — a STAR file containing the particle metadata, with _rlnImageName entries pointing to the new stack. An additional column _janas_csparc_rlnImageName records the original cryoSPARC _rlnImageName values for provenance.

- `outStack.mrcs` — a 2D particle stack assembled one particle at a time from the original cryoSPARC stack locations.

The resulting files can be used directly in JANAS or imported into other software (e.g. RELION).

(note: you can get the absolute path of your directory by using the bash command `pwd`)

# Import stack from cryoSPARC (legacy)

Before `csparc2star-stack` (v0.1.3.2 and earlier versions), the recommended procedure was to create a STAR file and rebuild a stack manually. You can still follow this approach if you prefer scripting.

Firstly, locate where your particle stack currently is. For example you need to import a full stack which is result from NU refinement job (e.g. J44), the file you are looking for is the:
J44_006_particles.cs, which contains the locations of the particles in the cryoSPARC directory. 

(1) **Convert .cs file to a STAR file**. Using pyem’s conversion script, run: `csparc2star.py J44_006_particles.cs J44_006_particles.star`. 
Alternatively, from JANAS v0.1.2 onwards you may use the built-in utility:
```
janas_utils csparc2star test/J44_006_particles.cs test/J44_006_particles.star
```

(2)  **Updating to reference .mrcs rather than .mrc.**
By default, cryoSPARC writes both single-particle volumes and image stacks with a .mrc extension.
RELION, however, requires that 2D particle stacks use the .mrcs extension to distinguish them from 3D volumes. This naming requirement is specific to RELION’s file-handling routines; most other cryo-EM packages will accept a stack named .mrc.
To satisfy RELION’s requirement, create a symbolic link to the .mrc particle stack by appending an “s” to the filename and update the STAR file so that each _rlnImageName entry refers to the .mrcs link. RELION will then recognise the file as a 2D stack without error and without increasing disk usage, since the symbolic link does not duplicate data.

You can do it by getting a symbolic link for all the files to point to mrcs:

```
src_dir="/source/absolute/path/J1235/extract"
dst_dir="/destination/absolute/path/J1235/extract"

mkdir -p "$dst_dir"
for src_path in "$src_dir"/*.mrc; do
    filename="${src_path##*/}"
    newname="${filename%.mrc}.mrcs"
    ln -s "$src_path" "$dst_dir/$newname"
done
```


and fix the star file, by opening the star file with any decent text editor (gedit on linux works fine) and add an s to all the mrc filenames in the _rlnImageName parameters (something like "edit -> replace all" in the text editor). Or using regular expressions just do:

```
sed -i 's/_particles\.mrc /_particles.mrcs /g' J1235_particles_exported.star
```

(3) **create the stack with relion using the command:**
Although not strictly necessary, it is beneficial to have a single stack of particles. This can be created with 
relion_stack_create --i J44_006_particles.star --o full_stack --one_by_one

# Update STAR file from cryoSPARC processing

You may wish to update an existing RELION STAR file with parameters from a cryoSPARC job (e.g. angles, shifts, CTF).
This preserves the original data_optics block and rewrites only the data_particles rows that are present in the input STAR header.
```
janas_utils update_from_csparc <input_cs> <input_star> <output_star>
```
where `<input_cs>` is the cryoSPARC .cs file, `<input_star>` is the STAR you want to update, and `<output_star>` is the path for the updated STAR.


Only columns that already exist in the data_particles header are modified. If a column is not present, it is ignored. The updater recognises and updates the following when available: `_rlnAngleRot`, `_rlnAngleTilt`, `_rlnAnglePsi`, `_rlnOriginXAngst`, `_rlnOriginYAngst`, `_rlnDefocusU`, `_rlnDefocusV`, `_rlnDefocusAngle`, `_rlnPhaseShift`, `_rlnCtfBfactor`, `_rlnOpticsGroup`, `_rlnRandomSubset`, and `_rlnClassNumber`.

The command will terminate with an error if the `.cs` array length and the number of particle rows in data_particles do not match. This check is strict by design to prevent row misalignment. If you encounter an off-by-one mismatch, inspect the end of the data_particles section and remove any trailing blank or comment lines so that only valid rows remain after the header; then rerun the command.

The `.cs` file can be downloaded from the cryoSPARC job that produced the parameters you wish to carry over (for example, Non-uniform refinement). In the job view, open the Particles output and export the particle metadata .cs file. The figure below shows the export location.

<img src="tutorial/tutorial_figures/export_cryosparc_file.png" alt="Export particle .cs from cryoSPARC job outputs" height="250"/>

The STAR and `.cs` files must refer to the same particle stack and be in the same order. The tool does not perform re-indexing or re-ordering; it updates row-by-row. Ensure that the number of particles is identical and that no filtering, re-sorting, or concatenation has changed row order between the two files.

On success, `<output_star>` contains the original data_optics section exactly as in `<input_star>` and a data_particles section whose existing columns have been replaced with values derived from `<input_cs>`, formatted with fixed numeric precision for reproducibility.
