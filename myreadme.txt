#steps to follow to reproduce the CO/H2 conversion factor maps shown in Sharda, Offner, Grudic et al. 2023
1.) Choose a starforge simulation and snapshot to analyze: https://docs.google.com/spreadsheets/d/1TGFJXnQW9GPp-9w6VFhyEW-iTb0W0i2KpsD01bt7QX4/edit#gid=750548296

2.) Get the snapshot from frontera on gadi (if you want to do the analysis on gadi)

3.) Run the jupyter notebook to get basic info of the snapshot. This uses both yt and meshoid.

4.) Run the ray tracing script in pytreegrav to obtain the column density of H2 for a starforge snapshot: https://github.com/mikegrudic/pytreegrav#ray-tracing
    Since this has to be done only once per snapshot, you can afford to use 100 rays with randomized directions so that you can decrease the noise and make it less correlated
    Example script (to be run on interactive node with 48 cores if using 100 rays): get_coldens_from_pytreegrav.py (in gizmo_carver/default_files)

5.) Initialize despotic cloud available in despotic/cloudfiles/basic_starforge.desp

6.) Now you want to find the abundance of CO based on different methods:
    a.) Simply assume CO abundance is 1e-4 everywhere
    b.) Use the CO abundance from Priestely+2023 (they postprocess a sim using UCLCHEM; include CO freezeout properly)
        Interpolate on their nh, NH to get CO
    c.) Use the column densities from pytreegrav in step 5, and run the NL99_GC chemical network in despotic to predict the CO abundance (starting from some initial abundance of C, O, H). Mimic CO freezeout by hand.
        Step c is very costly to do for each starforge particle. So, we create a 5D model grid in volume density, column density, ISRF, H2 abundance and gas temperature (these 5 matter the most for CO abundance)
        and interpolate across the grid to find the CO abundances. Works well if the grid is dense in volume and column densities (these 2 matter the most among these 5)
    d.) Same as step c but using UCLCHEM
    Steps c and d can be done by script do_chemistry.py (or, do_chemistry_parallel.py for more speed, but only run it on one node)  available in gizmo_carver/default_files
    Save the CO abundances as a txt file

7.) Load the CO abundances as a new field by creatting a new .HDF5 using h5py (yt doesn't do it for GIZMO datasets). This new field will be dimensionless in yt. 
    Example script: add_molecules_to_data.py (in gizmo_carver/default_files/)
    I name the modified file the same as the original one (if you want to keep the original, rename it) so that gizmo_carver scripts below do not need to be fiddled with

8.) Set the input parameters in the file input_gizmo_carver.py in gizmo_carver/src/. Checkout new_fields_list.py because thats where we add column density (in cm**-2) as a yt field.

9.) Run wrapper_gizmo_carver.py

10.) Step 9 will produce some output in the specified output directory. In particular, it will produce a file called numberdens_XXX.inp from the 3 different methods above in step 6.
     Make sure to rename the molecular numdens file you want to use to numberdens_co.inp as RADMC expects LAMDA Leiden database naming

11.) Run radmc3d to produce image.out using numberdens_co.inp takes 1.5 hrs, submit as a job)
     Example script: radmc3d image npix 256 loadlambda fluxcons doppcatch inclline linelist nostar writepop sizepc 5.0 phi 0 incl 0 | tee output.txt
     Parameter values '256' and '5.0' should be the same as box_dim and box_size in the input file in step 8

12.) Step 11 will produce an image.out file. Run radmc_moments.py (in gizmo_carver/default_files) on this file to produce txt files that store the mom0 information. This will write the 2D mom0 map as a txt file.
12.a) To produce mock observations, you can convolve the synthetic map with a beam and place the observe at distance > 0. Can do this via radmc_moments_obsv.py (in gizmo_carver/default_files/)

13.) To get an observational X_CO, you also need the column density along the direction you have I_CO (mom0 map of CO) from in step 12. Use the file recreate_cube.py to get the column density along the LOS
     by integrating the volume density of H2. save it as a txt file.
     NOTE: this column density is not the same as the column density used in step 6: that was local to the particle. This is along the LOS to the observer.

14.) Ratio of the mom0 map of CO 1-0 to H2 column density will give you the conversion factor X_CO

gizmo_carver original documentation below:

#########################################################
# gizmo_carver

[![Python](https://img.shields.io/badge/python-3.9-blue)](https://www.python.org/downloads/)
[![yt](https://img.shields.io/badge/yt-4.0.2-blue)](https://yt-project.org/)
[![yt_astro_analysis](https://yt-astro-analysis.readthedocs.io/en/latest/)](https://yt-project.org/)
[![radmc3dPy](https://img.shields.io/badge/radmc3dPy-0.30.2-blue)](https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_rmcpy/index.html)
[![matplotlib](https://img.shields.io/badge/matplotlib-3.5.0-blue)](https://matplotlib.org/)
[![numpy](https://img.shields.io/badge/numpy-4.0.2-blue)](https://numpy.org/)

Python scripts for generating synthetic observations from a GIZMO dataset
using yt and RADMC-3D

## File Contents

### globals_gizmo_carver.py

Contains constants and simple functions for use in the RADMC-3D carving routines.
Should not need to edit this file.

### inputs_gizmo_carver.py

Input file containing constants and parameters for RADMC-3D carving routines. 
This is the only file that requires editing during general use. 

### main_gizmo_carver.py

Driver file for RADMC carve routines. Generates necessary input files 
for RADMC-3D. Call this file to run the routine. Should not need to edit this file.

### writer_gizmo_carver.py

Writer class that compresses the relevant layers needed to create RADMC-3D 
amr, number density, line, and dust files. Contains custom covering_grid
function for the Gizmo file format. Should not need to edit this file.

### radmc_image_processing.py

Contains functions for plotting RADMC-3D output image files using RadMC3DPy 
and MatPlotLib. Edit this file to create custom plotting routines.

### image.py

Contains functions for handling image.out files, including calculating and plotting 
moment maps. Needs to be added to radmc3dPy python tools directory (then run python setup.py install)

### default_files

Contains input files for RADMC-3D that are not generated by the carver routine.
These files are still necessary for running RADMC-3D and are used by the carver routine, 
but do not depend on the contents of the dataset file. As such, modifications to 
these files may be needed for certain observation parameters.

## Getting Started

After cloning this repo to your local machine, modify the `inputs_gizmo_carver.py` file
to match the desired input and output parameters. Also ensure that the GIZMO dataset file
is present within the working directory. Copy your compiled `radmc3d` executable into the
default_files directory. Add your radmc submission command to `submit_script.sh`.
Then, simply run the `main_gizmo_carver.py` file to generate the complete set of 
RADMC-3D input files for your chosen parameters. 

After running RADMC-3D (`submit_script.sh`) on these files, the `radmc_image_processing.py` script can be used
to generate the moment 0, 1, and 2 maps of the output image file.

For more complete instructions, see the [full documentation](https://github.com/seafen7/gizmo_carver/blob/main/doc/FULL_DOC.md).

########################################################
