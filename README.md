# gizmo_carver

[![Python](https://img.shields.io/badge/python-3.9-blue)](https://www.python.org/downloads/)
[![yt](https://img.shields.io/badge/yt-4.0.2-blue)](https://yt-project.org/)
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
is present within the working directory. Then, simply run the `main_gizmo_carver.py` 
file to generate the complete set of RADMC-3D input files for your chosen parameters. 

After running RADMC-3D on these files, the `radmc_image_processing.py` script can be used
to generate the moment 0, 1, and 2 maps of the output image file.

For more complete instructions, see the [full documentation](https://github.com/seafen7/gizmo_carver/blob/main/doc/FULL_DOC.md).
