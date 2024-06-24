"""
   main_gizmo_carver.py

   Purpose:
        Driver file for RADMC carve routines. Call this function when you want
        to run the code. Should not ever need to edit this file.
        See inputs_gizmo_carver.py for parameter defaults

   Authors:
        Sean Feng, feng.sean01@utexas.edu
        Stella Offner
        Spring 2022
        
        Modified from: main_CarveOut.py, written by:
        Aaron T. Lee, aaron.t.lee@utexas.edu
        Spring 2018

   Written/Tested with Python 3.9, yt 4.0.2
"""

import yt
from writer_gizmo_carver import RadMC3DWriter_Gizmo
from globals_gizmo_carver import *
import inputs_gizmo_carver as inputs
from yt.units import *
import os
from datetime import datetime
import shutil
import numpy as np
import glob as glob
from new_fields_list import *
import sys
import h5py as h5py
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor, as_completed

# Overload box_center, snapshot number, and directory in case 
# we want to loop through snapshots or star locations:
def main_gizmo_carver(box_center=inputs.box_center, snap=inputs.snap, hdf5_dir=inputs.hdf5_dir, tag=inputs.tag):
    
    snap = str(snap).zfill(3)
    hdf5_file =  hdf5_dir+'snapshot_'+snap+'.hdf5'
    tag = 'sn'+snap+'_'+ np.str(np.int(np.imag(inputs.box_dim)))+'_'

    if inputs.mask_abundance == True:
        # Need to read and concatenate these files of accreted particles
        fns = glob.glob(hdf5_dir+"/blackhole_details/bhswallow*.txt")
      
        alllines = np.empty((0,19))
        for file in fns:
            lines = np.loadtxt(file, comments="#", delimiter=" ", unpack=False)
            #print(np.shape(alllines), np.shape(lines), file, lines)
            alllines=np.concatenate((alllines, lines), axis=0)

        #yt.add_field(('PartType0', 'Mask'), function=_Mask, units='g/cm**3', sampling_type='particle', force_override=True)

        # Must define this here so file information is in scope
        def _MaskedMolecularNumDensity(field, data):
            mask = data[('PartType0', 'ParticleIDs')]*0.0
            exist = np.where(np.isin(alllines[:,1], data[('PartType5','ParticleIDs')])== True)[0]  # Find all the rows of particles that exist here
            #print(" len of existing rows =", len(exist))
            ind = np.where(np.isin(data[('PartType0','ParticleIDs')], alllines[exist,6])== True)[0]
            #print(" len of accreted particles =", len(ind))
            mask[ind] = 1
            return data[('PartType0', 'H2NumDensity')]*inputs.molecular_abundance*( data[('PartType0','gas_temperature')] < yt.YTArray([3e3], "K"))*mask

        yt.add_field(('PartType0', 'MaskedMolecularNumDensity'), function=_MaskedMolecularNumDensity, units='cm**-3', sampling_type='particle', force_override=True)

#     def _H2NumDensity(field, data):
#         vel_x = data['PartType0','velocity_x'].to('km/s')
#         vel_y = data['PartType0','velocity_y'].to('km/s')
#         vel_z = data['PartType0','velocity_z'].to('km/s')
#         v_rms = (np.sqrt(vel_x**2 + vel_y**2 + vel_z**2))
#         return data[('PartType0', 'Density')]*data[('PartType0', 'MolecularMassFraction')]*data[('PartType0', 'NeutralHydrogenAbundance')]*(1-inputs.helium_mass_fraction)/(inputs.mol_hydrogen_ratio*mh)*(data[('PartType0','H2NumDensity')] < yt.YTArray([1e4], "cm**-3"))*( data[('PartType0','gas_temperature')] < yt.YTArray([3e3], "K"))*(v_rms < yt.YTArray([1e5], "cm/s"))

#     yt.add_field(('PartType0', 'H2MolecularNumDensity'), function=_H2NumDensity, units='cm**-3', sampling_type='particle', force_override=True)    ##updated by me
    
    # Loads file
    ds = yt.load(hdf5_file, unit_base=inputs.unit_base)
    print("domain = ", ds.domain_left_edge, ds.domain_right_edge)

    try:
        print("Loaded file " + str(ds)) 
    except NameError:
        assert False, "YT unable to properly load file!"
        
    now = datetime.now()
    dt_string = now.strftime("%m.%d.%y_%H.%M.%S")

    # Create working directory for this run
    current_dir = inputs.output_filepath + 'snapshot_' + str(int(snap))
    working_dir_name = current_dir#os.path.join(current_dir, 'RADMC_inputs_' + tag+dt_string)
    if not os.path.exists(working_dir_name):
        os.mkdir(working_dir_name)

    # Make a file to store I/O and setup parameters
    f = open(os.path.join(working_dir_name, inputs.out_makeinput), 'w')
    original_stdout = sys.stdout #Reset do: sys.stdout = original_stdout 
    sys.stdout = f
    print(">> INPUT PARAMETERS (inputs_gizmo_carver.py): <<")
    print("   dust_to_gas ", inputs.dust_to_gas)
    print("   molecular_abundance ", inputs.molecular_abundance)
    print("   mask abundance ", inputs.mask_abundance)
    print("   box_size, box_dim, box_center ", inputs.box_size, inputs.box_dim, box_center)
    print("   hdf5_file ", hdf5_file)
    f = h5py.File(hdf5_file, 'r') # Retrieve some basic information
    print("  Simulation time [code units] = ", f['Header'].attrs['Time'])
    if "PartType5" in f:
        part = f['PartType5']
        print("  Number of stars =", len(f['PartType5']['StellarFormationTime']))

    print("\n>> SETUP OUTPUT (main_gizmo_carver): <<")
    
    box_left = np.add(box_center, -inputs.box_size)
    box_right = np.add(box_center, inputs.box_size)
    box_left_cgs = [Convert(x, inputs.box_units, 'cm', 'cm') for x in box_left]
    box_right_cgs = [Convert(x, inputs.box_units, 'cm', 'cm') for x in box_right]
    box_left_cgs = unyt_array(box_left_cgs, 'cm')   # Format required by yt
    box_right_cgs = unyt_array(box_right_cgs, 'cm')

    print("\n Carving between Left = " + str(box_left_cgs))
    print("            to Right = " + str(box_right_cgs))
    print("       w/ Resolution = " + str(abs(inputs.box_dim)) + " x " + str(abs(inputs.box_dim)) + "\n")

    writer = RadMC3DWriter_Gizmo(ds, a_boxLeft=box_left_cgs, a_boxRight=box_right_cgs, a_boxDim=inputs.box_dim)
    velocity_fields = [("PartType0","velocity_x"), ("PartType0","velocity_y"), ("PartType0","velocity_z")]
    # Write the amr grid file (fast)
    print("1/7: Writing amr grid file (fast!)")
    writer.write_amr_grid(os.path.join(working_dir_name, inputs.out_afname))
    
    if inputs.mask_abundance == True:        
        tasks = [((writer.write_line_file),(('PartType0', 'MaskedMolecularNumDensity'), os.path.join(working_dir_name, inputs.out_nfname))),
            ((writer.write_dust_file),(('PartType0', 'DustDensity'), os.path.join(working_dir_name, inputs.out_ddfname))),
            ((writer.write_line_file),(("PartType0", "gas_temperature"), os.path.join(working_dir_name, inputs.out_tfname))),
            ((writer.write_dust_file),(("PartType0", "dust_temperature"), os.path.join(working_dir_name, inputs.out_dtfname))),
            ((writer.write_line_file),(velocity_fields, os.path.join(working_dir_name, inputs.out_vfname))),
            ((writer.write_line_file),(('PartType0', 'H2NumDensity'), os.path.join(working_dir_name, inputs.out_cfname)))]
    else: 
        tasks = [((writer.write_line_file),(('PartType0', 'MolecularNumDensity'), os.path.join(working_dir_name, inputs.out_nfname))),
            ((writer.write_dust_file),(('PartType0', 'DustDensity'), os.path.join(working_dir_name, inputs.out_ddfname))),
            ((writer.write_line_file),(("PartType0", "gas_temperature"), os.path.join(working_dir_name, inputs.out_tfname))),
            ((writer.write_dust_file),(("PartType0", "dust_temperature"), os.path.join(working_dir_name, inputs.out_dtfname))),
            ((writer.write_line_file),(velocity_fields, os.path.join(working_dir_name, inputs.out_vfname))),
            ((writer.write_line_file),(('PartType0', 'H2NumDensity'), os.path.join(working_dir_name, inputs.out_cfname)))]

    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(task, *args) for task, args in tasks]
    for future in as_completed(futures):
        try:
            future.result()  # Check for exceptions
        except Exception as exc:
            print(f'Generated an exception: {exc}')

    # Copy over existing files
    print('Copying default files...')
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_molname), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_wlmname), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_cwlname), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_dksname), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_dtpname), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_linname), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_rmcname), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_apfname), working_dir_name)
    # shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_execute), working_dir_name)
    # shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_subscript), working_dir_name)
        
    print('Done! Output files generated at: \n\n' + os.path.abspath(working_dir_name))

    return working_dir_name

main_gizmo_carver()
