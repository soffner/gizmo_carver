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
   Further modified by Piyush Sharda to include chemistry with despotic and uclchem (2023)
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

# Overload box_center, snapshot number, and directory in case 
# we want to loop through snapshots or star locations:

def main_gizmo_carver(box_center=inputs.box_center, snap=inputs.snap, hdf5_dir=inputs.hdf5_dir, tag=inputs.tag):

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
            return data[('PartType0', 'H2NumDensity')]*inputs.molecular_abundance*(data[('PartType0','H2NumDensity')] > yt.YTArray([1e2], "cm**-3"))*( data[('PartType0','gas_temperature')] < yt.YTArray([1e2], "K"))*mask

        yt.add_field(('PartType0', 'MaskedMolecularNumDensity'), function=_MaskedMolecularNumDensity, units='cm**-3', sampling_type='particle', force_override=True)

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
    current_dir = inputs.output_filepath
    working_dir_name = os.path.join(current_dir, 'RADMC_inputs_' + tag+dt_string)
    os.mkdir(working_dir_name)

    #new1_dir_name = os.path.join(working_dir_name, 'default')
    #os.mkdir(new1_dir_name)

    new2_dir_name = os.path.join(working_dir_name, 'despotic')
    os.mkdir(new2_dir_name)

    new3_dir_name = os.path.join(working_dir_name, 'uclchem')
    os.mkdir(new3_dir_name)

    lines_dir_name = os.path.join(working_dir_name, 'lines')
    os.mkdir(lines_dir_name)

    temps_dir_name = os.path.join(working_dir_name, 'temperatures')
    os.mkdir(temps_dir_name)

    CO_J10_dir_name = os.path.join(working_dir_name, 'CO_J10')
    os.mkdir(CO_J10_dir_name)

    CO_J21_dir_name = os.path.join(working_dir_name, 'CO_J21')
    os.mkdir(CO_J21_dir_name)

    CO_J32_dir_name = os.path.join(working_dir_name, 'CO_J32')
    os.mkdir(CO_J32_dir_name)

    CO_J43_dir_name = os.path.join(working_dir_name, 'CO_J43')
    os.mkdir(CO_J43_dir_name)

    CO_J54_dir_name = os.path.join(working_dir_name, 'CO_J54')
    os.mkdir(CO_J54_dir_name)

    CO_J65_dir_name = os.path.join(working_dir_name, 'CO_J65')
    os.mkdir(CO_J65_dir_name)

    CO_J76_dir_name = os.path.join(working_dir_name, 'CO_J76')
    os.mkdir(CO_J76_dir_name)

    CO_J87_dir_name = os.path.join(working_dir_name, 'CO_J87')
    os.mkdir(CO_J87_dir_name)

    CO_J98_dir_name = os.path.join(working_dir_name, 'CO_J98')
    os.mkdir(CO_J98_dir_name)

    CO_J109_dir_name = os.path.join(working_dir_name, 'CO_J109')
    os.mkdir(CO_J109_dir_name)

    CO13_J10_dir_name = os.path.join(working_dir_name, '13CO_J10')
    os.mkdir(CO13_J10_dir_name)

    CO13_J21_dir_name = os.path.join(working_dir_name, '13CO_J21')
    os.mkdir(CO13_J21_dir_name)

    CO13_J32_dir_name = os.path.join(working_dir_name, '13CO_J32')
    os.mkdir(CO13_J32_dir_name)

    CO13_J43_dir_name = os.path.join(working_dir_name, '13CO_J43')
    os.mkdir(CO13_J43_dir_name)

    CO13_J54_dir_name = os.path.join(working_dir_name, '13CO_J54')
    os.mkdir(CO13_J54_dir_name)

    CO13_J65_dir_name = os.path.join(working_dir_name, '13CO_J65')
    os.mkdir(CO13_J65_dir_name)

    CO13_J76_dir_name = os.path.join(working_dir_name, '13CO_J76')
    os.mkdir(CO13_J76_dir_name)

    CO13_J87_dir_name = os.path.join(working_dir_name, '13CO_J87')
    os.mkdir(CO13_J87_dir_name)

    CO13_J98_dir_name = os.path.join(working_dir_name, '13CO_J98')
    os.mkdir(CO13_J98_dir_name)

    CO13_J109_dir_name = os.path.join(working_dir_name, '13CO_J109')
    os.mkdir(CO13_J109_dir_name)

    C18O_J10_dir_name = os.path.join(working_dir_name, 'C18O_J10')
    os.mkdir(C18O_J10_dir_name)

    C18O_J21_dir_name = os.path.join(working_dir_name, 'C18O_J21')
    os.mkdir(C18O_J21_dir_name)

    C18O_J32_dir_name = os.path.join(working_dir_name, 'C18O_J32')
    os.mkdir(C18O_J32_dir_name)

    C18O_J43_dir_name = os.path.join(working_dir_name, 'C18O_J43')
    os.mkdir(C18O_J43_dir_name)

    C18O_J54_dir_name = os.path.join(working_dir_name, 'C18O_J54')
    os.mkdir(C18O_J54_dir_name)

    C18O_J65_dir_name = os.path.join(working_dir_name, 'C18O_J65')
    os.mkdir(C18O_J65_dir_name)

    C18O_J76_dir_name = os.path.join(working_dir_name, 'C18O_J76')
    os.mkdir(C18O_J76_dir_name)

    C18O_J87_dir_name = os.path.join(working_dir_name, 'C18O_J87')
    os.mkdir(C18O_J87_dir_name)

    C18O_J98_dir_name = os.path.join(working_dir_name, 'C18O_J98')
    os.mkdir(C18O_J98_dir_name)

    C18O_J109_dir_name = os.path.join(working_dir_name, 'C18O_J109')
    os.mkdir(C18O_J109_dir_name)

    HCOp_J10_dir_name = os.path.join(working_dir_name, 'HCOp_J10')
    os.mkdir(HCOp_J10_dir_name)

    HCOp_J21_dir_name = os.path.join(working_dir_name, 'HCOp_J21')
    os.mkdir(HCOp_J21_dir_name)

    HCOp_J32_dir_name = os.path.join(working_dir_name, 'HCOp_J32')
    os.mkdir(HCOp_J32_dir_name)

    HCOp_J43_dir_name = os.path.join(working_dir_name, 'HCOp_J43')
    os.mkdir(HCOp_J43_dir_name)

    HCOp_J54_dir_name = os.path.join(working_dir_name, 'HCOp_J54')
    os.mkdir(HCOp_J54_dir_name)

    HCOp_J65_dir_name = os.path.join(working_dir_name, 'HCOp_J65')
    os.mkdir(HCOp_J65_dir_name)

    HCOp_J76_dir_name = os.path.join(working_dir_name, 'HCOp_J76')
    os.mkdir(HCOp_J76_dir_name)

    HCOp_J87_dir_name = os.path.join(working_dir_name, 'HCOp_J87')
    os.mkdir(HCOp_J87_dir_name)

    HCOp_J98_dir_name = os.path.join(working_dir_name, 'HCOp_J98')
    os.mkdir(HCOp_J98_dir_name)

    HCOp_J109_dir_name = os.path.join(working_dir_name, 'HCOp_J109')
    os.mkdir(HCOp_J109_dir_name)

    H13COp_J10_dir_name = os.path.join(working_dir_name, 'H13COp_J10')
    os.mkdir(H13COp_J10_dir_name)

    H13COp_J21_dir_name = os.path.join(working_dir_name, 'H13COp_J21')
    os.mkdir(H13COp_J21_dir_name)

    H13COp_J32_dir_name = os.path.join(working_dir_name, 'H13COp_J32')
    os.mkdir(H13COp_J32_dir_name)

    H13COp_J43_dir_name = os.path.join(working_dir_name, 'H13COp_J43')
    os.mkdir(H13COp_J43_dir_name)

    H13COp_J54_dir_name = os.path.join(working_dir_name, 'H13COp_J54')
    os.mkdir(H13COp_J54_dir_name)

    H13COp_J65_dir_name = os.path.join(working_dir_name, 'H13COp_J65')
    os.mkdir(H13COp_J65_dir_name)

    H13COp_J76_dir_name = os.path.join(working_dir_name, 'H13COp_J76')
    os.mkdir(H13COp_J76_dir_name)

    H13COp_J87_dir_name = os.path.join(working_dir_name, 'H13COp_J87')
    os.mkdir(H13COp_J87_dir_name)

    H13COp_J98_dir_name = os.path.join(working_dir_name, 'H13COp_J98')
    os.mkdir(H13COp_J98_dir_name)

    H13COp_J109_dir_name = os.path.join(working_dir_name, 'H13COp_J109')
    os.mkdir(H13COp_J109_dir_name)

    HCN_J10_dir_name = os.path.join(working_dir_name, 'HCN_J10')
    os.mkdir(HCN_J10_dir_name)

    HCN_J21_dir_name = os.path.join(working_dir_name, 'HCN_J21')
    os.mkdir(HCN_J21_dir_name)

    HCN_J32_dir_name = os.path.join(working_dir_name, 'HCN_J32')
    os.mkdir(HCN_J32_dir_name)

    HCN_J43_dir_name = os.path.join(working_dir_name, 'HCN_J43')
    os.mkdir(HCN_J43_dir_name)

    HCN_J54_dir_name = os.path.join(working_dir_name, 'HCN_J54')
    os.mkdir(HCN_J54_dir_name)

    HCN_J65_dir_name = os.path.join(working_dir_name, 'HCN_J65')
    os.mkdir(HCN_J65_dir_name)

    HCN_J76_dir_name = os.path.join(working_dir_name, 'HCN_J76')
    os.mkdir(HCN_J76_dir_name)

    HCN_J87_dir_name = os.path.join(working_dir_name, 'HCN_J87')
    os.mkdir(HCN_J87_dir_name)

    HCN_J98_dir_name = os.path.join(working_dir_name, 'HCN_J98')
    os.mkdir(HCN_J98_dir_name)

    HCN_J109_dir_name = os.path.join(working_dir_name, 'HCN_J109')
    os.mkdir(HCN_J109_dir_name)

    HNC_J10_dir_name = os.path.join(working_dir_name, 'HNC_J10')
    os.mkdir(HNC_J10_dir_name)

    HNC_J21_dir_name = os.path.join(working_dir_name, 'HNC_J21')
    os.mkdir(HNC_J21_dir_name)

    HNC_J32_dir_name = os.path.join(working_dir_name, 'HNC_J32')
    os.mkdir(HNC_J32_dir_name)

    HNC_J43_dir_name = os.path.join(working_dir_name, 'HNC_J43')
    os.mkdir(HNC_J43_dir_name)

    HNC_J54_dir_name = os.path.join(working_dir_name, 'HNC_J54')
    os.mkdir(HNC_J54_dir_name)

    HNC_J65_dir_name = os.path.join(working_dir_name, 'HNC_J65')
    os.mkdir(HNC_J65_dir_name)

    HNC_J76_dir_name = os.path.join(working_dir_name, 'HNC_J76')
    os.mkdir(HNC_J76_dir_name)

    HNC_J87_dir_name = os.path.join(working_dir_name, 'HNC_J87')
    os.mkdir(HNC_J87_dir_name)

    HNC_J98_dir_name = os.path.join(working_dir_name, 'HNC_J98')
    os.mkdir(HNC_J98_dir_name)

    HNC_J109_dir_name = os.path.join(working_dir_name, 'HNC_J109')
    os.mkdir(HNC_J109_dir_name)


    # Make a file to store I/O and setup parameters
    f = open(os.path.join(working_dir_name, inputs.out_makeinput), 'w')
    original_stdout = sys.stdout #Reset do: sys.stdout = original_stdout 
    sys.stdout = f
    print(">> INPUT PARAMETERS (inputs_gizmo_carver.py): <<")
    print("   dust_to_gas ", inputs.dust_to_gas)
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

    # Write the amr grid file (fast)
    print("1/7: Writing amr grid file (fast!)")
    writer.write_amr_grid(os.path.join(working_dir_name, inputs.out_afname))

    # Write the number density file for species (slow)
    print("2/7: Writing number density file (slow.)")
    if inputs.mask_abundance == True:
        writer.write_line_file(('PartType0', 'MaskedMolecularNumDensity'), os.path.join(working_dir_name, inputs.out_nfname))
    else:
        #writer.write_line_file(('PartType0', 'CONumberDensityDefault'), os.path.join(new1_dir_name, inputs.out_codefaultname))
        writer.write_line_file(('PartType0', 'CONumberDensityDespotic'), os.path.join(new2_dir_name, inputs.out_codespname))
        writer.write_line_file(('PartType0', '13CONumberDensityDespotic'), os.path.join(new2_dir_name, inputs.out_13codespname))
        writer.write_line_file(('PartType0', 'C18ONumberDensityDespotic'), os.path.join(new2_dir_name, inputs.out_c18odespname))
        writer.write_line_file(('PartType0', 'HCOpNumberDensityDespotic'), os.path.join(new2_dir_name, inputs.out_hcopdespname))
        writer.write_line_file(('PartType0', 'H13COpNumberDensityDespotic'), os.path.join(new2_dir_name, inputs.out_h13copdespname))
        writer.write_line_file(('PartType0', 'CONumberDensityUCLCHEM'), os.path.join(new3_dir_name, inputs.out_couclname))
        writer.write_line_file(('PartType0', '13CONumberDensityUCLCHEM'), os.path.join(new3_dir_name, inputs.out_13couclname))
        writer.write_line_file(('PartType0', 'C18ONumberDensityUCLCHEM'), os.path.join(new3_dir_name, inputs.out_c18ouclname))
        writer.write_line_file(('PartType0', 'HCOpNumberDensityUCLCHEM'), os.path.join(new3_dir_name, inputs.out_hcopuclname))
        writer.write_line_file(('PartType0', 'H13COpNumberDensityUCLCHEM'), os.path.join(new3_dir_name, inputs.out_h13copuclname))
        writer.write_line_file(('PartType0', 'HCNNumberDensityUCLCHEM'), os.path.join(new3_dir_name, inputs.out_hcnuclname))
        writer.write_line_file(('PartType0', 'HNCNumberDensityUCLCHEM'), os.path.join(new3_dir_name, inputs.out_hncuclname))

    #writer.write_line_file(('PartType0', 'H2ColumnDensity'), os.path.join(working_dir_name, inputs.out_NHname))
    writer.write_line_file(('PartType0', 'H2NumDensity'), os.path.join(working_dir_name, inputs.out_nhname))

    # Write the dust density file for dust (slow)
    print("3/7: Writing dust density file (slow.)")
    writer.write_dust_file(('PartType0', 'DustDensity'), os.path.join(working_dir_name, inputs.out_ddfname))

    print("4.5/7: Writing microturbulence file (slow.)")
    writer.write_line_file(("PartType0", "microturbulence_speed"), os.path.join(working_dir_name, inputs.out_mtfname))

    # Write the temperature file for species or dust (slow)
    print("5/7: Writing temperature file (slower..)")
    #print(ds.find_max(("PartType0","gas_temperature_CO")))
    writer.write_line_file(("PartType0", "gas_temperature_CO"), os.path.join(temps_dir_name, inputs.out_tfname_CO))
    #print(ds.find_max(("PartType0","gas_temperature_HCOp")))
    writer.write_line_file(("PartType0", "gas_temperature_HCOp"), os.path.join(temps_dir_name, inputs.out_tfname_HCOp))
    #print(ds.find_max(("PartType0","gas_temperature_HCN")))
    writer.write_line_file(("PartType0", "gas_temperature_HCN"), os.path.join(temps_dir_name, inputs.out_tfname_HCN))
    writer.write_line_file(("PartType0", "gas_temperature_HNC"), os.path.join(temps_dir_name, inputs.out_tfname_HNC))

    # Assuming dust temperature is same as gas for now...
    print("6/7: Writing dust temperature file (slower..)")
    writer.write_dust_file(("PartType0", "dust_temperature"), os.path.join(working_dir_name, inputs.out_dtfname))

    # Write the gas velocity file (slow x 3)
    print("7/7: Writing velocity file (slowest...)")
    velocity_fields = [("PartType0","velocity_x"), ("PartType0", "velocity_y"), ("PartType0", "velocity_z")]
    writer.write_line_file(velocity_fields, os.path.join(working_dir_name, inputs.out_vfname))

    # Write the wavelength file
    if inputs.write_line_file == True:
        print("Writing line file ...")
        #file = os.path.join(working_dir_name, inputs.out_cwlname)
        #print("parameters =", inputs.vmax, inputs.nwav, inputs.restfreq, len(file) )
        writer.writeMolLineFile(os.path.join(CO_J10_dir_name, inputs.out_cwlname_CO_J10), inputs.vmax, inputs.nwav, inputs.restfreq_CO_J10)
        writer.writeMolLineFile(os.path.join(CO_J21_dir_name, inputs.out_cwlname_CO_J21), inputs.vmax, inputs.nwav, inputs.restfreq_CO_J21)
        writer.writeMolLineFile(os.path.join(CO_J32_dir_name, inputs.out_cwlname_CO_J32), inputs.vmax, inputs.nwav, inputs.restfreq_CO_J32)
        writer.writeMolLineFile(os.path.join(CO_J43_dir_name, inputs.out_cwlname_CO_J43), inputs.vmax, inputs.nwav, inputs.restfreq_CO_J43)
        writer.writeMolLineFile(os.path.join(CO_J54_dir_name, inputs.out_cwlname_CO_J54), inputs.vmax, inputs.nwav, inputs.restfreq_CO_J54)
        writer.writeMolLineFile(os.path.join(CO_J65_dir_name, inputs.out_cwlname_CO_J65), inputs.vmax, inputs.nwav, inputs.restfreq_CO_J65)
        writer.writeMolLineFile(os.path.join(CO_J76_dir_name, inputs.out_cwlname_CO_J76), inputs.vmax, inputs.nwav, inputs.restfreq_CO_J76)
        writer.writeMolLineFile(os.path.join(CO_J87_dir_name, inputs.out_cwlname_CO_J87), inputs.vmax, inputs.nwav, inputs.restfreq_CO_J87)
        writer.writeMolLineFile(os.path.join(CO_J98_dir_name, inputs.out_cwlname_CO_J98), inputs.vmax, inputs.nwav, inputs.restfreq_CO_J98)
        writer.writeMolLineFile(os.path.join(CO_J109_dir_name, inputs.out_cwlname_CO_J109), inputs.vmax, inputs.nwav, inputs.restfreq_CO_J109)

        writer.writeMolLineFile(os.path.join(CO13_J10_dir_name, inputs.out_cwlname_13CO_J10), inputs.vmax, inputs.nwav, inputs.restfreq_13CO_J10)
        writer.writeMolLineFile(os.path.join(CO13_J21_dir_name, inputs.out_cwlname_13CO_J21), inputs.vmax, inputs.nwav, inputs.restfreq_13CO_J21)
        writer.writeMolLineFile(os.path.join(CO13_J32_dir_name, inputs.out_cwlname_13CO_J32), inputs.vmax, inputs.nwav, inputs.restfreq_13CO_J32)
        writer.writeMolLineFile(os.path.join(CO13_J43_dir_name, inputs.out_cwlname_13CO_J43), inputs.vmax, inputs.nwav, inputs.restfreq_13CO_J43)
        writer.writeMolLineFile(os.path.join(CO13_J54_dir_name, inputs.out_cwlname_13CO_J54), inputs.vmax, inputs.nwav, inputs.restfreq_13CO_J54)
        writer.writeMolLineFile(os.path.join(CO13_J65_dir_name, inputs.out_cwlname_13CO_J65), inputs.vmax, inputs.nwav, inputs.restfreq_13CO_J65)
        writer.writeMolLineFile(os.path.join(CO13_J76_dir_name, inputs.out_cwlname_13CO_J76), inputs.vmax, inputs.nwav, inputs.restfreq_13CO_J76)
        writer.writeMolLineFile(os.path.join(CO13_J87_dir_name, inputs.out_cwlname_13CO_J87), inputs.vmax, inputs.nwav, inputs.restfreq_13CO_J87)
        writer.writeMolLineFile(os.path.join(CO13_J98_dir_name, inputs.out_cwlname_13CO_J98), inputs.vmax, inputs.nwav, inputs.restfreq_13CO_J98)
        writer.writeMolLineFile(os.path.join(CO13_J109_dir_name, inputs.out_cwlname_13CO_J109), inputs.vmax, inputs.nwav, inputs.restfreq_13CO_J109)

        writer.writeMolLineFile(os.path.join(C18O_J10_dir_name, inputs.out_cwlname_C18O_J10), inputs.vmax, inputs.nwav, inputs.restfreq_C18O_J10)
        writer.writeMolLineFile(os.path.join(C18O_J21_dir_name, inputs.out_cwlname_C18O_J21), inputs.vmax, inputs.nwav, inputs.restfreq_C18O_J21)
        writer.writeMolLineFile(os.path.join(C18O_J32_dir_name, inputs.out_cwlname_C18O_J32), inputs.vmax, inputs.nwav, inputs.restfreq_C18O_J32)
        writer.writeMolLineFile(os.path.join(C18O_J43_dir_name, inputs.out_cwlname_C18O_J43), inputs.vmax, inputs.nwav, inputs.restfreq_C18O_J43)
        writer.writeMolLineFile(os.path.join(C18O_J54_dir_name, inputs.out_cwlname_C18O_J54), inputs.vmax, inputs.nwav, inputs.restfreq_C18O_J54)
        writer.writeMolLineFile(os.path.join(C18O_J65_dir_name, inputs.out_cwlname_C18O_J65), inputs.vmax, inputs.nwav, inputs.restfreq_C18O_J65)
        writer.writeMolLineFile(os.path.join(C18O_J76_dir_name, inputs.out_cwlname_C18O_J76), inputs.vmax, inputs.nwav, inputs.restfreq_C18O_J76)
        writer.writeMolLineFile(os.path.join(C18O_J87_dir_name, inputs.out_cwlname_C18O_J87), inputs.vmax, inputs.nwav, inputs.restfreq_C18O_J87)
        writer.writeMolLineFile(os.path.join(C18O_J98_dir_name, inputs.out_cwlname_C18O_J98), inputs.vmax, inputs.nwav, inputs.restfreq_C18O_J98)
        writer.writeMolLineFile(os.path.join(C18O_J109_dir_name, inputs.out_cwlname_C18O_J109), inputs.vmax, inputs.nwav, inputs.restfreq_C18O_J109)

        writer.writeMolLineFile(os.path.join(HCOp_J10_dir_name, inputs.out_cwlname_HCOp_J10), inputs.vmax, inputs.nwav, inputs.restfreq_HCOp_J10)
        writer.writeMolLineFile(os.path.join(HCOp_J21_dir_name, inputs.out_cwlname_HCOp_J21), inputs.vmax, inputs.nwav, inputs.restfreq_HCOp_J21)
        writer.writeMolLineFile(os.path.join(HCOp_J32_dir_name, inputs.out_cwlname_HCOp_J32), inputs.vmax, inputs.nwav, inputs.restfreq_HCOp_J32)
        writer.writeMolLineFile(os.path.join(HCOp_J43_dir_name, inputs.out_cwlname_HCOp_J43), inputs.vmax, inputs.nwav, inputs.restfreq_HCOp_J43)
        writer.writeMolLineFile(os.path.join(HCOp_J54_dir_name, inputs.out_cwlname_HCOp_J54), inputs.vmax, inputs.nwav, inputs.restfreq_HCOp_J54)
        writer.writeMolLineFile(os.path.join(HCOp_J65_dir_name, inputs.out_cwlname_HCOp_J65), inputs.vmax, inputs.nwav, inputs.restfreq_HCOp_J65)
        writer.writeMolLineFile(os.path.join(HCOp_J76_dir_name, inputs.out_cwlname_HCOp_J76), inputs.vmax, inputs.nwav, inputs.restfreq_HCOp_J76)
        writer.writeMolLineFile(os.path.join(HCOp_J87_dir_name, inputs.out_cwlname_HCOp_J87), inputs.vmax, inputs.nwav, inputs.restfreq_HCOp_J87)
        writer.writeMolLineFile(os.path.join(HCOp_J98_dir_name, inputs.out_cwlname_HCOp_J98), inputs.vmax, inputs.nwav, inputs.restfreq_HCOp_J98)
        writer.writeMolLineFile(os.path.join(HCOp_J109_dir_name, inputs.out_cwlname_HCOp_J109), inputs.vmax, inputs.nwav, inputs.restfreq_HCOp_J109)

        writer.writeMolLineFile(os.path.join(H13COp_J10_dir_name, inputs.out_cwlname_H13COp_J10), inputs.vmax, inputs.nwav, inputs.restfreq_H13COp_J10)
        writer.writeMolLineFile(os.path.join(H13COp_J21_dir_name, inputs.out_cwlname_H13COp_J21), inputs.vmax, inputs.nwav, inputs.restfreq_H13COp_J21)
        writer.writeMolLineFile(os.path.join(H13COp_J32_dir_name, inputs.out_cwlname_H13COp_J32), inputs.vmax, inputs.nwav, inputs.restfreq_H13COp_J32)
        writer.writeMolLineFile(os.path.join(H13COp_J43_dir_name, inputs.out_cwlname_H13COp_J43), inputs.vmax, inputs.nwav, inputs.restfreq_H13COp_J43)
        writer.writeMolLineFile(os.path.join(H13COp_J54_dir_name, inputs.out_cwlname_H13COp_J54), inputs.vmax, inputs.nwav, inputs.restfreq_H13COp_J54)
        writer.writeMolLineFile(os.path.join(H13COp_J65_dir_name, inputs.out_cwlname_H13COp_J65), inputs.vmax, inputs.nwav, inputs.restfreq_H13COp_J65)
        writer.writeMolLineFile(os.path.join(H13COp_J76_dir_name, inputs.out_cwlname_H13COp_J76), inputs.vmax, inputs.nwav, inputs.restfreq_H13COp_J76)
        writer.writeMolLineFile(os.path.join(H13COp_J87_dir_name, inputs.out_cwlname_H13COp_J87), inputs.vmax, inputs.nwav, inputs.restfreq_H13COp_J87)
        writer.writeMolLineFile(os.path.join(H13COp_J98_dir_name, inputs.out_cwlname_H13COp_J98), inputs.vmax, inputs.nwav, inputs.restfreq_H13COp_J98)
        writer.writeMolLineFile(os.path.join(H13COp_J109_dir_name, inputs.out_cwlname_H13COp_J109), inputs.vmax, inputs.nwav, inputs.restfreq_H13COp_J109)

        writer.writeMolLineFile(os.path.join(HCN_J10_dir_name, inputs.out_cwlname_HCN_J10), inputs.vmax, inputs.nwav, inputs.restfreq_HCN_J10)
        writer.writeMolLineFile(os.path.join(HCN_J21_dir_name, inputs.out_cwlname_HCN_J21), inputs.vmax, inputs.nwav, inputs.restfreq_HCN_J21)
        writer.writeMolLineFile(os.path.join(HCN_J32_dir_name, inputs.out_cwlname_HCN_J32), inputs.vmax, inputs.nwav, inputs.restfreq_HCN_J32)
        writer.writeMolLineFile(os.path.join(HCN_J43_dir_name, inputs.out_cwlname_HCN_J43), inputs.vmax, inputs.nwav, inputs.restfreq_HCN_J43)
        writer.writeMolLineFile(os.path.join(HCN_J54_dir_name, inputs.out_cwlname_HCN_J54), inputs.vmax, inputs.nwav, inputs.restfreq_HCN_J54)
        writer.writeMolLineFile(os.path.join(HCN_J65_dir_name, inputs.out_cwlname_HCN_J65), inputs.vmax, inputs.nwav, inputs.restfreq_HCN_J65)
        writer.writeMolLineFile(os.path.join(HCN_J76_dir_name, inputs.out_cwlname_HCN_J76), inputs.vmax, inputs.nwav, inputs.restfreq_HCN_J76)
        writer.writeMolLineFile(os.path.join(HCN_J87_dir_name, inputs.out_cwlname_HCN_J87), inputs.vmax, inputs.nwav, inputs.restfreq_HCN_J87)
        writer.writeMolLineFile(os.path.join(HCN_J98_dir_name, inputs.out_cwlname_HCN_J98), inputs.vmax, inputs.nwav, inputs.restfreq_HCN_J98)
        writer.writeMolLineFile(os.path.join(HCN_J109_dir_name, inputs.out_cwlname_HCN_J109), inputs.vmax, inputs.nwav, inputs.restfreq_HCN_J109)

        writer.writeMolLineFile(os.path.join(HNC_J10_dir_name, inputs.out_cwlname_HNC_J10), inputs.vmax, inputs.nwav, inputs.restfreq_HNC_J10)
        writer.writeMolLineFile(os.path.join(HNC_J21_dir_name, inputs.out_cwlname_HNC_J21), inputs.vmax, inputs.nwav, inputs.restfreq_HNC_J21)
        writer.writeMolLineFile(os.path.join(HNC_J32_dir_name, inputs.out_cwlname_HNC_J32), inputs.vmax, inputs.nwav, inputs.restfreq_HNC_J32)
        writer.writeMolLineFile(os.path.join(HNC_J43_dir_name, inputs.out_cwlname_HNC_J43), inputs.vmax, inputs.nwav, inputs.restfreq_HNC_J43)
        writer.writeMolLineFile(os.path.join(HNC_J54_dir_name, inputs.out_cwlname_HNC_J54), inputs.vmax, inputs.nwav, inputs.restfreq_HNC_J54)
        writer.writeMolLineFile(os.path.join(HNC_J65_dir_name, inputs.out_cwlname_HNC_J65), inputs.vmax, inputs.nwav, inputs.restfreq_HNC_J65)
        writer.writeMolLineFile(os.path.join(HNC_J76_dir_name, inputs.out_cwlname_HNC_J76), inputs.vmax, inputs.nwav, inputs.restfreq_HNC_J76)
        writer.writeMolLineFile(os.path.join(HNC_J87_dir_name, inputs.out_cwlname_HNC_J87), inputs.vmax, inputs.nwav, inputs.restfreq_HNC_J87)
        writer.writeMolLineFile(os.path.join(HNC_J98_dir_name, inputs.out_cwlname_HNC_J98), inputs.vmax, inputs.nwav, inputs.restfreq_HNC_J98)
        writer.writeMolLineFile(os.path.join(HNC_J109_dir_name, inputs.out_cwlname_HNC_J109), inputs.vmax, inputs.nwav, inputs.restfreq_HNC_J109)

    else:
        print("Copying exisiting line file ...")
        shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_cwlname), working_dir_name)

    # Copy over existing files
    print('Copying default files...')
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_molname_CO), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_molname_13CO), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_molname_C18O), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_molname_HCOp), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_molname_H13COp), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_molname_HCN), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_molname_HNC), working_dir_name)

    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_linesname_CO), lines_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_linesname_13CO), lines_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_linesname_C18O), lines_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_linesname_HCOp), lines_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_linesname_H13COp), lines_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_linesname_HCN), lines_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_linesname_HNC), lines_dir_name)

    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_wlmname), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_dksname), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_dtpname), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_rmcname), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_execute), working_dir_name)

    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_jobscript_CO), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_jobscript_13CO), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_jobscript_C18O), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_jobscript_HCOp), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_jobscript_H13COp), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_jobscript_HCN), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_jobscript_HNC), working_dir_name)

    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_extra1), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_extra2), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_extra3), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_extra4), working_dir_name)
    shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_extra5), working_dir_name)

    print('Done! Output files generated at: \n\n' + os.path.abspath(working_dir_name))

    return working_dir_name
