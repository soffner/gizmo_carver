import sys
sys.path.append('/g/data1b/jh2/ps3459/gizmo_carver/src/')  # Add the directory to sys.path

import yt
import inputs_gizmo_carver as inputs
import numpy as np
import h5py
import astropy.constants as cons


#f = h5py.File('/g/data1b/jh2/ps3459/starforge_data/15/snapshot_2800.hdf5', 'a') #open in append mode so as to preserve all the other data
f = h5py.File(inputs.hdf5_dir + 'snapshot_'+inputs.snap+'.hdf5', 'a')
print('Opened file ', f)

#delete existing datasets of the same name:
if 'NumDensity' in f['/PartType0']:
	del f['/PartType0/NumDensity']
	print('Deleted existing dataset NumDensity')

if 'HColDensity' in f['/PartType0']:
	del f['/PartType0/HColDensity']
	print('Deleted existing dataset HColDensity')

if 'H2ColDensity' in f['/PartType0']:
	del f['/PartType0/H2ColDensity']
	print('Deleted existing dataset H2ColDensity')

if 'COAbund' in f['/PartType0']:
        del f['/PartType0/COAbund']
        print('Deleted existing dataset COAbund')

if 'COAbundDespotic' in f['/PartType0']:
        del f['/PartType0/COAbundDespotic']
        print('Deleted existing dataset COAbundDespotic')

if 'COAbundUCLCHEM' in f['/PartType0']:
        del f['/PartType0/COAbundUCLCHEM']
        print('Deleted existing dataset COAbundUCLCHEM')

if 'HCOpAbundDespotic' in f['/PartType0']:
           del f['/PartType0/HCOpAbundDespotic']
           print('Deleted existing dataset HCOpAbundDespotic')

if 'HCOpAbundUCLCHEM' in f['/PartType0']:
        del f['/PartType0/HCOpAbundUCLCHEM']
        print('Deleted existing dataset HCOpAbundUCLCHEM')

if 'HCNAbundUCLCHEM' in f['/PartType0']:
           del f['/PartType0/HCNAbundUCLCHEM']
           print('Deleted existing dataset HCNAbundUCLCHEM')

if 'HNCAbundUCLCHEM' in f['/PartType0']:
           del f['/PartType0/HNCAbundUCLCHEM']
           print('Deleted existing dataset HNCAbundUCLCHEM')

#we use the following unit base for starforge simulations
#unit_base = {'UnitMagneticField_in_gauss':  1e+4,
#                 'UnitLength_in_cm'         : cons.pc.cgs.value,
#                 'UnitMass_in_g'            : cons.M_sun.cgs.value,
#                 'UnitVelocity_in_cm_per_s' :      100}

nh_eff = np.loadtxt(inputs.hdf5_dir+'nh_eff_meshoid_10rays_'+inputs.snap+'.txt')
xH2 = f['PartType0']['MolecularMassFraction'][:]*f['PartType0']['NeutralHydrogenAbundance'][:]/inputs.mol_hydrogen_ratio
mean_mass_H = 1.0 + inputs.helium_mass_fraction
mean_mol_mass = mean_mass_H / (1.0 - xH2 + inputs.helium_mass_fraction/4.0)
mH = cons.m_p.cgs.value + cons.m_e.cgs.value
nh2 = nh_eff * ((cons.M_sun.cgs.value/cons.pc.cgs.value**2)/(mH*mean_mol_mass)) #in cm**-2; in new_fields_yt.py, we will create field H2ColumnDensity from this field (so the units are consistent)
f.create_dataset('/PartType0/H2ColDensity', data=nh2)
print('Created new dataset H2ColDensity')

nh = f['PartType0']['Density'][:]*(cons.M_sun.cgs.value/cons.pc.cgs.value**3)/(mH*mean_mol_mass) #in cm**-3; in new_fiels_yt.py we will add units
f.create_dataset('/PartType0/NumDensity', data=nh)
print('Created new dataset NumDensity')

co_desp = np.loadtxt(inputs.hdf5_dir+'co_abund_despotic_'+inputs.snap+'.txt')
f.create_dataset('/PartType0/COAbundDespotic', data=co_desp)
print('Created new dataset COAbundDespotic')


hcop_desp = np.loadtxt(inputs.hdf5_dir+'hcop_abund_despotic_'+inputs.snap+'.txt')
f.create_dataset('/PartType0/HCOpAbundDespotic', data=hcop_desp)
print('Created new dataset HCOpAbundDespotic')


co_ucl = np.loadtxt(inputs.hdf5_dir+'co_abund_uclchem_'+inputs.snap+'.txt')
f.create_dataset('/PartType0/COAbundUCLCHEM', data=co_ucl)
print('Created new dataset COAbundUCLCHEM')


co_ucl = np.loadtxt(inputs.hdf5_dir+'hcop_abund_uclchem_'+inputs.snap+'.txt')
f.create_dataset('/PartType0/HCOpAbundUCLCHEM', data=co_ucl)
print('Created new dataset HCOpAbundUCLCHEM')

co_ucl = np.loadtxt(inputs.hdf5_dir+'hcn_abund_uclchem_'+inputs.snap+'.txt')
f.create_dataset('/PartType0/HCNAbundUCLCHEM', data=co_ucl)
print('Created new dataset HCNAbundUCLCHEM')

co_ucl = np.loadtxt(inputs.hdf5_dir+'hnc_abund_uclchem_'+inputs.snap+'.txt')
f.create_dataset('/PartType0/HNCAbundUCLCHEM', data=co_ucl)
print('Created new dataset HNCAbundUCLCHEM')

f.close()

