import numpy as np
import h5py
import inputs_gizmo_carver as inputs
import astropy.constants as cons


#f = h5py.File('/g/data1b/jh2/ps3459/starforge_data/15/snapshot_2800.hdf5', 'a') #open in append mode so as to preserve all the other data
f = h5py.File(inputs.hdf5_dir + 'snapshot_'+inputs.snap+'.hdf5', 'a')
print('Opened file ', f)

#delete existing datasets of the same name:
if 'HColDensity' in f['/PartType0']:
	del f['/PartType0/HColDensity']
	print('Deleted existing dataset HColDensity')

if 'H2ColDensity' in f['/PartType0']:
	del f['/PartType0/H2ColDensity']
	print('Deleted existing dataset H2ColDensity')

if 'CONumDensity' in f['/PartType0']:
	del f['/PartType0/CONumDensity']
	print('Deleted existing dataset CONumDensity')

if 'CONumDensityDespotic' in f['/PartType0']:
        del f['/PartType0/CONumDensityDespotic']
        print('Deleted existing dataset CONumDensityDespotic')

if 'CONumDensityUCLCHEM' in f['/PartType0']:
        del f['/PartType0/CONumDensityUCLCHEM']
        print('Deleted existing dataset CONumDensityUCLCHEM')

if 'HCOpNumDensityDespotic' in f['/PartType0']:
           del f['/PartType0/HCOpNumDensityDespotic']
           print('Deleted existing dataset HCOpNumDensityDespotic')

if 'HCOpNumDensityUCLCHEM' in f['/PartType0']:
        del f['/PartType0/HCOpNumDensityUCLCHEM']
        print('Deleted existing dataset HCOpNumDensityUCLCHEM')

if 'HCNNumDensityUCLCHEM' in f['/PartType0']:
           del f['/PartType0/HCNNumDensityUCLCHEM']
           print('Deleted existing dataset HCNNumDensityUCLCHEM')


nh_eff = np.loadtxt('/g/data1b/jh2/ps3459/starforge_data/15/nh_eff_meshoid_100rays.txt')
xH2 = f['PartType0']['MolecularMassFraction'][:]*f['PartType0']['NeutralHydrogenAbundance'][:]/inputs.mol_hydrogen_ratio
mean_mass_H = 1.0 + inputs.helium_mass_fraction
mean_mol_mass = mean_mass_H / (1.0 - xH2 + inputs.helium_mass_fraction/4.0)
mH = cons.m_p.cgs.value + cons.m_e.cgs.value
nh2 = nh_eff * ((cons.M_sun.cgs.value/cons.pc.cgs.value**2)/(mH*mean_mol_mass)) #in cm**-2; in new_fields_yt.py, we will create field H2ColumnDensity from this field (so the units are consistent)
f.create_dataset('/PartType0/H2ColDensity', data=nh2)
print('Created new dataset H2ColDensity')

co_desp = np.loadtxt('/g/data1b/jh2/ps3459/starforge_data/15/co_abund_despotic_nov23_new.txt')
f.create_dataset('/PartType0/CONumDensityDespotic', data=co_desp)
print('Created new dataset CONumDensityDespotic')


hcop_desp = np.loadtxt('/g/data1b/jh2/ps3459/starforge_data/15/hcop_abund_despotic_nov23_new.txt')
f.create_dataset('/PartType0/HCOpNumDensityDespotic', data=hcop_desp)
print('Created new dataset HCOpNumDensityDespotic')


co_ucl = np.loadtxt('/g/data1b/jh2/ps3459/starforge_data/15/co_abund_uclchem_nov23_new.txt')
f.create_dataset('/PartType0/CONumDensityUCLCHEM', data=co_ucl)
print('Created new dataset CONumDensityUCLCHEM')


co_ucl = np.loadtxt('/g/data1b/jh2/ps3459/starforge_data/15/hcop_abund_uclchem_nov23_new.txt')
f.create_dataset('/PartType0/HCOpNumDensityUCLCHEM', data=co_ucl)
print('Created new dataset HCOpNumDensityUCLCHEM')

co_ucl = np.loadtxt('/g/data1b/jh2/ps3459/starforge_data/15/hcn_abund_uclchem_nov23_new.txt')
f.create_dataset('/PartType0/HCNNumDensityUCLCHEM', data=co_ucl)
print('Created new dataset HCNNumDensityUCLCHEM')


f.close()

