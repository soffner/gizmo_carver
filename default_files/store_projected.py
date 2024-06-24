import sys
sys.path.append('/g/data1b/jh2/ps3459/gizmo_carver/src/')  # Add the directory to sys.path

import yt
import inputs_gizmo_carver as inputs
from yt.units import *
import numpy as np
from new_fields_list import *

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
from meshoid import Meshoid
import h5py
import cmasher as cm
import astropy.constants as cons

fontsize = 12

hdf5_file = inputs.hdf5_dir+'snapshot_'+inputs.snap+'.hdf5'
F = h5py.File(hdf5_file, 'r')

rho = F["PartType0"]["Density"][:]
pdata = {}
for field in "Masses", "Coordinates", "SmoothingLength", "Velocities":
    pdata[field] = F["PartType0"][field][:]


pos = pdata["Coordinates"]
center = np.median(pos,axis=0) #np.array([15.95957649, 15.54566532, 15.19446488]) #np.median(pos,axis=0)
center[1]=center[1]-0.8125 #to get the dense clump in the center of the box
center[0]=center[0]-0.1875

rmax = inputs.box_size #float(input('Enter box size from center in pc '))
res = int(inputs.box_dim.imag)
X = Y = np.linspace(-rmax, rmax, res)
X, Y = np.meshgrid(X, Y)

pos -= center
pos, mass, hsml, v = pos, pdata["Masses"], pdata["SmoothingLength"], pdata["Velocities"]

M = Meshoid(pos, mass, hsml)

#projected average of number density
xH2 = F['PartType0']['MolecularMassFraction'][:]*F['PartType0']['NeutralHydrogenAbundance'][:]/inputs.mol_hydrogen_ratio #we need per H nuclei, so divide$
mean_mass_H = 1.0 + inputs.helium_mass_fraction
mean_mol_mass = mean_mass_H / (1.0 - xH2 + inputs.helium_mass_fraction/4.0) #para above section 2.2 of Sharda & Krumholz 2022
dens = F['PartType0']['Density'][:]
nh = dens*(cons.M_sun.cgs.value/cons.pc.cgs.value**3)/(mean_mol_mass*cons.m_p.cgs.value)
bb = M.ProjectedAverage(nh, center=np.array([0,0,0]), res=res, size=2*rmax)
np.savetxt(inputs.hdf5_dir + 'nh_proj_' + inputs.snap + '_' + format(rmax, '.1f') + '_' + str(res) + '.txt', bb, fmt='%0.1e')

#projected average of FUV ISRF
vol = mass/dens
isrf = F['PartType0']['PhotonEnergy'][:,1]/vol #photon energy is in msun (m/s)^2
isrf_norm = isrf*(cons.M_sun.cgs.value*1e4/cons.pc.cgs.value**3)/6.63e-14
bb = M.ProjectedAverage(isrf_norm, center=np.array([0,0,0]), res=res, size=2*rmax)
np.savetxt(inputs.hdf5_dir + 'isrf_proj_' + inputs.snap + '_' + format(rmax, '.1f') + '_' + str(res) + '.txt', bb, fmt='%0.1e')

#projected average of H2 column density - note this is the column locally seen by a particle (which we find from pytreegrav)
dens = F['PartType0']['H2ColDensity'][:]
bb = M.ProjectedAverage(dens, center=np.array([0,0,0]), res=res, size=2*rmax)
np.savetxt(inputs.hdf5_dir + 'pytreegrav_nh2_proj_' + inputs.snap + '_' + format(rmax, '.1f') + '_' + str(res) + '.txt', bb, fmt='%0.1e')

#projected average of gas temperature
dens = F['PartType0']['Temperature'][:]
bb = M.ProjectedAverage(dens, center=np.array([0,0,0]), res=res, size=2*rmax)
np.savetxt(inputs.hdf5_dir + 'gastemp_proj_' + inputs.snap + '_' + format(rmax, '.1f') + '_' + str(res) + '.txt', bb, fmt='%0.1e')

#projected average of dust temperature
dens = F['PartType0']['Dust_Temperature'][:]
bb = M.ProjectedAverage(dens, center=np.array([0,0,0]), res=res, size=2*rmax)
np.savetxt(inputs.hdf5_dir + 'dusttemp_proj_' + inputs.snap + '_' + format(rmax, '.1f') + '_' + str(res) + '.txt', bb, fmt='%0.1e')

print('All done!')
