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

fontsize = 12

def set_allticks(ax):
    ax.tick_params(axis='both', which='major', direction = 'in', top=True, right=True, labelsize=fontsize, length=7)
    ax.minorticks_on()
    ax.tick_params(axis='both', which='minor', direction = 'in', top=True, right=True, labelsize=fontsize, length=4)
    return None

def set_cbar(image,axis,cbar_label):
    #sets all the properties of the colorbar for any figure
    divider = make_axes_locatable(axis)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    cb=f.colorbar(image, cax=cax, orientation='vertical')
    cb.set_label(cbar_label, fontsize=fontsize)
    cb.ax.tick_params(which='major',direction='in',labelsize=fontsize,length=5)
    cb.ax.tick_params(which='minor',direction='in',labelsize=fontsize,length=3)
    return None



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

rmax = 20.0 #float(input('Enter box size from center in pc '))
res = 256
X = Y = np.linspace(-rmax, rmax, res)
X, Y = np.meshgrid(X, Y)

pos -= center
pos, mass, hsml, v = pos, pdata["Masses"], pdata["SmoothingLength"], pdata["Velocities"]

M = Meshoid(pos, mass, hsml)

f, ax = plt.subplots(figsize=(6,6))
sigma_gas_msun_pc2 = M.SurfaceDensity(M.m, center=np.array([0,0,0]),
                                      res=res, size=rmax)

np.savetxt(inputs.hdf5_dir + 'Sigma_gas_' + inputs.snap + '_' + str(rmax) + 'pc.txt', sigma_gas_msun_pc2, fmt='%0.1e')
p = ax.pcolormesh(X, Y, sigma_gas_msun_pc2, norm=colors.LogNorm(vmin=.1,vmax=1e3))
ax.set_aspect('equal')
set_cbar(p, ax, r"$\Sigma_{gas}$ $(\rm M_\odot\,pc^{-2})$")
xlims = ax.get_xlim()
ylims = ax.get_ylim()

#plot stars
if 'PartType5' in F.keys():
    masses = F["PartType5"]["Masses"][:]
    print('Number of stars formed ', len(masses))
    print('Min mean max stellar mass %0.3f %0.3f %0.3f'%(np.min(masses), np.mean(masses), np.max(masses)), ' Msun')
    sink_pos = F["PartType5"]["Coordinates"][:]
    sink_pos -= center
    sink_pos_x = sink_pos[:,0]
    sink_pos_y = sink_pos[:,1]

    massive = np.where(masses >= 8)[0]
    print('Massive stars are ', masses[massive], ' Msun')

    ax.scatter(sink_pos_x[massive], sink_pos_y[massive], marker='*', edgecolor='m', fc='white', s=100, alpha=0.5)
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
else:
    print('No stars formed yet.')

ax.set_xlabel("x (pc)")
ax.set_ylabel("y (pc)")
ax.grid()
f.savefig(inputs.hdf5_dir + 'Sigma_gas_' + inputs.snap + '_' + str(rmax) + 'pc.png', bbox_inches='tight')
plt.clf()
