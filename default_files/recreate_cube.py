import sys
sys.path.append('/g/data1b/jh2/ps3459/gizmo_carver/src/')  # Add the directory to sys.path

import yt
import inputs_gizmo_carver as inputs
from yt.units import *
import numpy as np
from new_fields_list import *

import matplotlib.pyplot as plt

hdf5_file = inputs.hdf5_dir+'snapshot_'+inputs.snap+'.hdf5'
ds = yt.load(hdf5_file, unit_base=inputs.unit_base)

print('using box size ', inputs.box_size, ' pc on either side of the center')

le = [inputs.box_center[0] - inputs.box_size, inputs.box_center[1] - inputs.box_size, inputs.box_center[2] - inputs.box_size]
re = [inputs.box_center[0] + inputs.box_size, inputs.box_center[1] + inputs.box_size, inputs.box_center[2] + inputs.box_size]
res = inputs.box_dim
ires = int(res.imag)

cs = ds.r[le[0]:re[0]:res, le[1]:re[1]:res, le[2]:re[2]:res]


numh2 = cs[('PartType0', 'H2NumDensity')]
cc = np.zeros((ires, ires))
#integrate along the z direction to find the column density
for i in range(0, ires):
	for j in range(0, ires):
		cc[i][j] = (2.0*inputs.box_size/res.imag)*3.086e18*np.sum(numh2[i][j])*yt.YTQuantity(1, 'cm**-2')

np.savetxt(inputs.hdf5_dir +'numh2_'+inputs.snap+'_'+format(inputs.box_size, '.1f')+'_'+str(ires)+'.txt', cc, fmt='%e')

#f, ax = plt.subplots(1,1)
#im = ax.imshow(np.log10(cc), origin='lower')
#f.colorbar(im)
#f.savefig('coldens_'+inputs.snap+'_'+str(inputs.box_size)+'_'+str(ires)+'.png', bbox_inches='tight')

numh2 = cs[('PartType0', 'CONumberDensityDespotic')]
cc = np.zeros((ires, ires))
#integrate along the z direction to find the column density
for i in range(0, ires):
           for j in range(0, ires):
                   cc[i][j] = (2.0*inputs.box_size/res.imag)*3.086e18*np.sum(numh2[i][j])*yt.YTQuantity(1, 'cm**-2')

np.savetxt(inputs.hdf5_dir +'numco_despotic_'+inputs.snap+'_'+format(inputs.box_size, '.1f')+'_'+str(ires)+'.txt', cc, fmt='%e')

numh2 = cs[('PartType0', 'CONumberDensityUCLCHEM')]
cc = np.zeros((ires, ires))
#integrate along the z direction to find the column density
for i in range(0, ires):
           for j in range(0, ires):
                   cc[i][j] = (2.0*inputs.box_size/res.imag)*3.086e18*np.sum(numh2[i][j])*yt.YTQuantity(1, 'cm**-2')

np.savetxt(inputs.hdf5_dir +'numco_uclchem_'+inputs.snap+'_'+format(inputs.box_size, '.1f')+'_'+str(ires)+'.txt', cc, fmt='%e')


numh2 = cs[('PartType0', 'HCNNumberDensityUCLCHEM')]
cc = np.zeros((ires, ires))
#integrate along the z direction to find the column density
for i in range(0, ires):
           for j in range(0, ires):
                   cc[i][j] = (2.0*inputs.box_size/res.imag)*3.086e18*np.sum(numh2[i][j])*yt.YTQuantity(1, 'cm**-2')

np.savetxt(inputs.hdf5_dir +'numhcn_uclchem_'+inputs.snap+'_'+format(inputs.box_size, '.1f')+'_'+str(ires)+'.txt', cc, fmt='%e')

numh2 = cs[('PartType0', 'HNCNumberDensityUCLCHEM')]
cc = np.zeros((ires, ires))
#integrate along the z direction to find the column density
for i in range(0, ires):
           for j in range(0, ires):
                   cc[i][j] = (2.0*inputs.box_size/res.imag)*3.086e18*np.sum(numh2[i][j])*yt.YTQuantity(1, 'cm**-2')

np.savetxt(inputs.hdf5_dir +'numhnc_uclchem_'+inputs.snap+'_'+format(inputs.box_size, '.1f')+'_'+str(ires)+'.txt', cc, fmt='%e')
