import yt
import inputs_gizmo_carver as inputs
from yt.units import *
import numpy as np
from new_fields_list import *

import matplotlib.pyplot as plt

hdf5_file = inputs.hdf5_dir+'snapshot_'+inputs.snap+'.hdf5'
ds = yt.load(hdf5_file, unit_base=inputs.unit_base)

le = [inputs.box_center[0] - inputs.box_size, inputs.box_center[1] - inputs.box_size, inputs.box_center[2] - inputs.box_size]
re = [inputs.box_center[0] + inputs.box_size, inputs.box_center[1] + inputs.box_size, inputs.box_center[2] + inputs.box_size]
res = inputs.box_dim
ires = int(res.imag)

cs = ds.r[le[0]:re[0]:res, le[1]:re[1]:res, le[2]:re[2]:res]

numh2 = cs[('PartType0', 'H2NumDensity')]

cc = np.zeros((ires, ires))

for i in range(0, ires):
	for j in range(0, ires):
		cc[i][j] = (inputs.box_size/res.imag)*3.086e18*np.sum(numh2[i][j])*yt.YTQuantity(1, 'cm**-2')

f, ax = plt.subplots(1,1)

im = ax.imshow(np.log10(cc), origin='lower')
f.colorbar(im)

f.savefig('coldens.png', bbox_inches='tight')
