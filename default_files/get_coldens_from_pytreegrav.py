import h5py
from meshoid import Meshoid
import pytreegrav

F = h5py.File('snapshot_000.hdf5', 'r')

rho = F["PartType0"]["Density"][:]
pdata = {}
for field in "Masses", "Coordinates", "SmoothingLength", "Velocities":
    pdata[field] = F["PartType0"][field][:]

pos = pdata["Coordinates"]
center = np.median(pos,axis=0)
pos -= center
pos, mass, hsml, v = pos, pdata["Masses"], pdata["SmoothingLength"], pdata["Velocities"]

vol = mass/rho
l_eff = (3 * vol / (4 * 3.14159))**(1./3.) #we model each gas particle as a sphere. This is also what pytreegrav models them as while doing RT.
kappa = 0.02 #opacity in code units. Accurate opacity really matters only if the column is very anisotropic.
sigma = mass*kappa

bb = pytreegrav.ColumnDensity(pos, sigma, l_eff, parallel=True, randomize_rays=True, rays=100)
bb_eff = -np.log(np.exp(-bb.clip(-300,300)).mean(axis=1)) # effective optical depth that would give the same radiation flux from a background; note clipping because overflow is not uncommon here
NH_eff = bb_eff / kappa # effective column density *for this opacity* in code mass/code length^2


#save the column densities in a txt file
np.savetxt('nh_eff_100rays.txt', NH_eff, fmt='%e')
