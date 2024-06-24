#0.1Solar case
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
from matplotlib.ticker import ScalarFormatter, MultipleLocator, LogLocator, Locator, NullFormatter
from matplotlib.transforms import Bbox
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Ellipse
from matplotlib.colors import Normalize
from matplotlib import colors

import pandas as pd
from pylab import genfromtxt
import pylab as P

import math
import itertools
import scipy
from scipy.interpolate import RegularGridInterpolator, griddata, NearestNDInterpolator


import sys
sys.path.append('/g/data1b/jh2/ps3459/despotic/')
from despotic import cloud, chemistry

import astropy.constants as cons
from astropy import units as u

import time
import yt
from joblib import Parallel, delayed, cpu_count
from meshoid import Meshoid
import h5py

import warnings
warnings.filterwarnings('ignore')

ncpus=48
path = '/g/data1b/jh2/ps3459/starforge_data/67/'
snap = '1400'
ucl_path = '/g/data1b/jh2/ps3459/UCLCHEM/mydata/'
despcloud_fname = '/g/data1b/jh2/ps3459/despotic/cloudfiles/basic_starforge.desp'
nheff_meshoid = np.loadtxt(path+'nh_eff_meshoid_10rays_'+snap+'.txt')

unit_base = {'UnitMagneticField_in_gauss':  1e+4,
             'UnitLength_in_cm'         : cons.pc.cgs.value,
             'UnitMass_in_g'            : cons.M_sun.cgs.value,
             'UnitVelocity_in_cm_per_s' :      100}


ds = yt.load(path+'snapshot_'+snap+'.hdf5', unit_base=unit_base)

def _NumberDensity(field, data):
    xH2 = data['PartType0', 'MolecularMassFraction']*data['PartType0', 'NeutralHydrogenAbundance']/2.0 #we need per H nuclei, so divide$
    mean_mass_H = 1.0 + 4*0.071
    mean_mol_mass = mean_mass_H / (1.0 - xH2 + 0.071) #para above section 2.2 of Sharda & Krumholz 2022
    nh = data['PartType0', 'Density']/(cons.m_p.cgs.value*mean_mol_mass*yt.YTQuantity(1, 'g'))
    return nh

def _H2ColumnDensity(field, data):
    bb = nheff_meshoid
    xH2 = data['PartType0', 'MolecularMassFraction']*data['PartType0', 'NeutralHydrogenAbundance']/2.0 #we need per H nuclei, so divide$
    mean_mass_H = 1.0 + 4*0.071
    mean_mol_mass = mean_mass_H / (1.0 - xH2 + 0.071) #para above section 2.2 of Sharda & Krumholz 2022
    nh = bb * ((cons.M_sun.cgs.value/cons.pc.cgs.value**2)/(cons.m_p.cgs.value*mean_mol_mass)) * \
         yt.YTQuantity(1, 'cm**-2')
    return nh

def _HColumnDensity(field, data):
    bb = nheff_meshoid
    nh = bb * ((cons.M_sun.cgs.value/cons.pc.cgs.value**2)/cons.m_p.cgs.value) * yt.YTQuantity(1, 'cm**-2')
    return nh

ds.add_field(('PartType0','NumberDensity'), function=_NumberDensity, units='cm**-3', 
             sampling_type='particle')
ds.add_field(('PartType0','H2ColumnDensity'), function=_H2ColumnDensity, units='cm**-2', 
             sampling_type='particle')
ds.add_field(('PartType0','HColumnDensity'), function=_HColumnDensity, units='cm**-2', 
             sampling_type='particle')

dd=ds.all_data()


vol = (dd['PartType0', 'Masses'].to('g') / dd['PartType0', 'Density'].to('g/cm**3'))

#photon energy is in msun (m/s)^2, but yt reads it in as a dimensionless value.
#isrf is photon energy divided by volume, expressed in units erg/cm^3.
isrf = (dd['PartType0', 'PhotonEnergy'][:,1])*yt.YTQuantity(cons.M_sun.cgs.value, 'g') * \
        yt.YTQuantity(1e2*1e2, 'cm/s')/vol
isrf_norm = isrf / (6.63e-14*yt.YTQuantity(1, 'g/cm**2/s')) #normalized to Solar from Habing

#effective length (needed to find approximate column density)
l_eff = (3 * vol / (4 * np.pi))**(1./3.)


#set initial abundances of C and O atoms per H nuclei (adopted from UCLCHEM default)
#!!VARY THIS DEPENDING ON THE METALLICITY!
info = {'xC':1.77e-5, 'xO':3.34e-5, 'xC+':0, 'xCO':0, 'xO+':0} #0.1 Solar


#import basic despotic cloud
c = cloud(fileName=despcloud_fname, verbose=False)
mean_mass_H = 1.0 + 4*c.comp.xHe #this is fixed for all clouds
nH = dd['PartType0', 'Density'].to('g/cm**3')/(mean_mass_H * yt.YTQuantity(cons.m_p.cgs.value, 'g')) #number density per H nuclei
xH2 = dd['PartType0', 'MolecularMassFraction']*dd['PartType0', 'NeutralHydrogenAbundance']/2.0 #we need per H nuclei, so divide by 2 

def get_interpolate(nh, NH, Tg, isrf, xH2, network=chemistry.NL99_GC, maxTime=3.154e12, info=info):
    c.nH = nh
    c.colDen = NH
        
    c.Tg = Tg
    c.Td = 27.5 # mean of 1000 randomly sampled points in starforge 2800 snapshot
    c.T_radDust = 30.0 # again from starforge 2800 snapshot
    c.chi = isrf
    c.comp.xH2 = xH2
    c.comp.xHplus = 1e-5
    c.comp.xHI = 1.0 - c.comp.xH2 - c.comp.xHplus
    
    init_xH2 = c.comp.xH2
    init_xHI = c.comp.xHI
    init_xHplus = c.comp.xHplus

    mean_mol_mass = mean_mass_H / (1.0 - c.comp.xH2 + c.comp.xHe) #para above section 2.2 of Sharda & Krumholz 2022
    c.sigmaNT = np.sqrt(cons.k_B.cgs.value * Tg / (mean_mol_mass * cons.m_p.cgs.value)) #1e5 #set to thermal dispersion. (lower limit). upper limit is bulk velocity dispersion of the gas

    try:
        #old way was to use c.setChemEq. Now we use c.chemEvol
        #chemEq = c.setChemEq(network=network, info=info, maxTime=maxTime, tEqGuess=1e11)
        #abund_CO = 0.5 * c.chemabundances['CO']   #this is per H nuclei, so divide by 2.0 if you want CO/H2
        tt, abund = c.chemEvol(maxTime, network=network, info=info, evolveTemp='fixed')
        if len(tt) != 101:
            raise ValueError('Something went wrong while doing c.chemEvol')
        abund_CO = abund['CO'][len(abund['CO'])-1]
        abund_HCOp = abund['HCO+'][len(abund['HCO+'])-1]
        
    except Exception as e:
        print('Despotic returned error: ', e)
        #chemEq = False
        abund_CO = np.NaN
        abund_HCOp = np.NaN

    return (abund_CO, abund_HCOp)

def clamp(value, minimum, maximum):
    return min(max(value, minimum), maximum)


# Modify the interpolate_5d function to accept a 2D array of points
def interpolate_5d(x1, x2, x3, x4, x5, data, points_set, nan_replacement=1e-30):

    interpolated_results = []

    #if some data is NAN, replace it by nan_replacement
    data = np.where(np.isnan(data), nan_replacement, data)
    
    for points in points_set:
        # Clamp each input coordinate to its respective bounds
        clamped_points = np.array([
            [clamp(points[0], x1.min(), x1.max()),
             clamp(points[1], x2.min(), x2.max()),
             clamp(points[2], x3.min(), x3.max()),
             clamp(points[3], x4.min(), x4.max()),
             clamp(points[4], x5.min(), x5.max())]
        ])

        # Create 1D interpolators for each dimension
        interpolators = [RegularGridInterpolator((x1, x2, x3, x4, x5), data, method='linear')]

        def interpolation_recursive(coordinates, interpolators):
            if len(interpolators) == 1:
                return interpolators[0](coordinates)
            else:
                return interpolation_recursive(coordinates, interpolators[:-1])(coordinates)

        interpolated_result = interpolation_recursive(clamped_points, interpolators)
        interpolated_results.append(interpolated_result)

    return np.array(interpolated_results)


def interpolate_5d_parallel(x1, x2, x3, x4, x5, data, points_set, nan_replacement=1e-30, ncpus=-1):
    def interpolate_single(points, x1, x2, x3, x4, x5, data, nan_replacement):
        # Clamp each input coordinate to its respective bounds
        clamped_points = np.array([
            [clamp(points[0], x1.min(), x1.max()),
             clamp(points[1], x2.min(), x2.max()),
             clamp(points[2], x3.min(), x3.max()),
             clamp(points[3], x4.min(), x4.max()),
             clamp(points[4], x5.min(), x5.max())]
        ])

        # Create 1D interpolators for each dimension
        interpolators = [RegularGridInterpolator((x1, x2, x3, x4, x5), data, method='linear')]

        def interpolation_recursive(coordinates, interpolators):
            if len(interpolators) == 1:
                return interpolators[0](coordinates)
            else:
                return interpolation_recursive(coordinates, interpolators[:-1])(coordinates)

        # Apply interpolation
        interpolated_result = interpolation_recursive(clamped_points, interpolators)
        return interpolated_result

    # Use joblib to parallelize interpolation for multiple points
    results = Parallel(n_jobs=ncpus)(
        delayed(interpolate_single)(
            points, x1, x2, x3, x4, x5, data, nan_replacement
        ) for points in points_set
    )

    return np.array(results)


#now, actual starforge data to be interpolated across
arr_1 = np.array(nH)
arr_2 = np.array(dd['PartType0', 'HColumnDensity'])
arr_3 = np.array(dd['PartType0', 'Temperature'])
arr_4 = np.array(isrf_norm)
arr_5 = np.array(dd['PartType0', 'MolecularMassFraction']*dd['PartType0', 'NeutralHydrogenAbundance'])/2.0


#preparing for despotic and uclchem analyses
nharray = np.sort(np.concatenate((np.logspace(2, 7, 15), np.logspace(1, 1.7, 3))))
NHarray = np.sort(np.concatenate((np.logspace(19, 23, 10), np.array([3.16e23, 1e24, 3.16e24, 1e25]))))
Tgarray = np.logspace(np.log10(5.0), np.log10(5e3), 4)
isrfarray = np.logspace(np.log10(1.0), np.log10(5e3), 4)
xH2array = np.logspace(-3, np.log10(0.5), 4)

despout_co = np.zeros((len(nharray), len(NHarray), len(Tgarray), len(isrfarray), len(xH2array)))
despout_hcop = np.zeros((len(nharray), len(NHarray), len(Tgarray), len(isrfarray), len(xH2array)))
#uclout_0 of the same shape will be imported from UCLCHEM output


# exclude cells where CO is gonna freeze out on dust grains
selrule1 = np.where(np.logical_and( nH > 1e4, dd['PartType0','Temperature'] < 15))[0]

# exclude cells where CO is negligible
selrule2 = np.where(nH < 1e1)[0]

# join the two regimes
badcells = np.concatenate((selrule1, selrule2))

# create an array of numbers from 0 to length of data
indices = np.arange(len(dd['PartType0','Density']))

goodindices_despotic = np.setdiff1d(indices, badcells) #this will return a sorted list of good indices


start_time = time.time()

for i in range(0, len(nharray)):
    nhi = nharray[i]
    for j in range(0, len(NHarray)):
        NHi = NHarray[j]
        for k in range(0, len(Tgarray)):
            Tgi = Tgarray[k]
            for l in range(0, len(isrfarray)):
                isrfi = isrfarray[l]
                for m in range(0, len(xH2array)):
                    xH2i = xH2array[m]
                    bb = get_interpolate(nhi, NHi, Tgi, isrfi, xH2i)
                    despout_co[i,j,k,l,m] = bb[0]
                    despout_hcop[i,j,k,l,m] = bb[1]
                
    print(i)
    
end_time = time.time()
elapsed_time = end_time - start_time

print("Elapsed time to create despotic grid: ", elapsed_time, "seconds")


#create array of inputs for interpolation
input_interp_despotic = [[arr_1[goodindices_despotic[i]], arr_2[goodindices_despotic[i]], 
                          arr_3[goodindices_despotic[i]], arr_4[goodindices_despotic[i]], 
                          arr_5[goodindices_despotic[i]]] for i in range(len(goodindices_despotic))]

#now interpolate for CO
start_time = time.time()
output_interp_despotic_co = interpolate_5d_parallel(nharray, NHarray, Tgarray, isrfarray, xH2array, despout_co, 
                                                    input_interp_despotic, ncpus=ncpus)
end_time = time.time()

elapsed_time = end_time - start_time

print("Elapsed time to interpolate and find CO abundances from despotic: ", elapsed_time, "seconds")

#now interpolate for HCO+
start_time = time.time()
output_interp_despotic_hcop = interpolate_5d_parallel(nharray, NHarray, Tgarray, isrfarray, xH2array, despout_hcop, 
                                                      input_interp_despotic, ncpus=ncpus)
end_time = time.time()

elapsed_time = end_time - start_time

print("Elapsed time to interpolate and find HCO+ abundances from despotic: ", elapsed_time, "seconds")


#create an array the same size as other simulation data
zig_co = np.zeros(len(dd['PartType0', 'Temperature']))
#fill it with the interpolated values where the indices are 'good'
zig_co[goodindices_despotic] = output_interp_despotic_co[:,0]
#save it
np.savetxt(path+'co_abund_despotic_'+snap+'.txt', zig_co, fmt='%e')

zig_hcop = np.zeros(len(dd['PartType0', 'Temperature']))
#fill it with the interpolated values where the indices are 'good'
zig_hcop[goodindices_despotic] = output_interp_despotic_hcop[:,0]
np.savetxt(path+'hcop_abund_despotic_'+snap+'.txt', zig_hcop, fmt='%e')


#reopen CO, HCO+ abundances and make them 0 if the temperature > max temperature available in
#lamda leiden database for collision rates 
bbco = np.loadtxt(path+'co_abund_despotic_'+snap+'.txt')
badtemp = np.where(dd['PartType0', 'Temperature'] >= 2e3)[0]
print(len(badtemp))
bbco[badtemp] = 0
np.savetxt(path+'co_abund_despotic_'+snap+'.txt', bbco, fmt='%e')

bbco = np.loadtxt(path+'hcop_abund_despotic_'+snap+'.txt')
badtemp = np.where(dd['PartType0', 'Temperature'] >= 2e2)[0]
print(len(badtemp))
bbco[badtemp] = 0
np.savetxt(path+'hcop_abund_despotic_'+snap+'.txt', bbco, fmt='%e')


#preparing for uclchem analyses
NHarray_new = np.sort(np.concatenate((NHarray, np.array([2.1e21, 2.2e21, 2.3e21]))))
Tgarray_new = np.sort(np.concatenate((Tgarray, np.array([10, 15, 20, 25, 35, 100, 300, 1000, 2500]))))
isrfarray_new = np.sort(np.concatenate((isrfarray, np.array([30, 80]))))


uclout_co   = np.zeros((len(nharray), len(NHarray_new), len(Tgarray_new), len(isrfarray_new), len(xH2array)))
uclout_hcop = np.zeros((len(nharray), len(NHarray_new), len(Tgarray_new), len(isrfarray_new), len(xH2array)))
uclout_hcn  = np.zeros((len(nharray), len(NHarray_new), len(Tgarray_new), len(isrfarray_new), len(xH2array)))
uclout_hnc  = np.zeros((len(nharray), len(NHarray_new), len(Tgarray_new), len(isrfarray_new), len(xH2array)))

#import uclchem output
uclout_co = np.loadtxt(ucl_path+'uclout_5D_CO_0.1solar.txt', delimiter=',')
uclout_hcop = np.loadtxt(ucl_path+'uclout_5D_HCOp_0.1solar.txt', delimiter=',')
uclout_hcn = np.loadtxt(ucl_path+'uclout_5D_HCN_0.1solar.txt', delimiter=',')
uclout_hnc = np.loadtxt(ucl_path+'uclout_5D_HNC_0.1solar.txt', delimiter=',')

uclout_co   = uclout_co.reshape((len(nharray), len(NHarray_new), len(Tgarray_new), len(isrfarray_new), 
                                 len(xH2array)))
uclout_hcop = uclout_hcop.reshape((len(nharray), len(NHarray_new), len(Tgarray_new), len(isrfarray_new), 
                                 len(xH2array)))
uclout_hcn  = uclout_hcn.reshape((len(nharray), len(NHarray_new), len(Tgarray_new), len(isrfarray_new), 
                                  len(xH2array)))
uclout_hnc  = uclout_hnc.reshape((len(nharray), len(NHarray_new), len(Tgarray_new), len(isrfarray_new),
                                  len(xH2array)))

# create an array of numbers from 0 to length of data
indices = np.arange(len(dd['PartType0','Density']))

# exclude cells where CO is negligible
goodindices_uclchem = np.setdiff1d(indices, selrule2) #this will return a sorted list of good indices

input_interp_uclchem = [[arr_1[goodindices_uclchem[i]], arr_2[goodindices_uclchem[i]], 
                         arr_3[goodindices_uclchem[i]], arr_4[goodindices_uclchem[i]], 
                         arr_5[goodindices_uclchem[i]]] for i in range(len(goodindices_uclchem))]


#now interpolate
start_time = time.time()
output_interp_uclchem_co = interpolate_5d_parallel(nharray, NHarray_new, Tgarray_new, isrfarray_new, xH2array, 
                                                   uclout_co, input_interp_uclchem, ncpus=ncpus)
end_time = time.time()

elapsed_time = end_time - start_time

print("Elapsed time to interpolate and find CO abundances from uclchem: ", elapsed_time, "seconds")

#now interpolate
start_time = time.time()
output_interp_uclchem_hcop = interpolate_5d_parallel(nharray, NHarray_new, Tgarray_new, isrfarray_new, xH2array, 
                                                     uclout_hcop, input_interp_uclchem, ncpus=ncpus)
end_time = time.time()

elapsed_time = end_time - start_time

print("Elapsed time to interpolate and find HCO+ abundances from uclchem: ", elapsed_time, "seconds")


#now interpolate
start_time = time.time()
output_interp_uclchem_hcn = interpolate_5d_parallel(nharray, NHarray_new, Tgarray_new, isrfarray_new, xH2array, 
                                                    uclout_hcn, input_interp_uclchem, ncpus=ncpus)
end_time = time.time()

elapsed_time = end_time - start_time

print("Elapsed time to interpolate and find HCN abundances from uclchem: ", elapsed_time, "seconds")


#now interpolate
start_time = time.time()
output_interp_uclchem_hnc = interpolate_5d_parallel(nharray, NHarray_new, Tgarray_new, isrfarray_new, xH2array,
                                                    uclout_hnc, input_interp_uclchem, ncpus=ncpus)
end_time = time.time()

elapsed_time = end_time - start_time

print("Elapsed time to interpolate and find HNC abundances from uclchem: ", elapsed_time, "seconds")


#create an array the same size as other simulation data
uzig_co = np.zeros(len(dd['PartType0', 'Temperature']))
#fill it with the interpolated values where the indices are 'good'
uzig_co[goodindices_uclchem] = output_interp_uclchem_co[:,0]
#save it
np.savetxt(path+'co_abund_uclchem_'+snap+'.txt', uzig_co, fmt='%e')

uzig_hcop = np.zeros(len(dd['PartType0', 'Temperature']))
uzig_hcop[goodindices_uclchem] = output_interp_uclchem_hcop[:,0]
np.savetxt(path+'hcop_abund_uclchem_'+snap+'.txt', uzig_hcop, fmt='%e')

uzig_hcn = np.zeros(len(dd['PartType0', 'Temperature']))
uzig_hcn[goodindices_uclchem] = output_interp_uclchem_hcn[:,0]
np.savetxt(path+'hcn_abund_uclchem_'+snap+'.txt', uzig_hcn, fmt='%e')


uzig_hnc = np.zeros(len(dd['PartType0', 'Temperature']))
uzig_hnc[goodindices_uclchem] = output_interp_uclchem_hnc[:,0]
np.savetxt(path+'hnc_abund_uclchem_'+snap+'.txt', uzig_hnc, fmt='%e')

#reopen CO, HCO+, HCN and HNC abundances and make them 0 if the temperature > max temperature available in
#lamda leiden database for collision rates 


bbco = np.loadtxt(path+'co_abund_uclchem_'+snap+'.txt')
badtemp = np.where(dd['PartType0', 'Temperature'] >= 2e3)[0]
print(len(badtemp))
bbco[badtemp] = 0
np.savetxt(path+'co_abund_uclchem_'+snap+'.txt', bbco, fmt='%e')

bbco = np.loadtxt(path+'hcop_abund_uclchem_'+snap+'.txt')
badtemp = np.where(dd['PartType0', 'Temperature'] >= 2e2)[0]
print(len(badtemp))
bbco[badtemp] = 0
np.savetxt(path+'hcop_abund_uclchem_'+snap+'.txt', bbco, fmt='%e')

bbco = np.loadtxt(path+'hcn_abund_uclchem_'+snap+'.txt')
badtemp = np.where(dd['PartType0', 'Temperature'] >= 5e2)[0]
print(len(badtemp))
bbco[badtemp] = 0
np.savetxt(path+'hcn_abund_uclchem_'+snap+'.txt', bbco, fmt='%e')


bbco = np.loadtxt(path+'hnc_abund_uclchem_'+snap+'.txt')
badtemp = np.where(dd['PartType0', 'Temperature'] >= 5e2)[0]
print(len(badtemp))
bbco[badtemp] = 0
np.savetxt(path+'hnc_abund_uclchem_'+snap+'.txt', bbco, fmt='%e')

print('All Done!')
