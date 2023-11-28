import yt
import inputs_gizmo_carver as inputs
from yt.units import *
import numpy as np

# Definition of the dust density. Uses dust to gas ratio from inputs file
def _DustDensity(field, data):
    return inputs.dust_to_gas * data[('PartType0', 'Density')].to('g/cm**3')

# Definition of the gas temperature. Uses info from inputs
def _gas_temp_CO(field, data):
    #y_helium = inputs.helium_mass_fraction/(4.0*(1-inputs.helium_mass_fraction))
    #mu = (1+4.0*y_helium)/(1+y_helium)

    ## In molecular gas electron abundance < 10^-5
    ##mu = (1+4.0*y_helium)/(1+y_helium+data[(‘PartType0’,‘ElectronAbundance’)]
    ##const = ((1.2*mh.in_mks()*(gamma-1))/kboltz.in_mks())

    #const = (mu*mh.in_mks()*(inputs.gamma-1))/(1e6*kboltz.in_mks())
    ## Note the LTE temperature tables are limited to 1e5K
    ## Non-LTE may require lower temperatures
    #return (data[('PartType0','InternalEnergy')]*const * ( data[('PartType0','InternalEnergy')] < yt.YTArray([300], "K")/const))
    return data[('PartType0', 'Temperature')] * yt.YTQuantity(1, 'K') * ( data[('PartType0','Temperature')] < yt.YTArray([2000], 'K') ) #2000 is the max temp available in CO lamda leiden collision rates file

def _gas_temp_HCOp(field, data):
    return data[('PartType0', 'Temperature')] * yt.YTQuantity(1, 'K') * ( data[('PartType0','Temperature')] < yt.YTArray([200], 'K') ) #200 is the max temp available in HCOp lamda leiden collision rates files

def _gas_temp_HCN(field, data):
    return data[('PartType0', 'Temperature')] * yt.YTQuantity(1, 'K') * ( data[('PartType0','Temperature')] < yt.YTArray([500], 'K') ) #500 is the max temp available in HCN lamda leiden collision rates files

# Definition of the target species field. Uses info from inputs
def _H2NumDensity(field, data):
    return data[('PartType0', 'Density')].to('g/cm**3')*data[('PartType0', 'MolecularMassFraction')]*data[('PartType0', 'NeutralHydrogenAbundance')]*(1-inputs.helium_mass_fraction)/(inputs.mol_hydrogen_ratio*yt.YTQuantity(mh, 'g'))

def _H2ColumnDensity(field, data):
    nh = data['PartType0', 'H2ColDensity'] * yt.YTQuantity(1, 'cm**-2')
    return nh

# Definition of the target species field. Uses info from inputs
#def _CONumberDensityDefault(field, data):
#    return data[('PartType0', 'H2NumDensity')].to('cm**-3')*inputs.molecular_abundance*(data[('PartType0','H2NumDensity')].to('cm**-3') > yt.YTArray([1e2], "cm**-3"))*( data[('PartType0','gas_temperature_CO')] < yt.YTArray([1e2], "K"))

# Definition of the target species field based on despotic chemistry. Uses info from inputs and a 9th degree polynomial (for CO)
# ONLY FOR CO!!!
#def _MolecularNumDensityPriestley(field, data):
    #old way - fit a poly to CO/H2
    #fit = np.array([-3.38719496e-04, 1.93973750e-02,  -4.54656223e-01, 5.82774352, -4.54162905e+01,
    #               2.24230736e+02,  -7.03639281e+02, 1.35711326e+03, -1.46211718e+03, 6.61468397e+02])
    #numdens = np.array(data[('PartType0', 'Density')].to('g/cm**3') / ((1.0 + inputs.helium_mass_fraction) * yt.YTQuantity(mh, 'g')))
    #freeze-out
    #selrule1 = np.where(np.logical_and( numdens > 1e4, data[('PartType0','Temperature')] < yt.YTArray([15.0], "K")))[0]
    #no CO at low densities
    #selrule2 = np.where(numdens < 1e2)[0]
    #badcells = np.concatenate((selrule1, selrule2))
    #numdens = np.log10(numdens)
    #logcoh2_numdens = fit[0]*numdens**9 + fit[1]*numdens**8 + fit[2]*numdens**7 + fit[3]*numdens**6 + \
    #             fit[4]*numdens**5 + fit[5]*numdens**4 + fit[6]*numdens**3 + fit[7]*numdens**2 + \
    #             fit[8]*numdens**1 + fit[9]
    #coh2_numdens = 10**logcoh2_numdens
    #co_numdens = coh2_numdens * data[('PartType0', 'H2NumDensity')].to('cm**-3')
    #co_numdens[badcells] = 0
#    return co_numdens

def _CONumberDensityDespotic(field, data):
    co_numdens = data[('PartType0', 'CONumDensityDespotic')] * yt.YTQuantity(1, 'cm**-3')
    return co_numdens

def _HCOpNumberDensityDespotic(field, data):
    hcop_numdens = data[('PartType0', 'HCOpNumDensityDespotic')] * yt.YTQuantity(1, 'cm**-3')
    return hcop_numdens

def _CONumberDensityUCLCHEM(field, data):
    co_numdens = data[('PartType0', 'CONumDensityUCLCHEM')] * yt.YTQuantity(1, 'cm**-3')
    return co_numdens

def _HCOpNumberDensityUCLCHEM(field, data):
    hcop_numdens = data[('PartType0', 'HCOpNumDensityUCLCHEM')] * yt.YTQuantity(1, 'cm**-3')
    return hcop_numdens

def _HCNNumberDensityUCLCHEM(field, data):
    hcn_numdens = data[('PartType0', 'HCNNumDensityUCLCHEM')] * yt.YTQuantity(1, 'cm**-3')
    return hcn_numdens

# Definition of the microturbulence at each point. Uses info from inputs
def _MicroTurb(field, data):
    turb = data['PartType0', 'velocity_x']
    #use the sizd linewidth relation: lw = 0.72 * (size/1pc)**0.56 km/s; size is of the regridded cell
    microturbulence_speed = 0.72 * ((2*inputs.box_size)/inputs.box_dim.imag)**0.56 * 1e5 #cgs
    print('Microturbulence speed is ', microturbulence_speed/1e5, ' km/s')
    turb[turb>=0] = yt.YTQuantity(microturbulence_speed, "cm/s")
    turb[turb<=0] = yt.YTQuantity(microturbulence_speed, "cm/s")
    return turb

# Definition of the dust temperature. Assumes same as gas temperature
def _DustTemperature(field, data):
    return data[('PartType0', 'Dust_Temperature')] * yt.YTQuantity(1, 'K')

# Create a mask based on accreted particles
# Need to update if want to exclude particles not yet formed
def _Mask(field, data):
    mask = data[('PartType0', 'Density')]*0.0
    ids = np.unique(alllines[:,1])
    # This is the slow way. See fast way below
    for id in ids:
        acclist = np.where(alllines[:,1]== id)[0] # All rows accreted by this star 
        ind = np.where(np.isin(data[('PartType0','ParticleIDs')], alllines[acclist,6])== True)[0]  # Particle ids of the accreted particles
        mask[ind] = 1
    return mask


# Add all the fields

yt.add_field(("PartType0", "gas_temperature_CO"), function=_gas_temp_CO, units="K", sampling_type='particle', force_override=True)
yt.add_field(("PartType0", "gas_temperature_HCOp"), function=_gas_temp_HCOp, units="K", sampling_type='particle', force_override=True)
yt.add_field(("PartType0", "gas_temperature_HCN"), function=_gas_temp_HCN, units="K", sampling_type='particle', force_override=True)

yt.add_field(("PartType0", "dust_temperature"), function=_DustTemperature, units="K", sampling_type='particle', force_override=True)

yt.add_field(('PartType0', 'DustDensity'), function=_DustDensity, units='g/cm**3', sampling_type='particle', force_override=True)

yt.add_field(('PartType0', 'H2ColumnDensity'), function=_H2ColumnDensity, units='cm**-2', sampling_type='particle', force_override=True)

yt.add_field(('PartType0', 'H2NumDensity'), function=_H2NumDensity, units='cm**-3', sampling_type='particle', force_override=True)

#yt.add_field(('PartType0', 'CONumberDensityDefault'), function=_CONumberDensityDefault, units='cm**-3', sampling_type='particle', force_override=True)

yt.add_field(('PartType0', 'CONumberDensityDespotic'), function=_CONumberDensityDespotic, units='cm**-3', sampling_type='particle', force_override=True)

yt.add_field(('PartType0', 'HCOpNumberDensityDespotic'), function=_HCOpNumberDensityDespotic, units='cm**-3', sampling_type='particle', force_override=True)

yt.add_field(('PartType0', 'CONumberDensityUCLCHEM'), function=_CONumberDensityUCLCHEM, units='cm**-3', sampling_type='particle', force_override=True)

yt.add_field(('PartType0', 'HCOpNumberDensityUCLCHEM'), function=_HCOpNumberDensityUCLCHEM, units='cm**-3', sampling_type='particle', force_override=True)

yt.add_field(('PartType0', 'HCNNumberDensityUCLCHEM'), function=_HCNNumberDensityUCLCHEM, units='cm**-3', sampling_type='particle', force_override=True)

yt.add_field(("PartType0", "microturbulence_speed"), function=_MicroTurb, units="cm/s", sampling_type='particle', force_override=True)

