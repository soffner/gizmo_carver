"""
   radmc_moments_obsv.py

   Purpose:
        Contains functions for plotting RADMC-3D output image files using MatPlotLib.
        Edit this file to create custom plotting routines.

   Author:
        Piyush Sharda (2024), sharda@strw.leidenuniv.nl

   Written/Tested with Python 3.9, radmc3dpy 0.30.2
"""

from radmc3dPy import *
import matplotlib.pylab as plb
import numpy as np
from radmc3dPy.image import radmc3dImage
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os
import sys

sys.path.insert(0, '/g/data1b/jh2/ps3459/gizmo_carver/src/')
import inputs_gizmo_carver as inputs

#CO
costr = 'restfreq_CO_J'
co13str = 'restfreq_13CO_J'
c18ostr = 'restfreq_C18O_J'
hcopstr = 'restfreq_HCOp_J'
h13copstr = 'restfreq_H13COp_J'
hcnstr = 'restfreq_HCN_J'
hncstr = 'restfreq_HNC_J'

Js = ['10', '21', '32'] #, '43', '54', '65', '76', '87', '98', '109']
fwhm_x = 5 #fwhm in arcsec in x direction
fwhm_y = 5 #fwhm in arcsec in y direction
dpc = 300 #distance to observer in pc
moment = [0, 1, 2]


#CO
for i in range(0, len(Js)):
    freq = getattr(inputs, costr + Js[i])
    print('frequency used over CO 1-0 frequency is ', freq/inputs.restfreq_CO_J10)

    m_image = radmc3dImage()
    m_image.readImage('despotic/image_CO_despotic_J'+Js[i]+'.out')
    bb = m_image.imConv(dpc=dpc, psfType='Gauss', pa=0., fwhm = [fwhm_x, fwhm_y])

    for j in range(0, len(moment)):
        vv = bb.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('despotic/despotic_co_mom_'+str(moment[j])+'_obsv_J'+Js[i]+'.txt', vv, fmt='%e')

    m_image = radmc3dImage()
    m_image.readImage('uclchem/image_CO_uclchem_J'+Js[i]+'.out')
    bb = m_image.imConv(dpc=dpc, psfType='Gauss', pa=0., fwhm = [fwhm_x, fwhm_y])

    for j in range(0, len(moment)):
        vv = bb.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('uclchem/uclchem_co_mom_'+str(moment[j])+'_obsv_J'+Js[i]+'.txt', vv, fmt='%e')

#13CO
for i in range(0, len(Js)):
    freq = getattr(inputs, co13str + Js[i])
    print('frequency used over 13CO 1-0 frequency is ', freq/inputs.restfreq_13CO_J10)

    m_image = radmc3dImage()
    m_image.readImage('despotic/image_13CO_despotic_J'+Js[i]+'.out')
    bb = m_image.imConv(dpc=dpc, psfType='Gauss', pa=0., fwhm = [fwhm_x, fwhm_y])

    for j in range(0, len(moment)):
        vv = bb.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('despotic/despotic_13co_mom_'+str(moment[j])+'_obsv_J'+Js[i]+'.txt', vv, fmt='%e')

    m_image = radmc3dImage()
    m_image.readImage('uclchem/image_13CO_uclchem_J'+Js[i]+'.out')
    bb = m_image.imConv(dpc=dpc, psfType='Gauss', pa=0., fwhm = [fwhm_x, fwhm_y])

    for j in range(0, len(moment)):
        vv = bb.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('uclchem/uclchem_13co_mom_'+str(moment[j])+'_obsv_J'+Js[i]+'.txt', vv, fmt='%e')


#C18O
for i in range(0, len(Js)):
    freq = getattr(inputs, c18ostr + Js[i])
    print('frequency used over C18O 1-0 frequency is ', freq/inputs.restfreq_C18O_J10)

    m_image = radmc3dImage()
    m_image.readImage('despotic/image_C18O_despotic_J'+Js[i]+'.out')
    bb = m_image.imConv(dpc=dpc, psfType='Gauss', pa=0., fwhm = [fwhm_x, fwhm_y])

    for j in range(0, len(moment)):
        vv = bb.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('despotic/despotic_c18o_mom_'+str(moment[j])+'_obsv_J'+Js[i]+'.txt', vv, fmt='%e')

    m_image = radmc3dImage()
    m_image.readImage('uclchem/image_C18O_uclchem_J'+Js[i]+'.out')
    bb = m_image.imConv(dpc=dpc, psfType='Gauss', pa=0., fwhm = [fwhm_x, fwhm_y])

    for j in range(0, len(moment)):
        vv = bb.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('uclchem/uclchem_c18o_mom_'+str(moment[j])+'_obsv_J'+Js[i]+'.txt', vv, fmt='%e')

#HCOp
for i in range(0, len(Js)):
    freq = getattr(inputs, hcopstr + Js[i])
    print('frequency used over HCOp 1-0 frequency is ', freq/inputs.restfreq_HCOp_J10)

    m_image = radmc3dImage()
    m_image.readImage('despotic/image_HCOp_despotic_J'+Js[i]+'.out')
    bb = m_image.imConv(dpc=dpc, psfType='Gauss', pa=0., fwhm = [fwhm_x, fwhm_y])

    for j in range(0, len(moment)):
        vv = bb.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('despotic/despotic_hcop_mom_'+str(moment[j])+'_obsv_J'+Js[i]+'.txt', vv, fmt='%e')

    m_image = radmc3dImage()
    m_image.readImage('uclchem/image_HCOp_uclchem_J'+Js[i]+'.out')
    bb = m_image.imConv(dpc=dpc, psfType='Gauss', pa=0., fwhm = [fwhm_x, fwhm_y])

    for j in range(0, len(moment)):
        vv = bb.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('uclchem/uclchem_hcop_mom_'+str(moment[j])+'_obsv_J'+Js[i]+'.txt', vv, fmt='%e')



#H13COp
for i in range(0, len(Js)):
    freq = getattr(inputs, h13copstr + Js[i])
    print('frequency used over H13COp 1-0 frequency is ', freq/inputs.restfreq_H13COp_J10)

    m_image = radmc3dImage()
    m_image.readImage('despotic/image_H13COp_despotic_J'+Js[i]+'.out')
    bb = m_image.imConv(dpc=dpc, psfType='Gauss', pa=0., fwhm = [fwhm_x, fwhm_y])
    for j in range(0, len(moment)):
        vv = bb.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('despotic/despotic_h13cop_mom_'+str(moment[j])+'_obsv_J'+Js[i]+'.txt', vv, fmt='%e')

    m_image = radmc3dImage()
    m_image.readImage('uclchem/image_H13COp_uclchem_J'+Js[i]+'.out')
    bb = m_image.imConv(dpc=dpc, psfType='Gauss', pa=0., fwhm = [fwhm_x, fwhm_y])
    for j in range(0, len(moment)):
        vv = bb.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('uclchem/uclchem_h13cop_mom_'+str(moment[j])+'_obsv_J'+Js[i]+'.txt', vv, fmt='%e')

        
#HCN
for i in range(0, len(Js)):
    freq = getattr(inputs, hcnstr + Js[i])
    print('frequency used over HCN 1-0 frequency is ', freq/inputs.restfreq_HCN_J10)

    m_image = radmc3dImage()
    m_image.readImage('uclchem/image_HCN_uclchem_J'+Js[i]+'.out')
    bb = m_image.imConv(dpc=dpc, psfType='Gauss', pa=0., fwhm = [fwhm_x, fwhm_y])

    for j in range(0, len(moment)):
        vv = bb.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('uclchem/uclchem_hcn_mom_'+str(moment[j])+'_obsv_J'+Js[i]+'.txt', vv, fmt='%e')


#HNC
for i in range(0, len(Js)):
    freq = getattr(inputs, hncstr + Js[i])
    print('frequency used over HNC 1-0 frequency is ', freq/inputs.restfreq_HNC_J10)

    m_image = radmc3dImage()
    m_image.readImage('uclchem/image_HNC_uclchem_J'+Js[i]+'.out')
    bb = m_image.imConv(dpc=dpc, psfType='Gauss', pa=0., fwhm = [fwhm_x, fwhm_y])
    for j in range(0, len(moment)):
        vv = bb.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('uclchem/uclchem_hnc_mom_'+str(moment[j])+'_obsv_J'+Js[i]+'.txt', vv, fmt='%e')


print('All done!')