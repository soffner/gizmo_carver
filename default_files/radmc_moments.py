"""
   radmc_moments.py

   Purpose:
        Contains functions for plotting RADMC-3D output image files using MatPlotLib.
        Edit this file to create custom plotting routines.

   Author:
        Sean Feng, feng.sean01@utexas.edu
        Spring 2022

   Written/Tested with Python 3.9, radmc3dpy 0.30.2
   Fully re-written by Piyush Sharda (2023)
"""

from radmc3dPy import *
import matplotlib.pylab as plb
from radmc3dPy import natconst as nc
import numpy as np
from radmc3dPy.image import radmc3dImage
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os
import sys

sys.path.insert(0, '/g/data1b/jh2/ps3459/gizmo_carver/src/')
import inputs_gizmo_carver as inputs

costr = 'restfreq_CO_J'
co13str = 'restfreq_13CO_J'
c18ostr = 'restfreq_C18O_J'
hcopstr = 'restfreq_HCOp_J'
h13copstr = 'restfreq_H13COp_J'
hcnstr = 'restfreq_HCN_J'
hncstr = 'restfreq_HNC_J'

Js = ['10', '21', '32', '43', '54', '65', '76', '87', '98', '109']
moment = [0] #, 1, 2]

#we calculate two different moment maps: one for X_CO, in K km/s, and another for CO SLED and line ratio r, in Jy/pix * km/s

#CO
for i in range(0, len(Js)):
    freq = getattr(inputs, costr + Js[i])
    print('frequency used over CO 1-0 frequency is ', freq/inputs.restfreq_CO_J10)

    m_image = radmc3dImage()
    m_image.readImage('despotic/image_CO_despotic_J'+Js[i]+'.out')
    #tb = m_image.compute_brightness_temperature(linear=False)
    #print('Min max Tb for despotic CO ' + Js[i] + 'is ', np.min(tb), np.max(tb), ' K')

    #moments in K (km/s): only for the first 5 J transitions because you are not in the RJ tail for higher J transitions
    if i < 5:
        for j in range(0, len(moment)):
            co_mom0 = m_image.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
            np.savetxt('despotic/despotic_co_mom_'+str(moment[j])+'_J'+Js[i]+'_Kkms.txt', co_mom0, fmt='%e')


    #moments in Jy/pix (km/s)
    jypix = m_image.imageJyppix
    #basically, copying the getMomentMap routine here and using imageJyppix to get the moment map
    v_kms = nc.cc * (freq - m_image.freq) / freq / 1e5
    vmap = np.zeros([m_image.nx, m_image.ny, m_image.nfreq])
    for ifreq in range(m_image.nfreq):
        vmap[:, :, ifreq] = v_kms[ifreq]

    for j in range(0, len(moment)):
        y = jypix * vmap**moment
        dum = np.abs(vmap[:, :, 1:] - vmap[:, :, :-1]) * (y[:, :, 1:] + y[:, :, :-1]) * 0.5
        co_mom0 = dum.sum(2)
        np.savetxt('despotic/despotic_co_mom_'+str(moment[j])+'_J'+Js[i]+'_Jykms.txt', co_mom0, fmt='%e')

    #UCLCHEM
    m_image = radmc3dImage()
    m_image.readImage('uclchem/image_CO_uclchem_J'+Js[i]+'.out')

    if i < 5:
        for j in range(0, len(moment)):
            co_mom0 = m_image.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
            np.savetxt('uclchem/uclchem_co_mom_'+str(moment[j])+'_J'+Js[i]+'_Kkms.txt', co_mom0, fmt='%e')

    jypix = m_image.imageJyppix
    #basically, copying the getMomentMap routine here and using imageJyppix to get the moment map
    v_kms = nc.cc * (freq - m_image.freq) / freq / 1e5
    vmap = np.zeros([m_image.nx, m_image.ny, m_image.nfreq])
    for ifreq in range(m_image.nfreq):
        vmap[:, :, ifreq] = v_kms[ifreq]

    for j in range(0, len(moment)):
        y = jypix * vmap**moment
        dum = np.abs(vmap[:, :, 1:] - vmap[:, :, :-1]) * (y[:, :, 1:] + y[:, :, :-1]) * 0.5
        co_mom0 = dum.sum(2)
        np.savetxt('uclchem/uclchem_co_mom_'+str(moment[j])+'_J'+Js[i]+'_Jykms.txt', co_mom0, fmt='%e')


    print('')
