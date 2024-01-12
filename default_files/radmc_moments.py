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
hcopstr = 'restfreq_HCOp_J'
h13copstr = 'restfreq_H13COp_J'
hcnstr = 'restfreq_HCN_J'
hncstr = 'restfreq_HNC_J'

Js = ['10', '21', '32'] #, '43', '54', '65', '76', '87', '98', '109']
moment = [0, 1, 2]

for j in range(0, len(moment)):

    for i in range(0, len(Js)):
        #CO
        freq = getattr(inputs, costr + Js[i])
        print('frequency used over CO 1-0 frequency is ', freq/inputs.restfreq_CO_J10)

        m_image = radmc3dImage()
        m_image.readImage('default/image_CO_default_J'+Js[i]+'.out')
        tb = m_image.compute_brightness_temperature(linear=True)
        print('Min max Tb for default CO ' + Js[i] + 'is ', np.min(tb), np.max(tb), ' K')
        co_mom0 = m_image.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('default/default_co_mom_'+str(moment[j])+'_J'+Js[i]+'.txt', co_mom0, fmt='%e')

        
        m_image = radmc3dImage()
        m_image.readImage('despotic/image_CO_despotic_J'+Js[i]+'.out')
        tb = m_image.compute_brightness_temperature(linear=True)
        print('Min max Tb for despotic CO ' + Js[i] + 'is ', np.min(tb), np.max(tb), ' K')
        co_mom0 = m_image.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('despotic/despotic_co_mom_'+str(moment[j])+'_J'+Js[i]+'.txt', co_mom0, fmt='%e')

        m_image = radmc3dImage()
        m_image.readImage('uclchem/image_CO_uclchem_J'+Js[i]+'.out')
        tb = m_image.compute_brightness_temperature(linear=True)
        print('Min max Tb for uclchem CO ' + Js[i] + 'is ', np.min(tb), np.max(tb), ' K')
        co_mom0 = m_image.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('uclchem/uclchem_co_mom_'+str(moment[j])+'_J'+Js[i]+'.txt', co_mom0, fmt='%e')
        
        #13CO
        freq = getattr(inputs, co13str + Js[i])
        print('frequency used over 13CO 1-0 frequency is ', freq/inputs.restfreq_13CO_J10)

        m_image = radmc3dImage()
        m_image.readImage('despotic/image_13CO_despotic_J'+Js[i]+'.out')
        tb = m_image.compute_brightness_temperature(linear=True)
        print('Min max Tb for despotic 13CO ' + Js[i] + 'is ', np.min(tb), np.max(tb), ' K')
        co_mom0 = m_image.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('despotic/despotic_13co_mom_'+str(moment[j])+'_J'+Js[i]+'.txt', co_mom0, fmt='%e')

        m_image = radmc3dImage()
        m_image.readImage('uclchem/image_13CO_uclchem_J'+Js[i]+'.out')
        tb = m_image.compute_brightness_temperature(linear=True)
        print('Min max Tb for uclchem 13CO ' + Js[i] + 'is ', np.min(tb), np.max(tb), ' K')
        co_mom0 = m_image.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('uclchem/uclchem_13co_mom_'+str(moment[j])+'_J'+Js[i]+'.txt', co_mom0, fmt='%e')

        
        #HCOp
        freq = getattr(inputs, hcopstr + Js[i])
        print('frequency used over HCO+ 1-0 frequency is ', freq/inputs.restfreq_HCOp_J10)

        m_image = radmc3dImage()
        m_image.readImage('despotic/image_HCOp_despotic_J'+Js[i]+'.out')
        tb = m_image.compute_brightness_temperature(linear=True)
        print('Min max Tb for despotic HCOp ' + Js[i] + 'is ', np.min(tb), np.max(tb), ' K')
        co_mom0 = m_image.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('despotic/despotic_hcop_mom_'+str(moment[j])+'_J'+Js[i]+'.txt', co_mom0, fmt='%e')

        m_image = radmc3dImage()
        m_image.readImage('uclchem/image_HCOp_uclchem_J'+Js[i]+'.out')
        tb = m_image.compute_brightness_temperature(linear=True)
        print('Min max Tb for uclchem HCOp ' + Js[i] + 'is ', np.min(tb), np.max(tb), ' K')
        co_mom0 = m_image.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('uclchem/uclchem_hcop_mom_'+str(moment[j])+'_J'+Js[i]+'.txt', co_mom0, fmt='%e')
        
        #H13COp
        freq = getattr(inputs, h13copstr + Js[i])
        print('frequency used over H13CO+ 1-0 frequency is ', freq/inputs.restfreq_H13COp_J10)

        m_image = radmc3dImage()
        m_image.readImage('despotic/image_H13COp_despotic_J'+Js[i]+'.out')
        tb = m_image.compute_brightness_temperature(linear=True)
        print('Min max Tb for despotic H13COp ' + Js[i] + 'is ', np.min(tb), np.max(tb), ' K')
        co_mom0 = m_image.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('despotic/despotic_h13cop_mom_'+str(moment[j])+'_J'+Js[i]+'.txt', co_mom0, fmt='%e')

        m_image = radmc3dImage()
        m_image.readImage('uclchem/image_H13COp_uclchem_J'+Js[i]+'.out')
        tb = m_image.compute_brightness_temperature(linear=True)
        print('Min max Tb for uclchem H13COp ' + Js[i] + 'is ', np.min(tb), np.max(tb), ' K')
        co_mom0 = m_image.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('uclchem/uclchem_h13cop_mom_'+str(moment[j])+'_J'+Js[i]+'.txt', co_mom0, fmt='%e')

        '''
        #HCN
        freq = getattr(inputs, hcnstr + Js[i])
        print('frequency used over HCN 1-0 frequency is ', freq/inputs.restfreq_HCN_J10)

        m_image = radmc3dImage()
        m_image.readImage('uclchem/image_HCN_uclchem_J'+Js[i]+'.out')
        tb = m_image.compute_brightness_temperature(linear=True)
        print('Min max Tb for uclchem HCN ' + Js[i] + 'is ', np.min(tb), np.max(tb), ' K')
        co_mom0 = m_image.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('uclchem/uclchem_hcn_mom_'+str(moment[j])+'_J'+Js[i]+'.txt', co_mom0, fmt='%e')
        
        #HNC
        freq = getattr(inputs, hncstr + Js[i])
        print('frequency used over HNC 1-0 frequency is ', freq/inputs.restfreq_HNC_J10)

        m_image = radmc3dImage()
        m_image.readImage('uclchem/image_HNC_uclchem_J'+Js[i]+'.out')
        tb = m_image.compute_brightness_temperature(linear=True)
        print('Min max Tb for uclchem HNC ' + Js[i] + 'is ', np.min(tb), np.max(tb), ' K')
        co_mom0 = m_image.getMomentMap(moment=moment[j], nu0=freq, Tb=True, linear=True)
        np.savetxt('uclchem/uclchem_hnc_mom_'+str(moment[j])+'_J'+Js[i]+'.txt', co_mom0, fmt='%e')
        '''

        print('')