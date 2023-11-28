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
hcopstr = 'restfreq_HCOp_J'
hcnstr = 'restfreq_HCN_J'

Js = ['10', '21', '32', '43', '54', '65', '76', '87', '98', '109']


for i in range(0, len(Js)):
    
    freq = getattr(inputs, costr + Js[i])
    print('frequency used over CO 1-0 frequency is ', freq/inputs.restfreq_CO_J10)

    
    m_image = radmc3dImage()
    m_image.readImage('despotic/image_CO_despotic_J'+Js[i]+'.out')
    co_mom0 = m_image.getMomentMap(0, freq, Tb=True)
    np.savetxt('despotic/despotic_co_mom0_J'+Js[i]+'.txt', co_mom0, fmt='%e')
    

    m_image = radmc3dImage()
    m_image.readImage('uclchem/image_CO_uclchem_J'+Js[i]+'.out')
    co_mom0 = m_image.getMomentMap(0, freq, Tb=True)
    np.savetxt('uclchem/uclchem_co_mom0_J'+Js[i]+'.txt', co_mom0, fmt='%e')
    
    
    freq = getattr(inputs, hcopstr + Js[i])
    print('frequency used over HCO+ 1-0 frequency is ', freq/inputs.restfreq_HCOp_J10)

    m_image = radmc3dImage()
    m_image.readImage('despotic/image_HCOp_despotic_J'+Js[i]+'.out')
    co_mom0 = m_image.getMomentMap(0, freq, Tb=True)
    np.savetxt('despotic/despotic_hcop_mom0_J'+Js[i]+'.txt', co_mom0, fmt='%e')

    m_image = radmc3dImage()
    m_image.readImage('uclchem/image_HCOp_uclchem_J'+Js[i]+'.out')
    co_mom0 = m_image.getMomentMap(0, freq, Tb=True)
    np.savetxt('uclchem/uclchem_hcop_mom0_J'+Js[i]+'.txt', co_mom0, fmt='%e')
    

    freq = getattr(inputs, hcnstr + Js[i])
    print('frequency used over HCN 1-0 frequency is ', freq/inputs.restfreq_HCN_J10)

    m_image = radmc3dImage()
    m_image.readImage('uclchem/image_HCN_uclchem_J'+Js[i]+'.out')
    co_mom0 = m_image.getMomentMap(0, freq, Tb=True)
    np.savetxt('uclchem/uclchem_hcn_mom0_J'+Js[i]+'.txt', co_mom0, fmt='%e')
    
