"""
   radmc_tau.py

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

Js = ['10'] #, '21', '32', '43', '54', '65', '76', '87', '98', '109']

#we calculate the max optical depth in a given spaxel by taking max of optical depths over aall spectral channels


#CO
for i in range(0, len(Js)):
    freq = getattr(inputs, costr + Js[i])
    print('frequency used over CO 1-0 frequency is ', freq/inputs.restfreq_CO_J10)

    m_image = radmc3dImage()
    m_image.readImage('despotic/tau_CO_despotic_J'+Js[i]+'.out')
    bb = m_image.image
    cc = np.max(bb, axis=2)
    np.savetxt('despotic/despotic_co_tau_J'+Js[i]+'.txt', cc, fmt='%e')

    #UCLCHEM
    m_image = radmc3dImage()
    m_image.readImage('uclchem/tau_CO_uclchem_J'+Js[i]+'.out')
    bb = m_image.image
    cc = np.max(bb, axis=2)
    np.savetxt('uclchem/uclchem_co_tau_J'+Js[i]+'.txt', cc, fmt='%e')

    print('')

print('All done!')
