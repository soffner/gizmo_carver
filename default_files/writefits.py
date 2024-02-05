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

Js = ['10'] #, '21', '32', '43', '54', '65', '76', '87', '98', '109']
moment = [0] #, 1, 2]

fwhm_x = 50 #in arcsec
fwhm_y = 50 #in arcsec
dpc = 300 # distance to source in pc

#we calculate two different moment maps: one for X_CO, in K km/s, and another for CO SLED and line ratio r, in Jy/pix * km/s


#CO
for i in range(0, len(Js)):
    freq = getattr(inputs, costr + Js[i])
    print('frequency used over CO 1-0 frequency is ', freq/inputs.restfreq_CO_J10)

    m_image = radmc3dImage()
    m_image.readImage('uclchem/image_CO_uclchem_J'+Js[i]+'.out')

    m_image.writeFits(fname='uclchem_co_J10_unconvolved.fits', dpc=dpc, nu0=freq, spectral_axis_vel=True)

    #convolve with beam
    bb=m_image.imConv(dpc=dpc, pa=0., fwhm = [fwhm_x, fwhm_y], tdiam_prim=None, tdiam_sec=None)
    print('convolved with beam size ', bb.fwhm, ' arcsec')
    bb.writeFits(fname='uclchem_co_J10_convolved_beam_'+format(fwhm_x, '.1f')+'.fits', dpc=dpc, nu0=freq, spectral_axis_vel=True)

