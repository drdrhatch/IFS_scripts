#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from fieldHelper import *
from momHelper import *
from parIOHelper import *
from geomHelper import *
from plotHelper import *
from windowFFT import *
import sys

suffix = sys.argv[1]

if not suffix =='.dat':
   suffix = '_'+suffix

tStart = float(sys.argv[2])
tEnd = float(sys.argv[3])

pars = init_read_parameters(suffix)
momen = momfile('mom_e'+suffix,pars)

zi = complex(0, 1)
nz = pars['nz0']
nx = pars['nx0']
ny = pars['nky0']
dz = 2.0/nz
zgrid = np.arange(nz)/float(nz-1)*(2.0-dz)-1.0
if 'lx_a' in pars:
    xgrid = np.arange(nx)/float(nx-1)*pars['lx_a']+pars['x0']-pars['lx_a']/2.0
else:
    xgrid = np.arange(nx)/float(nx-1)*pars['lx'] - pars['lx']/2.0

show_plots = True
plot_format = 'display'
#plot_format = 'ps'

kygrid = range(pars['nky0'])
zInd = nz / 2
xInd = nx * 5 / 8

nf = 200
lf = 10.

#for kyInd in kygrid:
if 1 == 1:
    kyInd = -1
    tgrid, dens_tky, tperp_tky = momen_tky(momen, \
                  zInd, \
                  kyInd, \
                  xInd, \
                  tStart, \
                  tEnd)
if 1 == 0:
    title = ' '
    filename = 'dens_tky01.ps'

    doublePlot2D(kygrid, tgrid, dens_tky, tperp_tky, 'dens_tky', 'tperp_tky', title, filename, 'ky', 't',plot_format)
if 1 == 1:
#    plot_format = 'ps'
    dens_fky = np.empty((nf,ny),dtype = 'complex128')
    tperp_fky = np.empty((nf,ny),dtype = 'complex128')
    show_plots = False
    for ky in kygrid:
        fgrid, dens_f = windowFFT(tgrid, dens_tky[:,ky], nf, lf, 'dens_ky=' + str(ky), show_plots, plot_format)
        fgrid, tperp_f = windowFFT(tgrid, tperp_tky[:,ky], nf, lf, 'tperp_ky=' + str(ky), show_plots, plot_format)
        dens_f = np.minimum(dens_f, 0.002)
        tperp_f = np.minimum(tperp_f, 0.002)
        #dens_f = np.log(dens_f)
        #tperp_f = np.log(tperp_f)
        dens_fky[:,ky] = dens_f.reshape(1, nf)
        tperp_fky[:,ky] = tperp_f.reshape(1, nf)
    #doublePlot2D(kygrid, fgrid, dens_fky, tperp_fky, 'dens_fky', 'tperp_fky', title, filename, 'ky', 'f',plot_format)
    filename = 'dens_fky01.ps'
    title = ' '
if 1 == 1:
    norm = np.amax(abs(dens_fky))
    dens_plot = np.log(1. + abs(dens_fky)/norm)
    f_kHz = np.array(fgrid) * 105.
    ky_cm = np.arange(0, ny) * pars['kymin'] * 0.338 / 0.063
    singlePlotABS2D(ky_cm, f_kHz, dens_plot, 'dens_fky', title, filename, 'ky (cm^-1)', 'f (kHz)',plot_format)

if 1 == 1:
    d0_m = 0.03
    R_m = 0.68
    n_grid = pars['n0_global'] * np.arange(0, pars['nky0']) 
    f_kHz = np.array(fgrid) * 105.
    filter_GPI = np.exp(- n_grid * d0_m / R_m)
    dens_fky_filtered = dens_fky * filter_GPI
    ky_cm = np.arange(0, ny) * pars['kymin'] * 0.338 / 0.063

    norm = np.amax(abs(dens_fky_filtered))
    dens_filtered_plot = np.log(1. + abs(dens_fky_filtered)/norm)
    singlePlotABS2D(ky_cm, f_kHz, dens_filtered_plot, 'dens_fky', title, filename, 'ky (cm^-1)', 'f (kHz)',plot_format)
