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
geom_type, geom_pars, geom_coeff = init_global_geometry(suffix, pars)
pref, cref, Omegaref, rhoref, rhorefStar, Aparref, Gamma_gb, Q_gb = otherRef(suffix, pars)
fref_kHz = cref / pars['Lref'] / 2. / np.pi / 1000.

zi = complex(0, 1)
nz = pars['nz0']
nx = pars['nx0']
dz = 2.0/nz
zgrid = np.arange(nz)/float(nz-1)*(2.0-dz)-1.0
if 'lx_a' in pars:
    xgrid = np.arange(nx)/float(nx-1)*pars['lx_a']+pars['x0']-pars['lx_a']/2.0
else:
    xgrid = np.arange(nx)/float(nx-1)*pars['lx'] - pars['lx']/2.0

show_plots = True
plot_format = 'display'
nf = 200
lf = 6.

kygrid = list(range(pars['nky0']))
zInd = nz/2
kyInd = -1

if 1 == 1:
    show_plots = False
    tgrid, dens_tx, tperp_tx = momen_tx(momen, \
                  geom_coeff, \
                  zgrid, \
                  kygrid, \
                  xgrid, \
                  zInd, \
                  tStart, \
                  tEnd, \
                  show_plots, \
                  plot_format)
if 1 == 1:
    f = open('radiometer.txt','w')
    sys.stdout = f
    print('start time ='+str(tStart)+', end time ='+str(tEnd))
    show_plots = False
    plot_format = 'ps'
    for xInd in range(nx / 4, nx * 3 / 4, 8):
        print('x=', xgrid[xInd])
        fgrid, dens_f = windowFFT(tgrid, dens_tx[:,xInd], nf, lf, 'dens_x=' + str(np.round(xgrid[xInd],4)), show_plots, plot_format)
        fgrid, tperp_f = windowFFT(tgrid, tperp_tx[:,xInd], nf, lf, 'tperp_' + str(np.round(xgrid[xInd],4)), show_plots, plot_format)
        wcm = radiometer(fgrid, tperp_f, 50., 300., fref_kHz)
        noise = radiometer(fgrid, tperp_f, 350., 600., fref_kHz)
        print('relative Tperp fluct level (50 ~ 300 kHz)='+str(np.round(wcm,5)))
        print('relative Tperp fluct level (350 ~ 600 kHz)='+str(np.round(noise,5)))
        print('relative Tperp fluct level='+str(np.round(wcm - noise,5)))
        wcm = radiometer(fgrid, dens_f, 50., 300., fref_kHz)
        noise = radiometer(fgrid, dens_f, 350., 600., fref_kHz)
        print('relative dens fluct level (50 ~ 300 kHz)='+str(np.round(wcm,5)))
        print('relative dens fluct level (350 ~ 600 kHz)='+str(np.round(noise,5)))
        print('relative dens fluct level='+str(np.round(wcm - noise,5)))
