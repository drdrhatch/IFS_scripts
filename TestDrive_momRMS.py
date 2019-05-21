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
#plot_format = 'ps'
nf = 200
lf = 10.

kygrid = range(pars['nky0'])
zInd = nz/2
kyInd = -1
xInd = nx * 5 / 8

if 1 == 1:
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
    plot_format = 'ps'
    for xInd in range(nx / 4, nx * 3 / 4 + 1, 8):
        #t_rms, tperp_rms = momen_rms(tgrid, tperp_tx, xInd, show_plots, plot_format)
        t_rms, dens_rms = momen_rms(tgrid, dens_tx, xInd, show_plots, plot_format)

