#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from fieldHelper import *
from momHelper import *
from geomHelper import *
from parIOHelper import *
from plotHelper import *
from windowFFT import *
import sys

suffix = sys.argv[1]

if not suffix =='.dat':
   suffix = '_'+suffix

zi = complex(0,1)
tStart = float(sys.argv[2])
tEnd = float(sys.argv[3])

pars = init_read_parameters(suffix)
field = fieldfile('field' + suffix, pars)

geom_type, geom_pars, geom_coeff = init_global_geometry(suffix, pars)
nz = pars['nz0']
nx = pars['nx0']
dz = 2.0/nz
zgrid = np.arange(nz)/float(nz-1)*(2.0-dz)-1.0
if 'lx_a' in pars:
    xgrid = np.arange(nx)/float(nx-1)*pars['lx_a']+pars['x0']-pars['lx_a']/2.0
else:
    xgrid = np.arange(nx)/float(nx-1)*pars['lx'] - pars['lx']/2.0

plot_format = 'display'
#plot_format = 'ps'

kygrid = list(range(pars['nky0']))

zInd = -1
xInd = -1

timeInd = np.argmin(abs(np.array(field.tfld) - tStart))
time, phi_xz, apar_xz = field_xz(field, geom_coeff, zgrid, kygrid, xgrid, timeInd, True, plot_format)
