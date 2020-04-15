#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from fieldHelper import *
from momHelper import *
from parIOHelper import *
from plotHelper import *
from windowFFT import *
import sys

suffix = sys.argv[1]

if not suffix =='.dat':
   suffix = '_'+suffix

tStart = float(sys.argv[2])
tEnd = float(sys.argv[3])

pars = init_read_parameters(suffix)
field = fieldfile('field' + suffix, pars)

nz = pars['nz0']
nx = pars['nx0']
dz = 2.0/nz
zgrid = np.arange(nz)/float(nz-1)*(2.0-dz)-1.0
if 'lx_a' in pars:
    xgrid = np.arange(nx)/float(nx-1)*pars['lx_a']+pars['x0']-pars['lx_a']/2.0
else:
    xgrid = np.arange(nx)/float(nx-1)*pars['lx'] - pars['lx']/2.0

#plot_format = 'display'
plot_format = 'ps'

kygrid = list(range(pars['nky0']))

zInd = -1
xInd = -1
timeInd = np.argmin(abs(np.array(field.tfld) - tStart))
for ky in kygrid:
    time, phi_xz, apar_xz = global_eigenfunctions(field, zInd, ky, xInd, timeInd)
    title = 'ky=' + str(ky)
    filename = 'n='+str(ky*6)+'_phi_apar_time='+str(np.round(time,4))+'.ps'
    doublePlot2D(xgrid, zgrid, phi_xz, apar_xz, 'phi_xz', 'apar_xz', title, filename, 'x', 'z', plot_format)


