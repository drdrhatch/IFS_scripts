#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import os
import optparse as op
from subprocess import call
from interp import *
from finite_differences import *
from read_write_geometry import *
from parIOHelper import *
#import matplotlib.pyplot as plt

parser=op.OptionParser(description='Calculates flux surface averaged grad rhot.')
options,args=parser.parse_args()
if len(args)!=2:
    exit("""
Please include parameters file name and GENE geometry file name.
    \n""")
parfile = args[0]
geomfile = args[1]

print parfile[10:]
pars = init_read_parameters(parfile[10:])
if 'x_local' in pars:
    gpars, geom = read_geometry_global(geomfile)
else:
    gpars, geom = read_geometry_local(geomfile)

print "gpars",gpars

print "gl_dxdR",geom['gl_dxdR']
print "gl_dxdZ",geom['gl_dxdZ']
abs_grad_rhot = (geom['gl_dxdR']**2 + geom['gl_dxdZ']**2)**0.5

print "abs(grad_rhot)",abs_grad_rhot
print "1/abs(grad_rhot)",1/abs_grad_rhot

zgrid = np.linspace(-np.pi,np.pi,num=pars['nz0'],endpoint = False)
print "zgrid",zgrid

plt.scatter(geom['gl_R'],geom['gl_z'])
plt.show()

plt.plot(zgrid,abs_grad_rhot)
plt.show()

#plt.plot(zgrid,1/abs_grad_rhot)
#plt.show()
