#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os
import matplotlib.pyplot as plt
from fieldlib import *
from ParIO import * 
import optparse as op
#from read_write_geometry import *

parser=op.OptionParser(description='Some infrastructure for reading in, manipulating, and plotting nonlinear field data.')
#parser.add_option('--plot_theta','-g',action='store_const',const=False,help = 'Plot global mode structures decomposed in poloidal m number.',default=True)
options,args=parser.parse_args()
print("options",options)
print("args",args)
if len(args)!=1:
    exit("""
Please include run number as argument (e.g., 0001)."
    \n""")
suffix = args[0]

if suffix in ['dat','.dat']:
     suffix = '.dat'
else:
     suffix = '_'+suffix

par = Parameters()
par.Read_Pars('parameters'+suffix)
pars = par.pardict

#field = fieldfile('field'+suffix,pars,False)
field = fieldfile('field'+suffix,pars)
#print "time0",time0
#print "field.tfld",field.tfld
time = np.array(field.tfld)
time0 = -1
if time0 == -1:
    itime = -1
    itime0 = len(time)-1
else: 
    itime = np.argmin(abs(time - time0))
    itime0 = itime

print("Looking at the spectra at time:",time[itime])
#This sets the time step you want to read in
field.set_time(time[itime])

ntot = field.nx*field.ny*field.nz

phi = np.zeros(ntot,dtype='complex128')
apar = np.zeros(ntot,dtype='complex128')
#print "ntot",field.nz*field.nx

if 'x_local' in pars:
    if pars['x_local']:
        x_local = True
    else:
        x_local = False 
else:
    x_local = True

if x_local:
    kxmin = 2.0*np.pi/pars['lx']
    kxgrid = np.linspace(-(pars['nx0']/2-1)*kxmin,pars['nx0']/2*kxmin,num=pars['nx0'])
    kxgrid = np.roll(kxgrid,int(pars['nx0']/2+1))
    print("kxgrid",kxgrid)
    kygrid = np.linspace(0,(pars['nky0']-1)*pars['kymin'],num=pars['nky0'])
    print("kygrid",kygrid)
    zgrid = np.linspace(-np.pi,np.pi,pars['nz0'],endpoint=False)
    print("zgrid",zgrid)
    phi=field.phi()[:,:,:]
    print('np.shape(phi)',np.shape(phi))
    apar=field.apar()[:,:,:]
    phi2 = abs(phi)**2
    print('zgrid[int(pars[nz0]/2)]',zgrid[int(pars['nz0']/2)])
    phi2_outboard = phi2[int(pars['nz0']/2),:,:]
    print('np.shape(phi2_outboard)',np.shape(phi2_outboard))
    phi2_ob_ky = np.sum(phi2_outboard,axis=1)
    plt.plot(kygrid,phi2_ob_ky,label=r'$\phi^2$')
    plt.title(r'$\phi^2$'+' at t= '+str(time[itime]))
    plt.xlabel(r'$k_y(\rho_s)$')
    plt.show()
else:  #x_local = False
    pass


