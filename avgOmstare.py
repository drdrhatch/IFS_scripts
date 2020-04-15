#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os
import matplotlib.pyplot as plt
from fieldlib import *
from ParIO import * 
from finite_differences import *
import optparse as op
from read_write_geometry import *
#from read_write_geometry_lilo import *
from subprocess import call
import sys
import math
from interp import *
from read_iterdb import *
#from calc_omega_from_field import *

parser=op.OptionParser(description='Plots mode structures and calculates various interesting quantities.')
parser.add_option('--plots','-p', action='store',dest = 'plots', help = 'Plot all plots.',default=False)
parser.add_option('--pfile','-f',type = 'str',action='store',dest="prof_file",help = 'profile file name.',default='empty')
parser.add_option('--w','-w',type = 'str',action='store',dest="weightFunc",help = 'weighting function.',default='empty')
options, args = parser.parse_args()
suffix = args[0]
plots = options.plots
#idb_file = options.idb_file
prof_file = options.prof_file
weightFunc = options.weightFunc

if suffix != '.dat':
    suffix = '_'+suffix

par = Parameters()
par.Read_Pars('parameters'+suffix)
pars = par.pardict

field = fieldfile('field'+suffix,pars)
time = np.array(field.tfld)
if 1 == 1:
    itime = -1
    itime0 = len(time)-1

field.set_time(time[itime])

if 'x_local' in pars:
    if pars['x_local']:
        x_local = True
    else:
        x_local = False 
else:
    x_local = True

if 'lilo' in pars:
    if pars['lilo']:
        lilo = True
    else:
        lilo = False
else:
    lilo = False

if 1 == 1:
    dz = 2.0/field.nz
    zgrid = np.arange(field.nz)/float(field.nz-1)*(2.0-dz)-1.0
    if 'lx_a' in pars:
        xgrid = np.arange(field.nx)/float(field.nx-1)*pars['lx_a']+pars['x0']-pars['lx_a']/2.0
    else:
        xgrid = np.arange(field.nx)/float(field.nx-1)*pars['lx'] - pars['lx']/2.0
    #gpars,geometry = read_geometry_global(pars['magn_geometry'][1:-1]+suffix)
    #qmin = np.min(geometry['q'])
    #qmax = np.max(geometry['q'])
    #mmin = math.ceil(qmin*pars['n0_global'])
    #mmax = math.floor(qmax*pars['n0_global'])
    #mnums = np.arange(mmin,mmax+1)
    #print "mnums",mnums
    #qrats = mnums/float(pars['n0_global'])
    
    phi = field.phi()[:,0,:]

    apar = field.apar()[:,0,:]


    if plots:
        plt.figure(figsize=(8.0,9.5))
        fig=plt.gcf()
        fig.subplots_adjust(right=0.9)
        fig.subplots_adjust(left=0.16)
        fig.subplots_adjust(hspace=0.35)
        plt.subplot(3,1,1)
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.title(r'$|\phi|$')
        plt.contourf(xgrid,zgrid,np.abs(phi),70)
        #for i in range(len(qrats)):
        #    ix = np.argmin(abs(geometry['q']-qrats[i]))
        #    plt.axvline(xgrid[ix],color='white')
        cb1=plt.colorbar()
        plt.subplot(3,1,2)
        plt.title(r'$Re[\phi]$')
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.contourf(xgrid,zgrid,np.real(phi),70)
        plt.colorbar()
        plt.subplot(3,1,3)
        plt.title(r'$Im[\phi]$')
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.contourf(xgrid,zgrid,np.imag(phi),70)
        plt.colorbar()
        plt.show()
    if plots:
        plt.figure(figsize=(8.0,9.5))
        fig=plt.gcf()
        fig.subplots_adjust(right=0.9)
        fig.subplots_adjust(left=0.16)
        fig.subplots_adjust(hspace=0.35)
        plt.subplot(3,1,1)
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.title(r'$|A_{||}|$')
        plt.contourf(xgrid,zgrid,np.abs(apar),70)
        #for i in range(len(qrats)):
        #    ix = np.argmin(abs(geometry['q']-qrats[i]))
        #    plt.axvline(xgrid[ix],color='white')
        cb1=plt.colorbar()
        plt.subplot(3,1,2)
        plt.title(r'$Re[A_{||}]$')
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.contourf(xgrid,zgrid,np.real(apar),70)
        plt.colorbar()
        plt.subplot(3,1,3)
        plt.title(r'$Im[A_{||}]$')
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.contourf(xgrid,zgrid,np.imag(apar),70)
        plt.colorbar()
        plt.show()


    if 1 == 1:
        if prof_file == 'empty':
            prof_file = input("Enter gene output profiles file name:\n")

        profData = np.genfromtxt(prof_file)
        rho = profData[:,0]
        omt = profData[:,4]
        omn = profData[:,5]
        omstar = pars['kymin'] * (omt + omn)

        weight = 0.
        norm = 0.
        totOmstar = 0.
        for i in range(pars['nx0']):
            phiX = np.sum(abs(phi[:,i]) ** 2)
            aparX = np.sum(abs(apar[:,i]) ** 2)
            if weightFunc == 'phi':
                weight = phiX
            elif weightFunc == 'apar':
                weight = aparX
            norm = norm + weight
            totOmstar = totOmstar + omstar[i] * weight
        print('average Om_star_e = ', totOmstar / norm)
