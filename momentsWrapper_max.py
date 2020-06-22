#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
#from fieldlib import *
from finite_differences import *
from momlib import *
#from finite_differences_x import *

def moments_from_mom_file(pars,suffix,center_only,plot,setTime=-1):

    #momen = momfile('mom_e'+suffix,pars)
    try:
        momen = momfile('mom_e'+suffix,pars)
       # print 'momen set'
    except:
        momen = momfile('mom_electrons'+suffix,pars)
    if (setTime == -1):
        momen.set_time(momen.tmom[setTime])
#        print 'Reading moments are at t = ', momen.tmom[setTime]
    else:
        isetTime = np.argmin(abs(np.array(momen.tmom)-setTime))
        momen.set_time(momen.tmom[isetTime])
#        print 'Reading moments are at t = ', momen.tmom[isetTime]

    nz = pars['nz0']
    nx = pars['nx0']

    if center_only:
        ikx_grid = [0]
        upar = np.zeros(nz,dtype='complex128')
        deln = np.zeros(nz,dtype='complex128')
        tpar = np.zeros(nz,dtype='complex128')
        tperp = np.zeros(nz,dtype='complex128')
        qpar = np.zeros(nz,dtype='complex128')
        qperp = np.zeros(nz,dtype='complex128')
    else:
        ikx_grid = np.arange(-nx/2+1,nx/2+1)
        upar = np.zeros(nx*nz,dtype='complex128')
        deln = np.zeros(nx*nz,dtype='complex128')
        tpar = np.zeros(nx*nz,dtype='complex128')
        tperp = np.zeros(nx*nz,dtype='complex128')
        qpar = np.zeros(nx*nz,dtype='complex128')
        qperp = np.zeros(nx*nz,dtype='complex128')

    if 'n0_global' in pars:
        phase_fac = -np.e**(-2.0*np.pi*(0.0+1.0J)*pars['n0_global']*pars['q0'])
    else:
        phase_fac = -1.0
    
    if pars['shat'] > 0.:
        for i in ikx_grid:
            this_upar = momen.upar()[:,0,i]*phase_fac**i
            upar[(i-ikx_grid[0])*momen.nz:(i-ikx_grid[0]+1)*momen.nz]=this_upar
            this_deln = momen.dens()[:,0,i]*phase_fac**i
            deln[(i-ikx_grid[0])*momen.nz:(i-ikx_grid[0]+1)*momen.nz]=this_deln
            this_tpar = momen.tpar()[:,0,i]*phase_fac**i
            tpar[(i-ikx_grid[0])*momen.nz:(i-ikx_grid[0]+1)*momen.nz]=this_tpar
            this_tperp = momen.tperp()[:,0,i]*phase_fac**i
            tperp[(i-ikx_grid[0])*momen.nz:(i-ikx_grid[0]+1)*momen.nz]=this_tperp
            this_qpar = momen.qpar()[:,0,i]*phase_fac**i
            qpar[(i-ikx_grid[0])*momen.nz:(i-ikx_grid[0]+1)*momen.nz]=this_qpar
            this_qperp = momen.qperp()[:,0,i]*phase_fac**i
            qperp[(i-ikx_grid[0])*momen.nz:(i-ikx_grid[0]+1)*momen.nz]=this_qperp
    else:
        for i in ikx_grid:
           this_upar = momen.upar()[:,0,-i]*phase_fac**i
           upar[(i-ikx_grid[0])*momen.nz:(i-ikx_grid[0]+1)*momen.nz]=this_upar
           this_deln = momen.dens()[:,0,-i]*phase_fac**i
           deln[(i-ikx_grid[0])*momen.nz:(i-ikx_grid[0]+1)*momen.nz]=this_deln
           this_tpar = momen.tpar()[:,0,-i]*phase_fac**i
           tpar[(i-ikx_grid[0])*momen.nz:(i-ikx_grid[0]+1)*momen.nz]=this_tpar
           this_tperp = momen.tperp()[:,0,-i]*phase_fac**i
           tperp[(i-ikx_grid[0])*momen.nz:(i-ikx_grid[0]+1)*momen.nz]=this_tperp
           this_qpar = momen.qpar()[:,0,-i]*phase_fac**i
           qpar[(i-ikx_grid[0])*momen.nz:(i-ikx_grid[0]+1)*momen.nz]=this_qpar
           this_qperp = momen.qperp()[:,0,-i]*phase_fac**i
           qperp[(i-ikx_grid[0])*momen.nz:(i-ikx_grid[0]+1)*momen.nz]=this_qperp

    if plot:
        if (setTime == -1):
            figTitle='t = '+ str(momen.tmom[setTime])
        else:
            figTitle='t = '+ str(momen.tmom[isetTime])
        if center_only:
            figTitle = figTitle+' center only'
        else:
            figTitle = figTitle+' entire simulation domain'
        plt.plot(np.real(deln),label='Re(deln)')
        plt.plot(np.imag(deln),label='Im(deln)')
        plt.plot(abs(deln),label='abs(deln)')
        plt.title(figTitle)
        plt.legend()
        plt.show()
        plt.plot(np.real(upar),label='Re(upar)')
        plt.plot(np.imag(upar),label='Im(upar)')
        plt.plot(abs(upar),label='abs(upar)')
        plt.title(figTitle)
        plt.legend()
        plt.show()
        if 1 == 0:
            plt.plot(np.real(tpar),label='Re(tpar)')
            plt.plot(np.imag(tpar),label='Im(tpar)')
            plt.plot(abs(tpar),label='abs(tpar)')
            plt.title(figTitle)
            plt.legend()
            plt.show()
   # print upar, deln, tpar, tperp, qpar, qperp
    return upar,deln,tpar,tperp,qpar,qperp

def LILO_moments_from_mom_file(pars,suffix,plot,setTime=-1):

    momen = momfile('mom_e'+suffix,pars)
    if (setTime == -1):
        momen.set_time(momen.tmom[setTime])
        print('Reading momentss are at t = ', momen.tmom[setTime])
    else:
        isetTime = np.argmin(abs(np.array(momen.tmom)-setTime))
        momen.set_time(momen.tmom[isetTime])
        print('Reading momentss are at t = ', momen.tmom[isetTime])

    nz = pars['nz0']
    nx = pars['nx0']

    if 1 == 0:
        upar = np.zeros(nz,dtype='complex128')
        deln = np.zeros(nz,dtype='complex128')
        tpar = np.zeros(nz,dtype='complex128')
        tperp = np.zeros(nz,dtype='complex128')
        qpar = np.zeros(nz,dtype='complex128')
        qperp = np.zeros(nz,dtype='complex128')
    else:
        upar = np.zeros(nx,dtype='complex128')
        deln = np.zeros(nx,dtype='complex128')
        tpar = np.zeros(nx,dtype='complex128')
        tperp = np.zeros(nx,dtype='complex128')
        qpar = np.zeros(nx,dtype='complex128')
        qperp = np.zeros(nx,dtype='complex128')

    
    upar = momen.upar()[int(nz/2),0,:]
    deln = momen.dens()[int(nz/2),0,:]
    deln_global = momen.dens()[:,:,:]

    if plot:
        if (setTime == -1):
            figTitle='t = '+ str(momen.tmom[setTime])
        else:
            figTitle='t = '+ str(momen.tmom[isetTime])
        plt.plot(np.real(upar),label='Re(upar)')
        plt.plot(np.imag(upar),label='Im(upar)')
        plt.plot(abs(upar),label='abs(upar)')
        plt.title(figTitle)
        plt.legend()
        plt.show()
    return upar,deln,deln_global

