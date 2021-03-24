#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from fieldlib import *
from geomHelper import *
from plotHelper import *

zi = complex(0,1)

def field_step_time(field, \
                    show_plots = True):
    field_time = field.tfld
    step_time = np.array(field_time[1:-1]) - np.array(field_time[0:-2])
    if show_plots:
        plt.plot(step_time)
        plt.title('Phi and Apar')
        plt.ylabel('step time (Lref / cref)')
        plt.show()
def global_eigenfunctions(field, \
                          zInd, \
                          kyInd, \
                          xInd, \
                          setTime = -1, \
                          show_plots = False, \
                          plot_format = 'display'):
    field.set_time(field.tfld[setTime])
    time = field.tfld[setTime]
    print( 'Reading eigenfunctions are at t = '+str( time))
    nz = field.nz
    ny = field.ny
    nx = field.nx
    dz = 2.0/nz
    zgrid = np.arange(nz)/float(nz-1)*(2.0-dz)-1.0
    if 'lx_a' in field.pars:
        xgrid = np.arange(nx)/float(nx-1)*field.pars['lx_a']+field.pars['x0']-field.pars['lx_a']/2.0
    else:
        xgrid = np.arange(nx)/float(nx-1)*field.pars['lx'] - field.pars['lx']/2.0
    if zInd == -1 and kyInd != -1 and xInd == -1:
        phi = field.phi()[0 : nz, kyInd, 0 : nx]
        apar = field.apar()[0 : nz, kyInd, 0 : nx]
    elif zInd != -1 and kyInd != -1 and xInd == -1:
        phi = field.phi()[zInd, kyInd, 0 : nx]
        apar = field.apar()[zInd, kyInd, 0 : nx]
    elif zInd != -1 and kyInd == -1 and xInd != -1:
        phi = field.phi()[zInd, 0:ny, xInd]
        apar = field.apar()[zInd, 0:ny, xInd]
    phi = phi*field.pars['rhostar']
    apar = apar*field.pars['rhostar']
    if show_plots:
        title = 'ky=' + str(kyInd)
        filename = 'n='+str(kyInd*6)+'_phi_apar_time='+str(np.round(time,4))+'.ps'
        doublePlot2D(xgrid, zgrid, phi, apar, 'phi', 'apar', title, filename, 'x', 'z', plot_format)
    return time, phi, apar
def field_xz(field, \
             geom_coeff, \
             zgrid, \
             kygrid, \
             xgrid, \
             timeInd = -1, \
             show_plots = False, \
             plot_format = 'display'):
    show_raw_plots = False
    q, Cy = q_Cy(geom_coeff)
    nGrid = np.array(kygrid)*field.pars['n0_global']
    thetaGrid = zgrid * np.pi
    thetaqMatrix = np.outer(thetaGrid, q)
    phi_xz = np.zeros((len(zgrid),len(q)), dtype = 'complex128')
    apar_xz = np.zeros((len(zgrid),len(q)), dtype = 'complex128')
    for ky in kygrid:
        time, this_phi, this_apar = global_eigenfunctions(field, -1, ky, -1, timeInd, show_raw_plots, plot_format)
        phi_xz += np.multiply(this_phi, np.exp(zi * nGrid[ky] * thetaqMatrix))
        apar_xz += np.multiply(this_apar, np.exp(zi * nGrid[ky] * thetaqMatrix))
        if ky != 0:
            phi_xz += np.multiply(np.conj(this_phi), np.exp(- zi * nGrid[ky] * thetaqMatrix))
            apar_xz += np.multiply(np.conj(this_apar), np.exp(- zi * nGrid[ky] * thetaqMatrix))
    if show_plots:
        title = 'time='+str(np.round(time,4))
        filename = 'phi_apar_time='+str(np.round(time,4))+'.ps'
#        singlePlot2D(xgrid, zgrid, phi_xz, 'dens_xz', title, filename, 'x', 'z', 'display')
        doublePlot2D(xgrid, zgrid, phi_xz, apar_xz, 'phi_xz', 'apar_xz', title, filename, 'x', 'z', plot_format)
    return time, phi_xz, apar_xz
            

def field_tx(field, \
             geom_coeff, \
             zgrid, \
             kygrid, \
             xgrid, \
             zInd, \
             tStart, \
             tEnd, \
             show_xz = False, \
             plot_format = 'display'):

    itStart = np.argmin(abs(np.array(field.tfld) - tStart))
    itEnd = np.argmin(abs(np.array(field.tfld) - tEnd))
    tsteps = itEnd - itStart + 1
    tgrid = []
    nz = field.nz
    nx = field.nx
    phi_tx = np.zeros((tsteps, nx), dtype='complex128')
    apar_tx = np.zeros((tsteps, nx), dtype='complex128')
    for timeInd in range(itStart, itEnd + 1):
        phi_x = np.zeros(field.nx, dtype='complex128')
        apar_x = np.zeros(field.nx, dtype='complex128')
        if show_xz:
            time, phi_xz, apar_xz = field_xz(field, geom_coeff, zgrid, kygrid, xgrid, timeInd, True, plot_format)
        else:
            time, phi_xz, apar_xz = field_xz(field, geom_coeff, zgrid, kygrid, xgrid, timeInd)
        phi_x = phi_xz[zInd,:]
        apar_x = apar_xz[zInd,:]
        phi_tx[timeInd - itStart, :] = phi_x.reshape(1, field.nx)
        apar_tx[timeInd - itStart, :] = apar_x.reshape(1, field.nx)
        tgrid.append(time)
    return tgrid, phi_tx, apar_tx

def local_eigenfunctions(pars, \
                         suffix, \
                         center_only = False, \
                         plot = True, \
                         setTime = -1):

    field = fieldfile('field'+suffix,pars)

    if (setTime == -1):
        field.set_time(field.tfld[setTime])
        print ('Reading eigenfunctions are at t = ', field.tfld[setTime])
    else:
        isetTime = np.argmin(abs(np.array(field.tfld)-setTime))
        field.set_time(field.tfld[isetTime])
        print ('Reading eigenfunctions are at t = ', field.tfld[isetTime])

    if center_only:
        ikx_grid = [0]
        phi = np.zeros(field.nz,dtype='complex128')
        apar = np.zeros(field.nz,dtype='complex128')
    else:
        ikx_grid = np.arange(-field.nx/2+1,field.nx/2+1)
        phi = np.zeros(field.nx*field.nz,dtype='complex128')
        apar = np.zeros(field.nx*field.nz,dtype='complex128')

    if 'n0_global' in pars:
        phase_fac = - np.e ** (- 2. * np.pi * \
                    zi * pars['n0_global'] * pars['q0'])
    else:
        phase_fac = -1.0
    
    if pars['shat'] > 0.:
        for i in ikx_grid:
            this_phi = field.phi()[:,0,i] * phase_fac ** i
            phi[(i - ikx_grid[0]) * field.nz: \
                (i - ikx_grid[0] + 1) * field.nz] = this_phi
            if pars['n_fields'] > 1 and pars['beta'] != 0:
                this_apar = field.apar()[:,0,i] * phase_fac ** i
                apar[(i - ikx_grid[0]) * field.nz: \
                     (i - ikx_grid[0] + 1) * field.nz] = \
                this_apar
    else:
        for i in ikx_grid:
            this_phi = field.phi()[:,0,-i] * phase_fac ** i
            phi[(i - ikx_grid[0]) * field.nz: \
                (i - ikx_grid[0] + 1) * field.nz] = this_phi
            if pars['n_fields'] > 1 and pars['beta'] != 0:
                this_apar = field.apar()[:,0,-i] * phase_fac ** i
                apar[(i - ikx_grid[0]) * field.nz: \
                     (i - ikx_grid[0] + 1) * field.nz] = \
                this_apar

    if plot:
        if (setTime == -1):
            figTitle='t = '+ str(field.tfld[setTime])
        else:
            figTitle='t = '+ str(field.tfld[isetTime])
        if center_only:
            figTitle = figTitle+' center only'
        else:
            figTitle = figTitle+' entire simulation domain'
        plt.plot(np.real(phi),label='Re(phi)')
        plt.plot(np.imag(phi),label='Im(phi)')
        plt.plot(abs(phi),label='abs(phi)')
        plt.title(figTitle)
        plt.legend()
        plt.show()
    if plot and pars['n_fields'] > 1 and pars['beta'] != 0:
        plt.plot(np.real(apar),label='Re(apar)')
        plt.plot(np.imag(apar),label='Im(apar)')
        plt.plot(abs(apar),label='abs(apar)')
        plt.title(figTitle)
        plt.legend()
        plt.show()
    return phi, apar

def calc_dphidz(pars, geom_coeff, suffix, show_plots = True):

    phi, apar = local_eigenfunctions(pars, suffix, False, False)
    zgrid, jacobian =  zGrid(geom_coeff, pars, False, False)
    if 1 == 1:
        nx = pars['nx0']
        nz = pars['nz0']
        if nx % 2 == 1:
            zgrid = np.linspace(- nx, nx, nx * nz, \
                         endpoint = False)
        else :
            zgrid = np.linspace(- (nx - 1), (nx + 1), \
                         nx * nz, endpoint = False)

    dphidz_mp = np.zeros(pars['nx0']*pars['nz0'],dtype='complex128')
    for i in np.arange(len(zgrid)-1):
        dphidz_mp[i] = (phi[i+1]-phi[i])/(zgrid[i+1]-zgrid[i]) * jacobian[i]

    if show_plots:
        plt.plot(zgrid, dphidz_mp, label = 'd phi / d z')
        plt.xlabel('z')
        plt.legend()
        plt.show()
    return zgrid, dphidz_mp, apar

def epar(pars, geom_coeff, suffix, show_plots = True):
    zgrid, dphidz_mp, apar = calc_dphidz(pars, geom_coeff, suffix)
    omega_filename = 'omega'+suffix
    omega_array = np.genfromtxt(omega_filename)
    #GENE has the convention of omega = -freq+i*gamma
    omega_complex = np.complex(-omega_array[2],omega_array[1])
    E_par_sum =sum(abs(-dphidz_mp+zi*omega_complex*apar))/(sum(abs(dphidz_mp))+sum(abs(zi*omega_complex*apar)))
    print ('E_par = ', E_par_sum)
    E_par =-dphidz_mp+zi*omega_complex*apar
    if show_plots:
        plt.plot(zgrid,np.abs(E_par),label=r'$|E_{\parallel}|$',color='black')
        plt.xlabel('z')
        plt.legend()
        plt.title('ky = '+str(pars['kymin']))
        plt.grid()
        plt.show()
