#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from momlib import *
from geomHelper import *
from plotHelper import *
from windowFFT import *
from parIOHelper import *

zi = complex(0, 1)

def momen_step_time(momen, \
                    show_plots = True):
    momen_time = momen.tmom
    step_time = np.array(momen_time[1:-1]) - np.array(momen_time[0:-2])
    if show_plots:
        plt.plot(step_time)
        plt.title('dens, Tperp')
        plt.ylabel('step time (Lref / cref)')
        plt.show()
def global_moments(momen, \
                   zInd, \
                   kyInd, \
                   xInd, \
                   setTime = - 1, \
                   show_plots = False, \
                   plot_format = 'display'):
    momen.set_time(momen.tmom[setTime])
    time = momen.tmom[setTime]
    print 'Reading moments are at t = ', time
    nz = momen.pars['nz0']
    nky = momen.pars['nky0']
    nx = momen.pars['nx0']
    dz = 2.0/nz
    zgrid = np.arange(nz)/float(nz-1)*(2.0-dz)-1.0
    if 'lx_a' in momen.pars:
        xgrid = np.arange(nx)/float(nx-1)*momen.pars['lx_a']+momen.pars['x0']-momen.pars['lx_a']/2.0
    else:
        xgrid = np.arange(nx)/float(nx-1)*momen.pars['lx'] - momen.pars['lx']/2.0
    if zInd == -1 and kyInd != -1 and xInd == -1: 
        deln = momen.dens()[0 : nz, kyInd, 0 : nx]
        tperp = momen.tperp()[0 : nz, kyInd, 0 : nx]
    elif zInd != -1 and kyInd != -1 and  xInd == -1:
        deln = momen.dens()[zInd, kyInd, 0 : nx]
        tperp = momen.tperp()[zInd, kyInd, 0 : nx]
    elif zInd != -1 and kyInd == -1 and  xInd != -1:
        deln = momen.dens()[zInd, 0 : nky, xInd]
        tperp = momen.tperp()[zInd, 0 : nky, xInd]
    elif zInd != -1 and kyInd != -1 and  xInd != -1:
        deln = momen.dens()[zInd, kyInd, xInd]
        tperp = momen.tperp()[zInd, kyInd, xInd]
    deln = deln * momen.pars['rhostar']
    tperp = tperp * momen.pars['rhostar']
    if show_plots:# and i == momen.pars['nky0'] - 1:
        title = 'ky=' + str(kyInd)
        filename = 'n='+str(kyInd*6)+'_dens_tperp_time='+str(np.round(time,4))+'.ps'
#        singlePlot2D(xgrid, zgrid, dens, 'dens_xz', title, filename, 'x', 'z', 'display')
        doublePlot2D(xgrid, zgrid, deln, tperp, 'dens', 'tperp', title, filename, 'x', 'z', plot_format)
    return time, deln, tperp
def momen_xz(momen, \
             geom_coeff, \
             zgrid, \
             kygrid, \
             xgrid, \
             timeInd = -1, \
             show_plots = False, \
             plot_format = 'display'):
    debug = False
    show_raw_plots = False
    q, Cy = q_Cy(geom_coeff)
    nGrid = np.array(range(momen.pars['nky0']))*momen.pars['n0_global']
    thetaGrid = zgrid * np.pi
    thetaqMatrix = np.outer(thetaGrid, q)
    if debug:
        print 'zi='+str(zi)
        print 'n0='+str(momen.pars['n0_global'])
        print q
        print nGrid
        print thetaGrid
        print thetaqMatrix
    dens_xz = np.zeros((len(zgrid), len(q)),dtype='complex128')
    tperp_xz = np.zeros((len(zgrid), len(q)),dtype='complex128')
    for ky in kygrid:
        time, this_dens, this_tperp = global_moments(momen, -1, ky, -1, timeInd, show_raw_plots, plot_format)
        dens_xz += np.multiply(this_dens, np.exp(zi * nGrid[ky] * thetaqMatrix))
        tperp_xz += np.multiply(this_tperp, np.exp(zi * nGrid[ky] * thetaqMatrix))
        if ky != 0:
            dens_xz += np.multiply(np.conj(this_dens), np.exp(- zi * nGrid[ky] * thetaqMatrix))
            tperp_xz += np.multiply(np.conj(this_tperp), np.exp(- zi * nGrid[ky] * thetaqMatrix))
    if show_plots:
        title = 'time='+str(np.round(time,4))
        filename = 'xz_dens_tperp_time='+str(np.round(time,4))+'.ps'
#        singlePlot2D(xgrid, zgrid, dens_xz, 'dens_xz', title, filename, 'x', 'z', 'display')
        doublePlot2D(xgrid, zgrid, dens_xz, tperp_xz, 'dens_xz', 'tperp_xz', title, filename, 'x', 'z', plot_format)
    return time, dens_xz, tperp_xz

def momen_tx(momen, \
                  geom_coeff, \
                  zgrid, \
                  kygrid, \
                  xgrid, \
                  zInd, \
                  tStart, \
                  tEnd, \
                  show_plots = False, \
                  plot_format = 'display'):
    show_xz = False
    itStart = np.argmin(abs(np.array(momen.tmom)-tStart))
    itEnd = np.argmin(abs(np.array(momen.tmom)-tEnd))
    tsteps = itEnd - itStart + 1
    tgrid = []
    nz = momen.pars['nz0']
    nx = momen.pars['nx0']
    deln_tx = np.zeros((tsteps, nx),dtype='complex128')
    tperp_tx = np.zeros((tsteps, nx),dtype='complex128')
    for timeInd in range(itStart, itEnd + 1):
        deln_x = np.zeros(nx,dtype='complex128')
        tperp_x = np.zeros(nx,dtype='complex128')
        time, dens_xz, tperp_xz = momen_xz(momen, geom_coeff, zgrid, kygrid, xgrid, timeInd, show_xz, plot_format)

        deln_x = dens_xz[zInd,:]
        tperp_x = tperp_xz[zInd,:]
        deln_tx[timeInd - itStart, :] = deln_x.reshape(1, nx)
        tperp_tx[timeInd - itStart, :] = tperp_x.reshape(1, nx)
        tgrid.append(time)
    if show_plots:
        title = ' '
        filename = 'tx_dens_tperp.ps'
        doublePlot2D(xgrid, tgrid, deln_tx, tperp_tx, 'dens_tx', 'tperp_tx', title, filename, 'x', 't',plot_format)
    return tgrid, deln_tx, tperp_tx
def momen_rms(tgrid, field_tx, xInd, xgrid, title, show_plots = False, plot_format='display'):
    momen_rms = []
    t_rms = []
    numerator = 0.
    denominator = 0.
    for i in range(len(tgrid) - 1):
        numerator += 0.5 * (abs(field_tx[i,xInd])**2 + \
                    abs(field_tx[i + 1,xInd])**2) * \
                    (tgrid[i + 1] - tgrid[i])
        denominator += tgrid[i + 1] - tgrid[i]
        if i > 10:
            momen_rms.append(np.sqrt(numerator / denominator))
            t_rms.append(0.5 * (tgrid[i] + tgrid[i + 1]))
    if show_plots:
        plt.figure()
        title = title + str(np.round(xgrid[xInd],4))
        #filename = 'rms_tperp_xInd='+str(xInd)+'.ps'
        filename = 'rms_dens_x='+str(np.round(xgrid[xInd],4))+'.ps'
        plt.plot(t_rms, momen_rms)
        plt.xlabel('t')
        plt.title(title)
        if plot_format == 'display':
            plt.show()
        elif plot_format == 'ps':
            fig=plt.gcf()
            fig.savefig(filename, format = 'ps', bbox_inches = 'tight')
    return t_rms, momen_rms
def momen_tky(momen, \
                  zInd, \
                  kyInd, \
                  xInd, \
                  tStart, \
                  tEnd):
    nky = momen.pars['nky0']
    itStart = np.argmin(abs(np.array(momen.tmom)-tStart))
    itEnd = np.argmin(abs(np.array(momen.tmom)-tEnd))
    tsteps = itEnd - itStart + 1
    tgrid = []
    deln_tky = np.zeros((tsteps, nky),dtype='complex128')
    tperp_tky = np.zeros((tsteps, nky),dtype='complex128')
    for timeInd in range(itStart, itEnd + 1):
        time, this_deln, this_tperp = global_moments(momen, zInd, kyInd, xInd, timeInd)
        deln_tky[timeInd - itStart, :] = this_deln.reshape(1, nky)
        tperp_tky[timeInd - itStart, :] = this_tperp.reshape(1, nky)
        tgrid.append(time)
    return tgrid, deln_tky, tperp_tky
def radiometer(fgrid, field_f, start, end, fref):
    startInd = np.argmin(abs(fgrid - start/fref))
    endInd = np.argmin(abs(fgrid - end/fref))
    #print start, end
    #print fgrid[startInd] * fref, fgrid[endInd]*fref
    numerator = 0.
    denominator = 0.
    for i in range(startInd, endInd):
        numerator += 0.5 * (abs(field_f[i])**2 + \
                     abs(field_f[i + 1])**2) * \
                     (fgrid[i + 1] - fgrid[i])
        denominator += fgrid[i + 1] - fgrid[i]
    return np.sqrt(numerator / denominator)
