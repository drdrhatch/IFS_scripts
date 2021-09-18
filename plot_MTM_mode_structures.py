#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from fieldlib import *
from ParIO import * 
from finite_differences import *
import optparse as op
from read_write_geometry import *
from subprocess import call
import sys
import math
from interp import *
from ParIO import Parameters 
from read_iterdb import *
from read_profiles import read_profile_file
#from calc_omega_from_field import *



name_list=['095','096','097','098','099','100','101','102','103','104','105']
stable_list=[1,    0,   0,    0,    1,    1,     1,    0,    0,   0,     1]
#name_list=name_list[:2]
#stable_list=stable_list[:2]

plot_real_part_apar=False #change to Ture if one wants to plot the real and imagainmart part of Apar

path_list=['']*len(name_list)
for i in range(len(path_list)):
    path_list[i]='/global/cscratch1/sd/maxcurie/D3D174819/GL_3560_q'+name_list[i]+'/'

path='/global/u1/m/maxcurie/genecode/prob_00APS_GL_3560_n6/'
profile_name_list=[path+'DIIID174864.iterdb']*len(path_list)
profile_type='ITERDB'



suffix='_1'
time0=-1 #-1 if plot the last time slice, float if plot the given time 



if plot_real_part_apar==True:
    fig, ax=plt.subplots(nrows=4,ncols=len(path_list),sharex=True)
    for i in range(4):
        for j in range(len(path_list)):
            #if i!=len(path_list)-1:
            #    ax[i,j].set_xticklabels([])
            if j!=0:
                ax[i,j].set_yticklabels([]) 
else: 
    fig, ax=plt.subplots(nrows=2,ncols=len(path_list),sharex=True)
    for i in range(2):
        for j in range(len(path_list)):
            #if i!=len(path_list)-1:
            #    ax[i,j].set_xticklabels([])
            if j!=0:
                ax[i,j].set_yticklabels([])



for i in range(len(path_list)):


    par = Parameters()
    par.Read_Pars(path_list[i]+'parameters'+suffix)
    pars = par.pardict
    try:
        gpars,geometry = read_geometry_local(path_list[i]+pars['magn_geometry'][1:-1]+suffix)
    except:
        gpars,geometry = read_geometry_global(path_list[i]+pars['magn_geometry'][1:-1]+suffix)
    
    

    #field = fieldfile('field'+suffix,pars,False)
    field = fieldfile(path_list[i]+'field'+suffix,pars)
    #print "time0",time0
    #print "field.tfld",field.tfld
    time = np.array(field.tfld)
    if time0 == -1:
        itime = -1
        itime0 = len(time)-1
    else: 
        itime = np.argmin(abs(time - time0))
        itime0 = itime
    
    print("Looking at the mode structure at time:",time[itime])
    #field.set_time(time[itime],itime0)
    field.set_time(time[itime])
    
    ntot = field.nz*field.nx
    
    phi = np.zeros(ntot,dtype='complex128')
    apar = np.zeros(ntot,dtype='complex128')
    #print "ntot",field.nz*field.nx

    dz = 2.0/field.nz
    zgrid = np.arange(field.nz)/float(field.nz-1)*(2.0-dz)-1.0
    zgrid_ext = np.arange(field.nz+4)/float(field.nz+4-1)*(2.0+3*dz)-(1.0+2.0*dz)
    #print zgrid
    #print zgrid_ext
    if 'lx_a' in pars:
        xgrid = np.arange(field.nx)/float(field.nx-1)*pars['lx_a']+pars['x0']-pars['lx_a']/2.0
    else:
        xgrid = np.arange(field.nx)/float(field.nx-1)*pars['lx'] - pars['lx']/2.0

    if 1==1:
        profile_name=profile_name_list[i]
        geomfile_name=path_list[i]+pars['magn_geometry'][1:-1]+suffix
        rhot0, rhop0, te0, ti0, ne0, ni0, nz0, vrot0 = read_profile_file(profile_type,profile_name,geomfile_name,suffix)
        
        uni_rhot = np.linspace(min(xgrid),max(xgrid),len(xgrid)*3)
        
        te_u = interp(rhot0,te0,uni_rhot)
        ne_u = interp(rhot0,ne0,uni_rhot)
        vrot_u = interp(rhot0,vrot0,uni_rhot)
        q      = interp(xgrid,geometry['q'],uni_rhot)
        tprime_e = -fd_d1_o4(te_u,uni_rhot)/te_u
        nprime_e = -fd_d1_o4(ne_u,uni_rhot)/ne_u

        x0_center = pars['x0']
        
        center_index = np.argmin(abs(uni_rhot-x0_center))
        B_gauss=10.**4              #1T=10^4Gauss
        qref = 1.6E-19              #in C
        c  = 1.                     #in 3*10^8m/s
        m_kg = 1.673E-27            #in kg
        Bref = pars['Bref']         #in Tesla
        Tref = pars['Tref']         #in keV
        nref = pars['nref']         #in 10^(19) /m^3
        Lref = pars['Lref']         #in m
        mref = pars['mref']         #in proton mass(kg)
        x0 = pars['x0']             #x/a, location
        kymin = pars['kymin']       #in rhoi
        nky0 = pars['nky0']         #in total number of ky
        n_step = pars['n0_global']  #in rhoi
        q0      = q[center_index]
        ne = ne_u[center_index]
        te = te_u[center_index] #it is in eV
        #Bref=float(Bref_Gauss)/10000
        m_SI = mref *1.6726*10**(-27)
        c  = 1
        qref = 1.6*10**(-19)
        nref = ne
        Tref = te * qref
        cref = np.sqrt(Tref / m_SI)
        Omegaref = qref * Bref / m_SI / c
        rhoref = cref / Omegaref 
        n0=1.
        m0 = n0*q0
        ky=n0*q0*rhoref/(Lref*x0_center)
        kymin = ky
        print("kymin="+str(kymin))
        n0_global = n0
        te_mid = te_u[center_index]
        kyGENE =kymin * (q/q0) * np.sqrt(te_u/te_mid) * (x0_center/uni_rhot) #Add the effect of the q varying
        #***Calculate omeage star********************************
        #from mtm_doppler
        omMTM = kyGENE*(tprime_e+nprime_e)
        gyroFreq = 9.79E3/np.sqrt(mref)*np.sqrt(te_u)/Lref
        #print("gyroFreq="+str(gyroFreq[center_index]))
        mtmFreq0 = omMTM*gyroFreq/2./np.pi/1000.
        omegaDoppler0 = vrot_u*n0_global/2./np.pi/1E3
        omega0 = mtmFreq0 + omegaDoppler0
    
    
    x,f,f_lab=uni_rhot,float(pars['n0_global'])*mtmFreq0,float(pars['n0_global'])*omega0

    xmin=np.argmin(abs(np.min(xgrid)-x))
    xmax=np.argmin(abs(np.max(xgrid)-x))
    x=x[xmin:xmax]
    f=f[xmin:xmax]
    
    #Find rational q surfaces
    qmin = np.min(geometry['q'])
    qmax = np.max(geometry['q'])
    mmin = math.ceil(qmin*pars['n0_global'])
    mmax = math.floor(qmax*pars['n0_global'])
    mnums = np.arange(mmin,mmax+1)
    print("mnums",mnums)
    qrats = mnums/float(pars['n0_global'])
    print("qrats",qrats)
    zgridm = np.arange(mmax*20)/float(mmax*20)*2.0-1.0
    nm = int(mmax*20)


    #field.set_time(field.tfld[-1],len(field.tfld)-1)
    field.set_time(field.tfld[itime])


    imax = np.unravel_index(np.argmax(abs(field.phi()[:,0,:])),(field.nz,field.nx))
    #print "imax",imax

    apar = field.apar()[:,0,:]
    apar_theta = np.zeros((nm,field.nx),dtype = 'complex128')     
    apar_m = np.zeros((nm,field.nx),dtype = 'complex128')     

    #10%
    x_min=0.96
    x_max=0.973
    

    ax[0,i].plot(xgrid,(geometry['q']-np.min(geometry['q'])*0.7)/np.max(geometry['q']-np.min(geometry['q'])*0.7),color='orange',label='safety factor')
    ax[0,i].plot(x,f/np.max(f),color='blue',label=r'$\omega_{*e}$')
    for j in range(len(qrats)):
        ix = np.argmin(abs(geometry['q']-qrats[j])) 
        if stable_list[i]==0 and x_min<=xgrid[ix] and xgrid[ix]<=x_max:
            ax[0,i].axvline(xgrid[ix],color='red')
        else:
            ax[0,i].axvline(xgrid[ix],color='green')
    
    if i==0:
        ax[0,i].set_ylabel(r'$a.u.$',fontsize=13)
        ax[1,i].set_ylabel(r'$|A_{||}|(z/\pi)$',fontsize=13)
        if plot_real_part_apar==True:
            ax[2,i].set_ylabel(r'$Re[A_{||}(z/\pi)]$',fontsize=13)
            ax[3,i].set_ylabel(r'$Im[A_{||}(z/\pi)]$',fontsize=13)

    ax[1,i].set_xlabel(r'$\rho_{tor}$',fontsize=13)
    ax[1,i].contourf(xgrid,zgrid,np.abs(apar),70)

    for j in range(len(qrats)):
        ix = np.argmin(abs(geometry['q']-qrats[j])) 
        ax[1,i].axvline(xgrid[ix],color='white')

    if plot_real_part_apar==True:
        ax[2,i].set_xlabel(r'$\rho_{tor}$',fontsize=13)
        ax[2,i].contourf(xgrid,zgrid,np.real(apar),70)
    
        for j in range(len(qrats)):
            ix = np.argmin(abs(geometry['q']-qrats[j])) 
            ax[2,i].axvline(xgrid[ix],color='white')
    
        ax[3,i].set_xlabel(r'$\rho_{tor}$',fontsize=13)
        ax[3,i].contourf(xgrid,zgrid,np.imag(apar),70)
    
        for j in range(len(qrats)):
            ix = np.argmin(abs(geometry['q']-qrats[j])) 
            ax[3,i].axvline(xgrid[ix],color='white')

    
    ax[0,i].set_title(r'$q=$'+name_list[i][0]+'.'+name_list[i][1:]+r'$q_0$')

plt.subplots_adjust(wspace=0, hspace=0)
#plt.tight_layout()
plt.show()
