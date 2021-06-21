#!/usr/bin/env python
# -*- coding: utf-8 -*-

import optparse as op
import math
import cmath
import sys
import numpy as np
import matplotlib.pyplot as plt
from fieldsWrapper import *
from parIOWrapper import init_read_parameters_file
from finite_differences import *
from fieldlib import *

#updated by Max Curie 06/21/2021

def omega_calc(suffix,plot=False):
    percent=40.#percentage of the time to calculate the omega
#    parser = op.OptionParser(description='')
#    parser.add_option('--show_plots','-p',action='store',dest='show_plots',help = 'Display the variation as a function of z',default=False)
    pars = init_read_parameters_file(suffix)
#    options = parser.parse_args()
#    show_plots = options.show_plots
    field = fieldfile('field'+suffix,pars)
    tend = field.tfld[-1]
    tstart = field.tfld[0]+(field.tfld[-1]-field.tfld[0])*(100.-percent)/100.
    #tstart = field.tfld[0]
    #print tend, tstart
    imax = np.unravel_index(np.argmax(abs(field.phi()[:,0,:])),(field.nz,field.nx))
    phi_t = []
    time = np.empty(0)
    #print np.array(field.tfld)
    istart = np.argmin(abs(np.array(field.tfld)-tstart))
    iend = np.argmin(abs(np.array(field.tfld)-tend))
    phi = np.empty(0,dtype='complex128')
    for i in range(istart,iend):
        field.set_time(field.tfld[i])
        phi, apar = eigenfunctions_from_field_file(pars,suffix,False,False,field.tfld[i],False,False)  
        phi_t.append(phi)
        time = np.append(time,field.tfld[i])
    phi_t = np.array(phi_t)
    
    #***********start of omega calculation***********

    function = np.mean(phi_t,axis=1)
    time=np.array(time)

    dt=time[1:]-time[:-1]
    dt_min=np.mean(dt)
    
    if abs(np.std(dt))>=np.min(dt)*0.01:
        print('time step is NOT uniform. interperlating')
        uni_time = np.linspace(min(time),max(time),int(abs((max(time)-min(time))/dt_min)*1.5)) #uniform time
        uni_function = np.interp(uni_time,time,function)
    else:
        uni_time=time
        uni_function=function
    
    phi_t=uni_function    
    print(np.shape(phi_t))

    time=uni_time
    
    delta_t = time[1]-time[0]
    omega_t=[]

    for t0 in range(len(time)-1):
        t=t0+1
        if phi_t[t-1]==0: continue
        elif abs(phi_t[t]/phi_t[t-1])==0: 
            omega_t.append(0)
        elif (phi_t[t]/phi_t[t-1]).imag ==0 and (phi_t[t]/phi_t[t-1]).imag < 0:
            omega_temp=complex(0,np.pi)*cmath.log(abs(phi_t[t]/phi_t[t-1]))/delta_t
            if abs(omega_temp)>10**3:
                continue
            else:
                omega_t.append(omega_temp)          
        else:
            omega_temp=cmath.log(phi_t[t]/phi_t[t-1])/delta_t
            if abs(omega_temp)>10**3:
                continue
            else:
                omega_t.append(omega_temp)

    omega_t=np.array(omega_t)
    gamma_avg=np.mean(omega_t.imag)
    gamma_std=np.std(omega_t.imag)
    omega_avg=np.mean(omega_t.real)
    omega_std=np.std(omega_t.real)
    
    if plot==True:
        plt.clf()
        plt.plot(omega_t.imag,label='imag')
        plt.plot(omega_t.real,label='real')
        plt.show()


    if abs(gamma_std/gamma_avg)>=0.3 or abs(omega_std/omega_avg)>=0.3:
        print('omega(cs/a)='+str(omega_avg)+'+-'+str(omega_std))
        print('gamma(cs/a)='+str(gamma_avg)+'+-'+str(gamma_std))
        print('No converged, please check')
        #gamma_avg=0
        #gamma_std=0
        #omega_avg=0
        #omega_std=0
    else:
        print('omega(cs/a)='+str(omega_avg)+'+-'+str(omega_std))
        print('gamma(cs/a)='+str(gamma_avg)+'+-'+str(gamma_std))
    
    
    try:
        from shutil import copyfile
        copyfile('omega'+suffix, 'old_omega'+suffix)
    except:
        pass
    
    
    f=open('omega'+suffix,'w')
    f.write('  '+str(pars['kymin'])+'    '+str(gamma_avg)+'  '+str(omega_avg))
    f.close()
    
    return gamma_avg,gamma_std,omega_avg,omega_std

    
    
#       if show_plots:
#           plt.plot((weight[:,t-1]*(gamma_avg-omega_diffs[:,t-1].real))**2/weight[:,t-1])
   
