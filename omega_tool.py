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
from get_nrg import *

def omega_calc(suffix,alg = 1):
    if alg == 1:
        omega_calc1(suffix)
    else:
        omega_calc2(suffix)

def omega_calc1(suffix):
#    parser = op.OptionParser(description='')
#    parser.add_option('--show_plots','-p',action='store',dest='show_plots',help = 'Display the variation as a function of z',default=False)
    pars = init_read_parameters_file(suffix)
#    options = parser.parse_args()
#    show_plots = options.show_plots
    field = fieldfile('field'+suffix,pars)
    tend = field.tfld[-1]
    tstart = field.tfld[-1]*0.9
    #print "tstart",tstart
    #print "tend",tend
    imax = np.unravel_index(np.argmax(abs(field.phi()[:,0,:])),(field.nz,field.nx))
    phi_t = []
    time = np.empty(0)
    istart = np.argmin(abs(np.array(field.tfld)-tstart))
    iend = np.argmin(abs(np.array(field.tfld)-tend))
    phi = np.empty(0,dtype='complex128')
    for i in range(istart,iend):
        #print "time",field.tfld[i]
        field.set_time(field.tfld[i])
        phi, apar = eigenfunctions_from_field_file(pars,suffix,False,False,field.tfld[i],False,False)   
        phi_t.append(phi)
        time = np.append(time,field.tfld[i])
    phi_t = np.array(phi_t)
    omega_diffs = []
    omega_t = []
    weight_t = []
    weight = []
    omega_avg = np.empty(0,dtype='complex128')
    delta_t = field.tfld[1]-field.tfld[0]
    zmax = len(phi_t[0])
    f=open('new_omega'+suffix,'w')
    f.write('   '+'t'+'        '+'gamma'+'        '+'omega'+'        '+'std_gamma'+'        '+'std_omega\n')
    for t in range(1,iend-istart):
        for z in range(zmax):
            omega_t.append(cmath.log((phi_t[t,z]/phi_t[t-1,z]))/(field.tfld[t]-field.tfld[t-1]))
            weight_t.append(abs(phi_t[t,z])+ abs(phi_t[t-1,z]))
        omega_diffs = np.array([omega_t[i*zmax:(i+1)*zmax] for i in range(len(omega_t)//zmax)],dtype='complex128')
        weight = np.array([weight_t[i*zmax:(i+1)*zmax] for i in range(len(weight_t)//zmax)],dtype='float128')
        omega_avg = np.append(omega_avg,np.sum(omega_diffs[:,t]*weight[:,t])/np.sum(weight[:,t]))
        gamma_avg = omega_avg[t-1].real
        omega_avg2 = omega_avg[t-1].imag
        delta_gamma2 = np.sum(weight[:,t-1]*(gamma_avg-omega_diffs[:,t-1].real)**2)/np.sum(weight[:,t-1])
        delta_omega2 = np.sum(weight[:,t-1]*(omega_avg2-omega_diffs[:,t-1].imag)**2)/np.sum(weight[:,t-1])
        f.write(str(field.tfld[t])+'    '+str(gamma_avg)+'    '+str(omega_avg2)+'    '+str(math.sqrt(delta_gamma2))+'    '+str(math.sqrt(delta_omega2))+'\n')
    f.close()


def omega_calc2(suffix):
    calc_from_apar=False
    #suffix = '_'+suffix

    #par = Parameters()
    #par.Read_Pars('parameters'+suffix)
    #pars = par.pardict

    pars = init_read_parameters_file(suffix)
    if pars['n_spec'] == 1:
        time, nrgi = get_nrg0(suffix,nspec=1)
    elif pars['n_spec'] == 2:
        time, nrgi, nrge = get_nrg0(suffix,nspec=2)
    elif pars['n_spec'] == 3:
        time, nrgi, nrge, nrg2 = get_nrg0(suffix,nspec=3)
    else:
        sys.exit("n_spec must be 1,2,3.")

    #Check for rescaling
    print "tmax",time[-1]
    print "tmin",time[0]
    for i in range(len(time)-1):
        if abs(nrgi[i,0] - nrgi[i+1,0])/(nrgi[i,0]+nrgi[i+1,0]) > 0.8:
            print "Rescaling at :",time[i]

    if calc_from_apar:
       print "Calculating growth rate from apar."
       #plt.semilogy(time,nrge[:,7])
       #plt.xlabel('time')
       #plt.show()
    else:
       print "Calculating growth rate from phi."
       #plt.semilogy(time,nrgi[:,6])
       #plt.title('QiES')
       #plt.xlabel('time')
       #plt.show()

    tstart = float(raw_input("Enter start time: "))
    tend = float(raw_input("Enter end time: "))

    field = fieldfile('field'+suffix,pars)
    istart = np.argmin(abs(np.array(field.tfld)-tstart))
    iend = np.argmin(abs(np.array(field.tfld)-tend))

    field.set_time(field.tfld[-1])
    imax = np.unravel_index(np.argmax(abs(field.phi()[:,0,:])),(field.nz,field.nx))
    phi = np.empty(0,dtype='complex128')
    if pars['n_fields'] > 1:
        imaxa = np.unravel_index(np.argmax(abs(field.apar()[:,0,:])),(field.nz,field.nx))
        apar = np.empty(0,dtype='complex128')

    time = np.empty(0)
    for i in range(istart,iend):
        #field.set_time(field.tfld[i],i)
        field.set_time(field.tfld[i])
        phi = np.append(phi,field.phi()[imax[0],0,imax[1]])
        if pars['n_fields'] > 1:
            apar = np.append(apar,field.apar()[imaxa[0],0,imaxa[1]])
        time = np.append(time,field.tfld[i])
    print "phi_0,phi_f",phi[0],phi[-1]     
    if len(phi) < 2.0:
        output_zeros = True
        omega = 0.0+0.0J
    else:
        output_zeros = False
        if calc_from_apar:
            print "Calculating omega from apar"
            if pars['n_fields'] < 2:
                stop
            omega = np.log(apar/np.roll(apar,1))
            dt = time - np.roll(time,1)
            omega /= dt
            omega = np.delete(omega,0)
            time = np.delete(time,0)
        else:
            omega = np.log(phi/np.roll(phi,1))
            dt = time - np.roll(time,1)
            omega /= dt
            print 'omega',omega
            omega = np.delete(omega,0)
            time = np.delete(time,0)
    
    gam_avg = np.average(np.real(omega))
    om_avg = np.average(np.imag(omega))
    
    if output_zeros:
        sys.exit( "Error: not enough time points in field selection.")
    else:
        print "Gamma:",gam_avg
        print "Omega:",om_avg
        f=open('new_omega'+suffix,'w')
        f.write('# omega_calc2\n   '+'kymin'+'        '+'gamma'+'        '+'omega\n'+'        '+'std_gamma'+'        '+'std_omega\n')
        f.write(str(pars['kymin'])+'    '+str(gam_avg)+'    '+str(om_avg)+'    '+str(np.nan)+ '    '+
    
str(np.nan)+'\n')
        f.close()
        
