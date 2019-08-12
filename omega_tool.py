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
def omega_calc(suffix):
#    parser = op.OptionParser(description='')
#    parser.add_option('--show_plots','-p',action='store',dest='show_plots',help = 'Display the variation as a function of z',default=False)
    pars = init_read_parameters_file(suffix)
#    options = parser.parse_args()
#    show_plots = options.show_plots
    field = fieldfile('field'+suffix,pars)
    tend = field.tfld[-1]
    tstart = field.tfld[-1]*0.9
    imax = np.unravel_index(np.argmax(abs(field.phi()[:,0,:])),(field.nz,field.nx))
    phi_t = []
    time = np.empty(0)
    istart = np.argmin(abs(np.array(field.tfld)-tstart))
    iend = np.argmin(abs(np.array(field.tfld)-tend))
    phi = np.empty(0,dtype='complex128')
    for i in range(istart,iend):
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
    zmax = len(phi_t[1])
    f=open('omega_new'+suffix,'w')
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
