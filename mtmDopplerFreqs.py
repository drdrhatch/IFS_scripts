#!/usr/bin/env python
# -*- coding: utf-8 -*-

from read_iterdb_file import *
from finite_differences import *
import matplotlib.pyplot as plt
from interp import *


iterdb_file_name = 'DIIID174082.iterdb'
aGENE_m = 0.778704682714
Bref_Gauss = 19304.
n0_global = 13
kymin = 0.13848289
x0 = 0.975

e = 1.6*10**(-19)
mref = 2.
M_kg = 3.3*10**(-27)

rhot0, te0, ti0, ne0, ni0, nz0, vrot0 = read_iterdb_file(iterdb_file_name)

uni_rhot = np.linspace(min(rhot0),max(rhot0),len(rhot0)*10.)

te_u = interp(rhot0,te0,uni_rhot)
ne_u = interp(rhot0,ne0,uni_rhot)
vrot_u = interp(rhot0,vrot0,uni_rhot)

tprime_e = -fd_d1_o4(te_u,uni_rhot)/te_u
nprime_e = -fd_d1_o4(ne_u,uni_rhot)/ne_u

x0Ind = np.argmin(abs(uni_rhot - x0))
te_mid = te_u[x0Ind]
kyGENE = kymin * np.sqrt(te_u/te_mid)
omMTM = kyGENE*(tprime_e+nprime_e)
gyroFreq = 9.79E3/np.sqrt(mref)*np.sqrt(te_u)/aGENE_m
mtmFreq = omMTM*gyroFreq/2./np.pi/1000.

omegaDoppler = vrot_u*n0_global/2./np.pi/1E3

if 1 == 1:
    plt.plot(uni_rhot,omegaDoppler,label='Doppler Shift')
    plt.plot(uni_rhot,mtmFreq,label='Electron Diamagnetic (MTM in plasma frame)')
    plt.plot(uni_rhot,mtmFreq + omegaDoppler,label='Diamagnetic plus Doppler (MTM in lab frame)')
    plt.axis([0.92,1.,-50.,300.])
    plt.xlabel('rhot')
    plt.ylabel('frequency (kHz)')
    plt.legend(loc = 2, prop = {'size':12})
    plt.show()
    #file_name = 'DIIID_Diallo_freq'
    #plt.savefig(file_name+'.pdf', format='pdf')

