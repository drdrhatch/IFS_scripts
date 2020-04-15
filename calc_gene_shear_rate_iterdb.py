#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculates the GENE shear rate from an iterdb file at a given value of rho_tor.  Arguments: rhot, lx_a, a, bfile, rtrpfile, iterdb file.
"""
#import matplotlib.pyplot as plt
import optparse as op
import numpy as np
from finite_differences import *
from interp import *
import matplotlib.pyplot as plt
from read_iterdb import *

plot_shear_rate = True

parser=op.OptionParser(description='Calculates the GENE shear rate from an iterdb file at a given value of rho_tor. Arguments: rhot, lx_a, a, bfile, rtrpfile, iterdb file.')
options,args=parser.parse_args()
if len(args)!=6:
    exit("""
Please include radial location (rho_tor), lx_a, minor radius a, bfile, rtrpfile, and iterdb file name."
    \n""")

e0 = 1.6e-19
mi = 2.014*1.66e-27

rt0=float(args[0])
lx_a=float(args[1])
a = float(args[2])
bfile = args[3]
rtrpfile = args[4]
idbfile= args[5]

rin = rt0 - lx_a/2.0
rout = rt0 + lx_a/2.0

print("Calculating average gene shear rate centered at ",rt0," with rhot_min = ",rin," and rhot_max = ",rout)

rhot, profiles, units = read_iterdb(idbfile)
rhotVR = rhot['VROT']
rhotTE = rhot['TE']
rho_tor = np.arange(4000)/3999.0

binfo = np.genfromtxt(bfile)
psi = binfo[:,1]
q = binfo[:,4]
rtrp = np.genfromtxt(rtrpfile)

q0 = full_interp(q,psi**0.5,rtrp[:,1],rtrp[:,0],rho_tor,verify_interp=True)

irbs = np.argmin(abs(rhotVR-rt0))
i0 = np.argmin(abs(rho_tor-rt0))

omtor = profiles['VROT']
omt0 = interp(rhotVR,omtor,rho_tor)
Te0 = interp(rhotTE,profiles['TE'],rho_tor)/1000.0

domt0 = fd_d1_o4(omt0,rho_tor)
gamma_ExB = rho_tor/q0*domt0

cs = np.sqrt(1000.0*Te0*e0/mi)

gamma_ExB_norm = gamma_ExB*a/cs
#plt.plot(rho_tor,gamma_ExB_norm)
#plt.show()
irin = np.argmin(abs(rho_tor-rin))
irout = np.argmin(abs(rho_tor-rout))
gam_avg_norm = np.sum(abs(gamma_ExB_norm[irin:irout]))/(irout-irin)
gam_avg_raw = np.sum(abs(gamma_ExB[irin:irout]))/(irout-irin)

print("Minor radius a:", a)
print("Sound speed cs at "+str(rt0)+": ",cs[i0])
print("GENE ExB shear rate at "+str(rt0)+": "+str(gamma_ExB_norm[i0]))
print("GENE ExB shear average from "+str(rin)+" to "+ str(rout) +": "+str(gam_avg_norm))
print("GENE ExB shear average (normalized to box center) from "+str(rin)+" to "+ str(rout) +": "+str(gam_avg_raw*a/cs[i0]))

if plot_shear_rate:
    plt.plot(rho_tor,gamma_ExB_norm)
    plt.xlabel(r'$\rho_{tor}$',size=18)
    plt.ylabel(r'$\gamma (c_s/a)$')
    plt.show()


