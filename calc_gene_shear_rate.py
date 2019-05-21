#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculates the GENE shear rate from the rbs file at a given value of rho_tor.
"""
#import matplotlib.pyplot as plt
import optparse as op
import numpy as np
from finite_differences import *
from interp import *
import matplotlib.pyplot as plt

plot_shear_rate = True

parser=op.OptionParser(description='Calculates the GENE shear rate from the rbs file at a given value of rho_tor.')
options,args=parser.parse_args()
if len(args)!=2:
    exit("""
Please include radial location (rho_tor) and lx_a."
    \n""")

def calc_a():
   f=open('rbsProfs','r')
   rbs = f.read()
   f.close()
   rbs = rbs.split('\n')  
   a_factor = float(rbs[1].split()[3])
   rbs = np.genfromtxt('rbsProfs')
   isep = np.argmin(abs(rbs[:,0]-1.0))
   a = rbs[isep,22]
   return a_factor*a

e0 = 1.6e-19
mi = 2.014*1.66e-27

rt0=float(args[0])
lx_a=float(args[1])

rin = rt0 - lx_a/2.0
rout = rt0 + lx_a/2.0

print "Calculating average gene shear rate centered at ",rt0," with rhot_min = ",rin," and rhot_max = ",rout

data = np.genfromtxt('rbsProfs')
rho_tor = np.arange(4000)/3999.0

irbs = np.argmin(abs(data[:,0]-rt0))
i0 = np.argmin(abs(rho_tor-rt0))

Er0 = interp(data[:,0],data[:,16],rho_tor)
Bp0 = interp(data[:,0],data[:,25],rho_tor)
R0 = interp(data[:,0],data[:,24],rho_tor)
q0 = interp(data[:,0],data[:,23],rho_tor)
Te0 = interp(data[:,0],data[:,5],rho_tor)

omt0 = Er0/(Bp0*R0)
domt0 = fd_d1_o4(omt0,rho_tor)
gamma_ExB = rho_tor/q0*domt0

cs = np.sqrt(1000.0*Te0*e0/mi)
a = calc_a()

gamma_ExB_norm = gamma_ExB*a/cs
#plt.plot(rho_tor,gamma_ExB_norm)
#plt.show()
irin = np.argmin(abs(rho_tor-rin))
irout = np.argmin(abs(rho_tor-rout))
gam_avg_norm = np.sum(abs(gamma_ExB_norm[irin:irout]))/(irout-irin)
gam_avg_raw = np.sum(abs(gamma_ExB[irin:irout]))/(irout-irin)

print "Minor radius a:", a
print "Sound speed cs at "+str(rt0)+": ",cs[i0]
print "GENE ExB shear rate at "+str(rt0)+": "+str(gamma_ExB_norm[i0])
print "GENE ExB shear average from "+str(rin)+" to "+ str(rout) +": "+str(gam_avg_norm)
print "GENE ExB shear average (normalized to box center) from "+str(rin)+" to "+ str(rout) +": "+str(gam_avg_raw*a/cs[i0])

if plot_shear_rate:
    plt.plot(rho_tor,gamma_ExB_norm)
    plt.xlabel(r'$\rho_{tor}$',size=18)
    plt.ylabel(r'$\gamma (c_s/a)$')
    plt.show()


