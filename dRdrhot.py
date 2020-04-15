#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import os
import optparse as op
from subprocess import call
from interp import *
from finite_differences import *
#import matplotlib.pyplot as plt

parser=op.OptionParser(description='Calculates dR / drhot, where R is outboard midplane major radius.')
options,args=parser.parse_args()
if len(args)!=2:
    exit("""
Please include EFIT file name and radial location (in rho_tor).
    \n""")
efit = args[0]
x0 = float(args[1])

Binfo_file = 'Binfo_'+efit
rtrp_file = 'rt_rp_'+efit
call(['my_efit_tools.py',efit,'-c','-n'])
call(['my_efit_tools.py',efit,'-p','-n'])
Binfo = np.genfromtxt(Binfo_file)
R0 = Binfo[:,0]
psi0 = Binfo[:,1]
rtrp = np.genfromtxt(rtrp_file)
rhot0 = rtrp[:,0]
rhop0 = rtrp[:,1]

#plt.plot(rhot0,rhop0)
#plt.ylabel('rhop')
#plt.xlabel('rhot')
#plt.show()

#plt.plot(psi0**0.5,R0)
#plt.ylabel('R')
#plt.xlabel('rhop')
#plt.show()

rhot = np.linspace(0.0,1.0,10000)

Rnew = full_interp(R0,psi0**0.5,rhop0,rhot0,rhot)

dRdrhot = fd_d1_o4(Rnew,rhot)

ix0 = np.argmin(abs(rhot-x0))

print("dR/drhot (m) at x0 ",x0)
print(dRdrhot[ix0])

