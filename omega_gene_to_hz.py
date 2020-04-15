#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import optparse as op
from ParIO import * 

parser=op.OptionParser(description='Converts from GENE frequencies to Hz.')
options,args=parser.parse_args()
if len(args)!=1:
    exit("""
Please include run number as argument (e.g., 0001)."
    \n""")
suffix = args[0]
suffix = '_'+suffix

omega = np.genfromtxt('omega'+suffix)
print("ky rhos",omega[0])
print("gamma (cs/a)",omega[1])
print("omega (cs/a)",omega[2])

par = Parameters()
par.Read_Pars('parameters'+suffix)
pars = par.pardict

print("Lref",pars['Lref'])
print("Tref",pars['Tref'])
print("mref",pars['mref'])

TJ = pars['Tref']*1000.0*1.602e-19
mi = pars['mref']*1.6726e-27
cs = np.sqrt(TJ/mi)
om_ref = cs/pars['Lref']

print("cs/Lref",om_ref)

print("omega (kHz)",omega[2]*om_ref/1000.0/2.0/np.pi)




