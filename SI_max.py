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


par = Parameters()
par.Read_Pars('parameters'+suffix)
pars = par.pardict
#spec = par.spec_nl
#print spec
coll = pars['coll']
Lref=pars['Lref']
Tref=pars['Tref']
nref=pars['nref']
mref=pars['mref']
#trpeps=pars['trpeps']
#trpeps=float(0.5)
#q0=pars['q0']
#R=pars['major_R']
#dens=spec['dens']
#dens=1
#Te=Tref #need further code of ParIO_max.py
#me=mref*(0.27244000*10**(-3))
#vte=np.sqrt(Te/me)
#print "dens", dens




#Reference from GENE manual page 29 this is in gauss unit
#coll_gauss=-2.3031*10**(-5)*Lref*nref*(24-np.log(np.sqrt(nref*10**13)/(10**3*Tref)))/(Tref)**2
coll_c=2.3031*10**(-5)*Lref*nref*(24-np.log(np.sqrt(nref*10**13)/(10**3*Tref)))/(Tref**2)
#coll_SI=coll_c*(unit_SI/unit_GENE)
coll_ei=pars['nu_ei']
#temp=Te
#print "checkpoint", trpeps

#coll_ei=(8/(3*np.sqrt(np.pi)))*(q0*(1**4)/(trpeps)**(3/2))*R*(dens/nref)*(Tref/Te)**2*coll_c
print(("ky rhos",omega[0]))
print(("gamma (cs/a)",omega[1]))
print(("omega (cs/a)",omega[2]))
print(("coll  (cs/a)", coll))
print(("nu_ei (cs/a)", coll_ei))

print(("Lref",pars['Lref']))
print(("Tref",pars['Tref']))
print(("mref",pars['mref']))

TJ = pars['Tref']*1000.0*1.602e-19
mi = pars['mref']*1.6726e-27
cs = np.sqrt(TJ/mi)
om_ref = cs/pars['Lref']

print(("cs/Lref",om_ref))
print(("coll_ei (kHz)", coll_ei*om_ref/1000.0))
print(("frequency (kHz)",omega[2]*om_ref/1000.0/(2.0*np.pi)))


