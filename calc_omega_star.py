#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import matplotlib.pyplot as plt
import optparse as op
import numpy as np
import matplotlib.pyplot as plt
from interp import *
from finite_differences import *


parser=op.OptionParser(description='Calculates gradient omega* from gene profile files normalized to values at selected rho_tor.  Arguments: profiles_e, profiles_i, rho_tor')
options,args=parser.parse_args()
if len(args)!=3:
    exit("""
Please include profiles_e and profiles_i file names."
    \n""")

f_ge = args[0]
f_gi = args[1]
rhot0 = float(args[2])

gene_e = np.genfromtxt(f_ge)
gene_i = np.genfromtxt(f_gi)

rhote = gene_e[:,0]
te = gene_e[:,2]
ne = gene_e[:,3]

rhoti = gene_i[:,0]
ti = gene_i[:,2]
ni = gene_i[:,3]

rhot1 = np.arange(1000)/999.0*(rhote[-1]-rhote[0])+rhote[0]
n1 = interp(rhote,ne,rhot1)
t1 = interp(rhote,te,rhot1)
rhot2 = np.arange(1000)/999.0*(rhoti[-1]-rhoti[0])+rhoti[0]
n2 = interp(rhoti,ni,rhot2)
t2 = interp(rhoti,ti,rhot2)

omt1 = abs(1.0/t1*fd_d1_o4(t1,rhot1))
omn1 = abs(1.0/n1*fd_d1_o4(n1,rhot1))
omt2 = abs(1.0/t2*fd_d1_o4(t2,rhot2))
omn2 = abs(1.0/n2*fd_d1_o4(n2,rhot2))

irhot1= np.argmin(abs(rhot1[:]  - rhot0))
irhot2= np.argmin(abs(rhot2[:]  - rhot0))

Te0 = t1[irhot1]
Ti0 = t2[irhot2]

omegastar_Te = t1**0.5*Te0**-0.5*omt1
omegastar_ne = t1**0.5*Te0**-0.5*omn1
omegastar_Ti = t2**0.5*Ti0**-0.5*omt2
omegastar_ni = t2**0.5*Ti0**-0.5*omn2

plt.plot(rhot1,omegastar_Te,label=r'$\omega_{*Te}$')
plt.plot(rhot1,omegastar_ne,label=r'$\omega_{*ne}$')
plt.xlabel(r'$\rho_{tor}$')
plt.legend()
plt.show()

plt.plot(rhot1,omegastar_Ti,label=r'$\omega_{*Ti}$')
plt.plot(rhot1,omegastar_ni,label=r'$\omega_{*ni}$')
plt.xlabel(r'$\rho_{tor}$')
plt.legend()
plt.show()

f = open('omega_star_info_'+f_ge,'w')
f.write('#1.rhot 2.omegastar_Te 3.omegastar_ne \n')
np.savetxt(f,np.column_stack((rhot1,omegastar_Te,omegastar_ne)))
f.close()

f = open('omega_star_info_'+f_gi,'w')
f.write('#1.rhot 2.omegastar_Ti 3.omegastar_ni \n')
np.savetxt(f,np.column_stack((rhot1,omegastar_Ti,omegastar_ni)))
f.close()

