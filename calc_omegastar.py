#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import matplotlib.pyplot as plt
import optparse as op
import numpy as np
from read_EFIT_file import *
from subprocess import call
from finite_differences import *

parser=op.OptionParser(description='Calculates omegastar profile from an EFIT file and gene profile files.  Arguments: efit_file_name, gene_profiles_file_name_e, gene_file_name_i')
options,args=parser.parse_args()
if len(args)!=3:
    exit("""
Please include efit_file_name, gene_profiles_file_name_e, gene_profiles_file_name_i."
    \n""")

efit_file_name = args[0]
gene_profiles_e = args[1]
gene_profiles_i = args[2]

mi = 1.673e-27
me = 9.109e-31
ee = 1.602e-19
mref = 2.0*mi

def smooth_q(psi,q):
    isep = np.argmin(abs(psi-1))
    if psi[isep] > 1:
        isep -+ 1
    q[isep:] = q[isep]
    return q

Lref, Bref, R_major, q0 = get_dimpar_pars(efit_file_name,0.9)
profe = np.genfromtxt(gene_profiles_e)
profi = np.genfromtxt(gene_profiles_i)

call(['my_efit_tools.py',efit_file_name,'-n','-c'])
call(['my_efit_tools.py',efit_file_name,'-n','-p'])
Binfo = np.genfromtxt('Binfo_'+efit_file_name)
rtrp = np.genfromtxt('rt_rp_'+efit_file_name)
q0 = Binfo[:,4]
psi = Binfo[:,1]
q0 = smooth_q(psi,q0)

rhop0 = Binfo[:,1]**0.5

rhote = profe[:,0]
#irhot99 = np.argmin(abs(rhote-0.99))
#irhot99 = len(rhote)-1
Te = profe[:,2]
Ti = profi[:,2]
ne = profe[:,3]
ni = profi[:,3]

#Interpolate q profile onto rhote grid
#qprof = full_interp_lin(q0,rhop0,rtrp[:,1],rtrp[:,0],rhote)
qprof = full_interp(q0,rhop0,rtrp[:,1],rtrp[:,0],rhote)
rhope = interp(rtrp[:,0],rtrp[:,1],rhote)
plt.plot(rhope,qprof,label='Interpolated')
plt.plot(rhop0,q0,label='Original')
plt.axis([0.0,1.0,0.0,10.0])
plt.xlabel('rhop')
plt.ylabel('q')
plt.title('Verify interpolated q profile')
plt.legend()
plt.show()

dummy = input("Assuming reference mass is deuterium (press any key to continue).\n")
Z = 1.0

crefSI = (Te*1000.0*ee/mref)**0.5
cref_gene = 9787.1518*np.sqrt((Te*1000.0)/2.0)
OmrefSI = ee*Bref/mref 
rhosSI = crefSI/OmrefSI
rhostar = rhosSI/Lref
#rhostar_gene = 3.2255E-3 * np.sqrt((2.0)*Te)/Bref/(minor_r*Lref)

dTe_drhot = -fd_d1_o4_uneven(Te,rhote)
dne_drhot = -fd_d1_o4_uneven(ne,rhote)
omte = dTe_drhot/Te
omne = dne_drhot/ne

omegastar = rhostar*crefSI/Lref*qprof/rhote*(omte+omne)
omegastarkHz = omegastar/2/np.pi/1000.0

plt.plot(rhote,omegastarkHz)
plt.xlabel('rhot')
plt.ylabel('omegastar (kHz)')
plt.title('omegastar for n=1')
plt.show()

print(len(rhote))
print(len(omegastarkHz))
f = open('omegastar.dat','w')
f.write('#1.rhot 2.omegastar(kHz)\n')
f.write('#n=1\n')
np.savetxt(f,np.column_stack((rhote,omegastarkHz)))
f.close()
