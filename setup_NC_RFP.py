#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import optparse as op
from finite_differences import *
from write_pfile import *
import os
from interp import *
from read_pfile import *

parser=op.OptionParser(description='Extracts the pressure profile from an EFIT file, constructs n and T profiles consistent with the pressure profile (according to user specs), and outputs to a pfile.')
parser.add_option('--T0','-t',type = 'float',action='store',dest="T0",help = 'Core T(keV).',default=1)

options,args=parser.parse_args()
if len(args)!=1:
    exit("""
Please include efit file name as argument.
    \n""")

gfile = args[0]
T0 = options.T0

e = 1.602e-19
os.system('efit_tools.py '+gfile+' --binfo -n')

Bfile = np.genfromtxt('Binfo_'+gfile)

psi0 = Bfile[:,1]
press0 = Bfile[:,5]
Bpol0 = Bfile[:,2]

entry_length = 300
psi = np.linspace(0,1,entry_length)
press = interp(psi0,press0,psi)
Bpol = interp(psi0,Bpol0,psi)

def setup_constant_T(psiN,T0kev,P_pasc):
    e = 1.602e-19
    TJ0 = T0kev*1000*e
    if len(psiN) != len(P_pasc):
        print("Error! psi_N and P_pasc must have same length")
        stop
    TJ = np.empty(len(P_pasc))
    TJ[:] = TJ0
    Tkev = TJ/1000/e
    n = 0.5*P_pasc / TJ
    return n,TJ,Tkev

ne, TeJ, Tekev = setup_constant_T(psi,T0,press)
#plt.plot(psi,ne,label='ne(m^-3)')
#plt.legend()
#plt.show()
#plt.plot(psi,Tekev,label='Te(KeV)')
#plt.legend()
#plt.show()
#plt.plot(psi,press,label='P from efit')
#plt.plot(psi,2*ne*Tekev*1000*e,'x',label='P from n,T')
#plt.legend()
#plt.show()
print("Setting up a C impurity species with very low density.  I can't get the NEO interface to work without it.")
nz = np.empty(len(ne))
nz[:] = 1.0e15
ni = ne - 6*nz
TiJ = TeJ
Tikev = Tekev

########Calculate beta_pol (RFP style)
#psi = np.linspace(0,1,entry_length)
#press = interp(psi0,press0,psi)

isep = np.argmin(abs(psi-1))
Bpol_sep = Bpol[isep]
print("Bpol_sep",Bpol_sep)
pint = np.sum(press)/len(press)
mu0 = 1.2566e-6
beta_pol = pint/(Bpol_sep**2/mu0/2.0)
print("beta_pol",beta_pol)






pdict = {}

pdict['entry_length'] = entry_length
pdict = add_to_pdict(pdict,psi,ne/1e20,'ne(10^20/m^3)')
pdict = add_to_pdict(pdict,psi,ni/1e20,'ni(10^20/m^3)')
pdict = add_to_pdict(pdict,psi,nz/1e20,'nz1(10^20/m^3)')
pdict = add_to_pdict(pdict,psi,Tikev,'ti(KeV)')
pdict = add_to_pdict(pdict,psi,Tekev,'te(KeV)')
ptot_kpa = (TeJ*ne + TiJ*ni)/1000.0
pdict = add_to_pdict(pdict,psi,ptot_kpa,'ptot(KPa)')
quants,grads = get_lists()
zeros = np.zeros(entry_length)
for i in quants:
    if i != 'ne(10^20/m^3)' and i != 'ni(10^20/m^3)' and i != 'ti(KeV)' and i != 'te(KeV)' and i != 'ptot(KPa)' and i != 'nz1(10^20/m^3)' and i != 'nb(10^20/m^3)' and i!= 'pb(KPa)':
        pdict = add_to_pdict(pdict,psi,zeros,i)
pdict['species info'] = '2 N Z A of ION SPECIES \n 6.000000   6.000000   12.00000\n 1.000000   1.000000   2.000000'
pdict['file_name'] = 'p'+gfile[1:]+'_T0_'+str(T0)



write_pfile(pdict)

pdict_in = read_pfile_direct(pdict['file_name']+'new')

fig = plt.figure(figsize = (6.0,8.0))
nstring = 'ne(10^20/m^3)'
Tstring = 'te(KeV)'
#plt.suptitle(r'$Q_{NL}/Q_{QLstd}$')
plt.subplot2grid((3,2),(0,0))
plt.plot(Bfile[:,1],Bfile[:,4])
plt.xlabel(r'$\Psi_N$')
plt.ylabel(r'$q$')
#ax = plt.axis()
#iq95 = np.argmin(abs(Bfile[:,1]-0.95))
#ymax = max(ax[3],
#plt.axis([0,1,ax[2],1.1*Bfile[iq95,4]])
plt.subplot2grid((3,2),(0,1))
plt.title(r'$\beta_{pol}$'+'(RFP) = '+str(beta_pol)[0:5])
plt.plot(Bfile[:,1],Bfile[:,5],label='P from efit')
plt.plot(pdict_in['psinorm_'+nstring],2*pdict_in[nstring]*1e20*pdict_in[Tstring]*1000*e,label='P from pfile (2x ne*Te)')
plt.xlabel(r'$\Psi_N$')
plt.ylabel(r'$P$'+'(pasc)')
plt.legend()
plt.subplot2grid((3,2),(1,0))
plt.plot(Bfile[:,1],Bfile[:,3])
plt.xlabel(r'$\Psi_N$')
plt.ylabel(r'$B_\phi(T)$')
plt.subplot2grid((3,2),(1,1))
plt.plot(Bfile[:,1],Bfile[:,2])
plt.xlabel(r'$\Psi_N$')
plt.ylabel(r'$B_\theta(T)$')
plt.xlabel(r'$\Psi_N$')
plt.subplot2grid((3,2),(2,0))
plt.plot(pdict_in['psinorm_'+nstring],pdict_in[nstring],label=nstring)
plt.xlabel(r'$\Psi_N$')
plt.ylabel(r'$n_e(10^{20}m^{-3})$')
plt.subplot2grid((3,2),(2,1))
plt.plot(pdict_in['psinorm_'+Tstring],pdict_in[Tstring],label=Tstring)
plt.xlabel(r'$\Psi_N$')
plt.ylabel(r'$T_e(KeV)$')
plt.tight_layout()
plt.show()

