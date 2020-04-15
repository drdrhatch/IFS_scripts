#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import matplotlib.pyplot as plt
import optparse as op
import numpy as np
from read_EFIT_file import *
from subprocess import call
from ParIO import * 
from read_iterdb_file import *

parser=op.OptionParser(description='Calculates nu_e, omega_*, beta_hat, etc. relevant for MTM physics. Arguments: suffix.  Note: profiles files must be named profiles_e_<suffix> and profiles_i_<suffix>.')
parser.add_option('--idb','-i',type = 'str',action='store',dest="idb_file",help = 'ITERDB file name.',default='')

options,args=parser.parse_args()
if len(args)!=1:
    exit("""
Please include suffix."
    \n""")

suffix = args[0]
idb_file = options.idb_file

par = Parameters()
par.Read_Pars('parameters_'+suffix)
pars = par.pardict
print("pars",pars)

for i in range(3):
   if 'name'+str(i+1) in pars:
      print("i",i,str(i+1),pars['name'+str(i+1)][1:-1])
      if pars['name'+str(i+1)][1:-1] == 'e':
         specnum = i+1


mi = 1.673e-27
me = 9.109e-31
ee = 1.602e-19
mref = 2.0*mi



dummy = input("Assuming ion charge is 1 and reference mass is deuterium (press any key to continue).\n")

Te = pars['Tref']
Ti = Te*pars['temp'+str(specnum)]
ne = pars['nref']
ni = ne*pars['dens'+str(specnum)]
Bref = pars['Bref']
Lref = pars['Lref']
omte = pars['omt'+str(specnum)]
omne = pars['omn'+str(specnum)]

beta = 403.0e-5*ne*Te/Bref**2  #From GENE documentation
crefSI = (Te*1000.0*ee/mref)**0.5
cref_gene = 9787.1518*np.sqrt((Te*1000.0)/2.0)
OmrefSI = ee*Bref/mref 
rhostar = crefSI/OmrefSI/Lref
#coll0 = 2.3031E-5*Lref*(ne)/(Te)**2*(24.0-log(sqrt(ne*1.0E13)/Te*0.001))
#coll0 = 2.3031e-5*(24.0-np.log(sqrt(ne*1.0E13)/Te/1000.0))*Lref*ne/Te**2
coll0 = 2.3031E-5*pars['Lref']*(pars['nref'])/(pars['Tref'])**2*(24.-np.log(np.sqrt(pars['nref']*1.0E13)/pars['Tref']*0.001))
coll = pars['coll']
nue = 4.0*coll*(mref/me)**0.5

shat = pars['shat']
n0 = pars['n0_global']

print("betae",beta)
print("coll (from GENE)",coll)
print("coll0",coll0)
print("coll0/coll",coll0/coll)
print("shat",shat)
print("a/LTe",omte)
print("beta_hat = betae*(Ls/LTe)^2",beta*(omte/shat)**2)
print("nu_e (normalized to cs/a)", nue)

kygrid = np.linspace(0.0,1.0,num=10)
plt.plot(kygrid,kygrid*(omte+omne),label='omega* = ky rhos a/LTe)')
plt.hlines(2.0*nue,0.0,1.0,color='black',label="2 * nu_e (cs/a)")
plt.hlines(10.0*nue,0.0,1.0,color='black',label="10 * nu_e (cs/a)")
plt.hlines(0.4*nue,0.0,1.0,color='black',label="0.4 * nu_e (cs/a)")
plt.xlabel('ky rhos')
plt.ylabel('omega a/cs')
plt.legend()
plt.show()

ngrid = np.arange(30)
plt.plot(ngrid,(ngrid/float(n0)*pars['kymin'])*(omte+omne),label='omega* = ky rhos( a/LTe+a/Ln)')
#plt.plot(kygrid,kygrid*(omte_max+omne_max),label='omega* max')
plt.hlines(2.0*nue,0.0,30.0,color='black',label="2 * nu_e (cs/a)")
plt.hlines(10.0*nue,0.0,30.0,color='black',label="10 * nu_e (cs/a)")
plt.hlines(0.4*nue,0.0,30.0,color='black',label="0.4 * nu_e (cs/a)")
plt.xlabel('toroidal n')
plt.ylabel('omega a/cs')
plt.legend()
plt.show()






