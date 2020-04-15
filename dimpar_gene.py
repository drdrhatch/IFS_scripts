#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import matplotlib.pyplot as plt
import optparse as op
import numpy as np
from ParIO import * 

parser=op.OptionParser(description='Calculates dimensionless parameters (beta, nu*, rho*) from GENE parameters.  Argument: suffix from output parameters file from a GENE simulation.')
options,args=parser.parse_args()
if len(args)!=1:
    exit("""
Please include gene parameters file suffix."
    \n""")

suffix = args[0]

mi = 1.673e-27
me = 9.109e-31
ee = 1.602e-19

#Lref, Bref, R_major, q0 = get_dimpar_pars(efit_file_name,rhot0)
#profe = np.genfromtxt(gene_profiles_e)
#profi = np.genfromtxt(gene_profiles_i)

#irhot0e = np.argmin(abs(profe[:,0]-rhot0))
#irhot0i = np.argmin(abs(profi[:,0]-rhot0))
if '_' in suffix:
   suffix = suffix
elif 'dat' in suffix:
   suffix = '.dat'
else:
   suffix = '_'+suffix

parfile = 'parameters'+suffix
par = Parameters()
par.Read_Pars(parfile)
pars = par.pardict



Tref = pars['Tref']
Ti = pars['temp1']*Tref
nref = pars['nref']
ni = pars['dens1']*nref
Bref = pars['Bref']
Lref = pars['Lref']
mref = pars['mref']
imass = mref*mi
if 'q0' in pars:
   q0 = pars['q0']
else:
   q0 = float(input("Enter q0:"))

print(( "Tref:",Tref))
print(( "nref:",nref))
print(( "Bref:",Bref))
print(( "Lref:",Lref))
minor_r = 1.0
major_R = pars['major_R']*Lref

trpeps = pars['x0']*Lref/major_R
Z = 1.0

beta = 403.0e-5*nref*Tref/Bref**2  #From GENE documentation
crefSI = abs(Tref*1000.0*ee/imass)**0.5
cref_gene = 9787.1518*np.sqrt((Tref*1000.0)/2.0)
OmrefSI = ee*Bref/imass 
rhostar = crefSI/OmrefSI/Lref
rhostar_gene = 3.2255E-3 * np.sqrt((2.0)*Tref)/Bref/(minor_r*Lref)
coll = 2.3031E-5*Lref*(nref)/(Tref)**2*(24.0-np.log(np.sqrt(nref*1.0E13)/Tref*0.001))

#nustar_i=8.0/3.0/np.pi**0.5*q0/trpeps**1.5*(major_R/Lref)*(ni/nref)*Z**4/(Tref/Ti)**2*coll

nustar_e=16.0/3.0/np.pi**0.5*q0/abs(trpeps)**1.5*(major_R/Lref)*Z**2*coll

print(( "x0:", pars['x0']))
print(( "Sound speed (m/s):", crefSI))
print(( "Sound gyroradius (m):", crefSI/OmrefSI))
print(( "Transit frequency (1/s):", crefSI/Lref))
print(( "QGBref = nref Tref vtref rhoref*^2 translated to (W/m^2):", nref*1.0e19*Tref*1000.0*ee*crefSI*rhostar**2))
print(( "GammaGBref = nref vtref rhoref*^2 translated to (particles / m^2 / s):", nref*1.0e19*crefSI*rhostar**2))
print(( "ChiGBref = Lref vtref rhoref*^2 translated to (m^2/s):", Lref*crefSI*rhostar**2))
print(( "QGBe = nref Tref vte rhoe*^2 translated to (W/m^2):", (nref*1.0e19*Tref*1000.0*ee*crefSI*rhostar**2)*(me/imass)**0.5))
print(( "coll (from GENE)",coll))

print( "*************************")
print( "Dimensionless parameters:")
print( "*************************")
print(( "rhostar = ",rhostar))
#print "nustar_i",nustar_i
print(( "nustar_e",nustar_e))
print(( "q0",q0))
print(( "beta = ",beta))
for i in range(pars['n_spec']):
    print(( "name"+str(i+1),pars['name'+str(i+1)]))
    print(( "omn"+str(i+1),pars['omn'+str(i+1)]))
    print(( "omt"+str(i+1),pars['omt'+str(i+1)]))

