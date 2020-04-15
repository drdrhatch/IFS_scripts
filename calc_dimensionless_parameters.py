#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import matplotlib.pyplot as plt
import optparse as op
import numpy as np
from read_EFIT_file import *

parser=op.OptionParser(description='Calculates dimensionless parameters (beta, nu*, rho*) from an EFIT file and gene profile files at a given value of rho_tor.  Arguments: efit_file_name, gene_profiles_file_name_e, gene_file_name_i, rhot0')
options,args=parser.parse_args()
if len(args)!=4:
    exit("""
Please include efit_file_name, gene_profiles_file_name_e, gene_profiles_file_name_i, rhot0."
    \n""")

efit_file_name = args[0]
gene_profiles_e = args[1]
gene_profiles_i = args[2]
rhot0 = float(args[3])

mi = 1.673e-27
me = 9.109e-31
ee = 1.602e-19
mref = 2.0*mi

Lref, Bref, R_major, q0 = get_dimpar_pars(efit_file_name,rhot0)
profe = np.genfromtxt(gene_profiles_e)
profi = np.genfromtxt(gene_profiles_i)

irhot0e = np.argmin(abs(profe[:,0]-rhot0))
irhot0i = np.argmin(abs(profi[:,0]-rhot0))

Te = profe[irhot0e,2]
Ti = profi[irhot0i,2]
ne = profe[irhot0e,3]
ni = profi[irhot0i,3]

print("Tref:",Te)
print("nref:",ne)
print("Bref:",Bref)
print("Lref:",Lref)
minor_r = 1.0

trpeps = rhot0*Lref/R_major

dummy = input("Assuming ion charge is 1 and reference mass is deuterium (press any key to continue).\n")
Z = 1.0

beta = 403.0e-5*ne*Te/Bref**2  #From GENE documentation
crefSI = (Te*1000.0*ee/mref)**0.5
cref_gene = 9787.1518*np.sqrt((Te*1000.0)/2.0)
OmrefSI = ee*Bref/mref 
rhostar = crefSI/OmrefSI/Lref
rhostar_gene = 3.2255E-3 * np.sqrt((2.0)*Te)/Bref/(minor_r*Lref)
coll = 2.3031E-5*Lref*(ne)/(Te)**2*(24.0-log(sqrt(ne*1.0E13)/Te*0.001))

nustar_i=8.0/3.0/pi**0.5*q0/trpeps**1.5*(R_major/Lref)*(ni/ne)*Z**4/(Te/Ti)**2*coll

nustar_e=16.0/3.0/pi**0.5*q0/trpeps**1.5*(R_major/Lref)*Z**2*coll

print("Sound speed (m/s):", crefSI)
print("Sound gyroradius (m):", crefSI/OmrefSI)
print("Transit frequency (1/s):", crefSI/Lref)
print("QGBref = nref Tref vtref rhoref*^2 translated to (W/m^2):", ne*1.0e19*Te*1000.0*ee*crefSI*rhostar**2)
print("GammaGBref = nref vtref rhoref*^2 translated to (particles / m^2 / s):", ne*1.0e19*crefSI*rhostar**2)
print("ChiGBref = Lref vtref rhoref*^2 translated to (m^2/s):", Lref*crefSI*rhostar**2)
print("QGBe = nref Tref vte rhoe*^2 translated to (W/m^2):", (ne*1.0e19*Te*1000.0*ee*crefSI*rhostar**2)*(me/mref)**0.5)
print("coll (from GENE)",coll)

print("*************************")
print("Dimensionless parameters:")
print("*************************")
print("rhostar = ",rhostar)
print("nustar_i",nustar_i)
print("nustar_e",nustar_e)
print("q0",q0)
print("beta = ",beta)

