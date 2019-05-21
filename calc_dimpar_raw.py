#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import matplotlib.pyplot as plt
import numpy as np

mi = 1.673e-27
me = 9.109e-31
ee = 1.602e-19

Te = 0.3 #keV
Ti = Te
ne = 3.0 #10^19 m^-3
ni = ne
mref = 2.0*mi
Bref = 1.2 #T
Lref = 0.07 #m minor radius
trpeps = 0.3
q0 = 2.0
R_major = Lref/trpeps
Z = 1

beta = 403.0e-5*ne*Te/Bref**2  #From GENE documentation
crefSI = (Te*1000.0*ee/mref)**0.5
cref_gene = 9787.1518*np.sqrt((Te*1000.0)/2.0)
OmrefSI = ee*Bref/mref 
rhostar = crefSI/OmrefSI/Lref
rhostar_gene = 3.2255E-3 * np.sqrt((2.0)*Te)/Bref/(Lref)
coll = 2.3031E-5*Lref*(ne)/(Te)**2*(24.0-np.log(np.sqrt(ne*1.0E13)/Te*0.001))

nustar_i=8.0/3.0/np.pi**0.5*q0/trpeps**1.5*(R_major/Lref)*(ni/ne)*Z**4/(Te/Ti)**2*coll

nustar_e=16.0/3.0/np.pi**0.5*q0/trpeps**1.5*(R_major/Lref)*Z**2*coll

print "Sound speed (m/s):", crefSI
print "Sound gyroradius (m):", crefSI/OmrefSI
print "Transit frequency (1/s):", crefSI/Lref
print "QGBref = nref Tref vtref rhoref*^2 (W/m^2):", ne*1.0e19*Te*1000.0*ee*crefSI*rhostar**2
print "ChiGBref = Lref vtref rhoref*^2 (m^2/s):", Lref*crefSI*rhostar**2
print "ChiBref = Lref vtref rhoref*^2 (m^2/s):", Lref*crefSI*rhostar
print "gyroBohm confinement time: Lref^2 / ChiGBref (s):", Lref**2/(Lref*crefSI*rhostar**2)
print "Bohm confinement time: Lref^2 / ChiBref (s):", Lref**2/(Lref*crefSI*rhostar)
print "QGBe = nref Tref vte rhoe*^2 (W/m^2):", (ne*1.0e19*Te*1000.0*ee*crefSI*rhostar**2)*(me/mref)**0.5

print "*************************"
print "Dimensionless parameters:"
print "*************************"
print "rhostar = ",rhostar
print "coll",coll
print "nustar_i",nustar_i
print "nustar_e",nustar_e
print "q0",q0
print "beta = ",beta

