#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import matplotlib.pyplot as plt
import optparse as op
import numpy as np
from read_EFIT_file import *
from subprocess import call

parser=op.OptionParser(description='Calculates nu_e, omega_*, beta_hat, etc. relevant for MTM physics. Arguments: efit_file_name, gene_profiles_file_name_e, gene_file_name_i, mid_ped rhot0, pedestal width')
options,args=parser.parse_args()
if len(args)!=5:
    exit("""
Please include efit_file_name, gene_profiles_file_name_e, gene_profiles_file_name_i, rhot0, ped_width."
    \n""")

efit_file_name = args[0]
gene_profiles_e = args[1]
gene_profiles_i = args[2]
rhot0 = float(args[3])
wped = float(args[4])
Binfo_file_name = 'Binfo_'+efit_file_name

mi = 1.673e-27
me = 9.109e-31
ee = 1.602e-19
mref = 2.0*mi

call(['my_efit_tools.py','-p','-n',efit_file_name])
binfo = np.genfromtxt(Binfo_file_name)
call(['calc_shat_from_efit.py',Binfo_file_name,gene_profiles_e,str(rhot0)])
shat_in = np.genfromtxt('shat.dat') #1.rhot, 2.q0, 3.shat
rhot_shat = shat_in[:,0]
shat = shat_in[:,2]
call(['calc_eta_from_gene_profiles.py',gene_profiles_e,gene_profiles_i,str(rhot0)])
profs = np.genfromtxt('profile_info_'+gene_profiles_e)
rhot_omte = profs[:,0]
omte = profs[:,4]

#Calculate average shat
imin_shat = np.argmin(abs(rhot_shat-(rhot0-wped/2.0)))
imax_shat = np.argmin(abs(rhot_shat-(rhot0+wped/2.0)))
shat_avg = np.sum(shat[imin_shat:imax_shat])/(imax_shat-imin_shat)
print("shat_avg",shat_avg)
plt.plot(rhot_shat[imin_shat:imax_shat],shat[imin_shat:imax_shat])
plt.title("shat_avg "+str(shat_avg))
plt.show()

#Calculate average omte
imin_omte = np.argmin(abs(rhot_omte-(rhot0-wped/2.0)))
imax_omte = np.argmin(abs(rhot_omte-(rhot0+wped/2.0)))
omte_avg = np.sum(omte[imin_omte:imax_omte])/(imax_omte-imin_omte)
print("omte_avg",omte_avg)
plt.plot(rhot_omte[imin_omte:imax_omte],omte[imin_omte:imax_omte])
plt.title("omte_avg "+str(omte_avg))
plt.show()

Lref, Bref, R_major, q0 = get_dimpar_pars(efit_file_name,rhot0)
profe = np.genfromtxt(gene_profiles_e)
profi = np.genfromtxt(gene_profiles_i)

irhot0e = np.argmin(abs(profe[:,0]-rhot0))
irhot0i = np.argmin(abs(profi[:,0]-rhot0))

Te = profe[irhot0e,2]
Ti = profi[irhot0i,2]
ne = profe[irhot0e,3]
ni = profi[irhot0i,3]

#print "Tref:",Te
#print "nref:",ne
#print "Bref:",Bref
#print "Lref:",Lref
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
nue = 4.0*coll*(mref/me)**0.5

#nustar_i=8.0/3.0/pi**0.5*q0/trpeps**1.5*(R_major/Lref)*(ni/ne)*Z**4/(Te/Ti)**2*coll
#nustar_e=16.0/3.0/pi**0.5*q0/trpeps**1.5*(R_major/Lref)*Z**2*coll

print("betae",beta)
print("coll",coll)
print("Avg a/Ls",shat_avg)
print("Avg a/LTe",omte_avg)
print("beta_hat = betae*(Ls/LTe)^2",beta*(omte_avg/shat_avg)**2)
print("nu_e (normalized to cs/a)", nue)

kygrid = np.linspace(0.0,1.0,num=10)
plt.plot(kygrid,kygrid*omte_avg,label='omega* = ky rhos a/LTe)')
plt.hlines(2.0*nue,0.0,1.0,color='black',label="2 * nu_e (cs/a)")
plt.hlines(10.0*nue,0.0,1.0,color='black',label="10 * nu_e (cs/a)")
plt.hlines(0.4*nue,0.0,1.0,color='black',label="0.4 * nu_e (cs/a)")
plt.xlabel('ky rhos')
plt.ylabel('omega a/cs')
plt.legend()
plt.show()




