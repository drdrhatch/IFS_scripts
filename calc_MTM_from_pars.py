#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import matplotlib.pyplot as plt
import optparse as op
import numpy as np
from read_EFIT_file import *
from subprocess import call
from ParIO import * 
from read_write_geometry import *
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
geomfile = pars['magn_geometry'][1:-1]+'_'+suffix

gpars, geom = read_geometry_global(geomfile)

mi = 1.673e-27
me = 9.109e-31
ee = 1.602e-19
mref = 2.0*mi

gene_profiles_e = 'profiles_e_'+suffix
gene_profiles_i = 'profiles_i_'+suffix

profe = np.genfromtxt(gene_profiles_e)
profi = np.genfromtxt(gene_profiles_i)

omte = profe[:,4]
omne = profe[:,5]
omti = profi[:,4]
omni = profi[:,5]
Teprof = profe[:,2]
neprof = profe[:,3]
xgrid = profe[:,0]

x0 = pars['x0']
lx_a = pars['lx_a']
qprof = geom['q']
dummy = raw_input("Assuming ions are first species, ion charge is 1 and reference mass is deuterium (press any key to continue).\n")

Te = pars['Tref']
Ti = Te*pars['temp1']
ne = pars['nref']
ni = ne*pars['dens1']
Bref = pars['Bref']
Lref = pars['Lref']


beta = 403.0e-5*ne*Te/Bref**2  #From GENE documentation
crefSI = (Te*1000.0*ee/mref)**0.5
cref_gene = 9787.1518*np.sqrt((Te*1000.0)/2.0)
OmrefSI = ee*Bref/mref 
rhostar = crefSI/OmrefSI/Lref
#coll0 = 2.3031E-5*Lref*(ne)/(Te)**2*(24.0-log(sqrt(ne*1.0E13)/Te*0.001))
#coll0 = 2.3031e-5*(24.0-np.log(sqrt(ne*1.0E13)/Te/1000.0))*Lref*ne/Te**2
coll0 = 2.3031E-5*pars['Lref']*(pars['nref'])/(pars['Tref'])**2*(24.-np.log(np.sqrt(pars['nref']*1.0E13)/pars['Tref']*0.001))

coll_prof = 2.3031E-5*Lref*(neprof)/(Teprof)**2*(24.0-np.log(sqrt(neprof*1.0E13)/Teprof*0.001))
coll = pars['coll']
nue = 4.0*coll*(mref/me)**0.5
nue_prof = 4.0*coll_prof*(mref/me)**0.5
print "nue",nue
print "nue_prof[x0]",nue_prof[int(pars['nx0']/2.0)]


ixmin = int(pars['nx0']*pars['l_buffer_size'])
ixmax = int(pars['nx0']-pars['nx0']*pars['u_buffer_size'])
print "ixmin",ixmin
print "ixmax",ixmax
omte_avg = np.sum(omte[ixmin:ixmax])/(ixmax-ixmin)
omne_avg = np.sum(omne[ixmin:ixmax])/(ixmax-ixmin)
omte_max = np.max(omte[ixmin:ixmax])
omne_max = np.max(omne[ixmin:ixmax])
shat = fd_d1_o4(qprof,xgrid)/qprof
shat_avg = np.sum(shat[ixmin:ixmax])/(ixmax-ixmin)
shat_min = np.min(shat[ixmin:ixmax])
n0 = pars['n0_global']
#Assuming negligible variation of Bref
rhostar_prof =  np.sqrt(1000.0*ee*Teprof/mref)*mref/ee/Bref/Lref
print "pars['rhostar']",pars['rhostar']
print "rhostar",rhostar_prof[int(pars['nx0']/2.0)]

#Calculate radial dependence of omega_star for given n0
omega_star_prof = n0*qprof/xgrid*rhostar_prof*(omte+omne)
omega_star_prof_kHz = omega_star_prof*np.sqrt(1000.0*ee*Teprof/mref)/Lref/1000.0/2.0/np.pi
omega_star_prof_locnorm = omega_star_prof_kHz*1000.0*Lref/crefSI*2.0*np.pi

betaprof = 403.0e-5*neprof*Teprof/Bref**2  
beta_avg = np.sum(betaprof[ixmin:ixmax])/(ixmax-ixmin)

if idb_file:
   if pars['n_spec'] > 2:
      rhotidb, teidb, tiidb, neidb, niidb, nzidb, vrot0 = read_iterdb_file(idb_file)
   else:
      rhotidb, teidb, tiidb, neidb, niidb, vrot0 = read_iterdb_file(idb_file)
   vrot_u = interp(rhotidb,vrot0,xgrid)
   omegaDoppler = vrot_u*n0/2.0/np.pi/1000.0
   omegaDoppler_norm = omegaDoppler*1000.0*Lref/crefSI*2.0*np.pi


plt.plot(xgrid[ixmin:ixmax],omega_star_prof[ixmin:ixmax])
plt.title('omega_star profile n0='+str(n0))
plt.show()
plt.plot(xgrid[ixmin:ixmax],omega_star_prof_kHz[ixmin:ixmax])
if idb_file:
   plt.plot(xgrid[ixmin:ixmax],omegaDoppler[ixmin:ixmax])
   plt.plot(xgrid[ixmin:ixmax],omegaDoppler[ixmin:ixmax]+omega_star_prof_kHz[ixmin:ixmax])
plt.title('omega_star profile kHz n0='+str(n0))
plt.show()
plt.plot(xgrid[ixmin:ixmax],omega_star_prof_locnorm[ixmin:ixmax])
if idb_file:
   plt.plot(xgrid[ixmin:ixmax],omegaDoppler_norm[ixmin:ixmax])
   plt.plot(xgrid[ixmin:ixmax],omegaDoppler_norm[ixmin:ixmax]+omega_star_prof_locnorm[ixmin:ixmax])
plt.title('omega_star profile locnorm n0='+str(n0))
plt.show()

ix_max_oms = np.argmax(omega_star_prof_locnorm)
omte_omne_max_oms = omega_star_prof_locnorm[ix_max_oms]/pars['kymin']
nue_max_oms = nue_prof[ix_max_oms]
shat_max_oms = shat[ix_max_oms]
beta_max_oms = betaprof[ix_max_oms]

plt.plot(xgrid[ixmin:ixmax],qprof[ixmin:ixmax])
plt.show()

plt.plot(xgrid[ixmin:ixmax],shat[ixmin:ixmax])
plt.title("shat_avg,min "+str(shat_avg)+' , '+str(shat_min))
plt.show()

plt.plot(xgrid[ixmin:ixmax],omte[ixmin:ixmax])
plt.title("omte_avg "+str(omte_avg))
plt.show()

print "betae",beta
print "betae max oms",beta_max_oms
print "coll (from GENE)",coll
print "coll0",coll0
print "coll0/coll",coll0/coll
print "Avg shat",shat_avg
print "Min shat",shat_min
print "Avg a/LTe",omte_avg
beta_hat_maxoms = beta_max_oms*(omte_omne_max_oms/shat_max_oms)**2
print "beta_hat avg = betae*(Ls/LTe)^2",beta_avg*(omte_avg/shat_avg)**2
print "beta_hat max oms = beta_max_oms*(omte_omne_max_oms/shat_max_oms)^2", beta_hat_maxoms
print "nu_e (normalized to cs/a)", nue

kygrid = np.linspace(0.0,1.0,num=10)
plt.plot(kygrid,kygrid*omte_avg,label='omega* = ky rhos a/LTe)')
#plt.plot(kygrid,kygrid*(omte_max+omne_max),label='omega* max')
plt.hlines(2.0*nue,0.0,1.0,color='black',label="2 * nu_e (cs/a)")
plt.hlines(10.0*nue,0.0,1.0,color='black',label="10 * nu_e (cs/a)")
plt.hlines(0.4*nue,0.0,1.0,color='black',label="0.4 * nu_e (cs/a)")
plt.xlabel('ky rhos')
plt.ylabel('omega a/cs')
plt.legend()
plt.show()

plt.plot(kygrid,kygrid*omte_omne_max_oms,label='omega* = ky rhos( a/LTe+a/Ln)')
#plt.plot(kygrid,kygrid*(omte_max+omne_max),label='omega* max')
plt.hlines(2.0*nue_max_oms,0.0,1.0,color='black',label="2 * nu_e (cs/a)")
plt.hlines(10.0*nue_max_oms,0.0,1.0,color='black',label="10 * nu_e (cs/a)")
plt.hlines(0.4*nue_max_oms,0.0,1.0,color='black',label="0.4 * nu_e (cs/a)")
plt.xlabel('ky rhos')
plt.ylabel('omega a/cs')
plt.title('quantities at maximum omega_star')
plt.legend()
plt.show()

ngrid = np.arange(30)
plt.plot(ngrid,(ngrid/float(n0)*pars['kymin'])*omte_omne_max_oms,label='omega* = ky rhos( a/LTe+a/Ln)')
#plt.plot(kygrid,kygrid*(omte_max+omne_max),label='omega* max')
plt.hlines(2.0*nue_max_oms,0.0,30.0,color='black',label="2 * nu_e (cs/a)")
plt.hlines(10.0*nue_max_oms,0.0,30.0,color='black',label="10 * nu_e (cs/a)")
plt.hlines(0.4*nue_max_oms,0.0,30.0,color='black',label="0.4 * nu_e (cs/a)")
plt.xlabel('toroidal n')
plt.ylabel('omega a/cs')
plt.title('quantities at maximum omega_star')
plt.legend()
plt.show()

f=open('MTM_nu_oms_info_'+suffix,'w')
f.write('#1.ky 2.n 3.Peak omega* \n')
f.write('#nue_max_oms: '+str(nue_max_oms)+'\n')
f.write("#beta_hat max oms = beta_max_oms*(omte_omne_max_oms/shat_max_oms)^2: "+str(beta_hat_maxoms) + '\n')
np.savetxt(f,np.column_stack((kygrid,(float(n0)/pars['kymin']*kygrid),kygrid*omte_omne_max_oms)))
f.close()




