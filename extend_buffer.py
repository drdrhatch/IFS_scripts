#!/usr/bin/env python
# -*- coding: utf-8 -*-
import optparse as op
import numpy as np
import matplotlib.pyplot as plt
from read_iterdb_file import *
from read_write_geometry import *
from write_iterdb import *
from ParIO import *
from interp import *
from finite_differences import *
from read_EFIT_file import *

parser=op.OptionParser(description='Extends profiles and geometry past separatrix.  Arguments: iterdb_file_name, geometry_file_name, parameters_file, efit_file_name.')
parser.add_option('--xbound','-x',type = 'float',action='store',dest="xbound",help = 'Rho_tor at which to begin smoothing (default: maximum rho_tor of previous simulation).',default=-1)
parser.add_option('--impurity','-i',action='store_const',const=1,help = 'Include impurity',default = False)
parser.add_option('--qtreatment','-q',action='store_const',const=1,help = '0: constant q, 1: constant gradient smoothed to 0',default = 1)
options,args=parser.parse_args()
if len(args)!=4:
    exit("""
Please include names of iterdb file, geometry file, parameters file, and EFIT file."
    \n""")

iterdb_file = args[0]
geom_file = args[1]
parameters_file = args[2]
efit_file = args[3]
xbound = options.xbound
impurity = options.impurity
qtreatment = options.qtreatment

geompars, geom = read_geometry_global(geom_file)
print("np.shape(geom['gzz'])",np.shape(geom['gzz']))

if impurity:
    rhot_te, te, ti, ne, ni, nb, vrot = read_iterdb_file(iterdb_file)
else:
    rhot_te, te, ti, ne, ni, vrot = read_iterdb_file(iterdb_file)

print("len(rhot_te)",len(rhot_te))
#for i in range(len(rhot_te)):
#    print "i, rhot",i,rhot_te[i]

par = Parameters()
par.Read_Pars(parameters_file)
pars = par.pardict

rhotor_min = pars['x0']-pars['lx_a']/2.0
rhotor_max = pars['x0']+pars['lx_a']/2.0
rhot0 = np.linspace(rhotor_min,rhotor_max,pars['nx0']) 
if xbound == -1:
   xbound = rhotor_max

print("This simulation is centered at x0="+str(pars['x0'])+" with lx_a ="+str(pars['lx_a'])+" and uses "+str(pars['nx0'])+" grid points.")
print("rho_tor = ["+str(pars['x0']-pars['lx_a']/2.0)+" , "+str(pars['x0']+pars['lx_a']/2.0)+"]")

extra_nx0 = int(float(input("How many gridpoints would you like to append to the simulation?\n")))
print(str(extra_nx0)+" additional gridpoints.")
drhotor = pars['lx_a']/(pars['nx0']-1)
lx_extend = drhotor*extra_nx0
erhotor_max = rhotor_max+lx_extend
elx_a = erhotor_max - rhotor_min
enx0 = int(pars['nx0']+extra_nx0)
print("New simulation range: rho_tor = [" + str(rhotor_min) + "," + str(rhotor_max+lx_extend) +"]")
print("New lx_a:", elx_a)
print("New nx0:", enx0)
ex0 = rhotor_min + (erhotor_max - rhotor_min)/2.0
print("New x0:", ex0)

print("Starting profile smoothing at xbound = ",xbound)

if impurity:
   rhot_te, te, ti, ne, ni, nb, vrot = read_iterdb_file(iterdb_file)
else:
   rhot_te, te, ti, ne, ni, vrot = read_iterdb_file(iterdb_file)

te = te/1000.0
ti = ti/1000.0
ne = ne/1.0e19
ni = ni/1.0e19
if impurity:
   nb = nb/1.0e19

erhot = np.linspace(rhotor_min,erhotor_max,enx0)
#itemp = np.argmin(abs(erhot-ex0))
#ex0 = erhot[itemp]
#print "New x0:", ex0

ete = interp(rhot_te,te,erhot)
eti = interp(rhot_te,ti,erhot)
ene = interp(rhot_te,ne,erhot)
eni = interp(rhot_te,ni,erhot)
if impurity:
   enb = interp(rhot_te,nb,erhot)
evrot = interp(rhot_te,vrot,erhot)

show_plots = False
if show_plots:
   plt.plot(rhot_te,te)
   plt.plot(erhot,ete)
   plt.show()
   plt.plot(rhot_te,ti)
   plt.plot(erhot,eti)
   plt.show()
   plt.plot(rhot_te,ne)
   plt.plot(erhot,ene)
   plt.show()
   plt.plot(rhot_te,ni)
   plt.plot(erhot,eni)
   plt.show()
   if impurity:
      plt.plot(rhot_te,nb)
      plt.plot(erhot,enb)
      plt.show()
   plt.plot(rhot_te,vrot)
   plt.plot(erhot,evrot)
   plt.show()

ixbound = np.argmin(abs(erhot-xbound))
print('ixbound =', ixbound)
lambdaT = 8.0/(erhotor_max-rhotor_max)
lambdaN = 6.0/(erhotor_max-rhotor_max)
print("lambdaT",lambdaT)
omte = -1.0/ete*fd_d1_o4(ete,erhot)
omti = -1.0/eti*fd_d1_o4(eti,erhot)
omne = -1.0/ene*fd_d1_o4(ene,erhot)
omni = -1.0/eni*fd_d1_o4(eni,erhot)
if impurity:
   omnb = -1.0/enb*fd_d1_o4(enb,erhot)
domega = -1.0*fd_d1_o4(evrot,erhot)

for i in range(enx0 - ixbound):
   #print ete[ixbound-1]*np.e**(-omte[ixbound]/lambdaT*np.e**(-lambdaT*(erhot[ixbound+i]-xbound)))
   #print "erhot[ixbound+i],xbound",erhot[ixbound+i],xbound
   #print "np.e**(-lambdaT*(erhot[ixbound+i]-xbound))",np.e**(-lambdaT*(erhot[ixbound+i]-xbound))
   ete[ixbound+i] = ete[ixbound]*np.e**(-omte[ixbound]/lambdaT)*np.e**(omte[ixbound]/lambdaT*np.e**(-lambdaT*(erhot[ixbound+i]-erhot[ixbound])))
   eti[ixbound+i] = eti[ixbound]*np.e**(-omti[ixbound]/lambdaT)*np.e**(omti[ixbound]/lambdaT*np.e**(-lambdaT*(erhot[ixbound+i]-erhot[ixbound])))
   ene[ixbound+i] = ene[ixbound]*np.e**(-omne[ixbound]/lambdaN)*np.e**(omne[ixbound]/lambdaN*np.e**(-lambdaN*(erhot[ixbound+i]-erhot[ixbound])))
   eni[ixbound+i] = eni[ixbound]*np.e**(-omni[ixbound]/lambdaN)*np.e**(omni[ixbound]/lambdaN*np.e**(-lambdaN*(erhot[ixbound+i]-erhot[ixbound])))
   if impurity:
      enb[ixbound+i] = enb[ixbound]*np.e**(-omnb[ixbound]/lambdaN)*np.e**(omnb[ixbound]/lambdaN*np.e**(-lambdaN*(erhot[ixbound+i]-erhot[ixbound])))
   evrot[ixbound+i] = domega[ixbound]/lambdaT*np.e**(-lambdaT*(erhot[ixbound+i]-erhot[ixbound])) + evrot[ixbound]-domega[ixbound]/lambdaT 

plt.plot(rhot_te,te,label='te')
plt.plot(erhot,ete,label='ete')
plt.axis((erhot[0],erhot[-1],0.0*min(ete),1.1*max(ete)))
plt.legend()
plt.show()

plt.plot(rhot_te,ti,label='ti')
plt.plot(erhot,eti,label='eti')
ax=plt.axis()
plt.axis((erhot[0],erhot[-1],0.0*min(eti),1.1*max(eti)))
plt.legend()
plt.show()

plt.plot(rhot_te,ne,label='ne')
plt.plot(erhot,ene,label='ene')
ax=plt.axis()
plt.axis((erhot[0],erhot[-1],0.0*min(ene),1.1*max(ene)))
plt.legend()
plt.show()

plt.plot(rhot_te,ni,label='ni')
plt.plot(erhot,eni,label='eni')
ax=plt.axis()
plt.axis((erhot[0],erhot[-1],0.0*min(eni),1.1*max(eni)))
plt.legend()
plt.show()

plt.plot(rhot_te,vrot,label='vrot')
plt.plot(erhot,evrot,label='evrot')
ax=plt.axis()
plt.axis((erhot[0],erhot[-1],0.0*min(evrot),1.1*max(evrot)))
plt.legend()
plt.show()

file_base =  iterdb_file[:-7]+"_extended"+"_nx0_"+str(enx0)+"_x0_"+str(ex0)+"_lx_a_"+str(elx_a)
if impurity:
   output_iterdb(erhot,erhot,ene,ete,eni,eti,file_base,'0','999.9',vrot=evrot,nimp=enb)
else:
   output_iterdb(erhot,erhot,ene,ete,eni,eti,file_base,'0','999.9',vrot=evrot)


Lref, Bref, R_major, q0, shat0 = get_geom_pars(efit_file,ex0)

iex0 = np.argmin(abs(erhot-ex0))
beta0 = 403.0e-5*ene[iex0]*ete[iex0]/Bref**2

geompars_new={}
geompars_new['gridpoints'] = geompars['gridpoints']
geompars_new['q0'] = q0
geompars_new['shat'] = shat0
geompars_new['s0'] = ex0**2
geompars_new['minor_r'] = geompars['minor_r']
geompars_new['major_R'] = geompars['major_R']
if geompars['trpeps'] == 0.0:
    print("Check trpeps!!")
    dummy = input("Press a key")
geompars_new['trpeps'] = geompars['trpeps']
geompars_new['beta'] = beta0
geompars_new['Lref'] = geompars['Lref']
geompars_new['Bref'] = geompars['Bref']
geompars_new['magn_geometry'] = geompars['magn_geometry']

geom_file_out = geom_file + "_extended"+"_enx0_"+str(enx0)+"_x0_"+str(ex0)

geom_new = {}
#geom_dummy = np.empty((geompars['gridpoints'],enx0))
geom_dummy1D = np.empty((enx0))

#print np.shape(geom_dummy)
#print np.shape(geom['gxx'])
#print "len(geom_dummy[:,0]",len(geom_dummy[:,0])
geompars['gridpoints'] = int(geompars['gridpoints'])
geom_new['gxx'] = np.empty((geompars['gridpoints'],enx0))
geom_new['gxx'][:,0:pars['nx0']] = geom['gxx']
for i in range(extra_nx0):
    geom_new['gxx'][:,pars['nx0']+i] = geom['gxx'][:,-1]

geom_new['gxy'] = np.empty((geompars['gridpoints'],enx0))
geom_new['gxy'][:,0:pars['nx0']] = geom['gxy']
for i in range(extra_nx0):
    geom_new['gxy'][:,pars['nx0']+i] = geom['gxy'][:,-1]

geom_new['gxz'] = np.empty((geompars['gridpoints'],enx0))
geom_new['gxz'][:,0:pars['nx0']] = geom['gxz']
for i in range(extra_nx0):
    geom_new['gxz'][:,pars['nx0']+i] = geom['gxz'][:,-1]

geom_new['gyy'] = np.empty((geompars['gridpoints'],enx0))
geom_new['gyy'][:,0:pars['nx0']] = geom['gyy']
for i in range(extra_nx0):
    geom_new['gyy'][:,pars['nx0']+i] = geom['gyy'][:,-1]

geom_new['gyz'] = np.empty((geompars['gridpoints'],enx0))
geom_new['gyz'][:,0:pars['nx0']] = geom['gyz']
for i in range(extra_nx0):
    geom_new['gyz'][:,pars['nx0']+i] = geom['gyz'][:,-1]

geom_new['gzz'] = np.empty((geompars['gridpoints'],enx0))
geom_new['gzz'][:,0:pars['nx0']] = geom['gzz']
for i in range(extra_nx0):
    geom_new['gzz'][:,pars['nx0']+i] = geom['gzz'][:,-1]

geom_new['Bfield'] = np.empty((geompars['gridpoints'],enx0))
geom_new['Bfield'][:,0:pars['nx0']] = geom['Bfield']
for i in range(extra_nx0):
    geom_new['Bfield'][:,pars['nx0']+i] = geom['Bfield'][:,-1]

geom_new['dBdy'] = np.empty((geompars['gridpoints'],enx0))
geom_new['dBdy'][:,0:pars['nx0']] = geom['dBdy']
for i in range(extra_nx0):
    geom_new['dBdy'][:,pars['nx0']+i] = geom['dBdy'][:,-1]

geom_new['dBdz'] = np.empty((geompars['gridpoints'],enx0))
geom_new['dBdz'][:,0:pars['nx0']] = geom['dBdz']
for i in range(extra_nx0):
    geom_new['dBdz'][:,pars['nx0']+i] = geom['dBdz'][:,-1]

geom_new['dBdx'] = np.empty((geompars['gridpoints'],enx0))
geom_new['dBdx'][:,0:pars['nx0']] = geom['dBdx']
for i in range(extra_nx0):
    geom_new['dBdx'][:,pars['nx0']+i] = geom['dBdx'][:,-1]

geom_new['jacobian'] = np.empty((geompars['gridpoints'],enx0))
geom_new['jacobian'][:,0:pars['nx0']] = geom['jacobian']
for i in range(extra_nx0):
    geom_new['jacobian'][:,pars['nx0']+i] = geom['jacobian'][:,-1]

C_y0 = ex0/q0
geom_new['C_y'] = np.empty((enx0))
geom_new['C_y'][0:pars['nx0']] = C_y0
for i in range(extra_nx0):
    geom_new['C_y'][pars['nx0']+i] = C_y0

geom_new['C_xy'] = np.empty((enx0))
geom_new['C_xy'][0:pars['nx0']] = geom['C_xy']
for i in range(extra_nx0):
    geom_new['C_xy'][pars['nx0']+i] = geom['C_xy'][-1]

geom_new['geo_R'] = np.empty((geompars['gridpoints'],enx0))
geom_new['geo_R'][:,0:pars['nx0']] = geom['geo_R']
for i in range(extra_nx0):
    geom_new['geo_R'][:,pars['nx0']+i] = geom['geo_R'][:,-1]

geom_new['geo_Z'] = np.empty((geompars['gridpoints'],enx0))
geom_new['geo_Z'][:,0:pars['nx0']] = geom['geo_Z']
for i in range(extra_nx0):
    geom_new['geo_Z'][:,pars['nx0']+i] = geom['geo_Z'][:,-1]

geom_new['geo_c1'] = np.empty((geompars['gridpoints'],enx0))
geom_new['geo_c1'][:,0:pars['nx0']] = geom['geo_c1']
for i in range(extra_nx0):
    geom_new['geo_c1'][:,pars['nx0']+i] = geom['geo_c1'][:,-1]

geom_new['geo_c2'] = np.empty((geompars['gridpoints'],enx0))
geom_new['geo_c2'][:,0:pars['nx0']] = geom['geo_c2']
for i in range(extra_nx0):
    geom_new['geo_c2'][:,pars['nx0']+i] = geom['geo_c2'][:,-1]

geom_new['dpdx_pm_arr'] = np.empty((enx0))
geom_new['dpdx_pm_arr'][0:pars['nx0']] = geom['dpdx_pm_arr']
for i in range(extra_nx0):
    geom_new['dpdx_pm_arr'][pars['nx0']+i] = geom['dpdx_pm_arr'][-1]

geom_new['q'] = np.empty((enx0))
geom_new['q'][0:pars['nx0']] = geom['q']
if qtreatment == 0:
    for i in range(extra_nx0):
        geom_new['q'][pars['nx0']+i] = geom['q'][-1]
elif qtreatment == 1:
    qprime = geom['q'][-1]-geom['q'][-2]
    for i in range(extra_nx0):
        geom_new['q'][pars['nx0']+i] = geom_new['q'][pars['nx0']+i-1] + qprime*np.e**(-float(i)/float(extra_nx0)*3.0)

plt.plot(geom_new['q'])
plt.title('q')
plt.show()

print("Writing new geometry file.")
write_tracer_efit_file(geompars_new,geom_new,geom_file_out)
print("Done.")

plt.contourf(geom_new['gxz'])
plt.title('gxz new (end)')
plt.show()

