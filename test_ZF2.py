import numpy as np
import matplotlib.pyplot as plt
from interp import *
from ParIO import *
from read_iterdb import *
from finite_differences import *
from read_EFIT_file import *

#parfile = 'parameters_4'
#genefile = 'gene_profiles_e'
#slicefile = 'slice_4.dat'
#prof_file = 'profiles_i_4'
#efit_file = 'g_temp'
#rbsfile = 'rbsProfs'

parfile = 'parameters_1'
genefile = 'gene_profiles_e'
slicefile = 'slice_1.dat'
prof_file = 'profiles_i_1'
efit_file = 'g_temp'
rbsfile = 'rbsProfs'



ee = 1.602177e-19
mi = 1.672622e-27

par = Parameters()
par.Read_Pars(parfile)
pars = par.pardict

#idbfile = pars['iterdb_file'][1:-1]
#rhotidb,profidb,unitsidb = read_iterdb(idbfile)
#rhotv = rhotidb['VROT']
#vrot = profidb['VROT']
#plt.plot(rhotv,vrot)
#plt.show()

gp = np.genfromtxt(genefile)
pr = np.genfromtxt(prof_file)

dat = np.genfromtxt(slicefile)
phi = dat[:,1]
rhot = pr[:,0]

rhostar = pars['rhostar']
Tref = pars['Tref']

rbs = np.genfromtxt(rbsfile)
rhot_rbs = rbs[:,0]
R_rbs = rbs[:,24]

Er_rbs = rbs[:,16]
Er_o_RBp_rbs = rbs[:,27]
Bp_rbs = rbs[:,25]

#plt.plot(rhot_rbs,R_rbs)
#plt.title('R')
#plt.show()


#plt.plot(rhot_rbs,Er_rbs)
#plt.title('Er')
#plt.show()


#plt.plot(rhot_rbs,Er_o_RBp_rbs)
#plt.plot(rhot_rbs,Er_rbs/R_rbs/Bp)
#plt.title('Er/RBpol')
#plt.show()

ilow = np.argmin(abs(rhot_rbs-rhot[0]))
R0 = np.linspace(R_rbs[ilow],R_rbs[-1],2000)
Er0 = interp(R_rbs,Er_rbs,R0)
Bp0 = interp(R_rbs,Bp_rbs,R0)
rhot0 = interp(R_rbs,rhot_rbs,R0)

phinorm = rhostar*Tref*1000.0
phiSI = phi*phinorm
phiSI0 = full_interp(phiSI,rhot,rhot_rbs[ilow:],R_rbs[ilow:],R0,verify_interp=True)

Er_ZF = fd_d1_o4(phiSI0,R0)
vrot_ZF = Er_ZF / R0/Bp0

plt.plot(rhot0,vrot_ZF)
plt.plot(rhot_rbs,Er_o_RBp_rbs)
plt.title('Er/RBp')
plt.show()

plt.plot(rhot0,Er_ZF)
plt.plot(rhot_rbs,Er_rbs)
plt.title('Er')
plt.show()

np.savetxt('vrot_ZF2.dat',np.column_stack((rhot0,vrot_ZF)))

Lref = pars['Lref']
Bref = pars['Bref']
nref = pars['nref']*1.0e19
Tnorm = rhostar*Tref*1000.0*ee
kynorm = 1.0/rhostar/Lref
QGB = phinorm*Tnorm*kynorm*nref/Bref
print("QGB (watts/m^2):",QGB)



