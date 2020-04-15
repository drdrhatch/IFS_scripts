import numpy as np
import matplotlib.pyplot as plt
from interp import *
from ParIO import *
from read_iterdb import *
from finite_differences import *
from read_EFIT_file import *

parfile = 'parameters_1'
genefile = 'gene_profiles_e'
slicefile = 'slice_1.dat'
prof_file = 'profiles_i_1'
efit_file = 'g_temp'

ee = 1.602177e-19
mi = 1.672622e-27

par = Parameters()
par.Read_Pars(parfile)
pars = par.pardict

idbfile = pars['iterdb_file'][1:-1]
rhotidb,profidb,unitsidb = read_iterdb(idbfile)
rhotv = rhotidb['VROT']
vrot = profidb['VROT']
#plt.plot(rhotv,vrot)
#plt.show()

gp = np.genfromtxt(genefile)
pr = np.genfromtxt(prof_file)

dat = np.genfromtxt(slicefile)
phi = dat[:,1]

rhostar = pars['rhostar']
Tref = pars['Tref']

rhot = pr[:,0]

rhop = interp(gp[:,0],gp[:,1],rhot)
psi = rhop**2

psi0 = np.linspace(psi[0],psi[-1],500)
phi0 = interp(psi,phi,psi0)
rhot0 = interp(gp[:,1],gp[:,0],psi0**0.5)
#plt.plot(psi,phi)
#plt.plot(psi0,phi0)
#plt.show()

psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw,psiax,psisep = read_EFIT_file(efit_file)

psisep0 = psisep-psiax
print("psisep0",psisep0)
psiSI = psi0*psisep0

phiSI = phi0*Tref*1000.0*rhostar
vrot_ZF = fd_d1_o4(phiSI,psiSI)
plt.plot(psi0,vrot_ZF)
plt.show()

plt.plot(rhot0,vrot_ZF)
plt.plot(rhotv,vrot)
plt.xlabel('rhot')
plt.ylabel('Vrot (rad/s)')
plt.show()
plt.show()




