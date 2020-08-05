from finite_differences import *
import matplotlib.pyplot as plt
from interp import *
import math
import csv
from read_EFIT import *
from read_EFIT_file import *
from read_iterdb_file import *

#**************Block for user******************************************
#**************Setting up*********************************************
iterdb_file_name = 'DIIID175823.iterdb' #name of the iterdb file
geomfile = 'g175823.04108_257x257'                     #name of the magnetic geometry file

n0=1  #Toroidal mode number
use_ky=True #True if user wants to define ky manually, which will override the n0
ky_define=0.18 #Defined ky for use_ky=True
center_define=0.98 #Define the center in case use_ky=True
mref = 2.        # mass of ion in proton mass

#**************End of Block for user******************************************

#*************Loading the data******************************************
rhot0, te0, ti0, ne0, ni0, nz0, vrot0 = read_iterdb_file(iterdb_file_name)
EFITdict = read_EFIT(geomfile)
xgrid = EFITdict['psipn']
q = EFITdict['qpsi']

uni_rhot = np.linspace(min(rhot0),max(rhot0),len(rhot0)*10.)

te_u = interp(rhot0,te0,uni_rhot)
ne_u = interp(rhot0,ne0,uni_rhot)
vrot_u = interp(rhot0,vrot0,uni_rhot)
q      = interp(xgrid,q,uni_rhot)
tprime_e = -fd_d1_o4(te_u,uni_rhot)/te_u
nprime_e = -fd_d1_o4(ne_u,uni_rhot)/ne_u

center_index = np.argmax((q * np.sqrt(te_u)*(tprime_e+nprime_e))[0:int(len(tprime_e)*0.99)])
x0_center=uni_rhot[center_index]

print('mid pedestal is at r/a = '+str(x0_center))

Lref, Bref, R_major, q0, shat0=get_geom_pars(geomfile,x0_center)

index_begin=np.argmin(abs(uni_rhot-x0_center+1.5*(1-x0_center)))

te_u = te_u[index_begin:len(uni_rhot)-1]
ne_u = ne_u[index_begin:len(uni_rhot)-1]
vrot_u = vrot_u[index_begin:len(uni_rhot)-1]
q      = q[index_begin:len(uni_rhot)-1]
tprime_e = tprime_e[index_begin:len(uni_rhot)-1]
nprime_e = nprime_e[index_begin:len(uni_rhot)-1]
uni_rhot = uni_rhot[index_begin:len(uni_rhot)-1]

center_index = np.argmax(abs(tprime_e*nprime_e))

#*************End of loading the data******************************************

#****************Start setting up ******************

q0      = q[center_index]
ne = ne_u[center_index]
te = te_u[center_index] #it is in eV
#Bref=float(Bref_Gauss)/10000
m_SI = mref *1.6726*10**(-27)
c  = 1
qref = 1.6*10**(-19)
nref = ne
Tref = te * qref
cref = np.sqrt(Tref / m_SI)
Omegaref = qref * Bref / m_SI / c
rhoref = cref / Omegaref 
#******************End setting up ****************

m0 = n0*q0
ky=n0*q0*rhoref/(Lref*x0_center)
kymin = ky
n0_global = n0
te_mid = te_u[center_index]
kyGENE =kymin * (q/q0) * np.sqrt(te_u/te_mid) * (x0_center/uni_rhot) #Add the effect of the q varying

if use_ky==True:
    define_center_index=np.argmin(abs(uni_rhot-center_define))
    kyGENE_center=kyGENE[define_center_index]
    factor=ky_define/kyGENE_center
    n_define=n0*factor
    n0_global=n_define
    print('The toroidal mode number will be '+str(n_define))
    kyGENE = kyGENE*factor


#***Calculate omeage star********************************
#from mtm_doppler
omMTM = kyGENE*(tprime_e+nprime_e)
gyroFreq = 9.79E3/np.sqrt(mref)*np.sqrt(te_u)/Lref
mtmFreq = omMTM*gyroFreq/2./np.pi/1000.
omegaDoppler = vrot_u*n0_global/2./np.pi/1E3
omega=mtmFreq + omegaDoppler


if 1 == 1:
    plt.plot(uni_rhot,omegaDoppler,label='Doppler Shift')
    plt.plot(uni_rhot,mtmFreq,label='Electron Diamagnetic (MTM in plasma frame)')
    plt.plot(uni_rhot,omega,label='Diamagnetic plus Doppler (MTM in lab frame)')
    #plt.axis([0.92,1.,-50.,300.])
    plt.xlabel('rhot')
    plt.ylabel('frequency (kHz)')
    plt.legend(loc = 2, prop = {'size':12})
    plt.show()
    #file_name = 'DIIID_Diallo_freq'
    #plt.savefig(file_name+'.pdf', format='pdf')

x=input("The location of interest:")
#print(de[:,0])
id1=np.argmin(abs(uni_rhot-float(x)))

print('The Doppler shift:')
print(str(omegaDoppler[id1])+"kHz")


