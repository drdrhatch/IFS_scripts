#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import matplotlib.pyplot as plt
import numpy as np
import optparse as op
import matplotlib.pyplot as plt
from finite_differences import *
from interp import *
from read_EFIT_file import *

parser=op.OptionParser(description='Calculates Er/(R B_pol) from a gene profiles file and an efit file.  Arguments: gene_profiles_file_name_i, efit_file_name.')
parser.add_option('--omtor','-o',type = 'str',action='store',dest="omtor_file",help = 'Include toroidal rotation (specify file name with format: (1)psi (2)omtor(rad/s).',default='empty')
options,args=parser.parse_args()
if len(args)!=2:
    exit("""
Please include gene_profiles_file_name_i and efit_file_name."
    \n""")

pfile = args[0]
efit_file_name = args[1]
omtor_file = options.omtor_file


data = np.genfromtxt(pfile) 

dummy = input("WARNING: profile file must have rho_poloidal as second column! Press any key to continue.\n")
dummy = input("Assuming ion charge is 1 and reference mass is deuterium (press any key to continue).\n")

rhot = data[:,0]
psi = data[:,1]**2
Ti = data[:,2]
ni = data[:,3]

mi = 1.673e-27
ee = 1.602e-19
mref = 2.0*mi
Z = 1.0

Lref, Bref, R_major, q0 = get_dimpar_pars(efit_file_name,0.9)

psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw,psiax,psisep = read_EFIT_file(efit_file_name)

print("R_major",R_major)
print("Bref",Bref)
print("psisep",psisep)
print("psiax",psiax)
psisep0 = psisep-psiax
print("psisep0",psisep0)

#Get everything on the same even psi grid
psi0 = np.linspace(0.0,1.0,num=3000)
#interopolate rho_tor, n, T, qpsi
ni0 = interp(psi,ni,psi0)
Ti0 = interp(psi,Ti,psi0)
Ti0J = Ti0*1000.0*ee
ni00 = ni0*1.0e19
pi0 = Ti0J*ni00
rhot0 = interp(psi,rhot,psi0)
qpsi0 = interp(psip_n,qpsi,psi0)
R0 = interp(psip_n,Rgrid,psi0)

trpeps = rhot0*Lref/R_major
vti = (0.5*Ti0*1000.0*ee/mref)**0.5
#nuii = 2.3031E-5*Lref*(ne)/(Te)**2*(24.0-log(sqrt(ne*1.0E13)/Te*0.001))
coll = 2.3031E-5*Lref*(ni0)/(Ti0)**2*(24.0-log(sqrt(ni0*1.0E13)/Ti0*0.001))
nustar_i=8.0/3.0/pi**0.5*qpsi0/trpeps**1.5*(R_major/Lref)*Z**4*coll

#Comments:
#The paper by Landreman and Ernst is probably the best. Presuming negligible externally driven toroidal rotation (our default assumption in a pedestal in a pedestal), the Er is given by setting eq(6) =0 (and ignoring differences between <B**2> and B**2).
#The quantity which L & E call k|| can be found in eq (15). (which is the negative of the quantity which Sauter calls alpha, but note that there is a typo in the Sauter paper formula for alpha, which was corrected in a later erratum by Sauter. The Ernst formula uses the correct version.)
#The exact definitions of nu* for use in this formula can be found in the original Suater paper. 
#One must have a trapped particle fraction f_t for these formulas. The formal definition of this is a particular integral which is roughly the same size as, but actually isn?t, what the name says. There are various approximations of varying degrees of complexity for this quantity. (The bootstrap code in Prashant?s loop uses the exact messy integral.) In the limit epsilon ->0, f_t = 1.46sqrt(epsilon) (epsilon = outboard minor radius/R ~ 0.3 for standard aspect ratio pedestals). However, for epsilon ~ 0.3, there are modest corrections to this. Without going into the details, from the literature, I derive an approximate formula with finite epsilon corrections, which is better than this at epsilon ~ 0.3:
#f_t = sqrt(epsilon)*(1. + 0.46* (1.-epsilon))
#which agrees with the asymptotic limit epsilon ->0, but goes to 1 at epsilon =1 as it should, and is also significantly closer to the right answer at epsilon = 0.3. 

ft = trpeps**0.5*(1.0+0.46*(1.0-trpeps))
fc = 1.0 - ft

a = 1.0/(1.0+0.5*nustar_i**0.5)
b = -1.17*fc/(1.0-0.22*ft-0.19*ft**2) + 0.25*(1-ft**2)*nustar_i**0.5
c =    0.315*nustar_i**2*ft**6
d = 1.0/(1.0+0.15*nustar_i**2*ft**6)
kpar = -(a*b+c)*d

plt.plot(psi0,kpar)
plt.title('kpar')
plt.show()

dTdpsi = fd_d1_o4(Ti0J,psi0)
dndpsi = fd_d1_o4(ni00,psi0)
dpdpsi = fd_d1_o4(pi0,psi0)
dRdpsi = fd_d1_o4(R0,psi0)

#This is -d Phi_0 / d psi = Er / (R B_theta) from from Landreman and Ernst PPCF 2012 (see comments below)
omegator = (1.0/psisep0/ee)*((1-kpar)*dTdpsi + Ti0J/ni00*dndpsi)
#This is -d Phi_0 / d psi the same thing without neoclassical--simply the pressure gradient.  Should be similar to the above
omegator0 = (1.0/psisep0/ee)*1/ni00*dpdpsi

outfile = 'omegator'+pfile
if omtor_file != 'empty':
   #This includes the effect of toroidal rotation on d Phi_0 / d psi--i.e. the V_{||} term in Eq. 6 of Landreman and Ernst
   print("Including toroidal rotation from file ",omtor_file)
   omtor = np.genfromtxt(omtor_file)
   omtor00 = interp(omtor[:,0],omtor[:,1],psi0)
   omegator += omtor00
   omegator0 += omtor00
   outfile += '_withVtor'

   f=open(outfile,'w')
   f.write('#1.rho_tor 2.psiN 3.Er/(R Bpol) 4.Er0(gradP)/(R Bpol) 5.Omtor(rad/s) \n')
   np.savetxt(f,np.column_stack((rhot0,psi0,omegator,omegator0,omtor00)))
   f.close()
else:
   f=open(outfile,'w')
   f.write('#1.rho_tor 2.psiN 3.Er/(R Bpol) 4.Er0(gradP)/(R Bpol) \n')
   np.savetxt(f,np.column_stack((rhot0,psi0,omegator,omegator0)))
   f.close()

