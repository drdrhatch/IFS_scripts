#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import matplotlib.pyplot as plt
import optparse as op
import numpy as np
import matplotlib.pyplot as plt
from interp import *
from finite_differences import *


parser=op.OptionParser(description='Calculates gradient scale lengths and eta from gene profile files (outputs full info to data files and outputs to screen info at selected rhotor).  Arguments: profiles_e, profiles_i, rho_tor')
parser.add_option('--imp','-i',action='store',type='string',help='Add an impurity profiles file to calculate also an impurity species.',dest='f_gz',default='')
options,args=parser.parse_args()
if len(args)!=3:
    exit("""
Please include profiles_e and profiles_i file names."
    \n""")

f_gz = False

if options.f_gz:
  f_gz = options.f_gz
  print("Impurity file:",f_gz)
f_ge = args[0]
f_gi = args[1]
rhot0 = float(args[2])

gene_e = np.genfromtxt(f_ge)
gene_i = np.genfromtxt(f_gi)
if f_gz:
  gene_z = np.genfromtxt(f_gz)

rhote = gene_e[:,0]
te = gene_e[:,2]
ne = gene_e[:,3]

rhoti = gene_i[:,0]
ti = gene_i[:,2]
ni = gene_i[:,3]

if f_gz:
   rhotz = gene_z[:,0]
   tz = gene_z[:,2]
   nz = gene_z[:,3]

rhot1 = np.arange(1000)/999.0*(rhote[-1]-rhote[0])+rhote[0]
n1 = interp(rhote,ne,rhot1)
t1 = interp(rhote,te,rhot1)
rhot2 = np.arange(1000)/999.0*(rhoti[-1]-rhoti[0])+rhoti[0]
n2 = interp(rhoti,ni,rhot2)
t2 = interp(rhoti,ti,rhot2)

if f_gz:
   rhot3 = np.arange(1000)/999.0*(rhotz[-1]-rhotz[0])+rhotz[0]
   n3 = interp(rhotz,nz,rhot2)
   t3 = interp(rhotz,tz,rhot2)

omt1 = abs(1.0/t1*fd_d1_o4(t1,rhot1))
omn1 = abs(1.0/n1*fd_d1_o4(n1,rhot1))
omt2 = abs(1.0/t2*fd_d1_o4(t2,rhot2))
omn2 = abs(1.0/n2*fd_d1_o4(n2,rhot2))
if f_gz:
   omt3 = abs(1.0/t3*fd_d1_o4(t3,rhot3))
   omn3 = abs(1.0/n3*fd_d1_o4(n3,rhot3))

#print len(rhot1)
#print len(n1)
#print len(omt1)
f = open('profile_info_'+f_ge,'w')
f.write('#1.rhot 2.dummy 3.Te 4.ne 5.omte 6.omne 7.etae \n')
np.savetxt(f,np.column_stack((rhot1,rhot1,t1,n1,omt1,omn1,omt1/omn1)))
f.close()


f = open('profile_info_'+f_gi,'w')
f.write('#1.rhot 2.dummy 3.Ti 4.ni 5.omti 6.omni 7.etai \n')
np.savetxt(f,np.column_stack((rhot2,rhot2,t2,n2,omt2,omn2,omt2/omn2)))
f.close()

if f_gz:
   f = open('profile_info_'+f_gz,'w')
   f.write('#1.rhot 2.dummy 3.Tz 4.nz 5.omtz 6.omnz 7.etaz \n')
   np.savetxt(f,np.column_stack((rhot3,rhot3,t3,n3,omt3,omn3,omt3/omn3)))
   f.close()

irhot1= np.argmin(abs(rhot1[:]  - rhot0))
irhot2= np.argmin(abs(rhot2[:]  - rhot0))
if f_gz:
   irhot3= np.argmin(abs(rhot3[:]  - rhot0))

n10 = n1[irhot1] 
t10 = t1[irhot1] 
n20 = n2[irhot2] 
t20 = t2[irhot2] 
if f_gz:
    n30 = n3[irhot3] 
    t30 = t3[irhot3] 


print('te at rho_tor = '+str(rhot0)+': ',t10)
print('ti at rho_tor = '+str(rhot0)+': ',t20)
print('ne at rho_tor = '+str(rhot0)+': ',n10)
print('ni at rho_tor = '+str(rhot0)+': ',n20)
if f_gz:
   print('tz at rho_tor = '+str(rhot0)+': ',tz[irhot3])
   print('nz at rho_tor = '+str(rhot0)+': ',nz[irhot3])


print('omte at rho_tor = '+str(rhot0)+': ',omt1[irhot1])
print('omti at rho_tor = '+str(rhot0)+': ',omt2[irhot2])
print('omne at rho_tor = '+str(rhot0)+': ',omn1[irhot1])
print('omni at rho_tor = '+str(rhot0)+': ',omn2[irhot2])
if f_gz:
   print('omtz at rho_tor = '+str(rhot0)+': ',omt3[irhot3])
   print('omnz at rho_tor = '+str(rhot0)+': ',omn3[irhot3])




