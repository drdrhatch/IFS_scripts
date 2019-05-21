#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import matplotlib.pyplot as plt
import optparse as op
import numpy as np
import matplotlib.pyplot as plt
from interp import *
from finite_differences import *

parser=op.OptionParser(description='Calculates shat from a file (first run my_efit_tools.py -p -n gfile_name to get it from an EFIT file).  Second argument is the file name of a gene profiles file (for rhotor and rhopol).  Outputs shat to file.  Arguments Binfo_filename, gene_profiles_filename, rho_tor.  ')
options,args=parser.parse_args()
if len(args)!=3:
    exit("""
Please include Binfo file name and rho_tor.
    \n""")

binfo = args[0]
gene = args[1]
rhot0 = float(args[2])

#efit_file_name = 'g030701.01187'

def calc_shat(q,rhot):
    rhot0 = np.arange(10000.0)/9999.0
    q0 = interp(rhot,q,rhot0)

    qprime = fd_d1_o4(q0,rhot0)
    shat = rhot0/q0*qprime
    return rhot0,q0,shat 

def calc_shat_wpsi(qin,psiin,rhot,rhop,rhot_range=[0.8,1.0]):
    drhot = rhot_range[1]-rhot_range[0]
    rhot0 = np.arange(10000.0)/9999.0*drhot+rhot_range[0]
    ind0 = np.argmin(abs(rhot[:]-rhot_range[0]))
    #print len(qin),len(psiin)
    #print len(rhop[ind0:]),len(rhot[ind0:])
    #print len(rhot0)
    #plt.plot(psiin,qin)
    #plt.show() 
    q0 = full_interp(qin,psiin,rhop[ind0:]**2,rhot[ind0:],rhot0,verify_interp = True)

    qprime = fd_d1_o4(q0,rhot0)
    shat = rhot0/q0*qprime
    plt.plot(rhot0,q0)
    plt.title('q')
    plt.show()
    plt.plot(rhot0,shat)
    plt.title('shat')
    plt.show()
    return rhot0,q0,shat 

    #Calculating shat
qfilein = np.genfromtxt(binfo)
q = qfilein[:,4]
psi = qfilein[:,1]
genefile = np.genfromtxt(gene)
ind8 = np.argmin(abs(genefile[:,0]-0.8))
rhot0,q0,shat = calc_shat_wpsi(q,psi,genefile[ind8:,0],genefile[ind8:,1])

#shat_out = shat[xind]
#print "Assuming shat = ",shat_out
#if j==0:
print "Saving shat.dat"
f = open('shat.dat','w')
f.write('#1.rhot0 2.q0 3.shat\n')
np.savetxt(f,np.column_stack((rhot0,q0,shat)))
f.close()





