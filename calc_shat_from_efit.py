#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import matplotlib.pyplot as plt
import optparse as op
import numpy as np
import matplotlib.pyplot as plt
from interp import *
from finite_differences import *
from subprocess import call

parser=op.OptionParser(description='Calculates shat from an efit file.  Argument: efit_file_name.  ')
options,args=parser.parse_args()
if len(args)!=1:
    exit("""
Please include efit file name as argument.
    \n""")

efit_file_name = args[0]

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
    q0 = full_interp(qin,psiin,rhop[ind0:]**2,rhot[ind0:],rhot0)

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

call(['my_efit_tools.py','-p','-n',efit_file_name])
call(['my_efit_tools.py','-c','-n',efit_file_name])

qfilein = np.genfromtxt('Binfo_'+efit_file_name)
q = qfilein[:,4]
psi = qfilein[:,1]

rtrp = np.genfromtxt('rt_rp_'+efit_file_name)
rt = rtrp[:,0]
rp = rtrp[:,1]

ind8 = np.argmin(abs(rt-0.8))
rhot0,q0,shat = calc_shat_wpsi(q,psi,rt[ind8:],rp[ind8:])

#shat_out = shat[xind]
#print "Assuming shat = ",shat_out
#if j==0:
print "Saving shat.dat"
f = open('shat.dat','w')
f.write('#1.rhot0 2.q0 3.shat\n')
np.savetxt(f,np.column_stack((rhot0,q0,shat)))
f.close()

