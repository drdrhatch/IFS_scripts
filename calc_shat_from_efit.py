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
parser.add_option('--show_plots','-p',action='store_true',help = 'Show plots.',default=False)
options,args=parser.parse_args()
if len(args)!=1:
    exit("""
Please include efit file name as argument.
    \n""")

efit_file_name = args[0]
show_plots = options.show_plots

#efit_file_name = 'g030701.01187'

def calc_shat(q,rhot):
    rhot0 = np.arange(10000.0)/9999.0
    q0 = interp(rhot,q,rhot0)

    qprime = fd_d1_o4_smoothend(q0,rhot0)
    shat = rhot0/q0*qprime
    return rhot0,q0,shat 

def calc_shat_wpsi(qin,psiin,rhot,rhop):
    #print len(qin),len(psiin)
    #print len(rhop[ind0:]),len(rhot[ind0:])
    #print len(rhot0)
    #plt.plot(psiin,qin)
    #plt.show() 
    #q0 = full_interp(qin,psiin,rhop[ind0:]**2,rhot[ind0:],rhot0)
    rhot0 = np.linspace(0,1,1000)
    psi0 = interp(rhot,rhop**2,rhot0)
    iend = np.argmin(abs(psi0-psiin[-2]))
    #psi0 = psi0[:iend]
    #rhot0 = rhot0[:iend]
    q0 = interp(psiin[:-2],qin[:-2],psi0[:iend])
    nmore = 1000-len(q0)
    for i in range(nmore):
        q0 = np.append(q0,q0[-1])
    #print('len(q0)',len(q0))
    #stop
    if show_plots:
        plt.plot(psiin,qin,label='original')
        plt.plot(psi0,q0,label='interoplated')
        plt.xlabel('psi')
        plt.ylabel('q')
        plt.legend()
        plt.show()

    qprime = fd_d1_o4_smoothend(q0,rhot0)
    shat = rhot0/q0*qprime
    if show_plots:
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
if show_plots:
    plt.plot(psi,q)
    plt.xlabel('psi')
    plt.ylabel('q')
    plt.show()

rtrp = np.genfromtxt('rt_rp_'+efit_file_name)
rt = rtrp[:,0]
rp = rtrp[:,1]

#ind8 = np.argmin(abs(rt-0.8))
rhot0,q0,shat = calc_shat_wpsi(q,psi,rt,rp)

#shat_out = shat[xind]
#print "Assuming shat = ",shat_out
#if j==0:
print( "Saving shat.dat")
f = open('shat.dat','w')
f.write('#1.rhot0 2.q0 3.shat\n')
np.savetxt(f,np.column_stack((rhot0,q0,shat)))
f.close()

