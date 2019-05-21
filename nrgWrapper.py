#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from get_nrg_x import *
import sys

def read_from_nrg_files(pars,suffix,plot,ncols=10):

    if pars['n_spec'] == 1:
        time, nrge = get_nrg0(suffix,nspec=1)
        #return time, nrge
    elif pars['n_spec'] == 2:
        time, nrgi, nrge = get_nrg0(suffix,nspec=2,ncols=10)
        #return time, nrgi, nrge
    elif pars['n_spec'] == 3:
        time, nrgi, nrge, nrgz = get_nrg0(suffix,nspec=3)
        #return time, nrgi, nrge, nrgz
    else:
        print "n_spec =", pars['n_spec']
        sys.exit("n_spec must be 1, 2 or 3.")

    if plot:
       # plt.semilogy(time,nrgi[:,0],label='nrgi')
        plt.semilogy(time,nrge[:,8],label='nrge')
        plt.legend()
        plt.xlabel('time')
        plt.ylabel(r'$|\delta n|^2$')
        plt.show()

    if pars['n_spec'] == 1:
        return time, nrge
    elif pars['n_spec'] == 2:
        return time, nrgi, nrge
    elif pars['n_spec'] == 3:
        return time, nrgi, nrge, nrgz

def read_Gamma_Q(time,nrgs,print_val,setTime=-1):
     
    if (setTime == -1):
        this_nrg = nrgs[setTime,:]
        print 'Reading nrg file are at t = ', time[setTime]
    else:
        isetTime = np.argmin(abs(time-setTime))
        this_nrg = nrgs[isetTime]
        print 'Reading nrg file are at t = ', time[isetTime]


    Gamma_es = this_nrg[4]
    Gamma_em = this_nrg[5]
    Qheat_es = this_nrg[6]
    Qheat_em = this_nrg[7]
    Pimom_es = this_nrg[8]
    Pimom_em = this_nrg[9]

    if print_val:
       #print "Gamma_es =", Gamma_es
       print ("Gamma_es = %12.4e" % Gamma_es)
       print ("Gamma_em = %12.4e" % Gamma_em)
       print ("Q_es = %12.4e" % Qheat_es)
       print ("Q_em = %12.4e" % Qheat_em)

    return Gamma_es, Gamma_em, Qheat_es, Qheat_em#, Pimom_es, Pimom_em
