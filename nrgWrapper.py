#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from get_nrg_x import *
from genetools import *
import sys

def species_order(suffix): 
    #*****Start of calcuation of diamagnetic frequency******
    paramfpath="parameters"+str(suffix)
    geneparam=read_parameters(paramfpath)
    #Determine the order of the species
    
    if 'species3' in geneparam: 
        species=['','','']
        for i in range(3):
            if geneparam['species'+str(i+1)]['charge']==-1:
                species[i]='e' #it is eletron
            elif geneparam['species'+str(i+1)]['charge']==1:
                species[i]='i' #it is ion
            elif geneparam['species'+str(i+1)]['charge']>1:
                species[i]='z' #it is impurity
    elif 'species3' not in geneparam: 
        species=['','']
        for i in range(2):
            if geneparam['species'+str(i+1)]['charge']==-1:
                species[i]='e' #it is eletron
            elif geneparam['species'+str(i+1)]['charge']==1:
                species[i]='i' #it is ion
    return species

def read_from_nrg_files(pars,suffix,plot,ncols=10):
    
    if pars['n_spec'] == 1:
        time, nrge = get_nrg0(suffix,nspec=1)
        #return time, nrge
    elif pars['n_spec'] == 2:
        species=species_order(suffix)
        if species[0]=='i' and species[1]=='e':
            time, nrgi, nrge = get_nrg0(suffix,nspec=2,ncols=10)
        elif species[0]=='e' and species[1]=='i':
            time, nrge, nrgi = get_nrg0(suffix,nspec=2,ncols=10)
        #return time, nrgi, nrge
    elif pars['n_spec'] == 3:
        species=species_order(suffix)
        if species[0]=='i' and species[1]=='e' and species[2]=='z':
              time,    nrgi,               nrge,               nrgz = get_nrg0(suffix,nspec=3)
        elif species[0]=='i' and species[1]=='z' and species[2]=='e':
              time,    nrgi,               nrgz,               nrge = get_nrg0(suffix,nspec=3)
        elif species[0]=='e' and species[1]=='i' and species[2]=='z':
              time,    nrge,               nrgi,               nrgz = get_nrg0(suffix,nspec=3)
        elif species[0]=='e' and species[1]=='z' and species[2]=='i':
              time,    nrge,               nrgz,               nrgi = get_nrg0(suffix,nspec=3)
        elif species[0]=='z' and species[1]=='e' and species[2]=='i':
              time,    nrgz,               nrge,               nrgi = get_nrg0(suffix,nspec=3)
        elif species[0]=='z' and species[1]=='i' and species[2]=='e':
              time,    nrgz,               nrgi,               nrge = get_nrg0(suffix,nspec=3)
        else:
            print('Error, species='+str(species))
        
        #return time, nrgi, nrge, nrgz
    else:
        print ("n_spec ="+str( pars['n_spec']))
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
        print ('Reading nrg file are at t = '+str(time[setTime]) )
    else:
        isetTime = np.argmin(abs(time-setTime))
        this_nrg = nrgs[isetTime]
        print ('Reading nrg file are at t = '+str(time[setTime]) )


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
