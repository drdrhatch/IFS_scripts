#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Created by Max T. Curie 09/07/2020
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ParIO import Parameters 
from LN_tools import get_suffix
from nrgWrapper import read_from_nrg_files


suffix_list=[1,2,3,6,7]
#average_transport=1.67
average_transport_EM=1.67
average_transport_ES=0.48
average_transport=average_transport_ES+average_transport_EM
mark_size=5

for i in range(len(suffix_list)):
    suffix=str(suffix_list[i])
    if suffix in ['dat','.dat']:
        suffix = '.dat'
    else:
        suffix = '_'+suffix
    suffix_list[i]=suffix

#Import the parameters from parameter file using ParIO
par = Parameters()
par.Read_Pars('parameters'+str(suffix_list[0]))
pars = par.pardict

def get_nrg(suffix,pars):  #suffix in the format of "_1" or ".dat"
    if pars['n_spec'] == 1:
        time, nrge = read_from_nrg_files(pars,suffix,False)
    elif pars['n_spec'] == 2:
        time, nrgi, nrge = read_from_nrg_files(pars,suffix,False)
    elif pars['n_spec'] == 3:
        time, nrgi, nrge, nrgz = read_from_nrg_files(pars,suffix,False)
    return time,nrge

def plot_trace(nrge_es,nrge_em):
    plt.clf()
    plt.plot(time,nrge_es,'.',ms=mark_size,label=r"$Q_{es}$ of electron")
    plt.plot(time,nrge_em,'.',ms=mark_size,label=r"$Q_{em}$ of electron")
    plt.title('nrg of electron')
    plt.xlabel('time')
    plt.legend()
    plt.show()
    
    scan_all = str(input("Scan all(Y/N):\n"))
    if scan_all=='n' or scan_all=='N':
        time_start = float(input("Start time:\n"))
        time_end = float(input("End time:\n"))
        time_start_index=np.argmin(abs(np.array(time)-float(time_start)))
        time_end_index=np.argmin(abs(np.array(time)-float(time_end)))
        time_start = time[time_start_index]
        time_end = time[time_end_index]
    elif scan_all=='y' or scan_all=='Y':
        time_start = time[0]
        time_end = time[-1]
        time_start_index=0
        time_end_index=len(time)-1
    else:
        print("Please respond with y, Y , n, N")
        sys.exit()
    
    '''
    plt.clf()
    plt.plot(time,nrge_es,'.',ms=mark_size,label=r"$Q_{es}$ of electron")
    plt.plot(time,nrge_em,'.',ms=mark_size,label=r"$Q_{em}$ of electron")
    plt.title('nrg of electron')
    plt.axvline(time_start,color='red',label="time start",alpha=1)
    plt.axvline(time_end,color='blue',label="time end",alpha=1)
    plt.xlabel('time(a/cs)')
    #plt.ylabel('Transport(MW)')
    plt.legend()
    plt.show()
    '''

    Q_em=np.mean(nrge_em[time_start_index:time_end_index-1])


    plt.clf()
    plt.plot(time,nrge_es/Q_em*average_transport_EM,'.',ms=mark_size,label=r"$Q_{es}$ of electron")
    plt.plot(time,nrge_em/Q_em*average_transport_EM,'.',ms=mark_size,label=r"$Q_{em}$ of electron")
    plt.title('Heat Transport of electron')
    plt.axhline(average_transport_EM,color='green',label=r"Average Heat Transport $Q_{em}$",alpha=1)
    plt.axhline(average_transport_ES,color='orange',label=r"Average Heat Transport $Q_{es}$",alpha=1)
    plt.axvline(time_start,color='red',label="time start",alpha=1)
    plt.axvline(time_end,color='blue',label="time end",alpha=1)
    plt.xlabel('time(a/cs)')
    plt.ylabel('Transport(MW)')
    plt.legend()
    plt.show()

    print('np.mean(nrge_es),np.std(nrge_es)'\
    	+str(np.mean(nrge_es[time_start_index:time_end_index-1]))+', '\
    	+str(np.std(nrge_es[time_start_index:time_end_index-1])) )

    print('np.mean(nrge_em),np.std(nrge_em)'\
    	+str(np.mean(nrge_em[time_start_index:time_end_index-1]))+', '\
    	+str(np.std(nrge_em[time_start_index:time_end_index-1])) )

    print('Q_es='\
    	+str(np.mean(nrge_es[time_start_index:time_end_index-1])/Q_em*average_transport_EM)+'+-'\
    	+str(np.std(nrge_es[time_start_index:time_end_index-1])/Q_em*average_transport_EM) )

    print('Q_em='\
    	+str(np.mean(nrge_em[time_start_index:time_end_index-1])/Q_em*average_transport_EM)+'+-'\
    	+str(np.std(nrge_em[time_start_index:time_end_index-1])/Q_em*average_transport_EM) )

    return time_start,time_end


nrge_sep_es=[]
nrge_sep_em=[]
time_list=[]
for suffix in suffix_list:
    time,nrge=get_nrg(suffix,pars)
    nrge_sep_es.append(np.ndarray.tolist(nrge[:,6]))
    nrge_sep_em.append(np.ndarray.tolist(nrge[:,7]))
    time_list.append(np.ndarray.tolist(time))

nrge_es=[]
nrge_em=[]
time=[]

for i in time_list:   
    time=time+i

for i in nrge_sep_es:
    nrge_es=nrge_es+i

for i in nrge_sep_em:
    nrge_em=nrge_em+i


plot_trace(nrge_es,nrge_em)