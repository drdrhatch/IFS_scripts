import sys
import numpy as np
import matplotlib.pyplot as plt
from D_chi_ratio import D_over_chi
from nrgWrapper import read_from_nrg_files
from parIOWrapper import read_species_gradients
from parIOWrapper import read_species_tempdens


def start_end_time(suffix,pars):  #suffix in the format of "_1" or ".dat"


    if pars['n_spec'] == 1:
        time, nrge = read_from_nrg_files(pars,suffix,False)
    elif pars['n_spec'] == 2:
        time, nrgi, nrge = read_from_nrg_files(pars,suffix,False)
    elif pars['n_spec'] == 3:
        time, nrgi, nrge, nrgz = read_from_nrg_files(pars,suffix,False)

    #print(str(np.shape(time)))
    #print(str(np.shape(nrge)))

    plt.clf()
    plt.plot(time,nrge[:,6],label="Q_es of electron")
    plt.plot(time,nrge[:,7],label="Q_em of electron")
    plt.title('nrg of electron')
    plt.xlabel('time')
    plt.legend()
    plt.show()
    
    scan_all = str(raw_input("Scan all(Y/N):\n"))
    if scan_all=='n' or scan_all=='N':
        time_start = float(raw_input("Start time:\n"))
        time_end = float(raw_input("End time:\n"))
        time_start_index=np.argmin(abs(time-float(time_start)))
        time_end_index=np.argmin(abs(time-float(time_end)))
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

    plt.clf()
    plt.plot(time,nrge[:,6],label="Q_es of electron")
    plt.plot(time,nrge[:,7],label="Q_em of electron")
    plt.title('nrg of electron')
    plt.axvline(time_start,color='red',label="time start",alpha=1)
    plt.axvline(time_end,color='blue',label="time end",alpha=1)
    plt.xlabel('time')
    plt.legend()
    plt.show()



    return time_start,time_end


    
