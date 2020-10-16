import matplotlib.pyplot as plt
import numpy as np
import csv
import os

from genetools import *

#Last edited by Max Curie 10/16/2020
#Criteria: D_chi, typical frequency, Epar

#**********Start of input setup**********************************************************************


def get_omega(suffix):
    evals = np.genfromtxt('eigenvalues_'+suffix) 
    
    #print(evals)
    gamma_list = []
    omega_list = []
   
    for line in evals:
        omega_list.append(line[1])
        gamma_list.append(line[0])

    indice_list = np.arange(len(omega_list))

    return gamma_list,omega_list,indice_list


def D_chi_e(suffix,index0):

    #from genetools.py
    paramfpath="parameters_"+str(suffix)
    geneparam=read_parameters(paramfpath)

    Tref=geneparam['units']['Tref']
    nref=geneparam['units']['nref']

    species=['e']
    Te=geneparam['species1']['temp']
    ne=geneparam['species1']['dens']
    omn_e=geneparam['species1']['omn']
    omt_e=geneparam['species1']['omt']
    
    #from genetools.py
    nrgfpath="nrg_"+str(suffix)
    nrgdata = read_nrg(nrgfpath)

    D_e=nrgdata[nrgfpath]['e']['PFluxes'][index0]+nrgdata[nrgfpath]['e']['PFluxem'][index0]
    D_e=D_e/ omn_e / ne
    Qes_e=nrgdata[nrgfpath]['e']['HFluxes'][index0]-3./2.*Te*nrgdata[nrgfpath]['e']['PFluxes'][index0]
    Qem_e=nrgdata[nrgfpath]['e']['HFluxem'][index0]-3./2.*Te*nrgdata[nrgfpath]['e']['PFluxem'][index0]
    Q_e = (Qes_e+Qem_e)
    chi_e = (Qes_e+Qem_e) / omt_e / ne / Te

    return Qes_e,Qem_e,Q_e,D_e,chi_e

def D_chi_e_judge(Qes_e,Qem_e,Q_e,D_e,chi_e):

    if abs(Qem_e/Qes_e) >= 1 and abs(D_e/chi_e)<0.6:
        mode="MTM"
    elif abs(Qem_e/Qes_e) >= 1 and abs(D_e/chi_e)>0.6:
        mode="KBM"
    elif abs(Qem_e/Qes_e) < 1 and abs(D_e/chi_e)<0.6:
        mode="ETG"
    else:
        mode="other"
    return mode


def scan_cases():
    csvfile_name='EV_log.csv'
    with open(csvfile_name, 'w', newline='') as csvfile:
        data = csv.writer(csvfile, delimiter=',')
        data.writerow(['Suffix','index','omega','gamma','mode','Qem_e/Qes_e','D_e/chi_e'])
        csvfile.close()

    cwd = os.getcwd()
    filelist = []
    for filename in os.listdir(cwd):
        if filename.startswith("field"):
            filelist.append(filename[-4:])
    filelist.sort()
    for suffix in filelist:
        print('*************reading'+suffix+'*************')
        gamma_list,omega_list,indice_list=get_omega(suffix)
        for index0 in indice_list:
            print('****'+str(index0)+'*****')
            Qes_e,Qem_e,Q_e,D_e,chi_e=D_chi_e(suffix,index0)
            mode=D_chi_e_judge(Qes_e,Qem_e,Q_e,D_e,chi_e)
            print('****'+str(mode)+'*****')
            with open(csvfile_name, 'a+', newline='') as csvfile:
                data = csv.writer(csvfile, delimiter=',')
                data.writerow([suffix,index0,omega_list[index0],gamma_list[index0],mode,Qem_e/Qes_e,D_e/chi_e])
                csvfile.close()


scan_cases()