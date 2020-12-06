#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Created by Max T. Curie 09/07/2020
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ParIO import Parameters 
from fieldlib import fieldfile
from geomWrapper import init_read_geometry_file
from read_write_geometry import read_geometry_local
from read_iterdb_file import read_iterdb_file

from LN_tools import get_suffix
from LN_tools import start_end_time
from LN_tools import B1_ky_f_spectrum_Z_sum_time_set
from LN_tools import n1_ky_f_spectrum_Z_sum_time_set

#*****************************************************************
#*******************Beginning of the User block*******************

Delta_Z=0.07  #7cm as bin for Z 
scan_all_Z=False #Change to True if one want to scan across the whole height
max_Z0=0.21   
min_Z0=-0.21
time_step=1     #read time with this step size

frequency_all=False      #Switch to True if one wants to sum over all fequency 

frequency_max=500.       #maximum frequency(Lab Frame) to sum over in kHz
frequency_min=150.       #minimum frequency(Lab Frame) in sum over in kHz

pic_path='pic'        #path one want to store the picture and video in
csv_path='csv'        #path one want to store the picture and video in
iterdb_file_name='/global/u1/m/maxcurie/max/Cases/DIIID175823_250k/DIIID175823.iterdb'

#*********************End of the User block***********************
#*****************************************************************

Outboard_mid_plane=False  #change to True if one wants to only want to look at outboard mid-plane
plot=False
show=False
csv_output=False

suffix=get_suffix()

#Import the parameters from parameter file using ParIO
par = Parameters()
par.Read_Pars('parameters'+suffix)
pars = par.pardict
#*********Get geometry from magn_geometry***********
#getting B field using read_write_geometry.py
gpars,geometry = read_geometry_local(pars['magn_geometry'][1:-1]+suffix)
#Get geom_coeff from ParIO Wrapper
geom_type, geom_pars, geom_coeff = init_read_geometry_file(suffix, pars)

real_R=geometry['gl_R'] #it is major radius in meter
real_Z=geometry['gl_z'] #it is height in meter, midland is 0, down is negative ,up is positive

n0_GENE = pars['nref']
B0_GENE =geometry['gBfield']




B_gauss=10.**4               #1T=10^4Gauss
qref = 1.6E-19              #in C
c  = 1.                     #in 3*10^8m/s
m_kg = 1.673E-27            #in kg
Bref = pars['Bref']         #in Tesla
Tref = pars['Tref']         #in keV
nref = pars['nref']         #in 10^(19) /m^3
Lref = pars['Lref']         #in m
mref = pars['mref']         #in proton mass(kg)
nref = nref * 1.E19         #in the unit of /m^3
Tref = Tref * qref * 1.E03  #in the unit of J
mref = mref * m_kg          #in the unit of kg
pref = nref * Tref          #in Pa*kB
cref = np.sqrt(Tref / mref) #in the unit of m/s
Omegaref = qref * Bref / mref / c  #in rad/s
rhoref = cref / Omegaref           #in m rho_i(ion gyroradii)
rhorefStar = rhoref / Lref         #Unitless




if scan_all_Z==True:
    min_Z0=min(real_Z)
    max_Z0=max(real_Z)

max_Z=max_Z0*1.00001    #Add a small number so it is even
min_Z=min_Z0

Z_grid=np.arange(min_Z,max_Z,Delta_Z)
Z_list = Z_grid[:-1]+Delta_Z/2.
print("Z_list: "+str(Z_list))
Z_list_cm=Z_list*100.

time_start,time_end=start_end_time(suffix,pars)

B1_list=np.zeros(len(Z_list))
B1_err=np.zeros(len(Z_list))
B0_list=np.zeros(len(Z_list))
B0_err=np.zeros(len(Z_list))
n1_list=np.zeros(len(Z_list))
n1_err=np.zeros(len(Z_list))
n0_list=np.zeros(len(Z_list))
n0_err=np.zeros(len(Z_list))
Ratio_list=np.zeros(len(Z_list))
Ratio_err=np.zeros(len(Z_list))

n0=nref #nref
B0=abs(np.mean(B0_GENE)*Bref*B_gauss)

for i_Z_list in range(len(Z_list)):

    max_Z0=Z_list[i_Z_list]+Delta_Z/2.    #in the unit of meter
    min_Z0=Z_list[i_Z_list]-Delta_Z/2.   #in the unit of meter

    n1_uni_freq,n1_amplitude_frequency_uni_ky_sum,n1_amplitude_frequency_uni=\
    n1_ky_f_spectrum_Z_sum_time_set(suffix,iterdb_file_name,\
    min_Z0,max_Z0,Outboard_mid_plane,time_step,time_start,time_end,\
    plot,show,csv_output,pic_path,csv_path)


    B1_uni_freq,B1_amplitude_frequency_uni_ky_sum,B1_amplitude_frequency_uni=\
        B1_ky_f_spectrum_Z_sum_time_set(suffix,iterdb_file_name,\
            min_Z0,max_Z0,Outboard_mid_plane,\
            time_step,time_start,time_end,\
            plot,show,csv_output,pic_path,csv_path)

    
    #uni_freq is 1D array
    #amplitude_frequency_uni_ky_sum is 1D array B1(f) length weighted average over Z, sum of ky
    #amplitude_frequency_uni is 2D array B1(ky,f) length weighted average over Z

    if frequency_all==True:
        B1=np.average(B1_amplitude_frequency_uni_ky_sum)/2.  #double count
        n1=np.average(n1_amplitude_frequency_uni_ky_sum)/2.  #double count
    else:
        beginning_index0=np.argmin(abs(B1_uni_freq-frequency_min))
        ending_index0=np.argmin(abs(B1_uni_freq-frequency_max))
        beginning_index=min(beginning_index0,ending_index0)
        ending_index=max(beginning_index0,ending_index0)
        B1_amplitude_frequency_sum_band=abs(B1_amplitude_frequency_uni_ky_sum[beginning_index:ending_index+1])
        B1=np.average(B1_amplitude_frequency_sum_band)
        n1_amplitude_frequency_sum_band=abs(n1_amplitude_frequency_uni_ky_sum[beginning_index:ending_index+1])
        n1=np.average(n1_amplitude_frequency_sum_band)
        print("B1="+str(B1)+"Gauss")
        print("n1="+str(n1)+"/m^3")

    
    B1_list[i_Z_list]=B1
    B0_list[i_Z_list]=B0
    n1_list[i_Z_list]=n1
    n0_list[i_Z_list]=n0

    Ratio_list[i_Z_list]=(B1/B0) / (n1/n0)



Z_err_cm=[Delta_Z*100./2.]*len(Z_list)
Ratio_err=[0.]*len(Z_list)
n1_err=[0.]*len(Z_list)
n0_err=[0.]*len(Z_list)
B1_err=[0.]*len(Z_list)
B0_err=[0.]*len(Z_list)


d = {'Z(cm)':Z_list_cm,'Z_err(cm)':Z_err_cm,\
    'Ratio':Ratio_list,'Ratio_err':Ratio_err,\
    'B_R(Gauss)':B1_list,'B_R_err(Gauss)':B1_err,\
    'B0(Gauss)':B0_list,'B0_err(Gauss)':B0_err,\
    'n1(m^-3)':n1_list,'n1_err(m^-3)':n1_err,\
    'n0(m^-3)':n1_list,'n0_err(m^-3)':n0_err }
df=pd.DataFrame(d, columns=['Z(cm)','Z_err(cm)',\
    'Ratio','Ratio_err',\
    'B_R(Gauss)','B_R_err(Gauss)',\
    'B0(Gauss)','B0_err(Gauss)',\
    'n1(m^-3)','n1_err(m^-3)',\
    'n0(m^-3)','n0_err(m^-3)'])
if frequency_all==True:
    df.to_csv('csv/RIP_from_all_freq.csv',index=False)
else:
    df.to_csv('csv/RIP_from '+str(round(frequency_min, 2))+'kHz to '+str(round(frequency_max, 2))+'kHz.csv',index=False)


plt.clf()
ax=df.plot(kind='scatter',x='Z(cm)',xerr='Z_err(cm)',y='B_R(Gauss)',yerr='B_R_err(Gauss)',grid=True,color='blue')
ax.set_xlabel(r'$Height(cm)$',fontsize=15)
ax.set_ylabel(r'$\bar{B}_r(Gauss)$',fontsize=15)
if frequency_all==True:
    plt.title(r'$\bar{B}_r$(Gauss) from _all_freq',fontsize=18)
else:
    plt.title(r'$\bar{B}_r$(Gauss) from '+str(round(frequency_min, 2))+'kHz to '+str(round(frequency_max, 2))+'kHz',fontsize=18)
#plt.xlim(min_Z0*100.,max_Z0*100.)
plt.savefig('pic/0summary_B1_freq_band.png')




plt.clf()
ax=df.plot(kind='scatter',x='Z(cm)',xerr='Z_err(cm)',y='n1(m^-3)',yerr='n1_err(m^-3)',grid=True,color='blue')
ax.set_xlabel(r'$Height(cm)$',fontsize=15)
ax.set_ylabel(r'$n_1(/m^3)$',fontsize=15)
if frequency_all==True:
    plt.title(r'$n_1(/m^3)$ from _all_freq',fontsize=18)
else:
    plt.title(r'$n_1(/m^3)$ from '+str(round(frequency_min, 2))+'kHz to '+str(round(frequency_max, 2))+'kHz',fontsize=18)
#plt.xlim(min_Z0*100.,max_Z0*100.)
plt.savefig('pic/0summary_n1_freq_band.png')


plt.clf()
ax=df.plot(kind='scatter',x='Z(cm)',xerr='Z_err(cm)',y='Ratio',yerr='Ratio_err',grid=True,color='blue')
ax.set_xlabel(r'$Height(cm)$',fontsize=15)
ax.set_ylabel(r'Ratio$\frac{B_1/B_0}{n_1/n_0}$',fontsize=15)
if frequency_all==True:
    plt.title(r'Ratio from _all_freq',fontsize=18)
else:
    plt.title(r'Ratio from '+str(round(frequency_min, 2))+'kHz to '+str(round(frequency_max, 2))+'kHz',fontsize=18)
#plt.xlim(min_Z0*100.,max_Z0*100.)
plt.savefig('pic/0summary_Ratio_freq_band.png')
plt.show()


