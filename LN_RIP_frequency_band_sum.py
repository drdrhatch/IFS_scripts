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

#*****************************************************************
#*******************Beginning of the User block*******************

Delta_Z=0.07  #7cm as bin for Z 
scan_all_Z=False #Change to True if one want to scan across the whole height
max_Z0=0.21   
min_Z0=-0.21
time_step=100     #read time with this step size

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

RIP_list=np.zeros(len(Z_list))
RIP_err=np.zeros(len(Z_list))

for i_Z_list in range(len(Z_list)):

    max_Z0=Z_list[i_Z_list]+Delta_Z/2.    #in the unit of meter
    min_Z0=Z_list[i_Z_list]-Delta_Z/2.   #in the unit of meter

    uni_freq,amplitude_frequency_uni_ky_sum,amplitude_frequency_uni=\
        B1_ky_f_spectrum_Z_sum_time_set(suffix,iterdb_file_name,\
            min_Z0,max_Z0,Outboard_mid_plane,\
            time_step,time_start,time_end,\
            plot,show,csv_output,pic_path,csv_path)
    #uni_freq is 1D array
    #amplitude_frequency_uni_ky_sum is 1D array B1(f) length weighted average over Z, sum of ky
    #amplitude_frequency_uni is 2D array B1(ky,f) length weighted average over Z

    if frequency_all==True:
        B1=np.average(amplitude_frequency_uni_ky_sum)/2.  #double count
    else:
        beginning_index0=np.argmin(abs(uni_freq-frequency_min))
        ending_index0=np.argmin(abs(uni_freq-frequency_max))
        beginning_index=min(beginning_index0,ending_index0)
        ending_index=max(beginning_index0,ending_index0)
        amplitude_frequency_sum_band=amplitude_frequency_uni_ky_sum[beginning_index:ending_index+1]
        B1=np.average(amplitude_frequency_sum_band)
        print("B1="+str(B1)+"Gauss")

    RIP_list[i_Z_list]=B1


Z_err_cm=[Delta_Z*100./2.]*len(Z_list)
RIP_err=[0.]*len(Z_list)

d = {'Z(cm)':Z_list_cm,'Z_err(cm)':Z_err_cm,'B_R(Gauss)':RIP_list,'B_R_err(Gauss)':RIP_err}
df=pd.DataFrame(d, columns=['Z(cm)','Z_err(cm)','B_R(Gauss)','B_R_err(Gauss)'])
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
plt.savefig('pic/0summary_RIP_freq_band.png')
plt.show()


