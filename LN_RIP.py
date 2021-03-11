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
from LN_tools import RIP_f_spectrum_FFT_then_sum
from LN_tools import RIP_f_spectrum_sum_then_FFT
from LN_tools import RIP_f_spectrum_sum_kx_then_FFT_then_sum_Z
from LN_tools import RIP_f_spectrum_density
from LN_tools import RIP_k_space_sum

#*****************************************************************
#*******************Beginning of the User block*******************
Outboard_mid_plane=False  #change to True if one wants to only want to look at outboard mid-plane
Delta_Z=0.07  #7cm as bin for Z 
scan_all_Z=False #Change to True if one want to scan across the whole height
max_Z0=0.035   
min_Z0=-0.035

time_step=1     #read time with this step size

frequency_all=True      #Switch to True if one wants to sum over all fequency 

frequency_max=-150.       #maximum frequency(Lab Frame) to sum over in kHz
frequency_min=-500.       #minimum frequency(Lab Frame) in sum over in kHz

run_mode=5	 			#change to 1, if one wants to sum over kx, Z and then do the FFT
						#change to 2, if one wants to do the FFT then sum over kx, Z 
						#change to 3, if one wants to do the sum over kx, FFT then sum over Z 
                        #change to 4, if one wants to have the spectral density (the sum over kx,Z then periodogram)
                        #change to 5, if one wants to sum over kx,ky,Z in k space

pic_path='RIP_pic'        #path one want to store the picture and video in
csv_path='RIP_csv'        #path one want to store the picture and video in
iterdb_file_name='/global/u1/m/maxcurie/max/Cases/DIIID175823_250k/DIIID175823.iterdb'
manual_Doppler=-6.	#if it is number, then manually put in the doppler shift in kHz for n=1, Type False if one to calculate the Doppler shift from ITERDB
#iterdb_file_name='DIIID164880.iterdb'

#*********************End of the User block***********************
#*****************************************************************

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

B0_Gauss=pars['Bref']*10000.

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
    if run_mode==1:
        f,B1_f=\
            RIP_f_spectrum_sum_then_FFT(suffix,iterdb_file_name,manual_Doppler,\
                min_Z0,max_Z0,Outboard_mid_plane,\
                time_step,time_start,time_end,\
                plot,show,csv_output,pic_path,csv_path)
    elif run_mode==2:
        f,B1_f=\
            RIP_f_spectrum_FFT_then_sum(suffix,iterdb_file_name,manual_Doppler,\
                min_Z0,max_Z0,Outboard_mid_plane,\
                time_step,time_start,time_end,\
                plot,show,csv_output,pic_path,csv_path)
    elif run_mode==3:
        f,B1_f=\
            RIP_f_spectrum_sum_kx_then_FFT_then_sum_Z(suffix,iterdb_file_name,manual_Doppler,\
                min_Z0,max_Z0,Outboard_mid_plane,\
                time_step,time_start,time_end,\
                plot,show,csv_output,pic_path,csv_path)
    elif run_mode==4:
        f,B1_f=\
            RIP_f_spectrum_density(suffix,iterdb_file_name,manual_Doppler,\
                min_Z0,max_Z0,Outboard_mid_plane,\
                time_step,time_start,time_end,\
                plot,show,csv_output,pic_path,csv_path)
    #f is 1D array with freqency
    #B1_f is 2D array B1(ky,f) length weighted average over Z

    if run_mode==5:
        B1,B1_error=\
            RIP_k_space_sum(suffix,iterdb_file_name,manual_Doppler,\
                min_Z0,max_Z0,Outboard_mid_plane,\
                time_step,time_start,time_end,\
                plot,show,csv_output,pic_path,csv_path)
    elif run_mode==4:
        B1_RIP0=0.
        B1_RIP_temp=0.
        if frequency_all==True:
            B1_RIP0=np.sum(abs(B1_f)**2.)*abs(f[1]-f[0])
        else:
            for i_f in range(len(f)):
                frequency_max=abs(frequency_max)
                frequency_min=abs(frequency_min)
                if frequency_max<frequency_min:
                    temp=frequency_max
                    frequency_max=frequency_min
                    frequency_min=temp
                if -frequency_max<=f[i_f] and f[i_f]<=-frequency_min:
                    print(B1_f[i_f])
                    B1_RIP0=B1_RIP0+abs(B1_f[i_f])**2.*abs(f[i_f]-f[i_f-1])
        B1=np.sqrt(B1_RIP0)
        B1_error=0.
    else:
        B1_RIP0=0.
        B1_RIP_temp=0.
        if frequency_all==True:
            B1=np.sum(abs(B1_f))
        else:
            for i_f in range(len(f)):
                if frequency_min<=f[i_f] and f[i_f]<=frequency_max:
                    print(B1_f[i_f])
                    B1_RIP0=B1_RIP0+abs(B1_f[i_f])
            B1=B1_RIP0
            
        B1_error=0.

    print("B1="+str(B1)+"Gauss")
    print("B1_erro="+str(B1_error)+"Gauss")
    RIP_list[i_Z_list]=B1
    RIP_err[i_Z_list]=B1_error



Z_err_cm=[Delta_Z*100./2.]*len(Z_list)
if run_mode!=5:
    RIP_err=[0.]*len(Z_list)
B0_list=[B0_Gauss]*len(Z_list)

d = {'Z(cm)':Z_list_cm,'Z_err(cm)':Z_err_cm,'B_R(Gauss)':RIP_list,'B_R_err(Gauss)':RIP_err,'B0(Gauss)':B0_list}
df=pd.DataFrame(d, columns=['Z(cm)','Z_err(cm)','B_R(Gauss)','B_R_err(Gauss)','B0(Gauss)'])
if frequency_all==True:
    df.to_csv(csv_path+'/0RIP_from_all_freq.csv',index=False)
else:
    df.to_csv(csv_path+'/0RIP_from '+str(round(frequency_min, 2))+'kHz to '+str(round(frequency_max, 2))+'kHz.csv',index=False)

plt.clf()
ax=df.plot(kind='scatter',x='Z(cm)',xerr='Z_err(cm)',y='B_R(Gauss)',yerr='B_R_err(Gauss)',grid=True,color='blue')
ax.set_xlabel(r'$Height(cm)$',fontsize=15)
ax.set_ylabel(r'$\bar{B}_r(Gauss)$',fontsize=15)
if frequency_all==True:
    plt.title(r'$\bar{B}_r$(Gauss) from _all_freq',fontsize=18)
else:
    plt.title(r'$\bar{B}_r$(Gauss) from '+str(round(frequency_min, 2))+'kHz to '+str(round(frequency_max, 2))+'kHz',fontsize=18)
#plt.xlim(min_Z0*100.,max_Z0*100.)
plt.savefig(pic_path+'/0summary_RIP_freq_band.png')
plt.show()


