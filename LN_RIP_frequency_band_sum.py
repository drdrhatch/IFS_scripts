#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Created by Max T. Curie 09/07/2020
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import imageio   #for making animation
import optparse as op
import sys
from fieldlib import fieldfile
from geomWrapper import ky
from geomWrapper import init_read_geometry_file
from read_write_geometry import read_geometry_local
from ParIO import Parameters 
from LN_tools import start_end_time
from LN_tools import LN_apar_frequency_nz

#!!!!!!!!!!!!!Needed to add the doppler shift effect!!!!!!!!!
#!!!!!!!!!!!!!Search 'Debatable' for the line that the author is not sure of.!!!!!!!!!!!

#*****************************************************************
#*******************Beginning of the User block*******************

Delta_Z=0.07  #7cm as bin for Z 
scan_all_Z=True #Change to True if one want to scan across the whole height
max_Z0=0.21    #Add a small number so it is even
min_Z0=-0.21

Detailed_report=True    #Swtich to true if one wants to add a detailed report

frequency_all=False      #Switch to True if one wants to sum over all fequency 
#!!!!!!!!!!If simulation has EXB shear on,  use lab frame.
#!!!!!!!!!!If simulation has EXB shear off, use plasma frame. 
#For the calculating the doppler shift, please use max_DopplerShift.py in the prob folder
frequency_max0=250.       #maximum frequency(Lab Frame) to sum over in kHz
frequency_min0=200.       #minimum frequency(Lab Frame) in sum over in kHz

EXB_on=False             #Switch to True is the simulation has the EXB shear
Doppler_shift = 34.55932120349099 #Doppler_shift in kHz

pic_path='pic'        #path one want to store the picture and video in
csv_path='csv'        #path one want to store the picture and video in

#*********************End of the User block***********************
#*****************************************************************
max_Z=max_Z0*1.00001    #Add a small number so it is even
min_Z=min_Z0

if EXB_on==True:
    Doppler_shift=0.

frequency_max=frequency_max0-Doppler_shift  #maximum frequency(Plasma Frame) to sum over in kHz
frequency_min=frequency_min0-Doppler_shift  #minimum frequency(Plasma Frame) in sum over in kHz

parser=op.OptionParser(description='Some infrastructure for reading in, manipulating, and plotting nonlinear field data.')
#parser.add_option('--plot_theta','-g',action='store_const',const=False,help = 'Plot global mode structures decomposed in poloidal m number.',default=True)
options,args=parser.parse_args()
print("options",options)
print("args",args)
if len(args)!=1:
    exit("""
Please include run number as argument (e.g., 0001)."
    \n""")
suffix = args[0]

if suffix in ['dat','.dat']:
     suffix = '.dat'
else:
     suffix = '_'+suffix

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
    min_Z=min(real_Z)
    max_Z=max(real_Z)
Z_grid=np.arange(min_Z,max_Z,Delta_Z)
Z_list = Z_grid[:-1]+Delta_Z/2
print("Z_list: "+str(Z_list))
Z_list_cm=Z_list*100.

#*********Get geometry from magn_geometry***********

time_start,time_end=start_end_time(suffix,pars)

RIP_list=np.zeros(len(Z_list))
RIP_err=np.zeros(len(Z_list))
for nZ_list in range(len(Z_list)):
    RIP_list_temp=[]
    for nZ in range(len(real_Z)):
        Z=real_Z[nZ]
        if Z_list[nZ_list]-Delta_Z/2.<Z and Z<=Z_list[nZ_list]+Delta_Z/2.:
            print('***************************')
            print('***************************')
            print('***************************')
            print('***************************')
            print('Looking at Z='+str(Z)+'m')
            #one can restrict the ky range corresprons to frequency 
            inz=nZ
            if Detailed_report==True:
                frequency_kHZ,amplitude_frequency_sum,amplitude_growth_sum=LN_apar_frequency_nz(suffix,inz,time_start,time_end,plot=True,pic_path='pic/nz_'+str(inz),csv_path='csv/nz_'+str(inz),output_csv=True,show=False)
            else:
                frequency_kHZ,amplitude_frequency_sum,amplitude_growth_sum=LN_apar_frequency_nz(suffix,inz,time_start,time_end,plot=False,pic_path='pic',csv_path='csv',output_csv=False,show=False)
            
            
            #!!!!!!!!!!Debatable
            #amplitude_frequency_sum=list(amplitude_frequency_sum)
            amplitude_frequency_sum=amplitude_frequency_sum[::-1]+amplitude_frequency_sum
            #amplitude_frequency_sum=Reverse+Itself
            #!!!!!!!!!!Debatable

            if frequency_all==True:
                B1=sum(amplitude_frequency_sum)/(2.*float(len(amplitude_frequency_sum)))  #double count
            else:
                beginning_index0=np.argmin(abs(frequency_kHZ-frequency_min))
                ending_index0=np.argmin(abs(frequency_kHZ-frequency_max))
                beginning_index=min(beginning_index0,ending_index0)
                ending_index=max(beginning_index0,ending_index0)
                amplitude_frequency_sum_band=amplitude_frequency_sum[beginning_index:ending_index+1]
                B1=sum(amplitude_frequency_sum_band)/(float(len(amplitude_frequency_sum_band)))
            print("B1="+str(B1)+"Gauss")
            

            #plt.clf()
            #plt.plot(frequency_kHZ,amplitude_frequency_sum)
            #plt.plot(frequency_kHZ)
            #plt.plot(frequency_kHZ[beginning_index:ending_index+1],amplitude_frequency_sum_band)
            #plt.axvline(frequency_kHZ[beginning_index],color='blue',label="freq min",alpha=1)
            #plt.axvline(frequency_kHZ[ending_index],color='red',label="freq max",alpha=1)
            #plt.show()
            
            
            RIP_list_temp.append(B1)
    RIP_list[nZ_list]=np.average(RIP_list_temp)
    RIP_err[nZ_list]=np.std(RIP_list_temp)
        
#*********************Output**************************

Z_err_cm=[Delta_Z*100./2.]*len(Z_list)
RIP_err=[0.]*len(Z_list)

d = {'Z(cm)':Z_list_cm,'Z_err(cm)':Z_err_cm,'B_R(Gauss)':RIP_list,'B_R_err(Gauss)':RIP_err}
df=pd.DataFrame(d, columns=['Z(cm)','Z_err(cm)','B_R(Gauss)','B_R_err(Gauss)'])
if frequency_all==True:
    df.to_csv('csv/RIP_from_all_freq.csv',index=False)
else:
    df.to_csv('csv/RIP_from '+str(round(frequency_min0, 2))+'kHz to '+str(round(frequency_max0, 2))+'kHz.csv',index=False)

plt.clf()
ax=df.plot(kind='scatter',x='Z(cm)',xerr='Z_err(cm)',y='B_R(Gauss)',yerr='B_R_err(Gauss)',grid=True,color='blue')
ax.set_xlabel(r'$Height(cm)$',fontsize=15)
ax.set_ylabel(r'$\bar{B}_r(Gauss)$',fontsize=15)
if frequency_all==True:
    plt.title(r'$\bar{B}_r$(Gauss) from _all_freq',fontsize=18)
else:
    plt.title(r'$\bar{B}_r$(Gauss) from '+str(round(frequency_min0, 2))+'kHz to '+str(round(frequency_max0, 2))+'kHz',fontsize=18)
#plt.xlim(min_Z0*100.,max_Z0*100.)
plt.savefig('pic/0summary_RIP_freq_band.png')
plt.show()
