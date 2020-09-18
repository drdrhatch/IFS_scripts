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

max_Z0=0.03    #
min_Z0=-0.03

Detailed_report=False    #Swtich to true if one wants to add a detailed report


EXB_on=False             #Switch to True is the simulation has the EXB shear
Doppler_shift = 34.55932120349099 #Doppler_shift in kHz

pic_path='pic'        #path one want to store the picture and video in
csv_path='csv'        #path one want to store the picture and video in

#*********************End of the User block***********************
#*****************************************************************

if EXB_on==True:
    Doppler_shift=0.

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

#*********Get geometry from magn_geometry***********

time_start,time_end=start_end_time(suffix,pars)


amplitude_frequency_sum_Z=0

nZ_count=0.

for nZ in range(len(real_Z)):
    Z=real_Z[nZ]
    if min_Z0<Z and Z<=max_Z0:
        print('***************************')
        print('***************************')
        print('***************************')
        print('***************************')
        print('Looking at Z='+str(Z)+'m')
        #one can restrict the ky range corresprons to frequency 
        nZ_count=nZ_count+1.
        inz=nZ

        frequency_kHZ,amplitude_frequency_sum,amplitude_growth_sum=LN_apar_frequency_nz(suffix,inz,time_start,time_end,plot=False,pic_path='pic',csv_path='csv',output_csv=False,show=False)
        amplitude_frequency_sum=amplitude_frequency_sum[::-1]+amplitude_frequency_sum
        amplitude_frequency_sum_Z = amplitude_frequency_sum_Z + amplitude_frequency_sum

frequency_kHZ=frequency_kHZ+Doppler_shift
amplitude_frequency_sum_Z=amplitude_frequency_sum_Z/nZ_count


d = {'f(kHz)':frequency_kHZ,'B_R(Gauss)':amplitude_frequency_sum_Z}
df=pd.DataFrame(d, columns=['f(kHz)','B_R(Gauss)'])
df.to_csv('0RIP_freq_spectrogram.csv',index=False)


plt.clf()
ax=df.plot(kind='scatter',x='f(kHz)',y='B_R(Gauss)',grid=True,color='blue')
ax.set_xlabel(r'$Frequency(kHz)$',fontsize=15)
ax.set_ylabel(r'$\bar{B}_r(Gauss)$',fontsize=15)
plt.title(r'$\bar{B}_r$(Gauss) spectrogram',fontsize=18)
plt.savefig('0RIP_freq_spectrogram.png')
plt.show()
