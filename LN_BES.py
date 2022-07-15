#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Created by Max T. Curie 09/07/2020
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ParIO import Parameters 
from LN_tools import get_suffix
from LN_tools import start_end_time
from LN_tools import BES_f_spectrum_density
from LN_tools import BES_f_spectrum_FFT
from FFT_general import FFT_sum
from FFT_general import spectral_density_sum
#*****************************************************************
#*******************Beginning of the User block*******************

Delta_Z=0.07  #7cm as bin for Z 

time_step=10     #read time with this step size
frequency_all=False      #Switch to True if one wants to sum over all fequency 

frequency_max=-249.5       #maximum frequency(Lab Frame) to sum over in kHz
frequency_min=-250.5       #minimum frequency(Lab Frame) in sum over in kHz

#frequency_max=-150       #maximum frequency(Lab Frame) to sum over in kHz
#frequency_min=-500       #minimum frequency(Lab Frame) in sum over in kHz

run_mode=2	 			#change to 1, if one wants to do the FFT and then sum over kx, Z 
                        #change to 2, if one wants to have the spectral density (the periodogram then sum over kx,Z )
                        

window_for_FFT='hann'     #Default is 'hann', other options
percent_window=1.        #enter 0(smooth,low resolution) to 1(not smooth,high resolution), 
						  #or 'Default'   

pic_path='BES_pic'        #path one want to store the picture and video in
csv_path='BES_csv'        #path one want to store the picture and video in
iterdb_file_name='/global/u1/m/maxcurie/max/Cases/DIIID175823_250k/DIIID175823.iterdb'
manual_Doppler=-7.	#if it is number, then manually put in the doppler shift in kHz for n=1, Type False if one to calculate the Doppler shift from ITERDB

#*********************End of the User block***********************
#*****************************************************************

max_Z0=Delta_Z/2.   
min_Z0=-Delta_Z/2.

Outboard_mid_plane=False  #change to True if one wants to only want to look at outboard mid-plane
plot=False
show=False
csv_output=False

suffix=get_suffix()
par = Parameters()
par.Read_Pars('parameters'+suffix)
pars = par.pardict
time_start,time_end=start_end_time(suffix,pars)



if run_mode==1:
    f,n1_f=\
        BES_f_spectrum_FFT(suffix,iterdb_file_name,manual_Doppler,min_Z0,max_Z0,\
                Outboard_mid_plane,time_step,time_start,time_end,\
                plot,show,csv_output,pic_path,csv_path)
    n1_BES0,n1_error=FFT_sum(f,n1_f,frequency_min,frequency_max,frequency_all)
elif run_mode==2:
    f,n1_f=\
        BES_f_spectrum_density(suffix,iterdb_file_name,manual_Doppler,min_Z0,max_Z0,\
                Outboard_mid_plane,time_step,time_start,time_end,percent_window,window_for_FFT,\
                plot,show,csv_output,pic_path,csv_path)
    n1_BES0,n1_error=spectral_density_sum(f,n1_f,frequency_min,frequency_max,frequency_all)


len_f=len(f)
n1_BES0=0.
n1_BES_temp=0.



if frequency_all==True:
    print(n1_f)
    n1_BES0=np.sum(abs(n1_f))*abs(f[1]-f[0])
    f_sum=abs(np.max(f)-np.min(f))
else:
    for i_f in range(len(f)):
        if frequency_min<=f[i_f] and f[i_f]<=frequency_max:
            print(n1_f[i_f])
            n1_BES0=n1_BES0+abs(n1_f[i_f])**2.*abs(f[i_f]-f[i_f-1])
n1_BES=n1_BES0**0.5
n1_error=0.


nref = pars['nref']         #in 10^(19) /m^3


file=open("0BES.txt","w")
file.write('n1_BES='+str(n1_BES)+'/m^3\n')
print('n1_BES='+str(n1_BES)+'/m^3')
file.write('n0='+str(nref)+'*10^19/m^3')
print('n0='+str(nref)+'*10^19/m^3')

