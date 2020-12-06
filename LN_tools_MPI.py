import os
import sys
import numpy as np
import optparse as op
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.pyplot as plt
import imageio   #for making animation

from nrgWrapper import read_from_nrg_files
from parIOWrapper import read_species_gradients
from parIOWrapper import read_species_tempdens
from ParIO import Parameters 
from fieldlib import fieldfile
from geomWrapper import ky
from geomWrapper import init_read_geometry_file
from read_write_geometry import read_geometry_local
from read_iterdb_file import read_iterdb_file
from interp import interp
from FFT_general import FFT_function_time
from max_stat_tool import sort_x_f

from mpi4py import MPI


def B1_ky_f_spectrum_Z_sum(suffix,iterdb_file_name,min_Z0=-0.03,max_Z0=0.03,\
    Outboard_mid_plane=False,plot=False,show=False,csv_output=False,pic_path='pic',csv_path='csv'):
    #Import the parameters from parameter file using ParIO
    par = Parameters()
    par.Read_Pars('parameters'+suffix)
    pars = par.pardict
    
    nky0 = int(pars['nky0'])
    
    #*********Start of Get geometry from magn_geometry***********
    #getting B field using read_write_geometry.py
    gpars,geometry = read_geometry_local(pars['magn_geometry'][1:-1]+suffix)
    #Get geom_coeff from ParIO Wrapper
    geom_type, geom_pars, geom_coeff = init_read_geometry_file(suffix, pars)
    
    real_R=geometry['gl_R'] #it is major radius in meter
    real_Z=geometry['gl_z'] #it is height in meter, midland is 0, down is negative ,up is positive
    
    #*********End of Get geometry from magn_geometry***********
    
    time_start,time_end=start_end_time(suffix,pars)

    
    #******************************************************************
    #********Start of Determin the length of the frequency array****************
    
    iky=nky0-1 #Scan the biggest ky, which has the largest doppler shift
    inz=int(len(real_Z)/2.)
    frequency_kHZ, amplitude_frequency, amplitude_growth=\
        LN_apar_frequency_nz_iky(suffix,inz,iky,time_start,time_end,\
        plot=False,pic_path='pic',csv_path='csv',output_csv=False)

    len_freq=len(frequency_kHZ)
    print("len_freq="+str(len_freq))
    

    #********End of Determin the length of the frequency array****************
    #******************************************************************
    
    #******************************************************************
    #********Start of Determin the length of the nZ****************
    Z_list=[]
    nZ_list=[]
    for nZ in range(len(real_Z)):
        Z=real_Z[nZ]
        if min_Z0<Z and Z<=max_Z0:
            Z_list.append(real_Z[nZ])
            nZ_list.append(nZ)

    if Outboard_mid_plane==True:
        nZ_list=[int(len(real_Z)/2)]
        Z_list=[real_Z[int(len(real_Z)/2)]]

    len_nZ=len(nZ_list)


    print('nZ_list: '+str(nZ_list))
    #********End of Determin the length of the nZ****************
    #******************************************************************
    
    #***********Start of FFT*******************************
    frequency_kHZ_ky=np.zeros((len_nZ,nky0,len_freq))
    amplitude_frequency_ky=np.zeros((len_nZ,nky0,len_freq))
    amplitude_growth_ky=np.zeros((len_nZ,nky0,len_freq))

    for i_Z in range(len_nZ):
        inZ=nZ_list[i_Z]
        Z=real_Z[inZ]
        print('***************************')
        print('***************************')
        print('***************************')
        print('***************************')
        print('Looking at Z='+str(Z)+'m')
        #one can restrict the ky range corresprons to frequency 

        for i_ky in range(nky0):
            amplitude_frequency, amplitude_frequency, amplitude_growth=\
                LN_apar_frequency_nz_iky(suffix,inZ,i_ky,time_start,time_end,\
                plot=False,pic_path='pic',csv_path='csv',output_csv=False)

            omegaDoppler_kHZ=Doppler_calc(suffix,i_ky,iterdb_file_name)

            new_frequency_kHZ_ky, new_amplitude_frequency, \
                new_amplitude_growth=frequency_Doppler(frequency_kHZ,\
                amplitude_frequency,amplitude_growth,omegaDoppler_kHZ)


            frequency_kHZ_ky[i_Z,i_ky,:]=new_frequency_kHZ_ky
            amplitude_frequency_ky[i_Z,i_ky,:]=new_amplitude_frequency
            amplitude_growth_ky[i_Z,i_ky,:]=new_amplitude_growth

    #***********End of FFT*******************************


    #***********Start of Sum of Z--Length-Weighted*************************

    frequency_kHZ_ky_sum_Z=np.zeros((nky0,len_freq))
    amplitude_frequency_ky_sum_Z=np.zeros((nky0,len_freq))
    amplitude_growth_ky_sum_Z=np.zeros((nky0,len_freq))



    for i_ky in range(nky0):
        sum_length_TEMP=0.
        for i_Z in range(len(nZ_list)):
            nZ=nZ_list[i_Z]
            length=np.sqrt( (real_R[nZ]-real_R[nZ-1])**2. \
                   +(real_Z[nZ]-real_Z[nZ-1])**2. )
            sum_length_TEMP=sum_length_TEMP+length
            frequency_kHZ_ky_sum_Z[i_ky,:]=frequency_kHZ_ky_sum_Z[i_ky,:]+frequency_kHZ_ky[i_Z,i_ky,:]*length
            amplitude_frequency_ky_sum_Z[i_ky,:]=amplitude_frequency_ky_sum_Z[i_ky,:]+amplitude_frequency_ky[i_Z,i_ky,:]*length
            amplitude_growth_ky_sum_Z[i_ky,:]=amplitude_growth_ky_sum_Z[i_ky,:]+amplitude_growth_ky[i_Z,i_ky,:]*length

        frequency_kHZ_ky_sum_Z[i_ky,:]=frequency_kHZ_ky_sum_Z[i_ky,:]/sum_length_TEMP
        amplitude_frequency_ky_sum_Z[i_ky,:]=amplitude_frequency_ky_sum_Z[i_ky,:]/sum_length_TEMP
        amplitude_growth_ky_sum_Z[i_ky,:]=amplitude_growth_ky_sum_Z[i_ky,:]/sum_length_TEMP

    #***********End of Sum of Z--Length-Weighted*************************

    #*************End of Interperlation******************
    df_min=min(abs(frequency_kHZ_ky_sum_Z[0,:-1]-frequency_kHZ_ky_sum_Z[0,1:]))
    print("df_min: "+str(df_min))
    print("min(frequency_kHZ_ky_sum_Z): "+str(np.amin(frequency_kHZ_ky_sum_Z)))
    uni_freq = np.linspace(np.amin(frequency_kHZ_ky_sum_Z),\
                           np.amax(frequency_kHZ_ky_sum_Z),\
                           num=int(abs((np.amax(frequency_kHZ_ky_sum_Z)-np.amin(frequency_kHZ_ky_sum_Z))/df_min)*10.))
    len_uni_freq=len(uni_freq)                     
    
    print("len_uni_freq: "+str(len_uni_freq))
    print("frequency_kHZ_ky_sum_Z[i_ky,:]: "+str(len(frequency_kHZ_ky_sum_Z[i_ky,:])))
    print("amplitude_frequency_ky_sum_Z[i_ky,:]: "+str(len(amplitude_frequency_ky_sum_Z[i_ky,:])))

    frequency_kHZ_uni=np.zeros((nky0,len_uni_freq))
    amplitude_frequency_uni=np.zeros((nky0,len_uni_freq))
    new_frequency_kHZ_uni=np.zeros((nky0,len_uni_freq))


    for i_ky in range(nky0):
        frequency_kHZ_uni[i_ky,:]=uni_freq
        amplitude_frequency_uni[i_ky,:]=np.interp(uni_freq,frequency_kHZ_ky_sum_Z[i_ky,:],amplitude_frequency_ky_sum_Z[i_ky,:])
        new_frequency_kHZ_uni[i_ky,:]=np.interp(uni_freq,frequency_kHZ_ky_sum_Z[i_ky,:],frequency_kHZ_ky_sum_Z[i_ky,:])
            

    #*************end of Interperlation******************

    n_list,ky_list=ky_list_calc(suffix)
    
    #****************start of Output***********************
    for i_ky in range(len(n_list)):
        d = {'f(kHz)':frequency_kHZ_ky_sum_Z[i_ky,:],'B_R(Gauss)':amplitude_frequency_ky_sum_Z[i_ky,:]}
        df=pd.DataFrame(d, columns=['f(kHz)','B_R(Gauss)'])
        df.to_csv(csv_path+'/B_r_freq_n_'+str(n_list[i_ky])+'.csv',index=False)
    
        plt.clf()
        plt.plot(frequency_kHZ_ky_sum_Z[i_ky,:],amplitude_frequency_ky_sum_Z[i_ky,:],label='Original')
        plt.xlabel('frequency (kHz)')
        plt.ylabel('B_r(Gauss)') 
        plt.title(r'$\bar{B}_r$(Gauss) spectrogram of n='+str(n_list[i_ky]),fontsize=18)
        plt.savefig(pic_path+'/B_r_freq_n_'+str(n_list[i_ky])+'.png')
        
        plt.clf()
        plt.plot(frequency_kHZ_ky_sum_Z[i_ky,:],amplitude_frequency_ky_sum_Z[i_ky,:],label='Original')
        plt.plot(frequency_kHZ_uni[i_ky,:],amplitude_frequency_uni[i_ky,:],label='Intper',alpha=0.5)
        plt.legend()
        plt.xlabel('frequency (kHz)')
        plt.ylabel('B_r(Gauss)') 
        plt.title(r'$\bar{B}_r$(Gauss) spectrogram of n='+str(n_list[i_ky]),fontsize=18)
        plt.savefig(pic_path+'/Interp_B_r_freq_n_'+str(n_list[i_ky])+'.png')

    
    B1_plot=amplitude_frequency_ky_sum_Z
    f_plot=frequency_kHZ_ky_sum_Z
    ky_plot=np.zeros(np.shape(frequency_kHZ_ky_sum_Z))
    
    for i_ky in range(nky0):
        ky_plot[i_ky,:]=  [ky_list[i_ky]] *len(frequency_kHZ_ky_sum_Z[i_ky,:])

    B1_plot=np.transpose(B1_plot)
    f_plot=np.transpose(f_plot)
    ky_plot=np.transpose(ky_plot)
    
    #print('shape'+str(np.shape(np.transpose(frequency_kHZ_ky_sum_Z))))

    
    if plot==True:
        plt.clf()
        plt.ylabel(r'$k_y \rho_i$',fontsize=10)
        plt.xlabel(r'$f(kHz)$',fontsize=10)
        plt.contourf(f_plot,ky_plot,np.log(B1_plot))#,level=[50,50,50])#,cmap='RdGy')
        for ky in ky_list:
            plt.axhline(ky,color='red',alpha=0.5)#alpha control the transparency, alpha=0 transparency, alpha=1 solid
        plt.axhline(ky_list[0],color='red',alpha=0.5,label='n starts from'+str(n_list[0]) )#alpha control the transparency, alpha=0 transparency, alpha=1 solid
        plt.legend()
        plt.colorbar()
        plt.title(r'log($B_r$) contour plot',fontsize=10)
        plt.savefig('0B_r_contour.png')
        #plt.show()


    #****************end of Output***********************
    


    #**********start of Sum over ky*********************
    
    amplitude_frequency_uni_ky_sum=np.sum(amplitude_frequency_uni,axis=0)
    #amplitude_growth_uni_ky_sum=np.sum(amplitude_growth_uni,axis=0)
    
    #**********end of Sum over ky*********************
    if csv_output==True:
        d = {'f(kHz)':uni_freq,'B_R(Gauss)':amplitude_frequency_uni_ky_sum}
        df=pd.DataFrame(d, columns=['f(kHz)','B_R(Gauss)'])
        df.to_csv('0B_r_from_'+str(max_Z0)+'_to_'+str(min_Z0)+'.csv',index=False)

    if plot==True:
        plt.clf()
        plt.xlabel(r'$Frequency(kHz)$',fontsize=15)
        plt.ylabel(r'$\bar{B}_r(Gauss)$',fontsize=15)
        plt.plot(uni_freq,amplitude_frequency_uni_ky_sum,label='Inpter')
        #plt.legend()
        plt.title(r'$\bar{B}_r$(Gauss) spectrogram',fontsize=18)
        plt.savefig('0RIP_freq_spectrogram.png')


    if plot==True: 
        plt.clf()
        plt.plot(uni_freq,amplitude_frequency_uni_ky_sum)
        plt.xlabel('frequency (kHz)')
        plt.ylabel('B_r(Gauss)') 
        plt.title(r'$\bar{B}_r$(Gauss) spectrogram',fontsize=18)
        plt.savefig('0RIP_freq_spectrogram.png')
        if show==True:
            plt.show()


    return uni_freq,amplitude_frequency_uni_ky_sum,amplitude_frequency_uni
    #uni_freq is 1D array
    #amplitude_frequency_uni_ky_sum is 1D array B1(f) length weighted average over Z, sum of ky
    #amplitude_frequency_uni is 2D array B1(ky,f) length weighted average over Z


#uni_freq,amplitude_frequency_uni_ky_sum,amplitude_frequency_uni=\
#        B1_ky_f_spectrum_Z_sum_time_set(suffix,iterdb_file_name,\
#            min_Z0,max_Z0,Outboard_mid_plane,\
#            time_start,time_end,\
#            plot,show,csv_output,pic_path,csv_path)
#uni_freq is 1D array
#amplitude_frequency_uni_ky_sum is 1D array B1(f) length weighted average over Z, sum of ky
#amplitude_frequency_uni is 2D array B1(ky,f) length weighted average over Z