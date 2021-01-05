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
from momentsWrapper_max import LILO_moments_from_mom_file

#input the suffix , plot nrg_es, em, return time_start,time_end
#from LN_tools import start_end_time
#time_start,time_end=start_end_time(suffix,pars)
def get_suffix():
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
    return suffix

def start_end_time(suffix,pars):  #suffix in the format of "_1" or ".dat"


    if pars['n_spec'] == 1:
        time, nrge = read_from_nrg_files(pars,suffix,False)
    elif pars['n_spec'] == 2:
        time, nrgi, nrge = read_from_nrg_files(pars,suffix,False)
    elif pars['n_spec'] == 3:
        time, nrgi, nrge, nrgz = read_from_nrg_files(pars,suffix,False)

    plt.clf()
    plt.plot(time,nrge[:,6],label="Q_es of electron")
    plt.plot(time,nrge[:,7],label="Q_em of electron")
    plt.title('nrg of electron')
    plt.xlabel('time')
    plt.legend()
    plt.show()
    
    scan_all = str(input("Scan all(Y/N):\n"))
    if scan_all=='n' or scan_all=='N':
        time_start = float(input("Start time:\n"))
        time_end = float(input("End time:\n"))
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

# n_list, ky_list=ky_list_calc(suffix)
def ky_list_calc(suffix):
    #Import the parameters from parameter file using ParIO
    par = Parameters()
    par.Read_Pars('parameters'+suffix)
    pars = par.pardict

    B_gauss=10.**4              #1T=10^4Gauss
    qref = 1.6E-19              #in C
    c  = 1.                     #in 3*10^8m/s
    m_kg = 1.673E-27            #in kg
    Bref = pars['Bref']         #in Tesla
    Tref = pars['Tref']         #in keV
    nref = pars['nref']         #in 10^(19) /m^3
    Lref = pars['Lref']         #in m
    mref = pars['mref']         #in proton mass(kg)
    q0 = pars['q0']              #unitless, safety factor/q95
    x0 = pars['x0']             #x/a, location
    kymin = pars['kymin']       #in rhoi
    nky0 = pars['nky0']         #in total number of ky
    n_step = pars['n0_global']  #in rhoi
    nref = nref * 1.E19         #in the unit of /m^3
    Tref = Tref * qref * 1.E03  #in the unit of J
    mref = mref * m_kg          #in the unit of kg
    pref = nref * Tref          #in Pa*k_{B}
    cref = np.sqrt(Tref / mref) #in the unit of m/s
    Omegaref = qref * Bref / mref / c  #in rad/s
    rhoref = cref / Omegaref           #in m rho_i(ion gyroradii)
    rhorefStar = rhoref / Lref         #Unitless
    gyroFreq= cref/Lref                 #the facor convert frequency from cs/a to rad/s


    #Determine the n for kymin
    ky_n1=1.*q0*rhoref/(Lref*x0)  #ky for n=1
    n_min=round(kymin/ky_n1)   #n for ky_min
    if n_min==0:
        n_min=1
    n_list=np.arange(int(n_min),int(n_min+nky0*n_step),int(n_step))
    ky_list=ky_n1*n_list
    return n_list, ky_list

#************Sample function line
#omegaDoppler=Doppler_calc(suffix,iky,iterdb_file_name)
def Doppler_calc(suffix,iky,iterdb_file_name):
    #Import the parameters from parameter file using ParIO
    par = Parameters()
    par.Read_Pars('parameters'+suffix)
    pars = par.pardict

    B_gauss=10.**4              #1T=10^4Gauss
    qref = 1.6E-19              #in C
    c  = 1.                     #in 3*10^8m/s
    m_kg = 1.673E-27            #in kg
    Bref = pars['Bref']         #in Tesla
    Tref = pars['Tref']         #in keV
    nref = pars['nref']         #in 10^(19) /m^3
    Lref = pars['Lref']         #in m
    mref = pars['mref']         #in proton mass(kg)
    q0 = pars['q0']              #unitless, safety factor/q95
    x0 = pars['x0']             #x/a, location
    kymin = pars['kymin']       #in rhoi
    nky0 = pars['nky0']         #in total number of ky
    n_step = pars['n0_global']  #in rhoi
    nref = nref * 1.E19         #in the unit of /m^3
    Tref = Tref * qref * 1.E03  #in the unit of J
    mref = mref * m_kg          #in the unit of kg
    pref = nref * Tref          #in Pa*k_{B}
    cref = np.sqrt(Tref / mref) #in the unit of m/s
    Omegaref = qref * Bref / mref / c  #in rad/s
    rhoref = cref / Omegaref           #in m rho_i(ion gyroradii)
    rhorefStar = rhoref / Lref         #Unitless
    gyroFreq= cref/Lref                 #the facor convert frequency from cs/a to rad/s
    
    #Determine the n for kymin
    ky_n1=1.*q0*rhoref/(Lref*x0)  #ky for n=1
    n_min=round(kymin/ky_n1)   #n for ky_min
    if n_min==0:
        n_min=1
    n_list=np.arange(int(n_min),int(n_min+nky0*n_step),int(n_step))

    #**********************Doppler shift**********************************************

    rhot0, te0, ti0, ne0, ni0, nz0, vrot0 = read_iterdb_file(iterdb_file_name)
    uni_rhot = np.linspace(min(rhot0),max(rhot0),int(len(rhot0)*10.))
    x0_index=np.argmin(abs(uni_rhot-x0))
    vrot_u = interp(rhot0,vrot0,uni_rhot)
    omegaDoppler_kHZ = vrot_u[x0_index]*n_list[iky]/(2.*np.pi*1000.)

    #**********************Doppler shift**********************************************

    return omegaDoppler_kHZ

#from LN_tools import frequency_Doppler
#new_frequency_kHZ, new_amplitude_frequency, new_amplitude_growth=frequency_Doppler(frequency_kHZ,amplitude_frequency,amplitude_growth,omegaDoppler_kHz)
def frequency_Doppler(frequency_kHZ,amplitude_frequency,amplitude_growth,omegaDoppler_kHZ):
    #!!!!!!!!!!Debatable
    #!!!!!!!!!!Debatable
    #!!!!!!!!!!Debatable
    #new_amplitude_frequency=amplitude_frequency[::-1]+amplitude_frequency
    new_amplitude_frequency=amplitude_frequency[::-1]
    new_amplitude_growth=amplitude_growth
    new_frequency_kHZ=frequency_kHZ+omegaDoppler_kHZ
    #!!!!!!!!!!Debatable
    #!!!!!!!!!!Debatable
    #!!!!!!!!!!Debatable
    return new_frequency_kHZ, new_amplitude_frequency, new_amplitude_growth
    
  

def BES_f_spectrum(suffix,iterdb_file_name,min_Z0,max_Z0,\
    Outboard_mid_plane,time_step,time_start,time_end,\
    plot,show,csv_output,pic_path,csv_path):
    #Where inz is the number of the element in nz

    #Import the parameters from parameter file using ParIO
    par = Parameters()
    par.Read_Pars('parameters'+suffix)
    pars = par.pardict
    #getting B field using read_write_geometry.py
    gpars,geometry = read_geometry_local(pars['magn_geometry'][1:-1]+suffix)
    #Get geom_coeff from ParIO Wrapper
    geom_type, geom_pars, geom_coeff = init_read_geometry_file(suffix, pars)
    real_Z=geometry['gl_z'] 
    real_R=geometry['gl_R'] 
    #Import the field file using fieldlib
    field = fieldfile('field'+suffix,pars)
    time = np.array(field.tfld)  #time stampes

    B_gauss=10.**4              #1T=10^4Gauss
    qref = 1.6E-19              #in C
    c  = 1.                     #in 3*10^8m/s
    m_kg = 1.673E-27            #in kg
    Bref = pars['Bref']         #in Tesla
    Tref = pars['Tref']         #in keV
    nref = pars['nref']         #in 10^(19) /m^3
    Lref = pars['Lref']         #in m
    mref = pars['mref']         #in proton mass(kg)
    q0 = pars['q0']              #unitless, safety factor/q95
    x0 = pars['x0']             #x/a, location
    kymin = pars['kymin']       #in rhoi
    nky0 = pars['nky0']         #in total number of ky
    n_step = pars['n0_global']  #in rhoi
    nref = nref * 1.E19         #in the unit of /m^3
    Tref = Tref * qref * 1.E03  #in the unit of J
    mref = mref * m_kg          #in the unit of kg
    pref = nref * Tref          #in Pa*k_{B}
    cref = np.sqrt(Tref / mref) #in the unit of m/s
    Omegaref = qref * Bref / mref / c  #in rad/s
    rhoref = cref / Omegaref           #in m rho_i(ion gyroradii)
    rhorefStar = rhoref / Lref         #Unitless
    gyroFreq= cref/Lref                 #the facor convert frequency from cs/a to rad/s


    #ky comes from geomWrapper.py
    ky_GENE_temp=ky(pars, geom_coeff,plot=False)     #ky for ky min for differen z
    #print('ky shape: '+str(np.shape(ky_GENE_temp)))

    n_list,ky_list=ky_list_calc(suffix)
    print('n0 list length: '+str(len(n_list)))
    print('n0 list: '+str(n_list))
    
    time_start_index=np.argmin(abs(time - time_start))
    time_end_index=np.argmin(abs(time - time_end))

    time_list_temp = time[time_start_index:time_end_index+1]
    time_list=[]

    for i in range(0,len(time_list_temp),time_step):
        time_list.append(time_list_temp[i])
    time_list=np.array(time_list)
    
    nky0=len(n_list)
    ntime=len(time_list)
    n1_ky_t=np.zeros((nky0,ntime))

    if os.path.isdir(csv_path):  #if path does not exist, then create 'csv'
        pass
    else:
        os.mkdir(csv_path) 
    if os.path.isdir(pic_path):  #if path does not exist, then create 'pic'
        pass
    else:
        os.mkdir(pic_path) 
    print("**********Scan starts, output in csv and pic***************")

    for i in range(len(time_list)):
        time0=time_list[i]
        itime = np.argmin(abs(time - time0))
        
        print("Looking at the spectra at time:"+str(time[itime]))
        #This sets the time step you want to read in
        #field.set_time(time[itime])
        
        if 'x_local' in pars:
            if pars['x_local']:
                x_local = True
            else:
                x_local = False 
        else:
            x_local = True

        if x_local:
            kxmin = 2.0*np.pi/pars['lx']
            kxgrid = np.linspace(-(pars['nx0']/2-1)*kxmin,pars['nx0']/2*kxmin,num=pars['nx0'])
            kxgrid = np.roll(kxgrid,int(pars['nx0']/2+1))
            #print("kxgrid"+str(kxgrid))
            
            #print("kygrid"+str(kygrid))
            zgrid = np.linspace(-np.pi,np.pi,pars['nz0'],endpoint=False)  

            

            upar,deln,deln_global= LILO_moments_from_mom_file(pars,suffix,False,setTime=itime)
            n1_GENE=abs(deln_global[:,:,:])
            n1_GENE_ky0=np.sum(n1_GENE[:,:,:],axis=2) 
            n1_GENE_ky= np.zeros(np.shape(n1_GENE_ky0[0,:]))

            #*****Sum over Z************
            sum_length_TEMP=0
            for nZ in range(len(real_Z)):
                if min_Z0<=real_Z[nZ] and real_Z[nZ]<=max_Z0:
                    length=np.sqrt( (real_R[nZ]-real_R[nZ-1])**2.\
                        +(real_Z[nZ]-real_Z[nZ-1])**2. )
                    sum_length_TEMP=sum_length_TEMP+length
                    n1_GENE_ky=n1_GENE_ky+n1_GENE_ky0[nZ,:]*length
            #*****Sum over Z************
            n1_GENE_ky=n1_GENE_ky/sum_length_TEMP
            n1_ky=n1_GENE_ky*(rhorefStar)*nref

        else:  #x_local = False
            print("Sorry, cannot handle Global Nonlinear yet...")
            pass
        #**Finished reading the Br
        #Recall B1_ky_t_inz=np.zeros((nky0,ntime))

        n1_ky_t[:,i]=n1_ky
    
    amplitude_frequency_sum=0
    amplitude_growth_sum=0
    #print(str(B1_ky_t_inz))
    if plot==True:
        ims_n1=[]

    n1_ky_f=[]
    growth_ky_f=[]
    f_ky_f=[]
    for iky in range(nky0):
        n1_inz_t=n1_ky_t[iky,:]
        frequency,amplitude_frequency,amplitude_growth=FFT_function_time(n1_inz_t,time_list,plot=False)
        frequency_kHZ=frequency*gyroFreq/(1000.)

        omegaDoppler_kHZ=Doppler_calc(suffix,iky,iterdb_file_name)

        new_frequency_kHZ_ky, new_amplitude_frequency, \
        new_amplitude_growth=frequency_Doppler(frequency_kHZ,\
        amplitude_frequency,amplitude_growth,omegaDoppler_kHZ)
        
        f_ky_f.append(new_frequency_kHZ_ky)
        n1_ky_f.append(new_amplitude_frequency)
        growth_ky_f.append(new_amplitude_growth)

    k_f_plot(f_ky_f,n1_ky_f,growth_ky_f,ky_list,n_list,pic_path,csv_path,name='n1')

    return f_ky_f,n1_ky_f,growth_ky_f


def RIP_f_spectrum(suffix,iterdb_file_name,min_Z0,max_Z0,\
    Outboard_mid_plane,time_step,time_start,time_end,\
    plot,show,csv_output,pic_path,csv_path):
    #Where inz is the number of the element in nz

    #Import the parameters from parameter file using ParIO
    par = Parameters()
    par.Read_Pars('parameters'+suffix)
    pars = par.pardict
    #getting B field using read_write_geometry.py
    gpars,geometry = read_geometry_local(pars['magn_geometry'][1:-1]+suffix)
    #Get geom_coeff from ParIO Wrapper
    geom_type, geom_pars, geom_coeff = init_read_geometry_file(suffix, pars)
    real_Z=geometry['gl_z'] 
    real_R=geometry['gl_R'] 
    #Import the field file using fieldlib
    field = fieldfile('field'+suffix,pars)
    time = np.array(field.tfld)  #time stampes

    B_gauss=10.**4.              #1T=10^4Gauss
    qref = 1.6E-19              #in C
    c  = 1.                     #in 3*10^8m/s
    m_kg = 1.673E-27            #in kg
    Bref = pars['Bref']         #in Tesla
    Tref = pars['Tref']         #in keV
    nref = pars['nref']         #in 10^(19) /m^3
    Lref = pars['Lref']         #in m
    mref = pars['mref']         #in proton mass(kg)
    q0 = pars['q0']              #unitless, safety factor/q95
    x0 = pars['x0']             #x/a, location
    kymin = pars['kymin']       #in rhoi
    nky0 = pars['nky0']         #in total number of ky
    n_step = pars['n0_global']  #in rhoi
    nref = nref * 1.E19         #in the unit of /m^3
    Tref = Tref * qref * 1.E03  #in the unit of J
    mref = mref * m_kg          #in the unit of kg
    pref = nref * Tref          #in Pa*k_{B}
    cref = np.sqrt(Tref / mref) #in the unit of m/s
    Omegaref = qref * Bref / mref / c  #in rad/s
    rhoref = cref / Omegaref           #in m rho_i(ion gyroradii)
    rhorefStar = rhoref / Lref         #Unitless
    gyroFreq= cref/Lref                 #the facor convert frequency from cs/a to rad/s


    #ky comes from geomWrapper.py
    ky_GENE_temp=ky(pars, geom_coeff,plot=False)     #ky for ky min for differen z
    #print('ky shape: '+str(np.shape(ky_GENE_temp)))

    #Determine the n for kymin
    ky_n1=1.*q0*rhoref/(Lref*x0)  #ky for n=1
    n_list,ky_list=ky_list_calc(suffix)
    print('n0 list length: '+str(len(n_list)))
    print('n0 list: '+str(n_list))
    n_min=np.min(n_list)
    ky_GENE_n1=ky_GENE_temp/float(n_min)
    ky_GENE_grid = np.outer(ky_GENE_n1,n_list) #outer product of the two vectors
    
    print("kygrid"+str(np.shape(ky_GENE_grid)))
    
    print('n0 list length: '+str(len(n_list)))
    print('n0 list: '+str(n_list))
    
    
    #B1=abs(np.mean(Apar_GENE[z,:])*len(Apar_GENE[z,:])*(ky_GENE_temp[z]/rhoref)*Bref*B_gauss*rhorefStar*rhoref)
    Apar_to_B1=abs((1./rhoref)*Bref*B_gauss*rhorefStar*rhoref)         #B1=Apar*ky_GENE_temp*Apar_to_B1

    
    time_start_index=np.argmin(abs(time - time_start))
    time_end_index=np.argmin(abs(time - time_end))

    time_list_temp = time[time_start_index:time_end_index+1]
    time_list=[]

    for i in range(0,len(time_list_temp),time_step):
        time_list.append(time_list_temp[i])
    time_list=np.array(time_list)
    
    nky0=len(n_list)
    ntime=len(time_list)
    B1_ky_t=np.zeros((nky0,ntime))

    if os.path.isdir(csv_path):  #if path does not exist, then create 'csv'
        pass
    else:
        os.mkdir(csv_path) 
    if os.path.isdir(pic_path):  #if path does not exist, then create 'pic'
        pass
    else:
        os.mkdir(pic_path) 
    print("**********Scan starts, output in csv and pic***************")

    for i in range(len(time_list)):
        time0=time_list[i]
        itime = np.argmin(abs(time - time0))
        field.set_time(time_list[i])
        print("Looking at the spectra at time:"+str(time[itime]))
        #This sets the time step you want to read in
        #field.set_time(time[itime])
        
        if 'x_local' in pars:
            if pars['x_local']:
                x_local = True
            else:
                x_local = False 
        else:
            x_local = True

        if x_local:
            kxmin = 2.0*np.pi/pars['lx']
            kxgrid = np.linspace(-(pars['nx0']/2-1)*kxmin,pars['nx0']/2*kxmin,num=pars['nx0'])
            kxgrid = np.roll(kxgrid,int(pars['nx0']/2+1))
            #print("kxgrid"+str(kxgrid))
            
            #print("kygrid"+str(kygrid))
            zgrid = np.linspace(-np.pi,np.pi,pars['nz0'],endpoint=False)  

            upar,deln,deln_global= LILO_moments_from_mom_file(pars,suffix,False,setTime=itime)
            apar=field.apar()[:,:,:]
            apar_ky=np.sum(apar[:,:,:],axis=2) 
            B1_GENE_ky0=ky_GENE_grid*apar_ky*Apar_to_B1
            B1_GENE_ky= np.zeros(np.shape(B1_GENE_ky0[0,:]))

            #*****Sum over Z************
            sum_length_TEMP=0
            for nZ in range(len(real_Z)):
                if min_Z0<=real_Z[nZ] and real_Z[nZ]<=max_Z0:
                    length=np.sqrt( (real_R[nZ]-real_R[nZ-1])**2.\
                        	+(real_Z[nZ]-real_Z[nZ-1])**2. )
                    sum_length_TEMP=sum_length_TEMP+length
                    B1_GENE_ky=B1_GENE_ky+B1_GENE_ky0[nZ,:]*length
            #*****Sum over Z************
            B1_GENE_ky=B1_GENE_ky/sum_length_TEMP
            B1_ky=B1_GENE_ky

        else:  #x_local = False
            print("Sorry, cannot handle Global Nonlinear yet...")
            pass
        #**Finished reading the Br
        #Recall B1_ky_t_inz=np.zeros((nky0,ntime))
        B1_ky_t[:,i]=B1_ky
    
    ky_GENE_inz = ky_GENE_grid[int(len(real_Z)/2),:]
    
    amplitude_frequency_sum=0
    amplitude_growth_sum=0
    #print(str(B1_ky_t_inz))
    if plot==True:
        ims_n1=[]

    B1_ky_f=[]
    growth_ky_f=[]
    f_ky_f=[]
    for iky in range(nky0):
        B1_inz_t=B1_ky_t[iky,:]
        frequency,amplitude_frequency,amplitude_growth=FFT_function_time(B1_inz_t,time_list,plot=False)
        frequency_kHZ=frequency*gyroFreq/(1000.)

        omegaDoppler_kHZ=Doppler_calc(suffix,iky,iterdb_file_name)

        new_frequency_kHZ_ky, new_amplitude_frequency, \
        new_amplitude_growth=frequency_Doppler(frequency_kHZ,\
        amplitude_frequency,amplitude_growth,omegaDoppler_kHZ)
        
        f_ky_f.append(new_frequency_kHZ_ky)
        B1_ky_f.append(new_amplitude_frequency)
        growth_ky_f.append(new_amplitude_growth)
    
    
    f_ky_f=np.array(f_ky_f)
    B1_ky_f=abs(np.array(B1_ky_f))
    growth_ky_f=np.array(growth_ky_f)

    f,amplitude_f=k_f_plot(f_ky_f,B1_ky_f,growth_ky_f,ky_list,n_list,pic_path,csv_path,name='B1')
    
    return f_ky_f,B1_ky_f,growth_ky_f


def k_f_plot(f_ky_f,amplitude_ky_f,growth_ky_f,ky_list,n_list,pic_path,csv_path,name='0'):
    f_ky_f=np.array(f_ky_f)
    amplitude_ky_f=np.array(amplitude_ky_f)
    growth_ky_f=np.array(growth_ky_f)

    (nky0,len_f)=np.shape(f_ky_f)

    df_min=np.min(abs(f_ky_f[0,:-1]-f_ky_f[0,1:]))
    uni_freq = np.linspace(np.min(f_ky_f), np.max(f_ky_f),\
                           num=int(abs((np.amax(f_ky_f)-np.amin(f_ky_f))/df_min)/10.)*100)
    print('len of uni_freq='+str(int(abs((np.amax(f_ky_f)-np.amin(f_ky_f))/df_min)/10.)*100))
    len_uni_freq=len(uni_freq)
    frequency_kHZ_uni=np.zeros((nky0,len_uni_freq))
    amplitude_frequency_uni=np.zeros((nky0,len_uni_freq))
    growth_frequency_kHZ_uni=np.zeros((nky0,len_uni_freq))
    ky_plot=np.zeros((nky0,len_uni_freq))
    for i_ky in range(nky0):
        frequency_kHZ_uni[i_ky,:]=uni_freq
        amplitude_frequency_uni[i_ky,:]=np.interp(uni_freq,f_ky_f[i_ky,:],amplitude_ky_f[i_ky,:])
        #growth_frequency_kHZ_uni[i_ky,:]=np.interp(uni_freq,f_ky_f[i_ky,:],growth_ky_f[i_ky,:])
        plt.clf()
        plt.xlabel(r'$f(kHz)$',fontsize=10)
        plt.ylabel(str(name),fontsize=10)
        plt.plot(frequency_kHZ_uni[i_ky,:],amplitude_frequency_uni[i_ky,:],label='intper')
        plt.plot(f_ky_f[i_ky,:],amplitude_ky_f[i_ky,:],label='original')
        plt.legend()
        plt.savefig(pic_path+'/0n='+str(n_list[i_ky])+'.png')
    ky_plot2=np.zeros((nky0,len_f))
    for i_ky in range(nky0):
        ky_plot[i_ky,:]=[ky_list[i_ky]]*len_uni_freq
        ky_plot2[i_ky,:]=[ky_list[i_ky]]*len_f

    plt.clf()
    plt.ylabel(r'$k_y \rho_s$',fontsize=10)
    plt.xlabel(r'$f(kHz)$',fontsize=10)
    plt.contourf(f_ky_f,ky_plot2,np.log(amplitude_ky_f))#,level=[50,50,50])#,cmap='RdGy')
    for ky in ky_list:
        plt.axhline(ky,color='red',alpha=0.5)#alpha control the transparency, alpha=0 transparency, alpha=1 solid
    plt.axhline(ky_list[0],color='red',alpha=0.5,label='n starts from'+str(n_list[0]) )#alpha control the transparency, alpha=0 transparency, alpha=1 solid
    plt.legend()
    plt.colorbar()
    plt.title('log('+str(name)+') contour plot',fontsize=10)
    plt.savefig(pic_path+'/'+str(name)+'_log_contour_plot.png')

    plt.clf()
    plt.ylabel(r'$k_y \rho_s$',fontsize=10)
    plt.xlabel(r'$f(kHz)$',fontsize=10)
    plt.contourf(f_ky_f,ky_plot2,amplitude_ky_f)#,level=[50,50,50])#,cmap='RdGy')
    for ky in ky_list:
        plt.axhline(ky,color='red',alpha=0.5)#alpha control the transparency, alpha=0 transparency, alpha=1 solid
    plt.axhline(ky_list[0],color='red',alpha=0.5,label='n starts from'+str(n_list[0]) )#alpha control the transparency, alpha=0 transparency, alpha=1 solid
    plt.legend()
    plt.colorbar()
    plt.title(str(name)+' contour plot',fontsize=10)
    plt.savefig(pic_path+'/'+str(name)+'_contour_plot.png')

    #plot the y(freq) 
    f=frequency_kHZ_uni[0,:]
    amplitude_f=np.sum(amplitude_frequency_uni,axis=0)

    plt.clf()
    plt.plot(f,amplitude_f)
    plt.ylabel(str(name),fontsize=10)
    plt.xlabel(r'$f(kHz)$',fontsize=10)
    plt.savefig(pic_path+'/0'+str(name)+'_spectrum.png')

    d = {'f(kHz)':f,'B_R(Gauss)':amplitude_f}
    df=pd.DataFrame(d, columns=['f(kHz)','B_R(Gauss)'])
    df.to_csv(csv_path+'/0RIP_spectrum_freq.csv',index=False)

    return f,amplitude_f

def LN_BES(frequency_min,frequency_max,iterdb_file_name,Delta_Z=0.02,time_step=1,pic_path='BES_pic',csv_path='BES_csv'):
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

    f_ky_f,n1_ky_f,growth_ky_f=\
    BES_f_spectrum(suffix,iterdb_file_name,min_Z0,max_Z0,\
        Outboard_mid_plane,time_step,time_start,time_end,\
        plot,show,csv_output,pic_path,csv_path)
    (n_ky,n_f)=np.shape(f_ky_f)
    print(np.shape(f_ky_f))
    print(np.max(f_ky_f))
    print(np.min(f_ky_f))
    print(n_ky)
    print(n_f)
    n1_BES0=0.
    n1_BES_temp=0.
    int_f=0.


    for iky in range(n_ky):
        for i_f in range(n_f):
            if frequency_min<=f_ky_f[iky][i_f] and f_ky_f[iky][i_f]<=frequency_max:
                int_f=int_f+abs(f_ky_f[iky][i_f]-f_ky_f[iky][i_f-1])
                n1_BES0=n1_BES0+abs(n1_ky_f[iky][i_f])*abs(f_ky_f[iky][i_f]-f_ky_f[iky][i_f-1])
                n1_BES_temp=n1_BES_temp+abs(n1_ky_f[iky][i_f])

    nref = pars['nref']         #in 10^(19) /m^3
    print('n1_BES_temp='+str(n1_BES_temp))
    print('nref='+str(nref))
    n1_BES=n1_BES0/int_f

    file=open("0BES.txt","w")
    file.write('n1_BES='+str(n1_BES)+'/m^3\n')
    print('n1_BES='+str(n1_BES)+'/m^3')
    file.write('n0='+str(nref)+'*10^19/m^3')
    print('n0='+str(nref)+'*10^19/m^3')

