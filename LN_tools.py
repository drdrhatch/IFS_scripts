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
from FFT_general import FFT_function_time


#input the suffix , plot nrg_es, em, return time_start,time_end
#from LN_tools import start_end_time
#time_start,time_end=start_end_time(suffix,pars)
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


#input the suffix and nz index, time_start and end, return the Br as function of frequency(frequency_kHZ,amplitude_frequency_sum,amplitude_growth_sum)
#from LN_tools import LN_apar_frequency_nz
#frequency_kHZ,amplitude_frequency_sum,amplitude_growth_sum=LN_apar_frequency(suffix,inz,time_start,time_end,plot=False,pic_path='pic',csv_path='csv',output_csv=False)
def LN_apar_frequency_nz(suffix,inz,time_start,time_end,plot=False,pic_path='pic',csv_path='csv',output_csv=False):
    #Where inz is the number of the element in nz

    #Import the parameters from parameter file using ParIO
    par = Parameters()
    par.Read_Pars('parameters'+suffix)
    pars = par.pardict
    #getting B field using read_write_geometry.py
    gpars,geometry = read_geometry_local(pars['magn_geometry'][1:-1]+suffix)
    #Get geom_coeff from ParIO Wrapper
    geom_type, geom_pars, geom_coeff = init_read_geometry_file(suffix, pars)
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

    #Determine the n for kymin
    ky_n1=1.*q0*rhoref/(Lref*x0)  #ky for n=1
    n_min=round(kymin/ky_n1)   #n for ky_min
    if n_min==0:
        n_min=1
    n_list=np.arange(int(n_min),int(n_min+nky0*n_step),int(n_step))

    ky_GENE_n1=ky_GENE_temp/float(n_min)
    ky_GENE_grid = np.outer(ky_GENE_n1,n_list) #outer product of the two vectors
    
    print("kygrid"+str(np.shape(ky_GENE_grid)))
    
    print('n0 list length: '+str(len(n_list)))
    print('n0 list: '+str(n_list))
    
    
    #B1=abs(np.mean(Apar_GENE[z,:])*len(Apar_GENE[z,:])*(ky_GENE_temp[z]/rhoref)*Bref*B_gauss*rhorefStar*rhoref)
    Apar_to_B1=abs((1./rhoref)*Bref*B_gauss*rhorefStar*rhoref)         #B1=Apar*ky_GENE_temp*Apar_to_B1

    
    time_start_index=np.argmin(abs(time - time_start))
    time_end_index=np.argmin(abs(time - time_end))
    time_list = time[time_start_index:time_end_index+1]
    
    nky0=len(n_list)
    ntime=len(time_list)
    B1_ky_t_inz=np.zeros((nky0,ntime))

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
        field.set_time(time[itime])
        
        
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
    
            apar=field.apar()[:,:,:]
            #print("apar"+str(np.shape(apar)))
            apar_ky = np.sum(apar,axis=2)         #sum over x_axis
            (nz0,nky0)=np.shape(apar_ky)
            B1_ky=np.zeros(np.shape(apar_ky))
            #print("apar_ky"+str(np.shape(apar_ky)))

            B1_ky=ky_GENE_grid*apar_ky*Apar_to_B1 #B1 in Gauss  (nz0,nky0)*(nz0,nky0)*scaler
            
        else:  #x_local = False
            print("Sorry, cannot handle Global Nonlinear yet...")
            pass
        #**Finished reading the Br
        #Recall B1_ky_t_inz=np.zeros((nky0,ntime))

        (nz,nky)=np.shape(B1_ky)
        B1_ky_t_inz[:,i]=B1_ky[inz,:]
    
    ky_GENE_inz = ky_GENE_grid[inz,:]
    
    amplitude_frequency_sum=0
    amplitude_growth_sum=0
    #print(str(B1_ky_t_inz))
    if plot==True:
        ims_B1=[]
    for iky in range(nky):
        B1_inz_t=B1_ky_t_inz[iky,:]
        frequency,amplitude_frequency,amplitude_growth=FFT_function_time(B1_inz_t,time_list,plot=False)
        frequency_kHZ=frequency*gyroFreq/(1000.)

        if plot==True:
            plt.clf()
            plt.plot(frequency_kHZ,amplitude_frequency,label='frequency')
            plt.plot(frequency_kHZ,amplitude_growth,label='growth')
            plt.title(r'B_r(Gauss)'+' at ky_out_board= '+str(ky_GENE_inz[iky]))
            plt.xlabel(r'$f(kHz)$')       
            plt.ylabel(r'$B_r(Gauss)$')
            plt.legend()
            plt.savefig('pic/B_r_ky_out_board='+str(ky_GENE_inz[iky])+'.png')
            

        if plot==True:
            ims_B1.append(imageio.imread('pic/B_r_ky_out_board='+str(ky_GENE_inz[iky])+'.png'))



        amplitude_frequency_sum=amplitude_frequency_sum+abs(amplitude_frequency)
        amplitude_growth_sum=amplitude_growth_sum+amplitude_growth

        if output_csv==True:
            d = {'ky_out_board':[ky_GENE_inz[iky]]*len(frequency_kHZ),'frequency(kHZ)':frequency_kHZ,'B_R(Gauss)amplitude_frequency':amplitude_frequency,'B_R(Gauss)amplitude_growth':amplitude_growth}
            df_k=pd.DataFrame(d, columns=['ky_out_board','frequency(kHZ)','B_R(Gauss)amplitude_frequency','B_R(Gauss)amplitude_growth'])
            df_k.to_csv('csv/B_r_inz='+str(inz)+'ky_out_board='+str(ky_GENE_inz[iky])+'.csv',index=False)
    if plot==True:
        imageio.mimwrite('pic/0_B_r_inz='+str(inz)+'ky_dynamic_images.gif', ims_B1)

        plt.clf()
        plt.plot(frequency_kHZ,amplitude_frequency_sum,label='frequency')
        plt.plot(frequency_kHZ,amplitude_growth_sum,label='growth')
        plt.title(r'$B_r(Gauss)$'+' sum over all mode numbers')
        plt.xlabel(r'$f(kHz)$')       
        plt.ylabel(r'$B_r(Gauss)$')
        plt.legend()
        plt.savefig('pic/0Sum_inz='+str(inz)+'B_r_ky_out_board.png') 
        plt.show()

    if output_csv==True:
        d = {'frequency(kHZ)':frequency_kHZ,'B_R(Gauss)amplitude_frequency':amplitude_frequency_sum,'B_R(Gauss)amplitude_growth':amplitude_growth_sum}
        df_sum=pd.DataFrame(d, columns=['frequency(kHZ)','B_R(Gauss)amplitude_frequency','B_R(Gauss)amplitude_growth'])
        df_sum.to_csv('csv/0Sum_B_r_ky_out_board=.csv',index=False)

    return frequency_kHZ,amplitude_frequency_sum,amplitude_growth_sum

    
