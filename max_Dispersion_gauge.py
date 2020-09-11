# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 16:56:45 2020

@author: jlara, maxcurie
"""
import numpy as np
import pandas as pd
from finite_differences import *
import matplotlib.pyplot as plt
from interp import *
import math
import csv
#The following function computes the growth rate of the slab MTM. The input are the physical parameters and it outputs the growth rate in units of omega_{*n}
#Parameters for function:
#nu: normalized to omega_{*n} 
#shat=L_n/L_s
#beta=plasma beta
#eta=L_n/L_T
#ky normalized to rho_i (ion gyroradius)

#Late Edited by Max Curie, 06/21/2020
#Way to run: python ~/max/scripts/MTMDispersion.py
#at x=0.98032 mid pedestal

n0=21. #toroidal mode number

(Re_1_min, Re_1_max)=(-10.,10.)
(Im_1_min, Im_1_max)=(-10.,10.)
(Re_2_min, Re_2_max)=(-10.,10.)
(Im_2_min, Im_2_max)=(-10.,10.)
(Re_3_min, Re_3_max)=(-10.,10.)
(Im_3_min, Im_3_max)=(-10.,10.)
alpha_scan_resolution=1.

plot=True 
csv_output=True

nu0=36.1471134158304
eta=2.14865959330967
shat=0.277548602587741
beta=0.000711463834687449
ky0=0.00772097143444933

factor_nu=0.317739083372979

ky=ky0*n0
nu=nu0/n0
len_list=50  #length of the list 
nu_max=0.5
nu_list=np.arange(0., nu_max/factor_nu, (nu_max/(factor_nu*float(len_list))))

eta_list=[eta] * len_list
shat_list=[shat] * len_list
beta_list=[beta] * len_list
ky_list=[ky] * len_list

nu_list_omega_star=nu_list*factor_nu



def Dispersion(nu,eta,shat,beta,ky,alpha1,alpha2,alpha3):
    factor=-alpha1*((ky*shat)/( alpha2*beta-eta*(eta+1.)/1836.))*(1/np.sqrt(1836.))+alpha3
    a=complex(0.,-15753.1)-32045.*factor
    b=complex(0.,15753.1)+complex(0,16426.6)*eta+11951.2*nu-complex(0,40594.8)*factor*nu
    c=-11951.2*nu-20267.7*nu*eta+9121.48*factor*nu**2
    gamma=(-b+np.sqrt(b**2-4.*a*c))/(2.*a)
    return gamma,factor
def Dispersion_list(nu,eta,shat,beta,ky,factor_nu,alpha1,alpha2,alpha3,n,plot):
    nx=len(nu)
    #print(nx)
    gamma_complex_temp=0
    factor_temp=0
    gamma_complex=[]
    factor=[]

    for i in range(nx):
        gamma_complex_temp,factor_temp=Dispersion(nu[i],eta[i],shat[i],beta[i],ky[i],alpha1,alpha2,alpha3)
        gamma_complex.append(gamma_complex_temp)
        factor.append(factor_temp)
    
    gamma_complex=np.asarray(gamma_complex)
    factor=np.asarray(factor)
    gamma=gamma_complex.real
    omega=gamma_complex.imag

    if plot==True:
        plt.clf()
        plt.xlabel('nu/omega*')
        plt.ylabel('gamma/omega*n') 
        plt.plot(nu*factor_nu,gamma)
        plt.savefig('Pic/'+str(n)+'.png')
    return gamma,omega,factor

def Dispersion_change(nu_list,eta_list,shat_list,beta_list,ky_list,factor_nu,alpha_scan_resolution,plot,csv_output):
    #Fit Parameters
    #def dispersion_gauge(iterdb_file_name,geomfile,scan_n0,plot_profile,plot_peak_scan,csv_profile,csv_peak_scan): 
    Re_1=np.arange(Re_1_min, Re_1_max, alpha_scan_resolution)
    Im_1=np.arange(Im_1_min, Im_1_max, alpha_scan_resolution)
    Re_2=np.arange(Re_2_min, Re_2_max, alpha_scan_resolution)
    Im_2=np.arange(Im_2_min, Im_2_max, alpha_scan_resolution)
    Re_3=np.arange(Re_3_min, Re_3_max, alpha_scan_resolution)
    Im_3=np.arange(Im_3_min, Im_3_max, alpha_scan_resolution)
    len_scan=len(Re_1)

    Re_1_list=[]
    Im_1_list=[]
    Re_2_list=[]
    Im_2_list=[]
    Re_3_list=[]
    Im_3_list=[]
    count_list=[]
    nu_peak=[]  #nu over omega*
    gamma_peak=[]

    with open('csv/0summery.csv', 'w') as csvfile:
        data = csv.writer(csvfile, delimiter=',')
        data.writerow(['n_count','nu_ei_omega','gamma_omega_n','Re_1','Im_1','Re_2','Im_2','Re_3','Im_3'])
    n=0
    for Ri_1 in range(len_scan):
        for Ii_1 in range(len_scan):
            for Ri_2 in range(len_scan):
                for Ii_2 in range(len_scan):
                    for Ri_3 in range(len_scan):
                        for Ii_3 in range(len_scan):
                            alpha1=complex(Re_1[Ri_1],Im_1[Ii_1])
                            alpha2=complex(Re_2[Ri_2],Im_2[Ii_2])
                            alpha3=complex(Re_3[Ri_3],Im_3[Ii_3])
                            Re_1_list.append(Re_1[Ri_1])
                            Im_1_list.append(Im_1[Ii_1])
                            Re_2_list.append(Re_2[Ri_2])
                            Im_2_list.append(Im_2[Ii_2])
                            Re_3_list.append(Re_3[Ri_3])
                            Im_3_list.append(Im_3[Ii_3])
                            gamma,omega,factor=Dispersion_list(nu_list,eta_list,shat_list,beta_list,ky_list,factor_nu,alpha1,alpha2,alpha3,n,plot)
                            index=np.argmax(gamma)
                            count_list.append(n)
                            nu_peak.append(nu_list[index]*factor_nu)
                            gamma_peak.append(gamma[index])

                            if csv_output==True:
                                d = {'nu_ei_omega':nu_list_omega_star,'gamma_omega_n':gamma,'nu_omegan':nu_list,'eta':eta_list,'shat':shat_list,'beta':beta_list,'ky':ky_list}
                                df=pd.DataFrame(d, columns=['nu_ei_omega','gamma_omega_n','nu_omegan','eta','shat','beta','ky'])
                                df.to_csv('csv/'+str(n)+'.csv',index=False)

                            with open('csv/0summery.csv', 'a') as csvfile:
                                data = csv.writer(csvfile, delimiter=',')
                                data.writerow([n,nu_list[index]*factor_nu,gamma[index],Re_1[Ri_1],Im_1[Ii_1],Re_2[Ri_2],Im_2[Ii_2],Re_3[Ri_3],Im_3[Ii_3]])

                            n=n+1

    if csv_output==True:
        d1 = {'n_count':count_list,'nu_ei_omega':nu_peak,'gamma_omega_n':gamma_peak,'Re_1':Re_1_list,'Im_1':Im_1_list,'Re_2':Re_2_list,'Im_2':Im_2_list,'Re_3':Re_3_list,'Im_3':Im_3_list}
        df1=pd.DataFrame(d1, columns=['n_count','nu_ei_omega','gamma_omega_n','Re_1','Im_1','Re_2','Im_2','Re_3','Im_3'])
        df1.to_csv('csv/summary.csv',index=False)
    
    #return gamma,factor
Dispersion_change(nu_list,eta_list,shat_list,beta_list,ky_list,factor_nu,alpha_scan_resolution,plot,csv_output)

