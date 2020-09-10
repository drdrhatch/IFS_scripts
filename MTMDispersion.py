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
from read_EFIT import *
from read_EFIT_file import *
from read_iterdb_file import *
from max_pedestal_finder import find_pedestal
from max_stat_tool import *
#The following function computes the growth rate of the slab MTM. The input are the physical parameters and it outputs the growth rate in units of omega_{*n}
#Parameters for function:
#nu: normalized to omega_{*n} 
#shat=L_n/L_s
#beta=plasma beta
#eta=L_n/L_T
#ky normalized to rho_i (ion gyroradius)

#Late Edited by Max Curie, 06/21/2020
#Way to run: python ~/max/scripts/MTMDispersion.py


#**************Block for user******************************************
#**************Setting up*********************************************

#iterdb_file_name='DIIID175823.iterdb'  #name of the iterdb file
#geomfile='g175823.04108_257x257'           #name of the magnetic geometry file
iterdb_file_name='DIIID174864.iterdb'  #name of the iterdb file
geomfile='g174864'           #name of the magnetic geometry file

#iterdb_file_name = 'efit_Dial_Nmod5_Zp2_48_new.iterdb' #name of the iterdb file
#iterdb_file_name = 'efit_Dial_N_Zp2_48_new.iterdb' #name of the iterdb file
#geomfile = 'efit_Dial_N_Zp2_48'                     #name of the magnetic geometry file
#geomfile = 'efit_temp'                     #name of the magnetic geometry file

run_mode_finder=True        #Change to True if one want to run mode finder 
run_nu_scan=False           #Change to True if one want to run collisionality scan 


omega_percent=10.                       #choose the omega within the top that percent defined in(0,100)
n_min=1                                #minmum mode number (include) that finder will cover
n_max=30                               #maximum mode number (include) that finder will cover
bins=800                               #sizes of bins to smooth the function
plot_profile=True                    #Set to True is user want to have the plot of the profile
plot_n_scan=False                      #Set to True is user want to have the plot of the gamma over n
csv_profile=False                    #Set to True is user want to have the csv file "profile_output.csv" of the profile
csv_n_scan=True                       #Set to True is user want to have the csv file "MTM_dispersion_n_scan.csv" of the gamma over n
plot_spectrogram=True
peak_of_plasma_frame=False             #Set to True if one want to look around the peak of omega*e in plasam frame
manual_ped=1
mid_ped0=0.958

#******For scaning********
scan_n0=3.
choose_location=True    #Change to True if one wants to change the location manually 
location=0.984139203080616
plot_peak_scan=True
csv_peak_scan=True
nu_percent=10  #about the nu0 for x% 1=100%
#**************End of Setting up*********************************************
#**************End of Block for user******************************************

def Dispersion_V1(nu,eta,shat,beta,ky):
  #Fit Parameters
  alpha1=-10.-1.j
  alpha2=-10.-4.j
  alpha3=-7. - 3.j
  factor=-alpha1*((ky*shat)/( alpha2*beta-eta*(eta+1.)/1836.))*(1/np.sqrt(1836))+alpha3
  a=complex(0,-15753.1)-32045.*factor
  b=complex(0,15753.1)+complex(0,16426.6)*eta+11951.2*nu-complex(0,40594.8)*factor*nu
  c=-11951.2*nu-20267.7*nu*eta+9121.48*factor*nu**2
  gamma=(-b+np.sqrt(b**2-4*a*c))/(2.*a)
  return gamma,factor


#Version 2.0
def Dispersion(nu,eta,shat,beta,ky):
    #Fit Parameters
    factor=-((ky*shat/beta))*(1/np.sqrt(1836))
    a=complex(0,-15753.1)-32045.*factor
    b=complex(0,15753.1)+complex(0,16426.6)*eta+11951.2*nu-complex(0,40594.8)*factor*nu
    c=-11951.2*nu-20267.7*nu*eta+9121.48*factor*nu**2
    gamma=(-b+np.sqrt(b**2-4*a*c))/(2.*a)
    return gamma,factor
 

#return nu,ky for the case n_tor=1 for the given location(default to be pedestal)
def Parameter_reader(iterdb_file_name,geomfile,manual_ped,mid_ped0,plot,output_csv):
    n0=1.
    mref = 2.        # mass of ion in proton mass
    rhot0, te0, ti0, ne0, ni0, nz0, vrot0 = read_iterdb_file(iterdb_file_name)
    EFITdict = read_EFIT(geomfile)
    xgrid = EFITdict['psipn']
    q = EFITdict['qpsi']

    uni_rhot = np.linspace(min(rhot0),max(rhot0),len(rhot0)*10.)

    te_u = interp(rhot0,te0,uni_rhot)
    ne_u = interp(rhot0,ne0,uni_rhot)
    ni_u = interp(rhot0,ni0,uni_rhot)
    vrot_u = interp(rhot0,vrot0,uni_rhot)
    q      = interp(xgrid,q,uni_rhot)
    tprime_e = -fd_d1_o4(te_u,uni_rhot)/te_u
    nprime_e = -fd_d1_o4(ne_u,uni_rhot)/ne_u
    qprime = fd_d1_o4(q,uni_rhot)/q


    #center_index = np.argmax((tprime_e*te_u+nprime_e*ne_u)[0:int(len(tprime_e)*0.99)])
    
    midped, topped=find_pedestal(file_name=geomfile, path_name='', plot=False)
    x0_center = midped
    if manual_ped==1:
        x0_center=mid_ped0

    print('mid pedestal is at r/a = '+str(x0_center))
    
    Lref, Bref, R_major, q0, shat0=get_geom_pars(geomfile,x0_center)

    index_begin=np.argmin(abs(uni_rhot-topped+1.*(1.-x0_center)))

    te_u = te_u[index_begin:len(uni_rhot)-1]
    ne_u = ne_u[index_begin:len(uni_rhot)-1]
    ni_u = ni_u[index_begin:len(uni_rhot)-1]
    vrot_u = vrot_u[index_begin:len(uni_rhot)-1]
    q      = q[index_begin:len(uni_rhot)-1]
    tprime_e = tprime_e[index_begin:len(uni_rhot)-1]
    nprime_e = nprime_e[index_begin:len(uni_rhot)-1]
    qprime   = qprime[index_begin:len(uni_rhot)-1]
    uni_rhot = uni_rhot[index_begin:len(uni_rhot)-1]

    Lt=1./tprime_e
    Ln=1./nprime_e
    Lq=1./qprime
    
    center_index = np.argmin(abs(uni_rhot-x0_center))

    ne=ne_u/(10**19)      # in 10^19 /m^3
    ni=ni_u/(10**19)      # in 10^19 /m^3
    te=te_u/1000.          #in keV
    m_SI = mref *1.6726*10**(-27)
    me_SI = 9.11*10**(-31)
    c  = 1.
    qref = 1.6*10**(-19)
    #refes to GENE manual
    coll_c=2.3031*10**(-5)*Lref*ne/(te)**2*(24-np.log(np.sqrt(ne*10**13)/(te*1000)))
    coll_ei=4*(ni/ne)*coll_c*np.sqrt(te*1000.*qref/me_SI)/Lref
    nuei=coll_ei
    beta=403.*10**(-5)*ne*te/Bref**2.

    nref=ne_u[center_index]
    te_mid=te_u[center_index]
    Tref=te_u[center_index] * qref

    cref = np.sqrt(Tref / m_SI)
    Omegaref = qref * Bref / m_SI / c
    rhoref = cref / Omegaref 
    kymin=n0*q0*rhoref/(Lref*x0_center)
    kyGENE =kymin * (q/q0) * np.sqrt(te_u/te_mid) * (x0_center/uni_rhot) #Add the effect of the q varying
    #from mtm_doppler
    omMTM = kyGENE*(tprime_e+nprime_e)
    gyroFreq = 9.79E3/np.sqrt(mref)*np.sqrt(te_u)/Lref
    mtmFreq = omMTM*gyroFreq/(2.*np.pi*1000.)
    omegaDoppler = vrot_u*n0/(2.*np.pi*1000.)
    omega=mtmFreq + omegaDoppler

    omega_n_GENE=kyGENE*(nprime_e)       #in cs/a
    omega_n=omega_n_GENE*gyroFreq/(2.*np.pi*1000.)  #in kHz

    coll_ei=coll_ei/(1000.)  #in kHz

    shat=Ln/Lq
    eta=Ln/Lt
    ky=kyGENE
    nu=(coll_ei)/(2.*np.pi*omega_n) 


    if plot==True:
        plt.clf()
        plt.xlabel('r/a')
        plt.ylabel('eta') 
        plt.plot(uni_rhot,eta,label='eta')
        plt.show()

        plt.clf()
        #plt.title('mode number finder')
        plt.xlabel('r/a')
        plt.ylabel('omega*(Lab), kHz') 
        plt.plot(uni_rhot,omega,label='omega*(Lab)')
        plt.show()

        plt.clf()
        plt.xlabel('r/a')
        plt.ylabel('eta') 
        plt.plot(uni_rhot,eta,label='eta')
        plt.show()

        plt.clf()
        plt.xlabel('r/a')
        plt.ylabel('shat') 
        plt.plot(uni_rhot,shat,label='shat')
        plt.show()


        plt.clf()
        plt.xlabel('r/a')
        plt.ylabel('beta') 
        plt.plot(uni_rhot,beta,label='beta')
        plt.show()
        
        plt.clf()
        plt.xlabel('r/a')
        plt.ylabel('ky rhoi') 
        plt.plot(uni_rhot,ky,label='ky')
        plt.show()

    
    if output_csv==True:
        with open('profile_output.csv','w') as csvfile:
            data = csv.writer(csvfile, delimiter=',')
            data.writerow(['x/a','nu_ei(kHz)','omega*n(kHz)','omega* plasma(kHz)','Doppler shift(kHz)','nu/omega*n','eta','shat','beta','ky rhoi(for n=1)'])
            for i in range(len(uni_rhot)):
                data.writerow([uni_rhot[i],coll_ei[i],omega_n[i],mtmFreq[i],omegaDoppler[i],nu[i],eta[i],shat[i],beta[i],ky[i]])
        csvfile.close()
    
    return uni_rhot,nu,eta,shat,beta,ky,q,mtmFreq,omegaDoppler,omega_n,omega_n_GENE

#scan the Dispersion for the given location(default to be pedestal)
def Dispersion_list(uni_rhot,nu,eta,shat,beta,ky,plot):
    nx=len(uni_rhot)
    #print(nx)
    gamma_complex_temp=0
    factor_temp=0
    gamma_complex=[]
    factor=[]

    for i in range(nx):
        gamma_complex_temp,factor_temp=Dispersion(nu[i],eta[i],shat[i],beta[i],ky[i])
        gamma_complex.append(gamma_complex_temp)
        factor.append(factor_temp)
    
    gamma_complex=np.asarray(gamma_complex)
    factor=np.asarray(factor)
    gamma=gamma_complex.imag
    omega=gamma_complex.real

    if plot==True:
        plt.clf()
        plt.xlabel('r/a')
        plt.ylabel('gamma/omega*n') 
        plt.plot(uni_rhot,gamma,label='gamma')
        plt.show()
    return gamma,omega,factor

#this function takes the q and n, returns the locations of the rational surfaces
def Rational_surface(uni_rhot,q,n0):
    x_list=[]
    m_list=[]
    qmin = np.min(q)
    qmax = np.max(q)

    #print("q_min is: "+str(qmin))
    #print("q_max is: "+str(qmax))
    m_min = math.ceil(qmin*n0)
    m_max = math.floor(qmax*n0)
    mnums = np.arange(m_min,m_max+1)

    print(m_min,m_max)
    #or m in range(m_min,m_max+1):
    for m in mnums:
    	#print(m)
        q0=float(m)/float(n0)
        index0=np.argmin(abs(q-q0))
        if abs(q[index0]-q0)<0.1:
            x_list.append(uni_rhot[index0])
            m_list.append(m)

    #print(x_list)
    #x_list=np.asarray(x_list)
    #m_list=np.asarray(m_list)
    return x_list, m_list

#this function finds the a peak of the 
def Peak_of_drive(uni_rhot,mtmFreq,omegaDoppler,omega_percent):
    x_peak_range=[]
    x_range_ind=[]
    #omega=mtmFreq+omegaDoppler
    if peak_of_plasma_frame==True:
        omega=mtmFreq
    else:
        omega=mtmFreq+omegaDoppler
    #if manual_ped==1:
    #    mid_ped0
    #mid_ped0=0.958
    omega_max=np.max(omega)
    print(mtmFreq)
    print(omega_max)
    omega_min=omega_max*(100.-omega_percent)/100.
    for i in range(len(uni_rhot)):
        if omega[i] >= omega_min:
            x_peak_range.append(uni_rhot[i])
            x_range_ind.append(i)
    print("x="+str(min(x_peak_range))+"~"+str(max(x_peak_range)))
    return x_peak_range, x_range_ind

#scan toroidial mode number
def Dispersion_n_scan(uni_rhot,nu,eta,shat,beta,ky,q,omega_n,omega_n_GENE,mtmFreq,omegaDoppler,x_peak_range,x_range_ind,n_min,n_max,plot,output_csv):
    ind_min  =min(x_range_ind)
    ind_max  =max(x_range_ind)
    uni_rhot_full=uni_rhot
    nu_full=nu
    eta_full=eta
    shat_full=shat
    beta_full=beta
    ky_full=ky
    q_full=q
    omega_n_full=omega_n
    omega_n_GENE_full=omega_n_GENE
    omegaDoppler_full=omegaDoppler
    mtmFreq_full=mtmFreq

    uni_rhot_top=uni_rhot[ind_min:ind_max]
    print("uni_rhot="+str(min(uni_rhot_top))+"~"+str(max(uni_rhot_top)))
    nu_top=nu[ind_min:ind_max]
    eta_top=eta[ind_min:ind_max]
    shat_top=shat[ind_min:ind_max]
    beta_top=beta[ind_min:ind_max]
    ky_top=ky[ind_min:ind_max]
    q_top=q[ind_min:ind_max]
    omega_n_top=omega_n[ind_min:ind_max]
    omega_n_GENE_top=omega_n_GENE[ind_min:ind_max]
    omegaDoppler_top=omegaDoppler[ind_min:ind_max]
    mtmFreq_top=mtmFreq[ind_min:ind_max]
    omega_star_max=max(mtmFreq_top)

    n_list=[]
    m_list=[]
    x_list=[]
    gamma_list=[]
    omega_list=[]
    factor_list=[]
    gamma_list_kHz=[]
    omega_list_kHz=[]
    omega_list_Lab_kHz=[]
    omega_star_list_kHz=[]
    omega_star_list_Lab_kHz=[]
    nu_ei_omega_plasma_list=[]
    omega_omega_peak_list=[]
    kHz_to_cs_a_list=[]
    
    if plot==True:
        plt.clf()
        plt.title('mode number scan')
        plt.xlabel('r/a')
        plt.ylabel('gamma/omega*n') 
        
        

    for n0 in range(n_min,n_max+1):
        print("************n="+str(n0)+"************")

        gamma,omega,factor=Dispersion_list(uni_rhot_full,nu_full/float(n0),eta_full,shat_full,beta_full,ky_full*float(n0),plot=False)
        x0_list, m0_list=Rational_surface(uni_rhot_top,q_top,n0)
        #print(x_list)
        if plot==True and max(gamma)>0:
            plt.plot(uni_rhot,gamma)   #,label='n='+str(n0))

        for i in range(len(x0_list)):
            x=x0_list[i]
            #print(x)
            m=m0_list[i]
            #print(m)

            x_index=np.argmin(abs(x-uni_rhot_top))
            gamma,factor=Dispersion(nu_top[x_index]/float(n0),eta_top[x_index],shat_top[x_index],beta_top[x_index],ky_top[x_index]*float(n0))
            
            gamma_complex=gamma
            gamma=gamma_complex.imag
            omega=gamma_complex.real

            
            kHz_to_cs_a=omega_n_GENE_top[x_index]/omega_n_top[x_index]   # kHz * (kHz_to_cs_a)= cs/a unit
            omega_n_temp=omega_n_top[x_index]*float(n0)
            gamma_kHz=gamma*omega_n_temp
            omega_kHz=omega*omega_n_temp
            omega_Lab_kHz=omega_kHz-omegaDoppler_top[x_index]*float(n0)
            omega_star_kHz=mtmFreq_top[x_index]*float(n0)
            omega_star_Lab_kHz=mtmFreq_top[x_index]*float(n0)+omegaDoppler_top[x_index]*float(n0)
            nu_ei_omega_plasma=nu_top[x_index]/float(n0)*omega_n_top[x_index]*float(n0)/(2.*np.pi*omega_star_kHz)

            gamma_list_kHz.append(gamma_kHz)
            kHz_to_cs_a_list.append(kHz_to_cs_a)
            omega_list_kHz.append(omega_kHz)
            omega_list_Lab_kHz.append(omega_Lab_kHz)
            omega_star_list_kHz.append(omega_star_kHz)
            omega_star_list_Lab_kHz.append(omega_star_Lab_kHz)
            nu_ei_omega_plasma_list.append(nu_ei_omega_plasma)
            gamma_list.append(gamma)
            omega_list.append(omega)
            factor_list.append(factor)
            n_list.append(n0)
            m_list.append(m)
            x_list.append(x)
            omega_omega_peak_list.append(omega_star_kHz/omega_star_max/float(n0))

            #print("x="+str(x)+", gamma(kHz)="+str(gamma_kHz))
            if plot==True and gamma>0:
                plt.axvline(x,color='red',alpha=0.05)

    if plot==True:
        plt.axvline(min(uni_rhot_top),color='green',label="near the peak of omega*")
        plt.axvline(max(uni_rhot_top),color='green',label="near the peak of omega*")
        plt.axhline(0,color='red',label="gamma=0")
        plt.legend()
        plt.show()
    
    #print(len(x_list),len(n_list),len(m_list))
    if output_csv==True:
        with open('MTM_dispersion_n_scan.csv','w') as csvfile:
            data = csv.writer(csvfile, delimiter=',')
            data.writerow([   'x/a',     'n',      'm',      'gamma(kHz)',      'gamma(cs/a)',    'nu_ei/omega*(plasma)','omega*/omega*_max' ,'omega_plasma(kHz)','omega_lab(kHz)',   'omega_star_plasma(kHz)','omega_star_lab(kHz)',    'factor'])
            for i in range(len(x_list)):
                data.writerow([x_list[i],n_list[i],m_list[i],gamma_list_kHz[i],gamma_list_kHz[i]*kHz_to_cs_a_list[i],nu_ei_omega_plasma_list[i],omega_omega_peak_list[i],omega_list_kHz[i],omega_list_Lab_kHz[i],omega_star_list_kHz[i],omega_star_list_Lab_kHz[i],factor_list[i]])
        csvfile.close()

    return x_list,n_list,m_list,gamma_list,omega_list,factor_list,gamma_list_kHz,omega_list_kHz,omega_list_Lab_kHz,omega_star_list_kHz,omega_star_list_Lab_kHz

#define a normal distribusion, input the x axis as x_list, output norm(x_list)
def normal(x_list,x0,mu,sigma):
    #return x0*1./(sigma*sqrt(2.*np.pi))*exp(-1./2.*((x-mu)/sigma)**2.)
    return x0*exp(-1./2.*((x_list-mu)/sigma)**2.)

#Calculate the gamma as function of frequency
#input the omega*n for n=1  
def Spectrogram(gamma_list_kHz,omega_list_kHz):

    omega_min=min(omega_list_kHz)
    omega_max=max(omega_list_kHz)
    #print(omega_min)
    #print(omega_max)
    f=np.arange(omega_min,omega_max,0.01)
    gamma_f=np.zeros(len(np.arange(omega_min,omega_max,0.01)))
    #print(gamma_f)
    for i in range(len(gamma_list_kHz)):
        if gamma_list_kHz[i]>0:
            x0=gamma_list_kHz[i]
        else:
            x0=0
        mu=omega_list_kHz[i]
        sigma=0.2
        gamma_f=gamma_f+normal(f,x0,mu,sigma)
    return f,gamma_f

def Spectrogram_2_frames(gamma_list_kHz,omega_list_kHz,omega_list_Lab_kHz,bins,plot):
    
    f_lab,gamma_f_lab=Spectrogram(gamma_list_kHz,omega_list_Lab_kHz)
    #f_lab_pd_frame=pd.DataFrame(f_lab, columns=['a']) 
    #f_lab_smooth=f_lab_pd_frame.diff().rolling(window=dt).mean()
    #gamma_f_lab_smooth=gamma_f_lab.diff().rolling(window=dt).mean()
    if plot==True:
        plt.clf()
        plt.title('Spectrogram in lab frame')
        plt.xlabel('f(kHz)')
        plt.ylabel('gamma(a.u.)') 
        plt.plot(f_lab,gamma_f_lab)
        plt.show()

        plt.clf()
        plt.title('Spectrogram in lab frame(rolling average)')
        plt.xlabel('f(kHz)')
        plt.ylabel('gamma(a.u.)') 
        plt.plot(smooth(f_lab,bins)[0],smooth(gamma_f_lab,bins)[0])
        plt.show()


    f_plasma,gamma_f_plasma=Spectrogram(gamma_list_kHz,omega_list_kHz)
    if plot==True:
        plt.clf()
        plt.title('Spectrogram in plasma frame')
        plt.xlabel('f(kHz)')
        plt.ylabel('gamma(a.u.)') 
        plt.plot(f_plasma,gamma_f_plasma)
        plt.show()

        plt.clf()
        plt.title('Spectrogram in plasma frame(rolling average)')
        plt.xlabel('f(kHz)')
        plt.ylabel('gamma(a.u.)') 
        plt.plot(smooth(f_plasma,bins)[0],smooth(gamma_f_plasma,bins)[0])
        plt.show()
    return f_lab,gamma_f_lab,f_plasma,gamma_f_plasma

def peak_scan(uni_rhot,nu,eta,shat,beta,ky,q,omega_n,omega_n_GENE,mtmFreq,omegaDoppler,n0,plot_peak_scan,csv_peak_scan):
    center_index=np.argmax(mtmFreq)
    x0=uni_rhot[center_index]
    if choose_location==True:
        x0=location
        center_index=np.argmin(abs(uni_rhot-x0))
    #x0=0.98
    #center_index=np.argmin(abs(uni_rhot-x0))
    #print(x0)
    nu0=nu[center_index]/float(n0)
    eta0=eta[center_index]
    shat0=shat[center_index]
    beta0=beta[center_index]
    ky0=ky[center_index]*float(n0)

    omega_n_0=omega_n[center_index]*float(n0)
    omega_n_GENE_0=omega_n_GENE[center_index]*float(n0)
    mtmFreq0=mtmFreq[center_index]*float(n0)
    omegaDoppler0=omegaDoppler[center_index]*float(n0)
    
    uni_rhot_list=[]
    nu_list=[]
    eta_list=[]
    shat_list=[]
    beta_list=[]
    ky_list=[]
    nu_list_omega_star=[]

    #print(uni_rhot[center_index])
    nmax=400
    

    nu_min=nu0*(1.-nu_percent/2.)
    d_nu=nu0*nu_percent/float(nmax)

    for i in range(nmax):
        nu_list.append(nu_min+float(i)*d_nu)
        uni_rhot_list.append(x0)
        eta_list.append(eta0)
        shat_list.append(shat0)
        beta_list.append(beta0)
        ky_list.append(ky0)
        nu_list_omega_star.append((nu_min+float(i)*d_nu)*omega_n_0/(mtmFreq0))
    
    gamma,omega,factor=Dispersion_list(uni_rhot_list,nu_list,eta_list,shat_list,beta_list,ky_list,plot=False)
    gamma_cs_a=gamma*omega_n_GENE_0
    
    if csv_peak_scan==True:
        d = {'nu_omegan':nu_list,'eta':eta_list,'shat':shat_list,'beta':beta_list,'ky':ky_list,'nu_ei_omega_plasma':nu_list_omega_star,'gamma_cs_a':gamma_cs_a}
        df=pd.DataFrame(d, columns=['nu_omegan','eta','shat','beta','ky','nu_ei_omega_plasma','gamma_cs_a'])
        df.to_csv('nu_scan.csv',index=False)
    #print(nu_list)
    #nu_list_omega_star=nu_list*omega_n_0/(mtmFreq0+omegaDoppler0)

    if plot_peak_scan==True:
        plt.clf()
        plt.title('Nu scan')
        plt.xlabel('nu/omega*')
        plt.ylabel('gamma(cs/a)') 
        plt.plot(nu_list_omega_star,gamma_cs_a)
        plt.show()


def MTM_scan(iterdb_file_name,geomfile,omega_percent,bins,n_min,n_max,plot_profile,plot_n_scan,csv_profile,csv_n_scan): 
    #return nu,ky for the case n_tor=1 for the given location(default to be pedestal)
    uni_rhot,nu,eta,shat,beta,ky,q,mtmFreq,omegaDoppler,omega_n,omega_n_GENE=Parameter_reader(iterdb_file_name,geomfile,manual_ped,mid_ped0,plot=plot_profile,output_csv=csv_profile)
    print(mtmFreq)
    print(omegaDoppler)
    x_peak_range, x_range_ind=Peak_of_drive(uni_rhot,mtmFreq,omegaDoppler,omega_percent)

    #scan toroidial mode number
    x_list,n_list,m_list,gamma_list,\
    omega_list,factor_list,gamma_list_kHz,\
    omega_list_kHz,omega_list_Lab_kHz,\
    omega_star_list_kHz,omega_star_list_Lab_kHz\
    =Dispersion_n_scan(uni_rhot,nu,eta,shat,beta,ky,q,\
    omega_n,omega_n_GENE,mtmFreq,omegaDoppler,x_peak_range,x_range_ind,\
    n_min,n_max,plot=plot_n_scan,output_csv=csv_n_scan)

    f_lab,gamma_f_lab,f_plasma,gamma_f_plasma=Spectrogram_2_frames(gamma_list_kHz,omega_star_list_kHz,omega_star_list_Lab_kHz,bins,plot=plot_spectrogram)
    return x_list,n_list,m_list,gamma_list,omega_list,factor_list


def coll_scan(iterdb_file_name,geomfile,n0,plot_profile,plot_peak_scan,csv_profile,csv_peak_scan): 
    uni_rhot,nu,eta,shat,beta,ky,q,mtmFreq,omegaDoppler,omega_n,omega_n_GENE=Parameter_reader(iterdb_file_name,geomfile,manual_ped,mid_ped0,plot=plot_profile,output_csv=csv_profile)
    peak_scan(uni_rhot,nu,eta,shat,beta,ky,q,omega_n,omega_n_GENE,mtmFreq,omegaDoppler,n0,plot_peak_scan,csv_peak_scan)

if run_mode_finder==True:
    x_list,n_list,m_list,gamma_list,omega_list,factor_list=MTM_scan(iterdb_file_name,geomfile,omega_percent,bins,n_min,n_max,plot_profile,plot_n_scan,csv_profile,csv_n_scan)
if run_nu_scan==True:
    coll_scan(iterdb_file_name,geomfile,scan_n0,plot_profile,plot_peak_scan,csv_profile,csv_peak_scan)
