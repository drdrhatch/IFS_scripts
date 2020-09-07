#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

Delta_Z=0.07  #7cm as bin for Z 
scan_all_Z=False #Change to True if one want to scan across the whole height
max_Z0=0.21    #Add a small number so it is even
min_Z0=-0.21

max_Z=max_Z0*1.00001    #Add a small number so it is even
min_Z=min_Z0

pic_path='pic'        #path one want to store the picture and video in
csv_path='csv'        #path one want to store the picture and video in


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

#Import the field file using fieldlib
field = fieldfile('field'+suffix,pars)
time = np.array(field.tfld)  #time stampes

#print('shape of field:'+str(np.shape(field)))
#print('shape of real_R:'+str(np.shape(real_R)))
#print('shape of real_Z:'+str(np.shape(real_Z)))


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

time_start,time_end=start_end_time(suffix,pars)
time_start_index=np.argmin(abs(time - time_start))
time_end_index=np.argmin(abs(time - time_end))
time_list = time[time_start_index:time_end_index+1]

if os.path.isdir(csv_path):  #if path does not exist, then create 'csv'
    pass
else:
    os.mkdir(csv_path) 
if os.path.isdir(pic_path):  #if path does not exist, then create 'pic'
    pass
else:
    os.mkdir(pic_path) 
print("**********Scan starts, output in csv and pic***************")


RIP_time_list=np.zeros((len(time_list),len(Z_list)))

for time0 in time_list:
    
    itime = np.argmin(abs(time - time0))
    itime0 = itime
    print("Looking at the spectra at time:"+str(time[itime]))
    #This sets the time step you want to read in
    field.set_time(time[itime])
    
    ntot = field.nx*field.ny*field.nz
    
    phi = np.zeros(ntot,dtype='complex128')
    apar = np.zeros(ntot,dtype='complex128')
    #print "ntot",field.nz*field.nx
    
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

        apar=abs(field.apar()[:,:,:])
        #print("apar"+str(np.shape(apar)))
        apar_ky = np.sum(apar,axis=2)
        (nz0,nky0)=np.shape(apar_ky)
        B1_ky=np.zeros(np.shape(apar_ky))
        #print("apar_ky"+str(np.shape(apar_ky)))

        B1_ky=ky_GENE_grid*apar_ky*Apar_to_B1 #B1 in Gauss


        RIP_list=np.zeros(len(Z_list))
        for nZ_list in range(len(Z_list)):
            RIP_list_temp=0
            for nZ in range(len(real_Z)):
                Z=real_Z[nZ]
                if Z_list[nZ_list]-Delta_Z/2<Z and Z<=Z_list[nZ_list]+Delta_Z/2:
                    #one can restrict the ky range corresprons to frequency 
                    B1=np.sum(B1_ky[nZ,:])
                    RIP_list_temp=RIP_list_temp+B1
            RIP_list[nZ_list]=RIP_list_temp
        
        
        d = {'Z(cm)':Z_list_cm,'B_R(Gauss)':RIP_list}
        df=pd.DataFrame(d, columns=['Z(cm)','B_R(Gauss)'])
        df.to_csv('csv/RIP_t='+str(time[itime])+'.csv',index=False)

    else:  #x_local = False
        pass

#****Read, plot and analyze the data*********
RIP_summary_avg=np.zeros(len(Z_list))
RIP_summary_err=np.zeros(len(Z_list))
Z_err=[Delta_Z/2.]*len(Z_list)
#Z_err_cm=Z_err * 100.
Z_err_cm=[Delta_Z*100./2.]*len(Z_list)

for iZ in range(len(Z_list)):
    RIP_temp_list=[]
    for time0 in time_list:
        itime = np.argmin(abs(time - time0))
        itime0 = itime
        df = pd.read_csv('csv/RIP_t='+str(time[itime])+'.csv')
        #print("df['B_R(Gauss)'][iZ]:"+str(df['B_R(Gauss)'][iZ]))
        RIP_temp_list.append(df['B_R(Gauss)'][iZ])
    RIP_summary_avg[iZ]=np.average(RIP_temp_list)
    RIP_summary_err[iZ]=np.std(RIP_temp_list)
    #print('err:' +str(RIP_summary_err[iZ]))
    #print('avg:' +str(RIP_summary_avg[iZ]))




d = {'Z(cm)':Z_list_cm,'Z_err(cm)':Z_err_cm,'B_R(Gauss)':RIP_summary_avg,'B_R_err(Gauss)':RIP_summary_err}
df_summary=pd.DataFrame(d, columns=['Z(cm)','Z_err(cm)','B_R(Gauss)','B_R_err(Gauss)'])
df_summary.to_csv('csv/0summary_RIP.csv',index=False)


plt.clf()
ax=df_summary.plot(kind='scatter',x='Z(cm)',xerr='Z_err(cm)',y='B_R(Gauss)',yerr='B_R_err(Gauss)',grid=True,label='Average',color='blue')
ax.set_xlabel(r'$Height(cm)$',fontsize=15)
ax.set_ylabel(r'$\bar{B}_r(Gauss)$',fontsize=15)
#plt.xlim(min_Z0*100.,max_Z0*100.)
plt.savefig('pic/0summary_RIP.png')


#*******Makting animation!!!*****
ims_RIP=[]
for time0 in time_list:
    itime = np.argmin(abs(time - time0))
    itime0 = itime
    print("Making animation at time:"+str(time[itime]))
    
    df = pd.read_csv('csv/RIP_t='+str(time[itime])+'.csv')
    RIP_temp_list.append(df['B_R(Gauss)'])
    '''
    plt.clf()
    
    #print(df_summary)
    ax2=df_summary.plot(kind='scatter',x='Z(cm)',xerr='Z_err(cm)',y='B_R(Gauss)',yerr='B_R_err(Gauss)',grid=True,label='Average',color='red',lw=3)
    df.plot.bar(x='Z(cm)',y='B_R(Gauss)',label='t= '+str(time[itime]),color='blue',ax=ax2)
    ax2.set_xlabel(r'$Height(cm)$',fontsize=15)
    ax2.set_ylabel(r'$\bar{B}_r(Gauss)$',fontsize=15)
    plt.title(r'$\bar{B}_r$'+' at t='+str(time[itime]))
    #plt.legend()
    #plt.xlim(min_Z0*100.,max_Z0*100.)
    plt.savefig('pic/RIP_t='+str(time[itime])+'.png')
'''
    plt.clf()
    plt.errorbar(df_summary['Z(cm)'],df_summary['B_R(Gauss)'],xerr=df_summary['Z_err(cm)'],yerr=df_summary['B_R_err(Gauss)'],color='red',label='Average')
    plt.bar(df['Z(cm)'], df['B_R(Gauss)'], color ='blue', width = 4,label='t= '+str(time[itime])) 
    plt.xlabel('Height(m)',fontsize=10)
    plt.ylabel(r'$\bar{B}_r(Gauss)$',fontsize=10)
    plt.title(r'$\bar{B}_r$'+' at t='+str(time[itime]))
    plt.legend()
    plt.xlim(min_Z0*100.,max_Z0*100.)
    plt.savefig('pic/RIP_t='+str(time[itime])+'.png')

    file_name='pic/RIP_t='+str(time[itime])+'.png'
    ims_RIP.append(imageio.imread(file_name))


imageio.mimwrite('pic/0LN_RIP_dynamic_images.gif', ims_RIP)

print("Finished, Everything is in the directory csv and pic")

