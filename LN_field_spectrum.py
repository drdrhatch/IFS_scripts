#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os
import matplotlib.pyplot as plt
import imageio   #for making animation
from fieldlib import fieldfile
from geomWrapper import ky
from geomWrapper import init_read_geometry_file
from read_write_geometry import read_geometry_local
from ParIO import Parameters 
import optparse as op
import sys



#from read_write_geometry import *

time_scan=True    #Change to True if one want to scan trough time. 
#time_start
#time_stop
path='pic'        #path one want to store the picture and video in


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
nky0 = pars['nky0']       #in rhoi
n_step = pars['n0_global']       #in rhoi
nref = nref * 1.E19         #in the unit of /m^3
Tref = Tref * qref * 1.E03  #in the unit of J
mref = mref * m_kg          #in the unit of kg
pref = nref * Tref          #in Pa*k_{B}
cref = np.sqrt(Tref / mref) #in the unit of m/s
Omegaref = qref * Bref / mref / c  #in rad/s
rhoref = cref / Omegaref           #in m rho_i(ion gyroradii)
rhorefStar = rhoref / Lref         #Unitless


#ky comes from geomWrapper.py
ky_GENE_temp=ky(pars, geom_coeff,plot=False)     #ky for ky min
#print('ky shape: '+str(np.shape(ky_GENE_temp)))

#Determine the n for kymin
ky_n1=1.*q0*rhoref/(Lref*x0)  #ky for n=1
n_min=round(kymin/ky_n1)   #n for ky_min
if n_min==0:
    n_min=1
n_list=np.arange(int(n_min),int(n_min+nky0*n_step),int(n_step))
ky_GENE_n1=ky_GENE_temp/float(n_min)
kygrid = n_list*ky_n1

ky_GENE_ob_n1=ky_GENE_n1[int(pars['nz0']/2)]
kygrid_ob=ky_GENE_ob_n1*n_list

print('n0 list length: '+str(len(n_list)))
print('n0 list: '+str(n_list))




#B1=abs(np.mean(Apar_GENE[z,:])*len(Apar_GENE[z,:])*(ky_GENE_temp[z]/rhoref)*Bref*B_gauss*rhorefStar*rhoref)
Apar_to_B1=abs((1./rhoref)*Bref*B_gauss*rhorefStar*rhoref)         #B1=Apar*ky_GENE_temp*Apar_to_B1


if time_scan==True:
    time_list = np.array(field.tfld)
    if os.path.isdir(path):  #if path does not exist, then create 'pic'
        pass
    else:
        print("*************************")
        os.mkdir(path) 
    #for test 
    #time_list=time_list[int(99*len(time_list)/100):-1]
else:
    time_list = [-1]       #Just look at the last second


for time0 in time_list:
    if time0 == -1:
        itime = -1
        itime0 = len(time)-1
    else: 
        itime = np.argmin(abs(time - time0))
        itime0 = itime
    print("Looking at the spectra at time:"+str(time[itime])+'ms')
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
        #print("zgrid"+str(zgrid))
        phi=field.phi()[:,:,:]
        print('np.shape(phi)'+str(np.shape(phi)))
        phi2 = abs(phi)**2
        print('zgrid[int(pars[nz0]/2)]'+str(zgrid[int(pars['nz0']/2)]))
        phi2_outboard = phi2[int(pars['nz0']/2),:,:]
        print('np.shape(phi2_outboard)'+str(np.shape(phi2_outboard)))
        phi2_ob_ky = np.sum(phi2_outboard,axis=1)
        
        apar=field.apar()[:,:,:]

        apar=abs(apar)
        apar_outboard = apar[int(pars['nz0']/2),:,:]
        apar_ob_ky = np.sum(apar_outboard,axis=1)

        apar2=abs(apar)**2
        apar2_outboard = apar2[int(pars['nz0']/2),:,:]
        apar2_ob_ky = np.sum(apar2_outboard,axis=1)

        B1_ob_ky=kygrid_ob*apar_ob_ky*Apar_to_B1*B_gauss  #B1 in Gauss



        plt.clf()
        plt.plot(kygrid,phi2_ob_ky,label=r'$\phi^2$')
        plt.title(r'$Outboard\ \phi^2$'+' at t= '+str(time[itime])+'ms')
        plt.xlabel(r'$k_y(\rho_s)$')
        plt.legend()
        if time_scan==True:
            plt.savefig('pic/Phi_t='+str(time[itime])+'ms.png')
        else:
            plt.show()

        
        
    
        plt.clf()
        plt.plot(kygrid,apar2_ob_ky,label=r'$A_{||}^2$')
        plt.title(r'$Outboard\ A_{||}^2$'+' at t= '+str(time[itime])+'ms')
        plt.xlabel(r'$k_y(\rho_s)$')       
        plt.legend()
        if time_scan==True:
            plt.savefig('pic/Apar_t='+str(time[itime])+'ms.png')
        else:
            plt.show()

        plt.clf()
        plt.plot(kygrid,B1_ob_ky,label=r'$B_{R}$')
        plt.title(r'$Outboard\ B_{R}$'+' at t= '+str(time[itime])+'ms')
        plt.xlabel(r'$k_y(\rho_s)$')       
        plt.ylabel(r'$B_R(Gauss)$')
        plt.legend()
        if time_scan==True:
            plt.savefig('pic/B1_t='+str(time[itime])+'ms.png')
        else:
            plt.show()

    
    else:  #x_local = False
        pass


if time_scan==False:
    sys.exit()

#*******Makting animation!!!*****
ims_phi=[]
ims_apar=[]
ims_B1=[]
for time0 in time_list:
    if time0 == -1:
        itime = -1
        itime0 = len(time)-1
    else: 
        itime = np.argmin(abs(time - time0))
        itime0 = itime
    print("Looking at the spectra at time:"+str(time[itime])+'ms')

    

    file_name='pic/Phi_t='+str(time[itime])+'ms.png'
    ims_phi.append(imageio.imread(file_name))

    file_name='pic/Apar_t='+str(time[itime])+'ms.png'
    ims_apar.append(imageio.imread(file_name))

    file_name='pic/B1_t='+str(time[itime])+'ms.png'
    ims_B1.append(imageio.imread(file_name))


imageio.mimwrite('pic/0_phi_dynamic_images.gif', ims_phi)
imageio.mimwrite('pic/0_apar_dynamic_images.gif', ims_apar)
imageio.mimwrite('pic/0_B1_dynamic_images.gif', ims_B1)

print("Everything is in the directory pic")

