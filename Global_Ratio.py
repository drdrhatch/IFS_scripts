#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import h5py
import optparse as op
import math
import cmath
import sys
import numpy as np
#import scipy
from scipy import interpolate
import matplotlib.pyplot as plt
from fieldsWrapper import *
from parIOWrapper import init_read_parameters_file
from finite_differences import *
from fieldlib import *
from max_stat_tool import *
from momlib import *
import sys
from nrgWrapper import *
from momentsWrapper_max import *
from read_write_geometry import *
from read_pfile import *
from SI_Gauss_GENE_unit import *
from fieldHelper import *
from max_profile_reader import *
import csv

#Developed by Max Curie on 06/30/2020
#testing path:    /global/cscratch1/sd/maxcurie/global_scan/n0_10
#short-cut:      RIP_global
#short-cut for testing:      RIP_global_test 
#This is the V3 that area weighted

#********************************************************************
#*****************Block for the user*********************************
Outboard_mid_plane=True  #change to True if one wants to only want to look at outboard mid-plane
Delta_Z=0.07  #7cm as bin for Z 
scan_all_Z=False #Change to True if one want to scan across the whole height
max_Z0=0.035
min_Z0=-0.035
#***************End of Block for the user****************************
#********************************************************************

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

def Read_parameter(suffix,plot=False):
    #Initiate the momlib and fieldlib
    pars = init_read_parameters_file(suffix)
    field = fieldfile('field'+suffix,pars)
    moms = momfile('mom_e'+suffix,pars)
    #getting B field using read_write_geometry.py
    gpars,geometry = read_geometry_global(pars['magn_geometry'][1:-1]+suffix)
    #Get geom_coeff from geomWrapper.py
    geom_type, geom_pars, geom_coeff = init_global_geometry(suffix, pars)
    
    real_R=geometry['geo_R'] #it is major radius in meter
    real_Z=geometry['geo_Z'] #it is height in meter, midland is 0, down is negative ,up is positive

    J=geometry['jacobian'] #Jacobian
    
    gxx=geometry['gxx']
    gxy=geometry['gxy']
    gxz=geometry['gxz']
    gyy=geometry['gyy']
    gyz=geometry['gyz']
    gzz=geometry['gzz']

    #from max_profile_reader.py
    x_a,x_rho_ref,T,n0,omt,omn  = profile_e_info(suffix)

    #From plot mode structures

    #Setup the field file
    #************************Setting up the time*****************
    time0=float(options.time0)
    time = np.array(field.tfld)
    timemom = np.array(moms.tmom)
    if time0 == -1:
        itime = -1
        itime0 = len(time)-1
    else: 
        itime = np.argmin(abs(time - time0))
        itime0 = itime
    print(("Looking at the RIP at time:",time[itime]))
    #field.set_time(time[itime],itime0)
    field.set_time(time[itime])
    moms.set_time(timemom[itime])

    upar,deln,deln_global= LILO_moments_from_mom_file(pars,suffix,False,setTime=-1)

    #********************************************

    dz = 2.0/field.nz
    zgrid = np.arange(field.nz)/float(field.nz-1)*(2.0-dz)-1.0
    zgrid_ext = np.arange(field.nz+4)/float(field.nz+4-1)*(2.0+3*dz)-(1.0+2.0*dz)
    #print zgrid
    #print zgrid_ext
    if 'lx_a' in pars:
        xgrid = np.arange(field.nx)/float(field.nx-1)*pars['lx_a']+pars['x0']-pars['lx_a']/2.0
    else:
        xgrid = np.arange(field.nx)/float(field.nx-1)*pars['lx'] - pars['lx']/2.0

    nx=len(xgrid) #number of grid in x axis
    nz=len(zgrid) #number of grid in z axis

    if nx < 10 or nz < 10:
        print("Please increase the resolution")
        return 0
 

    #density
    n0_GENE = np.tile(n0, (nz, 1))
    #print(n0_GENE)

    #B field
    B0_GENE=geometry['Bfield']

    #delta density, check on page 57
    upar,deln,deln_global= LILO_moments_from_mom_file(pars,suffix,False,setTime=-1)
    n1_GENE=abs(deln_global[:,0,:])

    #print geometry['geo_R']
    #print np.shape(geometry['geo_R'])

    #getting phi averaged apar and delta n
    (i1,i2,i3)=np.shape(field.apar())
    #print((np.shape(deln_global)))
    #print(np.shape(field.apar()))
    #Apar_GENE = np.zeros((i1,i3))
    #for i in range(i1):
        #Apar_GENE = Apar_GENE + field.apar()[i,0,:]
        #n1_GENE   = n1_GENE   + deln_global[]
    #Apar_GENE=Apar_GENE/(i2+1)
    Apar_GENE = abs(field.apar()[:,0,:])
    
    #*****************************************************************
    #************Normalized the density and magnetic field************
    B0=np.zeros((nz,nx))
    B1=np.zeros((nz,nx))
    n0=np.zeros((nz,nx))
    n1=np.zeros((nz,nx))
    #ky_GENE=np.zeros(np.shape(deln_global))

    B_gauss=10.**4               #1T=10^4Gauss
    qref = 1.6E-19              #in C
    c  = 1.                     #in 3*10^8m/s
    m_kg = 1.673E-27            #in kg
    Bref = pars['Bref']         #in Tesla
    Tref = pars['Tref']         #in keV
    nref = pars['nref']         #in 10^(19) /m^3
    Lref = pars['Lref']         #in m
    mref = pars['mref']         #in proton mass(kg)
    nref = nref * 1.E19         #in the unit of /m^3
    Tref = Tref * qref * 1000.   #in the unit of J
    mref = mref * m_kg          #in the unit of kg
    pref = nref * Tref          #in Pa*kB
    cref = np.sqrt(Tref / mref) #Speed of sound in the unit of m/s
    Omegaref = qref * Bref / (mref * c)  #Gyrofrequency in rad/s
    rhoref = cref / Omegaref           #gyroradius in m
    rhorefStar = rhoref / Lref         #gyroradius/minor radius in Unitless
    
    for x in range(0,nx):
        #ky_global comes from geomWrapper.py
        ky_GENE_temp=ky_global(pars, geom_coeff, x)
        
        #print("calc Br")
        #print(x)
        for z in range(0,nz):
            B0[z,x]=abs(B0_GENE[z,x]*Bref*B_gauss)
            B1[z,x]=abs(Apar_GENE[z,x]*ky_GENE_temp[z]*Bref*B_gauss*rhorefStar)
            n0[z,x]=abs(n0_GENE[z,x]*nref)
            n1[z,x]=abs(n1_GENE[z,x]*(n0_GENE[z,x]*rhorefStar)*nref)
            #print("rhoi="+str(rhoref)+"m")
            #print("Bref="+str(Bref))
            #print("B_gauss="+str(B_gauss))
            #print("Factor="+str(Bref*B_gauss*rhorefStar))
    #************End of Normalized the density and magnetic field************
    if plot==True:
        plt.clf()
        plt.ylabel(r'$Height(m)$',fontsize=10)
        plt.xlabel(r'$Major\ raduis(m)$',fontsize=10)
        plt.figure(
        figsize=(4*(np.max(real_R)-np.min(real_R)), 4*(np.max(real_Z)-np.min(real_Z))),
        dpi=96)
        plt.contourf(real_R,real_Z,B0)
        plt.title('B0 in real space',fontsize=10)
        plt.savefig('B0.png')

        plt.clf()
        plt.ylabel(r'$Height(m)$',fontsize=10)
        plt.xlabel(r'$Major\ raduis(m)$',fontsize=10)
        plt.figure(
        figsize=(4*(np.max(real_R)-np.min(real_R)), 4*(np.max(real_Z)-np.min(real_Z))),
        dpi=96)
        plt.contourf(real_R,real_Z,B1)
        plt.title('B1 in real space',fontsize=10)
        plt.savefig('B1.png')

        plt.clf()
        plt.ylabel(r'$Height(m)$',fontsize=10)
        plt.xlabel(r'$Major\ raduis(m)$',fontsize=10)
        plt.figure(
        figsize=(4*(np.max(real_R)-np.min(real_R)), 4*(np.max(real_Z)-np.min(real_Z))),
        dpi=96)
        plt.contourf(real_R,real_Z,n0)
        plt.title('n0 in real space',fontsize=10)
        plt.savefig('n0.png')

        plt.clf()
        plt.ylabel(r'$Height(m)$',fontsize=10)
        plt.xlabel(r'$Major\ raduis(m)$',fontsize=10)
        plt.figure(
        figsize=(4*(np.max(real_R)-np.min(real_R)), 4*(np.max(real_Z)-np.min(real_Z))),
        dpi=96)
        plt.contourf(real_R,real_Z,n1)
        plt.title('n1 in real space',fontsize=10)
        plt.savefig('n1.png')

    return J,real_R,real_Z,xgrid,zgrid,B0,B1,n0,n1,gxx,gxy,gyy,gyz,gzz

def Ratio_calc(J,real_R,real_Z,min_Z0,max_Z0,B0,B1,n0,n1,Outboard_mid_plane=False):
    (nz,nx)=np.shape(J)

    if Outboard_mid_plan==True:
        n1_mean=np.mean( (n1*J)[int(nz/2),:] )/np.mean(J[int(nz/2),:])
        n0_mean=np.mean( (n0*J)[int(nz/2),:] )/np.mean(J[int(nz/2),:])
        B1_mean=np.mean( (B1*J)[int(nz/2),:] )/np.mean(J[int(nz/2),:])
        B0_mean=np.mean( (B0*J)[int(nz/2),:] )/np.mean(J[int(nz/2),:])
        BES_highlight_Z=real_Z[int(nz/2),:]
        BES_highlight_R=real_R[int(nz/2),:]
    else:
        n1_list=[]
        n0_list=[]
        B1_list=[]
        B0_list=[]
        J_list=[]
        for i in range(nz):
            for j in range(nx):
                if real_Z[n,m]>=RIP_Z_location-0.5*beam_width and real_Z[n,m]<=RIP_Z_location+0.5*beam_width:
                    #area=real_R[i,j]
                    n1_list.append(n1[i,j]*J[i,j])
                    n0_list.append(n0[i,j]*J[i,j])
                    B1_list.append(B1[i,j]*J[i,j])
                    B0_list.append(B0[i,j]*J[i,j])
                    J_list.append(J[i,j])
        n1_mean=np.mean(n1_list)/np.mean(J_list)
        n0_mean=np.mean(n0_list)/np.mean(J_list)
        B1_mean=np.mean(B1_list)/np.mean(J_list)
        B0_mean=np.mean(B0_list)/np.mean(J_list)

    return n1_mean,n0_mean,B1_mean,B0_mean

suffix=get_suffix()

J,real_R,real_Z,xgrid,zgrid,B0,B1,n0,n1,gxx,gxy,gyy,gyz,gzz=\
    Read_parameter(suffix)

if scan_all_Z==True:
    min_Z0=min(real_Z)
    max_Z0=max(real_Z)

max_Z=max_Z0*1.00001    #Add a small number so it is even
min_Z=min_Z0

Z_grid=np.arange(min_Z,max_Z,Delta_Z)
Z_list = Z_grid[:-1]+Delta_Z/2.
print("Z_list: "+str(Z_list))


with open('Ratio'+suffix+'.csv', 'w') as csvfile:     #clear all and then write a row
    data = csv.writer(csvfile, delimiter=',')
    data.writerow(['Height(m)','n1','n0','B1','B0','(B1/B0)/(n1/n0)'])
csvfile.close()

if Outboard_mid_plane==True:
    n1,n0,B1,B0 =\
         Ratio_calc(J,real_R,real_Z,min_Z0,max_Z0,B0,B1,n0,n1,Outboard_mid_plane)

    with open('Ratio'+suffix+'.csv', 'a') as csvfile:     #clear all and then write a row
        data = csv.writer(csvfile, delimiter=',')
        data.writerow(['mid_plane',n1,n0,B1,B0,(B1/B0)/(n1/n0)])
    csvfile.close()
else:
    for i_Z_list in range(len(Z_list)):
        max_Z0=Z_list[i_Z_list]+Delta_Z/2.    #in the unit of meter
        min_Z0=Z_list[i_Z_list]-Delta_Z/2.   #in the unit of meter
        n1,n0,B1,B0 =\
            Ratio_calc(J,real_R,real_Z,min_Z0,max_Z0,B0,B1,n0,n1,Outboard_mid_plane)

        with open('Ratio'+suffix+'.csv', 'a') as csvfile:     #clear all and then write a row
            data = csv.writer(csvfile, delimiter=',')
            data.writerow([Z_list[i_Z_list],n1,n0,B1,B0,(B1/B0)/(n1/n0)])
        csvfile.close()