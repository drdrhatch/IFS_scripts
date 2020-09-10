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

#Developed by Max Curie on 01/27/2020
#testing path:    /global/cscratch1/sd/maxcurie/global_test/n0_10
#short-cut:      RIP_global
#short-cut for testing:      RIP_global_test 
#V1 is the one does not have the interperlation

def l_RIP(suffix):

# Define parameter
    del_x= 0.135 #13.5cm for DIIID RIP
    beam_width=0.02 #2cm for the laser beam width
    real_grid= 0.02 # unit: meter, resolution of the the line integral on height 


    ratio_location=0. #mid-plan
    aveage_delta=0.010  #average around the z=ratio_location\pm aveage_delta (beam width) in meter
    #minor radius
    #read-in radial location

#Initiate the momlib and fieldlib
    pars = init_read_parameters_file(suffix)
    field = fieldfile('field'+suffix,pars)
    moms = momfile('mom_e'+suffix,pars)
    #getting B field using read_write_geometry.py
    gpars,geometry = read_geometry_local(pars['magn_geometry'][1:-1]+suffix)
    #Get geom_coeff from geomWrapper.py
    geom_type, geom_pars, geom_coeff = init_read_geometry_file(suffix, pars)

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
    print("Looking at the RIP at time:",time[itime])
    #field.set_time(time[itime],itime0)
    field.set_time(time[itime])
    moms.set_time(timemom[itime])

#********************************************
    dz = 2.0/field.nz
    zgrid = np.arange(field.nz)/float(field.nz-1)*(2.0-dz)-1.0

    nz=len(zgrid) #number of grid in z axis
    #density n0
    n0_GENE = pars['nref']
    #from momentsWrapper.py
    #delta density, check on page 57
    upar,deln,deln_global= LILO_moments_from_mom_file(pars,suffix,False,setTime=-1)
    n1_GENE=abs(deln_global[:,0,:]) 
    #B field
    B0_GENE=geometry['gBfield']

    #From fieldlib     
    Apar_GENE = abs(field.apar()[:,0,:])

    #print('**********************')
    #print(f'the size of the n0 is {np.shape(n0_GENE)}')
    #print('**********************')
    #print(f'the size of the deln is {np.shape(n1_GENE)}')
    #print('**********************')
    #print(f'the size of the B0 is {np.shape(B0_GENE)}')
    #print('**********************')
    #print(f'the size of the apar is {np.shape(Apar_GENE)}')   
    
#**************************************************
#**************Calculating RIP globally**************
#************Normalized the density and magnetic field************
    B0=np.zeros(nz)
    B1=np.zeros(nz)
    n1=np.zeros(nz)
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
    Tref = Tref * qref * 1.E03  #in the unit of J
    mref = mref * m_kg          #in the unit of kg
    pref = nref * Tref          #in Pa*kB
    cref = np.sqrt(Tref / mref) #in the unit of m/s
    Omegaref = qref * Bref / mref / c  #in rad/s
    rhoref = cref / Omegaref           #in m rho_i(ion gyroradii)
    rhorefStar = rhoref / Lref         #Unitless
    
    #ky comes from geomWrapper.py
    ky_GENE_temp=ky(pars, geom_coeff)

    n0=nref #nref

    for z in range(nz):
        B0[z]=abs(B0_GENE[z]*Bref*B_gauss)
        B1[z]=abs(np.mean(Apar_GENE[z,:])*len(Apar_GENE[z,:])*(ky_GENE_temp[z]/rhoref)*Bref*B_gauss*rhorefStar*rhoref)
        n1[z]=abs(np.mean(n1_GENE[z,:])*len(n1_GENE[z,:])*(rhorefStar)*nref)

    print('******n0*********')
    print(n0)
    print('*******n1********')
    print(n1)
    print('*******B0********')
    print(B0)
    print('********B1*******')
    print(B1)

    z_loop, B0_loop = loop(zgrid,B0,-0.5,0.5)
    z_loop, B1_loop = loop(zgrid,B1,-0.5,0.5)
    z_loop, n1_loop = loop(zgrid,n1,-0.5,0.5)
    print('*******n0 Loop********')
    print(n0)
    print('*******n1 Loop********')
    print(n1_loop)
    print('*******B0 Loop********')
    print(B0_loop)
    print('********B1 Loop*******')
    print(B1_loop)

    Ratio=np.zeros(len(z_loop))
    B_ratio=np.zeros(len(z_loop))
    n_ratio=np.zeros(len(z_loop))
    RIP_gauss=np.zeros(len(z_loop))
    for z in range(len(z_loop)):
        Ratio[z]=B1_loop[z]*n0/(B0_loop[z]*n1_loop[z])
        RIP_gauss[z]=B1_loop[z]
        B_ratio[z]=B1_loop[z]/B0_loop[z]
        n_ratio[z]=n1_loop[z]/n0


    x0=z_loop
    y0=Ratio
    location_index_max=np.argmin(abs(x0-ratio_location-aveage_delta))
    location_index_min=np.argmin(abs(x0-ratio_location+aveage_delta))
    local_temp=0
    if location_index_max > location_index_min: 
        location_index_max = location_index_max
        location_index_min = location_index_min
    elif location_index_min > location_index_max:
        local_temp = location_index_max
        location_index_max = location_index_min
        location_index_min = local_temp
    else:
        local_temp=np.argmin(abs(x0-ratio_location))
        location_index_max=local_temp
        location_index_min=local_temp-1
        
    print(location_index_min)
    print(location_index_max)
    RIP_0=np.mean(y0[location_index_min:location_index_max])
    print("RIP ratio average at Z="+str(ratio_location)+":  "+str(RIP_0))
    
    #dB_0=np.mean(B_ratio[location_index_min:location_index_max])
    #print("dB ratio average at Z="+str(ratio_location)+":  "+str(dB_0))

    #dn_0=np.mean(n_ratio[location_index_min:location_index_max])
    #print("dn ratio average at Z="+str(ratio_location)+":  "+str(dn_0))


    print("Start ploting")

    plt.clf()
    plt.title(r'$Magnetic\ fluctuation\ ratio$')
    plt.xlabel(r'$z/\pi$',fontsize=13)
    plt.ylabel(r'$\frac{ n_e \delta B_r}{ \delta n_e B_0}$',fontsize=10)
    plt.plot(z_loop,Ratio)
    plt.savefig('Ratio.png')
    #plt.show()

    plt.clf()
    plt.title(r'$Magnetic\ fluctuation$')
    plt.xlabel(r'$z/\pi$',fontsize=13)
    plt.ylabel(r'$\delta B_r$',fontsize=10)
    plt.plot(z_loop,RIP_gauss)
    plt.savefig('RIP_gauss.png')
    #plt.show()
    

    plt.clf()
    plt.title(r'$Magnetic\ fluctuation\ Ratio$')
    plt.xlabel(r'$z/\pi$',fontsize=13)
    plt.ylabel(r'$\delta B_r$',fontsize=10)
    plt.plot(z_loop,B_ratio)
    plt.savefig('dB_ratio.png')
    #plt.show()

    plt.clf()
    plt.title(r'$Density\ fluctuation\ Ratio$')
    plt.xlabel(r'$z/\pi$',fontsize=13)
    plt.ylabel(r'$\delta B_r$',fontsize=10)
    plt.plot(z_loop,n_ratio)
    plt.savefig('dn_ratio.png')
    #plt.show()

    print("End ploting")

parser=op.OptionParser(description='Plots mode structures and calculates various interesting quantities.')
parser.add_option('--plot_theta','-g',action='store_const',const=False,help = 'Plot all plots.',default=True)
parser.add_option('--plot_ballooning','-b',action='store_const',const=False,help = 'Plot all plots.',default=True)
parser.add_option('--plot_all','-a',action='store_const',const=1,help = 'Plot all plots.',default=False)
parser.add_option('--time','-t',type = 'float',action='store',dest="time0",help = 'Time to plot mode structure.',default=-1)
parser.add_option('--idb','-i',type = 'str',action='store',dest="idb_file",help = 'ITERDB file name.',default='empty')
parser.add_option('--pfile','-f',type = 'str',action='store',dest="prof_file",help = 'ITERDB file name.',default='empty')
options,args=parser.parse_args()
    
if len(args)!=1:
    exit("""
Please include run number as argument (e.g., 0001)."
    \n""")
suffix = args[0]
suffix = '_'+suffix
#suffix = args[0]
#A_par('_0002')

l_RIP(suffix)
