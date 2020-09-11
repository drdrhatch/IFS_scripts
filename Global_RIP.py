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
r_BES=0.98    #x/a of the location of BES
z_BES=0.       #Height of the location of BES, Z=0 being the outboard midplan
BES_radius=0.02  #BES measurement radius in m
#Choose a location to 
RIP_Z_location=0. #mid-plan
beam_width=0.02 #2cm for the laser beam width

Parameter_demo_plot=False  #Change to True if one wants to see the plot of equilibrium
BES_demo_plot=False  #Change to True if one wants to see the plot of the location where the BES data is collected
RIP_demo_plot=True   #Change to True if one wants to see the plot of the RIP and the location where the RIP data is collected

scan_along_z=False  #Change to True if one to scan the RIP result across different height
bin_smooth=10       #size of the bin for rolling average to smooth the data
#
#del_x= 0.135 #13.5cm for DIIID RIP
#real_grid= 0.02 # unit: meter, resolution of the the line integral on height

#***************End of Block for the user****************************
#********************************************************************


def Read_parameter(suffix,plot):

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

    return real_R,real_Z,xgrid,zgrid,B0,B1,n0,n1,gxx,gxy,gyy,gyz,gzz

def BES_real_space(z_BES,r_BES,BES_radius,xgrid,zgrid,real_R,real_Z,gxx,gxy,gyy,gyz,gzz,B0,B1,n0,n1,plot):
    
    (nz,nx)=np.shape(real_R)
    index_x=np.argmin(abs(xgrid-r_BES))
    index_z=np.argmin(abs(real_Z[int(nz/4):int(3*nz/4),index_x]-z_BES))+int(nz/4)
    
    BES_highlight_R=[]
    BES_highlight_Z=[]

    R_BES=real_R[index_z,index_x]
    Z_BES=real_Z[index_z,index_x]
    
    n1_BES=0.
    n0_BES=0.
    B0_BES=0.

    for i in range(nz):
        for j in range(nx):
            if ((real_R[i,j]-R_BES)**2.+(real_Z[i,j]-Z_BES)**2.)**(1./2.)<=BES_radius:
                #area=real_R[i,j]
                n1_BES=n1_BES+n1[i,j]
                n0_BES=n0_BES+n0[i,j]
                B0_BES=B0_BES+B0[i,j]
                BES_highlight_Z.append(real_Z[i,j])
                BES_highlight_R.append(real_R[i,j])

    if plot==True:
        plt.clf()
        plt.ylabel(r'$Height(m)$',fontsize=10)
        plt.xlabel(r'$Major\ raduis(m)$',fontsize=10)
        plt.figure(
        figsize=(4*(np.max(real_R)-np.min(real_R)), 4*(np.max(real_Z)-np.min(real_Z))),
        dpi=960)
        plt.contourf(real_R,real_Z,n1)
        for i in range(len(BES_highlight_R)):
            plt.plot(BES_highlight_R[i], BES_highlight_Z[i],'bo')
        plt.title('BES_highlight',fontsize=10)
        #plt.xlim(1.0,1.3)
        #plt.ylim(-0.03,0.03)
        plt.savefig('BES_highlight.png')

        plt.clf()
        plt.ylabel(r'$Height(m)$',fontsize=10)
        plt.xlabel(r'$Major\ raduis(m)$',fontsize=10)
        plt.figure(
        figsize=(4*(np.max(real_R)-np.min(real_R)), 4*(np.max(real_Z)-np.min(real_Z))),
        dpi=960)
        plt.contourf(real_R,real_Z,n1)
        for i in range(len(BES_highlight_R)):
            plt.plot(BES_highlight_R[i], BES_highlight_Z[i],'bo')
        plt.title('BES_highlight',fontsize=10)
        plt.xlim(1.9,2.3)
        plt.ylim(-0.1,0.1)
        plt.savefig('BES_highlight_zoom.png')

    return n0_BES,n1_BES,B0_BES

def run_BES(z_BES,r_BES,BES_radius,xgrid,zgrid,real_R,real_Z,gxx,gxy,gyy,gyz,gzz,B0,B1,n0,n1,plot):
    for i in range(10):
        n0_BES,n1_BES,B0_BES=BES_real_space(z_BES,r_BES,float(i+1)*BES_radius,xgrid,zgrid,gxx,gxy,gyy,gyz,gzz,real_R,real_Z,B0,B1,n0,n1,plot)
        if n0_BES>0.0001:
            break
    return n0_BES,n1_BES,B0_BES

#corners = [(2.0, 1.0), (4.0, 5.0), (7.0, 8.0)]
def PolygonArea(corners):
    n = len(corners) # of corners
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area

def g_RIP_single_location(RIP_Z_location,beam_width,real_R,real_Z,gxx,gxy,gyy,gyz,gzz,B0,B1,n0,n1,plot):
    integrate_B1_n0_dR=0.
    integrate_n0_dR=0.
    (nz,nx)=np.shape(B0)
    RIP_highlight_R=[]
    RIP_highlight_Z=[]
    area_list=[]
    area_tot=0
    
    for n in range(nz):
        for m in range(nx):
            if real_Z[n,m]>=RIP_Z_location-0.5*beam_width and real_Z[n,m]<=RIP_Z_location+0.5*beam_width:
                n_temp=n
                m_temp=m
                if n >= nz-1:
                    n_temp=nz-2
                if m >= nx-1:
                    m_temp=nx-2
                corners=[(real_R[n_temp,m_temp], real_Z[n_temp,m_temp]),\
                         (real_R[n_temp+1,m_temp], real_Z[n_temp+1,m_temp]),\
                         (real_R[n_temp+1,m_temp+1], real_Z[n_temp+1,m_temp+1]),\
                         (real_R[n_temp,m_temp+1], real_Z[n_temp,m_temp+1])]
                area=PolygonArea(corners)

                print("("+str(real_R[n,m])+", "+str(real_Z[n,m])+")")
                print("area="+str(area))

                area_tot=area_tot+area
                integrate_B1_n0_dR=integrate_B1_n0_dR+B1[n,m]*n0[n,m]*area
                integrate_n0_dR=integrate_n0_dR+n0[n,m]*area
                RIP_highlight_Z.append(real_Z[n,m])
                RIP_highlight_R.append(real_R[n,m])
                area_list.append(area)

    integrate_B1_n0_dR=integrate_B1_n0_dR/area_tot
    integrate_n0_dR=integrate_n0_dR/area_tot

    if plot==True:
        plt.clf()
        plt.ylabel(r'$Height(m)$',fontsize=10)
        plt.xlabel(r'$Major\ raduis(m)$',fontsize=10)
        plt.figure(
        figsize=(4*(np.max(real_R)-np.min(real_R)), 4*(np.max(real_Z)-np.min(real_Z))),
        dpi=960)
        plt.contourf(real_R,real_Z,B1)
        for i in range(len(RIP_highlight_R)):
            plt.plot(RIP_highlight_R[i], RIP_highlight_Z[i],'bo')#,alpha=area_list[i]*30.)
        plt.title('RIP_highlight',fontsize=10)
        plt.xlim(1.0,1.3)
        plt.ylim(-beam_width,beam_width)
        plt.savefig('RIP_highlight_inboard.png')
     
        plt.clf()
        plt.ylabel(r'$Height(m)$',fontsize=10)
        plt.xlabel(r'$Major\ raduis(m)$',fontsize=10)
        plt.figure(
        figsize=(4*(np.max(real_R)-np.min(real_R)), 4*(np.max(real_Z)-np.min(real_Z))),
        dpi=960)
        plt.contourf(real_R,real_Z,B1)
        for i in range(len(RIP_highlight_R)):
            plt.plot(RIP_highlight_R[i], RIP_highlight_Z[i],'bo')#,alpha=area_list[i]*30.)
        plt.title('RIP_highlight',fontsize=10)
        plt.xlim(1.9,2.3)
        plt.ylim(-beam_width,beam_width)
        plt.savefig('RIP_highlight_outboard.png')
        #plt.show()

    return integrate_B1_n0_dR,integrate_n0_dR

def run_RIP_BES_Z_location(suffix,RIP_Z_location,beam_width,z_BES,r_BES,BES_radius,\
                           Parameter_demo_plot,BES_demo_plot,RIP_demo_plot):
    real_R,real_Z,xgrid,zgrid,B0,B1,n0,n1,gxx,gxy,gyy,gyz,gzz=Read_parameter(suffix,plot=Parameter_demo_plot)
    n0_BES,n1_BES,B0_BES=run_BES(z_BES,r_BES,BES_radius,xgrid,zgrid,real_R,real_Z,gxx,gxy,gyy,gyz,gzz,B0,B1,n0,n1,plot=BES_demo_plot)
    integrate_B1_n0_dR,integrate_n0_dR=g_RIP_single_location(RIP_Z_location,beam_width,real_R,real_Z,gxx,gxy,gyy,gyz,gzz,B0,B1,n0,n1,plot=RIP_demo_plot)
    print(n0_BES,n1_BES,B0_BES,integrate_B1_n0_dR,integrate_n0_dR)
    ratio=integrate_B1_n0_dR/(integrate_n0_dR*B0_BES)  *  (n0_BES/n1_BES)
    print("***********************")
    print("***********************")
    print("***********************")
    print("(dB/B)/(dn/n)="+str(ratio))

    return integrate_B1_n0_dR,integrate_n0_dR,ratio,n0_BES,n1_BES,B0_BES


def run_RIP_BES_Z_scan(suffix,RIP_Z_location,beam_width,z_BES,r_BES,BES_radius,\
                           Parameter_demo_plot,BES_demo_plot,RIP_demo_plot):
    real_R,real_Z,xgrid,zgrid,B0,B1,n0,n1,gxx,gxy,gyy,gyz,gzz=Read_parameter(suffix,plot=Parameter_demo_plot)
    n0_BES,n1_BES,B0_BES=run_BES(z_BES,r_BES,BES_radius,xgrid,zgrid,real_R,real_Z,gxx,gxy,gyy,gyz,gzz,B0,B1,n0,n1,plot=BES_demo_plot)
    
    Z_grid=np.arange(np.min(real_Z), np.max(real_Z), beam_width/2.)
    #Z_grid=np.arange(-0.1, 0.1, beam_width/2.)
    integrate_B1_n0_dR_list=[]
    integrate_n0_dR_list=[]
    ratio_list=[]
    with open('RIP_ratio.csv', 'w') as csvfile:
        RIP_data = csv.writer(csvfile, delimiter=',')
        RIP_data.writerow(['Z(m)','Ratio','integrate_B1_n0_dR','integrate_n0_dR','n0_BES','n1_BES','B0_BES'])

        for RIP_location in Z_grid:
            integrate_B1_n0_dR,integrate_n0_dR=g_RIP_single_location(RIP_location,beam_width,real_R,real_Z,gxx,gxy,gyy,gyz,gzz,B0,B1,n0,n1,plot=RIP_demo_plot)
            integrate_B1_n0_dR_list.append(integrate_B1_n0_dR)
            integrate_n0_dR_list.append(integrate_n0_dR)
            ratio=integrate_B1_n0_dR/(integrate_n0_dR*B0_BES)  *  (n0_BES/n1_BES)
            ratio_list.append(ratio)
            RIP_data.writerow([RIP_location,ratio,integrate_B1_n0_dR,integrate_n0_dR,n0_BES,n1_BES,B0_BES])

    csvfile.close()


    print(n0_BES,n1_BES,B0_BES,integrate_B1_n0_dR,integrate_n0_dR)
    
    print("***********************")
    print("***********************")
    print("***********************")
    print("(dB/B)/(dn/n)="+str(ratio))

    return integrate_B1_n0_dR_list,integrate_n0_dR_list,ratio_list,n0_BES,n1_BES,B0_BES,Z_grid



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

if scan_along_z==True:
    integrate_B1_n0_dR_list,integrate_n0_dR_list,ratio_list,n0_BES,n1_BES,B0_BES,Z_grid\
    =run_RIP_BES_Z_scan(suffix,RIP_Z_location,beam_width,z_BES,r_BES,BES_radius,\
                        Parameter_demo_plot,BES_demo_plot,RIP_demo_plot)
    plt.clf()
    plt.ylabel('(dB/B)/(dn/n)',fontsize=10)
    plt.xlabel('Height(m)',fontsize=10)
    plt.plot(Z_grid, ratio_list)
    plt.title(r'$\frac{\int n_e \delta B_r dR}{B_0\int n_e dR} \frac{n_e}{\delta n_e}$',fontsize=10)
    plt.savefig('Ratio.png')
    

    plt.clf()
    plt.ylabel('(dB/B)/(dn/n)',fontsize=10)
    plt.xlabel('Height(m)',fontsize=10)
    plt.plot(smooth(Z_grid,bin_smooth)[0], smooth(ratio_list,bin_smooth)[0])
    plt.title(r'$\frac{\int n_e \delta B_r dR}{B_0\int n_e dR} \frac{n_e}{\delta n_e}$',fontsize=10)
    plt.savefig('Ratio_smooth.png')

integrate_B1_n0_dR,integrate_n0_dR,ratio,n0_BES,n1_BES,B0_BES\
=run_RIP_BES_Z_location(suffix,RIP_Z_location,beam_width,z_BES,r_BES,BES_radius,\
                        Parameter_demo_plot,BES_demo_plot,RIP_demo_plot)
