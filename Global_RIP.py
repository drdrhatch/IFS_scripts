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
from momentsWrapper import *
from read_write_geometry import *
from read_pfile import *
from SI_Gauss_GENE_unit import *
from fieldHelper import *
from max_profile_reader import *
import csv

#Developed by Max Curie on 06/07/2020
#testing path:    /global/cscratch1/sd/maxcurie/global_scan/n0_10
#short-cut:      RIP_global
#short-cut for testing:      RIP_global_test 
#This is the V2 that include the real BES ratio

def g_RIP(suffix):

# Define parameter
    del_x= 0.135 #13.5cm for DIIID RIP
    beam_width=0.02 #2cm for the laser beam width
    r_BES=0.98    #x/a of the location of BES
    z_BES=0       #theta/pi of the location of BES
    grides_BES=2  #Take average tot_grid=(1+grides_BES)^2 grids around (z_BES, R_BES), 
    real_grid= 0.02 # unit: meter, resolution of the the line integral on height
    #minor radius
    #read-in radial location

#Initiate the momlib and fieldlib
    pars = init_read_parameters_file(suffix)
    field = fieldfile('field'+suffix,pars)
    moms = momfile('mom_e'+suffix,pars)
    #getting B field using read_write_geometry.py
    gpars,geometry = read_geometry_global(pars['magn_geometry'][1:-1]+suffix)
    geom_type, geom_pars, geom_coeff = init_global_geometry(suffix, pars)
    
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
    print(n0_GENE)

    #B field
    B0_GENE=geometry['Bfield']

    #delta density, check on page 57
    upar,deln,deln_global= LILO_moments_from_mom_file(pars,suffix,False,setTime=-1)
    n1_GENE=abs(deln_global[:,0,:])

    #print geometry['geo_R']
    #print np.shape(geometry['geo_R'])

    #getting phi averaged apar and delta n
    (i1,i2,i3)=np.shape(field.apar())
    print((np.shape(deln_global)))
    print(np.shape(field.apar()))
    #Apar_GENE = np.zeros((i1,i3))
    #for i in range(i1):
        #Apar_GENE = Apar_GENE + field.apar()[i,0,:]
        #n1_GENE   = n1_GENE   + deln_global[]
    #Apar_GENE=Apar_GENE/(i2+1)
    Apar_GENE = abs(field.apar()[:,0,:])
    
#**************************************************
#**************Calculating RIP globally**************
#************Normalized the density and magnetic field************
    B0=np.zeros((nz,nx))
    B1=np.zeros((nz,nx))
    n0=np.zeros((nz,nx))
    n1=np.zeros((nz,nx))
    #ky_GENE=np.zeros(np.shape(deln_global))


    qref = 1.6E-19              #in C
    c  = 1.                     #in 3*10^8m/s
    m_kg = 1.673E-27            #in kg
    Bref = pars['Bref']         #in Tesla
    Tref = pars['Tref']         #in keV
    nref = pars['nref']         #in 10^(19) /m^3
    Lref = pars['Lref']         #in m
    mref = pars['mref']         #in proton mass(kg)
    nref = nref * 1.E19         #in the unit of /m^3
    Tref = Tref * qref * 1000   #in the unit of J
    mref = mref * m_kg          #in the unit of kg
    pref = nref * Tref          #in Pa*kB
    cref = np.sqrt(Tref / mref) #in the unit of m/s
    Omegaref = qref * Bref / mref / c  #in rad/s
    rhoref = cref / Omegaref           #in m
    rhorefStar = rhoref / Lref         #Unitless
    
    for x in range(0,nx):
        #ky_global comes from geomWrapper.py
        ky_GENE_temp=ky_global(pars, geom_coeff, x)
        #print ky_GENE_temp
        print("calc Br")
        print(x)
        for z in range(0,nz):
            B0[z,x]=abs(B0_GENE[z,x]*Bref*B_gauss)
            B1[z,x]=abs(Apar_GENE[z,x]*ky_GENE_temp[z]*Bref*B_gauss*rhorefStar)
            n0[z,x]=abs(n0_GENE[z,x]*nref)
            n1[z,x]=abs(n1_GENE[z,x]*(n0_GENE[z,x]*rhorefStar)*nref)



#**************** Finding boundary of real space*************
    #Get to the real space from read_write_geometry.py Use geo_R and geo_z
    print("Dividing to 4 sections starts")
    
    real_R=geometry['geo_R'] #it is major radius in meter
    real_Z=geometry['geo_Z'] #it is height in meter, midland is 0, down is negative ,up is positive

    inner_Z=np.zeros(nz) # inner boundary of the ring
    inner_R=np.zeros(nz) # inner boundary of the ring
    outer_Z=np.zeros(nz) # outer boundary of the ring
    outer_R=np.zeros(nz) # outer boundary of the ring

    for i in range(0,nz):
        inner_Z[i]=real_Z[i,0]
        inner_R[i]=real_R[i,0]
        outer_Z[i]=real_Z[i,nx-1]
        outer_R[i]=real_R[i,nx-1]

    print("Dividing to 4 sections ends")


#***********************Dividing grid for the calculation***********

    #Divid the ring into four pieces/ 3 section - Higher, Middle, Lower/ 8 line segments
    nR_real=nx
    dR=(outer_R[int(nz/2-1)]-inner_R[int(nz/2-1)])/nx
    #Z=np.arange(np.min(outer_Z),np.max(outer_Z),(np.max(outer_Z)-np.min(outer_Z))/nz)
    dZ=(np.max(outer_Z)-np.min(outer_Z))/nz

    o_max=np.argmax(outer_Z)
    o_min=np.argmin(outer_Z)
    i_max=np.argmax(inner_Z)
    i_min=np.argmin(inner_Z)


    #x = real_Z
    #y = real_R
    #xx, yy = np.meshgrid(x, y)
    #B1 = np.sin(xx**2+yy**2)
    #B0 = np.sin(xx**2+yy**2)
    #n1 = np.sin(xx**2+yy**2)
    #n0 = np.sin(xx**2+yy**2)


# *****************line integrate the real space********************
    #print "Start 2D interpolate"
    #B1_r=interpolate.interp2d(real_Z,real_R,B1,kind='cubic')
    #B1_r=interpolate.interp2d(real_Z,real_R,B1,kind='linear')
    #B0_r=interpolate.interp2d(real_Z,real_R,B0,kind='linear')
    #n1_r=interpolate.interp2d(real_Z,real_R,n1,kind='linear')
    #n0_r=interpolate.interp2d(real_Z,real_R,n0,kind='linear')
    #print "End 2D interpolate"

    print("Start interpolate boundary")

    #B0_cood=np.zeros(np.shape(deln_global)) 
    #B1_cood=np.zeros(np.shape(deln_global)) 
    #n0_cood=np.zeros(np.shape(deln_global)) 
    #n1_cood=np.zeros(np.shape(deln_global)) 

    # Fitting the curve

    outer_outboard_Z=outer_Z[o_min:o_max]
    outer_outboard_R=outer_R[o_min:o_max]

    outer_inboard_Z=list(outer_Z[0:o_min+1])+list(outer_Z[o_max-1:nz-1])
    outer_inboard_R=list(outer_R[0:o_min+1])+list(outer_R[o_max-1:nz-1])

    inner_outboard_Z=inner_Z[i_min:i_max]
    inner_outboard_R=inner_R[i_min:i_max]
    inner_inboard_Z=list(inner_Z[0:i_min+1])+list(inner_Z[i_max:nz-1])
    inner_inboard_R=list(inner_R[0:i_min+1])+list(inner_R[i_max:nz-1])

    outer_outboard=interpolate.interp1d(outer_outboard_Z,outer_outboard_R,kind='cubic')
    outer_inboard =interpolate.interp1d(outer_inboard_Z, outer_inboard_R, kind='cubic')
    inner_outboard=interpolate.interp1d(inner_outboard_Z,inner_outboard_R,kind='cubic')
    inner_inboard =interpolate.interp1d(inner_inboard_Z, inner_inboard_R, kind='cubic')

    print("End interpolate boundary")




    Z=np.arange(np.min(outer_Z), np.max(outer_Z), dZ)
    Z_in=np.arange(np.min(inner_Z), np.max(inner_Z), dZ)

    outer_outboard_cood=np.zeros(len(Z))
    outer_inboard_cood =np.zeros(len(Z))
    inner_outboard_cood=np.zeros(len(Z_in))
    inner_inboard_cood =np.zeros(len(Z_in))
    top_r=np.zeros(len(Z))
    bottom_r=np.zeros(len(Z))
    dens_r=np.zeros(len(Z))
    Ratio_r=np.zeros(len(Z))
    RIP_gauss_r=np.zeros(len(Z))
    Ratio_BES_r=np.zeros(len(Z))
    dn_list=np.zeros(len(Z))
    n_list=np.zeros(len(Z))
    dB_list=np.zeros(len(Z))
    B_list=np.zeros(len(Z))

    j=0
    k=0

    #print len(Z)

    print("Start Calc RIP")


    nz_BES=(np.abs(zgrid - z_BES)).argmin()
    nx_BES=(np.abs(xgrid - r_BES)).argmin()
    n1_BES=np.mean(n1[nz_BES-grides_BES:nx_BES+grides_BES,nx_BES-grides_BES:nx_BES+grides_BES])
    n0_BES=np.mean(n0[nz_BES-grides_BES:nx_BES+grides_BES,nx_BES-grides_BES:nx_BES+grides_BES])
    B0_BES=np.mean(B0[nz_BES-grides_BES:nx_BES+grides_BES,nx_BES-grides_BES:nx_BES+grides_BES])

    for Zi in Z: 
        print(j)

        outer_outboard_cood[j]=outer_outboard(Zi)
        outer_inboard_cood[j]=outer_inboard(Zi)

        temp_top=0
        temp_bottom=0
        temp_dens=0
        temp_bottom_BES=0
        dn_temp=0
        n_temp=0
        B_temp=0
        dB_temp=0
        l=0
        temp_n=0
#**********the Lower section***************
        if Zi < np.min(inner_Z):
            for Ri in np.arange(outer_inboard(Zi), outer_outboard(Zi), dR):

                
                for m in range(0,nx):
                    for n in range(0,nz):
                        if real_R[n,m]>Ri and real_R[n,m]<Ri+dR and real_Z[n,m]>Zi and real_Z[n,m]<Zi+dZ:
                            temp_top=temp_top+B1[n,m]*n0[n,m]
                            temp_bottom=temp_bottom+B0[n,m]*n1[n,m]
                            temp_bottom_BES=temp_bottom_BES+B0[n,m]*n1_BES
                            temp_dens=temp_dens+n0[n,m]
                            dn_temp=dn_temp+n1[n,m]
                            n_temp=n_temp+n0[n,m]
                            B_temp=B_temp+B0[n,m]
                            dB_temp=dB_temp+B1[n,m]
                            temp_n=temp_n+1


                #B0_cood[j,l]=B0_r(Zi,Ri)
                #B1_cood[j,l]=B1_r(Zi,Ri)
                #n0_cood[j,l]=n0_r(Zi,Ri)
                #n1_cood[j,l]=n1_r(Zi,Ri)

            top_r[j]=temp_top
            bottom_r[j]=temp_bottom
            dens_r[j]=temp_dens
            if temp_n==0:
                dn_list[j]=0
                n_list[j]=0
                dB_list[j]=0
                B_list[j]=0
            else:
                dn_list[j]=dn_temp/temp_n
                n_list[j]=n_temp/temp_n
                dB_list[j]=dB_temp/temp_n
                B_list[j]=B_temp/temp_n
            

            l=l+1
            if temp_dens!=0:
                RIP_gauss_r[j]=temp_top/temp_dens
            else:
                RIP_gauss_r[j]=0


            if temp_bottom!=0:
                Ratio_r[j]=temp_top/temp_bottom
            else:
                Ratio_r[j]=0

            if temp_bottom_BES!=0:
                Ratio_BES_r[j]=temp_top/temp_bottom_BES
            else:
                Ratio_BES_r[j]=0

#*********the middle section****************

        elif Zi < np.max(inner_Z):

            inner_outboard_cood[k]=inner_outboard(Zi)
            inner_inboard_cood[k]=inner_inboard(Zi)


            for Ri in np.arange(outer_inboard(Zi), inner_inboard(Zi), dR):
                #B0_cood[j,l]=B0_r(Zi,Ri)
                #B1_cood[j,l]=B1_r(Zi,Ri)
                #n0_cood[j,l]=n0_r(Zi,Ri)
                #n1_cood[j,l]=n1_r(Zi,Ri)

                for m in range(0,nx):
                    for n in range(0,nz):
                        if real_R[n,m]>Ri and real_R[n,m]<Ri+dR and real_Z[n,m]>Zi and real_Z[n,m]<Zi+dZ:
                            temp_top=temp_top+B1[n,m]*n0[n,m]
                            temp_bottom=temp_bottom+B0[n,m]*n1[n,m]
                            temp_bottom_BES=temp_bottom_BES+B0[n,m]*n1_BES
                            temp_dens=temp_dens+n0[n,m]
                            dn_temp=dn_temp+n1[n,m]
                            n_temp=n_temp+n0[n,m]
                            B_temp=B_temp+B0[n,m]
                            dB_temp=dB_temp+B1[n,m]
                            temp_n=temp_n+1

                l=l+1

            for Ri in np.arange(inner_outboard(Zi), outer_outboard(Zi), dR):
                #B0_cood[j,l]=B0_r(Zi,Ri)
                #B1_cood[j,l]=B1_r(Zi,Ri)
                #n0_cood[j,l]=n0_r(Zi,Ri)
                #n1_cood[j,l]=n1_r(Zi,Ri)

                for m in range(0,nx):
                    for n in range(0,nz):
                        if real_R[n,m]>Ri and real_R[n,m]<Ri+dR and real_Z[n,m]>Zi and real_Z[n,m]<Zi+dZ:
                            temp_top=temp_top+B1[n,m]*n0[n,m]
                            temp_bottom=temp_bottom+B0[n,m]*n1[n,m]
                            temp_bottom_BES=temp_bottom_BES+B0[n,m]*n1_BES
                            temp_dens=temp_dens+n0[n,m]
                            dn_temp=dn_temp+n1[n,m]
                            n_temp=n_temp+n0[n,m]
                            B_temp=B_temp+B0[n,m]
                            dB_temp=dB_temp+B1[n,m]
                            temp_n=temp_n+1

                l=l+1

            top_r[j]=temp_top
            bottom_r[j]=temp_bottom
            dens_r[j]=temp_dens
            if temp_n==0:
                dn_list[j]=0
                n_list[j]=0
                dB_list[j]=0
                B_list[j]=0
            else:
                dn_list[j]=dn_temp/temp_n
                n_list[j]=n_temp/temp_n
                dB_list[j]=dB_temp/temp_n
                B_list[j]=B_temp/temp_n

            if temp_dens!=0:
                RIP_gauss_r[j]=temp_top/temp_dens
            else:
                RIP_gauss_r[j]=0


            if temp_bottom!=0:
                Ratio_r[j]=temp_top/temp_bottom
            else:
                Ratio_r[j]=0

            if temp_bottom_BES!=0:
                Ratio_BES_r[j]=temp_top/temp_bottom_BES
            else:
                Ratio_BES_r[j]=0
 
            k=k+1


#********the upper section******************
        else:
            for Ri in np.arange(outer_inboard(Zi), outer_outboard(Zi), dR):
                #B0_cood[j,l]=B0_r(Zi,Ri)
                #B1_cood[j,l]=B1_r(Zi,Ri)
                #n0_cood[j,l]=n0_r(Zi,Ri)
                #n1_cood[j,l]=n1_r(Zi,Ri)

                for m in range(0,nx):
                    for n in range(0,nz):
                        if real_R[n,m]>Ri and real_R[n,m]<Ri+dR and real_Z[n,m]>Zi and real_Z[n,m]<Zi+dZ:
                            temp_top=temp_top+B1[n,m]*n0[n,m]
                            temp_bottom=temp_bottom+B0[n,m]*n1[n,m]
                            temp_bottom_BES=temp_bottom_BES+B0[n,m]*n1_BES
                            temp_dens=temp_dens+n0[n,m]
                            dn_temp=dn_temp+n1[n,m]
                            n_temp=n_temp+n0[n,m]
                            B_temp=B_temp+B0[n,m]
                            dB_temp=dB_temp+B1[n,m]
                            temp_n=temp_n+1

                l=l+1

            top_r[j]=temp_top
            bottom_r[j]=temp_bottom
            dens_r[j]=temp_dens
            if temp_n==0:
                dn_list[j]=0
                n_list[j]=0
                dB_list[j]=0
                B_list[j]=0
            else:
                dn_list[j]=dn_temp/temp_n
                n_list[j]=n_temp/temp_n
                dB_list[j]=dB_temp/temp_n
                B_list[j]=B_temp/temp_n

            if temp_dens!=0:
                RIP_gauss_r[j]=temp_top/temp_dens
            else:
                RIP_gauss_r[j]=0


            if temp_bottom!=0:
                Ratio_r[j]=temp_top/temp_bottom
            else:
                Ratio_r[j]=0

            if temp_bottom_BES!=0:
                Ratio_BES_r[j]=temp_top/temp_bottom_BES
            else:
                Ratio_BES_r[j]=0

        print((RIP_gauss_r[j]))
        j=j+1

    print("End Calc RIP")


    print("Start ploting")
        
    plt.clf()
    plt.ylabel(r'$Height(m)$',fontsize=10)
    plt.xlabel(r'$Major\ raduis(m)$',fontsize=10)
    plt.figure(
    figsize=(4*(np.max(outer_R)-np.min(outer_R)), 4*(np.max(outer_Z)-np.min(outer_Z))),
    dpi=96)
    plt.plot(inner_R,inner_Z)
    plt.plot(outer_R,outer_Z)
    plt.title(r'$boundary\ of\ the\ simulation$',fontsize=10)
    plt.savefig('boundary.png')


    
    plt.clf()
    plt.ylabel(r'$Height(m)$',fontsize=10)
    plt.xlabel(r'$Major\ raduis(m)$',fontsize=10)
    plt.figure(
    figsize=(4*(np.max(outer_R)-np.min(outer_R)), 4*(np.max(outer_Z)-np.min(outer_Z))),
    dpi=96)
    plt.plot(outer_outboard_cood,Z)
    plt.plot(outer_inboard_cood,Z)
    plt.plot(inner_outboard_cood,Z_in)
    plt.plot(inner_inboard_cood,Z_in)
    plt.title(r'$4\ sections\ of\ boundary\ of\ the\ simulation$',fontsize=10)
    plt.savefig('boundary_4_sections.png')

    plt.clf()
    plt.title(r'$Line\ integrated\ magnetic\ fluctuation$')
    plt.xlabel(r'$\frac{\int n_e \delta B_r dR}{\int n_e dR} (Arbitrary unit)$',fontsize=10)
    plt.ylabel(r'$Height (m)$',fontsize=13)
    plt.plot(norm(RIP_gauss_r),Z)
    plt.savefig('RIP_gauss_r.png')
    #plt.show()

    plt.clf()
    plt.title(r'$Line\ integrated\ magnetic\ fluctuation$')
    plt.xlabel(r'$\frac{\int n_e \delta B_r dR}{\int \delta n_e B_0 dR}$',fontsize=10)
    plt.ylabel(r'$Height (m)$',fontsize=13)
    plt.plot(Ratio_r,Z)
    plt.savefig('Ratio_r.png')
    #plt.show()
    
    Ratio_BES_r=RIP_gauss_r*n0_BES/(B0_BES*n1_BES)

    #*****Start save the data********* 
    with open('RIP_gauss.csv', 'w') as csvfile:
        nz=len(Z)
        RIP_data = csv.writer(csvfile, delimiter=',')
        RIP_data.writerow(['Z','Ratio_r','Ratio_BES','RIP_gauss','dB','B','dn','n'])

        for j in range(nz):
            RIP_data.writerow([Z[j],Ratio_r[j],Ratio_BES_r[j],RIP_gauss_r[j],dB_list[j],B_list[j],dn_list[j],n_list[j]])
    csvfile.close()
    #************End save the data**************


    plt.clf()
    plt.title(r'$Line\ integrated\ magnetic\ fluctuation$')
    plt.xlabel(r'$\frac{\int n_e \delta B_r dR}{\int n_e dR} (Arbitrary unit)$',fontsize=10)
    plt.ylabel(r'$Height (m)$',fontsize=13)
    plt.plot(norm(smooth(RIP_gauss_r,15)[0]),smooth(Z,15)[0])
    plt.savefig('RIP_gauss_r_smooth.png')
    #plt.show()

    plt.clf()
    plt.title(r'$Line\ integrated\ magnetic\ fluctuation$')
    plt.xlabel(r'$\frac{\int n_e \delta B_r dR}{\int n_e dR} (Arbitrary unit)$',fontsize=10)
    plt.ylabel(r'$Height (m)$',fontsize=13)
    plt.plot(smooth(RIP_gauss_r,15)[0],smooth(Z,15)[0])
    plt.plot(RIP_gauss_r,Z)
    plt.savefig('RIP_gauss_r_compare.png')
    #plt.show()

    plt.clf()
    plt.title(r'$Line\ integrated\ magnetic\ fluctuation$')
    plt.xlabel(r'$\frac{\int n_e \delta B_r dR}{\int \delta n_e B_0 dR}$',fontsize=10)
    plt.ylabel(r'$Height (m)$',fontsize=13)
    plt.plot(smooth(Ratio_r,15)[0],smooth(Z,15)[0])
    plt.savefig('Ratio_r_smooth.png')
    #plt.show()

    plt.clf()
    plt.title(r'$fluctuation Ratio\ Based\ On\ RIP\ and\ BES$')
    plt.xlabel(r'$\frac{\int n_e \delta B_r dR}{B_0\int n_e dR} \frac{n_e}{\delta n_e}$',fontsize=10)
    plt.ylabel(r'$Height (m)$',fontsize=13)
    plt.plot(smooth(Ratio_BES_r,15)[0],smooth(Z,15)[0])
    plt.savefig('Ratio_BES_r_smooth.png')
    #plt.show()

    Z_zoom, RIP_gauss_r_zoom= zoom1D(smooth(Z,15)[0],smooth(RIP_gauss_r,15)[0],-0.25,0.25)
    plt.clf()
    plt.title(r'$Line\ integrated\ magnetic\ fluctuation$')
    plt.xlabel(r'$\frac{\int n_e \delta B_r dR}{\int n_e dR} (Arbitrary\ unit)$',fontsize=10)
    plt.ylabel(r'$Height (m)$',fontsize=13)
    plt.plot(RIP_gauss_r_zoom,Z_zoom)
    plt.savefig('RIP_gauss_r_smooth_zoom.png')
    #plt.show()
    

    Z_zoom, Ratio_BES_r_zoom= zoom1D(smooth(Z,15)[0],smooth(Ratio_BES_r,15)[0],-0.25,0.25)
    plt.clf()
    plt.title(r'$fluctuation Ratio\ Based\ On\ RIP\ and\ BES$')
    plt.xlabel(r'$\frac{\int n_e \delta B_r dR}{B_0\int n_e dR} \frac{n_e}{\delta n_e}$',fontsize=10)
    plt.ylabel(r'$Height (m)$',fontsize=13)
    plt.plot(Ratio_BES_r_zoom,Z_zoom)
    plt.savefig('Ratio_BES_r_smooth_zoom.png')
    #plt.show()

    print("End ploting")



e_SI= 1.6*10**(-19)
kB_SI=1.3807*10**(-23)
c_SI=2.9979*10**8 #Speed of light

#Regarding SI is unity

#*********************
#*******Gaussian unit to 1 SI unit
Charge_gauss=2.998*10**9
Charge_density_gauss=2.998*10**3
E_gauss=1/Charge_gauss
D_gauss=4*np.pi*2.998*10**5
B_gauss=10**4
H_gauss=4*np.pi*10**(-3)
phi_gauss=1/(2.998*10**2)
mass_gauss=10**3
Energy_gauss=10**7
Length_gauss=10**2
Density_gauss=10**(-6)
kD_gauss=1.3807*10**(-17)
T_gauss=1
   

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

g_RIP(suffix)
