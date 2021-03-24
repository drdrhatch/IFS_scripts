#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import h5py
import optparse as op
import math
import cmath
import sys
import numpy as np
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
from geomWrapper import *
from ParIO import *

def Br_Gauss(apar,suffix,ky_GENE):#From GENE Apar to Gauss B_r

    #***********1 GENE unit to SI unit 
    #**** For more info, check on GENE manual page 72
    #par = Parameters()
    #par.Read_Pars('parameters'+suffix)
    #pars = par.pardict
    pars = init_read_parameters_file(suffix)
    #field = fieldfile('field'+suffix,pars)
    #***********geom file
    #gpars,geometry = read_geometry_global(pars['magn_geometry'][1:-1]+suffix)
    #geom_type, geom_pars, geom_coeff = init_global_geometry(suffix, pars)
    qref = 1.6E-19
    c  = 1.
    m_kg = 1.673E-27
    Bref = pars['Bref']
    Tref = pars['Tref']
    nref = pars['nref']
    Lref = pars['Lref']
    mref = pars['mref']
    nref = nref * 1.E19
    Tref = Tref * qref * 1.E03
    mref = mref * m_kg
    pref = nref * Tref
    cref = np.sqrt(Tref / mref)
    Omegaref = qref * Bref / mref / c
    rhoref = cref / Omegaref
    rhorefStar = rhoref / Lref

    Apar_norm = Bref * Lref * rhorefStar **2

    #print rho_ref
    #print rho_ref_star

    A_SI=apar*Apar_norm
    #print'*********SI test***********'
    #print apar
    #print A_SI
    #print '************************'

    #ky_SI=ky_global(pars, geom_coeff, x)[z]
    #print ky_SI
    #print ky_SI
    ky_SI=ky_GENE/rhoref
    B_r=ky_SI*A_SI*B_gauss
    #print B_r

    return B_r

def Br_Gauss_group(apar,suffix,ky_GENE):#From GENE Apar to Gauss B_r

    #***********1 GENE unit to SI unit 
    #**** For more info, check on GENE manual page 72
    #par = Parameters()
    #par.Read_Pars('parameters'+suffix)
    #pars = par.pardict
    pars = init_read_parameters_file(suffix)
    #field = fieldfile('field'+suffix,pars)
    #***********geom file
    #gpars,geometry = read_geometry_global(pars['magn_geometry'][1:-1]+suffix)
    #geom_type, geom_pars, geom_coeff = init_global_geometry(suffix, pars)
    qref = 1.6E-19
    c  = 1.
    m_kg = 1.673E-27
    Bref = pars['Bref']
    Tref = pars['Tref']
    nref = pars['nref']
    Lref = pars['Lref']
    mref = pars['mref']
    nref = nref * 1.E19
    Tref = Tref * qref * 1.E03
    mref = mref * m_kg
    pref = nref * Tref
    cref = np.sqrt(Tref / mref)
    Omegaref = qref * Bref / mref / c
    rhoref = cref / Omegaref
    rhorefStar = rhoref / Lref

    Apar_norm = Bref * Lref * rhorefStar **2

    #print rho_ref
    #print rho_ref_star

#***************For global group

    #for x in range(0,nx):
        #ky_global comes from geomWrapper.py
        #ky_GENE_temp=ky_global(pars, geom_coeff, x)

#***************For global group

    A_SI=apar*Apar_norm
    #print'*********SI test***********'
    #print apar
    #print A_SI
    #print '************************'

    #ky_SI=ky_global(pars, geom_coeff, x)[z]
    #print ky_SI
    #print ky_SI
    ky_SI=ky_GENE/rhoref
    B_r=ky_SI*A_SI*B_gauss
    #print B_r

    return B_r
def Br_Gauss_local_group(apar,suffix,ky_GENE):#From GENE Apar to Gauss B_r

    #***********1 GENE unit to SI unit 
    #**** For more info, check on GENE manual page 72
    #par = Parameters()
    #par.Read_Pars('parameters'+suffix)
    #pars = par.pardict
    pars = init_read_parameters_file(suffix)
    #field = fieldfile('field'+suffix,pars)
    #***********geom file
    #gpars,geometry = read_geometry_global(pars['magn_geometry'][1:-1]+suffix)
    #geom_type, geom_pars, geom_coeff = init_global_geometry(suffix, pars)
    qref = 1.6E-19
    c  = 1.
    m_kg = 1.673E-27
    Bref = pars['Bref']
    Tref = pars['Tref']
    nref = pars['nref']
    Lref = pars['Lref']
    mref = pars['mref']
    nref = nref * 1.E19
    Tref = Tref * qref * 1.E03
    mref = mref * m_kg
    pref = nref * Tref
    cref = np.sqrt(Tref / mref)
    Omegaref = qref * Bref / mref / c
    rhoref = cref / Omegaref
    rhorefStar = rhoref / Lref

    Apar_norm = Bref * Lref * rhorefStar **2

    #print rho_ref
    #print rho_ref_star

    A_SI=apar*Apar_norm
    #print'*********SI test***********'
    #print apar
    #print A_SI
    #print '************************'

    #ky_SI=ky_global(pars, geom_coeff, x)[z]
    #print ky_SI
    #print ky_SI
    ky_SI=ky_GENE/rhoref
    B_r=ky_SI*A_SI*B_gauss
    #print B_r

    return B_r

def Br_Gauss_local_factor(apar,suffix,ky_GENE):#From GENE Apar to Gauss B_r

    #***********1 GENE unit to SI unit 
    #**** For more info, check on GENE manual page 72
    #par = Parameters()
    #par.Read_Pars('parameters'+suffix)
    #pars = par.pardict
    pars = init_read_parameters_file(suffix)
    #field = fieldfile('field'+suffix,pars)
    #***********geom file
    #gpars,geometry = read_geometry_global(pars['magn_geometry'][1:-1]+suffix)
    #geom_type, geom_pars, geom_coeff = init_global_geometry(suffix, pars)
    qref = 1.6E-19
    c  = 1.
    m_kg = 1.673E-27
    Bref = pars['Bref']
    Tref = pars['Tref']
    nref = pars['nref']
    Lref = pars['Lref']
    mref = pars['mref']
    nref = nref * 1.E19
    Tref = Tref * qref * 1.E03
    mref = mref * m_kg
    pref = nref * Tref
    cref = np.sqrt(Tref / mref)
    Omegaref = qref * Bref / mref / c
    rhoref = cref / Omegaref
    rhorefStar = rhoref / Lref

    Apar_norm = Bref * Lref * rhorefStar **2

    #factor=ky_GENE/rhoref*Apar_norm*B_gauss
    factor=ky_GENE*Bref*B_gauss*rhorefStar

    return factor

def B0_Gauss(B0,suffix):#From GENE Apar to Gauss B_r
    #***********1 GENE unit to SI unit 
    #**** For more info, check on GENE manual page 72
    #pars = init_read_parameters_file(suffix)
    par = Parameters()
    par.Read_Pars('parameters'+suffix)
    pars = par.pardict
    B0_SI = B0*(pars['Bref'])*B_gauss
    return B0_SI

def n1_SI(n1,n0,suffix):#From GENE Apar to Gauss B_r
    #***********1 GENE unit to SI unit 
    #**** For more info, check on GENE manual page 72
    #pars = init_read_parameters_file(suffix)
    par = Parameters()
    par.Read_Pars('parameters'+suffix)
    pars = par.pardict
    nref=pars['nref']*10**(19)
    mref=pars['mref']*(1.6726*10**(-27)) #proton mass
    Lref=pars['Lref']
    Tref=pars['Tref']*1000*e_SI  #in keV, Assume kB=1
    Bref=pars['Bref']
    qref=e_SI #in electron charge
    cref=np.sqrt(Tref/mref) #Speed of sound in plasma cs=cref
    Omega_ref=qref*Bref/mref
    rho_ref=cref/Omega_ref
    rho_ref_star=rho_ref/Lref

    #******************
    #ky_GENE= k_y rho_ref
    #rho_ref = c_ref/Omega_ref
    n1_SI=n1*(n0*rho_ref_star)*nref

    return n1_SI

def n0_SI(n0,suffix):#From GENE Apar to Gauss B_r
    #***********1 GENE unit to SI unit 
    #**** For more info, check on GENE manual page 72
    #pars = init_read_parameters_file(suffix)
    par = Parameters()
    par.Read_Pars('parameters'+suffix)
    pars = par.pardict
    nref=pars['nref']*10**(19)
    
    n0_SI=n0*nref

    return n0_SI

def coll_SI(suffix):#From GENE Apar to Gauss B_r
    #***********1 GENE unit to SI unit 
    #**** For more info, check on GENE manual page 72
    #pars = init_read_parameters_file(suffix)
    par = Parameters()
    par.Read_Pars('parameters'+suffix)
    pars = par.pardict
    print('GENE coll')

    print(pars['coll'])

    factor=10**19*(e_SI*1000/kB_SI)**(-2)*10**(-3)
    
    coll_Hz=pars['coll']*factor

    return coll_Hz
def coll_ei_SI(suffix):#From GENE Apar to Gauss B_r
    #***********1 GENE unit to SI unit 
    #**** For more info, check on GENE manual page 72
    #pars = init_read_parameters_file(suffix)
    par = Parameters()
    par.Read_Pars('parameters'+suffix)
    pars = par.pardict

    factor=10**19*(e_SI*1000/kB_SI)**(-2)*10**(-3)
    
    coll_Hz=4*pars['nu_ei']*factor

    return coll_Hz
def omn_SI(suffix):#From GENE Apar to Gauss B_r
    #***********1 GENE unit to SI unit 
    #**** For more info, check on GENE manual page 72
    #pars = init_read_parameters_file(suffix)
    par = Parameters()
    par.Read_Pars('parameters'+suffix)
    pars = par.pardict

    TJ = pars['Tref']*1000.0*1.602e-19
    mi = pars['mref']*1.6726e-27
    cs = np.sqrt(TJ/mi)
    om_ref = cs/pars['Lref']
    omn_Hz=17*om_ref/1000.0/(2.0*np.pi)

    return omn_Hz


#********************************************************
#***********Define the function**************************
#Function average:
#Use the weighted average to "smooth" out the bump and rescaling
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
#*********************


    

#*********Calculation

# Gaussian
#e_gauss=4.8032*10**(-10)
#Gauss_to_SI=np.pi*e_gauss**4*1/(2**(3/2))*1/(Charge_gauss**4*Density_gauss*Length_gauss/(T_gauss)**2)
#print Gauss_to_SI
#SI_to_GENE=Gauss_to_SI*qref**4*nref*Lref*(kB_SI/10**3)**(-2)
#SI_to_GENE_e_SI=Gauss_to_SI*nref*Lref*(kB_SI/(10**3*e_SI))**(-2)
#SI_to_GENE_e_1=Gauss_to_SI*qref**4*nref*Lref*(kB_SI/(10**3))**(-2)
#SI_to_GENE_wrong=Gauss_to_SI*qref**4*nref*Lref*(kB_SI/(10**3*e_SI))**(-2)

#print SI_to_GENE_e_SI
#print SI_to_GENE_e_1
#print SI_to_GENE_wrong


#print('This script is transform between SI unit, Gaussian unit and GENE unit')
#print('write down the formular in terms of Charge, Charge_density, E, D, B, H, phi, mass, Engergy, Length, Density, and make sure it is in python form, for instance 2^10=2**10')
#formular = str(raw_input("The formular: "))
#print('Which way to transform? SI_to_Gauss, SI_to_GENE, Gauss_GENE, GENE_to_SI, Gauss_to_SI')
#transform = str(raw_input("Way of transform: "))


#if transform='SI_to_Gauss':
#   SI_to_Gauss(SI)
#elif transform='SI_to_GENE':
#    SI_to_GENE(SI)
#elif transform='SI_to_Gauss':
#    SI_to_Gauss(SI)
#elif transform='SI_to_Gauss':
#    SI_to_Gauss(SI)
#elif transform='SI_to_Gauss':
#    SI_to_Gauss(SI)
#elif transform='SI_to_Gauss':
#    SI_to_Gauss(SI)

#print('Use the following ')

#def SI_to_Gauss(SI):


