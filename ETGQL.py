#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from scipy import integrate
from fieldlib import *
import numpy
import matplotlib.pyplot as plt
import re
import optparse as op
from ParIO import * 
from interp import *
from os.path import exists
from read_write_geometry import read_geometry_local
from get_nrg import get_nrg0

parser=op.OptionParser(description='Constructs the quasilinear estimate used in final publication.  First argument is nonlinear run number (0 if none), second argumetn is scanfiles directory for linear sims.')

parser.add_option("-n", "--no_bessel",
                  action="store_true", dest="no_bessel", default=False,
                  help="Don't include the Bessel approximation in the averaging.")

parser.add_option("-z", "--select_zrange",
                  action="store_true", dest="select_zrange", default=False,
                  help="Select the zrange for averaging.")

options,args=parser.parse_args()
if len(args)!=2:
    exit("""
Please include nonlinear run number as argument (e.g., 0001) and scanfiles suffix as second argument.  If no nonlinear enter 0 in first argument."
    \n""")
suffix = args[0]
sfsuffix = args[1]

no_bessel = options.no_bessel
select_zrange = options.select_zrange
#print("no_bessel",no_bessel)
#stop

if 'dat' in suffix:
   suffix = '.dat'
elif '_' not in suffix:
   suffix = '_' + suffix

def construct_extended_ballooning(pars,field):
    ntot = pars['nx0']*pars['nz0']
    dz = float(2.0)/float(pars['nz0'])
    zgrid = np.pi*(np.arange(ntot)/float(ntot-1)*(2*pars['nx0']-dz)-pars['nx0'])
    field_ext = np.zeros(ntot,dtype='complex128')

    if 'n0_global' in pars:
        phase_fac = -np.e**(-2.0*np.pi*(0.0+1.0J)*pars['n0_global']*pars['q0'])
    else:
        phase_fac = -1.0
    #print "phase_fac",phase_fac

    if pars['shat'] < 0.0:
        for i in range(int(pars['nx0']/2)+1):
            field_ext[(i+int(pars['nx0']/2))*pars['nz0']:(i+int(pars['nx0']/2)+1)*pars['nz0']]=field[:,0,-i]*phase_fac**i
            if i < int(pars['nx0']/2):
                field_ext[(int(pars['nx0']/2)-i-1)*pars['nz0'] : (int(pars['nx0']/2)-i)*pars['nz0'] ]=field[:,0,i+1]*phase_fac**(-(i+1))
    else:
        for i in range(int(pars['nx0']/2)):
            #print("phase_fac**i",phase_fac**i)
            #print("phase_fac**(-(i+1))",phase_fac**(-(i+1)))
            field_ext[(i+int(pars['nx0']/2))*pars['nz0']:(i+int(pars['nx0']/2)+1)*pars['nz0']]=field[:,0,i]*phase_fac**i
            if i < int(pars['nx0']/2):
                field_ext[(int(pars['nx0']/2)-i-1)*pars['nz0'] : (int(pars['nx0']/2)-i)*pars['nz0'] ]=field[:,0,-1-i]*phase_fac**(-(i+1))

    return zgrid,field_ext

def get_jacxB_extended(pars,zgrid,geom):
    nx0 = pars['nx0']
    nz0 = pars['nz0']
    jacxB = geom['gjacobian']*geom['gBfield']
    jacxB_extended = np.empty(len(zgrid))
    for i in range(nx0):
        jacxB_extended[i*nz0:(i+1)*nz0] = jacxB[:]
    return jacxB_extended

#def calc_kperp(pars,geom_coeff):
#
#    nx = int(pars['nx0'])
#    ikx_grid = np.arange(- nx // 2 + 1, nx // 2 + 1)
#    nz = int(pars['nz0'])
#    lx = float(pars['lx'])
#    ky = float(pars['kymin'])
#    dkx = 2. * np.pi * float(pars['shat']) * float(ky)
#
#    if 'kx_center' in pars:
#        kx_center = float(pars['kx_center'])
#    else:
#        kx_center = 0.
#
#    ggxx = geom_coeff['ggxx'].astype(float)
#    ggxy = geom_coeff['ggxy'].astype(float)
#    ggyy = geom_coeff['ggyy'].astype(float)
#    gBfield = geom_coeff['gBfield'].astype(float)
#
#    kperp = np.zeros(nx*nz,dtype='float128') # changed to longdouble ..
#
#    for i in ikx_grid:
#        kx = i*dkx+kx_center
#        this_kperp = np.sqrt(ggxx*kx**2+2.*ggxy*kx*ky+ggyy*ky**2)
#
#        kperp[(i-ikx_grid[0])*nz:(i-ikx_grid[0]+1)*nz]=this_kperp
#    return kperp

def get_kperp(pars,geom):
    nx = pars['nx0']
    ikx_grid = np.arange(- nx // 2 + 1, nx // 2 + 1)
    ky = float(pars['kymin'])
    dkx = 2. * np.pi * float(pars['shat']) * float(ky)

    nz = int(pars['nz0'])
    lx = float(pars['lx'])
    dkx = 2. * np.pi * float(pars['shat']) * float(ky)

    if 'kx_center' in pars:
        kx_center = pars['kx_center']
    else:
        kx_center = 0.

    kperp = np.zeros(nx*nz,dtype='float128')

    for i in ikx_grid:
        kx = i*dkx+kx_center
        this_kperp = np.sqrt(geom['ggxx']*kx**2+2.*geom['ggxy']*kx*ky+geom['ggyy']*ky**2)
        kperp[(i-ikx_grid[0])*nz:(i-ikx_grid[0]+1)*nz]=this_kperp
    return kperp


#def eigenfunction_average_bessel(kperp,quantity,field,pars,geometry,mass_ratio = 1):
#
#    zgrid, field_ext = construct_extended_ballooning(pars,field)
#
#    jacxBpi = geometry['gjacobian']*geometry['gBfield']*np.pi
#    jacxBpi_extended = get_jacxBpi_extended(zgrid,jacxBpi)
#    #mass_ratio = 9.1094e-31/(pars['mref']*1.6726e-27)
#    kperp_bessel = kperp*np.sqrt(mass_ratio)
#
#    alpha = 2./3.
#    bessel_factor = 1. / np.sqrt(1. + 2. * (kperp_bessel**2 + np.pi * alpha * kperp_bessel**4) /
#                    (1. + alpha * kperp_bessel**2))
#
#    ave_quant = 0.
#    denom = 0.
#
#    for i in np.arange(len(field_ext)-1):
#        ave_quant = ave_quant + (quantity[i]*abs(field_ext[i])**2 * bessel_factor[i]+\
#            quantity[i+1]*abs(field_ext[i+1])**2 * bessel_factor[i+1])/2.*\
#            (zgrid[i+1]-zgrid[i])*jacxBpi_extended[i]
#        denom = denom + (abs(field_ext[i])**2 * bessel_factor[i] \
#            +abs(field_ext[i+1])**2 * bessel_factor[i+1])/2.*\
#            (zgrid[i+1]-zgrid[i])*jacxBpi_extended[i]
#
#    ave_quant = ave_quant/denom
#
#    return ave_quant

def eigenfunction_average_bessel(pars,geom,kperp,quantity,field,mass_ratio = 1,charge = 1, field_weighted = True, zstart = 0, zend = -1):

    ntot = pars['nx0']*pars['nz0']
    dz = float(2.0)/float(pars['nz0'])
    zgrid = np.pi*(np.arange(ntot)/float(ntot-1)*(2*pars['nx0']-dz)-pars['nx0'])

    jacxB_ext =  get_jacxB_extended(pars,zgrid,geom)

    kperp = kperp*np.sqrt(mass_ratio)/abs(charge)
    alpha = 2./3.
    bessel_factor = 1. / np.sqrt(1. + 2. * (kperp**2 + np.pi * alpha * kperp**4) /
                (1. + alpha * kperp**2))
    if no_bessel:
        bessel_factor[:] = 1.0

    show_plots = False
    if show_plots:
        plt.plot(bessel_factor)
        plt.title('bessel_factor')
        plt.show()
        plt.plot(jacxB_ext)
        plt.title('jacxB_ext')
        plt.show()
        np.savetxt('bessel_jacxB.dat',np.column_stack((zgrid,bessel_factor,jacxB_ext)))
    sum_ddz = 0
    denom = 0
    if zend == -1:
        zend = len(zgrid) - 1
    if field_weighted:
        weight = field
    else:
        weight = np.ones_like(field)

    sum_ddz = 0.5*integrate.simps(quantity[zstart:zend] *abs(weight[zstart:zend])**2* bessel_factor[zstart:zend]*jacxB_ext[zstart:zend],zgrid[zstart:zend])
    denom =  0.5*integrate.simps(abs(field[zstart:zend])**2 * bessel_factor[zstart:zend] *jacxB_ext[zstart:zend],zgrid[zstart:zend])
    #print('sum_ddz',sum_ddz)
    #print('denom',denom)
    avg_quantity = sum_ddz/denom
    return avg_quantity

par = Parameters()
par.Read_Pars('parameters'+suffix)
pars = par.pardict
print(pars)
print(pars['n_spec'])
for i in range(pars['n_spec']):
    if pars['charge'+str(i+1)] == -1:
        enum = i+1
        ename = pars['name'+str(i+1)]

print('electron species name:', ename)
print('electron species number:', enum)

omt = pars['omt'+str(enum)]
omn = pars['omn'+str(enum)]
mass_ratio = pars['mass'+str(enum)]

N=int(np.floor(pars['nx0']/2))+1
#f=numpy.loadtxt('fluxspectrae'+suffix+'.dat')
#ky1=[]
#Qes=[]
#for i in range(N,len(f)):
#    ky1.append(f[i][0])
#    Qes.append(f[i][2])

#parlin = Parameters()
#parlin.Read_Pars('scanfiles'+sfsuffix+'/parameters')
#parslin = parlin.pardict

kylin = []
kperp_arr = []
gamma = []
kx_center = []
chi_mixl = []
chi_ml_max = 0
Q_over_G = []

first_time = True
for i in range(200):   #large enough range to cover all possible in linear scan
    lin_suffix = '0000'+str(i)
    lin_suffix = lin_suffix[-4:]
    parlin_path = 'scanfiles'+sfsuffix+'/parameters_'+lin_suffix
    izstart = 0
    izend = -1
    if exists(parlin_path):
        parlin = Parameters()
        parlin.Read_Pars(parlin_path)
        parslin = parlin.pardict
        if parslin['n_spec'] == 1:
            time,nrg = get_nrg0('_'+lin_suffix,nspec=1,path='scanfiles'+sfsuffix+'/')
        elif parslin['n_spec'] == 2:
            time,nrg1,nrg2 = get_nrg0('_'+lin_suffix,nspec=2,path='scanfiles'+sfsuffix+'/')
            if enum == 1:
                nrg = nrg1
            else:
                nrg = nrg2
        elif parslin['n_spec'] == 3:
            time,nrg1,nrg2,nrg3 = get_nrg0('_'+lin_suffix,nspec=3,path='scanfiles'+sfsuffix+'/')
            if enum == 1:
                nrg = nrg1
            elif enum == 2:
                nrg = nrg2
            else:
                nrg = nrg3
        print('shape of nrg',np.shape(nrg))
        Q_over_G.append(nrg[-1,6]/nrg[-1,4])
        field = fieldfile('scanfiles'+sfsuffix+'/field_'+lin_suffix,parslin)
        time = np.array(field.tfld)
        itime = len(time)-1
        field.set_time(time[itime])
        phi = field.phi()
        kylin.append(parslin['kymin'])
        geomfile = 'scanfiles'+sfsuffix+'/'+parslin['magn_geometry'][1:-1]+'_'+lin_suffix
        gpars,geom = read_geometry_local(geomfile)
        #print(geom.keys())
        kperp = get_kperp(parslin,geom)
        if 'kx_center' in parslin:
            kx_center.append(parslin['kx_center'])
        else:
            kx_center.append(0)
        zgrid_ext,phi_ext = construct_extended_ballooning(parslin,phi)
        if first_time:
            first_time = False
            if select_zrange:
                print("Select z range now.")
                print("Present z range (in units of pi):",zgrid_ext[0]/np.pi,zgrid_ext[-1]/np.pi)
                zrange = float(input("Enter a value for z to define the range (units of pi):"))
                izstart = np.argmin(abs(zgrid_ext-(-zrange*np.pi)))            
                izend = np.argmin(abs(zgrid_ext-(zrange*np.pi)))            
                print("Integrating from z/pi = ",-zrange," to z/pi = ",zrange)
                print("izstart,izend",izstart,izend)
                print("len(zgrid_ext)",len(zgrid_ext))
                dummy = input("press any key")
            
        kperp_avg = eigenfunction_average_bessel(parslin,geom,kperp,kperp,phi_ext,mass_ratio=mass_ratio,zstart = izstart,zend = izend)
        kperp_arr.append(kperp_avg)
        f = open('scanfiles'+sfsuffix+'/'+'omega_'+lin_suffix,'r')
        data = f.read().split()
        #print('data',data)
        gamma.append(float(data[1]))
        chi_mixl.append(gamma[-1]/kperp_arr[-1]**2)
        if kx_center[-1] ==0 and chi_mixl[-1] > chi_ml_max:
            chi_ml_max = chi_mixl[-1]
        #print('ky,gamma,kx_center,kperp',kylin[-1],gamma[-1],kx_center[-1],kperp_arr[-1])
        #plt.plot(zgrid_ext,abs(phi_ext)/np.max(abs(phi_ext))*np.max(kperp),label='abs(phi)')
        #plt.plot(zgrid_ext,kperp,label='kperp')
        #plt.show()

ipeak = np.argmax(np.array(chi_mixl)/np.array(Q_over_G))
Q_over_G_peak = Q_over_G[ipeak]


QQL = chi_ml_max*(omt/(1+omn))**2*0.877*omt
GQL = QQL/Q_over_G_peak*(omt/(1+omn))
print("Q: Prediction of quasilinear model",QQL)
print("G: Prediction of quasilinear model",GQL)
    




