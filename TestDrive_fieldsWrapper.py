#!/usr/bin/env python
# -*- coding: utf-8 -*-

import optparse as op
import matplotlib.pyplot as plt
import numpy as np
from parIOWrapper import init_read_parameters_file
from geomWrapper import *
from fieldsWrapper import *
from momentsWrapper import *

parser=op.OptionParser(description='')
parser.add_option('--avg','-a', action='store',dest = 'avg', help = 'Eigenfunction average kperp, omega', default=False)
options,args=parser.parse_args()
suffix = args[0]
avg = options.avg

if not suffix =='.dat':
   suffix = '_'+suffix

pars = init_read_parameters_file(suffix)

lilo = False
xlocal = True

if lilo:
    plot = True
    phi, apar = LILO_eigenfunctions_from_field_file(pars,suffix,plot,setTime=-1)
    upar, dens = LILO_moments_from_mom_file(pars,suffix,plot,setTime=-1)
elif xlocal:

    geom_type, geom_pars, geom_coeff = init_read_geometry_file(suffix,pars)

    #phi = np.zeros(pars['nx0']*pars['nz0'],dtype='complex128')
    #apar = np.zeros(pars['nx0']*pars['nz0'],dtype='complex128')
    print('Eigenfunction at outboard midplane')
    zgrid, jacobian = reconstruct_zgrid(geom_coeff, pars, False, False,0.)
    kperp, omd_curv, omd_gradB = calc_kperp_omd(geom_type,geom_coeff,pars,False,False)
    print('Eigenfunction total')
    phi, apar = eigenfunctions_from_field_file(pars,suffix,False,False,-1)
    if avg:
        ave_kperp, ave_omd =  eigenfunction_average(zgrid, jacobian,kperp,omd_curv,phi,'phi')
        #ave_kperp, ave_omd =  eigenfunction_average(zgrid, jacobian,kperp,omd_curv,apar,'apar')
        print('Input to SKiM kperp =', ave_kperp / float(pars['kymin']))
        print('Input to SKiM omd =', ave_omd / pars['kymin'])
        print('calculate average kz using FT')
        phi1 = list(phi)
        apar1 = list(apar)
        #field_kz, kz_grid = fourierTrans(pars, zgrid, jacobian, apar1, True, 'apar')
        print('calculate average kz using d field/ dz')
#        ave_kz = kz_from_dfielddz(zgrid,jacobian,apar,True, 'apar')
        ave_kz = kz_from_dfielddz(zgrid,jacobian,phi,True, 'phi')
        print('ave_kz', ave_kz)

    if 'kx_center' in pars:
        kx_center = float(pars['kx_center'])
    else:
        kx_center = 0.
    theta_0 = kx_center/np.pi/float(pars['shat'])\
              /float(pars['kymin'])
    theta_0 = np.around(theta_0,decimals=3)

    if 1 == 1:
        plt.plot(zgrid,abs(phi),label='abs phi',color='black')
        plt.plot(zgrid,np.real(phi),label='real phi',color='blue')
        plt.plot(zgrid,np.imag(phi),label='imag phi',color='red')
        plt.plot(zgrid,kperp,label='kperp',color='lawngreen')
        plt.plot(zgrid,omd_curv,label='omd curv',color='aqua')
        plt.title('ky = '+str(pars['kymin'])+', kx_center = '+str(theta_0)+r'$\pi$')
        plt.grid()
        plt.legend()
        plt.show()
    if 1 == 1:
        plt.plot(zgrid,abs(apar),label='abs apar',color='black')
        plt.plot(zgrid,np.real(apar),label='real apar',color='blue')
        plt.plot(zgrid,np.imag(apar),label='imag apar',color='red')
        plt.plot(zgrid,kperp,label='kperp',color='lawngreen')
        plt.plot(zgrid,omd_curv,label='omd curv',color='aqua')
        plt.title('ky = '+str(pars['kymin'])+', kx_center = '+str(theta_0)+r'$\pi$')
        plt.grid()
        plt.legend()
        plt.show()

    upar,deln,tpar,tperp,qpar,qperp = moments_from_mom_file(pars,suffix,False,False,setTime=-1)

#    if avg:
#        ave_kz = kz_from_dfielddz(zgrid,jacobian,upar,True, 'upar')
#        print 'ave_kz', ave_kz

    if 1 == 1:
        plt.plot(zgrid,abs(upar),label='abs upar',color='black')
        plt.plot(zgrid,np.real(upar),label='real upar',color='blue')
        plt.plot(zgrid,np.imag(upar),label='imag upar',color='red')
        plt.plot(zgrid,kperp,label='kperp',color='lawngreen')
        plt.plot(zgrid,omd_curv,label='omd curv',color='aqua')
        plt.title('ky = '+str(pars['kymin'])+', kx_center = '+str(theta_0)+r'$\pi$')
        plt.grid()
        plt.legend()
        plt.show()
