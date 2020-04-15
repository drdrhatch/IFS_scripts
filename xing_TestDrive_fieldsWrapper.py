#!/usr/bin/env python
# -*- coding: utf-8 -*-

import optparse as op
import matplotlib.pyplot as plt
import numpy as np
from parIOWrapper import init_read_parameters_file
from geomWrapper import *
from fieldsWrapper import *
from momentsWrapper import *
import sys


parser=op.OptionParser(description='')
parser.add_option('--show_plots','-p', action='store',dest = 'show_plots', help = 'Eigenfunction average kperp, omega', default=False)
parser.add_option('--avg','-a', action='store',dest = 'avg', help = 'Eigenfunction average kperp, omega', default=False)
parser.add_option('--g','-g', action='store',dest = 'gauss', help = 'Average kperp, omega with Gaussian width', default=0.)
parser.add_option('--smooth_field','-s', action='store',dest = 'smooth_field', help = 'smooth field with phi_i = 0.5phi_i + 0.25(phi_i-1 + phi_i+1)', default=False)
parser.add_option('--c','-c', action='store',dest = 'central', help = 'Only look at central', default=0.)
options,args=parser.parse_args()
suffix = args[0]
avg = options.avg
gauss = float(options.gauss)
show_plots = options.show_plots
smooth_field = options.smooth_field
central = options.central
print(('gauss =', gauss))
if not suffix =='.dat':
   suffix = '_'+suffix

pars = init_read_parameters_file(suffix)

om = np.genfromtxt('omega'+suffix)
#gamma.append(om[1])
#omega.append(om[2])

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
    zgrid, jacobian = reconstruct_zgrid(geom_coeff, pars, True, False,0.)
    kperp, omd_curv, omd_gradB = calc_kperp_omd(geom_type,geom_coeff,pars,True,False)
    ggxx,ggxy,ggxz,ggyy,ggyz,ggzz,gdBdx,gdBdy,gdBdz,gBfield,gjacobian,gl_R,gl_phi,gl_z,gl_dxdR,gl_dxdZ = read_geom_coeff_raw(geom_type,geom_coeff,True)

    f = open('geom_coeff' + suffix + '.txt','w')
    f.write('#1.zgrid 2.gxx 3.gxy 4.gyy \n')
    np.savetxt(f, np.column_stack((zgrid, ggxx, ggxy, ggyy)))
    f.close()

    if not central:
        zgrid, jacobian = reconstruct_zgrid(geom_coeff, pars, False, False, 0.)
        kperp, omd_curv, omd_gradB = calc_kperp_omd(geom_type,geom_coeff,pars,False,False)
        f = open('kperp_vector' + suffix + '.txt','w')
        f.write('#1.zgrid 2.kperp 3.omd_curv 4.omd_gradB\n')
        np.savetxt(f, np.column_stack((zgrid, kperp, omd_curv, omd_gradB)))
        f.close()
#    phi = np.exp(-(zgrid / gauss)**2)
    #if not gauss == 0.:
    #    phi = np.exp(-(zgrid / gauss)**2)
    #else:
    phi, apar = eigenfunctions_from_field_file(pars,suffix,False,False,-1, smooth_field)
    #upar,deln,tpar,tperp,qpar,qperp = moments_from_mom_file(pars,suffix,False,False,setTime=-1)


    if 1 == 0: #using j_par to weight k_perp
        upar_smoothed01 = field_smoother(upar)
        tmp0, tmp1 =  eigenfunction_average(zgrid, jacobian, kperp, omd_curv, upar,'upar')
        upar_avg_kperp, tmp =  eigenfunction_average(zgrid, jacobian, kperp, omd_curv, upar_smoothed01,'upar')
        f = open('upar_avg_kperp' + suffix,'w')
        f.write('#1.ky 2.upar_avg_kperp 3. ratio \n')
        ky = float(pars['kymin'])
        np.savetxt(f, np.column_stack((ky, upar_avg_kperp, upar_avg_kperp/ky)), fmt='%4.4f', delimiter = "  ")
        f.close()
    if avg:
        # compute kz, kperp, omd using phi
        if not central:
            phi, apar = eigenfunctions_from_field_file(pars,suffix,False,False,-1)
        else:
            phi, apar = eigenfunctions_from_field_file(pars,suffix,True,False,-1)
        phi_smoothed01 = field_smoother(phi)
        phi_smoothed02 = field_smoother(phi_smoothed01)

        if not central:
            ave_kz, zstart, zend = kz_from_dfielddz(zgrid,jacobian, phi, True, 'phi')
        else:
            ave_kz, zstart, zend = kz_from_dfielddz(zgrid,jacobian, phi, False, 'phi', zstart = -.95, zend = .95)
        ave_kperp, ave_omd =  eigenfunction_average(zgrid, jacobian, kperp, omd_curv, phi,'phi')
        ave_sq_int, ave_int_sq = eigenfunction_squared(zgrid, jacobian, phi)

        ave_kz_smoothed01, zstart01, zend01  = kz_from_dfielddz(zgrid, jacobian, phi_smoothed01, False, 'phi smoothed once', zstart = zstart, zend = zend)
        ave_kperp_smoothed01, ave_omd_smoothed01 =  eigenfunction_average(zgrid, jacobian, kperp, omd_curv, phi_smoothed01,'phi smoothed once')
        ave_sq_int01, ave_int_sq01 = eigenfunction_squared(zgrid, jacobian, phi_smoothed01)

        ave_kz_smoothed02, zstart02, zend02  = kz_from_dfielddz(zgrid, jacobian, phi_smoothed02, False, 'phi smoothed twice', zstart = zstart, zend = zend)
        ave_kperp_smoothed02, ave_omd_smoothed02 =  eigenfunction_average(zgrid, jacobian, kperp, omd_curv, phi_smoothed02, 'phi smoothed twice')
        ave_sq_int02, ave_int_sq02 = eigenfunction_squared(zgrid, jacobian, phi_smoothed02)
        
        f = open('averaged' + suffix + str(gauss).replace(".",""), 'w')
        sys.stdout = f
        print('Using phi from GENE run')
        print('kz = ', ave_kz)
        print('kperp =', ave_kperp)
        print('omd =', ave_omd)
        print('gamma/kperp**2 =', float(om[1]) / ave_kperp**2)
        print('')
        print('')
        print('Input to SKiM')
        print('kperp =', ave_kperp / float(pars['kymin']))
        print('omd =', ave_omd / float(pars['kymin']))
        print('')
        print('')
        print('kz * (int phi)**2 / int(phi**2)/ sqrt(2pi) = ', ave_kz/np.sqrt(2. * np.pi) *ave_int_sq / ave_sq_int)
        print('Using once smoothed phi from GENE run')
        print('kz = ', ave_kz_smoothed01)
        print('kperp =', ave_kperp_smoothed01)
        print('omd =', ave_omd_smoothed01)
        print('gamma/kperp**2 =', float(om[1]) / ave_kperp_smoothed01**2)
        print('')
        print('')
        print('Input to SKiM')
        print('kperp =', ave_kperp_smoothed01 / float(pars['kymin']))
        print('omd =', ave_omd_smoothed01 / float(pars['kymin']))
        print('')
        print('')
        print('kz * (int phi)**2 / int(phi**2)/ sqrt(2pi) = ', ave_kz_smoothed01/np.sqrt(2. * np.pi) *ave_int_sq01 / ave_sq_int01)
        print('Using twice smoothed phi from GENE run')
        print('kz = ', ave_kz_smoothed02)
        print('kperp =', ave_kperp_smoothed02)
        print('omd =', ave_omd_smoothed02)
        print('gamma/kperp**2 =', float(om[1]) / ave_kperp_smoothed02**2)
        print('')
        print('')
        print('Input to SKiM')
        print('kperp =', ave_kperp_smoothed02 / float(pars['kymin']))
        print('omd =', ave_omd_smoothed02 / float(pars['kymin']))
        print('')
        print('')
        print('kz * (int phi)**2 / int(phi**2)/ sqrt(2pi) = ', ave_kz_smoothed02/np.sqrt(2. * np.pi) *ave_int_sq02 / ave_sq_int02)
        f.close()


        f = open('avg_kz_kperp_omd' + suffix,'w')
        f.write('#1.ky 2.kz 3.kperp 4.omd\n')
        ky = float(pars['kymin'])
        np.savetxt(f, np.column_stack((ky, ave_kz_smoothed02, ave_kperp_smoothed02, ave_omd_smoothed02)), fmt='%4.4f', delimiter = "  ")
        f.close()

    if 'kx_center' in pars:
        kx_center = float(pars['kx_center'])
    else:
        kx_center = 0.
    theta_0 = kx_center/np.pi/float(pars['shat'])/float(pars['kymin'])
    theta_0 = np.around(theta_0,decimals=1)

    if show_plots:
        plt.plot(zgrid,abs(phi),label='abs phi',color='black')
        plt.plot(zgrid,np.real(phi),label='real phi',color='blue')
        plt.plot(zgrid,np.imag(phi),label='imag phi',color='red')
        plt.plot(zgrid,kperp,label='kperp',color='lawngreen')
        plt.plot(zgrid,omd_curv,label='omd curv',color='aqua')
       #plt.title('ky = '+str(pars['kymin'])+', kx_center = '+str(theta_0)+r'$\pi$')
        plt.xlabel('z(pi)')
        plt.grid()
        plt.legend()
        plt.show()
    if 1 == 0:
        plt.plot(zgrid,abs(apar),label='abs apar',color='black')
        plt.plot(zgrid,np.real(apar),label='real apar',color='blue')
        plt.plot(zgrid,np.imag(apar),label='imag apar',color='red')
        plt.plot(zgrid,kperp,label='kperp',color='lawngreen')
        plt.plot(zgrid,omd_curv,label='omd curv',color='aqua')
        plt.title('ky = '+str(pars['kymin'])+', kx_center = '+str(theta_0)+r'$\pi$')
        plt.grid()
        plt.legend()
        plt.show()

    if 1 == 0:
        print('calculate average kz using d field/ dz')
        ave_kz = kz_from_dfielddz(zgrid,jacobian, phi, True, 'phi')
        f = open('averaged' + suffix + str(gauss).replace(".",""), 'w')
        sys.stdout = f
        if gauss == 0.:
            print('Using field from GENE run')
        else:
            print('Gaussian width z = ', gauss)
        print('Input to SKiM kz = ', ave_kz)
        ave_kperp, ave_omd =  eigenfunction_average(zgrid, jacobian,kperp,omd_curv,phi,'phi')
        print('Input to SKiM kperp =', ave_kperp / pars['kymin'])
        print('Input to SKiM omd =', ave_omd / pars['kymin'])
#        print 'calculate average kz using FT'
#        phi1 = list(phi)
#        field_kz, kz_grid = fourierTrans(pars, zgrid, jacobian, phi1, True, 'phi')
