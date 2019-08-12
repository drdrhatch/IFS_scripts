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
from get_nrg import *
from finite_differences import *
from interp import *

parser=op.OptionParser(description='')
parser.add_option('--show_plots','-p', action='store_true',dest = 'show_plots', help = 'Eigenfunction average kperp, omega', default=False)
##parser.add_option("-a", "--avg", action="store_true", dest="avg_flag", default=False, help="Eigenfunction average kperp, omega")

options,args=parser.parse_args()
suffix = args[0]
avg = True
smooth_field = True
show_plots = options.show_plots
if not suffix =='.dat':
   suffix = '_'+suffix

def eigenfunction_average_bessel(z_grid,jacobian,kperp,omega_d,field,name):

    alpha = 2./3.
    bessel_factor = 1. / np.sqrt(1. + 2. * (kperp**2 + np.pi * alpha * kperp**4) /
                    (1. + alpha * kperp**2))

    ave_kperp2 = 0.
    denom = 0.
    for i in np.arange(len(field)-1):
        ave_kperp2 = ave_kperp2 + (kperp[i]**2*abs(field[i])**2 +\
            kperp[i+1]**2*abs(field[i+1])**2)/2.*\
            (z_grid[i+1]-z_grid[i])/jacobian[i]
        denom = denom + (abs(field[i])**2 +abs(field[i+1])**2)/2.*\
            (z_grid[i+1]-z_grid[i])/jacobian[i]
    ave_kperp2 = ave_kperp2/denom
    ave_kperp = np.sqrt(ave_kperp2)
    #print name + ' weighted k_perp^2 =', ave_kperp2
    print name + ' weighted k_perp =', ave_kperp

    ave_omegad = 0.
    denom = 0.
    for i in np.arange(len(field)-1):
        ave_omegad = ave_omegad + (omega_d[i]*abs(field[i])**2 * bessel_factor[i]+\
            omega_d[i+1]*abs(field[i+1])**2 * bessel_factor[i+1])/2.*\
            (z_grid[i+1]-z_grid[i])/jacobian[i]
        denom = denom + (abs(field[i])**2 * bessel_factor[i] \
            +abs(field[i+1])**2 * bessel_factor[i+1])/2.*\
            (z_grid[i+1]-z_grid[i])/jacobian[i]

    ave_omegad = ave_omegad/denom
    print name + ' weighted omega_d =', ave_omegad

    return ave_kperp, ave_omegad

def kz_from_dfielddz_bessel(kperp, zgrid, jacobian, field, plot, name):
    alpha = 2./3.
    bessel_factor = 1. / np.sqrt(1. + 2. * (kperp**2 + np.pi * alpha * kperp**4) /
                    (1. + alpha * kperp**2))
    if 1 == 0:
        plt.plot(zgrid, bessel_factor, label = 'bessel_factor')
        plt.legend()
        plt.show()

    dfielddz = np.empty(len(field),dtype='complex128')
    for i in range(len(field)-1):
        dfielddz[i] = (field[i+1]-field[i])/\
            (zgrid[i+1]-zgrid[i])*jacobian[i]
    if plot:
        plt.plot(zgrid[:-1], np.abs(dfielddz[:-1]), label = 'abs d'+name+'/dz')
        plt.plot(zgrid[:-1], np.real(dfielddz[:-1]), label = 'real d'+name+'/dz')
        plt.plot(zgrid[:-1], np.imag(dfielddz[:-1]), label = 'imag d'+name+'/dz')
        plt.legend()
        plt.xlabel('z')
        plt.show()
    sum_ddz = 0.
    denom = 0.
    
    zstart = 5
    zend = len(zgrid)-5
    #startInd = np.argmin(abs(zgrid - zstart))
    #endInd = np.argmin(abs(zgrid - zend))
    #for i in (startInd, endInd + 1):
    for i in range(zstart,zend):
        sum_ddz = sum_ddz + 0.5*(abs(dfielddz[i])**2 * bessel_factor[i] \
                  + abs(dfielddz[i+1])**2 * bessel_factor[i+1])*\
                  (zgrid[i+1]-zgrid[i])/jacobian[i]
        denom = denom + 0.5*(abs(field[i])**2 * bessel_factor[i] \
                + abs(field[i+1])**2* bessel_factor[i+1])*\
                (zgrid[i+1]-zgrid[i])/jacobian[i]
    ave_kz = np.sqrt(sum_ddz/denom)
    print name + ' averaged kz = ', ave_kz
    print 'Input to SKiM kz = ', ave_kz
    return ave_kz

pars = init_read_parameters_file(suffix)
if pars['comp_type'][1:-1] == 'EV':
    omfile = 'eigenvalues'+suffix
else:
    omfile = 'omega'+suffix

om = np.genfromtxt(omfile)

if 'kx_center' in pars:
    kx_center = pars['kx_center']
else:
    kx_center = 0.0

if pars['n_spec'] == 1:
    time,nrg = get_nrg0(suffix,nspec = 1)
elif pars['n_spec'] == 2:
    time,nrg1,nrg2 = get_nrg0(suffix,nspec = 2)
    if 'e' in pars['name2']:
        nrg = nrg2
    elif 'e' in pars['name1']:
        nrg = nrg1
    else:
        print "Error! Can't fine electron species number in parameters file."
        sys.exit()
else:
    print "Error! This script can only handle 2 species at this time."
    sys.exit()

if pars['comp_type'][1:-1] == 'EV':
    EVrun = True
    NEV = len(om[:,0])
    gamma = om[:,0]
    omega = om[:,1]
else:
    EVrun = False
    NEV = 1
    gamma = [om[1]]
    omega = [om[2]]
#gamma.append(om[1])
#omega.append(om[2])

xlocal = True
geom_type, geom_pars, geom_coeff = init_read_geometry_file(suffix,pars)

zgrid, jacobian = reconstruct_zgrid(geom_coeff, pars, True, False,0.)
kperp, omd_curv, omd_gradB = calc_kperp_omd(geom_type,geom_coeff,pars,True,False)
ggxx,ggxy,ggxz,ggyy,ggyz,ggzz,gdBdx,gdBdy,gdBdz,gBfield,gjacobian,gl_R,gl_phi,gl_z,gl_dxdR,gl_dxdZ = read_geom_coeff_raw(geom_type,geom_coeff,True)

f = open('geom_coeff' + suffix + '.txt','w')
f.write('#1.zgrid 2.gxx 3.gxy 4.gyy 5.gBfield\n')
np.savetxt(f, np.column_stack((zgrid, ggxx, ggxy, ggyy, gBfield)))
f.close()

zgrid, jacobian = reconstruct_zgrid(geom_coeff, pars, False, False, 0.)
kperp, omd_curv, omd_gradB = calc_kperp_omd(geom_type,geom_coeff,pars,False,False)
f = open('kperp_vector3' + suffix + '.txt','w')
f.write('#1.zgrid 2.kperp 3.omd_curv 4.omd_gradB\n')
np.savetxt(f, np.column_stack((zgrid, kperp, omd_curv, omd_gradB)))
f.close()

for i in range(NEV):

   if EVrun:
       itime = i+1
   else:
       itime = -1

   print("\n################ Eigenvalue "+str(i+1)+" ################\n" )
   Qn2 = nrg[itime-1,6]/nrg[itime-1,0]

   #phi, apar = eigenfunctions_from_field_file(pars,suffix,False,False,setTime = itime, smooth_field = True)
    #upar,deln,tpar,tperp,qpar,qperp = moments_from_mom_file(pars,suffix,False,False,setTime=-1)
   #abs_phi = abs(phi)
   #imax = np.argmax(abs_phi)
   #phi00 = phi[imax]
   #print('phi00 = ', phi00)
   #phi = phi / phi00
   #phi_plot = phi.copy()
   #phi = np.real(phi)

   #if 1 == 0: #using j_par to weight k_perp
   #    upar_smoothed01 = field_smoother(upar)
   #    tmp0, tmp1 =  eigenfunction_average_bessel(zgrid, jacobian, kperp, omd_curv, upar,'upar')
   #    upar_avg_kperp, tmp =  eigenfunction_average_bessel(zgrid, jacobian, kperp, omd_curv, upar_smoothed01,'upar')
   #    f = open('upar_avg_kperp3_ev'+str(i+1) + suffix,'w')
   #    f.write('#1.ky 2.upar_avg_kperp 3. ratio \n')
   #    ky = float(pars['kymin'])
   #    np.savetxt(f, np.column_stack((ky, upar_avg_kperp, upar_avg_kperp/ky)), fmt='%4.4f', delimiter = "  ")
   #    f.close()
   if avg:
       # compute kz, kperp, omd using phi
       phi, apar = eigenfunctions_from_field_file(pars,suffix,False,False,setTime = itime)
       abs_phi = abs(phi)
       imax = np.argmax(abs_phi)
       phi00 = phi[imax]
       print('phi00 = ', phi00)
       phi = phi / phi00

       #Calculate number of zeros
       nzfine = len(zgrid)*10
       dzfine = (zgrid[-1]-zgrid[0])/float(nzfine)
       zgrid_fine = np.arange(zgrid[0],zgrid[-1],dzfine)
       phi_fine_i = interp(zgrid,np.imag(phi),zgrid_fine)
       phi_fine_r = interp(zgrid,np.real(phi),zgrid_fine)
       phi_fine = phi_fine_r + 1.0J*phi_fine_i
       dphidz = fd_d1_o4(abs(phi_fine),zgrid_fine)
       d2phidz2 = fd_d1_o4(dphidz,zgrid_fine)
       nzeros = 0
       zero_loc = np.empty(0)
       #plt.plot(zgrid_fine,abs(phi_fine))
       #plt.plot(zgrid,abs(phi))
       for k in range(len(zgrid_fine)):
           if abs(dphidz[k]) < abs(dphidz[k-1]) and abs(dphidz[k]) < abs(dphidz[k+1]) and abs(dphidz[k])/abs(phi_fine)[k] < 0.2 and abs(phi_fine[k]) > 0.2 and d2phidz2[k] < 0 :
               zero_loc = np.append(zero_loc,zgrid_fine[k])
               nzeros += 1
               #plt.axvline(x=zgrid_fine[k],color = 'gray')
       #plt.title('nzeros = '+str(nzeros))
       #plt.show()

       #plt.plot(zgrid_fine,abs(dphidz),'x-')
       #plt.title('nzeros = '+str(nzeros))
       #plt.show()

       #plt.plot(zgrid_fine,d2phidz2)
       #plt.title('nzeros = '+str(nzeros))
       #plt.show()


       phi_plot = phi.copy()
       phi = np.real(phi)
       print('phi = ', phi[0])
       phi_smoothed01 = field_smoother(phi)
       print('phi = ', phi[0])
       phi_smoothed02 = field_smoother(phi_smoothed01)
   
       ave_kz = kz_from_dfielddz_bessel(kperp, zgrid,jacobian, phi, False, 'phi')
       ave_kperp, ave_omd = eigenfunction_average_bessel(zgrid, jacobian, kperp, omd_curv, phi,'phi')
       ave_sq_int, ave_int_sq = eigenfunction_squared(zgrid, jacobian, phi)
   
       ave_kz_smoothed01  = kz_from_dfielddz_bessel(kperp, zgrid, jacobian, phi_smoothed01, False, 'phi smoothed once' )
       ave_kperp_smoothed01, ave_omd_smoothed01 =  eigenfunction_average_bessel(zgrid, jacobian, kperp, omd_curv, phi_smoothed01,'phi smoothed once')
       ave_sq_int01, ave_int_sq01 = eigenfunction_squared(zgrid, jacobian, phi_smoothed01)
   
       ave_kz_smoothed02  = kz_from_dfielddz_bessel(kperp, zgrid, jacobian, phi_smoothed02, False, 'phi smoothed twice' )
       ave_kperp_smoothed02, ave_omd_smoothed02 =  eigenfunction_average_bessel(zgrid, jacobian, kperp, omd_curv, phi_smoothed02, 'phi smoothed twice')
       ave_sq_int02, ave_int_sq02 = eigenfunction_squared(zgrid, jacobian, phi_smoothed02)
       
       f = open('averaged3_ev' +str(i+1) + suffix , 'w')
       #sys.stdout = f
       f.write( '')
       f.write( '')
       f.write( 'Using phi from GENE run')
       f.write( 'kz = '+str( ave_kz))
       f.write( 'kperp ='+str( ave_kperp))
       f.write( 'omd ='+str( ave_omd))
       f.write( 'gamma/kperp**2 ='+str( float(gamma[i]) / ave_kperp**2))
       f.write( 'kz**2 ='+str( ave_kz**2))
       f.write( 'omd * abs(omega) ='+str( ave_omd * np.sqrt(om[2]**2 + om[1]**2)))
       f.write( 'kz * (int phi)**2 / int(phi**2)/ sqrt(2pi) = '+str( ave_kz/np.sqrt(2. * np.pi) *ave_int_sq / ave_sq_int))
       f.write( 'Input to SKiM')
       f.write( 'kperp ='+str( ave_kperp / float(pars['kymin'])))
       f.write( 'omd ='+str( ave_omd / float(pars['kymin'])))
       f.write( '')
       f.write( '')
       f.write( 'Using once smoothed phi from GENE run')
       f.write( 'kz = '+str( ave_kz_smoothed01))
       f.write( 'kperp ='+str( ave_kperp_smoothed01))
       f.write( 'omd ='+str( ave_omd_smoothed01))
       f.write( 'gamma/kperp**2 ='+str( float(gamma[i]) / ave_kperp_smoothed01**2))
       f.write( 'kz**2 ='+str( ave_kz_smoothed01**2))
       f.write( 'omd * abs(omega) ='+str( ave_omd_smoothed01 * np.sqrt(omega[i]**2 + gamma[i]**2)))
       f.write( 'kz * (int phi)**2 / int(phi**2)/ sqrt(2pi) = '+str( ave_kz_smoothed01/np.sqrt(2. * np.pi) *ave_int_sq01 / ave_sq_int01))
       f.write( 'Input to SKiM')
       f.write( 'kperp ='+str( ave_kperp_smoothed01 / float(pars['kymin'])))
       f.write( 'omd ='+str( ave_omd_smoothed01 / float(pars['kymin'])))
       f.write( '')
       f.write( '')
       f.write( 'Using twice smoothed phi from GENE run')
       f.write( 'kz = '+str( ave_kz_smoothed02))
       f.write( 'kperp ='+str( ave_kperp_smoothed02))
       f.write( 'omd ='+str( ave_omd_smoothed02))
       f.write( 'gamma/kperp**2 ='+str( float(gamma[i]) / ave_kperp_smoothed02**2))
       f.write( 'kz**2 ='+str( ave_kz_smoothed02**2))
       f.write( 'omd * abs(omega) ='+str( ave_omd_smoothed02 * np.sqrt(om[2]**2 + om[1]**2)))
       f.write( 'kz * (int phi)**2 / int(phi**2)/ sqrt(2pi) = '+str( ave_kz_smoothed02/np.sqrt(2. * np.pi) *ave_int_sq02 / ave_sq_int02))
       f.write( 'Input to SKiM')
       f.write( 'kperp ='+str( ave_kperp_smoothed02 / float(pars['kymin'])))
       f.write( 'omd ='+str( ave_omd_smoothed02 / float(pars['kymin'])))
       f.close()
   
       f = open('mode_info_'+str(i+1) + suffix,'w')
       f.write('#1.ky 2.kz 3.kperp 4.omd 5.gamma 6.omega 7.gam/kperp^2 8.Q/n^2 9.num_zeros 10.kx_center \n')
       ky = float(pars['kymin'])
       np.savetxt(f, np.column_stack((ky, 
                  ave_kz_smoothed01, 
                  ave_kperp_smoothed01, 
                  ave_omd_smoothed01,
                  float(gamma[i]),
                  float(omega[i]),
                  float(gamma[i]) / ave_kperp_smoothed01**2,
                  float(Qn2),nzeros,kx_center)),fmt='%4.4f', delimiter = "  ")
       f.close()
   
       if 'kx_center' in pars:
           kx_center = float(pars['kx_center'])
       else:
           kx_center = 0.
       theta_0 = kx_center/np.pi/float(pars['shat'])/float(pars['kymin'])
       theta_0 = np.around(theta_0,decimals=1)
   
       if show_plots:
           plt.plot(zgrid,abs(phi_plot),label='abs phi',color='black')
           plt.plot(zgrid,np.real(phi_plot),label='real phi',color='blue')
           plt.plot(zgrid,np.imag(phi_plot),label='imag phi',color='red')
           plt.plot(zgrid,kperp,label='kperp',color='lawngreen')
           plt.plot(zgrid,omd_curv,label='omd curv',color='aqua')
          #plt.title('ky = '+str(pars['kymin'])+', kx_center = '+str(theta_0)+r'$\pi$')
           plt.xlabel('z(pi)')
           plt.grid()
           plt.legend()
           plt.show()

