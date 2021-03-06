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
parser.add_option('--nzeros','-n', action='store_true',dest = 'nzeros', help = 'Calculates the number of bumps in the eigenfunction.', default=False)
##parser.add_option("-a", "--avg", action="store_true", dest="avg_flag", default=False, help="Eigenfunction average kperp, omega")

options,args=parser.parse_args()
suffix = args[0]
avg = True
smooth_field = True
show_plots = options.show_plots
nzeros = options.nzeros
if not suffix =='.dat':
   suffix = '_'+suffix

def eigenfunction_width(z_grid,jacobian,field,name):

    field2_integral = 0.
    dz = zgrid[1]-zgrid[0]
    denom = np.max(abs(field)**2)
    #for i in np.arange(len(field)-1):
    #    field2_integral +=  dz*abs(field[i])**2/jacobian[i] 

    for i in np.arange(len(field)-1):
        field2_integral +=  (abs(field[i])**2 +\
            abs(field[i+1])**2)/2.*\
            (zgrid[i+1]-zgrid[i])/jacobian[i]

    width = field2_integral / denom
    #print (name + ' width =', width)
    return width

def eigenfunction_average_bessel(z_grid,jacobian,kperp,omega_d,field,name):

    alpha = 2./3.
    bessel_factor = 1. / np.sqrt(1. + 2. * (kperp**2 + np.pi * alpha * kperp**4) /
                    (1. + alpha * kperp**2))

    #plt.plot(bessel_factor)
    #plt.title('bessel')
    #plt.show()
    #plt.plot(1/jacobian)
    #plt.title('J')
    #plt.show()
    #plt.plot(kperp)
    #plt.title('kperp')
    #plt.show()
    #plt.plot(omega_d)
    #plt.title('omd')
    #plt.show()
    #plt.plot(np.real(field))
    #plt.plot(np.imag(field))
    #plt.title('field')
    #plt.show()
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
    #print (name + ' weighted k_perp =', ave_kperp)

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
    #print (name + ' weighted omega_d =', ave_omegad)

    return ave_kperp, ave_omegad

def eigenfunction_average_kx(z_grid,jacobian,kx_ext,field,name):

    ave_kx2 = 0.
    denom = 0.
    for i in np.arange(len(field)-1):
        ave_kx2 = ave_kx2 + (kx_ext[i]**2*abs(field[i])**2 +\
            kx_ext[i+1]**2*abs(field[i+1])**2)/2.*\
            (z_grid[i+1]-z_grid[i])/jacobian[i]
        denom = denom + (abs(field[i])**2 +abs(field[i+1])**2)/2.*\
            (z_grid[i+1]-z_grid[i])/jacobian[i]
    ave_kx2 = ave_kx2/denom
    ave_kx = np.sqrt(ave_kx2)
    #print name + ' weighted k_perp^2 =', ave_kperp2
    #print (name + ' weighted kx=', ave_kx)

    return ave_kx

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
    #print (name + ' averaged kz = ', ave_kz)
    #print ('Input to SKiM kz = ', ave_kz)
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
        print ("Error! Can't find electron species number in parameters file.")
        sys.exit()
elif pars['n_spec'] == 3:
    time,nrg1,nrg2,nrg3 = get_nrg0(suffix,nspec = 3)
    if 'e' in pars['name2']:
        nrg = nrg2
    elif 'e' in pars['name1']:
        nrg = nrg1
    elif 'e' in pars['name3']:
        nrg = nrg3
    else:
        print ("Error! Can't find electron species number in parameters file.")
        sys.exit()

#else:
#    print ("Error! This script can only handle 2 species at this time.")
#    sys.exit()

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

#zgrid, jacobian = reconstruct_zgrid(geom_coeff, pars, True, False,0.)
#kperp, omd_curv, omd_gradB = calc_kperp_omd(geom_type,geom_coeff,pars,True,False)
ggxx,ggxy,ggxz,ggyy,ggyz,ggzz,gdBdx,gdBdy,gdBdz,gBfield,gjacobian,gl_R,gl_phi,gl_z,gl_dxdR,gl_dxdZ = read_geom_coeff_raw(geom_type,geom_coeff,True)

#f = open('geom_coeff' + suffix + '.txt','w')
#f.write('#1.zgrid 2.gxx 3.gxy 4.gyy 5.gBfield\n')
#np.savetxt(f, np.column_stack((zgrid, ggxx, ggxy, ggyy, gBfield)))
#f.close()

zgrid, jacobian = reconstruct_zgrid(geom_coeff, pars, False, False, 0.)
kperp, omd_curv, omd_gradB = calc_kperp_omd(geom_type,geom_coeff,pars,False,False)
kx_ext = calc_kx_extended(pars,False)
#f = open('kperp_vector3' + suffix + '.txt','w')
#f.write('#1.zgrid 2.kperp 3.omd_curv 4.omd_gradB\n')
#np.savetxt(f, np.column_stack((zgrid, kperp, omd_curv, omd_gradB)))
#f.close()

for i in range(NEV):

   if EVrun:
       itime = i+1
   else:
       itime = -1

   print("\n################ Eigenvalue "+str(i+1)+" ################\n" )
   Qn2 = nrg[itime-1,6]/nrg[itime-1,0]

   if avg:
       # compute kz, kperp, omd using phi
       phi, apar = eigenfunctions_from_field_file(pars,suffix,False,False,setTime = itime)
       abs_phi = abs(phi)
       imax = np.argmax(abs_phi)
       phi00 = phi[imax]
       #print('phi00 = ', phi00)
       phi = phi / phi00
      
       f=open('phi_apar_'+str(i)+suffix,'w')
       f.write('#1.z 2.Re(phi) 3.Im(phi) 4.Re(apar) 5.Im(apar)\n')
       np.savetxt(f,np.column_stack((zgrid,np.real(phi),np.imag(phi),np.real(apar),np.imag(apar))))
       f.close()

       #print("len(kperp)",len(kperp))
       #print("len(phi)",len(phi))
       #plt.plot(zgrid,np.abs(phi))
       #plt.plot(zgrid,kperp/np.max(kperp))
       #plt.show()


       phi_plot = phi.copy()
       width = eigenfunction_width(zgrid,jacobian, phi, 'phi')
       #phi = np.real(phi)
       #print('phi = ', phi[0])
       #phi_smoothed01 = field_smoother(phi)
       #print('phi = ', phi[0])
       #phi_smoothed02 = field_smoother(phi_smoothed01)

   
       ave_kz = kz_from_dfielddz_bessel(kperp, zgrid,jacobian, phi, False, 'phi')
       ave_kperp, ave_omd = eigenfunction_average_bessel(zgrid, jacobian, kperp, omd_curv, phi,'phi')
       ave_kx = eigenfunction_average_kx(zgrid, jacobian, kx_ext, phi,'phi')
   
       print ('')
       print ('')
       print ('Using phi from GENE run')
       print ('kz = '+str( ave_kz))
       print ('kperp ='+str( ave_kperp))
       print ('kx ='+str( ave_kx))
       print ('omd ='+str( ave_omd))
       print ('gamma/kperp**2 ='+str( float(gamma[i]) / ave_kperp**2))
       print ('kz**2 ='+str( ave_kz**2))
       print ('omd * abs(omega) ='+str( ave_omd * np.sqrt(om[2]**2 + om[1]**2)))
       print ('width ='+str( width))
       print ('Input to SKiM')
       #print ('kperp ='+str( ave_kperp / float(pars['kymin'])))
       print ('omd ='+str( ave_omd / float(pars['kymin'])))
       print ('')
       print ('')
   
       f = open('mode_info_'+str(i+1) + suffix,'w')
       f.write('#1.ky 2.<kz> 3.kperp 4.omd 5.gamma 6.omega 7.gam/kperp^2 8.Q/n^2 9.num_zeros 10.kx_center 11.<kx>  12.width\n')
       ky = float(pars['kymin'])
       np.savetxt(f, np.column_stack((ky, 
                  ave_kz, 
                  ave_kperp, 
                  ave_omd,
                  float(gamma[i]),
                  float(omega[i]),
                  float(gamma[i]) / ave_kperp**2,
                  float(Qn2),np.nan,kx_center, ave_kx,width)),fmt='%4.4f', delimiter = "  ")
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

