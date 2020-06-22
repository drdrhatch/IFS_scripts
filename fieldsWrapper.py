#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from fieldlib import *
from finite_differences import *
#from finite_differences_x import *

def field_smoother(field):
    field_smooth = np.zeros(len(field), dtype = 'complex128')
    field_tmp = np.zeros(len(field), dtype = 'complex128')
    field_tmp[0] = field[0]
    field_tmp[len(field) - 1] = field[len(field) - 1]
    for i in range(1, len(field) - 1):
        field_tmp[i] = 0.5 * field[i] + 0.25 * (field[i - 1] + \
                          field[i + 1])
    field_smooth[0] = field_tmp[0]
    field_smooth[len(field) - 1] = field_tmp[len(field) - 1]
    for i in range(1, len(field) - 1):
        field_smooth[i] = 0.5 * field_tmp[i] + 0.25 * (field_tmp[i - 1] + \
                          field_tmp[i + 1])
    if 1 == 0:
        plt.plot(np.real(field), label = 're before')
        plt.plot(np.real(field_smooth), label = 're after')
        plt.legend()
        plt.show()
        plt.plot(np.imag(field), label = 'im before')
        plt.plot(np.imag(field_smooth), label = 'im after')
        plt.legend()
        plt.show()
    return field_smooth

def eigenfunctions_from_field_file(pars,suffix,center_only,plot,setTime=-1,smooth_field = False, scale_field = True):

    field = fieldfile('field'+suffix,pars)
    nz = int(field.nz)
    nx = int(field.nx)

    if (setTime == -1):
        field.set_time(field.tfld[setTime])
#        print 'Reading eigenfunctions are at t = ', field.tfld[setTime]
    else:
        isetTime = np.argmin(abs(np.array(field.tfld)-setTime))
        field.set_time(field.tfld[isetTime])
#        print 'Reading eigenfunctions are at t = ', field.tfld[isetTime]

    if center_only:
        ikx_grid = [0]
        phi = np.zeros(nz,dtype='complex128')
        apar = np.zeros(nz,dtype='complex128')
    else:
        print("nx",nx,"nx/2",nx/2,"floor(nx/2)",np.floor(nx/2))
        ikxmin = int( -np.ceil(nx/2)+1)
        ikx_grid = np.arange(nx)+ikxmin
        #ikx_grid = np.arange(-int(np.floor(nx/2)+1),int(np.floor(nx/2)+1))
        print("ikx_grid",ikx_grid)
        phi = np.zeros(nx*nz,dtype='complex128')
        apar = np.zeros(nx*nz,dtype='complex128')

    if 'n0_global' in pars:
        phase_fac = -np.e**(-2.0*np.pi*(0.0+1.0J)*int(pars['n0_global']) * float(pars['q0']))
    else:
        phase_fac = -1.0
    
    if float(pars['shat']) > 0.:
        for i in ikx_grid:
            this_phi = field.phi()[:,0,i]*phase_fac**i
            phi[(i-ikx_grid[0])*nz:(i-ikx_grid[0]+1)*nz]=this_phi
            if int(pars['n_fields']) > 1 and float(pars['beta']) !=0:
                this_apar = field.apar()[:,0,i]*phase_fac**i
                apar[(i-ikx_grid[0])*nz:(i-ikx_grid[0]+1)*nz]=\
                                                                 this_apar
    else:
        for i in ikx_grid:
            print("ikx_grid",ikx_grid)
            this_phi = field.phi()[:,0,-int(i)]*phase_fac**i
            print("len(this_phi)",len(this_phi))
            print("(i-ikx_grid[0])*nz,(i-ikx_grid[0]+1)*nz",(i-ikx_grid[0])*nz,(i-ikx_grid[0]+1)*nz)
            print("len(phi)",len(phi))
            phi[(i-ikx_grid[0])*nz:(i-ikx_grid[0]+1)*nz]=this_phi
            if pars['n_fields'] > 1 and pars['beta'] !=0:
                this_apar = field.apar()[:,0,-i]*phase_fac**i
                apar[(i-ikx_grid[0])*nz:(i-ikx_grid[0]+1)*nz]=\
                                                                 this_apar

    # Normalize phi and apar by highest value so that the peak abs val = 1
    if scale_field:
        phi = phi/np.max(abs(field.phi()[:,0,:]))
        if int(pars['n_fields']) > 1 and float(pars['beta']) !=0:
            apar = apar/np.max(abs(field.apar()[:,0,:]))
    if plot:
        if (setTime == -1):
            figTitle='t = '+ str(field.tfld[setTime])
        else:
            figTitle='t = '+ str(field.tfld[isetTime])
        if center_only:
            figTitle = figTitle+' center only'
        else:
            figTitle = figTitle+' entire simulation domain'
        plt.plot(np.real(phi),label='Re(phi)')
        plt.plot(np.imag(phi),label='Im(phi)')
        plt.plot(abs(phi),label='abs(phi)')
        plt.title(figTitle)
        plt.legend()
        plt.show()
    if plot and pars['n_fields'] > 1 and pars['beta'] !=0:
        plt.plot(np.real(apar),label='Re(apar)')
        plt.plot(np.imag(apar),label='Im(apar)')
        plt.plot(abs(apar),label='abs(apar)')
        plt.title(figTitle)
        plt.legend()
        plt.show()
    return phi, apar

def eigenfunction_average(z_grid,jacobian,kperp,omega_d,field,name):
    ave_kperp2 = 0.
    ave_omegad = 0.
    denom = 0.
    for i in np.arange(len(field)-1):
        ave_kperp2 = ave_kperp2 + (kperp[i]**2*abs(field[i])**2 +\
            kperp[i+1]**2*abs(field[i+1])**2)/2.*\
            (z_grid[i+1]-z_grid[i])/jacobian[i]
        ave_omegad = ave_omegad + (omega_d[i]*abs(field[i])**2 +\
            omega_d[i+1]*abs(field[i+1])**2)/2.*\
            (z_grid[i+1]-z_grid[i])/jacobian[i]
        denom = denom + (abs(field[i])**2 +abs(field[i+1])**2)/2.*\
            (z_grid[i+1]-z_grid[i])/jacobian[i]
    ave_kperp2 = ave_kperp2/denom
    ave_kperp = np.sqrt(ave_kperp2)
    #print name + ' weighted k_perp^2 =', ave_kperp2
    print( name + ' weighted k_perp =', ave_kperp)

    ave_omegad = ave_omegad/denom
    print( name + ' weighted omega_d =', ave_omegad)

    return ave_kperp, ave_omegad

def eigenfunction_squared(z_grid,jacobian,field):
    ave_sq_int = 0.
    ave_int_sq = 0.
    for i in np.arange(len(field)-1):
        ave_sq_int = ave_sq_int + ((field[i])**2 +\
            (field[i+1])**2)/2.*\
            (z_grid[i+1]-z_grid[i])/jacobian[i]
        ave_int_sq = ave_int_sq + (field[i] +\
            field[i+1])/2.*\
            (z_grid[i+1]-z_grid[i])/jacobian[i]

    ave_sq_int = abs(ave_sq_int)
    print( 'int phi**2 =', ave_sq_int)
    ave_int_sq= abs(ave_int_sq)**2
    print( ' (int phi)**2 =', ave_int_sq)

    return ave_sq_int, ave_int_sq

def kz_from_dfielddz(zgrid, jacobian, field, plot, name, zstart = 0., zend = 0.):
    dfielddz = np.empty(len(field),dtype='complex128')
    for i in range(len(field)-1):
        dfielddz[i] = (field[i+1]-field[i])/\
            (zgrid[i+1]-zgrid[i])*jacobian[i]
    if plot:
        plt.plot(zgrid, np.abs(dfielddz), label = 'abs d'+name+'/dz')
        plt.plot(zgrid, np.real(dfielddz), label = 'real d'+name+'/dz')
        plt.plot(zgrid, np.imag(dfielddz), label = 'imag d'+name+'/dz')
        plt.legend()
        plt.xlabel('z')
        plt.show()
    sum_ddz = 0.
    denom = 0.
    
    ##if not (zstart == 0. and zend == 0.):
    if zstart == zend:
        zstart = float(raw_input("Enter start z: "))
        zend = float(raw_input("Enter end z: "))
    startInd = np.argmin(abs(zgrid - zstart))
    endInd = np.argmin(abs(zgrid - zend))
    for i in range(startInd, endInd + 1):
        sum_ddz = sum_ddz + 0.5*(abs(dfielddz[i])**2+\
                  abs(dfielddz[i+1])**2)*\
                  (zgrid[i+1]-zgrid[i])/jacobian[i]
        denom = denom + 0.5*(abs(field[i])**2 + abs(field[i+1])**2)*\
                  (zgrid[i+1]-zgrid[i])/jacobian[i]
    ave_kz = np.sqrt(sum_ddz/denom)
    print( name + ' averaged kz = ', ave_kz)
    print( 'Input to SKiM kz = ', ave_kz)
    return ave_kz, zstart, zend

def fourierTrans(pars,zgrid,jacobian,field,plot,name):
    zi=complex(0,1)
    field_kz = np.empty(0,dtype='complex128')
    nkz = 100
    lkz = pars['nz0']/2
    kz_grid = np.linspace(-lkz,lkz,nkz,endpoint=False)
    for k in np.arange(len(kz_grid)):
        this_field_kz = 0.
        for i in np.arange(len(zgrid)-1):
            this_field_kz = this_field_kz + 0.5*(field[i]*\
                np.exp(zi*kz_grid[k]*zgrid[i])\
                +field[i+1]*np.exp(zi*kz_grid[k]*zgrid[i+1]))*\
                (zgrid[i+1]-zgrid[i])/jacobian[i]
        field_kz = np.append(field_kz,this_field_kz)
    if plot:
        plt.plot(kz_grid,np.abs(field_kz),label='abs('+name+'_kz)')
        plt.plot(kz_grid,np.real(field_kz),label='real('+name+'_kz)')
        plt.plot(kz_grid,np.imag(field_kz),label='imag('+name+'_kz)')
        plt.xlabel('kz')
        plt.title('ky = '+str(pars['kymin']))
        plt.legend()
        plt.show()

    kzstart = float(raw_input("Enter start kz: "))
    kzend = float(raw_input("Enter end kz: "))
    startInd = np.argmin(abs(kz_grid - kzstart))
    endInd = np.argmin(abs(kz_grid - kzend))
    sum_kz2 = 0.
    denom = 0.
    for i in range(len(kz_grid)-1):
        sum_kz2 = sum_kz2 + 0.5*(kz_grid[i]**2*abs(field_kz[i])**2+\
                  kz_grid[i+1]**2*abs(field_kz[i+1])**2)*\
                  (kz_grid[i+1]-kz_grid[i])
        denom = denom + 0.5*(abs(field_kz[i])**2 + \
                abs(field_kz[i+1])**2)*(kz_grid[i+1]-kz_grid[i])
    ave_kz = np.sqrt(sum_kz2/denom)
    print( name + ' averaged kz = ', ave_kz)
    #print 'input to SKiM averaged kz = ', ave_kz/np.pi/pars['q0']/pars['major_R']

    return field_kz, kz_grid

def LILO_eigenfunctions_from_field_file(pars,suffix,plot,setTime=-1):

    field = fieldfile('field'+suffix,pars)
    nz = int(field.nz)
    nx = int(field.nx)

    if (setTime == -1):
        field.set_time(field.tfld[setTime])
        print( 'Reading eigenfunctions are at t = ', field.tfld[setTime])
    else:
        isetTime = np.argmin(abs(np.array(field.tfld)-setTime))
        field.set_time(field.tfld[isetTime])
        print( 'Reading eigenfunctions are at t = ', field.tfld[isetTime])

    if 1 == 0:
        phi = np.zeros(nz,dtype='complex128')
        apar = np.zeros(nz,dtype='complex128')
    else:
        phi = np.zeros(field.nx,dtype='complex128')
        apar = np.zeros(field.nx,dtype='complex128')

    if 1 == 1:
        phi = field.phi()[nz/2,0,:]
        apar = field.apar()[nz/2,0,:]

    # Normalize phi and apar by highest value so that the peak abs val = 1
    phi = phi/np.max(abs(field.phi()[nz/2,0,:]))
    apar = apar/np.max(abs(field.apar()[nz/2,0,:]))
    if plot:
        if (setTime == -1):
            figTitle='t = '+ str(field.tfld[setTime])
        else:
            figTitle='t = '+ str(field.tfld[isetTime])
        plt.plot(np.real(phi),label='Re(phi)')
        plt.plot(np.imag(phi),label='Im(phi)')
        plt.plot(abs(phi),label='abs(phi)')
        plt.title(figTitle)
        plt.legend()
        plt.show()
    if plot and pars['n_fields'] > 1 and pars['beta'] !=0:
        plt.plot(np.real(apar),label='Re(apar)')
        plt.plot(np.imag(apar),label='Im(apar)')
        plt.plot(abs(apar),label='abs(apar)')
        plt.title(figTitle)
        plt.legend()
        plt.show()
    return phi, apar
