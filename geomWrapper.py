#!/usr/bin/env python
# -*- coding: utf-8 -*-

from read_write_geometry import * 
import matplotlib.pyplot as plt
import sys
import optparse as op
from parIOWrapper import init_read_parameters_file
from finite_differences import *

def init_read_geometry_file(suffix, pars):

    geom_type = pars['magn_geometry'][1:-1]
    geom_file = geom_type + suffix
    geom_pars, geom_coeff = read_geometry_local(geom_file)

    return geom_type, geom_pars, geom_coeff
   
def init_read_geometry_file_glob(suffix, pars):

    geom_type = pars['magn_geometry'][1:-1]
    geom_file = geom_type + suffix
    geom_pars, geom_coeff = read_geometry_global(geom_file)

    return geom_type, geom_pars, geom_coeff

def plot_grid_points(geom_coeff, plot = True):

    Z = geom_coeff['gl_z']
    R = geom_coeff['gl_R']
    
    if plot:
        plt.scatter(R,Z)
        plt.xlabel('R')
        plt.ylabel('Z')
        plt.axis('equal')
        plt.title('simulation grid points')
        plt.show()

    return R, Z

def read_geom_coeff_raw(geom_type, geom_coeff, plot = False):

    ggxx = geom_coeff['ggxx']
    ggxy = geom_coeff['ggxy']
    ggxz = geom_coeff['ggxz']
    ggyy = geom_coeff['ggyy']
    ggyz = geom_coeff['ggyz']
    ggzz = geom_coeff['ggzz']

    gdBdx = geom_coeff['gdBdx']
    gdBdy = geom_coeff['gdBdy']
    gdBdz = geom_coeff['gdBdz']

    gBfield = geom_coeff['gBfield']
    gjacobian = geom_coeff['gjacobian']

    gl_R = geom_coeff['gl_R']
    gl_phi = geom_coeff['gl_phi']
    gl_z = geom_coeff['gl_z']
    gl_dxdR = geom_coeff['gl_dxdR']
    gl_dxdZ = geom_coeff['gl_dxdZ']

    return ggxx,ggxy,ggxz,ggyy,ggyz,ggzz,gdBdx,gdBdy,gdBdz,gBfield,gjacobian,\
           gl_R,gl_phi,gl_z,gl_dxdR,gl_dxdZ

def read_curv_coeff(geom_type, geom_coeff, plot = False):

    ggxx = geom_coeff['ggxx']
    ggxy = geom_coeff['ggxy']
    ggxz = geom_coeff['ggxz']
    ggyy = geom_coeff['ggyy']
    ggyz = geom_coeff['ggyz']
    ggzz = geom_coeff['ggzz']

    gamma1 = ggxx * ggyy - ggxy ** 2
    gamma2 = ggxx * ggyz - ggxy * ggxz
    gamma3 = ggxy * ggyz - ggyy * ggxz

    gdBdx = geom_coeff['gdBdx']
    gdBdy = geom_coeff['gdBdy']
    gdBdz = geom_coeff['gdBdz']

    gBfield = geom_coeff['gBfield']

    if 1 == 0:
        plt.plot(gBfield,label='Bfield')
        plt.legend()
        plt.show()

    Kx = - gdBdy - gamma2 / gamma1 * gdBdz
    Ky = gdBdx - gamma3 / gamma1 * gdBdz

    if (geom_type == 's_alpha'):
        Kx = Kx / gBfield
        Ky = Ky / gBfield

    if plot:
        plt.plot(Kx,label='K_x')
        plt.plot(Ky,label='K_y')
        plt.title('curvature coefficients')
        plt.legend()
        plt.show()

    return Kx, Ky

def reconstruct_zgrid(geom_coeff, pars, center_only, plot = True, edge_opt = -1):

    nx = int(pars['nx0'])
    nz = int(pars['nz0'])
    if 'gBfield' in geom_coeff:
        gBfield = geom_coeff['gBfield']
    else:
        gBfield = geom_coeff['Bfield']
    if 'gjacobian' in geom_coeff:
        gjacobian = geom_coeff['gjacobian']
    else:
        gjacobian = geom_coeff['jacobian']

    if center_only:
        ikx_grid = [0]
    else:
        ikx_grid = np.arange(- nx // 2 + 1, nx // 2 + 1)

    zgrid_even_center = np.linspace(-1., 1., nz, endpoint = False)
    #print 'Debug: ', 'dz_even =', \
    #      zgrid_even_center[1] - zgrid_even_center[0], \
    #      'should be equal to 2/nz =', 2./nz

    if not center_only:
        if nx % 2 == 1:
            zgrid_even = np.linspace(- nx, nx, nx * nz, endpoint = False)
        else :
            zgrid_even = np.linspace(- (nx - 1), (nx + 1), nx * nz, \
                                                        endpoint = False)
    
    if 'edge_opt' in pars:
      if edge_opt == -1:
        edge_opt = pars['edge_opt']
      else:
        edge_opt = float(edge_opt)
    else:
      edge_opt = 0.

    if edge_opt != 0:
        #sys.exit('edge_opt ~= 0 code for zgrid is not ready.')
        zgrid_edge = np.zeros(nx * nz, dtype = 'float128')

        
        N = np.arcsinh(edge_opt*zgrid_even_center[0]*np.pi)/\
            zgrid_even_center[0]/np.pi
        zgrid_edge_center = 1./edge_opt*np.sinh(N*zgrid_even_center*np.pi)/np.pi

        dz = np.zeros(nz, dtype = 'float128')
        for i in np.arange(int(nz / 2 + 1), nz):
            dz[i] = zgrid_edge_center[i]-zgrid_edge_center[i - 1]
        for i in np.arange(int(nz / 2 - 1), - 1, - 1):
            dz[i] = zgrid_edge_center[i + 1]-zgrid_edge_center[i]

        for i in ikx_grid:
            this_zgrid_edge = i * 2. + zgrid_edge_center
            zgrid_edge[(i - ikx_grid[0]) * nz: (i - ikx_grid[0] + 1) * nz] \
                      = this_zgrid_edge
        if center_only:
            zgrid = zgrid_edge_center
        else:
            zgrid = zgrid_edge
    else:
        if 'edge_opt' in pars and pars['edge_opt'] != 0:
            print ('Warning:edge_opt ~= 0 code for zgrid is not ready.')
        if center_only:
            zgrid = zgrid_even_center
        else:
            zgrid = zgrid_even

        dz = np.ones(nz, dtype = 'float128')*2./nz

    jacobian_center = 1./np.pi/gjacobian/gBfield
    
    if not center_only:
        jacobian = np.zeros(nx * nz, dtype = 'float128')
        for i in ikx_grid:
            jacobian[(i-ikx_grid[0])*nz:(i-ikx_grid[0]+1)*nz]=jacobian_center
    if center_only:
        jacobian = jacobian_center

    if plot:
        plt.plot(zgrid, label = 'zgrid')
        plt.legend(loc=2)
        plt.show()
        plt.plot(jacobian, label = 'jacobian')
        plt.legend(loc=2)
        plt.show()

    return zgrid, jacobian

def calc_kperp_omd(geom_type, geom_coeff,pars,center_only,plot, ky =-1):

    nx = int(pars['nx0'])
    if center_only:
        ikx_grid = [0]
    else:
        ikx_grid = np.arange(- nx // 2 + 1, nx // 2 + 1)
    nz = int(pars['nz0'])
    lx = float(pars['lx'])
    if ky == -1:
        ky = float(pars['kymin'])
    print ('ky = ', ky)
    dkx = 2. * np.pi * float(pars['shat']) * float(ky)

    dpdx_tot = float(pars['beta']) * \
               (float(pars['omn1']) + float(pars['omt1']))
    if  int(pars['n_spec']) > 1:
        dpdx_tot = dpdx_tot + float(pars['beta']) * \
                   (float(pars['omn2']) + float(pars['omt2']))
        if int(pars['n_spec']) > 2:
           dpdx_tot = dpdx_tot + float(pars['beta']) * \
                      (float(pars['omn3']) + float(pars['omt3']))
    #print 'Debug: ', 'dpdx_pm =', pars['dpdx_pm'], \
    #      'should be equal to beta*(omni+omti+omne+omte) =', dpdx_tot

    if 'kx_center' in pars:
        kx_center = float(pars['kx_center'])
    else:
        kx_center = 0.

    Kx, Ky = read_curv_coeff(geom_type, geom_coeff, False)
    ggxx = geom_coeff['ggxx'].astype(float)
    ggxy = geom_coeff['ggxy'].astype(float)
    ggyy = geom_coeff['ggyy'].astype(float)
    gBfield = geom_coeff['gBfield'].astype(float)

    if center_only:
        kperp = np.zeros(nz,dtype='float128')
        omd_curv = np.zeros(nz,dtype='float128')
        omd_gradB = np.zeros(nz,dtype='float128')
    else:
        kperp = np.zeros(nx*nz,dtype='float128')
        omd_curv = np.zeros(nx*nz,dtype='float128')
        omd_gradB = np.zeros(nx*nz,dtype='float128')

    for i in ikx_grid:
        kx = i*dkx+kx_center
        this_kperp = np.sqrt(ggxx*kx**2+2.*ggxy*kx*ky+ggyy*ky**2)
        this_omegad_gradB = -(Kx*kx+Ky*ky)/gBfield
        this_omegad_curv = this_omegad_gradB + \
                           ky * float(pars['dpdx_pm'])/gBfield**2/2.
        #this_omegad_curv = 2.*this_omegad

        kperp[(i-ikx_grid[0])*nz:(i-ikx_grid[0]+1)*nz]=this_kperp
        omd_curv[(i-ikx_grid[0])*nz:(i-ikx_grid[0]+1)*nz]=this_omegad_curv
        omd_gradB[(i-ikx_grid[0])*nz:(i-ikx_grid[0]+1)*nz]=this_omegad_gradB

    if geom_type == 's_alpha' or geom_type == 'slab':
        if 'amhd' in pars:
            amhd = pars['amhd']
        else:
            amhd = 0.
        z_grid = np.linspace(-1.,1., nz, endpoint = False)
        Kx0 = -np.sin(z_grid*np.pi)/pars['major_R']
        Ky0 = -(np.cos(z_grid*np.pi)+np.sin(z_grid*np.pi)*\
              (pars['shat']*z_grid*np.pi-amhd*\
              np.sin(z_grid*np.pi)))/pars['major_R']
        omega_d0 = -(Kx0*kx_center+Ky0*ky)
        omega_d00 = omega_d0+amhd/pars['q0']**2/pars['major_R']/2.*ky/gBfield**2
        gxx0 = 1.
        gxy0 = pars['shat']*z_grid*np.pi-amhd*np.sin(z_grid*np.pi)
        gyy0 = 1+(pars['shat']*z_grid*np.pi-amhd*np.sin(z_grid*np.pi))**2
        kperp0 = np.sqrt(gxx0*kx_center**2+2.*gxy0*kx_center*ky+gyy0*ky**2)
    if plot:
        plt.plot(kperp,label='kperp')
        plt.title('entire simulation domain')
        if geom_type == 's_alpha' and plot and center_only:
            plt.plot(kperp0,label='check')
            plt.title('center only')
        plt.legend()
        plt.show()
        plt.plot(omd_curv,label='omd_curv')
        plt.plot(omd_gradB,label='omd_gradB')
        plt.title('entire simulation domain')
        if geom_type == 's_alpha' and plot and center_only:
            plt.plot(omega_d0,label='check')
            plt.plot(omega_d00,label='check')
            plt.title('center only')
        plt.legend()
        plt.show()

    return kperp, omd_curv, omd_gradB


def calc_kx_extended(pars,plot, ky =-1):

    nx = int(pars['nx0'])
    ikx_grid = np.arange(- nx // 2 + 1, nx // 2 + 1)
    nz = int(pars['nz0'])
    lx = float(pars['lx'])
    if ky == -1:
        ky = float(pars['kymin'])
    print ('ky = ', ky)
    dkx = 2. * np.pi * float(pars['shat']) * float(ky)

    dpdx_tot = float(pars['beta']) * \
               (float(pars['omn1']) + float(pars['omt1']))
    if  int(pars['n_spec']) > 1:
        dpdx_tot = dpdx_tot + float(pars['beta']) * \
                   (float(pars['omn2']) + float(pars['omt2']))
        if int(pars['n_spec']) > 2:
           dpdx_tot = dpdx_tot + float(pars['beta']) * \
                      (float(pars['omn3']) + float(pars['omt3']))

    if 'kx_center' in pars:
        kx_center = float(pars['kx_center'])
    else:
        kx_center = 0.

    kx_ext = np.zeros(nx*nz,dtype='float128')

    for i in ikx_grid:
        kx = i*dkx+kx_center

        kx_ext[(i-ikx_grid[0])*nz:(i-ikx_grid[0]+1)*nz]=kx

    if plot:
        plt.plot(kx_ext,label='kperp')
        plt.title('entire simulation domain')
        plt.legend()
        plt.show()

    return kx_ext



def bounce_averaged_omd(suffix,pars,geom_coeff,omega_d1,omega_d2,z_grid,ky=-1):
    nl = 1024
    me = 0.27240000E-03
    lx = pars['lx']
    nz = pars['nz0']
    nx = pars['nx0']
    if ky == -1:
        ky = pars['kymin']
    print ('ky = ', ky)
    nx = 1

    gBfield = geom_coeff['gBfield']
    gdBdz = geom_coeff['gdBdz']

    dBdz_c = np.zeros(nz-1,dtype='float128')
    Bfield_c = gBfield
    for k in range(nz-1):
        dBdz_c[k] = (gBfield[k+1]-gBfield[k])/\
                    (z_grid[k+1]-z_grid[k])
    if 1 == 0:
        plt.plot(dBdz_c,'.',label='dBdz_check')
        plt.legend()
        plt.show()
    if 1 == 0:
        plt.plot(gBfield,'.',label='Bfield_gene')
        plt.legend()
        plt.show()

    #l_grid_long = np.linspace(1./max(Bfield_c),1./min(Bfield_c),nl+1,endpoint = False)
    #l_grid = l_grid_long[1:]
    l_grid = np.linspace(1./max(Bfield_c),1./min(Bfield_c),nl,endpoint = False)

    J = np.zeros(nl,dtype='float128')
    t_b = np.zeros(nl,dtype='float128')
    bnc_avg_omd1 = np.zeros(nl,dtype='float128')
    bnc_avg_omd2 = np.zeros(nl,dtype='float128')
    bnc_avg_omd3 = np.zeros(nl,dtype='float128')
    v_parallel = np.zeros((nl,nz),dtype='float128') 
    v_perp = np.zeros((nl,nz),dtype='float128') 
    for j in np.arange(nl):
        for k in np.arange(nz):
             if 1. - l_grid[j] * Bfield_c[k] < 0.:
                 v_parallel[j,k] = 0.
                 v_perp[j,k] = np.sqrt(1./me)
             else:
                 v_parallel[j,k] = np.sqrt(1. - \
                                 l_grid[j] * Bfield_c[k]) * np.sqrt(1./me)
                 v_perp[j,k] = np.sqrt(l_grid[j]*Bfield_c[k]) * np.sqrt(1./me)
        if 1 == 0:
            plt.plot(v_parallel[j,:],'*',label='v_parallel')
            plt.plot(v_perp[j,:],'.',label='v_perp')
            plt.plot(np.sqrt(v_parallel[j,:]**2+v_perp[j,:]**2),'.',label='v_tot')
            plt.legend()
            plt.show()

    if 1 == 0:
        for j in np.arange(nl):
            for k in np.arange(nz-1):
                J[j] = J[j] + (v_parallel[j,k] + v_parallel[j,k+1]) / 2. *\
                          (z_grid[k+1] - z_grid[k])
        plt.plot(l_grid,J,'.',label='invariant J')
        plt.legend()
        plt.show()

    this_t_b = np.zeros(nz-1,dtype='float128')
    this_ba_omd1 = np.zeros(nz-1,dtype='float128')
    this_ba_omd2 = np.zeros(nz-1,dtype='float128')
    this_ba_omd3 = np.zeros(nz-1,dtype='float128')
    for j in np.arange(nl):
        for k in np.arange(nz-1):
            if dBdz_c[k]==0:
                this_t_b[k] = (z_grid[k+1]-z_grid[k])*\
                              2./(v_parallel[j,k]+v_parallel[j,k+1])
            else:
                this_t_b[k] = 1./l_grid[j]/dBdz_c[k]*\
                              (v_parallel[j,k]-v_parallel[j,k+1])
            this_ba_omd1[k] = (omega_d1[k]+omega_d1[k+1])/2.*this_t_b[k]
            this_ba_omd2[k] = (omega_d2[k]+omega_d2[k+1])/2.*this_t_b[k]
            this_ba_omd3[k] = ((omega_d1[k]*v_parallel[j,k]**2 + \
                              omega_d1[k+1]*v_parallel[j,k+1]**2)\
                              +(omega_d2[k]*v_perp[j,k]**2/2. + \
                              omega_d2[k+1]*v_perp[j,k+1]**2/2.))/2.*this_t_b[k]
            #this_ba_omd3[k] = ((omega_d1[k]*v_parallel[j,k]**2 + \
            #                  omega_d1[k+1]*v_parallel[j,k+1]**2)\
            #                  /2.)*this_t_b[k]
        if 1 == 0:
            #plt.plot(this_t_b,'*',label='time in each interval')
            plt.plot(v_parallel[j,:]**2*me,label='v_parallel**2')
            plt.plot(v_perp[j,:]**2*me,label='v_perp**2')
            plt.plot(this_ba_omd1,'.',label='omega_d1(lambda,z)')
            plt.plot(this_ba_omd2,'.',label='omega_d2(lambda,z)')
            plt.plot(this_ba_omd3*me,'*',label='omega_d3(lambda,z)')
            plt.legend()
            plt.show()
 
        t_b[j] = sum(this_t_b)
        bnc_avg_omd1[j] = sum(this_ba_omd1)/t_b[j]
        bnc_avg_omd2[j] = sum(this_ba_omd2)/t_b[j]
        bnc_avg_omd3[j] = sum(this_ba_omd3)/t_b[j]

    if 1 == 0:
        plt.plot(l_grid,t_b,'.',label='half bounce period')
        plt.legend(loc=2)
        plt.show()
    plt.plot(l_grid,bnc_avg_omd1, '.', \
    label='bnc avg omega_drift v_parallel')
    plt.plot(l_grid,bnc_avg_omd2, '.', \
    label='bnc avg omega_drift v_perp')
    plt.plot(l_grid,bnc_avg_omd3*me, '.', \
    label='bnc avg omega_drift total')
    #plt.legend()
    plt.show()
    plt.plot(l_grid,bnc_avg_omd3/bnc_avg_omd3[-1], '.', \
    label='bnc avg omega_drift total')
    #plt.legend()
    plt.show()

    print ('bounce averaged omd = ', np.mean(bnc_avg_omd3)*me)
    print ('normalized <omd> = ', np.mean(bnc_avg_omd3)/bnc_avg_omd1[-1]*me)

def calc_shatloc(geom_coeff, z_grid, plot = False):

    temp = geom_coeff['ggxy']/geom_coeff['ggxx']
    shatLoc = fd_d1_o4(temp, z_grid)
    if plot:
        plt.plot(z_grid, shatLoc,label='shat loc')
        plt.legend()
        plt.show()   
    return shatLoc

def smoothWdiff(geomInput,s=0.3,nmax=3):

#    plt.plot(geomInput,label='input')
    n = 0
    newArray = np.empty(len(geomInput),dtype='complex128')
    while n < nmax:
        for i in range(len(geomInput)):
            if i == 0 or i == len(geomInput)-1:
                newArray[i] = geomInput[i]
            else:
                newArray[i] = geomInput[i] + s*(geomInput[i+1]+geomInput[i-1]\
                              - 2.*geomInput[i])
#        plt.plot(newArray,label='after '+str(n+1)+' iterations')
        geomInput = newArray
        n = n + 1
#    plt.legend()
#    plt.show()
    return geomInput
               
def smoothWhypdiff(geomInput,nmax=3):

#    plt.plot(geomInput,label='input')
    n = 0
    newArray = np.empty(len(geomInput),dtype='complex128')
    while n < nmax:
        for i in range(len(geomInput)):
            if i == 0 or i == len(geomInput)-1:
                newArray[i] = geomInput[i]
            elif i == 1 or i == len(geomInput)-2:
                newArray[i] = geomInput[i] + 0.25*(geomInput[i+1]+geomInput[i-1]\
                              - 2.*geomInput[i])
            else:
                newArray[i] = -geomInput[i+2] / 12. +\
                              geomInput[i+1] / 3. + geomInput[i] / 2. + \
                              geomInput[i-1] / 3. - geomInput[i-2] / 12.
#        plt.plot(newArray,label='after '+str(n+1)+' iterations')
        geomInput = newArray
        n = n + 1
#    plt.legend()
#    plt.show()
    return geomInput

def ktheta_factor(geom_pars, geom_coeff, plot = False):

    ggyy = geom_coeff['ggyy']
    ggyz = geom_coeff['ggyz']
    ggzz = geom_coeff['ggzz']

    q = geom_pars['q0']
    Cy = geom_pars['Cy']
    ktheta = q * Cy / np.sqrt(q**2 * Cy**2 * ggyy + \
             2. * q * Cy * ggyz + ggzz)

    if plot:
        plt.plot(ktheta,label='ktheta factor')
        plt.legend()
        plt.show()

    print ('ktheta factor at z = 0', ktheta[geom_pars['gridpoints']/2])

def k2_factor(geom_type, geom_coeff, plot = False):

    ggxx = geom_coeff['ggxx']
    ggxy = geom_coeff['ggxy']
    ggyy = geom_coeff['ggyy']

    gamma1 = ggxx * ggyy - ggxy ** 2

    k2 = np.sqrt(gamma1/ggxx)

    if plot:
        plt.plot(k2,label='k2 factor')
        plt.legend()
        plt.show()

    print('k2 factor min =', np.min(k2))

    return np.min(k2)

def k2_factor_global(geom_type, geom_coeff, xInd, plot = False):

    ggxx = geom_coeff['gxx'][:, xInd]
    ggxy = geom_coeff['gxy'][:, xInd]
    ggyy = geom_coeff['gyy'][:, xInd]

    gamma1 = ggxx * ggyy - ggxy ** 2

    k2 = np.sqrt(gamma1/ggxx)

    if plot:
        plt.plot(k2,label='k2 factor')
        plt.legend()
        plt.show()

    print ('k2 factor min =', np.min(k2))

    return np.min(k2)
def ky(pars, geom_coeff, plot):

    ggxx = geom_coeff['ggxx']
    ggxy = geom_coeff['ggxy']
    ggyy = geom_coeff['ggyy']

    gamma1 = ggxx * ggyy - ggxy ** 2

    kymin = pars['kymin']
    ky = np.sqrt(gamma1/ggxx)*kymin

    R, Z = plot_grid_points(geom_coeff, False)
    if plot:
        plt.plot(R,label='R')
        plt.plot(Z,label='Z')
        plt.legend()
        plt.show()
        plt.plot(ky,label='ky')
        plt.legend()
        plt.show()

    return ky

def ky_global(pars, geom_coeff, xInd):

    ggxx = geom_coeff['gxx'][:, xInd]
    ggxy = geom_coeff['gxy'][:, xInd]
    ggyy = geom_coeff['gyy'][:, xInd]

    gamma1 = ggxx * ggyy - ggxy ** 2

    kymin = pars['kymin']
    ky = np.sqrt(gamma1/ggxx)*kymin

    return ky

def kthetaConversion(ktheta_cm, pars):

    T_ev = float(pars['Tref']) * 1.E03
    B_Gauss = float(pars['Bref']) * 1.E04
    rho_cm = 102. * np.sqrt(float(pars['mref'])) * np.sqrt(T_ev) / B_Gauss
    ktheta_GENE = ktheta_cm * rho_cm

    return ktheta_GENE
