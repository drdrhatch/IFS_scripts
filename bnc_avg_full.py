#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import sys
from parIOWrapper import *
from geomWrapper import *
from fieldsWrapper import *
import optparse as op

parser=op.OptionParser(description='')
parser.add_option('--show_plots','-p', action='store',dest = 'show_plots', help = 'Eigenfunction average kperp, omega', default=False)
options,args=parser.parse_args()
suffix = args[0]
show_plots = options.show_plots
#peak_period = int(sys.argv[2])
#width = float(sys.argv[3])
#nl = int(sys.argv[2])
nl = 20


if not suffix=='.dat':
   suffix = '_'+suffix

me = 0.27240000E-03

def int_dz(f_z, zgrid, jacobian):
    sum_dz = 0.
    for i in range(nz - 1):
        sum_dz += (f_z[i] + f_z[i + 1]) /2. *\
                 (zgrid[i + 1] - zgrid[i]) / jacobian[i]
    return sum_dz

def cosine_theta(lambda_i, B_z):
    descrim = 1. - lambda_i * B_z
    if descrim >= 0.:
        return np.sqrt(descrim)
    else:
        return 0.

def N_avg(N_z, Q_z, zgrid, jacobian):
    dz_N_Q = int_dz(N_z * Q_z, zgrid, jacobian)
    dz_N = int_dz(N_z, zgrid, jacobian)
    return dz_N_Q / dz_N, dz_N

def lambda_int(f_lambda):
    return sum(f_lambda)


pars = init_read_parameters_file(suffix)
nz = int(pars['nz0'])
nx = int(pars['nx0'])

geom_type, geom_pars, geom_coeff = init_read_geometry_file(suffix, pars)
Bfield = geom_coeff['gBfield']

zgrid, jacobian = \
    reconstruct_zgrid(geom_coeff, pars, True, False, 0.)
zgrid_full, jacobian_full = \
    reconstruct_zgrid(geom_coeff, pars, False, False, 0.)

kperp, omd_curv, omd_gradB = \
    calc_kperp_omd(geom_type, geom_coeff, pars, True, False)
kperp_full, omd_curv_full, omd_gradB_full = \
    calc_kperp_omd(geom_type, geom_coeff, pars, False, False)

phi, apar = \
    eigenfunctions_from_field_file(pars, suffix, True, False, -1, False)
phi_full, apar_full = \
    eigenfunctions_from_field_file(pars, suffix, False, False, -1, False)

abs_phi = abs(phi)
imax = np.argmax(abs_phi)
phi00 = phi[imax]
phi = phi / phi00
phi = np.real(phi)

abs_phi_full = abs(phi_full)
imax = np.argmax(abs_phi_full)
phi_full00 = phi_full[imax]
phi_full = phi_full / phi_full00
phi_full = np.real(phi_full)

#if peak_period == 0:
#    phi = np.exp(-zgrid**2/width**2) #* np.cos(kz * zgrid)
#else:
#    phi = np.zeros(nz)
# construct lambda grid to be uniform or more points around 1/Bmax, 1/Bmin
uniform_lambda_grid = False
if uniform_lambda_grid:
    lambda_grid = np.linspace(1./max(Bfield), 1./min(Bfield), nl + 1, endpoint = True)
else:
    lambda_grid = [1./max(Bfield) + \
                  (1./min(Bfield) - 1./max(Bfield)) * 0.5 * (1. - \
                   np.cos(np.pi * i / nl)) \
                  for i in range(nl + 1)]

N_trapped = np.zeros((nl, nz), dtype = 'float128')
omd_grid = np.zeros((nl, nz), dtype = 'float128')
for i in range(nl):
    for j in range(nz):
        N_trapped[i, j] = cosine_theta(lambda_grid[i], Bfield[j]) - \
                          cosine_theta(lambda_grid[i + 1], Bfield[j])
        omd_grid[i, j] = omd_gradB[j] * (1. - (lambda_grid[i] + \
                          lambda_grid[i + 1])/2. * Bfield[j] / 2.)

N_avg_omd = []
int_dz_N = []
for i in range(nl):
    tmp1, tmp2 = N_avg(N_trapped[i,:], omd_grid[i,:],
          zgrid, jacobian)
    N_avg_omd.append(tmp1)
    int_dz_N.append(tmp2)

if show_plots:
    plt.plot(zgrid,abs(phi),label='abs phi',color='black')
    plt.plot(zgrid,np.real(phi),label='real phi',color='blue')
    plt.plot(zgrid,np.imag(phi),label='imag phi',color='red')
    plt.xlabel('z(pi)')
    plt.grid()
    plt.legend()
    plt.show()

N_avg_phi = []
for i in range(nl):
    tmp1, tmp2 = N_avg(N_trapped[i,:], phi,
          zgrid, jacobian)
    N_avg_phi.append(tmp1)

lambda_int_N = []
for i in range(nz):
    lambda_int_N.append(lambda_int(N_trapped[:,i]))
N_avg_phi = np.array(N_avg_phi)
int_dz_N = np.array(int_dz_N)
int_dz_phi2 = int_dz(abs(phi)**2, zgrid, jacobian)

frac_trap = lambda_int(abs(N_avg_phi)**2 * int_dz_N) \
            / int_dz_phi2

bnc_avg_omd = lambda_int(N_avg_omd *abs(N_avg_phi)**2 * int_dz_N) \
              /lambda_int(abs(N_avg_phi)**2 * int_dz_N)

f = open('bnc2_avg_full' + suffix, 'w')
sys.stdout = f
print('[central period only] frac_trap = {}'.format(np.round(frac_trap, 4)))
print('[central period only] bounce averaged omega_d = {}'.format(np.round(bnc_avg_omd, 4)))
print('[central period only] bounce averaged for SKiM omega_d = {}'.format(np.round(bnc_avg_omd / pars['kymin'], 4)))


if show_plots:
    plt.plot(zgrid_full, abs(phi_full), label='abs phi full',color='black')
    plt.plot(zgrid_full, np.real(phi_full),label='real phi full',color='blue')
    plt.plot(zgrid_full, np.imag(phi_full),label='imag phi full',color='red')
    plt.xlabel('z(pi)')
    plt.grid()
    plt.legend()
    plt.show()

frac_trap_numer = []
frac_trap_denom = []
bnc_avg_omd_numer = []
bnc_avg_omd_denom = []

ikx_grid_full = np.arange(-nx/2+1,nx/2+1)
for i in ikx_grid_full:
    this_phi = phi_full[(i-ikx_grid_full[0])*nz:(i-ikx_grid_full[0]+1)*nz]

    this_zgrid = zgrid_full[(i-ikx_grid_full[0])*nz:(i-ikx_grid_full[0]+1)*nz]
    this_jacobian = jacobian_full[(i-ikx_grid_full[0])*nz:(i-ikx_grid_full[0]+1)*nz]
#    if i == peak_period:
#        this_phi = np.exp(-(this_zgrid - i*2)**2/width**2) #* np.cos(kz * zgrid)
#    else:
#        this_phi = np.zeros(nz)
#    plt.plot(this_zgrid, this_phi)
    
    N_avg_this_phi = []
    for j in range(nl):
        tmp1, tmp2 = N_avg(N_trapped[j, :], this_phi, this_zgrid, this_jacobian)
        N_avg_this_phi.append(tmp1)
    N_avg_this_phi = np.array(N_avg_this_phi)
    int_dz_this_phi2 = int_dz(abs(this_phi)**2, this_zgrid, this_jacobian)
    
    frac_trap_numer.append(lambda_int(abs(N_avg_this_phi)**2 * int_dz_N))
    frac_trap_denom.append(int_dz_this_phi2)
    bnc_avg_omd_numer.append(lambda_int(N_avg_omd * abs(N_avg_this_phi)**2 * \
                             int_dz_N))
    bnc_avg_omd_denom.append(lambda_int(abs(N_avg_this_phi)**2 * int_dz_N))

frac_trap_full = sum(frac_trap_numer) / sum(frac_trap_denom)
bnc_avg_omd_full = sum(bnc_avg_omd_numer) / sum(bnc_avg_omd_denom)

print('[full simulation periods] frac_trap = {}'.format(np.round(frac_trap_full, 4)))
print('[full simulation periods] bounce averaged omega_d = {}'.format(np.round(bnc_avg_omd_full, 4)))
print('[full simulation periods] bounce averaged for SKiM omega_d = {}'.format(np.round(bnc_avg_omd_full / pars['kymin'], 4)))

f.close()
