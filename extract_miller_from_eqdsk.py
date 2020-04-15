#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 09:32:25 2013

@author: dtold
"""

import numpy as np
import matplotlib.pyplot as plt
from sys import exit, stdout
import optparse as op

nr = 120
ntheta = 160
parser = op.OptionParser(description='Extract Miller shaping parameters from EQDSK files.')
parser.add_option('--rovera', '-r', action='store_const', const=1)
parser.add_option('--conv', '-c', action='store_const', const=1)
parser.add_option('--noplots', '-n', action='store_const', const=1)
options, args = parser.parse_args()
use_r_a = options.rovera
show_plots = (not options.noplots)
write_rhotor_rhopol_conversion = options.conv
if not write_rhotor_rhopol_conversion:
    if len(args) != 2:
        exit("""
Please give two arguments: <EQDSK filename> <Position in rho_tor>
optional: -r <Position in r/a>
          -c: Write rhotor/rhopol conversion to file
          -n: Suppress plots\n""")

filename = args[0]
file = open(filename, 'r')


def find(val, arr):
    return np.argmin(np.abs(arr - val))


eqdsk = file.readlines()
print(('Header: {0:s}'.format(eqdsk[0])))
# set resolutions
nw = np.int(eqdsk[0].split()[-2])
nh = np.int(eqdsk[0].split()[-1])
pw = np.int((nw/8/2)*2)  # psi-width, number of flux surfaces around position of interest
print(('Resolution: {0:4d} x {1:4d}'.format(nw, nh)))

entrylength = 16
try:
    rdim, zdim, rctr, rmin, zmid = [np.float(eqdsk[1][j*entrylength:(j + 1)*entrylength]) for j in
                                    range(len(eqdsk[1])//entrylength)]
except:
    entrylength = 15
    try:
        rdim, zdim, rctr, rmin, zmid = [np.float(eqdsk[1][j*entrylength:(j + 1)*entrylength]) for j
                                        in range(len(eqdsk[1])//entrylength)]
    except:
        exit('Error reading EQDSK file, please check format!')
rmag, zmag, psiax, psisep, Bctr = [np.float(eqdsk[2][j*entrylength:(j + 1)*entrylength]) for j in
                                   range(len(eqdsk[2])//entrylength)]
_, psiax2, _, rmag2, _ = [np.float(eqdsk[3][j*entrylength:(j + 1)*entrylength]) for j in
                          range(len(eqdsk[3])//entrylength)]
zmag2, _, psisep2, _, _ = [np.float(eqdsk[4][j*entrylength:(j + 1)*entrylength]) for j in
                           range(len(eqdsk[4])//entrylength)]
if rmag != rmag2:
    np.sys.exit('Inconsistent rmag: %7.4g, %7.4g'%(rmag, rmag2))
if psiax2 != psiax:
    np.sys.exit('Inconsistent psiax: %7.4g, %7.4g'%(psiax, psiax2))
if zmag != zmag2:
    np.sys.exit('Inconsistent zmag: %7.4g, %7.4g'%(zmag, zmag2))
if psisep2 != psisep:
    np.sys.exit('Inconsistent psisep: %7.4g, %7.4g'%(psisep, psisep2))
F = np.empty(nw, dtype=np.float)
p = np.empty(nw, dtype=np.float)
ffprime = np.empty(nw, dtype=np.float)
pprime = np.empty(nw, dtype=np.float)
qpsi = np.empty(nw, dtype=np.float)
psirz_1d = np.empty(nw*nh, dtype=np.float)
start_line = 5
lines = list(range(nw//5))
if nw%5 != 0:
    lines = list(range(nw//5 + 1))
for i in lines:
    n_entries = len(eqdsk[i + start_line])//entrylength
    F[i*5:i*5 + n_entries] = [np.float(eqdsk[i + start_line][j*entrylength:(j + 1)*entrylength]) for
                              j in range(n_entries)]
start_line = i + start_line + 1

for i in lines:
    n_entries = len(eqdsk[i + start_line])//entrylength
    p[i*5:i*5 + n_entries] = [np.float(eqdsk[i + start_line][j*entrylength:(j + 1)*entrylength]) for
                              j in range(n_entries)]
start_line = i + start_line + 1

for i in lines:
    n_entries = len(eqdsk[i + start_line])//entrylength
    ffprime[i*5:i*5 + n_entries] = [
        np.float(eqdsk[i + start_line][j*entrylength:(j + 1)*entrylength]) for j in
        range(n_entries)]
start_line = i + start_line + 1

for i in lines:
    n_entries = len(eqdsk[i + start_line])//entrylength
    pprime[i*5:i*5 + n_entries] = [
        np.float(eqdsk[i + start_line][j*entrylength:(j + 1)*entrylength]) for j in
        range(n_entries)]
start_line = i + start_line + 1

lines_twod = list(range(nw*nh//5))
if nw*nh%5 != 0:
    lines_twod = list(range(nw*nh//5 + 1))
for i in lines_twod:
    n_entries = len(eqdsk[i + start_line])//entrylength
    psirz_1d[i*5:i*5 + n_entries] = [
        np.float(eqdsk[i + start_line][j*entrylength:(j + 1)*entrylength]) for j in
        range(n_entries)]
start_line = i + start_line + 1
psirz = psirz_1d.reshape(nh, nw)

for i in lines:
    n_entries = len(eqdsk[i + start_line])//entrylength
    qpsi[i*5:i*5 + n_entries] = [np.float(eqdsk[i + start_line][j*entrylength:(j + 1)*entrylength])
                                 for j in range(n_entries)]
start_line = i + start_line + 1

# invert sign of psi if necessary to guarantee increasing values for interpolation
if psisep < psiax:
    psirz = -psirz
    ffprime = -ffprime
    pprime = -pprime
    psiax *= -1
    psisep *= -1

# ignore limiter data etc. for the moment
dw = rdim/(nw - 1)
dh = zdim/(nh - 1)
rgrid = np.array([rmin + i*dw for i in range(nw)])
zgrid = np.array([zmid - zdim/2. + i*dh for i in range(nh)])
# contourf(rgrid,zgrid,psirz,70);gca().set_aspect('equal')
# show()

# create 5th order 2D spline representation of Psi(R,Z)
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline

interpol_order = 3
psi_spl = RectBivariateSpline(zgrid, rgrid, psirz, kx=interpol_order, ky=interpol_order)

# linear grid of psi, on which all 1D fields are defined
linpsi = np.linspace(psiax, psisep, nw)
# create rho_tor grid
x_fine = np.linspace(psiax, psisep, nw*10)
phi_fine = np.empty((nw*10), dtype=np.float)

phi_fine[0] = 0.
# spline of q for rhotor grid
q_spl_psi = UnivariateSpline(linpsi, qpsi, k=interpol_order, s=1e-5)
q_fine = q_spl_psi(x_fine)

for i in range(1, nw*10):
    x = x_fine[:i + 1]
    y = q_spl_psi(x)
    phi_fine[i] = np.trapz(y, x)
rho_tor_fine = np.sqrt(phi_fine/phi_fine[-1])
rho_tor_spl = UnivariateSpline(x_fine, rho_tor_fine, k=interpol_order, s=1e-5)
rho_tor = np.empty(nw, dtype=np.float)
for i in range(nw):
    rho_tor[i] = rho_tor_spl(linpsi[i])

if write_rhotor_rhopol_conversion:
    rt_rp_filename = 'rt_rp_%s'%filename
    rt_rp_file = open(rt_rp_filename, 'w')
    rt_rp_file.write('# rho_tor          rho_pol      q\n')
    for i in range(len(x_fine)):
        rho_pol = np.sqrt((x_fine[i] - psiax)/(psisep - psiax))
        rt_rp_file.write('%16.8e %16.8e %16.8e\n'%(rho_tor_fine[i], rho_pol, q_fine[i]))
    rt_rp_file.close()
    exit('\nWrote rhotor/rhopol conversion to %s.'%rt_rp_filename)

t1 = np.arctan2(zmid - zdim/2. - zmag, rmin - rmag)
t2 = np.arctan2(zmid - zdim/2. - zmag, rmin + rdim - rmag)
t3 = np.arctan2(zmid + zdim/2. - zmag, rmin + rdim - rmag)
t4 = np.arctan2(zmid + zdim/2. - zmag, rmin - rmag)

theta_arr = np.linspace(-np.pi, np.pi, ntheta)
# for i in range(nw):
#    curr_psi=linpsi[i]

print('Finding flux surface shapes...')
R = np.empty((nw, ntheta), dtype=np.float)
Z = np.empty((nw, ntheta), dtype=np.float)
dr = rdim*np.cos(theta_arr)
dz = rdim*np.sin(theta_arr)

for j in range(len(theta_arr)):
    stdout.write('\r Finished %4.1f%%.'%(j*100./(ntheta - 1)))
    stdout.flush()
    theta = theta_arr[j]
    r_pol = np.linspace(rmag, rmag + dr[j], nr)  # array([rmag+i*dr for i in range(nr)])
    z_pol = np.linspace(zmag, zmag + dz[j], nr)  # array([zmag+i*dz for i in range(nr)])
    psi_rad = psi_spl.ev(z_pol, r_pol)
    psi_rad_sav = psi_rad
    psi_rad[0] = psiax
    # must restrict interpolation range because of non-monotonic psi around coils
    cutoff = 0
    for i in range(1, len(psi_rad)):
        if psi_rad[i] < psi_rad[i - 1]:
            cutoff = i
            break
    psi_rad = psi_rad[:i]
    end_ind = np.argmin(np.abs(psi_rad - psisep))
    end_ind += (1 if (psi_rad[end_ind] < psisep) else 0)
    indsep = end_ind + 1

    R_int = interp1d(psi_rad[:indsep], r_pol[:indsep], kind=interpol_order)
    R[:, j] = R_int(linpsi)
    Z_int = interp1d(psi_rad[:indsep], z_pol[:indsep], kind=interpol_order)
    Z[:, j] = Z_int(linpsi)

print('\nFinding flux surface centers...')
# find average elevation for all flux surfaces
Z_avg = np.empty(nw, dtype=np.float)
ds = np.empty(ntheta, dtype=np.float)
for i in range(1, nw):
    ds[1:ntheta - 1] = 0.5*np.sqrt(
            (R[i, 2:ntheta] - R[i, 0:ntheta - 2]) ** 2 + (Z[i, 2:ntheta] - Z[i, 0:ntheta - 2]) ** 2)
    ds[0] = 0.5*np.sqrt((R[i, 1] - R[i, -1]) ** 2 + (Z[i, 1] - Z[i, -1]) ** 2)
    ds[-1] = 0.5*np.sqrt((R[i, 0] - R[i, -2]) ** 2 + (Z[i, 0] - Z[i, -2]) ** 2)
    Z_avg[i] = np.average(Z[i, :], weights=ds)
# Treat the magnetic axis separately as no average is required and ds==0
Z_avg[0] = Z[0, 0]

# find R0,Z0 for all flux surfaces
R0 = np.empty(nw, dtype=np.float)
R0[0] = rmag
r_avg = np.empty(nw, dtype=np.float)
r_avg[0] = 0.
r_maxmin = np.empty(nw, dtype=np.float)
r_maxmin[0] = 0.
for i in range(1, nw):
    stdout.write('\r Finished %4.1f%%.'%(i*100./(nw - 1)))
    stdout.flush()
    R_array = R[i, ntheta//4:3*ntheta//4]
    Z_array = Z[i, ntheta//4:3*ntheta//4]
    # low field side
    Z_int = interp1d(Z_array, list(range(ntheta//2)), kind=interpol_order)
    ind_Zavg = Z_int(Z_avg[i])
    R_int = interp1d(list(range(ntheta//2)), R_array, kind=interpol_order)
    R_out = R_int(ind_Zavg)
    R_max = np.amax(R_array)
    # high field side
    R_array = np.roll(R[i, :-1], ntheta//2)[ntheta//4:3*ntheta//4]
    Z_array = np.roll(Z[i, :-1], ntheta//2)[ntheta//4:3*ntheta//4]

    # have to use negative Z_array here to have increasing order
    Z_int = interp1d(-Z_array, list(range(ntheta//2)), kind=interpol_order)
    # again negative
    ind_Zavg = Z_int(-Z_avg[i])

    R_int = interp1d(list(range(ntheta//2)), R_array, kind=interpol_order)
    R_in = R_int(ind_Zavg)
    R_min = np.amin(R_array)
    R0[i] = 0.5*(R_out + R_in)
    r_avg[i] = 0.5*(R_out - R_in)
    r_maxmin[i] = 0.5*(R_max - R_min)

radpos = np.float(args[1])
if use_r_a:
    r_a = radpos
    # find psi index of interest (for the specified r/a position)
    poi_ind = find(radpos, r_avg/r_avg[-1])
    print(('\nExamine {0:3d} flux surfaces around position rho_tor={1:7.4g}...'.format(pw, radpos)))
else:
    # find psi index of interest (for the specified rho_tor position)
    poi_ind = find(radpos, rho_tor)
    ravg_rho_spl = UnivariateSpline(rho_tor, r_avg/r_avg[-1], k=interpol_order, s=1e-5)
    r_a = np.float(ravg_rho_spl(radpos))
    print(('\nExamine {0:3d} flux surfaces around position r/a={1:7.4g}...'.format(pw, r_a)))

# modified theta grid for each flux surface
# arrays equidistant on modified theta grid are marked by 'tm' index!!!
linpsi_spl = UnivariateSpline(r_avg, linpsi, k=interpol_order, s=1e-5)
ravg_spl = UnivariateSpline(linpsi, r_avg, k=interpol_order, s=1e-5)

rmaxmin_spl = UnivariateSpline(linpsi, r_maxmin, k=interpol_order, s=1e-5)
q_spl = UnivariateSpline(r_avg, qpsi, k=interpol_order, s=1e-5)
R0_spl = UnivariateSpline(r_avg, R0, k=interpol_order, s=1e-5)
Z0_spl = UnivariateSpline(r_avg, Z_avg, k=interpol_order, s=1e-5)
F_spl = UnivariateSpline(r_avg, F, k=interpol_order, s=1e-5)
r = r_a*r_avg[-1]
psi = np.float(linpsi_spl(r))
psi_N = (psi - psiax)/(psisep - psiax)
R0_pos = np.float(R0_spl(r))
Z0_pos = np.float(Z0_spl(r))
F_pos = np.float(F_spl(r))
Bref_miller = F_pos/R0_pos

print(('Coordinates: r={0:8.5g}, psi={1:8.5g}, psi_N={2:8.5g}, r/R0={3:8.5g}, rho_tor={4:8.5g}, '
      'r_maxmin={5:8.5g}'.format(r, psi, psi_N, r/R0_pos, np.float(rho_tor_spl(psi)),
                                 np.float(rmaxmin_spl(psi)))))
psi_stencil = list(range(poi_ind - pw//2, poi_ind + pw//2))
if psi_stencil[0] < 1:
    psi_stencil = [psi_stencil[i] + 1 - psi_stencil[0] for i in range(len(psi_stencil))]
if psi_stencil[-1] > nw - 1:
    psi_stencil = [psi_stencil[i] - (psi_stencil[-1] - nw + 1) for i in range(len(psi_stencil))]
R_tm = np.empty((pw, ntheta), dtype=np.float)
Z_tm = np.empty((pw, ntheta), dtype=np.float)
R_extended = np.empty(2*ntheta - 1, dtype=np.float)
Z_extended = np.empty(2*ntheta - 1, dtype=np.float)
# theta_mod[0]=theta_arr
# R_tm[0]=R[0]
# Z_tm[0]=Z[0]
theta_tmp = np.linspace(-2.*np.pi, 2*np.pi, 2*ntheta - 1)

print('Interpolating to flux-surface dependent (proper) theta grid...')
for i in psi_stencil:
    stdout.write('\r Finished %4.1f%%.'%(psi_stencil.index(i)*100./(len(psi_stencil) - 1)))
    stdout.flush()
    imod = i - psi_stencil[0]
    # print('Finished {0:4.1f}%%.'.format(float(i)/(pw-1)*100))
    R_extended[0:(ntheta - 1)//2] = R[i, (ntheta + 1)//2:-1]
    R_extended[(ntheta - 1)//2:(3*ntheta - 3)//2] = R[i, :-1]
    R_extended[(3*ntheta - 3)//2:] = R[i, 0:(ntheta + 3)//2]
    Z_extended[0:(ntheta - 1)//2] = Z[i, (ntheta + 1)//2:-1]
    Z_extended[(ntheta - 1)//2:(3*ntheta - 3)//2] = Z[i, :-1]
    Z_extended[(3*ntheta - 3)//2:] = Z[i, 0:(ntheta + 3)//2]
    # for j in range(ntheta):
    theta_mod_ext = np.arctan2(Z_extended - Z_avg[i], R_extended - R0[i])
    # introduce 2pi shifts to theta_mod_ext
    for ind in range(ntheta):
        if theta_mod_ext[ind + 1] < 0. and theta_mod_ext[ind] > 0. and np.abs(
                theta_mod_ext[ind + 1] - theta_mod_ext[ind]) > np.pi:
            lshift_ind = ind
        if theta_mod_ext[-ind - 1] > 0. and theta_mod_ext[-ind] < 0. and np.abs(
                theta_mod_ext[-ind - 1] - theta_mod_ext[-ind]) > np.pi:
            rshift_ind = ind
    theta_mod_ext[-rshift_ind:] += 2.*np.pi
    theta_mod_ext[:lshift_ind + 1] -= 2.*np.pi
    # print theta_mod, theta_arr
    #    plot(theta_mod_ext)
    #    plot(theta_tmp)
    #    show()
    theta_int = interp1d(theta_mod_ext, theta_tmp, kind=interpol_order)
    theta_orig_tm = theta_int(theta_arr)
    R_int = interp1d(theta_mod_ext, R_extended, kind=interpol_order)
    Z_int = interp1d(theta_mod_ext, Z_extended, kind=interpol_order)
    R_tm[imod] = R_int(theta_arr)
    Z_tm[imod] = Z_int(theta_arr)
#    plot(R_tm[imod],Z_tm[imod])
# gca().set_aspect('equal')

# now we have the flux surfaces on a symmetric grid in theta (with reference to R0(r), Z0(r))
# symmetrize flux surfaces
# figure()
R_sym = np.empty((pw, ntheta), dtype=np.float)
Z_sym = np.empty((pw, ntheta), dtype=np.float)
for i in psi_stencil:
    imod = i - psi_stencil[0]
    Z_sym[imod, :] = 0.5*(Z_tm[imod, :] - Z_tm[imod, ::-1]) + Z_avg[i]
    R_sym[imod, :] = 0.5*(R_tm[imod, :] + R_tm[imod, ::-1])
#    plot(R_sym[imod],Z_sym[imod])
# gca().set_aspect('equal')
# show()

dq_dr_avg = np.empty(pw, dtype=np.float)
dq_dpsi = np.empty(pw, dtype=np.float)
drR = np.empty(pw, dtype=np.float)
drZ = np.empty(pw, dtype=np.float)
kappa = np.empty(pw, dtype=np.float)
delta = np.empty(pw, dtype=np.float)
s_kappa = np.empty(pw, dtype=np.float)
s_delta = np.empty(pw, dtype=np.float)
delta_upper = np.empty(pw, dtype=np.float)
delta_lower = np.empty(pw, dtype=np.float)
zeta_arr = np.empty((pw, 4), dtype=np.float)
zeta = np.empty(pw, dtype=np.float)
s_zeta = np.empty(pw, dtype=np.float)
for i in psi_stencil:
    imod = i - psi_stencil[0]
    # calculate delta
    stencil_width = ntheta//10
    for o in range(2):
        if o:
            ind = np.argmax(Z_sym[imod])
            section = list(range(ind + stencil_width//2, ind - stencil_width//2, -1))
        else:
            ind = np.argmin(Z_sym[imod])
            section = list(range(ind - stencil_width//2, ind + stencil_width//2))
        x = R_sym[imod, section]
        y = Z_sym[imod, section]
        y_int = interp1d(x, y, kind=interpol_order)
        x_fine = np.linspace(np.amin(x), np.amax(x), stencil_width*100)
        y_fine = y_int(x_fine)
        if o:
            x_at_extremum = x_fine[np.argmax(y_fine)]
            delta_upper[imod] = (R0[i] - x_at_extremum)/r_avg[i]
            Z_max = np.amax(y_fine)
        else:
            x_at_extremum = x_fine[np.argmin(y_fine)]
            delta_lower[imod] = (R0[i] - x_at_extremum)/r_avg[i]
            Z_min = np.amin(y_fine)
    # calculate kappa
    kappa[imod] = (Z_max - Z_min)/2./r_avg[i]

# linear extrapolation (in psi) for axis values
# delta_upper[0]=2*delta_upper[1]-delta_upper[2]
# delta_lower[0]=2*delta_lower[1]-delta_lower[2]
# kappa[0]=2*kappa[1]-kappa[2]
# zeta[0]=2*zeta[1]-zeta[2]
delta = 0.5*(delta_upper + delta_lower)

# calculate zeta
for i in psi_stencil:
    imod = i - psi_stencil[0]
    x = np.arcsin(delta[imod])
    # find the points that correspond to Miller-theta=+-pi/4,+-3/4*pi and extract zeta from those
    for o in range(4):
        if o == 0:
            val = np.pi/4.
            searchval = np.cos(val + x/np.sqrt(2))
            searcharr = (R_sym[imod] - R0[i])/r_avg[i]
        elif o == 1:
            val = 3.*np.pi/4
            searchval = np.cos(val + x/np.sqrt(2))
            searcharr = (R_sym[imod] - R0[i])/r_avg[i]
        elif o == 2:
            val = -np.pi/4.
            searchval = np.cos(val - x/np.sqrt(2))
            searcharr = (R_sym[imod] - R0[i])/r_avg[i]
        elif o == 3:
            val = -3.*np.pi/4
            searchval = np.cos(val - x/np.sqrt(2))
            searcharr = (R_sym[imod] - R0[i])/r_avg[i]
        if o in [0, 1]:
            searcharr2 = searcharr[ntheta//2:]
            ind = find(searchval, searcharr2) + ntheta//2
        else:
            searcharr2 = searcharr[0:ntheta//2]
            ind = find(searchval, searcharr2)
        #        print o,ind
        section = list(range(ind - stencil_width//2, ind + stencil_width//2))
        theta_sec = theta_arr[section]
        if o in [0, 1]:
            theta_int = interp1d(-searcharr[section], theta_sec, kind=interpol_order)
            theta_of_interest = theta_int(-searchval)
        else:
            theta_int = interp1d(searcharr[section], theta_sec, kind=interpol_order)
            theta_of_interest = theta_int(searchval)
        Z_sec = Z_sym[imod, section]
        Z_sec_int = interp1d(theta_sec, Z_sec, kind=interpol_order)
        #        print searchval,val, theta_sec
        Z_val = Z_sec_int(theta_of_interest)
        zeta_arr[imod, o] = np.arcsin((Z_val - Z_avg[i])/kappa[imod]/r_avg[i])
    zeta_arr[imod, 1] = np.pi - zeta_arr[imod, 1]
    zeta_arr[imod, 3] = -np.pi - zeta_arr[imod, 3]
    #    print zeta_arr[i]
    zeta[imod] = 0.25*(
            np.pi + zeta_arr[imod, 0] - zeta_arr[imod, 1] - zeta_arr[imod, 2] + zeta_arr[imod, 3])

Bref_efit = np.abs(F[0]/R0[0])
Lref_efit = np.sqrt(2*np.abs(phi_fine[-1])/Bref_efit)

kappa_spl = UnivariateSpline(r_avg[psi_stencil], kappa, k=interpol_order, s=1e-5)
delta_spl = UnivariateSpline(r_avg[psi_stencil], delta, k=interpol_order, s=1e-5)
zeta_spl = UnivariateSpline(r_avg[psi_stencil], zeta, k=interpol_order, s=1e-5)
amhd = np.empty(pw, dtype=np.float)
amhd_Miller = np.empty(pw, dtype=np.float)
Vprime = np.empty(pw, dtype=np.float)
dV_dr = np.empty(pw, dtype=np.float)
V = np.empty(pw, dtype=np.float)
V_manual = np.empty(pw, dtype=np.float)
r_FS = np.empty(pw, dtype=np.float)
for i in psi_stencil:
    imod = i - psi_stencil[0]
    Vprime[imod] = np.abs(np.sum(qpsi[i]*R_sym[imod] ** 2/F[i])*4*np.pi ** 2/ntheta)
    dV_dr[imod] = np.abs(np.sum(qpsi[i]*R_sym[imod] ** 2/F[i])*4*np.pi ** 2/ntheta)/ \
                  ravg_spl.derivatives(linpsi[i])[1]
    #    V[imod]=trapz(Vprime[:imod+1],linpsi[psi_stencil])
    r_FS[imod] = np.average(np.sqrt((R_sym[imod] - R0[i]) ** 2 + (Z_sym[imod] - Z_avg[i]) ** 2),
                            weights=qpsi[i]*R_sym[imod] ** 2/F[i])
    amhd[imod] = -qpsi[i] ** 2*R0[i]*pprime[i]*8*np.pi*1e-7/Bref_miller ** 2/ \
                 ravg_spl.derivatives(linpsi[i])[1]
    #    amhd_Miller[imod]=-2*Vprime[imod]/(2*pi)**2*(V[imod]/2/pi**2/R0[i])**0.5*4e-7*pi*pprime[i]
    dq_dr_avg[imod] = q_spl.derivatives(r_avg[i])[1]
    dq_dpsi[imod] = q_spl_psi.derivatives(linpsi[i])[1]
    drR[imod] = R0_spl.derivatives(r_avg[i])[1]
    drZ[imod] = Z0_spl.derivatives(r_avg[i])[1]
    s_kappa[imod] = kappa_spl.derivatives(r_avg[i])[1]*r_avg[i]/kappa[imod]
    s_delta[imod] = delta_spl.derivatives(r_avg[i])[1]*r_avg[i]/np.sqrt(1 - delta[imod] ** 2)
    s_zeta[imod] = zeta_spl.derivatives(r_avg[i])[1]*r_avg[i]
amhd_spl = UnivariateSpline(r_avg[psi_stencil], amhd, k=interpol_order, s=1e-5)
rFS_spl = UnivariateSpline(r_avg[psi_stencil], r_FS, k=interpol_order, s=1e-5)
drR_spl = UnivariateSpline(r_avg[psi_stencil], drR, k=interpol_order, s=1e-5)
drZ_spl = UnivariateSpline(r_avg[psi_stencil], drZ, k=interpol_order, s=1e-5)
Zavg_spl = UnivariateSpline(r_avg, Z_avg, k=interpol_order, s=1e-5)

# plot(r_avg[psi_stencil],V)
# figure()
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(r_avg[psi_stencil], kappa)
ax.set_title('Elongation')
ax.set_xlabel(r'$r_{avg}$', fontsize=14)
ax.set_ylabel(r'$\kappa$', fontsize=14)
ax.axvline(r, 0, 1, ls='--', color='k', lw=2)

fig2 = plt.figure()
ax2 = fig2.add_subplot(1, 1, 1)
ax2.plot(r_avg[psi_stencil], s_kappa)
ax2.set_title('Elongation (Derivative)')
ax2.set_xlabel(r'$r_{avg}$', fontsize=14)
ax2.set_ylabel(r'$s_\kappa$', fontsize=14)
ax2.axvline(r, 0, 1, ls='--', color='k', lw=2)

fig3 = plt.figure()
ax3 = fig3.add_subplot(1, 1, 1)
ax3.plot(r_avg[psi_stencil], delta)
ax3.set_title('Triangularity')
ax3.set_xlabel(r'$r_{avg}$', fontsize=14)
ax3.set_ylabel(r'$\delta$', fontsize=14)
ax3.axvline(r, 0, 1, ls='--', color='k', lw=2)

fig4 = plt.figure()
ax4 = fig4.add_subplot(1, 1, 1)
ax4.plot(r_avg[psi_stencil], s_delta)
ax4.set_title('Triangularity (Derivative)')
ax4.set_xlabel(r'$r_{avg}$', fontsize=14)
ax4.set_ylabel(r'$s_\delta$', fontsize=14)
ax4.axvline(r, 0, 1, ls='--', color='k', lw=2)
# figure()
# plot(r_avg[psi_stencil],delta_diff)
# xlabel(r'$r_{avg}$',fontsize=14)
# ylabel(r'$\delta_u-\delta_l$',fontsize=14)

fig5 = plt.figure()
ax5 = fig5.add_subplot(1, 1, 1)
ax5.plot(r_avg[psi_stencil], zeta)
ax5.set_title('Squareness')
ax5.set_xlabel(r'$r_{avg}$', fontsize=14)
ax5.set_ylabel(r'$\zeta$', fontsize=14)
ax5.axvline(r, 0, 1, ls='--', color='k', lw=2)

fig6 = plt.figure()
ax6 = fig6.add_subplot(1, 1, 1)
ax6.plot(r_avg[psi_stencil], s_zeta)
ax6.set_title('Squareness (Derivative)')
ax6.set_xlabel(r'$r_{avg}$', fontsize=14)
ax6.set_ylabel(r'$s_\zeta$', fontsize=14)
ax6.axvline(r, 0, 1, ls='--', color='k', lw=2)

# select a given flux surface
ind = poi_ind - psi_stencil[0]  # find(r_a,r_avg/r_avg[-1])
Lref = R0_pos
print(('\n\nShaping parameters for flux surface r={0:9.5g}, r/a={1:9.5g}:'.format(r, r_a)))
# print 'r_FS= %9.5g (flux-surface averaged radius)\n' %rFS_spl(r)
print('Copy the following block into a GENE parameters file:\n')
print(('trpeps  = {0:9.5g}'.format(r/R0_pos)))
print(('q0      = {0:9.5g}'.format(np.float(q_spl(r)))))
print(('shat    = {0:9.5g} !(defined as r/q*dq_dr)'.format(
        r/np.float(q_spl(r))*q_spl.derivatives(r)[1])))
# print 'shat=%9.5g (defined as (psi-psiax)/q*dq_dpsi)' %((psi-linpsi[0])/q_spl(
# r)*q_spl_psi.derivatives(psi)[1])
print(('amhd    = {0:9.5g}'.format(np.float(amhd_spl(r)))))
# print 'amhd_Miller=%9.5g' %amhd_Miller[ind]
# print test[ind]
print(('drR     = {0:9.5g}'.format(np.float(drR_spl(r)))))
print(('drZ     = {0:9.5g}'.format(np.float(drZ_spl(r)))))
print(('kappa   = {0:9.5g}'.format(np.float(kappa_spl(r)))))
print(('s_kappa = {0:9.5g}'.format(kappa_spl.derivatives(r)[1]*r/np.float(kappa_spl(r)))))
print(('delta   = {0:9.5g}'.format(np.float(delta_spl(r)))))
print(('s_delta = {0:9.5g}'.format(
        delta_spl.derivatives(r)[1]*r/np.sqrt(1 - np.float(delta_spl(r)) ** 2))))
print(('zeta    = {0:9.5g}'.format(np.float(zeta_spl(r)))))
print(('s_zeta  = {0:9.5g}'.format(zeta_spl.derivatives(r)[1]*r)))
print(('minor_r = {0:9.5g}'.format(1.0)))
print(('major_R = {0:9.5g}'.format(R0_pos/r*r_a)))
print('\nFor normalization to major radius, set instead:')
print(('minor_r = {0:9.5g}'.format(r/R0_pos/r_a)))
print(('major_R = {0:9.5g}'.format(1.0)))
print('The same conversion factor must be considered for input frequencies and gradients.')
print('\nAdditional information:')
print(('Lref        = {0:9.5g} !for Lref=a convention'.format(r_avg[-1])))
print(('Lref        = {0:9.5g} !for Lref=R0 convention'.format(R0_pos)))
print(('Bref        = {0:9.5g}'.format(Bref_miller)))
# minor radius at average flux surface elevation (GYRO definition)
print(('\na (avg elev)= {0:9.5g}'.format(r_avg[-1])))
print(('R0          = {0:9.5g}'.format(R0_pos)))
print(('Z0          = {0:9.5g}'.format(np.float(Zavg_spl(r)))))
print(('Lref_efit   = {0:9.5g}'.format(Lref_efit)))
print(('Bref_efit   = {0:9.5g}'.format(Bref_efit)))
print(('B_unit(GYRO)= {0:9.5g}'.format(q_spl(r)/r/ravg_spl.derivatives(psi)[1])))
# print 'Vprime=%9.5g; dV_dr=%9.5g' %(Vprime[ind],dV_dr[ind])
# print 'V=%9.5g' %V[ind]
# minor radius defined by 0.5(R_max-R_min), where R_max and R_min can have any elevation
# print 'a_maxmin    = %9.5g' %r_maxmin[-1]
print(('dpsi/dr     = {0:9.5g}'.format(1./ravg_spl.derivatives(psi)[1])))
# print 'dr_maxmin/dr     = %9.5g' %(rmaxmin_spl.derivatives(psi)[1]/ravg_spl.derivatives(psi)[1])
print(('drho_tor/dr = {0:9.5g}'.format(rho_tor_spl.derivatives(psi)[1]/ravg_spl.derivatives(psi)[1])))
print(('Gradient conversion omt(rho_tor) -> a/LT; factor = {0:9.5g}'.format(
        r_avg[-1]*(rho_tor_spl.derivatives(psi)[1]/ravg_spl.derivatives(psi)[1]))))
print(('Gradient conversion omt(rho_tor) -> R/LT; factor = {0:9.5g}'.format(
        R0_pos*(rho_tor_spl.derivatives(psi)[1]/ravg_spl.derivatives(psi)[1]))))

fig7 = plt.figure()
ax7 = fig7.add_subplot(1, 1, 1)
ax7.plot(R_tm[ind], Z_tm[ind], 'k-', lw=2, label='original')
ax7.plot(R_sym[ind], Z_sym[ind], 'r-', lw=2, label='symmetrized')
# for v in range(-10,10,4):
# zeta_test=-0.011086#zeta[ind]
ax7.plot(R0_pos + r*np.cos(theta_arr + np.arcsin(delta_spl(r))*np.sin(theta_arr)),
         Zavg_spl(r) + kappa_spl(r)*r*np.sin(theta_arr + zeta_spl(r)*np.sin(2*theta_arr)),
         label='Miller')
ax7.set_title('Flux surface shapes')
ax7.set_xlabel('$R$/m', fontsize=14)
ax7.set_ylabel('$Z$/m', fontsize=14)
ax7.set_aspect('equal')
ax7.legend(loc=10, prop={'size': 10})

provide_conversions = 0
if provide_conversions:
    f = open('rho_tor-r_a_conv', 'w')
    f.write('#r/a    rho_tor\n')
    for i in range(0, nw):
        f.write('%16.8e %16.8e\n'%(ravg_spl(linpsi[i])/r_avg[-1], rho_tor_spl(linpsi[i])))

if show_plots:
    plt.show()

# #Fourier series representation of flux surfaces
# RZ_c=R+1j*Z
# for i in range(0,nw,4):
#    RZ_c_fft=fft(RZ_c[i])
#    #plot(RZ_c_fft)
#    trunc_order=2
#    RZ_c_fft[trunc_order:nw/2]=0.;RZ_c_fft[nw/2+1:-trunc_order+1]=0.
#    RZ_trunc=ifft(RZ_c_fft)
#    R_trunc=real(RZ_trunc)
#    Z_trunc=imag(RZ_trunc)
#    plot(R_trunc,Z_trunc,lw=2,ls='--')
# contour(rgrid,zgrid,psirz,20)
# gca().set_aspect('equal')
# show()
