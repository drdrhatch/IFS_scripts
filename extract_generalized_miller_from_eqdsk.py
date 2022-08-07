#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 02 17:51:25 2020

@authors: tbg, partially based on extract_miller_from_eqdsk.py by dtold
"""

import numpy as np
import matplotlib.pyplot as plt
from sys import exit, stdout
import optparse as op
from scipy import interpolate
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from create_toroidal_flux_grid import create_toroidal_flux


nr = 256 #120
ntheta = 256 #160
parser = op.OptionParser(description='Extract generalized Miller shaping parameters from EQDSK files.')
parser.add_option('--rovera', '-r', action='store_const', const=1)
parser.add_option('--noplots', '-n', action='store_const', const=1)
parser.add_option('--conv', '-c', action='store_const', const=1)
parser.add_option('--method', '-m', action='store', type='int', dest="method",default=0)
options, args = parser.parse_args()
use_r_a = options.rovera
show_plots = (not options.noplots)
write_rhotor_rhopol_conversion = options.conv
tfmethod=options.method
if not write_rhotor_rhopol_conversion:
    if len(args) != 2:
        exit("""
Please give two arguments: <EQDSK filename> <Position in rho_tor>
optional: -r <Position in r/a>
          -c: Write rhotor/rhopol conversion to file
          -n: Suppress plots\n
          -m 0|1: Use 1d/2d method for creating rho_tor grid\n""")

filename = args[0]
radpos = float(args[1])

file = open(filename, 'r')


def get_flux_surface_distance(R0,Z0,R,Z):
    nz0 = len(R)
    #poloidal theta (COCOS=2/12); theta is running counter-clockwise
    #in a poloidal plain to the right of the magnetic axis
    theta = np.arctan((Z-Z0)/(R-R0))
    for iz in range(int(nz0/2),nz0-1):
        thetadiff = (theta[iz+1]-theta[iz])
        if thetadiff>0.5*np.pi:
            theta[iz+1]-=np.pi
        elif thetadiff<-0.5*np.pi:
            theta[iz+1]+=np.pi
        thetadiff = (theta[nz0-1-iz-1]-theta[nz0-1-iz])
        if thetadiff>0.5*np.pi:
            theta[nz0-1-iz-1]-=np.pi
        elif thetadiff<-0.5*np.pi:
            theta[nz0-1-iz-1]+=np.pi

    #sorted theta (starting at th=0)
    theta=np.where(abs(theta)<1E-12,0.0,theta)
#        print ('theta = ', theta)
    theta_sort = (theta+2.0*np.pi)%(2.0*np.pi)
    sort_ind = np.argsort(theta_sort)
    #        print('sort_ind = ', sort_ind)
    theta_sort = theta_sort[sort_ind]
    #        print('theta_sort = ', theta_sort)
    aN = np.sqrt((R[sort_ind]-R0)**2+(Z[sort_ind]-Z0)**2) #/R0

    return {"theta" : theta,          #poloidal angle for original data
            "sort_ind" : sort_ind,    #sort index for original data
            "aN" : aN,                #sorted flux surface radius funct.
            "theta_sort" : theta_sort #sorted theta
    }

def interpolate_periodic(theta_in,theta_out,var_in):
    #axis and data needs to be sorted
    #now extend var and theta_sort to definitely cover >= [0.0,2*pi]
    theta_ext = theta_in
    var_ext = var_in
    shift = 0
    while theta_ext[0] > 0.0:
        theta_ext = np.append(theta_ext[-1]-2.0*np.pi,theta_ext)
        var_ext = np.append(var_ext[-1],var_ext)
        shift += 1
    while theta_ext[-1] < 2.0*np.pi:
        theta_ext = np.append(theta_ext,theta_ext[shift]+2.0*np.pi)
        var_ext = np.append(var_ext,var_ext[shift])
    #interpolate to [0,2*pi]:
    fc = interpolate.interp1d(theta_ext,var_ext,kind='cubic') #,fill_value='extrapolate')
    return fc(theta_out)

def find(val, arr):
    return np.argmin(np.abs(arr - val))

####################### Main program ###############################

###first read g-eqdsk file
eqdsk = file.readlines()
print('Header: {0:s}'.format(eqdsk[0]))
# set resolutions
nw = np.int(eqdsk[0].split()[-2])
nh = np.int(eqdsk[0].split()[-1])
pw = max([np.int((nw/8/2)*2),6])  # psi-width, number of flux surfaces around position of interest
print('Resolution: {0:4d} x {1:4d}'.format(nw, nh))

entrylength = 16
try:
    rdim, zdim, rctr, rmin, zmid = [float(eqdsk[1][j*entrylength:(j + 1)*entrylength]) for j in
                                    range(len(eqdsk[1])//entrylength)]
except:
    entrylength = 15
    try:
        rdim, zdim, rctr, rmin, zmid = [float(eqdsk[1][j*entrylength:(j + 1)*entrylength]) for j
                                        in range(len(eqdsk[1])//entrylength)]
    except:
        exit('Error reading EQDSK file, please check format!')
rmag, zmag, psiax, psisep, Bctr = [float(eqdsk[2][j*entrylength:(j + 1)*entrylength]) for j in
                                   range(len(eqdsk[2])//entrylength)]
_, psiax2, _, rmag2, _ = [float(eqdsk[3][j*entrylength:(j + 1)*entrylength]) for j in
                          range(len(eqdsk[3])//entrylength)]
zmag2, _, psisep2, _, _ = [float(eqdsk[4][j*entrylength:(j + 1)*entrylength]) for j in
                           range(len(eqdsk[4])//entrylength)]
if rmag != rmag2:
    np.sys.exit('Inconsistent rmag: %7.4g, %7.4g'%(rmag, rmag2))
if psiax2 != psiax:
    np.sys.exit('Inconsistent psiax: %7.4g, %7.4g'%(psiax, psiax2))
if zmag != zmag2:
    np.sys.exit('Inconsistent zmag: %7.4g, %7.4g'%(zmag, zmag2))
if psisep2 != psisep:
    np.sys.exit('Inconsistent psisep: %7.4g, %7.4g'%(psisep, psisep2))
F = np.empty(nw, dtype=float)
p = np.empty(nw, dtype=float)
ffprime = np.empty(nw, dtype=float)
pprime = np.empty(nw, dtype=float)
qpsi = np.empty(nw, dtype=float)
psirz_1d = np.empty(nw*nh, dtype=float)
start_line = 5
lines = range(nw//5)
if nw%5 != 0:
    lines = range(nw//5 + 1)
for i in lines:
    n_entries = len(eqdsk[i + start_line])//entrylength
    F[i*5:i*5 + n_entries] = [float(eqdsk[i + start_line][j*entrylength:(j + 1)*entrylength]) for
                              j in range(n_entries)]
start_line = i + start_line + 1

for i in lines:
    n_entries = len(eqdsk[i + start_line])//entrylength
    p[i*5:i*5 + n_entries] = [float(eqdsk[i + start_line][j*entrylength:(j + 1)*entrylength]) for
                              j in range(n_entries)]
start_line = i + start_line + 1

for i in lines:
    n_entries = len(eqdsk[i + start_line])//entrylength
    ffprime[i*5:i*5 + n_entries] = [
        float(eqdsk[i + start_line][j*entrylength:(j + 1)*entrylength]) for j in
        range(n_entries)]
start_line = i + start_line + 1

for i in lines:
    n_entries = len(eqdsk[i + start_line])//entrylength
    pprime[i*5:i*5 + n_entries] = [
        float(eqdsk[i + start_line][j*entrylength:(j + 1)*entrylength]) for j in
        range(n_entries)]
start_line = i + start_line + 1

lines_twod = range(nw*nh//5)
if nw*nh%5 != 0:
    lines_twod = range(nw*nh//5 + 1)
for i in lines_twod:
    n_entries = len(eqdsk[i + start_line])//entrylength
    psirz_1d[i*5:i*5 + n_entries] = [
        float(eqdsk[i + start_line][j*entrylength:(j + 1)*entrylength]) for j in
        range(n_entries)]
start_line = i + start_line + 1
psirz = psirz_1d.reshape(nh, nw)

for i in lines:
    n_entries = len(eqdsk[i + start_line])//entrylength
    qpsi[i*5:i*5 + n_entries] = [float(eqdsk[i + start_line][j*entrylength:(j + 1)*entrylength])
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

# create 3rd order 2D spline representation of Psi(R,Z)
interpol_order = 3
psi_spl = RectBivariateSpline(zgrid, rgrid, psirz, kx=interpol_order, ky=interpol_order)

# linear grid of psi, on which all 1D fields are defined
linpsi = np.linspace(psiax, psisep, nw)

### compute rho_tor (rho_pol) grid
# create rho_tor grid
# method=0 uses a 1d integration of the qprofile
# method=1 uses a 2d integration of the toroidal field (more reliable, but slow for now)
psi_fine,phi_fine,rho_tor_fine=create_toroidal_flux(linpsi,F,psi_spl,qpsi,nw,nh,rmag,zmag,rmin,rdim,zmid,zdim,tfmethod)


# spline of q(linpsi) for rhotor grid
q_spl_psi = UnivariateSpline(linpsi, qpsi, k=interpol_order, s=1e-5)
q_fine = q_spl_psi(psi_fine)

rho_tor_spl = UnivariateSpline(psi_fine, rho_tor_fine, k=interpol_order, s=1e-5)
rho_tor = np.empty(nw, dtype=float)
for i in range(nw):
    rho_tor[i] = rho_tor_spl(linpsi[i])

if write_rhotor_rhopol_conversion:
    rt_rp_filename = 'rt_rp_%s'%filename
    rt_rp_file = open(rt_rp_filename, 'w')
    rt_rp_file.write('# rho_tor          rho_pol      q\n')
    for i in range(len(psi_fine)):
        rho_pol = np.sqrt((psi_fine[i] - psiax)/(psisep - psiax))
        rt_rp_file.write('%16.8e %16.8e %16.8e\n'%(rho_tor_fine[i], rho_pol, q_fine[i]))
    rt_rp_file.close()
    exit('\nWrote rhotor/rhopol conversion to %s.'%rt_rp_filename)


theta_arr = np.linspace(-np.pi, np.pi, ntheta)

### Analyze flux surfaces for linpsi values
print('Finding flux surface shapes...')
R = np.empty((nw, ntheta), dtype=float)
Z = np.empty((nw, ntheta), dtype=float)
dr = rdim*np.cos(theta_arr)
dz = rdim*np.sin(theta_arr)
for j in range(len(theta_arr)):
    stdout.write('\r Finished %4.1f%%.'%(j*100./(ntheta - 1)))
    stdout.flush()
    theta = theta_arr[j]
    r_pol = np.linspace(rmag, rmag + dr[j], nr)
    z_pol = np.linspace(zmag, zmag + dz[j], nr)
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
# find average elevation for all flux surfaces (needed to determine ravg)
Z_avg = np.empty(nw, dtype=float)
ds = np.empty(ntheta, dtype=float)
for i in range(1, nw):
    ds[1:ntheta - 1] = 0.5*np.sqrt(
            (R[i, 2:ntheta] - R[i, 0:ntheta - 2]) ** 2 + (Z[i, 2:ntheta] - Z[i, 0:ntheta - 2]) ** 2)
    ds[0] = 0.5*np.sqrt((R[i, 1] - R[i, -1]) ** 2 + (Z[i, 1] - Z[i, -1]) ** 2)
    ds[-1] = 0.5*np.sqrt((R[i, 0] - R[i, -2]) ** 2 + (Z[i, 0] - Z[i, -2]) ** 2)
    Z_avg[i] = np.average(Z[i, :], weights=ds)
# Treat the magnetic axis separately as no average is required and ds==0
Z_avg[0] = Z[0, 0]

# find R0,Z0 for all flux surfaces
R0 = np.empty(nw, dtype=float)
R0_maxmin = np.empty(nw, dtype=float)
R0[0] = rmag
R0_maxmin[0] = rmag
r_avg = np.empty(nw, dtype=float)
r_avg[0] = 0.
r_maxmin = np.empty(nw, dtype=float)
r_maxmin[0] = 0.
for i in range(1, nw):
    stdout.write('\r Finished %4.1f%%.'%(i*100./(nw - 1)))
    stdout.flush()
    R_array = R[i, ntheta//4:3*ntheta//4]
    Z_array = Z[i, ntheta//4:3*ntheta//4]
    # low field side
    Z_int = interp1d(Z_array, range(ntheta//2), kind=interpol_order)
    ind_Zavg = Z_int(Z_avg[i])
    R_int = interp1d(range(ntheta//2), R_array, kind=interpol_order)
    R_out = R_int(ind_Zavg)
    R_max = np.amax(R_array)
    # high field side
    R_array = np.roll(R[i, :-1], ntheta//2)[ntheta//4:3*ntheta//4]
    Z_array = np.roll(Z[i, :-1], ntheta//2)[ntheta//4:3*ntheta//4]

    # have to use negative Z_array here to have increasing order
    Z_int = interp1d(-Z_array, range(ntheta//2), kind=interpol_order)
    # again negative
    ind_Zavg = Z_int(-Z_avg[i])

    R_int = interp1d(range(ntheta//2), R_array, kind=interpol_order)
    R_in = R_int(ind_Zavg)
    R_min = np.amin(R_array)
    R0[i] = 0.5*(R_out + R_in)
    r_avg[i] = 0.5*(R_out - R_in)
    r_maxmin[i] = 0.5*(R_max - R_min)
    R0_maxmin[i] = 0.5*(R_max + R_min)

# GENE x coord (either ravg or r_maxmin)
xgene = r_maxmin

if use_r_a:
    print('\nExamine {0:3d} flux surfaces around position r/a={1:7.4g}...'.format(pw, radpos))
    r_a = radpos
    # find psi index of interest (for the specified r/a position)
    poi_ind = find(radpos, xgene/xgene[-1])
    linpsi_spl = UnivariateSpline(xgene, linpsi, k=interpol_order, s=0)
    psipos = linpsi_spl(radpos*xgene[-1])
else:
    print('\nExamine {0:3d} flux surfaces around position rho_tor={1:7.4g}...'.format(pw, radpos))
    # find nearest psi index of interest (for the specified rho_tor position)
    poi_ind = find(radpos, rho_tor)
    #determine psi value for user-defined rho_tor value
    psi_rho_spl = UnivariateSpline(rho_tor_fine, psi_fine, k=interpol_order, s=1e-8)
    psipos = psi_rho_spl(radpos)

#stencil involving the psi points on the original linpsi grid next to user-defined psipos
psi_stencil = range(poi_ind - pw//2, poi_ind + pw//2)
if psi_stencil[0] < 1:
    psi_stencil = [psi_stencil[i] + 1 - psi_stencil[0] for i in range(len(psi_stencil))]
if psi_stencil[-1] > nw - 1:
    psi_stencil = [psi_stencil[i] - (psi_stencil[-1] - nw + 1) for i in range(len(psi_stencil))]


#Determine flux surface at psipos
Rtarget = np.empty(ntheta, dtype=float)
Ztarget = np.empty(ntheta, dtype=float)
for j in range(len(theta_arr)):
    theta = theta_arr[j]
    r_pol = np.linspace(rmag, rmag + dr[j], nr)
    z_pol = np.linspace(zmag, zmag + dz[j], nr)
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
    Z_int = interp1d(psi_rad[:indsep], z_pol[:indsep], kind=interpol_order)

    Rtarget[j] = R_int(psipos)
    Ztarget[j] = Z_int(psipos)

#target radius
r_maxmin_target = 0.5*(np.amax(Rtarget)-np.amin(Rtarget))

#target center
Rtarget_maxmin = 0.5 * (np.amin(Rtarget) + np.amax(Rtarget))
Ztarget_maxmin = 0.5 * (np.amin(Ztarget) + np.amax(Ztarget))

r = r_maxmin_target
r_a = r/xgene[-1]

xgene_spl = UnivariateSpline(linpsi, xgene, k=interpol_order, s=1e-5)
ravg_spl = UnivariateSpline(linpsi, r_avg, k=interpol_order, s=1e-5)
rmaxmin_spl = UnivariateSpline(linpsi, r_maxmin, k=interpol_order, s=1e-5)

q_spl = UnivariateSpline(xgene, qpsi, k=interpol_order, s=1e-5)
F_spl = UnivariateSpline(xgene, F, k=interpol_order, s=1e-8)
psi_N = (psipos - psiax)/(psisep - psiax)
#R0_pos = float(R0_spl(r))
#Z0_pos = float(Z0_spl(r))
F_pos = float(F_spl(r))



print('Coordinates: r={0:8.5g}, r/R0={1:8.5g}, r/a={2:8.5g}, psi={3:8.5g}, psi_N={4:8.5g}, rho_tor={5:8.5g}'
      .format(r, r/Rtarget_maxmin, r_a, psipos, psi_N, float(rho_tor_spl(psipos))))

Bref_efit = np.abs(F[0]/R0_maxmin[0])
Lref_efit = np.sqrt(2*np.abs(phi_fine[-1])/Bref_efit)

Bref_arbmiller = float(F_spl(r))/Rtarget_maxmin
Lref = Rtarget_maxmin


amhd = np.empty(pw, dtype=float)
for i in psi_stencil:
    imod = i - psi_stencil[0]
    amhd[imod] = -qpsi[i]**2*R0_maxmin[i]*pprime[i]*8.0*np.pi*1e-7/Bref_arbmiller**2/ \
                 xgene_spl.derivatives(linpsi[i])[1]

amhd_spl = UnivariateSpline(xgene[psi_stencil], amhd, k=interpol_order, s=1e-8)
Zavg_spl = UnivariateSpline(r_avg, Z_avg, k=interpol_order, s=1e-8)


# COS/SIN series representation of flux surfaces using r=r_maxmin definition
# @todo: need to interpolate to exact flux surface - otherwise we may parametrize
# a neighbouring one instead

theta_grid = (np.arange(0.0,ntheta,1.0))*2.0*np.pi/ntheta

#flux surface distance
aN_thgrid = np.empty((pw, ntheta), dtype=float)
#print(np.shape(theta_grid),np.shape(aN_thgrid))
for i in psi_stencil:
    imod = i - psi_stencil[0]
    res = get_flux_surface_distance(Rtarget_maxmin,Ztarget_maxmin,R[i,0:ntheta-1],Z[i,0:ntheta-1])
    aN_nonequidist = res["aN"]/Lref
    aN_thgrid[imod,:] = interpolate_periodic(res["theta_sort"],theta_grid,aN_nonequidist)

res = get_flux_surface_distance(Rtarget_maxmin,Ztarget_maxmin,Rtarget[0:ntheta-1],Ztarget[0:ntheta-1])
aN_nonequidist = res["aN"]/Lref
aN_thgrid_target = interpolate_periodic(res["theta_sort"],theta_grid,aN_nonequidist)

aN_intrp = np.empty(ntheta, dtype=float)
aN_dr = np.empty(ntheta, dtype=float)
for iz in range(0,ntheta):
    aN_spl = UnivariateSpline(r_maxmin[psi_stencil], aN_thgrid[:,iz], k=interpol_order, s=1e-5)
    aN_intrp[iz] = aN_spl(r_maxmin_target)
    aN_dr[iz] = aN_spl.derivatives(r_maxmin_target)[1]*Lref

fft_aN = np.fft.fft(aN_thgrid_target)
fft_aN_intrp = np.fft.fft(aN_intrp)
fft_aN_dr = np.fft.fft(aN_dr)

#keep up to nth_kept modes (typically around 30)
nth_kept = 30
cN = 2.0*np.real(fft_aN)[0:nth_kept]/ntheta
cN[0] *= 0.5
sN = -2.0*np.imag(fft_aN)[0:nth_kept]/ntheta
sN[0] *= 0.5

cN2 = 2.0*np.real(fft_aN)[0:nth_kept]/ntheta
cN2[0] *= 0.5
sN2 = -2.0*np.imag(fft_aN)[0:nth_kept]/ntheta
sN2[0] *= 0.5


cN_dr = 2.0*np.real(fft_aN_dr)[0:nth_kept]/ntheta
cN_dr[0] *= 0.5
sN_dr = -2.0*np.imag(fft_aN_dr)[0:nth_kept]/ntheta
sN_dr[0] *= 0.5

#Write information to screen
print('\n\nShaping parameters for flux surface r={0}, r/a={1}:'.format(r, r_a))
print('Copy the following block into a GENE parameters file:\n')

print("&geometry")
print("magn_geometry = 'miller_general'")
print('q0      = {0}'.format(float(q_spl(r))))
print('shat    = {0} !(defined as r/q*dq/dr)'.format(
        r/float(q_spl(r))*q_spl.derivative()(r)))
print('trpeps  = {0}'.format(r/Rtarget_maxmin))
print('minor_r = {0}'.format(xgene[-1]/Lref))
print("major_R = {0}".format(Rtarget_maxmin/Lref))
print("major_Z = {0}".format(Ztarget_maxmin/Lref))
print('amhd    = {0}'.format(float(amhd_spl(r))))
print("")
print("cN_m = ",end="")
for mode in cN:
    print ("{0:15.8e},".format(mode),end="")
print("")
print ("sN_m = ",end="")
for mode in sN:
    print ("{0:15.8e},".format(mode),end="")
print("")
print ("cNdr_m = ",end="")
for mode in cN_dr:
    print ("{0:15.8e},".format(mode),end="")
print("")
print ("sNdr_m = ",end="")
for mode in sN_dr:
    print ("{0:15.8e},".format(mode),end="")
print("")
print("")
print("!adjust the following signs:")
print("sign_Bt_CW = 1")
print("sign_Ip_CW = 1")
print("")

print("rhostar = -1")
print('/\n')
print('\nAdditional information:')
print('&units')
print('Lref        = {0}'.format(Lref))
print('Bref        = {0}'.format(Bref_arbmiller))
print('/\n')
# minor radius defined by 0.5(R_max-R_min), where R_max and R_min can have any elevation
print('a_maxmin    = {0}'.format(r_maxmin[-1]))
# minor radius at average flux surface elevation (GYRO definition)
print('a (avg elev)= {0}'.format(r_avg[-1]))
#print('R0          = {0}'.format(R0_pos))
#print('Z0          = {0}'.format(float(Zavg_spl(r))))
print('Lref_efit   = {0}'.format(Lref_efit))
print('Bref_efit   = {0}'.format(Bref_efit))
print('B_unit(GYRO)= {0}'.format(q_spl(r)/r/xgene_spl.derivatives(psipos)[1]))
print('dpsi/dr     = {0}'.format(1./xgene_spl.derivatives(psipos)[1]))
print('dr_maxmin/dr = {0}'.format(rmaxmin_spl.derivatives(psipos)[1]/xgene_spl.derivatives(psipos)[1]))
print('dr_avg/dr    = {0}'.format(ravg_spl.derivatives(psipos)[1]/xgene_spl.derivatives(psipos)[1]))
print('drho_tor/dr  = {0}'.format(rho_tor_spl.derivatives(psipos)[1]/xgene_spl.derivatives(psipos)[1]))
print('Gradient conversion omt(rho_tor) -> Lref/LT; factor = {0:9.5g}'.format(
        Lref*(rho_tor_spl.derivatives(psipos)[1]/xgene_spl.derivatives(psipos)[1])))



fig = plt.figure()
plt.rcParams.update({'axes.titlesize': 'small',
                     'axes.labelsize': 'small',
                     'xtick.labelsize': 'small',
                     'ytick.labelsize': 'small',
                     'legend.fontsize': 'small'})

########## Plot (R,Z) contours ##############

ax1 = fig.add_subplot(1, 2, 1)
ax1.set_title('Flux surface shape (selected $\Psi$ + stencil)')
for i in psi_stencil:
    ax1.plot(R[i,:],Z[i,:],color='lightgray',linestyle='solid',linewidth=1)

#ax1.plot(Rtarget,Ztarget, label="exact")
#plt.plot(R_sym[ind], Z_sym[ind], 'r-', lw=2, label='symmetrized')
#flux surface from Fourier coefficients
nind = np.arange(0.0,nth_kept,1.0)
a_check = []
a_check2 = []  #aN interpolated from neighbouring psi points
a_dr_check = []
ntheta_check = 500
theta_check = (np.arange(0.0,ntheta_check,1.0))*2.0*np.pi/ntheta_check

for iz in range(0,ntheta_check):
    a_check += [np.sum(cN*np.cos(nind*theta_check[iz])+
                       sN*np.sin(nind*theta_check[iz]))]
    a_check2 += [np.sum(cN2*np.cos(nind*theta_check[iz])+
                       sN2*np.sin(nind*theta_check[iz]))]
theta_check = np.append(theta_check,2.0*np.pi)
a_check = np.append(a_check,a_check[0])
R_check = Rtarget_maxmin+a_check*np.cos(theta_check)*Lref
Z_check = Ztarget_maxmin+a_check*np.sin(theta_check)*Lref
ax1.plot(R_check,Z_check,label="$\Psi$={0:8.5g}".format(psipos))

a_check2 = np.append(a_check2,a_check2[0])
R_check2 = Rtarget_maxmin+a_check2*np.cos(theta_check)*Lref
Z_check2 = Ztarget_maxmin+a_check2*np.sin(theta_check)*Lref
#ax1.plot(R_check2,Z_check2,label="aN_intrp",linestyle='dashed')

ax1.plot(Rtarget_maxmin,Ztarget_maxmin,marker='.',color='black')

ax1.set_xlabel('$R$/m')
ax1.set_ylabel('$Z$/m')
ax1.set_aspect('equal')
ax1.legend(loc=10)

############ plot q profile ###################
ax2 = fig.add_subplot(2, 2, 2)
ax2.plot(xgene[psi_stencil]/xgene[-1], qpsi[psi_stencil],marker='o',linestyle='none',color='lightgray')
xgene_finesten = xgene[psi_stencil[0]]+np.arange(0,pw*5+1,1)*(xgene[psi_stencil[-1]]-xgene[psi_stencil[0]])/(5*pw)
ax2.plot(xgene_finesten/xgene[-1],q_spl(xgene_finesten))
ax2.plot(r_a,q_spl(r),marker='x',color='black')
ax2.set_title('Safety factor')
ax2.set_xlabel(r'$r_{maxmin}/a$')
ax2.set_ylabel(r'$q$')
ax2.axvline(r_a, 0, 1, ls='--', color='k',lw=1)
#ax2a = ax2.twiny()
#ax2a.set_xlim(ax2.get_xlim())
#print(ax2.get_xticks())
#ax2a.set_xticks(ax2.get_xticks())
#ax2a.set_xticklabels(ax2.get_xticks()/xgene[-1])
#ax2a.plot(rho_tor[psi_stencil], qpsi[psi_stencil])
#ax2a.cla()

############ plot shat profile ###################
ax3 = fig.add_subplot(2, 2, 4)
shat = xgene[psi_stencil]/qpsi[psi_stencil]*q_spl.derivative()(xgene[psi_stencil])
ax3.plot(xgene[psi_stencil]/xgene[-1],shat,marker='o',color='lightgray',linestyle='none')
shat = xgene_finesten/q_spl(xgene_finesten)*q_spl.derivative()(xgene_finesten)
ax3.plot(xgene_finesten/xgene[-1],shat)
ax3.plot(r_a,r/float(q_spl(r))*q_spl.derivative()(r),marker='x',color='black')
ax3.set_title('Magnetic shear')
ax3.set_xlabel(r'$r_{maxmin}/a$')
ax3.set_ylabel(r'$\hat{s}$')
ax3.axvline(r_a, 0, 1, ls='--', color='k',lw=1)

if show_plots:
    fig.tight_layout()
    plt.show()
