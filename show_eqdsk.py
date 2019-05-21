#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to visualize and print g-eqdsk contents
Created on Thu Jul 30 12:55:00 2015

@author: tbg, dtold
"""

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser=argparse.ArgumentParser(description='Plot EQDSK contents.')
parser.add_argument("filename", help="EQDSK filename")
args=parser.parse_args()

try:
    with open(args.filename,'r') as file:
        eqdsk=file.readlines()
except IOError:
    sys.exit("EQDSK file not found")

endline=len(eqdsk)
print('Header: {0:s}'.format(eqdsk[0]))
#set resolutions
nw=int(eqdsk[0].split()[-2])
nh=int(eqdsk[0].split()[-1])
print('Resolution: {0:4d} x {1:4d}'.format(nw,nh))

entrylength=16
try:
    rdim,zdim,rctr,rmin,zmid=([float(eqdsk[1][j*entrylength:(j+1)*entrylength])
                               for j in range(len(eqdsk[1])//entrylength)])
except:
    entrylength=15
    try:
        rdim,zdim,rctr,rmin,zmid=([float(eqdsk[1][j*entrylength:(j+1)*entrylength])
                                   for j in range(len(eqdsk[1])//entrylength)])
    except:
        sys.exit('Error reading EQDSK file, please check format!')

rmag,zmag,psiax,psisep,Bctr=([float(eqdsk[2][j*entrylength:(j+1)*entrylength])
                              for j in range(len(eqdsk[2])//entrylength)])
current,psiax2,dum,rmag2,dum=([float(eqdsk[3][j*entrylength:(j+1)*entrylength])
                               for j in range(len(eqdsk[3])//entrylength)])
zmag2,dum,psisep2,dum,dum=([float(eqdsk[4][j*entrylength:(j+1)*entrylength])
                            for j in range(len(eqdsk[4])//entrylength)])

if rmag!=rmag2:
    sys.exit('Inconsistent rmag: {:7.4g}, {:7.4g}'.format(rmag,rmag2))
if psiax2!=psiax:
    sys.exit('Inconsistent psiax: {:7.4g}, {:7.4g}'.format(psiax,psiax2))
if zmag!=zmag2:
    sys.exit('Inconsistent zmag: {:7.4g}, {:7.4g}'.format(zmag,zmag2))
if psisep2!=psisep:
    sys.exit('Inconsistent psisep: {:7.4g}, {:7.4g}'.format(psisep,psisep2))

print('R of magn. axis [m]:                 {0:6.3f}'.format(rmag))
print('Z of magn. axis [m]:                 {0:6.3f}'.format(zmag))
print('Pol. flux at magn. axis [Weber/rad]: {0:6.3f}'.format(psiax))
print('Pol. flux at separatrix [Weber/rad]: {0:6.3f}'.format(psisep))
print('R of vacuum tor. magn. field [m]:    {0:6.3f}'.format(rctr))
print('Vacuum toroidal magnetic field [T]:  {0:6.3f}'.format(Bctr))
print('Plasma current [MA]:                 {0:6.3f}'.format(current/1E6))

fpol=np.empty(nw,dtype=float)
pres=np.empty(nw,dtype=float)
ffprime=np.empty(nw,dtype=float)
pprime=np.empty(nw,dtype=float)
qpsi=np.empty(nw,dtype=float)
psirz_1d=np.empty(nw*nh,dtype=float)
start_line=5
lines=range(nw//5)
if nw%5!=0:
    lines=range(nw//5+1)
for i in lines:
    n_entries=len(eqdsk[i+start_line])//entrylength
    fpol[i*5:i*5+n_entries]=([float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength])
                              for j in range(n_entries)])
start_line=i+start_line+1
print('Toroidal magn. field on axis [T]:    {0:6.3f}'.format(fpol[0]/rmag))

for i in lines:
    n_entries=len(eqdsk[i+start_line])//entrylength
    pres[i*5:i*5+n_entries]=([float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength])
                              for j in range(n_entries)])
start_line=i+start_line+1

for i in lines:
    n_entries=len(eqdsk[i+start_line])//entrylength
    ffprime[i*5:i*5+n_entries]=([float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength])
                                 for j in range(n_entries)])
start_line=i+start_line+1

for i in lines:
    n_entries=len(eqdsk[i+start_line])//entrylength
    pprime[i*5:i*5+n_entries]=([float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength])
                                for j in range(n_entries)])
start_line=i+start_line+1

lines_twod=range(nw*nh//5)
if nw*nh%5!=0:
    lines_twod=range(nw*nh//5+1)
for i in lines_twod:
    n_entries=len(eqdsk[i+start_line])//entrylength
    psirz_1d[i*5:i*5+n_entries]=([float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength])
                                  for j in range(n_entries)])
start_line=i+start_line+1
psirz=psirz_1d.reshape(nh,nw)

for i in lines:
    n_entries=len(eqdsk[i+start_line])//entrylength
    qpsi[i*5:i*5+n_entries]=([float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength])
                              for j in range(n_entries)])
start_line=i+start_line+1
print('Safety factor q on axis:             {0:6.3f}'.format(qpsi[0]))
print('Safety factor q at separatrix:       {0:6.3f}'.format(qpsi[nw-1]))

#invert sign of psi if necessary to guarantee increasing values for interpolation
if psisep<psiax:
    psirz=-psirz
    ffprime=-ffprime
    pprime=-pprime
    psiax*=-1
    psisep*=-1

#boundary and limiter data if available
if start_line<endline:
    nbbbs=int(eqdsk[start_line].split()[0])
    limitr=int(eqdsk[start_line].split()[1])
    print('Boundary grid points: {0:4d}'.format(nbbbs))
    print('Limiter grid points:  {0:4d}'.format(limitr))
    start_line = start_line+1

    rzbbbs=np.empty(2*nbbbs,dtype=float)
    rbbbs=np.empty(nbbbs+1,dtype=float)
    zbbbs=np.empty(nbbbs+1,dtype=float)

    if (2*nbbbs)%5!=0:
        lines=range((2*nbbbs)//5+1)
    else: lines=range((2*nbbbs)//5)
    for i in lines:
        n_entries=len(eqdsk[i+start_line])//entrylength
        rzbbbs[i*5:i*5+n_entries]=([float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength])
                                    for j in range(n_entries)])
    start_line=i+start_line+1

    rbbbs[0:nbbbs]=rzbbbs[0::2]
    zbbbs[0:nbbbs]=rzbbbs[1::2]
    rbbbs[nbbbs]=rbbbs[0]
    zbbbs[nbbbs]=zbbbs[0]

    rzlim=np.empty(2*limitr,dtype=float)
    rlim=np.empty(limitr+1,dtype=float)
    zlim=np.empty(limitr+1,dtype=float)
    if (2*limitr)%5!=0:
        lines=range((2*limitr)//5+1)
    else:
        lines=range((2*limitr)//5)
    for i in lines:
        n_entries=len(eqdsk[i+start_line])//entrylength
        rzlim[i*5:i*5+n_entries]=([float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength])
                                   for j in range(n_entries)])
    start_line=i+start_line+1
    rlim[0:limitr]=rzlim[0::2]
    zlim[0:limitr]=rzlim[1::2]
    rlim[limitr]=rlim[0]
    zlim[limitr]=zlim[0]
else:
    nbbbs=0
    limitr=0

#further processing
dw=rdim/(nw-1)
dh=zdim/(nh-1)
rgrid=np.array([rmin+i*dw for i in range(nw)])
zgrid=np.array([zmid-zdim/2.+i*dh for i in range(nh)])
#contourf(rgrid,zgrid,psirz,70);gca().set_aspect('equal')
#show()

#create 5th order 2D spline representation of Psi(R,Z)
from scipy.interpolate import RectBivariateSpline as RBS
from scipy.interpolate import UnivariateSpline as US
interpol_order=3
psi_spl=RBS(zgrid,rgrid,psirz,kx=interpol_order,ky=interpol_order)

#linear grid of psi, on which all 1D fields are defined
linpsi=np.linspace(psiax,psisep,nw)
#create rho_tor grid
x_fine=np.linspace(psiax,psisep,nw*10)
phi_fine=np.empty((nw*10),dtype=float)
phi_fine[0]=0.
#spline of q for rhotor grid
q_spl_psi=US(linpsi,qpsi,k=interpol_order,s=1e-5)

for i in range(1,nw*10):
    x=x_fine[:i+1]
    y=q_spl_psi(x)
    phi_fine[i]=np.trapz(y,x)
rho_tor_fine=np.sqrt(phi_fine/phi_fine[-1])
rho_tor_spl=US(x_fine,rho_tor_fine,k=interpol_order,s=1e-5)
rho_tor=np.empty(nw,dtype=float)
for i in range(nw):
    rho_tor[i]=rho_tor_spl(linpsi[i])

nlev=80
fig1=plt.figure()
ax1=fig1.add_subplot(111, aspect='equal')
cf=ax1.contourf(rgrid,zgrid,psirz,nlev)#, cmap='gist_heat_r')
if nbbbs>0:
    ax1.plot(rbbbs,zbbbs)
if limitr>0:
    ax1.plot(rlim,zlim)
ax1.plot(rmag,zmag,'p')
fig1.colorbar(cf)
ax1.set_xlabel(r'$R / m$')
ax1.set_ylabel(r'$Z / m$')
ax1.set_title(r'$\Psi(R,Z)$')

fig2=plt.figure()
ax2_1=fig2.add_subplot(3,1,1)
ax2_1.plot(linpsi,qpsi)
ax2_1.set_xlabel(r'$\Psi$')
ax2_1.set_ylabel(r'$q$')

ax2_2=fig2.add_subplot(3,1,2)
ax2_2.plot(linpsi,pres,label=r'$p$')
ax2_2.plot(linpsi,pprime,label=r"$p^\prime$")
ax2_2.legend()
ax2_2.set_xlabel(r'$\Psi$')

ax2_3=fig2.add_subplot(3,1,3)
ax2_3.plot(linpsi,fpol,label=r'$f_{pol}$')
ax2_3.plot(linpsi,ffprime,label=r"$ff^\prime$")
ax2_3.legend()
ax2_3.set_xlabel(r'$\Psi$')

provide_file=0
if provide_file:
    f=open('rho_tor-p','w')
    f.write('#rho_tor              p[kPa]        q\n')
    for i in range(0,nw):
        f.write('%16.8e %16.8e %16.8e\n' %(rho_tor[i],pres[i],qpsi[i]))
    exit(0)

plt.show()
