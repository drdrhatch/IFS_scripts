#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import re
import ParIO as pario
import fieldlib
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("suffix", help="run number or .dat suffix of output data")
parser.add_argument(
    "--stime", "-s", action="store", default=0, type=float, help="start time for POD"
)
parser.add_argument("--etime", "-e", action="store", help="end time for POD")
parser.add_argument(
    "--time", "-t", dest="time0", action="store", help="end time for POD"
)
args = parser.parse_args()

if re.search("dat$", args.suffix):
    suffix = ".dat"
elif re.search("[0-9]{1,4}$", args.suffix):
    match = re.search(r"([0-9]{1,4})$", args.suffix)
    suffix = "_" + match.group(0).zfill(4)
    print(suffix)
else:
    print("Please enter a valid run number, e.g. .dat or 0123")

par = pario.Parameters()
par.Read_Pars("parameters" + suffix)
pars = par.pardict

field = fieldlib.fieldfile("field" + suffix, pars)
min_time, max_time = field.get_minmaxtime()
stime = args.stime

if args.stime == None:
    stime = max_time
else:
    stime = float(args.stime)

if args.etime == None:
    etime = max_time
else:
    etime = float(args.etime)

print("start time = ", stime)
print("end time = ", etime)

if args.time0 == None:
    time0 = max_time
else:
    time0 = float(args.time0)

time = np.array(field.tfld)
print(time)
print(type(field.tfld))
print(type(time))
print(type(time0))
itime = np.argmin(abs(time - time0))
field.set_time(time[itime])

nx = field.nx
ny = field.ny
nz = field.nz
nexc = pars["nexc"]

print('(nx, ny, nz)  = ',nx,ny,nz)

ntot = nx * nz

phase_fac = (-1) ** nexc

ndz = 2.0 / float(nz)
# zgrid = np.linspace(-nx, nx-dz, ntot)

# for j in ky_list:
j = 1


def kxrange(ky, nx, N):
    hmodes = np.arange(0, nx / 2, N * ky, dtype=np.intc)
    lmodes = np.arange(0, -nx / 2, -N * ky, dtype=np.intc)
    modes = np.union1d(lmodes, hmodes)
    return modes


def phase(ky, N):
    phase = (-1) ** (ky * N)
    return phase


def zrange(kx, nz):
    ncon = (kx.size - 1) // 2
    ncon1 = ncon+1
    zgrid = np.pi * np.linspace(-(ncon+1), ncon+1, (ncon+1) * nz, endpoint=False)
    return zgrid


kx = kxrange(j, nx, nexc)
for i in kxrange(j, nx, nexc):
    phi = np.zeros(kx.size*nz, dtype="complex128")
    # istart1 = (i+nx//2)*nz
    # iend1   = (i+nx//2+1)*nz
    # phi[istart1:iend1] = field.phi()[:,j,i]*phase_fac**(j)

    # istart2 = (nx//2-i-1)*nz
    # iend2   = (nx//2-i)*nz
    # phi[istart2:iend2] = field.phi()[:,j,-(i+1)]*phase_fac**(j)
    # print(type(field.phi()))
    print(i)
    print(field.phi()[:, j, i].size,field.phi3d.size)
    phi[i * nz : (i + 1) * nz] = field.phi()[:, j, i] * phase(j, nexc)
    print(i * nz, (i + 1) * nz, phi.size,phi,)

print(kx)
zgrid = zrange(kx, nz)
print(zgrid.size,phi.size,kx.size)
plt.title(r"$\phi$")
plt.plot(zgrid, np.real(phi), color="red", label=r"$Re[\phi]$")
plt.plot(zgrid, np.imag(phi), color="blue", label=r"$Im[\phi]$")
plt.plot(zgrid, np.abs(phi), color="black", label=r"$|\phi|$")
ax = plt.axis()
plt.legend()
plt.xlabel(r"$z/\pi$", size=18)
plt.show()
