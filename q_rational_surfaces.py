#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import optparse as op
import matplotlib.pyplot as plt
from read_iterdb import *
from read_write_geometry import *
import math

parser=op.OptionParser(description='Plots rational surfaces for a given toroidal mode number n')

options,args = parser.parse_args()
if len(args)!=5:
    exit("""
Please include toroidal mode number, iterdb file, and geomfile (tracer_efit or gene, etc), lx_a, and x0.
    \n""")

n0 = int(float(args[0]))
idb_file = args[1]
geomfile = args[2]
x0 = float(args[3])
lx_a = float(args[4])

rhotmin = x0-lx_a/2.0
rhotmax = x0+lx_a/2.0

rhot_idb,profs_idb,units_idb = read_iterdb(idb_file)
gpars,geometry = read_geometry_global(geomfile)

rhot = np.linspace(rhotmin,rhotmax,len(geometry['q']))

#plt.plot(rhot_idb['TE'],profs_idb['TE'])
#plt.xlabel('rhot')
#plt.ylabel('Te')
#plt.show()

#plt.plot(rhot,geometry['q'])
#plt.ylabel('q')
#plt.show()

#Find rational q surfaces
qmin = np.min(geometry['q'])
qmax = np.max(geometry['q'])
mmin = math.ceil(qmin*n0)
mmax = math.floor(qmax*n0)
mnums = np.arange(mmin,mmax+1)
print "mnums",mnums
qrats = mnums/float(n0)
print "qrats",qrats

plt.figure(figsize=(6,6))
plt.subplot(2,1,1)
plt.plot(rhot_idb['TE'],profs_idb['TE'])
for i in range(len(qrats)):
    ix = np.argmin(abs(geometry['q']-qrats[i])) 
    plt.axvline(rhot[ix],color='black')
plt.xlabel('rhot')
plt.ylabel('Te')
ix = np.argmin(abs(rhot_idb['TE']-rhotmin))
plt.axis([rhotmin,rhotmax,0,profs_idb['TE'][ix]*1.2])
plt.subplot(2,1,2)
plt.plot(rhot,geometry['q'])
plt.ylabel('q')
plt.xlabel('rhot')
plt.axis([rhot[0],rhot[-1],geometry['q'][0]*0.9,geometry['q'][-1]*1.1])
for i in range(len(qrats)):
    ix = np.argmin(abs(geometry['q']-qrats[i])) 
    plt.axvline(rhot[ix],color='black')
plt.show()


