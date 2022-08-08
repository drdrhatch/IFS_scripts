#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import optparse as op

parser=op.OptionParser(description='Plots neoclassical.')
parser.add_option('--log','-l',action='store_true',dest="plot_log",help = 'Plot on log scale (semilog).',default=False)
options,args=parser.parse_args()


if len(args)!=1:
    exit("""
Please enter file name\n""")

filename = args[0]
plot_log = options.plot_log
print("plot_log",plot_log)

f=open(filename,'r')
data = f.read()

lines = data.split('\n')


keep_going = True
i=2
while keep_going:
    ls = lines[i].split()
    if len(ls) == 1:
        num_spec = i-2
        keep_going = False
    i += 1

print("Number of species:",num_spec)

time = np.empty(0)
# delete header
lines = np.delete(lines,0)
temp = lines[0::num_spec+1]
for i in range(len(temp)):
    if temp[i]:
        time = np.append(time,float(temp[i]))
ntime = len(time)
fluxes = np.empty((ntime,4,num_spec))
for i in range(num_spec):
    temp  = lines[i+1::num_spec+1]
    for j in range(len(temp)):
        if temp[j]:
            ls = np.array(temp[j].split())
            ls = ls.astype(float)
            fluxes[j,:,i] = ls[:]


fig = plt.figure(figsize = (6.0,6.0))
#plt.suptitle(r'$Q_{NL}/Q_{QLstd}$')
plt.subplot2grid((2,2),(0,0))
for i in range(num_spec):
    if plot_log:
        plt.semilogy(time,fluxes[:,0,i],label='species '+str(i))
    else:
        plt.plot(time,fluxes[:,0,i],label='species '+str(i))
plt.legend(loc = 'upper left')
plt.ylabel(r'$\Gamma / \Gamma_{GB}$',size = 14)
plt.xlabel(r'$t(L_{ref}/c_{ref})$',size = 14)
plt.subplot2grid((2,2),(0,1))
for i in range(num_spec):
    if plot_log:
        plt.semilogy(time,fluxes[:,1,i],label='species '+str(i))
    else:
        plt.plot(time,fluxes[:,1,i],label='species '+str(i))
plt.legend(loc = 'upper left')
plt.ylabel(r'$Q / Q_{GB}$',size = 14)
plt.xlabel(r'$t(L_{ref}/c_{ref})$',size = 14)
plt.subplot2grid((2,2),(1,0))
for i in range(num_spec):
    if plot_log:
        plt.semilogy(time,fluxes[:,2,i],label='species '+str(i))
    else:
        plt.plot(time,fluxes[:,2,i],label='species '+str(i))
plt.legend(loc = 'upper left')
plt.ylabel(r'$\Pi / \Pi_{GB}$',size = 14)
plt.xlabel(r'$t(L_{ref}/c_{ref})$',size = 14)
plt.subplot2grid((2,2),(1,1))
for i in range(num_spec):
    if plot_log:
        plt.semilogy(time,fluxes[:,3,i],label='species '+str(i))
    else:
        plt.plot(time,fluxes[:,3,i],label='species '+str(i))
plt.legend(loc = 'upper left')
plt.ylabel(r'$J/(n_{ref}c_{ref}B_{ref}\rho^*_{ref})$',size = 14)
plt.xlabel(r'$t(L_{ref}/c_{ref})$',size = 14)
plt.tight_layout()
plt.show()
