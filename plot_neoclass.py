#!/usr/bin/env python

# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import optparse as op

parser=op.OptionParser(description='Plot neoclassical.')
options,args=parser.parse_args()

if len(args)!=1:
    exit("""
Please enter file name\n""")

filename = args[0]

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

print ("Number of species:",num_spec )
    
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
            ls = ls.astype(np.float)
            fluxes[j,:,i] = ls[:]

for i in range(num_spec):
    plt.semilogy(time,fluxes[:,0,i],label='species '+str(i),basey=10)
plt.legend(loc = 'upper left')
plt.ylabel(r'$\Gamma / \Gamma_{GB}$',size = 14)
plt.xlabel(r'$t(L_{ref}/c_{ref})$',size = 14)
plt.show()

for i in range(num_spec):
    plt.semilogy(time,fluxes[:,1,i],label='species '+str(i),basey=10)
plt.legend(loc = 'upper left')
plt.ylabel(r'$Q / Q_{GB}$',size = 14)
plt.xlabel(r'$t(L_{ref}/c_{ref})$',size = 14)
plt.show()

for i in range(num_spec):
    plt.semilogy(time,fluxes[:,1,i],label='species '+str(i),basey=10)
plt.legend(loc = 'upper left')
plt.ylabel(r'$\Pi / \Pi_{GB}$',size = 14)
plt.xlabel(r'$t(L_{ref}/c_{ref})$',size = 14)
plt.show()

for i in range(num_spec):
    plt.semilogy(time,fluxes[:,1,i],label='species '+str(i),basey=10)
plt.legend(loc = 'upper left')
plt.ylabel(r'$J/(n_{ref}c_{ref}B_{ref}\rho^*_{ref})$',size = 14)
plt.xlabel(r'$t(L_{ref}/c_{ref})$',size = 14)
plt.show()




