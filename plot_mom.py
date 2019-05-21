#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plots: Several columns from nrg file along with useful information from parameters file and outputs figure.
"""

import matplotlib.pyplot as plt
import numpy as np
from get_nrg_x import get_nrg0
from subprocess import call
import optparse as op
from sys import path
GENE_path="/global/homes/h/halfmoon/knl-refactored/tools/python"
path.append(GENE_path)
script_path="/global/homes/h/halfmoon/scripts"
path.append(script_path)
from ParIO import Parameters

parser=op.OptionParser(description='Useful nrg plots.')
parser.add_option('--tmax','-t',action='store',type = 'float')
#parser.add_option('--conv','-c',action='store_const',const=1,help = 'Output rho_tor vs rho_pol')
#parser.add_option('--binfo','-p',action='store_const',const=1,help = 'Ouptut R, psi_pol, B_pol, B_tor')
#parser.add_option('--noplot','-n',action='store_const',const=1,help = 'Suppress plots')
options,args=parser.parse_args()
if len(args)<2:
    exit("""
Arguments must be 1.Case name 2.Run number (e.g. .dat or _1)
    \n""")
case_name = args[0]
run_number = args[1]

print "run_number", run_number
tmax = options.tmax
#print "tmax",tmax

par=Parameters()
par.Read_Pars('parameters'+run_number)
n_spec = int(float(par.pardict['n_spec']))
print 'n_spec',n_spec
nrgcols = int(float(par.pardict['nrgcols']))
rhostar = float(par.pardict['rhostar'])
print 'ncols',nrgcols

time,nrgi,nrge=get_nrg0(run_number,nspec=n_spec,ncols = nrgcols)

if not tmax:
    print "Identify end time for plots."

    plt.plot(time,nrgi[:,0],'x-')
    plt.xlabel(r'$t(L_{ref}/v_{ref})$',size=18)
    plt.ylabel(r'$n^2/(n_0\rho^*)^2$',size=18)
    plt.show()

    tmax = float(raw_input("Enter end time for plot:"))

end_index = np.argmin(abs(time-tmax))
print "time[end_index]",time[end_index]

plt.figure(figsize=(14,8))
ax=plt.subplot(231)
plt.subplots_adjust(wspace = 0.3)
plt.subplots_adjust(hspace = 0.3)
plt.subplots_adjust(right = 0.8)
plt.subplots_adjust(left = 0.1)
plt.title(r'$n^2/(n_0 \rho^*)^2$')
plt.plot(time[:end_index],nrgi[:end_index,0],label='ni')
plt.plot(time[:end_index],nrge[:end_index,0],label='ne')
ax=plt.subplot(232)
plt.title(r'$u_{||}^2/(c_s \rho^*)^2$')
plt.plot(time[:end_index],nrgi[:end_index,1],label='ui')
plt.plot(time[:end_index],nrge[:end_index,1],label='ue')
ax=plt.subplot(233)
plt.title(r'$T_{||}^2/(T_0 \rho^*)^2$')
plt.plot(time[:end_index],nrgi[:end_index,2],label='Tpari')
plt.plot(time[:end_index],nrge[:end_index,2],label='Tpare')
plt.annotate('Case: '+case_name,[0.82,0.9],xycoords='figure fraction')
plt.annotate('nrg'+run_number,[0.82,0.87],xycoords='figure fraction')
nx0 = str(par.pardict['nx0'])
print "nx0",nx0
print "type(nx0)",type(nx0)
nky0 = str(par.pardict['nky0'])
nz0 = str(par.pardict['nz0'])
nv0 = str(par.pardict['nv0'])
nw0 = str(par.pardict['nw0'])
x0 = str(par.pardict['x0'])
lx_a = str(par.pardict['lx_a'])
lx = str(par.pardict['lx'])
rhostar = str(par.pardict['rhostar'])
beta = str(par.pardict['beta'])
coll = str(par.pardict['coll'])
hyp_x = str(par.pardict['hyp_x'])
hyp_y = str(par.pardict['hyp_y'])
hyp_z = str(par.pardict['hyp_z'])
plt.annotate('nx0: '+nx0,[0.82,0.84],xycoords='figure fraction')
plt.annotate('nky0: '+nky0,[0.82,0.81],xycoords='figure fraction')
plt.annotate('nz0: '+nz0,[0.82,0.78],xycoords='figure fraction')
plt.annotate('x0: '+x0,[0.82,0.75],xycoords='figure fraction')
plt.annotate('lx_a: '+lx_a,[0.82,0.72],xycoords='figure fraction')
plt.annotate('lx: '+lx,[0.82,0.69],xycoords='figure fraction')
plt.annotate('rhostar: '+rhostar,[0.82,0.66],xycoords='figure fraction')
plt.annotate('beta: '+beta,[0.82,0.63],xycoords='figure fraction')
plt.annotate('coll: '+coll,[0.82,0.6],xycoords='figure fraction')
plt.annotate('hyp_x: '+hyp_x,[0.82,0.57],xycoords='figure fraction')
plt.annotate('hyp_y: '+hyp_y,[0.82,0.54],xycoords='figure fraction')
plt.annotate('hyp_z: '+hyp_z,[0.82,0.51],xycoords='figure fraction')
diagdir = par.pardict['diagdir']
plt.annotate('diagdir: ',[0.82,0.48],xycoords='figure fraction')
plt.annotate(diagdir[1:-1],[0.78,0.45],xycoords='figure fraction')
ax=plt.subplot(234)
plt.title(r'$T_{\perp}^2/(T_0 \rho^*)^2$')
plt.plot(time[:end_index],nrgi[:end_index,3],label='Tperpi')
plt.plot(time[:end_index],nrge[:end_index,3],label='Tperpe')
ax=plt.subplot(235)
plt.title(r'$Q_{ES}/Q_{GB}$')
plt.plot(time[:end_index],nrgi[:end_index,6],label='QESi')
plt.plot(time[:end_index],nrge[:end_index,6],label='QESe')
plt.xlabel(r'$t(a/c_s)$',size=18)
ax=plt.subplot(236)
plt.title(r'$Q_{EM}/Q_{GB}$')
plt.plot(time[:end_index],nrgi[:end_index,7],label='QEMi')
plt.plot(time[:end_index],nrge[:end_index,7],label='QEMe')
plt.savefig(case_name+run_number+'.png')
plt.close()





