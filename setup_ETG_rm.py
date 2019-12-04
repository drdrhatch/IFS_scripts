#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from subprocess import call
import os
from interp import *
from finite_differences import *
import sys
from ParIO import * 
from get_nrg import *
import optparse as op
import csv
import matplotlib.pyplot as plt

#####Setup
ky_scan_string = "10, 20, 30, 40, 50, 60, 80, 120, 160, 240"
num_kxcenter = 3
template_dir = '/global/u2/d/drhatch/gene-dev/prob_ETG_template'
GENE_dir = '/global/u2/d/drhatch/gene-dev/'

parser=op.OptionParser(description='Synthesizes data from nonlinear ETG simulation and sets up corresponding linear ETG run.  Run this script in the ouput directory of the NL simulation.')
#parser.add_option('--time','-t',type = 'float',action='store',dest="time0",help = 'Time to plot mode structure.',default=-1)
options,args=parser.parse_args()
if len(args)!=1:
    exit("""
Please include run number as argument (e.g., 0001)."
    \n""")
suffix = args[0]

if 'dat' in suffix:
   suffix = '.dat'
elif '_' not in suffix:
   suffix = '_'+suffix

par = Parameters()
par.Read_Pars('parameters'+suffix)
pars = par.pardict

print pars
if pars['n_spec'] == 1:
   time,nrg = get_nrg0(suffix,nspec=1)
elif pars['n_spec'] == 2:
   time,nrg_ion,nrg = get_nrg0(suffix,nspec=2)
   dummy = raw_input("Assuming electrons are second species (press any key to continue).")
else:
   print "Can only handle 1 or 2 species right now."
   exit()
   

plt.plot(time,nrg[:,6])
plt.xlabel('t(cs/a)')
plt.ylabel('Q/QGB')
plt.show()

start_time = float(raw_input('Enter start time for average:'))
start_index = np.argmin(abs(time-start_time))

Qavg = np.sum(nrg[start_index:,6])/(len(time)-start_index-1)
Q_CV = np.std(nrg[start_index:,6])/Qavg
print "Qavg",Qavg
print "Q_CV",Q_CV

summary = {}
summary['Qavg'] = Qavg
summary['Q_CV'] = Q_CV
summary['omt'] = pars['omt1']
summary['omn'] = pars['omn1']
summary['nx0'] = pars['nx0']
summary['nky0'] = pars['nky0']
summary['x0'] = pars['x0']
summary['lx'] = pars['lx']
summary['kymin'] = pars['kymin']
summary['beta'] = pars['beta']
summary['debye2'] = pars['debye2']
summary['n_spec'] = pars['n_spec']
if 'ExBrate' in pars:
   summary['ExBrate'] = pars['ExBrate']
else:
   summary['ExBrate'] = 0.0
summary['q0'] = pars['q0']
summary['shat'] = pars['shat']
summary['geomfile'] = pars['geomfile']
if 'edge_opt' in pars:
    summary['edge_opt'] = pars['edge_opt']
else:
    summary['edge_opt'] = 0
summary['minor_r'] = pars['minor_r']
summary['major_R'] = pars['major_R']
summary['Bref'] = pars['Bref']
summary['Tref'] = pars['Tref']
summary['nref'] = pars['nref']
summary['Lref'] = pars['Lref']
summary['mref'] = pars['mref']
summary['nu_ei'] = pars['nu_ei']
if 'nustar_e' in pars:
    summary['nustar_e'] = pars['nustar_e']

w = csv.writer(open('summary'+str(suffix)+'.csv','w'))
for key, val in summary.items():
   w.writerow([key, val])

cwd = os.getcwd()
this_dir = cwd.split('/')[-1]


kx_center_scan_string = '  !scanlist: 0.0 '
for i in range(num_kxcenter-1):
    kx_center_scan_string += ', '+str(0.5*(i+1)/float(num_kxcenter-1)*2*np.pi)+'*'+str(pars['shat'])+'*kymin(1)'

print "kx_center_scan_string",kx_center_scan_string

f=open('parameters'+suffix,'r')
parfile=f.read()
f.close()
parfile_split = parfile.split('\n')
for i in range(len(parfile_split)):
    if 'diagdir' in parfile_split[i]:
        parfile_split[i] = 'diagdir = \''+cwd+'\''
    if 'kymin' in parfile_split[i]: 
        parfile_split[i] = 'kymin = 0.05  !scanlist: '+ky_scan_string \
          +'\nkx_center = 0.0 '+kx_center_scan_string
        print "parfile_split[i]",parfile_split[i]
    #if 'kx_center' in parfile_split[i]: 
    #    parfile_split[i] = ''
    if 'n_procs_v' in parfile_split[i]:
        parfile_split[i] = 'n_procs_v = 4' 
    if 'n_procs_sim' in parfile_split[i]:
        parfile_split[i] = 'n_procs_sim = 256' 
    if 'n_procs_w' in parfile_split[i]:
        parfile_split[i] = 'n_procs_w = 8' 
    if 'n_procs_z' in parfile_split[i]:
        parfile_split[i] = 'n_procs_z = 8' 
    if 'nx0' in parfile_split[i]:
        parfile_split[i] = 'nx0 = 7' 
    if 'nky0' in parfile_split[i]:
        parfile_split[i] = 'nky0 = 1' 
    if 'adapt_lx' in parfile_split[i]:
        parfile_split[i] = 'adapt_lx = T'
    if 'n0_global' in parfile_split[i]:
        parfile_split[i] = ''
    if 'read_checkpoint' in parfile_split[i]:
        parfile_split[i] = 'read_checkpoint = F'
    if 'nonlinear' in parfile_split[i]:
        parfile_split[i] = 'nonlinear = F'
    if 'nblocks' in parfile_split[i]:
        parfile_split[i] = ''
    if 'general' in parfile_split[i]:
        parfile_split[i] = '&general \ncalc_dt = T \nomega_prec = 1.0e-1'
    if 'simtimelim' in parfile_split[i]:
        parfile_split[i] = 'simtimelim = 1.5'
    if 'ExBrate' in parfile_split[i]:
        parfile_split[i] = '!'+parfile_split[i]
    if 'nexc' in parfile_split[i]:
        parfile_split[i] = ''
    if 'magn_geometry' in parfile_split[i]:
        print 'magn_geometry',pars['magn_geometry']
        if pars['magn_geometry'] == '\'tracer_efit\'':
            parfile_split[i] = 'magn_geometry = \'gene\''
    if 'geomfile' in parfile_split[i]:
        if pars['magn_geometry'] == '\'tracer_efit\'':
            parfile_split[i] = 'geomfile = \'tracer_efit'+suffix+'\''
            geomfile = 'tracer_efit'+suffix
        elif pars['magn_geometry'] == '\'gene\'':
            parfile_split[i] = 'geomfile = \'gene'+suffix+'\''
            geomfile = 'gene'+suffix
        
parfile_out='\n'.join(parfile_split)
#print 'parfile_out',parfile_out
f=open('parameters'+suffix+'_linear','w')
f.write(parfile_out)
f.close()

probdir = GENE_dir+'prob_'+this_dir+'_ll'+suffix
call(['cp','-r',template_dir,probdir])
call(['cp','parameters'+suffix+'_linear',probdir+'/parameters'])
call(['cp',geomfile,probdir])
##call(['cp',pars['iterdbfile'],probdir])


