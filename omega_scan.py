#!/usr/bin/env python
# -*- coding: utf-8 -*-

from omega_tool import omega_calc
import os
from ParIO import *
import sys
import optparse as op
import numpy as np

parser=op.OptionParser(description='Calculates growth rate and frequncy from linear scan.')
parser.add_option('--alg','-a',action='store',type = int,help = 'Choose 1 for eigenmode average method and 2 for maximum value.',default = 2)
parser.add_option('--single','-s',action='store',type = str,help = 'Analyze only one run in the scan (e.g. 0001).',default = 0)
options,args=parser.parse_args()

alg = options.alg
single = int(options.single)
print( "alg",alg)
print( "single",single)

cwd = os.getcwd()
filelist = []
if float(single) != 0.0:
   suffix = '0000'+str(single)
   suffix = suffix[-4:]
   filelist.append('_'+suffix)
else:
    for filename in os.listdir(cwd):
        if filename.startswith("omega_0"):
            filelist.append(filename[-5:])
    filelist.sort()

for suffix in filelist:
    print( "suffix",suffix)
    omega_calc(suffix,alg = alg)


fn = open('new_omega'+suffix,'r')
fl = fn.readlines()
vals_f = fl[-1].split()
fn.close()

omega_gene_file = 'omega_'+suffix
if os.path.isfile(omega_gene_file):
    omega_gene = np.genfromtxt(omega_gene_file)
    print( "GENE calculation:")
    print( '  ' + 'kymin' + '      ' + 'gamma' + '    ' + 'omega' )
    print( omega_gene)

par = Parameters()
par.Read_Pars('parameters'+suffix)
pars = par.pardict
print( "New calculation:")
print( '  ' + 'kymin' + '      ' + 'gamma' + '    ' + 'omega' + 'std_gamma' + 'std_omega')
print( str(pars['kymin']) + '    ' + vals_f[1] + '    ' + vals_f[2] + '    ' + vals_f[3] + '    ' + vals_f[4])
selection = int(float(input('How to proceed:\n1. Accept calculation \n2. Manually enter gamma and omega \n3. Don\'t output anything\n')))
if selection == 1:
    f=open('omega'+suffix,'w')
    f.write(str(pars['kymin'])+'    '+vals_f[1]+'    '+vals_f[2]+'\n')
    f.close()
elif selection == 2:
#    gam_avg = float(input('Enter gamma: '))
#    om_avg = float(input('Enter omega: '))
#    f=open('omega'+suffix,'w')
#    f.write(str(pars['kymin'])+'    '+str(gam_avg)+'    '+str(om_avg)+'\n')
    f.close()

