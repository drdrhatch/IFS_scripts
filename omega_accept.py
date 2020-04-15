#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from ParIO import *
import os
import numpy as np

suffix = sys.argv[1]
fn = open('new_omega_'+suffix,'r')
fl = fn.readlines()
vals_f = fl[-1].split()
fn.close()

omega_gene_file = 'omega_'+suffix
if os.path.isfile(omega_gene_file):
    omega_gene = np.genfromtxt(omega_gene_file)
    print( "GENE calculation:")
    print(( '  ' + 'kymin' + '      ' + 'gamma' + '    ' + 'omega' ))
    print( omega_gene)

par = Parameters()
par.Read_Pars('parameters_'+suffix)
pars = par.pardict
print( "New calculation:")
print(( '  ' + 'kymin' + '      ' + 'gamma' + '    ' + 'omega' + 'std_gamma' + 'std_omega'))
print(( str(pars['kymin']) + '    ' + vals_f[1] + '    ' + vals_f[2] + '    ' + vals_f[3] + '    ' + vals_f[4]))
selection = int(float(eval(input('How to proceed:\n1. Accept calculation \n2. Manually enter gamma and omega \n3. Don\'t output anything\n'))))
if selection == 1:
    f=open('omega_'+suffix,'w')
    f.write(str(pars['kymin'])+'    '+vals_f[1]+'    '+vals_f[2]+'\n')
    f.close()
elif selection == 2:
#    gam_avg = float(input('Enter gamma: '))
#    om_avg = float(input('Enter omega: '))
#    f=open('omega'+suffix,'w')
#    f.write(str(pars['kymin'])+'    '+str(gam_avg)+'    '+str(om_avg)+'\n')
    f.close()

