#!/usr/bin/env python
# -*- coding: utf-8 -*-

from omega_tool import omega_calc
import os
import optparse as op

parser=op.OptionParser(description='Calculates growth rate and frequncy from linear scan.')
parser.add_option('--alg','-a',action='store',type = int,help = 'Choose 1 for eigenmode average method and 2 for maximum value.',default = 2)
parser.add_option('--single','-s',action='store',type = str,help = 'Analyze only one run in the scan (e.g. 0001).',default = 0)
options,args=parser.parse_args()

alg = options.alg
single = int(options.single)
print(( "alg",alg))
print(( "single",single))

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
    print(( "suffix",suffix))
    omega_calc(suffix,alg = alg)
