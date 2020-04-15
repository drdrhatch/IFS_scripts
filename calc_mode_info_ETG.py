#!/usr/bin/env python
# -*- coding: utf-8 -*-

from subprocess import call
import os

for i in range(100):
    num = '0000'+str(i)
    num = num[-4:]
    print(num)
    if os.path.exists('./parameters_'+num):
        call(['calc_kpar_kperp_omd.py',num])

files = []

base = 'mode_info_'
cat_string = 'cat '
for i in range(100):  
    for j in range(50):
        suffix = '0000'+str(i)
        suffix = suffix[-4:]
        efile = base+str(j)+'_'+suffix
        if os.path.exists('./'+efile):
            print(efile)
            files.append(efile)
            cat_string += efile+' '
            print(efile)

cat_string += '> '+'mode_info_all'

call(cat_string,shell=True)

