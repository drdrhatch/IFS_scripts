#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Finds the rhot value corresponding to a given T from a gene parameters file.
"""
#import matplotlib.pyplot as plt
import optparse as op
import numpy as np
from read_EFIT_file import *

parser=op.OptionParser(description='Finds the rhot value corresponding to a given T from a gene parameters files.  Arguments: Te, gene_profiles_file_name_e' )
options,args=parser.parse_args()
if len(args)!=2:
    exit("""
Please include Te and gene_profiles_file_name_e."
    \n""")

Te = float(args[0])
gene_profiles_e = args[1]
data = np.genfromtxt(gene_profiles_e)
irhot = np.argmin(abs(data[:,2]-Te))
rhot = data[irhot,0]
print("rho_toroidal with Te = "+str(Te)+":\n",rhot)
print("density at = "+str(rhot)+":\n",data[irhot,3])

