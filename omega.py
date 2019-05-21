#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import os
import numpy as np
import sys
from parIOWrapper import *

#title = sys.argv[1]
filelist = []
X = []
Y1 = []
Y2 =[]
cwd = os.getcwd()
for filename in os.listdir(cwd):
    if filename.startswith("omega_0"):
        filelist.append(filename)
filelist.sort()
for fname in filelist:
    data=np.loadtxt(fname)
#    print(data)
    if data[0] > 4:
     X.append(data[0])
     Y1.append(data[1]/60)
    else:
     X.append(data[0])
     Y1.append(data[1])
    Y2.append(data[2])   
 
#plt.title()
plt.plot(X,Y1,'-bo',label= r'$\gamma/\Omega_i$')
#plt.plot(X,Y2,'-ro',label= r'$\omega/\Omega_i$')
plt.legend(loc='upper left')
plt.xlabel(r'$\rho_{tor}$')
#plt.ylabel(r'$\gamma/\Omega_i$')
plt.show()
plt.plot(X,Y2,'-ro',label= r'$\omega/\Omega_i$')
plt.xlabel(r'$\rho_{tor}$')
plt.show()
