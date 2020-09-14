#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from read_iterdb_file import *
import optparse as op


parser=op.OptionParser(description='Converts from iterdb file to GENE profile files.')

options,args=parser.parse_args()
if len(args)!=1:
    exit("""
Please include the name of the iterdb file as argument."
    \n""")

fiterdb = args[0]

rhot, te, ti, ne, ni, nb, vrot = read_iterdb_file(fiterdb)

plt.plot(rhot,te)
plt.xlabel('rhot')
plt.ylabel('Te(ev)')
plt.show()
plt.plot(rhot,ti)
plt.xlabel('rhot')
plt.ylabel('Ti(ev)')
plt.show()
plt.plot(rhot,ne)
plt.xlabel('rhot')
plt.ylabel('ne(m^-3)')
plt.show()
plt.plot(rhot,ni)
plt.xlabel('rhot')
plt.ylabel('ni(m^-3)')
plt.show()
if nb[0]:
   plt.plot(rhot,nb)
   plt.xlabel('rhot')
   plt.ylabel('nz(m^-3)')
   plt.show()
plt.plot(rhot,vrot)
plt.xlabel('rhot')
plt.ylabel('vrot(rad/s)')
plt.show()

zeros = np.zeros(len(rhot))

f = open('profiles_e','w')
f.write('# 1.rhot, 2.zeros 3.Te(keV) 4.ne(10^19 m^-3)\n#\n')
np.savetxt(f,np.column_stack((rhot,zeros,te/1000.0,ne/1e19)))
f.close()
        
f = open('profiles_i','w')
f.write('# 1.rhot, 2.zeros 3.Ti(keV) 4.ni(10^19 m^-3)\n#\n')
np.savetxt(f,np.column_stack((rhot,zeros,ti/1000.0,ni/1e19)))
f.close()

if nb[0]:
   Z = ((ni-ne)/nb)
   Z = abs(int(round(Z[int(len(ni)/2)])))
   print('Z for impurity:',Z)
   f = open('profiles_z','w')
   f.write('# 1.rhot, 2.zeros 3.Tz(keV) 4.nz(10^19 m^-3)\n#Z = '+str(Z)+'\n')
   np.savetxt(f,np.column_stack((rhot,zeros,ti/1000.0,nb/1e19)))
   f.close()


