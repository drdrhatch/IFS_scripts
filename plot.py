#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import os
import numpy as np
import sys

f_i = open('profiles_i','r')
i_lines = f_i.readlines()
f_z = open('profiles_z','r')
z_lines = f_z.readlines()
f_i.close()
f_z.close()
X_i = []
Y_i = []
Y_z = []

for line in i_lines:
    X_i.append(line.split(' ')[0])
    Y_i.append(line.split(' ')[3])

for line in z_lines:
    Y_z.append(line.split(' ')[3])

i=0
while i < len(Y_i):
    print(float(Y_i[i])/float(Y_z[i]))
    i += 1 

plt.title('Ion and Impurity Density Profiles')
plt.plot(X_i,Y_i,':ro')
plt.plot(X_i,Y_z)

plt.show()
