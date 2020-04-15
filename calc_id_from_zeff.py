import numpy as np
import matplotlib.pyplot as plt

Zeff = float(input("Enter Zeff:\n"))
Z = float(input("Enter Z:\n"))
ne = 1.0

#zeff = (ni+Z**2*nz)/ne
#nz = (zeff*ne - ni)/Z**2 
#ni = ne - Z*(zeff*ne - ni)/Z**2
#ni - ni/Z = ne - zeff*ne/Z

ni = (ne - Zeff*ne/Z)/(1-1.0/Z)
nz = (Zeff*ne - ni)/Z**2 

print("Zeff = ",Zeff)
print("Z = ", Z)
print("ni/ne", ni)
print("nz/ne", nz)

print("Quasineutrality (should be zero):",ni-ne+Z*nz)
print("Zeff",(ni+Z**2*nz)/ne)

