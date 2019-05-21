import numpy as np
import sys

def first_derivative(f_in,x_in):
    x = x_in.flatten()
    f = f_in.flatten()
    dx = x[1]-x[0]
    dx1 = x[2] - x[1]
    if abs(dx - dx1) > 1.E-10:
        print(dx, dx1)
        sys.exit("x grid must be uniform")
    dfdx = np.empty(len(f))
    for i in range(len(f)):
        if i == 0:
            dfdx[i] = (f[i+1] - f[i])/dx
        elif i == 1 or i == len(f)-2:
            dfdx[i] = (f[i+1]-f[i-1])/2./dx
        elif i == len(f)-1:
            dfdx[i] = (f[i] - f[i-1])/dx
        else:
            dfdx[i] = (-f[i+2]+8.*f[i+1]-8.*f[i-1]+f[i-2])/12./dx
	
    return dfdx

