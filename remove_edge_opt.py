import numpy as np
from interp import *

def remove_edge_opt_complex(array,edge_opt):
    nz = len(array)
    dz = 2.0*np.pi/nz
    zgrid_even = np.arange(nz)/float(nz-1)*(2.0*np.pi-dz)-np.pi
    N = np.arcsinh(edge_opt*zgrid_even[-1])/zgrid_even[-1]
    zprime_even = zgrid_even
    z_of_zprime_even = 1.0/edge_opt*np.sinh(N*zprime_even) 
    Rarray = interp(z_of_zprime_even,np.real(array),zgrid_even)
    Iarray = interp(z_of_zprime_even,np.imag(array),zgrid_even)
    array_out = Rarray + 1.0J*Iarray
    return  array_out

def remove_edge_opt(array,edge_opt):
    nz = len(array)
    dz = 2.0*np.pi/nz
    zgrid_even = np.arange(nz)/float(nz-1)*(2.0*np.pi-dz)-np.pi
    N = np.arcsinh(edge_opt*zgrid_even[-1])/zgrid_even[-1]
    #zprime = 1/N*np.arcsinh(edge_opt*zgrid_even)
    zprime_even = zgrid_even
    z_of_zprime_even = 1.0/edge_opt*np.sinh(N*zprime_even) 
    array_out = interp(z_of_zprime_even,array,zgrid_even)
    return  array_out

