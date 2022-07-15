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

def get_xgrid(geomfile_without_edge_opt,geomfile_with_edge_opt,lx_a,x0,plot=False):
    from read_write_geometry import read_geometry_global
    parameters, geometry  = read_geometry_global(geomfile_without_edge_opt)
    nx_norm=256
    q_norm = geometry['q']
    xgrid_norm = np.arange(nx_norm)/float(nx_norm-1)*lx_a+x0-lx_a/2.0

    
    parameters, geometry  = read_geometry_global(geomfile_with_edge_opt)
    nx_opt=256
    q_opt = geometry['q']
    
    xgrid_opt=np.interp(q_opt,q_norm,xgrid_norm)    #Interprolate f with new_x 
    
    if plot==True:
        plt.clf()
        plt.plot(xgrid_norm,q_norm,label='q0(norm)')
        plt.plot(xgrid_opt,q_opt,label='q0(edge opt)')
        plt.legend()
        plt.show()

        plt.clf()
        plt.plot(xgrid_opt)
        plt.show()

    return xgrid_opt


