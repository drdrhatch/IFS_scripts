#from parIOHelper import init_read_parameters
from read_write_geometry import *
import matplotlib.pyplot as plt
import numpy as np

def init_read_geometry(suffix, pars):

    geom_type = pars['magn_geometry'][1:-1]
    geom_file = geom_type + suffix
    geom_pars, geom_coeff = read_geometry_local(geom_file)

    return geom_type, geom_pars, geom_coeff

def local_grid_points(geom_coeff, show_plots = False):

    Z = geom_coeff['gl_z']
    R = geom_coeff['gl_R']

    if show_plots:
        plt.scatter(R,Z)
        plt.xlabel('R')
        plt.ylabel('Z')
        plt.axis('equal')
        plt.title('simulation grid points')
        plt.show()

    return R, Z

def ky(pars, geom_coeff, ktheta_cm = -1., show_plots = False):

    ggxx = geom_coeff['ggxx']
    ggxy = geom_coeff['ggxy']
    ggyy = geom_coeff['ggyy']

    gamma1 = ggxx * ggyy - ggxy ** 2

    kymin = pars['kymin']
    ky = np.sqrt(gamma1/ggxx)*kymin
    if ktheta_cm != -1.:
        ky /= ky[pars['nz0']/2]
        ky *= ktheta_cm * 100.


    R, Z = local_grid_points(geom_coeff, show_plots)

    if show_plots:
        plt.plot(ky,label='ky (1/m)')
        plt.legend()
        plt.show()

    return ky

def init_global_geometry(suffix, pars):

    geom_type = pars['magn_geometry'][1:-1]
    geom_file = geom_type + suffix
    geom_pars, geom_coeff = read_geometry_global(geom_file)

    return geom_type, geom_pars, geom_coeff

def q_Cy(geom_coeff):
    return geom_coeff['q'], geom_coeff['C_y']

def zGrid(geom_coeff, \
          pars, \
          center_only = False, \
          plot = True, \
          edge_opt = -1):

    nx = pars['nx0']
    nz = pars['nz0']
    if 'gBfield' in geom_coeff:
        gBfield = geom_coeff['gBfield']
    else:
        gBfield = geom_coeff['Bfield']
    if 'gjacobian' in geom_coeff:
        gjacobian = geom_coeff['gjacobian']
    else:
        gjacobian = geom_coeff['jacobian']

    if center_only:
        ikx_grid = [0]
    else:
        ikx_grid = np.arange(- nx / 2 + 1, nx / 2 + 1)

    zgrid_even_center = np.linspace(-1., 1., nz, \
                        endpoint = False)
    if not center_only:
        if nx % 2 == 1:
            zgrid_even = np.linspace(- nx, nx, nx * nz, \
                         endpoint = False)
        else :
            zgrid_even = np.linspace(- (nx - 1), (nx + 1), \
                         nx * nz, endpoint = False)

    if 'edge_opt' in pars:
      if edge_opt == -1:
        edge_opt = pars['edge_opt']
      else:
        edge_opt = float(edge_opt)
    else:
      edge_opt = 0.

    if edge_opt != 0:
        zgrid_edge = np.zeros(nx * nz, dtype = 'float128')
        
        N = np.arcsinh(edge_opt * zgrid_even_center[0] * \
            np.pi) / zgrid_even_center[0] / np.pi

        zgrid_edge_center = 1. / edge_opt * \
                            np.sinh(N * zgrid_even_center * \
                            np.pi) / np.pi

        dz = np.zeros(nz, dtype = 'float128')
        for i in np.arange(nz / 2 + 1, nz):
            dz[i] = zgrid_edge_center[i] - \
                    zgrid_edge_center[i - 1]
        for i in np.arange(nz / 2 - 1, - 1, - 1):
            dz[i] = zgrid_edge_center[i + 1] - \
                    zgrid_edge_center[i]

        for i in ikx_grid:
            this_zgrid_edge = i * 2. + zgrid_edge_center
            zgrid_edge[(i - ikx_grid[0]) * \
            nz : (i - ikx_grid[0] + 1) * nz] \
            = this_zgrid_edge
        if center_only:
            zgrid = zgrid_edge_center
        else:
            zgrid = zgrid_edge
    else:
        if center_only:
            zgrid = zgrid_even_center
        else:
            zgrid = zgrid_even

        dz = np.ones(nz, dtype = 'float128') * 2. / nz

    jacobian_center = 1. / np.pi / gjacobian / gBfield
    
    if not center_only:
        jacobian = np.zeros(nx * nz, dtype = 'float128')
        for i in ikx_grid:
            jacobian[(i-ikx_grid[0])*nz:(i-ikx_grid[0]+1)*nz]=\
                jacobian_center
    if center_only:
        jacobian = jacobian_center

    if plot:
        plt.plot(zgrid, label = 'zgrid')
        plt.legend(loc=2)
        plt.show()
        plt.plot(jacobian, label = 'jacobian')
        plt.legend(loc=2)
        plt.show()

    return zgrid, jacobian
