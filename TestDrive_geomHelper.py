import matplotlib.pyplot as plt
import numpy as np
from geomHelper import *
from fieldHelper import *
from parIOHelper import *
import sys

suffix = sys.argv[1]

if not suffix =='.dat':
   suffix = '_'+suffix

pars = init_read_parameters(suffix)

geom_type, geom_pars, geom_coeff = init_read_geometry(suffix, pars)

zgrid, jacobian = zGrid(geom_coeff, pars, False, False)

phi, apar = local_eigenfunctions(pars, suffix, False, False)

if 'kx_center' in pars:
    kx_center = pars['kx_center']
else:
    kx_center = 0.
theta_0 = kx_center/np.pi/pars['shat']/pars['kymin']
theta_0 = np.around(theta_0,decimals=1)
show_plots = True
if show_plots:
    plt.plot(zgrid,abs(phi),label='abs phi',color='black')
    plt.plot(zgrid,np.real(phi),label='real phi',color='blue')
    plt.plot(zgrid,np.imag(phi),label='imag phi',color='red')
    plt.title('ky = '+str(pars['kymin'])+', kx_center = '+str(theta_0)+r'$\pi$')
    plt.grid()
    plt.legend()
    plt.show()
if show_plots:
    plt.plot(zgrid,abs(apar),label='abs apar',color='black')
    plt.plot(zgrid,np.real(apar),label='real apar',color='blue')
    plt.plot(zgrid,np.imag(apar),label='imag apar',color='red')
    plt.title('ky = '+str(pars['kymin'])+', kx_center = '+str(theta_0)+r'$\pi$')
    plt.grid()
    plt.legend()
    plt.show()

epar(pars, geom_coeff, suffix)
