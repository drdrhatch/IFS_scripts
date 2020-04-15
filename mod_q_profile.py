import numpy as np
import matplotlib.pyplot as plt
from read_write_geometry import *

geomfile = 'tracer_efit_0001_J78697_nx0128_nz080'
plot_change = False

qmult = 0.958  #multiply q profile by this factor
qoff = 0.0  #offset q profile by this amount

parameters, geometry  = read_geometry_global(geomfile)
q_original = geometry['q']*1.0

geometry['q'] = geometry['q']*qmult
geometry['q'] = geometry['q'][:] + qoff

parameters['q0'] = parameters['q0']*qmult
parameters['q0'] = parameters['q0'] + qoff

write_tracer_efit_file(parameters,geometry,geomfile+'_qmult'+str(qmult)+'_qoff'+str(qoff))

if plot_change:
   pfile = input('Enter profile file name:\n')
   prof = np.genfromtxt(pfile)
   plt.plot(prof[:,0],geometry['q'],label='new')
   plt.plot(prof[:,0],q_original,label='original')
   ax = plt.axis()
   plt.axis((ax[0],ax[1],0.0,ax[3]))
   plt.legend(loc='lower right')
   plt.xlabel('rhot')
   plt.title('q')
   plt.show()

