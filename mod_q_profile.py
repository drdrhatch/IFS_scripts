import numpy as np
import matplotlib.pyplot as plt
from read_write_geometry import *

geomfile = 'gene_0001_qmult0.958_hager_78697_nx0320_nz060'   #name of the geometry file output from GENE
plot_change = True         #change to True if one wants to plot(global) or see the info(local)
local=False     #Change to True if one wants to write local tracer efit file
qmult = 1.04384133611691   #multiply q profile by this factor
qoff = 0.0     #offset q profile by this amount

if local==True:
   with open(geomfile, 'r') as file:
      # read a list of lines into data
      data = file.readlines()
   q0=float(data[2].split()[2])
   q_mod=q0*qmult+qoff
   data[2]='q0 = '+str(q_mod)+'\n'

   with open(geomfile+'_qmult'+str(qmult)+'_qoff'+str(qoff), 'w') as file:
      file.writelines( data )
   if plot_change:
      print('q0='+str(q0))
      print('q_modified='+str(q_mod))
else: 
   parameters, geometry  = read_geometry_global(geomfile)

   q_original =  geometry['q']
   geometry['q'] = geometry['q']*qmult
   geometry['q'] = geometry['q'][:] + qoff
   
   parameters['q0'] = parameters['q0']*qmult
   parameters['q0'] = parameters['q0'] + qoff

   write_tracer_efit_file(parameters,geometry,geomfile+'_qmult'+str(qmult)+'_qoff'+str(qoff))

   if plot_change:
      #pfile = input('Enter profile file name:\n')
      #prof = np.genfromtxt(pfile)
      #plt.plot(prof[:,0], geometry['q'],label='new')
      #plt.plot(prof[:,0], q_original,label='original')
      #ax = plt.axis()
      #plt.axis((ax[0],ax[1],0.0,ax[3]))
      #plt.xlabel('rhot')

      plt.plot(geometry['q'],label='new')
      plt.plot(q_original,label='original')
      plt.legend(loc='lower right')
      plt.title('q')
      plt.show()

