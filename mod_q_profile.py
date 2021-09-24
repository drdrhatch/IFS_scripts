import numpy as np
import matplotlib.pyplot as plt
from read_write_geometry import *

geomfile = 'tracer_efit'   #name of the geometry file output from GENE
plot_change = True         #change to True if one wants to plot(global) or see the info(local)
local=True     #Change to True if one wants to write local tracer efit file
qmult = 1.05   #multiply q profile by this factor
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

   geometry['q0'] = geometry['q0']*qmult
   geometry['q0'] = geometry['q0'][:] + qoff
   
   parameters['q0'] = parameters['q0']*qmult
   parameters['q0'] = parameters['q0'] + qoff

   write_tracer_efit_file(parameters,geometry,geomfile+'_qmult'+str(qmult)+'_qoff'+str(qoff))

   if plot_change:
      #pfile = input('Enter profile file name:\n')
      #prof = np.genfromtxt(pfile)
      plt.plot(geometry['q'],label='new')
      plt.plot(q_original,label='original')
      #ax = plt.axis()
      #plt.axis((ax[0],ax[1],0.0,ax[3]))
      plt.legend(loc='lower right')
      plt.xlabel('rhot')
      plt.title('q')
      plt.show()

