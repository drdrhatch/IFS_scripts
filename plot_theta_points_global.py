import numpy as np
from read_write_geometry import *

geomfile = 'tracer_efit_0001'

parameters, geom = read_geometry_global(geomfile)

R = geom['geo_R']
Z = geom['geo_Z']
nx = np.size(R)/parameters['gridpoints'] 

for i in range(int(parameters['gridpoints'])):
    plt.scatter(R[i,nx/2],Z[i,nx/2])
plt.xlabel('R(m)')
plt.ylabel('Z(m)')
plt.show()


