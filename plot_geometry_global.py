import numpy as np
from read_write_geometry import *
from mpl_toolkits import mplot3d

geomfile = 'tracer_efit_1'
kymin = 0.1511
ny0 = 10

pars, geom = read_geometry_global(geomfile)

#plt.scatter(geom['gl_z'],geom['gl_phi'])
#plt.xlabel('z')
#plt.ylabel('phi')
#plt.show()
#
#plt.plot(geom['gl_z'],geom['gl_phi'])
#plt.xlabel('z')
#plt.ylabel('phi')
#plt.show()

dz = float(2.0)/float(pars['gridpoints'])
zgrid = np.arange(pars['gridpoints'])/float(pars['gridpoints']-1)*(2.0-dz)-1.0
zgrid *= np.pi
print "zgrid",zgrid

Ly = 2.0*np.pi/kymin
dy = Ly/ny0
ygrid = np.arange(-Ly/2,Ly/2,dy)
print "ygrid",ygrid
print "Ly",Ly

print np.shape(geom['geo_R'])
print np.shape(geom['geo_Z'])
print np.shape(geom['C_y'])
print np.shape(geom['q'])

nx0 = len(geom['geo_R'][0,:])
nz0 = pars['gridpoints']
R0 = np.empty(0)
Z0 = np.empty(0)
phi0 = np.empty(0)

for i in range(nx0):
    for j in range(ny0):
        for k in range(nz0):
            R0 = np.append(R0,geom['geo_R'][k,i])
            Z0 = np.append(Z0,geom['geo_Z'][k,i])
            #note: phi = -y/Cy +q theta
            phi0 = np.append(phi0,-ygrid[j]/geom['C_y'][i] + geom['q'][i]*zgrid[k])

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(R0*np.sin(phi0),R0*np.cos(phi0),Z0,marker='.')
plt.xlabel('X(m)')
plt.ylabel('Y(m)')
plt.show()


#fig = plt.figure()
#ax = plt.axes(projection='3d')
#ax.plot3D(geom['gl_R']*np.sin(geom['gl_phi']),geom['gl_R']*np.cos(geom['gl_phi']),geom['gl_z'])
#plt.xlabel('X(m)')
#plt.ylabel('Y(m)')
#plt.show()


