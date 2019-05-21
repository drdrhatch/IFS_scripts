import numpy as np
import matplotlib.pyplot as plt
from ParIO import *

pfile = 'parameters_1'
efile = 'energy3d_1'
#pfile = 'parameters.dat'
#efile = 'energy3d.dat'
time = 10.0

par = Parameters()
par.Read_Pars(pfile)
pars = par.pardict
ntot = pars['nx0']*pars['nky0']*pars['nz0']
n_entries = 6
mem_single = ntot*8
mem_tot = n_entries*ntot*8  #6 entries, real data (8)
nint = 4
nreal = 8

f = open(efile,'rb')
f.seek(4)
time0 = np.fromfile(f,dtype='float64',count=1)
#print "Time step 1"
#print time0
#print "Time step 2"
#for i in range(200):
#   f.seek(4+8+mem_tot-100+i)
#   time0 = np.fromfile(f,dtype='float64',count=1)
#   if abs(41-time0) < 5.0:
#      print i, time0
#print "Time step 3"
#for i in range(400):
#   f.seek(4+2*8+2*mem_tot-200+i)
#   time0 = np.fromfile(f,dtype='float64',count=1)
#   if abs(41-time0) < 10.0:
#      print i, time0

#for i in range(30):
#   f.seek(4+i*8+i*mem_tot+i*56)
#   time0 = np.fromfile(f,dtype='float64',count=1)
#   print time0
#   f.seek(4+i*8+i*mem_tot+i*56+8+8)
#   test = np.fromfile(f,dtype='float64',count=1)
#   print test
#   num = 3
#   f.seek(4+i*8+i*mem_tot+i*56+8+num*8+num*mem_single+8)
#   test = np.fromfile(f,dtype='float64',count=1)
#   print test

zgrid = np.empty(pars['nz0'])
xgrid = np.linspace(pars['x0'] - pars['lx_a']/2.0 , pars['x0'] + pars['lx_a']/2.0, pars['nx0'])
zgrid = np.linspace(-np.pi,np.pi, pars['nz0'],endpoint = False)
print 'xgrid',xgrid

continue_read = True
i=0
while continue_read:
   f.seek(4+i*8+i*mem_tot+i*56)
   tin = np.fromfile(f,dtype='float64',count=1)
   if tin == 0.0 or tin:
      time = np.append(time,tin)
      i += 1
   else:
      continue_read = False

i -= 1
print "i",i
print "time",time

f.seek(4+i*8+i*mem_tot+i*56+16)
ein = np.fromfile(f,dtype='float64',count=ntot)
ein = np.reshape(ein,(pars['nky0'],pars['nx0'],pars['nz0']),order='F')
print "np.shape(ein)",np.shape(ein)
plt.contourf(xgrid,zgrid,np.transpose(ein[0,:,:]),50)
plt.colorbar()
plt.show()
plt.contourf(xgrid,zgrid,np.transpose(ein[1,:,:]),50)
plt.colorbar()
plt.show()

for i in range(5):
   print pars['nx0']/(i+2)
   plt.plot(zgrid,ein[1,pars['nx0']/(i+2),:])
   #plt.plot(zgrid,ein[1,190,:])
plt.show()

ein_kysum = np.sum(ein,axis=0)
ein_kyzsum = np.sum(ein_kysum,axis=1)
ein_kyxsum = np.sum(ein_kysum,axis=0)

plt.contourf(xgrid,zgrid,np.transpose(ein_kysum[:,:]),50)
plt.colorbar()
plt.show()

plt.plot(xgrid,ein_kyzsum)
plt.show()

plt.plot(zgrid,ein_kyxsum)
plt.show()

plt.contourf(xgrid,zgrid,np.transpose(np.log(ein_kysum[:,:])),50)
plt.colorbar()
plt.show()

plt.semilogy(xgrid,ein_kyzsum)
plt.show()

plt.semilogy(zgrid,ein_kyxsum)
plt.show()



#for i in range(5):
#   print "###############"
#   f.seek(4+i*8+i*mem_tot+i*56)
#   time0 = np.fromfile(f,dtype='float64',count=1)
#   print time0
#   f.seek(4+i*8+i*mem_tot+i*56 + 16 )
#   time0 = np.fromfile(f,dtype='float64',count=1)
#   print time0
#   f.seek(4+i*8+i*mem_tot+i*56 + 16 + mem_single + 8 )
#   time0 = np.fromfile(f,dtype='float64',count=1)
#   print time0
#   f.seek(4+i*8+i*mem_tot+i*56 + 16 + mem_single + 8 +mem_single + 8 )
#   time0 = np.fromfile(f,dtype='float64',count=1)
#   print time0
#   f.seek(4+i*8+i*mem_tot+i*56 + 16 + mem_single + 8 +mem_single + 8  + mem_single + 8)
#   time0 = np.fromfile(f,dtype='float64',count=1)
#   print time0
#   f.seek(4+i*8+i*mem_tot+i*56 + 16 + mem_single + 8 +mem_single + 8  + mem_single + 8 + mem_single + 8)
#   time0 = np.fromfile(f,dtype='float64',count=1)
#   print time0
#   f.seek(4+i*8+i*mem_tot+i*56 + 16 + mem_single + 8 +mem_single + 8  + mem_single + 8 + mem_single + 8 + mem_single + 8)
#   time0 = np.fromfile(f,dtype='float64',count=1)
#   print time0
#   f.seek(4+i*8+i*mem_tot+i*56 + 16 + mem_single + 8 +mem_single + 8  + mem_single + 8 + mem_single + 8 + mem_single + 8 + mem_single + 8  )
#   time0 = np.fromfile(f,dtype='float64',count=1)
#   print time0
   

