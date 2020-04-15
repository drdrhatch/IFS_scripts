#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from ParIO import *
import optparse as op
import sys
import os

parser=op.OptionParser(description='Eigenvector projection routine.  Required arguments: ky index and kx_center index for df file.')
parser.add_option('--start_time','-s',type='float',dest='start_time',default=0.0,help = 'Start time [default: 0.0]')
parser.add_option('--ev_dir','-e',type='str',dest='ev_dir',default='./',help = 'Eigenvector directory [default: \'./\']')
parser.add_option('--df_dir','-d',type='str',dest='df_dir',default='./',help = 'Directory for df files [default: \'./\']')
parser.add_option('--nev','-n',type='int',dest='nev',default=20,help = 'Number of eigenvalues for projection [default: 20]')
parser.add_option('--test_ortho','-o',action='store_const',const=1,help = 'Verify orthogonality of selected eigenvectors')
options,args=parser.parse_args()
start_time=options.start_time
ev_dir = options.ev_dir
df_dir = options.df_dir
nev = options.nev
test_ortho = options.test_ortho

#use_r_a=options.rovera
if len(args)!=2:
    exit("""
Please include ky_index and kx_index as arguments."
    \n""")

ky_index = int(float(args[0]))
kx_index = int(float(args[1]))

ky_index_string = '%04d' % ky_index
kx_index_string = '%04d' % kx_index
df_file_name = df_dir + 'df_ky'+ky_index_string+'kx'+kx_index_string+'.dat'
print("df_file_name",df_file_name)

####### Start functions ########
####### Start functions ########
####### Start functions ########

def read_evec(file_name,which_ev,par,swap_endian=False):
   """Reads an eigenvector from an eigenvector: file.  Eigenvector determined by \'which_ev\'"""
   f = open(file_name,'r')
   ntot = par['nx0']*par['nz0']*par['nv0']*par['nw0']*par['n_spec']
   mem_tot = ntot*16
   #gt0 = np.empty((nkx_keep,par['nz0'],par['nv0'],par['nw0'],par['n_spec']))
   gt0 = np.empty(ntot)
   #f.seek(4+which_ev*(20+mem_tot))
   f.seek(4+which_ev*(20+mem_tot))
   evnum = np.fromfile(f,dtype='int32',count=1)
   #f.seek(4+which_ev*(20+mem_tot)+8+8)
   f.seek(4+which_ev*(20+mem_tot)+12)
   gt0 = np.fromfile(f,dtype='complex128',count=ntot) 
   gt0 = np.reshape(gt0,(par['nx0'],par['nz0'],par['nv0'],par['nw0'],par['n_spec']),order='F')
   return evnum,gt0

def scalar_product(g1,g2):
   return np.sum(np.conj(g1)*g2)

def get_sorted_evecs(ev_dir,par,nev):
   evals = np.genfromtxt(ev_dir+'/eigenvalues.dat') 
   gammas = np.empty(0)
   omegas = np.empty(0)
   indices = np.empty(0,dtype='int')
   evs = evals*1.0
   levecs = np.empty((par['nx0'],par['nz0'],par['nv0'],par['nw0'],par['n_spec'],nev),dtype='complex128')
   revecs = np.empty((par['nx0'],par['nz0'],par['nv0'],par['nw0'],par['n_spec'],nev),dtype='complex128')
   for i in range(nev):
      indices=np.append(indices,np.argmax(evs[:,0]))
      gammas=np.append(gammas,evs[indices[i],0])
      omegas=np.append(omegas,evs[indices[i],1])
      #print gammas[i],omegas[i],indices[i]
      evs[indices[i],0] = -1.0e20
      #Get left eigenvectors  
      evnum, levecs[:,:,:,:,:,i] = read_evec(ev_dir+'/eigenvectors_l.dat',indices[i],par)
      evnum, revecs[:,:,:,:,:,i] = read_evec(ev_dir+'/eigenvectors_r.dat',indices[i],par)
   return gammas,omegas,indices,levecs,revecs 

def get_time_from_dffile(file_name,ntot):
        
    print("Reading file",file_name)
    file_exists=os.path.isfile(file_name)
    if file_exists:
        pass
    else:
        print("File does not exist:",file_name)
        sys.exit()

    f=open(file_name,'r')
    
    mem_tot=ntot*16
    time=np.empty(0)
    continue_read=1
    i=0
    while (continue_read): 
      f.seek(4+i*(mem_tot+24))
      i=i+1
      input=np.fromfile(f,dtype='float64',count=1)
      if input==0 or input:
          time = np.append(time,input)
          #print input
      else:
          continue_read=0

    f.close()
    return time

def read_dfout_file(file_name,num_ky):
   f = open(file_name,'r')
   data = f.read()
   lines = data.split('\n')
   ky_indices = np.empty(num_ky,dtype = 'int')
   kx_center = np.empty(num_ky, dtype = 'int')
   kx_shift = np.empty(num_ky, dtype = 'int')
   nkx_keep = np.empty(num_ky, dtype = 'int')
   kx_modes = {} 
  
   nline = 2
   for i in range(num_ky):
      ky_indices[i] = int(float(lines[nline].split()[-1]))
      nline += 1
      kx_center[i] = int(float(lines[nline].split()[-1]))
      nline += 1
      kx_shift[i] = int(float(lines[nline].split()[-1]))
      nline += 1
      nkx_keep[i] = int(float(lines[nline].split()[-1]))
      nline += 2
      kx_m = np.empty(nkx_keep[i],dtype = 'int')
      for j in range(nkx_keep[i]):
         kx_m[j] = int(float(lines[nline]))
         nline += 1
      kx_modes[i] = kx_m      
      nline += 2

   #print ky_indices
   #print kx_center
   #print kx_shift
   #print nkx_keep
   #print kx_modes
   return ky_indices, kx_center, kx_shift, nkx_keep, kx_modes

def read_time_step_df(file_name,which_itime,nkx_keep,par,swap_endian=False):
   """Reads a time step from dfout file.  Time step determined by \'which_itime\'"""
   f = open(file_name,'r')
   ntot = nkx_keep*par['nz0']*par['nv0']*par['nw0']*par['n_spec']
   mem_tot = ntot*16
   #gt0 = np.empty((nkx_keep,par['nz0'],par['nv0'],par['nw0'],par['n_spec']))
   gt0 = np.empty(ntot)
   f.seek(4+which_itime*(24+mem_tot))
   time = np.fromfile(f,dtype='float64',count=1)
   f.seek(4+which_itime*(24+mem_tot)+8+8)
   gt0 = np.fromfile(f,dtype='complex128',count=ntot) 
   gt0 = np.reshape(gt0,(nkx_keep,par['nz0'],par['nv0'],par['nw0'],par['n_spec']),order='F')
   return time,gt0


def crop_df(dft,dfnkx,evnkx):
   df0 = np.empty( (evnkx,np.shape(dft)[1], \
                          np.shape(dft)[2], \
                          np.shape(dft)[3], \
                          np.shape(dft)[4]) \
                 ,dtype='complex128')
   df0[0,:,:,:,:] = dft[0,:,:,:,:]
   for i in range((evnkx-1)/2):
      #print "crop_dft i",i
      #print "crop_dft evnkx-(i+1)",evnkx-(i+1)
      #print "crop_dft dfnkx-(i+1)",dfnkx-(i+1)
      df0[i+1,:,:,:,:] = dft[i+1,:,:,:,:]
      df0[evnkx-(i+1),:,:,:,:] = dft[dfnkx-(i+1),:,:,:,:]
   return df0

####### Start script ########
####### Start script ########
####### Start script ########



####### Set up df stuff ########
####### Set up df stuff ########
####### Set up df stuff ########

par = Parameters()

par.Read_Pars(df_dir+'/parameters.dat')

dfpars = par.pardict

ky_indices, kx_center, kx_shift, nkx_keep, kx_modes = read_dfout_file(df_dir+'/dfout.info',dfpars['num_ky_modes'])

df_index = -999
for i in range(dfpars['num_ky_modes']):
   if kx_index == kx_center[i] and ky_index == ky_indices[i]:
      df_index = i
if df_index == -999:
   print("ERROR! invalid kx_index or ky_index.")
   sys.exit()

nkx = nkx_keep[df_index]
nv0 = dfpars['nv0']
nw0 = dfpars['nw0']
nz0 = dfpars['nz0']
n_spec = dfpars['n_spec']

time = get_time_from_dffile(df_file_name,nv0*nw0*nz0*nkx*n_spec)
#print time

start_index = np.argmin(abs(time - start_time))
print("Starting at t = ",start_time)
print("start_index",start_index)


####### Set up ev stuff ########
####### Set up ev stuff ########
####### Set up ev stuff ########

parev = Parameters()
parev.Read_Pars(ev_dir+'/parameters.dat')
evpars = parev.pardict

print("Checking consistency of df parameters and ev parameters.")
if evpars['nz0'] != dfpars['nz0']:
   print("Error! Mismatch between dfout and ev nz0!")
   sys.exit()
if evpars['nv0'] != dfpars['nv0']:
   print("Error! Mismatch between dfout and ev nv0!")
   sys.exit()
if evpars['nw0'] != dfpars['nw0']:
   print("Error! Mismatch between dfout and ev nw0!")
   sys.exit()
if evpars['n_spec'] != dfpars['n_spec']:
   print("Error! Mismatch between dfout and ev n_spec!")
   sys.exit()
if evpars['nx0'] % 2 != 1:
   print("Error! This routine can only handle nx0=odd.  Fix if necessary.")
   sys.exit()
if evpars['nx0'] > nkx:
   print("Error! This routine expects nx0 for eigenvectors to be <= nx0 for dfout data.  Can be fixed if necessary.")
   sys.exit()

print("Sorting eigenvalues and eigenvectors.")
gammas, omegas, indices, levecs, revecs = get_sorted_evecs(ev_dir,evpars,nev)
plt.scatter(gammas,omegas)
plt.xlabel(r'$\gamma (v_{ref}/L_{ref})$',size=18)
plt.ylabel(r'$\omega (v_{ref}/L_{ref})$',size=18)
plt.show()

if test_ortho:
   print("Testing orthogonality of eigenvectors.")
   sps = np.empty((nev,nev))
   for i in range(nev):
      for j in range(nev):
         sps[i,j] = np.abs(scalar_product(revecs[:,:,:,:,:,i],levecs[:,:,:,:,:,j]))
         if i==j and sps[i,j] < 0.999999:
            print("Orthogonality failure at ",i,j)
            print("Scalar product:",sps[i,j])
            sys.exit()
         if i!=j and sps[i,j] > 1.0e-13:
            print("Orthogonality Failure at ",i,j)
            print("Scalar product:",sps[i,j])
            sys.exit()
   plt.contourf(sps,50)
   plt.colorbar()
   plt.show()      
    
####### Do projection ########
####### Do projection ########
####### Do projection ########

ct = np.empty((len(time)-start_index,nev),dtype='complex128')
print("Doing eigenvector projection.")

for i in range(start_index,len(time)):
   t0,dft = read_time_step_df(df_file_name,i,nkx,dfpars)   
   if nkx > evpars['nx0']:
      df0 = crop_df(dft,nkx,evpars['nx0'])
   for j in range(nev):
      ct[i-start_index,j] = scalar_product(revecs[:,:,:,:,:,j],df0) 
  
for i in range(nev):
   plt.plot(time[start_index:],abs(ct[:,i]),label=str(i+1))
plt.legend()
plt.show()   

for i in range(nev):
   plt.semilogy(time[start_index:],abs(ct[:,i]),label=str(i+1))
plt.legend()
plt.show()   

