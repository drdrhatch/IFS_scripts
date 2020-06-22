import numpy as np

def get_psi0(filename):
   f=open(filename,'r')
   rbs = f.read()
   f.close()
   rbs = rbs.split('\n')  
   psi0 = float(rbs[1].split()[2])
   print( "psi0",psi0)
   return psi0 

def calc_a(filename):
   f=open(filename,'r')
   rbs = f.read()
   f.close()
   rbs = rbs.split('\n')  
   a_factor = float(rbs[1].split()[3])
   rbs = np.genfromtxt(filename)
   isep = np.argmin(abs(rbs[:,0]-1.0))
   a = rbs[isep,22]
   return a_factor*a

