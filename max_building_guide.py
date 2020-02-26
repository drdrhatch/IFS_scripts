import optparse as op
import math
import cmath
import sys
import numpy as np
from fieldlib import *
from ParIO import *
from read_write_geometry import *
#from max_shat import *

def initi_with_input(): #this function initalized the parameter and time

#this requires field file present.

    parser=op.OptionParser(description='Plots mode structures and calculates various interesting quantities.')
    parser.add_option('--plot_theta','-g',action='store_const',const=False,help = 'Plot all plots.',default=True)
    parser.add_option('--plot_ballooning','-b',action='store_const',const=False,help = 'Plot all plots.',default=True)
    parser.add_option('--plot_all','-a',action='store_const',const=1,help = 'Plot all plots.',default=False)
    parser.add_option('--time','-t',type = 'float',action='store',dest="time0",help = 'Time to plot mode structure.',default=-1)
    parser.add_option('--idb','-i',type = 'str',action='store',dest="idb_file",help = 'ITERDB file name.',default='empty')
    parser.add_option('--pfile','-f',type = 'str',action='store',dest="prof_file",help = 'ITERDB file name.',default='empty')
    options,args=parser.parse_args()
    print(("options",options))
    print(("args",args))
    if len(args)!=1:
        exit("""
    Please include run number as argument (e.g., 0001)."
        \n""")
    suffix = args[0]
    plot_all=options.plot_all
    print(("options.plot_ballooning", options.plot_ballooning))
    plot_ballooning=options.plot_ballooning
    plot_theta=options.plot_theta
    idb_file = options.idb_file
    prof_file = options.prof_file
    time0=float(options.time0)
    
    suffix = '_'+suffix
    
    par = Parameters()
    par.Read_Pars('parameters'+suffix)
    pars = par.pardict

    #field = fieldfile('field'+suffix,pars,False)
    field = fieldfile('field'+suffix,pars)
    #print "time0",time0
    #print "field.tfld",field.tfld
    time = np.array(field.tfld)
    if time0 == -1:
        itime = -1
        itime0 = len(time)-1
    else: 
        itime = np.argmin(abs(time - time0))
        itime0 = itime

    print(("Looking at the mode structure at time:",time[itime]))
    #field.set_time(time[itime],itime0)
    field.set_time(time[itime])
    
    

    return time, suffix, pars, field

def phi_and_apar(suffix,field,pars): #The for loop is adding the phase correction into the field output

    if 'x_local' in pars:
        if pars['x_local']:
            x_local = True
        else:
            x_local = False 
    else:
        x_local = True

    dz = float(2.0)/float(field.nz)
    ntot = field.nz*field.nx
    
    #print "ntot",field.nz*field.nx

    #Adding phase correction to the field
    if 'n0_global' in pars:
        if 'q0' in pars:
            q0=float(pars['q0'])
        else:
            q0=0.5
        phase_fac = -np.e**(-2.0*np.pi*(0.0+1.0J)*pars['n0_global']*q0)
    else:
        phase_fac = -1.0
        
    #print "phase_fac",phase_fac
    if 'shat' not in pars:
        print('shat is not in pars')
        if x_local:
            parameters, geometry = read_geometry_local(pars['magn_geometry'][1:-1]+suffix)
        else:
            parameters, geometry = read_geometry_global(pars['magn_geometry'][1:-1]+suffix)
        shat = parameters['shat']
        #shat = local_shat(suffix, pars['x0'])
    else:
        shat=pars['shat']
    #print('shat={shat}')
    
    if x_local:
        phi = np.zeros(ntot,dtype='complex128')
        apar = np.zeros(ntot,dtype='complex128')

        if float(shat) < 0.0:
            for i in range(int(field.nx/2)+1):
                phi[(i+int(field.nx/2))*field.nz:(i+int(field.nx/2)+1)*field.nz]=field.phi()[:,0,-i]*phase_fac**i
                if i < int(field.nx/2):
                    phi[(int(field.nx/2)-i-1)*field.nz : (int(field.nx/2)-i)*field.nz ]=field.phi()[:,0,i+1]*phase_fac**(-(i+1))
                if pars['n_fields']>1:
                    apar[(i+int(field.nx/2))*field.nz:(i+int(field.nx/2)+1)*field.nz]=field.apar()[:,0,-i]*phase_fac**i
                    if i < int(field.nx/2):
                        apar[(int(field.nx/2)-i-1)*field.nz : (int(field.nx/2)-i)*field.nz ]=field.apar()[:,0,i+1]*phase_fac**(-(i+1))
        else:
            for i in range(int(field.nx/2)):
                #print phi[(i+int(field.nx/2))*field.nz:(i+int(field.nx/2)+1)*field.nz]
                #print (i+int(field.nx/2))*field.nz
                #print (i+int(field.nx/2)+1)*field.nz
                phi[(i+int(field.nx/2))*field.nz:(i+int(field.nx/2)+1)*field.nz]=field.phi()[:,0,i]*phase_fac**i
                if i < int(field.nx/2):
                    phi[(int(field.nx/2)-i-1)*field.nz : (int(field.nx/2)-i)*field.nz ]=field.phi()[:,0,-1-i]*phase_fac**(-(i+1))
                if pars['n_fields']>1:
                    apar[(i+int(field.nx/2))*field.nz:(i+int(field.nx/2)+1)*field.nz]=field.apar()[:,0,i]*phase_fac**i
                    if i < int(field.nx/2):
                        apar[(int(field.nx/2)-i-1)*field.nz : (int(field.nx/2)-i)*field.nz ]=field.apar()[:,0,-1-i]*phase_fac**(-(i+1))
    
        #phi = phi/field.phi()[int(field.nz/2),0,0]
        #apar = apar/field.phi()[int(field.nz/2),0,0]
            zgrid = np.arange(ntot)/float(ntot-1)*(2*field.nx-dz)-field.nx
            return phi,apar,zgrid
            #********Return of location simultion


    else:
        apar = field.apar()[:,0,:]
        phi  = field.phi() [:,0,:]

        #gpars,geometry = read_geometry_global(pars['magn_geometry'][1:-1]+suffix)
        #phi_bnd = np.zeros((field.nz+4,field.ny,field.nx),dtype = 'complex128')
        #print(geometry['q'])
        #phase = (0.0+1.0J)*pars['n0_global']*2.0*np.pi*geometry['q']
        #for i in range(field.nx):
        #    phi_bnd[-2,:,i] = phi_bnd[2,:,i]*np.e**(-phase[i])
        #    phi_bnd[-1,:,i] = phi_bnd[3,:,i]*np.e**(-phase[i])
        #    phi_bnd[0,:,i]  = phi_bnd[-4,:,i]*np.e**(phase[i])
        #    phi_bnd[1,:,i]  = phi_bnd[-3,:,i]*np.e**(phase[i])

        #print(np.shape(phi_bnd))

        dz = 2.0/field.nz
        zgrid = np.arange(field.nz)/float(field.nz-1)*(2.0-dz)-1.0

        if 'lx_a' in pars:
            xgrid = np.arange(field.nx)/float(field.nx-1)*pars['lx_a']+pars['x0']-pars['lx_a']/2.0
        else:
            xgrid = np.arange(field.nx)/float(field.nx-1)*pars['lx'] - pars['lx']/2.0
        return phi, apar, xgrid, zgrid


