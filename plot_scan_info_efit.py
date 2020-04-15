#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ParIO import *
#from matplotlib.mlab import griddata
from scipy.interpolate import griddata
import optparse as op
from interp import *
#"""
#Plots data from scan_info.dat (created by scan_info.py).
#Assumes a scan in ky and one of the following: kx_center or x0.
#The routine can handle a ky-dependent kx_center scan.  
#An x0 scan is assumed to be ky-independent.
#"""

parser=op.OptionParser(description='Plots data from scan_info.dat.')
parser.add_option('--get_global','-g',action='store_true',help='Filter out local modes.')
options,args=parser.parse_args()
if len(args)!=1:
    exit("""
Please include scan number as argument (e.g., 0001)."
    \n""")
suffix = args[0]
get_global = options.get_global
os.chdir('scanfiles'+suffix)

data = np.genfromtxt('scan_info.dat')

par = Parameters()
par.Read_Pars('parameters')
pars = par.pardict

if type(pars['scan_dims']) == str:
    scan_dims = pars['scan_dims'].split()
    numscan_tot = 1
    for i in scan_dims:
        numscan_tot *= int(i)
else:
    numscan_tot = pars['scan_dims']

def filter_numerical_modes(data_in,pars):
    ntot = len(data_in[:,0])
    dz = 2.0/float(pars['nz0'])
    for i in range(ntot):
        if data_in[i,7] < 2.0*dz:
            data_in[i,4] = np.nan
            data_in[i,5] = np.nan
            data_in[i,10] = np.nan
        elif data_in[i,4] < -1.0:
            data_in[i,4] = np.nan 
            data_in[i,5] = np.nan
            data_in[i,10] = np.nan
    return data_in

def filter_local_modes(data_in,pars):
    #Modes with global_factor < ? cannot fit in the scale range of the drive (rough estimate)
    ntot = len(data_in[:,0])
    for i in range(ntot):
        if data_in[i,13] < 2.0:
            data_in[i,4] = np.nan
            data_in[i,5] = np.nan
            data_in[i,10] = np.nan
    return data_in

def classify_modes(data_in,pars):
    if 'ExBrate' in pars and pars['ExBrate'] != 0.0:
        wExBshear = True
    else:
        wExBshear = False
    ntot = len(data_in[:,0])
    modes = []
    epar_threshold = 0.2
    epar_threshold2 = 0.6
    if 'x_local' in pars and not pars['x_local']:
        epar_threshold = 0.4
        epar_threshold2 = 0.7
    dz = 2.0/float(pars['nz0'])
    print(( "epar_threshold",epar_threshold))
    print(( "epar_threshold2",epar_threshold2))
    for i in range(ntot):
        #Test for numerical modes (i.e. grid-scale z)
        if data_in[i,7] <= 2*dz:
            if data_in[i,8]/data_in[i,9] > 1.0:
                #Numerical microtearing
                modes.append('NMTM') 
            else:
                #Numerical ballooning parity
                modes.append('NBP')  
        #Test for MTM
        #Negative frequency
        #Mostly tearing parity
        #Mostly EM transport
        if data_in[i,10] > 0.5 and \
           data_in[i,8] / data_in[i,9] >= 0.8 and \
           len(modes) == i:
            if not wExBshear:
                if data_in[i,5] < 0.0:
                    modes.append('MTM')  
            else:
                modes.append('MTM')  
              
        #Test for KBM
        #Positive frequency
        #Peaked mode structure
        #Mostly ballooning parity
        #Epar cancelation 
        if data_in[i,8] / data_in[i,9] < 0.1 and \
           data_in[i,11] < epar_threshold and \
           len(modes) == i:
            if not wExBshear:
                if data_in[i,5] > 0.0:
                    modes.append('KBM')  
            else:
                modes.append('KBM')  
        if data_in[i,11] < epar_threshold2 and \
           len(modes) == i:
            modes.append('ID')
        if len(modes) == i:
            modes.append('other')
    print(( "len(modes)",len(modes)))
    for i in range(len(modes)):
        print(( i+1,modes[i]))
    return modes

def get_grids(data_in):
    numscan_tot = len(data_in[:,0])
    num_col = len(data_in[0,:])
    ky_array = np.empty(0)
    x0_array = np.empty(0)
    for i in range(numscan_tot):
        if data[i,0] not in ky_array:
            ky_array = np.append(ky_array,data[i,0])
        if data[i,1] not in x0_array and not np.isnan(data[i,1]):
            x0_array = np.append(x0_array,data[i,1])
 
    if len(x0_array) == 0:
        x0_array = np.append(x0_array,0.0)
    nkxc = int(numscan_tot/(len(ky_array)*len(x0_array)) )
    kxc_array = np.empty((len(ky_array),nkxc))
    kxc_array[:,:] = -1.0
    fill_index = np.empty(len(ky_array),dtype='int')
    fill_index[:] = -1
    for i in range(numscan_tot):
        kyind = np.argmin(abs(data_in[i,0]-ky_array))
        if data_in[i,2] not in kxc_array[kyind,:]:
            fill_index[kyind] += 1
            print(( 'kyind',kyind))
            kxc_array[kyind,fill_index[kyind]] = data_in[i,2]

    #print 'fill_index',fill_index
    return ky_array, x0_array, kxc_array

def get_kxc_array(data_in,ky_array,x0_array):
    numscan_tot = len(data_in[:,0])
    num_col = len(data_in[0,:])
    

def get_ky_arrays(data_in):
    numscan_tot = len(data_in[:,0])
    num_col = len(data_in[0,:])
    ky_array = np.empty(0)
    for i in range(numscan_tot):
        if data[i,0] not in ky_array:
            ky_array = np.append(ky_array,data[i,0])
    data_ky = np.empty((len(ky_array),int(numscan_tot/len(ky_array)),num_col))
    index_ky = np.empty(len(ky_array),dtype='int')
    index_ky[:] = -1
    for i in range(numscan_tot):
        index = np.argmin(abs(ky_array-data[i,0]))
        index_ky[index] += 1
        data_ky[index,index_ky[index],:] = data[i,:]
    
    return data_ky 

#def calc_pedestal_average_gHB(rbs):
#    rmin = float(raw_input("Enter xmin for pedestal average:\n"))
#    rmax = float(raw_input("Enter xmax for pedestal average:\n"))
#    rhot = np.arange(2000)/1999.0
#    ghb = interp(rbs[:,0],rbs[:,9],rhot)
#    imin = np.argmin(abs(rhot-rmin))
#    imax = np.argmin(abs(rhot-rmax))
#    ghb_avg = np.sum(abs(ghb[imin:imax]))/float(imax-imin)   
#    return ghb_avg,rmin,rmax


#def get_x0_arrays(data_in):
#    numscan_tot = len(data_in[:,0])
#    num_col = len(data_in[0,:])
#    x0_array = np.empty(0)
#    for i in range(numscan_tot):
#        if data[i,1] not in x0_array:
#            x0_array = np.append(x0_array,data[i,1])
#    data_x0 = np.empty((len(x0_array),numscan_tot/len(x0_array),num_col))
#    index_x0 = np.empty(len(x0_array),dtype='int')
#    index_x0[:] = -1
#    for i in range(numscan_tot):
#        index = np.argmin(abs(x0_array-data[i,1]))
#        index_x0[index] += 1
#        data_x0[index,index_x0[index],:] = data[i,:]
#    return x0_array, data_x0 


ky_array, x0_array, kxc_array = get_grids(data)
data_ky = get_ky_arrays(data)

modes = classify_modes(data,pars)

if len(x0_array) > 1:
    plt.contourf(x0_array,ky_array,data_ky[:,:,4],50)
    plt.colorbar()
    for i in range(len(data[:,0])):
        if modes[i] == 'MTM':
            mark = 'x'
            col = 'red'
        elif modes[i] == 'KBM':
            mark = 'o'
            col = 'black'
        elif modes[i] == 'ID':
            mark = 'o'
            col = 'magenta'
        elif modes[i] == 'NMTM' or modes[i] == 'NBP':
            mark = '+'
            col = 'orange'
        else:
            mark = '+'
            col = 'green'
        plt.scatter(data[i,1],data[i,0],marker=mark,color=col)
    plt.show()

    plt.contourf(x0_array,ky_array,data_ky[:,:,11],50)
    plt.colorbar()
    plt.show()

    fig=plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    for i in range(len(ky_array)):
        ax.plot(x0_array[:],data_ky[i,:,0],data_ky[i,:,4],color = 'red',marker = 'x')
    #for i in range(len(x0_array)):
    #    ax.plot(data_x0[i,:,1],ky_array[:],data_x0[i,:,3],color = 'blue',marker = 'x')
    ax.set_xlabel(r'$\rho_{tor}$')
    ax.set_ylabel(r'$k_y \rho_s$')
    ax.set_zlabel(r'$\gamma (c_s/a)$')
    plt.show()

    fig=plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    for i in range(len(ky_array)):
        ax.plot(x0_array[:],data_ky[i,:,0],data_ky[i,:,11],color = 'red',marker = 'x')
    #for i in range(len(x0_array)):
    #    ax.plot(data_x0[i,:,1],ky_array[:],data_x0[i,:,10],color = 'blue',marker = 'x')
    ax.set_xlabel(r'$\rho_{tor}$')
    ax.set_ylabel(r'$k_y \rho_s$')
    #ax.set_zlabel(r'$\gamma (c_s/a)$')
    plt.show()

dz = 2.0/float(pars['nz0'])
if len(kxc_array[0,:]) > 1:

    data = filter_numerical_modes(data,pars)    
    if get_global:
       data = filter_local_modes(data,pars)

    kxcmax = kxc_array[-1,-1]
    kymax = ky_array[-1]
    print( 'kxcmax')
    #kxgrid = np.arange(100)/99.0*kxcmax
    #kygrid = np.arange(100)/99.0*kymax
    #zi = griddata((data[:,2],data[:,0]),data[:,4],(kxgrid,kygrid))#,interp='linear')
    #z2 = griddata((data[:,2],data[:,0]),data[:,10],(kxgrid,kygrid))#,interp='linear')
    #print(np.shape(zi))
    #plt.contourf(kxgrid,kygrid,zi,50,cmap=plt.cm.RdYlBu)
    print((data[:,4]))
    grs = np.nan_to_num(data[:,4],nan=0.0)
    print(grs)
    plt.tricontourf(data[:,2],data[:,0],grs,levels=50)
    plt.colorbar()
    for i in range(len(data[:,0])):
        if modes[i] == 'MTM':
            mark = 'x'
            col = 'red'
        elif modes[i] == 'KBM':
            mark = 'o'
            col = 'black'
        elif modes[i] == 'ID':
            mark = 'o'
            col = 'magenta'
        elif modes[i] == 'NMTM' or modes[i] == 'NBP':
            mark = '+'
            col = 'orange'
        else:
            mark = '+'
            col = 'green'
        plt.scatter(data[i,2],data[i,0],marker=mark,color=col)
    plt.xlabel(r'$k_{x,c} \rho_s$',size=18)
    plt.ylabel(r'$k_{y} \rho_s$',size=18)

    plt.show()

    #Find maximum gr per ky
    max_gr_ky = np.zeros(len(ky_array))
    max_gr_ky_fr = np.zeros(len(ky_array))
    modes_ky = [None]*len(ky_array)
    for i in range(len(modes_ky)):
        modes_ky[i] = 'NBP'
    for i in range(len(data[:,0])):
        kyind = np.argmin(abs(ky_array-data[i,0]))
        if data[i,4] > max_gr_ky[kyind]:
            max_gr_ky[kyind] = data[i,4]
            max_gr_ky_fr[kyind] = data[i,5]
            modes_ky[kyind] = modes[i]
    glob_string = ''
    fglob_string = ''
    if get_global:
       glob_string += '   (Filtered for global)'
       fglob_string += '_gfiltered'
    
    f=open('gamma_omega_ky'+fglob_string+'.dat','w')
    f.write('# x0 = '+str(pars['x0'])+glob_string+'\n')
    f.write('# diagdir ='+pars['diagdir']+'\n')
    np.savetxt(f,np.column_stack((ky_array,max_gr_ky,max_gr_ky_fr)))
    f.close()
    
    plt.plot(ky_array,max_gr_ky)
    f=open('modes_ky'+fglob_string+'.dat','w')
    f.write('# x0 = '+str(pars['x0'])+glob_string+'\n')
    f.write('# diagdir ='+pars['diagdir']+'\n')
    for i in range(len(ky_array)):
        if modes_ky[i] == 'MTM':
            mark = 'x'
            col = 'red'
        elif modes_ky[i] == 'KBM':
            mark = 'o'
            col = 'black'
        elif modes_ky[i] == 'ID':
            mark = (5,2,0)
            col = 'magenta'
        elif modes_ky[i] == 'NMTM' or modes_ky[i] == 'NBP':
            mark = '+'
            col = 'orange'
        else:
            mark = '+'
            col = 'green'
        plt.plot(ky_array[i],max_gr_ky[i],color=col,marker=mark,markeredgewidth=2,markersize=5)
        #plt.plot(ky_array[i],max_gr_ky_fr[i],color=col,marker=mark,markeredgewidth=2,markersize=5)
        print(( "modes_ky[i]",modes_ky[i]))
        f.write(str(ky_array[i])+'\t'+modes_ky[i]+'\n')
      
    f.close()
    plt.xlabel(r'$k_y \rho_s$',size=18)
    plt.ylabel(r'$\gamma (c_s/a)$',size=18)
    plt.show()
    

if 'x_local' not in pars or pars['x_local']:
    print( "Local scan detected.")
    data = filter_numerical_modes(data,pars)    
    
    f=open('modes_ky.dat','w')
    f.write('# x0 = '+str(pars['x0'])+'\n')
    f.write('# diagdir ='+pars['diagdir']+'\n')
    ky_array = data[:,0]
    for i in range(len(data[:,0])):
        if modes[i] == 'MTM':
            mark = 'x'
            col = 'red'
        elif modes[i] == 'KBM':
            mark = 'o'
            col = 'black'
        elif modes[i] == 'ID':
            mark = (5,2,0)
            col = 'magenta'
        elif modes[i] == 'NMTM' or modes[i] == 'NBP':
            mark = '+'
            col = 'orange'
        else:
            mark = '+'
            col = 'green'
        plt.plot(ky_array[i],data[i,4],color=col,marker=mark,markeredgewidth=2,markersize=5)
        print(( "modes[i]",modes[i]))
        f.write(str(ky_array[i])+'\t'+modes[i]+'\n')
      
    f.close()
    plt.xlabel(r'$k_y \rho_s$',size=18)
    plt.ylabel(r'$\gamma (c_s/a)$',size=18)
    plt.show()
else:

    homedir = os.path.expanduser("~")
    thisdir = os.getcwd()
    case = thisdir.split('/')[-2]
    #rbs = np.genfromtxt(homedir+'/pmv_eqs/'+case+'/rbsProfs')
    #ghbpa, rmin, rmax = calc_pedestal_average_gHB(rbs) 
    #ind_x0 = np.argmin(abs(rbs[:,0]-data[0,1]))
    #gamma_HB_x0 = abs(rbs[ind_x0,9])
    #print "Local HB shear rate:",gamma_HB_x0
    #print "Pedestal averaged HB shear rate:",ghbpa


