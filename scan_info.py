#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from ParIO import *
import os
from fieldlib import *
from get_nrg import get_nrg0
import sys
from calc_gr import *
from finite_differences import *
from read_write_geometry import *
import optparse as op
from subprocess import call
from get_abs_psi_prime import *
from remove_edge_opt import *

parser=op.OptionParser(description='Calculates mode information and synthesizes scan info.')
options,args=parser.parse_args()
if len(args)!=1:
    exit("""
Please include scan number as argument (e.g., 0001)."
    \n""")
suffix = args[0]

calc_grs = True
os.chdir('scanfiles'+suffix)


def my_corr_func_complex(v1,v2,time,show_plot=False,v1eqv2=True):
    #print "len(time)",len(time)
    #print "len(v1)",len(v1)
    #print "len(v2)",len(v2)
    dt=time[1]-time[0]
    #print "dt:", dt
    N=len(time)
    cfunc=np.zeros(N,dtype='complex')
    for i in range(N):
        i0=i+1
        cfunc[-i0]=np.sum(np.conj(v1[-i0:])*v2[:i0])
    tau=np.arange(N)
    tau=tau*dt
    if v1eqv2:
        cfunc=np.real(cfunc)
    max_corr=max(np.abs(cfunc))
    corr_time=0.0
    i=0
    while corr_time==0.0:
        if (abs(cfunc[i])-max_corr/np.e) > 0.0 and \
           (abs(cfunc[i+1])-max_corr/np.e) <= 0.0:
            slope=(cfunc[i+1]-cfunc[i])/(tau[i+1]-tau[i])
            zero=cfunc[i]-slope*tau[i]
            corr_time=(max_corr/np.e-zero)/slope
        i+=1
    neg_loc = 10000.0
    i=0
    while neg_loc==10000.0 and i < N:
        if cfunc[i] < 0.0:
            neg_loc = tau[i]
        i+=1

    if neg_loc < corr_time:
        print "WARNING: neg_loc < corr_time"
        corr_time = neg_loc

    if show_plot:
        plt.plot(tau,cfunc,'x-')
        ax=plt.axis()
        plt.vlines(corr_time,ax[2],ax[3])
        plt.show()
    return cfunc,tau,corr_time

par = Parameters()
par.Read_Pars('parameters')
pars = par.pardict
edge_opt = pars['edge_opt']
print "edge_opt = ",edge_opt
#dummy = raw_input("Press any key to continue:\n")

print type(pars['scan_dims'])
print "scan_dims",pars['scan_dims']
if type(pars['scan_dims']) == str:
    scan_dims = pars['scan_dims'].split()
    numscan_tot = 1
    for i in scan_dims:
        numscan_tot *= int(i)
else:
    numscan_tot = pars['scan_dims']

print "Total number of runs: ", numscan_tot 

if pars['n_spec'] == 3:
    dummy = raw_input("Assuming electrons are third species (press any key).\n")

#Test if global scan
if 'x_local' in pars and not pars['x_local']:
    scan_info = np.zeros((numscan_tot,14),dtype='float64')
    for i in range(numscan_tot):
        par0 = Parameters()
        scan_num = '000'+str(i+1)
        scan_num = scan_num[-4:]
        print "Analyzing ",scan_num
        if os.path.isfile('parameters_'+scan_num):
            par0.Read_Pars('parameters_'+scan_num)
            pars0 = par0.pardict
            nspec = pars0['n_spec']
            print pars0['kymin']
            scan_info[i,0] = pars0['kymin']
            if 'x0' in pars0:
                scan_info[i,1] = pars0['x0']
            else:
                scan_info[i,1] = pars['x0']
            scan_info[i,2] = 0.0
            if 'n0_global' in pars0:
                scan_info[i,3] = pars0['n0_global']
            else:
                scan_info[i,3] = np.nan
        else:
            par0.Read_Pars('in_par/parameters_'+scan_num)
            pars0 = par0.pardict
            nspec = pars0['n_spec']
            scan_info[i,0] = float(str(pars0['kymin']).split()[0])
            if 'n0_global' in pars0:
                scan_info[i,3] = float(str(pars0['n0_global']).split()[0])
            else:
                scan_info[i,3] = np.nan
            scan_info[i,1] = pars0['x0']
            scan_info[i,2] = 0.0
        if os.path.isfile('omega_'+scan_num):
            omega0 = np.genfromtxt('omega_'+scan_num)
            if omega0.any() and omega0[1] != 0.0:
                scan_info[i,4]=omega0[1]
                scan_info[i,5]=omega0[2]
            elif calc_grs:
                call(['calc_omega_from_field.py',scan_num])
                om = np.genfromtxt('omega_'+scan_num)
                scan_info[i,4]=om[1]
                scan_info[i,5]=om[2]
            else:
                scan_info[i,4]=np.nan
                scan_info[i,5]=np.nan
        elif calc_grs and os.path.isfile('field_'+scan_num):
            call(['calc_omega_from_field.py',scan_num])
            om = np.genfromtxt('omega_'+scan_num)
            scan_info[i,4]=om[1]
            scan_info[i,5]=om[2]
        else:
            scan_info[i,4]=np.nan
            scan_info[i,5]=np.nan
        if os.path.isfile('field_'+scan_num):
            field = fieldfile('field_'+scan_num,pars0)
            #field.set_time(field.tfld[-1],len(field.tfld[-1])-1)
            field.set_time(field.tfld[-1])
            fntot = field.nz*field.nx
    
            dz = float(2.0)/float(field.nz)
            zgrid = np.arange(field.nz)/float(field.nz-1)*(2.0-dz)-1.0
            zgrid_ext = np.arange(field.nz+4)/float(field.nz+4-1)*(2.0+3*dz)-(1.0+2.0*dz)
            #field.set_time(field.tfld[-1],len(field.tfld[-1])-1)
            field.set_time(field.tfld[-1])
       
            imax = np.unravel_index(np.argmax(abs(field.phi()[:,0,:])),(field.nz,field.nx))
            imaxa = np.unravel_index(np.argmax(abs(field.apar()[:,0,:])),(field.nz,field.nx))
            phi = field.phi()[:,0,:]
            apar = field.apar()[:,0,:] 
            cfunc,zed,corr_len=my_corr_func_complex(phi[:,imax[1]],phi[:,imax[1]],zgrid,show_plot=False)
            scan_info[i,6] = np.nan
            scan_info[i,7] = corr_len
            parity_factor_apar = np.abs(np.sum(apar[:,imaxa[1]]))/np.sum(np.abs(apar[:,imaxa[1]]))
            scan_info[i,8] = parity_factor_apar
            parity_factor_phi = np.abs(np.sum(phi[:,imax[1]]))/np.sum(np.abs(phi[:,imax[1]]))
            scan_info[i,9] = parity_factor_phi
    
            #KBM test with E||

            dz = 2.0/field.nz
            zgrid_ext = np.arange(field.nz+4)/float(field.nz+4-1)*(2.0+3*dz)-(1.0+2.0*dz)
            xgrid = np.arange(field.nx)/float(field.nx-1)*pars0['lx_a']+pars0['x0']-pars0['lx_a']/2.0
            gpars,geometry = read_geometry_global(pars0['magn_geometry'][1:-1]+'_'+scan_num)
            phase = (0.0+1.0J)*pars0['n0_global']*2.0*np.pi*geometry['q']
            phi_bnd = np.zeros((field.nz+4,field.nx),dtype = 'complex128')     
            gradphi= np.zeros((field.nz+4,field.nx),dtype = 'complex128')     
            phi_bnd[2:-2,:] = phi
            for j in range(field.nx):
                phi_bnd[-2,j] = phi_bnd[2,j]*np.e**(-phase[j])
                phi_bnd[-1,j] = phi_bnd[3,j]*np.e**(-phase[j])
                phi_bnd[0,j] = phi_bnd[-4,j]*np.e**(phase[j])
                phi_bnd[1,j] = phi_bnd[-3,j]*np.e**(phase[j])
                gradphi[:,j] = fd_d1_o4(phi_bnd[:,j],zgrid_ext)
                gradphi[2:-2,j] = gradphi[2:-2,j]/np.pi/(geometry['jacobian'][:,j]*geometry['Bfield'][:,j])
        
            #field.set_time(field.tfld[-1],len(field.tfld[-1])-1)
            field.set_time(field.tfld[-1])
        
            #apar = field.apar()[:,:]
        
            omega_complex = (scan_info[i,5]*(0.0+1.0J) + scan_info[i,4])
        
            #fig=plt.gcf()
            #fig.subplots_adjust(right=0.9)
            #fig.subplots_adjust(left=0.16)
            #fig.subplots_adjust(hspace=0.35)
            #plt.subplot(3,1,1)
            #plt.plot(zgrid,(np.real(gradphi[2:-2,0,field.nx/4])),color = 'black')
            #plt.plot(zgrid,(np.real(-omega_complex*field.apar()[:,0,field.nx/4])),color = 'red')
            #plt.plot(zgrid,(np.imag(gradphi[2:-2,0,field.nx/4])),'-.',color = 'black')
            #plt.plot(zgrid,(np.imag(-omega_complex*field.apar()[:,0,field.nx/4])),'-.',color = 'red')
            #plt.subplot(3,1,2)
            #plt.plot(zgrid,(np.real(gradphi[2:-2,0,field.nx/2])),color = 'black')
            #plt.plot(zgrid,(np.real(-omega_complex*field.apar()[:,0,field.nx/2])),color = 'red')
            #plt.plot(zgrid,(np.imag(gradphi[2:-2,0,field.nx/2])),'-.',color = 'black')
            #plt.plot(zgrid,(np.imag(-omega_complex*field.apar()[:,0,field.nx/2])),'-.',color = 'red')
            #plt.subplot(3,1,3)
            #plt.plot(zgrid,(np.real(gradphi[2:-2,0,3*field.nx/4])),color = 'black')
            #plt.plot(zgrid,(np.real(-omega_complex*field.apar()[:,0,3*field.nx/4])),color = 'red')
            #plt.plot(zgrid,(np.imag(gradphi[2:-2,0,3*field.nx/4])),'-.',color = 'black')
            #plt.plot(zgrid,(np.imag(-omega_complex*field.apar()[:,0,3*field.nx/4])),'-.',color = 'red')
            #plt.show()
            print "omega_complex",omega_complex
            
            diff = np.sum(np.abs(gradphi[2:-2,:] + omega_complex*apar[:,:]))
            phi_cont = np.sum(np.abs(gradphi[2:-2,:]))
            apar_cont = np.sum(np.abs(omega_complex*apar[:,:]))
            print "diff",diff
            print "phi_cont",phi_cont
            print "apar_cont",apar_cont
            print "diff/abs",diff/(phi_cont+apar_cont)
            scan_info[i,11] = diff/(phi_cont+apar_cont)
            scan_info[i,12] = np.nan           
   
        else:
            scan_info[i,6] = np.nan
            scan_info[i,7] = np.nan
            scan_info[i,8] = np.nan
            scan_info[i,9] = np.nan
            scan_info[i,11] = np.nan
            scan_info[i,12] = np.nan
            scan_info[i,13] = np.nan
    
        if os.path.isfile('nrg_'+scan_num):
            if nspec==1:
                tn,nrg1=get_nrg0('_'+scan_num,nspec=nspec)
                scan_info[i,10]=nrg1[-1,7]/abs(nrg1[-1,6])
            elif nspec==2:
                tn,nrg1,nrg2=get_nrg0('_'+scan_num,nspec=nspec)
                scan_info[i,10]=nrg2[-1,7]/(abs(nrg2[-1,6])+abs(nrg1[-1,6]))
            elif nspec==3:
                tn,nrg1,nrg2,nrg3=get_nrg0('_'+scan_num,nspec=nspec)
                scan_info[i,10]=nrg3[-1,7]/(abs(nrg3[-1,6])+abs(nrg1[-1,6]))
            else:
                sys.exit("Not ready for nspec>2")
        else:
            scan_info[i,10] = np.nan
    
    f=open('scan_info.dat','w')
    p = open('parameters','r')
    lines = p.read()
    p.close()
    lines = lines.split('\n')
    for line in lines:
        if 'scan' in line:
            f.write('#'+line+'\n')
    f.write('#1.kymin 2.x0 3.kx_center 4.n0_global 5.gamma(cs/a) 6.omega(cs/a) 7.<z> 8.lambda_z 9.parity(apar) 10.parity(phi) 11.QEM/QES 12.Epar cancelation \n')
    np.savetxt(f,scan_info)
    f.close()
#local scan
else:

    scan_info = np.zeros((numscan_tot,14),dtype='float64')
    
    for i in range(numscan_tot):
        par0 = Parameters()
        scan_num = '000'+str(i+1)
        scan_num = scan_num[-4:]
        print "Analyzing ",scan_num
        if os.path.isfile('parameters_'+scan_num):
            par0.Read_Pars('parameters_'+scan_num)
            pars0 = par0.pardict
            nspec = pars0['n_spec']
            print pars0['kymin']
            scan_info[i,0] = pars0['kymin']
            if 'x0' in pars0:
                scan_info[i,1] = pars0['x0']
            else:
                scan_info[i,1] = pars['x0']
            if 'kx_center' in pars0:
                scan_info[i,2] = pars0['kx_center']
            else:
                scan_info[i,2] = 0.0
            if 'n0_global' in pars0:
                scan_info[i,3] = pars0['n0_global']
            else:
                scan_info[i,3] = np.nan
        else:
            par0.Read_Pars('in_par/parameters_'+scan_num)
            pars0 = par0.pardict
            nspec = pars0['n_spec']
            scan_info[i,0] = float(str(pars0['kymin']).split()[0])
            scan_info[i,1] = float(str(pars0['x0']).split()[0])
            if 'kx_center' in pars0:
                scan_info[i,2] = float(str(pars0['kx_center']).split()[0])
            else:
                scan_info[i,2] = 0.0
            if 'n0_global' in pars0:
                scan_info[i,3] = pars0['n0_global']
            else:
                scan_info[i,3] = 0.0
        if os.path.isfile('omega_'+scan_num):
            omega0 = np.genfromtxt('omega_'+scan_num)
            if omega0.any() and omega0[1] != 0.0:
                scan_info[i,4]=omega0[1]
                scan_info[i,5]=omega0[2]
            elif calc_grs:
                scan_info[i,4]=calc_gr2('_'+scan_num,nspec=nspec)
                scan_info[i,5]= 0.0
                np.savetxt('omega_'+scan_num,[scan_info[i,0],scan_info[i,4],np.nan])
            else:
                scan_info[i,4]=np.nan
                scan_info[i,5]=np.nan
        elif calc_grs and os.path.isfile('nrg_'+scan_num):
            scan_info[i,4]=calc_gr2('_'+scan_num,nspec=nspec)
            scan_info[i,5] = 0.0
            np.savetxt('omega_'+scan_num,[scan_info[i,0],scan_info[i,4],np.nan])
        else:
            scan_info[i,4]=np.nan
            scan_info[i,5]=np.nan
        
        if os.path.isfile('field_'+scan_num):
            field = fieldfile('field_'+scan_num,pars0)
            #field.set_time(field.tfld[-1],len(field.tfld)-1)
            field.set_time(field.tfld[-1])
            fntot = field.nz*field.nx
    
            dz = float(2.0)/float(field.nz)
            zgrid = np.arange(fntot)/float(fntot-1)*(2*field.nx-dz)-field.nx
            zgrid0 = np.arange(field.nz)/float(field.nz-1)*(2.0-dz)-1.0
            phi = np.zeros(fntot,dtype='complex128')
            apar = np.zeros(fntot,dtype='complex128')
            phikx = field.phi()[:,0,:]
            aparkx = field.phi()[:,0,:]
            phase_fac = -np.e**(-2.0*np.pi*(0.0+1.0J)*pars0['n0_global']*pars0['q0'])
            for j in range(field.nx/2+1):
                phi[(j+field.nx/2)*field.nz:(j+field.nx/2+1)*field.nz]=field.phi()[:,0,j]*phase_fac**j
                if j < field.nx/2:
                    phi[(field.nx/2-j-1)*field.nz : (field.nx/2-j)*field.nz ]=field.phi()[:,0,-1-j]*phase_fac**(-(j+1))
                if pars0['n_fields']>1:
                    apar[(j+field.nx/2)*field.nz:(j+field.nx/2+1)*field.nz]=field.apar()[:,0,j]*phase_fac**j
                    if j < field.nx/2:
                        apar[(field.nx/2-j-1)*field.nz : (field.nx/2-j)*field.nz ]=field.apar()[:,0,-1-j]*phase_fac**(-(j+1))
        
            zavg=np.sum(np.abs(phi)*np.abs(zgrid))/np.sum(np.abs(phi))
            scan_info[i,6] = zavg
            cfunc,zed,corr_len=my_corr_func_complex(phi,phi,zgrid,show_plot=False)
            scan_info[i,7] = corr_len
            parity_factor_apar = np.abs(np.sum(apar))/np.sum(np.abs(apar))
            scan_info[i,8] = parity_factor_apar
            parity_factor_phi = np.abs(np.sum(phi))/np.sum(np.abs(phi))
            scan_info[i,9] = parity_factor_phi
    
            #KBM test with E||
            gpars,geometry = read_geometry_local(pars0['magn_geometry'][1:-1]+'_'+scan_num)
            jacxB = geometry['gjacobian']*geometry['gBfield']
            if scan_info[i,5] == scan_info[i,5]:
                omega_complex = (scan_info[i,5]*(0.0+1.0J) + scan_info[i,4])
                gradphi = fd_d1_o4(phi,zgrid)
                for j in range(pars0['nx0']):
                    gradphi[pars0['nz0']*j:pars0['nz0']*(j+1)] = gradphi[pars0['nz0']*j:pars0['nz0']*(j+1)]/jacxB[:]/np.pi
                #plt.plot(gradphi)
                #plt.plot(omega_complex*apar)
                #plt.show()
                diff = np.sum(np.abs(gradphi + omega_complex*apar))
                phi_cont = np.sum(np.abs(gradphi))
                apar_cont = np.sum(np.abs(omega_complex*apar))
                scan_info[i,11] = diff/(phi_cont+apar_cont)
                #print 'diff/abs',scan_info[i,10]
            else:
                scan_info[i,11] = np.nan
            phi0 = np.empty(np.shape(phikx),dtype = 'complex') 
            apar0 = np.empty(np.shape(aparkx),dtype = 'complex') 
            #print "Shape of phikx",np.shape(phikx)
            #print "Shape of zgrid",np.shape(zgrid)
            #dummy = raw_input("Press any key:\n")
            #if edge_opt > 0.0: 
            #    for ix in range(len(phikx[0,:])):
            #        phi0[:,ix] = remove_edge_opt_complex(phikx[:,ix],edge_opt)
            #        apar0[:,ix] = remove_edge_opt_complex(aparkx[:,ix],edge_opt)
            phi0 = phikx
            apar0 = aparkx
            #Calculate <gamma_HB> / gamma
            geomfile = pars0['magn_geometry'][1:-1]+'_'+scan_num
            print "geomfile",geomfile
            zgrid_pp, Btheta_R, prefactor = get_abs_psi_prime(geomfile,'../rbsProfs',pars['x0'])
            rbs = np.genfromtxt('../rbsProfs')
            ind_rbs_x0 = np.argmin(abs(rbs[:,0]-pars['x0'])) 
            gamma_HB_norm_x0 = rbs[ind_rbs_x0,9]
            ind_z0 = np.argmin(abs(zgrid_pp)) 
            prefactor_norm = prefactor/prefactor[ind_z0]
            gamma_HB_theta = abs(gamma_HB_norm_x0*prefactor_norm)
            #plt.plot(zgrid_pp,gamma_HB_theta)
            #plt.xlabel(r'$z/\pi$',size=18)
            #plt.ylabel(r'$\gamma_{HB}(c_s/a)$',size=18)
            #plt.show()
            gamma_HB_sum = 0.0
            phi_sum = 0.0
            for ix in range(len(phi0[0,:])):
                gamma_HB_sum += np.sum(abs(phi0[:,ix])**2*gamma_HB_theta*geometry['gjacobian'])
                phi_sum += np.sum(abs(phi0[:,ix])**2*geometry['gjacobian'])
            gamma_HB_avg = gamma_HB_sum / phi_sum
            scan_info[i,12] = gamma_HB_avg
            #gamma_HB_theta = abs(gamma_HB_norm_x0*prefactor_norm)
            #theta0 = (scan_info[i,2]/(pars0['shat']*pars0['kymin']*np.pi))
            #if theta0 > 1.0:
            #    theta0 -= 2.0
            #ind_theta0 = np.argmin(abs(zgrid_pp-theta0)) 
            #print 'kx_center',scan_info[i,2]
            #print 'shat',pars0['shat']
            #print 'kymin',pars0['kymin']
            #print 'theta0',theta0
            #print 'ind_theta0',ind_theta0
 	    scan_info[i,13] = np.min(gamma_HB_theta)
        else:
            scan_info[i,6] = np.nan
            scan_info[i,7] = np.nan
            scan_info[i,8] = np.nan
            scan_info[i,9] = np.nan
            scan_info[i,11] = np.nan
            scan_info[i,12] = np.nan
            scan_info[i,13] = np.nan
    
        if os.path.isfile('nrg_'+scan_num):
            if nspec==1:
                tn,nrg1=get_nrg0('_'+scan_num,nspec=nspec)
                scan_info[i,10]=nrg1[-1,7]/abs(nrg1[-1,6])
            elif nspec==2:
                tn,nrg1,nrg2=get_nrg0('_'+scan_num,nspec=nspec)
                scan_info[i,10]=nrg2[-1,7]/(abs(nrg2[-1,6])+abs(nrg1[-1,6]))
            elif nspec==3:
                tn,nrg1,nrg2,nrg3=get_nrg0('_'+scan_num,nspec=nspec)
                scan_info[i,10]=nrg3[-1,7]/(abs(nrg3[-1,6])+abs(nrg1[-1,6]))
            else:
                sys.exit("Not ready for nspec>2")
        else:
            scan_info[i,10] = np.nan
    
    f=open('scan_info.dat','w')
    p = open('parameters','r')
    lines = p.read()
    p.close()
    lines = lines.split('\n')
    for line in lines:
        if 'scan' in line:
            f.write('#'+line+'\n')
    f.write('#1.kymin 2.x0 3.kx_center 4.n0_global 5.gamma(cs/a) 6.omega(cs/a) 7.<z> 8.lambda_z 9.parity(apar) 10.parity(phi) 11.QEM/QES 12.Epar cancelation 13.gamma_HB_avg 14.gamma_HB_min \n')
    np.savetxt(f,scan_info)
    f.close()

os.chdir('../')
call(['./plot_scan_info.py',suffix])


