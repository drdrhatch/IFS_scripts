import numpy as np
from read_write_geometry import *
from interp import *
from finite_differences import *
from remove_edge_opt import *
from rbs_tools import *

def get_abs_psi_prime(geomfile,rbsfile,x0,show_plots = False, do_tests = True, remove_edge_opt = False, edge_opt = 0.0):

    print( "Calculating absolution value of psi_prime (i.e. Bpol R).")
    print( "x0 (rhot):",x0)
    #Note edge_opt is not removed--must include jacobian in integrals

    rbs = np.genfromtxt(rbsfile)
    
    parameters, geom = read_geometry_local(geomfile)
    psi0 = abs(get_psi0(rbsfile))
    
    nz0 = float(parameters['gridpoints'])
    dz = 2.0/nz0
    zgrid_even = np.arange(nz0)/float(nz0-1)*(2.0-dz)-1.0
    # zgrid is even zgrid of length nz0
    # gl_dxdR = drhot / dR on zgrid un-normalized
    # gl_dxdZ = drhot / dZ on zgrid un-normalized
    # gl_R = R on zgrid 
    # gB = Total B on zgrid un-normalized
    if remove_edge_opt:
        if edge_opt <= 0.0:
            edge_opt = float(raw_input("Enter edge_opt:\n"))    
        print( "Removing edge_opt.")
        gl_dxdR = remove_edge_opt(geom['gl_dxdR'],edge_opt)
        gl_dxdZ = remove_edge_opt(geom['gl_dxdZ'],edge_opt)
        gl_R = remove_edge_opt(geom['gl_R'],edge_opt)
        gl_z = remove_edge_opt(geom['gl_z'],edge_opt)
        gB = remove_edge_opt(geom['gBfield'],edge_opt)
    else:
        gl_dxdR = geom['gl_dxdR']
        gl_dxdZ = geom['gl_dxdZ']
        gl_R = geom['gl_R']
        gl_z = geom['gl_z']
        gB = geom['gBfield']
    
    zgrid = zgrid_even
    gB = gB*parameters['Bref']
    gl_dxdZ = gl_dxdZ/parameters['Lref']
    gl_dxdR = gl_dxdR/parameters['Lref']
    
    #drhot_dr = drhot / dr (where dr is radial perpendicular to flux surface) on zgrid
    drhot_dr = (gl_dxdR**2+gl_dxdZ**2)**0.5
    
    rhot_even = np.arange(1000)/999.0
    
    psi_rhot = interp(rbs[:,0],rbs[:,1],rhot_even)
    #dpsi_drhot = dpsi_drhot on even rhot grid
    dpsi_drhot = fd_d1_o4(psi_rhot,rhot_even)
    
    #ind = rhot index of radial position (on even rhot grid)
    ind = np.argmin(abs(rhot_even - x0))
    #dpsi_drhot0 = dpsi_drhot at x0
    dpsi_drhot0 = dpsi_drhot[ind]
    
    #dpsi_dr = R B_theta on zgrid 
    dpsi_dr = drhot_dr*dpsi_drhot0*psi0
    if show_plots:
        plt.plot(zgrid,dpsi_dr,'x-')
        plt.xlabel('z(rad)')
        plt.ylabel(r'$\frac{d \psi}{ d r}$',size=18)
        plt.show()
        np.savetxt('dpsidr.dat',np.column_stack((zgrid,dpsi_dr)))
    
    if do_tests:
        #####Tests#####
        #####Tests#####
        #####Tests#####
    
        #ind0 = index of outboard midplane on zgrid
        ind0 = np.argmin(abs(zgrid-0))
        #ind_gexb = index of x0 on rbs grid
        x0_ind_rbs = np.argmin(abs(rbs[:,0]-x0))
        #prefactor = should be B_theta**2 R**2 / B
        prefactor = dpsi_dr**2/gB
        #prefactor_norm = prefactor normalized to outboard midplane value
        prefactor_norm = prefactor/prefactor[ind0]
        
        #Test 1:
        #Calculate R**2 B_theta**2 / B
        #Compare outboard midplane vs rbs
        #print "\n\n"
        #print "gB(geom) at outboard",gB[ind0]
        #print "B(rbs) at x0",(rbs[x0_ind_rbs,25]**2+rbs[x0_ind_rbs,26]**2)**0.5
        error = abs(gB[ind0]-(rbs[x0_ind_rbs,25]**2+rbs[x0_ind_rbs,26]**2)**0.5)/abs(gB[ind0])
        if error>0.02:
            print( "Bfield error=",error)
            print( "Error=",error)
            dummy = raw_input("Press any key to continue")
        #print "error",error
        
        #print "\n\n"
        #print "R(geom) at outboard",gl_R[ind0]
        #print "R(rbs) at x0",rbs[x0_ind_rbs,24]
        error = abs(gl_R[ind0]-rbs[x0_ind_rbs,24])/abs(gl_R[ind0])
        if error>0.02:
            print( "R error=",error)
            dummy = raw_input("Press any key to continue")
        #print "error",error
        #print "\n\n"
        
        Rrbs = rbs[:,24]
        R_even = np.arange(1000)/999.0*(Rrbs[-1]-Rrbs[0])+Rrbs[0]
        rhot_R = interp(Rrbs,rbs[:,0],R_even)
        psi_R = interp(Rrbs,rbs[:,1],R_even)
        ind_x0_R_even = np.argmin(abs(rhot_R-x0))
        drhot_dR_R = fd_d1_o4(rhot_R,R_even)
        #print "drhot_dr (geom) at outboard", drhot_dr[ind0]
        #print "drhot_dr (rbs) at x0", drhot_dR_R[ind_x0_R_even]
        error = abs(drhot_dr[ind0]-drhot_dR_R[ind_x0_R_even])/drhot_dr[ind0]
        if error>0.02:
            print( "drhot_dr error=",error)
            dummy = raw_input("Press any key to continue")
        #print "error",error
        
        #print "\n\n"
        
        dpsi_dR_R = fd_d1_o4(psi_R,R_even)*psi0
        #print "dpsi_dr (geom) at outboard", dpsi_dr[ind0]
        #print "dpsi_dr (rbs) at x0", dpsi_dR_R[ind_x0_R_even]
        #print "Btheta*R (rbs) at x0", rbs[x0_ind_rbs,25]*rbs[x0_ind_rbs,24]
        error = abs(dpsi_dr[ind0]-dpsi_dR_R[ind_x0_R_even])/dpsi_dr[ind0]
        if error>0.02:
            print( "dpsi_dr error=",error)
            dummy = raw_input("Press any key to continue")
        #print "error",error
        
        #print "\n\n"
        
        #print "prefactor[ind0]",prefactor[ind0]
        prefactor_rbs = rbs[:,25]**2*rbs[:,24]**2/(rbs[:,25]**2+rbs[:,26]**2)**0.5
        #print "B_theta**2 R**2 / B", prefactor_rbs[x0_ind_rbs] 
        error = abs(prefactor[ind0]-prefactor_rbs[x0_ind_rbs])/abs(prefactor[ind0])
        if error>0.02:
            print( "prefactor error=",error)
            dummy = raw_input("Press any key to continue")
        #print "error",error
        
        #print "\n\n"
        
        #####Tests#####
        #####Tests#####
        #####Tests#####
    return  zgrid, dpsi_dr, prefactor
        
        
        
        
