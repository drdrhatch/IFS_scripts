#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plots: q(psi)
       P(psi)
       psi(R,Z)
       psi(R,Z) for psi = 0.9 and 1.0 along with grid points
-c: outputs rho_tor and rho_pol
-p: outputs R, psi, B_pol
-n: suppresses plots

Modified from extract_miller_from_eqdsk.py
"""

from pylab import *
from sys import argv,exit,stdout
import optparse as op
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as US
from scipy import interpolate

parser=op.OptionParser(description='Tools for extracting information from an EFIT file.')
#parser.add_option('--rovera','-r',action='store_const',const=1)
parser.add_option('--conv','-c',action='store_const',const=1,help = 'Output rho_tor vs rho_pol')
parser.add_option('--binfo','-p',action='store_const',const=1,help = 'Ouptut R, psi_pol, B_pol, B_tor')
parser.add_option('--noplot','-n',action='store_const',const=1,help = 'Suppress plots')
options,args=parser.parse_args()
write_rhotor_rhopol_conversion=options.conv
write_binfo=options.binfo
show_plots = not options.noplot
#use_r_a=options.rovera
if len(args)!=1:
    exit("""
Please include efit file name as argument."
    \n""")

filename=args[0]
file=open(filename,'r')

def find(val,arr):
    return argmin(abs(arr-val))

def fd_d1_o4(var,grid,mat=False):
    """Centered finite difference, first derivative, 4th order.
    var: quantity to be differentiated.
    grid: grid for var
    mat: matrix for the finite-differencing operator. if mat=False then it is created"""

    if not mat:
        mat=get_mat_fd_d1_o4(len(var),grid[1]-grid[0])

    dvar=np.dot(mat,var)
    dvar[0]=0.0
    dvar[1]=0.0
    #dvar[2]=0.0
    dvar[-1]=0.0
    dvar[-2]=0.0
    #dvar[-3]=0.0
    return dvar

def get_mat_fd_d1_o4(size,dx,plot_matrix=False):
    """Creates matrix for centered finite difference, first derivative, 4th order.
    size: size of (number of elements in) quantity to be differentiated
    dx: grid spacing (for constant grid)."""

    prefactor=1.0/(12.0*dx)
    mat=np.zeros((size,size),dtype='float')
    for i in range(size):
        if i-1 >= 0:
            mat[i,i-1]=-8
        if i-2 >= 0:
            mat[i,i-2]=1
        if i+1 <= size-1:
            mat[i,i+1]=8
        if i+2 <= size-1:
            mat[i,i+2]=-1

    mat=prefactor*mat

    if plot_matrix:
        plt.contourf(mat,50)
        plt.colorbar()
        plt.show()

    return mat

def interp(xin,yin,xnew):
    """
    xin: x variable input
    yin: y variable input
    xnew: new x grid on which to interpolate
    yout: new y interpolated on xnew
    """

    #splrep returns a knots and coefficients for cubic spline
    rho_tck = interpolate.splrep(xin,yin)
    #Use these knots and coefficients to get new y
    yout = interpolate.splev(xnew,rho_tck,der=0)

    return yout

def full_interp(func_xin,xin,xconv,yconv,yout):
    """
    Takes function func_xin on grid xin and outputs the function on yout grid
    func_xin: function to interpolate
    xin: grid corresponding to func_xin
    xconv: xgrid for conversion
    yconv: ygrid for conversion
    yout: output grid
    """

    #If necessary truncate func_xin onto correct range
    if xin[0] < xconv[0]:
        low_index = np.argmin(abs(xconv-xin[0]))
    else:
        low_index = 0
    if xin[-1] > xconv[-1]:
        high_index = np.argmin(abs(xconv-xin[-1]))
    else:
        high_index = -1

    if high_index == -1:
        func_xin = func_xin[low_index:]
    else:
        func_xin = func_xin[low_index:high_index]

    func_xconv = interp(xin,func_xin,xconv)
    func_yout = interp(yconv,func_xconv,yout)

    return func_yout

eqdsk=file.readlines()
print 'Header: %s' %eqdsk[0]
#set resolutions
nw=int(eqdsk[0].split()[-2]);nh=int(eqdsk[0].split()[-1])
pw=(nw/8/2)*2 #psi-width, number of flux surfaces around position of interest
print 'Resolution: %d x %d' %(nw,nh)

entrylength=16
#note: here rmin is rleft from EFIT
try:
    rdim,zdim,rctr,rmin,zmid=[float(eqdsk[1][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk[1])/entrylength)]
except:
    entrylength=15
    try:
        rdim,zdim,rctr,rmin,zmid=[float(eqdsk[1][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk[1])/entrylength)]
    except:
        exit('Error reading EQDSK file, please check format!')

rmag,zmag,psiax,psisep,Bctr=[float(eqdsk[2][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk[2])/entrylength)]
dum,psiax2,dum,rmag2,dum=[float(eqdsk[3][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk[3])/entrylength)]
zmag2,dum,psisep2,dum,dum=[float(eqdsk[4][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk[4])/entrylength)]
if rmag!=rmag2: sys.exit('Inconsistent rmag: %7.4g, %7.4g' %(rmag,rmag2))
if psiax2!=psiax: sys.exit('Inconsistent psiax: %7.4g, %7.4g' %(psiax,psiax2))
if zmag!=zmag2: sys.exit('Inconsistent zmag: %7.4g, %7.4g' %(zmag,zmag2) )
if psisep2!=psisep: sys.exit('Inconsistent psisep: %7.4g, %7.4g' %(psisep,psisep2))

print "rmag", rmag
print "zmag", zmag
print "psiax", psiax
print "psisep", psisep
print "Bctr", Bctr
Rgrid = np.arange(nw)/float(nw-1)*rdim+rmin
print "rdim",rdim
print "rmin",rmin
#print "Rgrid",Rgrid
Zgrid = np.arange(nh)/float(nh-1)*zdim+(zmid-zdim/2.0)
print "zdim",zdim
print "zmid",zmid
#print "Zgrid",Zgrid

F=empty(nw,dtype=float)
p=empty(nw,dtype=float)
ffprime=empty(nw,dtype=float)
pprime=empty(nw,dtype=float)
qpsi=empty(nw,dtype=float)
psirz_1d=empty(nw*nh,dtype=float)
start_line=5
lines=range(nw/5)
if nw%5!=0: lines=range(nw/5+1)
for i in lines:
    n_entries=len(eqdsk[i+start_line])/entrylength
    F[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1

for i in lines:
    n_entries=len(eqdsk[i+start_line])/entrylength
    p[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1

for i in lines:
    n_entries=len(eqdsk[i+start_line])/entrylength
    ffprime[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1

for i in lines:
    n_entries=len(eqdsk[i+start_line])/entrylength
    pprime[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1

lines_twod=range(nw*nh/5)
if nw*nh%5!=0: lines_twod=range(nw*nh/5+1)
for i in lines_twod:
    n_entries=len(eqdsk[i+start_line])/entrylength
    psirz_1d[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1
psirz=psirz_1d.reshape(nh,nw)

for i in lines:
    n_entries=len(eqdsk[i+start_line])/entrylength
    qpsi[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1

#linear grid of psi, on which all 1D fields are defined
linpsi=linspace(psiax,psisep,nw)
#create rho_tor grid
x_fine=linspace(psiax,psisep,nw*10)
phi_fine=empty((nw*10),dtype=float)
phi_fine[0]=0.
#spline of q for rhotor grid
interpol_order=3
q_spl_psi=US(linpsi,qpsi,k=interpol_order,s=1e-5)

for i in range(1,nw*10):
    x=x_fine[:i+1]
    y=q_spl_psi(x)
    phi_fine[i]=trapz(y,x)
rho_tor_fine=sqrt(phi_fine/phi_fine[-1])
rho_tor_spl=US(x_fine,rho_tor_fine,k=interpol_order,s=1e-5)
rho_tor=empty(nw,dtype=float)
for i in range(nw):
    rho_tor[i]=rho_tor_spl(linpsi[i])

if write_rhotor_rhopol_conversion:
    rt_rp_filename='rt_rp_%s' %filename
    rt_rp_file=open(rt_rp_filename,'w')
    rt_rp_file.write('# rho_tor          rho_pol\n')
    for i in range(len(x_fine)):
        rho_pol=sqrt((x_fine[i]-psiax)/(psisep-psiax))
        rt_rp_file.write('%16.8e %16.8e\n' %(rho_tor_fine[i],rho_pol))
    rt_rp_file.close()
    exit('\nWrote rhotor/rhopol conversion to %s.' %rt_rp_filename)



Z0_ind = np.argmin(np.abs(Zgrid-zmid))
psi=np.arange(nw)/float(nw-1)*(psisep-psiax)
#print "psi",psi
psi_midplane = psirz[Z0_ind,:]

#plt.plot(psi_midplane)
#plt.title(r'$\Psi$'+' midplane')
#plt.show()

#plt.plot(psi_midplane-psiax)
#plt.title(r'$\Psi$'+' midplane(adjusted)')
#plt.show()

if show_plots:
    plt.plot(psi/(psisep-psiax),qpsi,'-')
    plt.xlabel(r'$\Psi/\Psi_0$')
    plt.title(r'$q$')
    plt.show()

    plt.plot(psi/(psisep-psiax),p)
    plt.xlabel(r'$\Psi/\Psi_0$')
    plt.title(r'$P(Nm^{-2})$')
    plt.show()

#plt.plot(psi,F)
#plt.xlabel(r'$\Psi$')
#plt.title(r'$RB_\phi(mT)$')
#plt.show()

    plt.contour(Rgrid,Zgrid,(psirz[:]-psiax)/(psisep-psiax),50)
    plt.colorbar()
    plt.contour(Rgrid,Zgrid,(psirz[:]-psiax)/(psisep-psiax),levels = [1.0],color = 'black',linewidth = 2)
    plt.contour(Rgrid,Zgrid,(psirz[:]-psiax)/(psisep-psiax),levels = [0.9],color = 'black',linewidth = 2)
    plt.xlabel('R(m)')
    plt.ylabel('Z(m)')
    plt.show()

#plt.contour(Rgrid,Zgrid,(psirz[:]-psiax)/(psisep-psiax),50)
#plt.contour(Rgrid,Zgrid,(psirz[:]-psiax)/(psisep-psiax),levels = [1.0],color = 'black', linewidth = 2)
    plt.contour(Rgrid,Zgrid,(psirz[:]-psiax)/(psisep-psiax),levels = [0.9,1.0],color = 'black', linewidth = 4)
    plt.colorbar()
    Zplot = np.zeros(nw)
    for j in range(len(Zgrid)):
        Zplot[:] = Zgrid[j]
        plt.plot(Rgrid,Zplot,'.',color = 'black', markersize = 1)
    plt.xlabel('R(m)')
    plt.ylabel('Z(m)')
    plt.show()

if write_binfo:
    #Calculate midplane B_pol
    B_pol = fd_d1_o4(psi_midplane,Rgrid)/Rgrid
    psi_norm_out = (psi_midplane-psiax)/(psisep-psiax)
    F_out = interp(psi/(psisep-psiax),F,psi_norm_out)
    q_out = interp(psi/(psisep-psiax),qpsi,psi_norm_out)
    f=open('Binfo_'+filename,'w')
    f.write('# Outer Midplane')
    f.write('# 1.R(m) 2.psi_norm 3.B_pol(T) 4.B_tor(T) 5.q\n')
    f.write('# R at magnetic axis = '+str(rmag)+'\n')
    f.write('# psisep - psiax = '+str(psisep-psiax)+'\n')
    Rmag_ind = np.argmin(abs(Rgrid - rmag))
    print "rmag",rmag
    print "Rmag_ind",Rmag_ind
    print "Rgrid[Rmag_ind]",Rgrid[Rmag_ind]
    temp = psi_norm_out
    temp[0:Rmag_ind] = 0
    psi_ind_sep = np.argmin(abs(temp-1.05))
    print "psi_ind_sep",psi_ind_sep
    B_tor = F_out / Rgrid
    np.savetxt(f,np.column_stack((Rgrid[Rmag_ind:psi_ind_sep],psi_norm_out[Rmag_ind:psi_ind_sep],B_pol[Rmag_ind:psi_ind_sep],B_tor[Rmag_ind:psi_ind_sep],q_out[Rmag_ind:psi_ind_sep])))
    f.close()
    #plt.plot(psi_norm_out[Rmag_ind:psi_ind_sep],q_out[Rmag_ind:psi_ind_sep])
    #plt.plot(psi/(psisep-psiax),qpsi)
    #plt.plot(psi_norm_out[Rmag_ind:psi_ind_sep],(Rgrid[Rmag_ind:psi_ind_sep]-rmag)*B_tor[Rmag_ind:psi_ind_sep]/(Rgrid[Rmag_ind:psi_ind_sep]*B_pol[Rmag_ind:psi_ind_sep]))
    #plt.show()

    if show_plots:
        plt.plot(Rgrid[Rmag_ind:psi_ind_sep],B_pol[Rmag_ind:psi_ind_sep])
        plt.title(r'$B_\theta$')
        plt.xlabel('R(m)')
        plt.show()

        plt.plot(psi_norm_out[Rmag_ind:psi_ind_sep],B_pol[Rmag_ind:psi_ind_sep])
        plt.title(r'$B_\theta$')
        plt.xlabel(r'$\Psi$')
        plt.show()
