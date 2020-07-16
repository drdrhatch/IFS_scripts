from pylab import *
from sys import argv,exit,stdout
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as US
from scipy import interpolate
import numpy as np
from finite_differences import *
from interp import *

def read_EFIT_file(efit_file_name):

    f = open(efit_file_name,'r')
    eqdsk=f.readlines()
    line1=eqdsk[0].split()
    nw=int(eqdsk[0].split()[-2])
    nh=int(eqdsk[0].split()[-1])
    print ('EFIT file Resolution: %d x %d' %(nw,nh))

    entrylength=16
    #note: here rmin is rleft from EFIT
    try:
        rdim,zdim,rctr,rmin,zmid=[float(eqdsk[1][j*entrylength:(j+1)*entrylength]) for j in range(int(len(eqdsk[1])/entrylength))]
    except:
        entrylength=15
        try:
            rdim,zdim,rctr,rmin,zmid=[float(eqdsk[1][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk[1])/entrylength)]
        except:
            exit('Error reading EQDSK file, please check format!')

    rmag,zmag,psiax,psisep,Bctr=[float(eqdsk[2][j*entrylength:(j+1)*entrylength]) for j in range(int(len(eqdsk[2])/entrylength))]
    dum,psiax2,dum,rmag2,dum=[float(eqdsk[3][j*entrylength:(j+1)*entrylength]) for j in range(int(len(eqdsk[3])/entrylength))]
    zmag2,dum,psisep2,dum,dum=[float(eqdsk[4][j*entrylength:(j+1)*entrylength]) for j in range(int(len(eqdsk[4])/entrylength))]
    if rmag!=rmag2: sys.exit('Inconsistent rmag: %7.4g, %7.4g' %(rmag,rmag2))
    if psiax2!=psiax: sys.exit('Inconsistent psiax: %7.4g, %7.4g' %(psiax,psiax2))
    if zmag!=zmag2: sys.exit('Inconsistent zmag: %7.4g, %7.4g' %(zmag,zmag2) )
    if psisep2!=psisep: sys.exit('Inconsistent psisep: %7.4g, %7.4g' %(psisep,psisep2))

    ###print "rmag", rmag
    ###print "zmag", zmag
    ###print "psiax", psiax
    ###print "psisep", psisep
    ###print "Bctr", Bctr
    # (R,Z) grid on which psi_pol is written
    Rgrid = np.arange(nw)/float(nw-1)*rdim+rmin
    ###print "rdim",rdim
    ###print "rmin",rmin
    ###print "first few Rgrid points", Rgrid[0:6]
    ###print "last few Rgrid points", Rgrid[-7:-1]
    Zgrid = np.arange(nh)/float(nh-1)*zdim+(zmid-zdim/2.0)
    ###print "zdim",zdim
    ###print "zmid",zmid
    ###print "first few Zgrid points", Zgrid[0:6]
    ###print "last few Zgrid points", Zgrid[-7:-1]

    # F, p, ffprime, pprime, q are written on uniform psi_pol grid
    # uniform grid of psi_pol~[psiax,psisep], resolution=nw
    F=empty(nw,dtype=float)
    p=empty(nw,dtype=float)
    ffprime=empty(nw,dtype=float)
    pprime=empty(nw,dtype=float)
    qpsi=empty(nw,dtype=float)
    # psi_pol is written on uniform (R,Z) grid (res=nw(R)*nh(Z))
    psirz_1d=empty(nw*nh,dtype=float)
    
    start_line=5
    lines=range(int(nw/5))
    if nw%5!=0: lines=range(int(nw/5)+1)
    for i in lines:
        n_entries=int(len(eqdsk[i+start_line])/entrylength)
        F[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
    start_line=i+start_line+1

    for i in lines:
        n_entries=int(len(eqdsk[i+start_line])/entrylength)
        p[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
    start_line=i+start_line+1

    for i in lines:
        n_entries=int(len(eqdsk[i+start_line])/entrylength)
        ffprime[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
    start_line=i+start_line+1

    for i in lines:
        n_entries=int(len(eqdsk[i+start_line])/entrylength)
        pprime[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
    start_line=i+start_line+1

    lines_twod=range(int(nw*nh/5))
    if nw*nh%5!=0: lines_twod=range(int(nw*nh/5)+1)
    for i in lines_twod:
        n_entries=int(len(eqdsk[i+start_line])/entrylength)
        psirz_1d[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
    start_line=i+start_line+1
    psirz=psirz_1d.reshape(nh,nw)

    for i in lines:
        n_entries=int(len(eqdsk[i+start_line])/entrylength)
        qpsi[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
    start_line=i+start_line+1

    # even grid of psi_pol, on which all 1D fields are defined
    psip_n = np.linspace(0.0,1.0,nw)
    # return data read from efit file 
    # psip_n: uniform flux grid from magnetic axis to separatrix
    # F, p, ffprime, pprime, qpsi are on psip_n
    # uniform (R,Z) grid, psirz is on this grid
    return psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw,psiax,psisep

def calc_rho_tor(psip_n, psiax, psisep, qpsi, nw,psip_n_max=0.999):

    #create rho_tor grid on even psi_pol grid
    interpol_order = 3 
    psi_pol = np.empty(len(psip_n))
    for i in range(len(psip_n)):
        psi_pol[i] = psiax+psip_n[i]*(psisep-psiax)

    q_spl_psi = US(psi_pol, qpsi, k=interpol_order, s=1e-5)
    psi_pol_fine = linspace(psi_pol[0], psi_pol[-1], nw*10)
    psi_tor_fine = empty((nw*10),dtype=float)
    psi_tor_fine[0] = 0.
    qfine = np.empty(nw*10)
    ##################
    ##################
    qnew = q_spl_psi(psi_pol_fine)
    psip_n_fine = (psi_pol_fine-psiax)/(psisep-psiax)
    #plt.plot(psip_n_fine,qnew)
    #plt.plot(psip_n,qpsi,'x-')
    #plt.show()
    ##################
    ##################
    for i in range(1,nw*10):
        x=psi_pol_fine[:i+1]
        y=q_spl_psi(x)
        psi_tor_fine[i]=np.trapz(y,x)

    phi_edge = psi_tor_fine[-1]
    rhot_n_fine=np.sqrt(psi_tor_fine/(psi_tor_fine[-1]-psi_tor_fine[0]))
    rho_tor_spl=US(psi_pol_fine, rhot_n_fine, k=interpol_order, s=1e-5)
    # rhot_n grid (not uniform, on even grid of psi_pol) of resolution=nw 
    rhot_n=rho_tor_spl(psi_pol)
    
    # rho_tor_spl takes value of psi_pol (not normalized) and convert into rhot_n
    return rho_tor_spl, rhot_n, phi_edge

def calc_B_fields(Rgrid, rmag, Zgrid, zmag, psirz, psiax, psisep, F, nw, psip_n):
    
    # Z0_ind is the index of Zgrid of midplane
    Z0_ind = np.argmin(np.abs(Zgrid-zmag))
    # psi_midplane is psi_pol at midplane on even Rgrid
    psi_pol_mp = psirz[Z0_ind,:]
    # Rmag_ind is index of unif_R at rmag
    Rmag_ind = np.argmin(np.abs(Rgrid - rmag))
    ###print "rmag",rmag
    ###print "Rmag_ind",Rmag_ind
    ###print "Rgrid[Rmag_ind]~rmag", Rgrid[Rmag_ind]
    ###print "psi_pol_mp[Rmag_ind]~psiax", psi_pol_mp[Rmag_ind]
    psi_pol_obmp = psi_pol_mp[Rmag_ind:].copy()
    #normalize psi_pol_obmp to psip_n_temp
    psip_n_temp = np.empty(len(psi_pol_obmp))
    for i in range(len(psi_pol_obmp)):
        psip_n_temp[i] = (psi_pol_obmp[i]-psiax)/(psisep-psiax)
    unif_R = np.linspace(Rgrid[Rmag_ind],Rgrid[-1],nw*10)
#    unif_R = np.linspace(Rgrid[Rmag_ind],Rgrid[-1],nw)
    psip_n_unifR = interp(Rgrid[Rmag_ind:],psip_n_temp,unif_R)
    psisep_ind = np.argmin(abs(psip_n_unifR-1.02))
    ###print "psisep_ind", psisep_ind
    ###print "psip_n_temp[psisep_ind]~1", psip_n_unifR[psisep_ind]
    #print "we have a problem here because uniform R grid doesn't have enough resolution near separatrix"
    psip_n_obmp = psip_n_unifR[:psisep_ind].copy()
    ###print "psip_n_obmp[0]~0", psip_n_obmp[0]
    ###print "psip_n_obmp[-1]~1", psip_n_obmp[-1]
    #plt.plot(psi_pol_obmp)
    #plt.show()
    R_obmp = unif_R[:psisep_ind].copy()
    # B_pol is d psi_pol/ d R * (1/R)
    #B_pol = fd_d1_o4(psi_pol_obmp, unif_R[Rmag_ind:Rmag_ind+psisep_ind])/unif_R[Rmag_ind:Rmag_ind+psisep_ind]
    B_pol = fd_d1_o4(psip_n_obmp*(psisep-psiax)+psiax,R_obmp)/R_obmp
    # convert F(on even psi_pol grid) to F(on even R grid)
    F_obmp = interp(psip_n, F, psip_n_obmp)
    # B_tor = F/R
    B_tor = F_obmp/R_obmp

    # psip_n_obmp is normalized psi_pol at outboard midplane on uniform unif_R
    # B_tor and B_pol are on uniform unif_R as well
    # psip_n_obmp is unlike psip_n ([0,1]), it goes from 0 to 1.06 here
    return psip_n_obmp, R_obmp, B_pol, B_tor

def read_EFIT_parameters(efit_file_name):

    f = open(efit_file_name,'r')
    eqdsk=f.readlines()
    line1=eqdsk[0].split()
    nw=int(eqdsk[0].split()[-2])
    nh=int(eqdsk[0].split()[-1])

    entrylength=16
    #note: here rmin is rleft from EFIT
    try:
        rdim,zdim,rctr,rmin,zmid=[float(eqdsk[1][j*entrylength:(j+1)*entrylength]) for j in range(int(len(eqdsk[1])/entrylength))]
    except:
        entrylength=15
        try:
            rdim,zdim,rctr,rmin,zmid=[float(eqdsk[1][j*entrylength:(j+1)*entrylength]) for j in range(int(len(eqdsk[1])/entrylength))]
        except:
            exit('Error reading EQDSK file, please check format!')

    rmag,zmag,psiax,psisep,Bctr=[float(eqdsk[2][j*entrylength:(j+1)*entrylength]) for j in range(int(len(eqdsk[2])/entrylength))]
    curr,psiax2,dum,rmag2,dum=[float(eqdsk[3][j*entrylength:(j+1)*entrylength]) for j in range(int(len(eqdsk[3])/entrylength))]
    zmag2,dum,psisep2,dum,dum=[float(eqdsk[4][j*entrylength:(j+1)*entrylength]) for j in range(int(len(eqdsk[4])/entrylength))] 
    if rmag!=rmag2: sys.exit('Inconsistent rmag: %7.4g, %7.4g' %(rmag,rmag2))
    if psiax2!=psiax: sys.exit('Inconsistent psiax: %7.4g, %7.4g' %(psiax,psiax2))
    if zmag!=zmag2: sys.exit('Inconsistent zmag: %7.4g, %7.4g' %(zmag,zmag2) )
    if psisep2!=psisep: sys.exit('Inconsistent psisep: %7.4g, %7.4g' %(psisep,psisep2))

    return rdim,zdim,rctr,rmin,zmid,Bctr,curr,nh

def get_dimpar_pars(efit_file_name,rhot0):

    psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw,psiax,psisep = read_EFIT_file(efit_file_name)
    #plt.plot(Rgrid)
    #plt.title('Rgrid (rmag = '+str(rmag)+')')
    #plt.show()
    R_major = rmag
    dummy, rhot_n, phi_edge = calc_rho_tor(psip_n, psiax, psisep, qpsi, nw,psip_n_max=0.999)
    psip_n_obmp, R_obmp, B_pol, B_tor = calc_B_fields(Rgrid, rmag, Zgrid, zmag, psirz, psiax, psisep, F, nw, psip_n)
    Bref = abs(B_tor[0])
    Lref = np.sqrt(2.0*abs(phi_edge/Bref))
    #plt.plot(rhot_n,psip_n)
    #plt.xlabel('rhot_n')
    #plt.ylabel('psi_n')
    #plt.show()
    irhot_n = np.argmin(abs(rhot_n-rhot0))
    q0 = qpsi[irhot_n]
    return Lref, Bref, R_major, q0

def get_current_density(efit_file_name):

    psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw,psiax,psisep = read_EFIT_file(efit_file_name)
    Jtot = Rgrid*pprime+ffprime/Rgrid
    return psip_n,Rgrid,Jtot


def get_geom_pars(efit_file_name,rhot0):

    psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw,psiax,psisep = read_EFIT_file(efit_file_name)
    R_major = rmag
    dummy, rhot_n, phi_edge = calc_rho_tor(psip_n, psiax, psisep, qpsi, nw,psip_n_max=0.999)
    psip_n_obmp, R_obmp, B_pol, B_tor = calc_B_fields(Rgrid, rmag, Zgrid, zmag, psirz, psiax, psisep, F, nw, psip_n)
    Bref = abs(B_tor[0])
    Lref = np.sqrt(2.0*abs(phi_edge/Bref))
    irhot_n = np.argmin(abs(rhot_n-rhot0))
    q0 = qpsi[irhot_n]


    rhot_new = np.linspace(rhot_n[0],rhot_n[-1],4.0*len(rhot_n))
    qpsi_new = interp(rhot_n,qpsi,rhot_new)  
    shat = rhot_new/qpsi_new*fd_d1_o4(qpsi_new,rhot_new)
    irhot_new = np.argmin(abs(rhot_new-rhot0))
    shat0 = shat[irhot_new]
    return Lref, Bref, R_major, q0, shat0



