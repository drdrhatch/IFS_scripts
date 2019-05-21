#! /usr/bin/python

from pylab import *
from sys import argv,exit,stdout
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy import interpolate
import numpy as np
from finite_differences_x import *
from interp import *

def read_EFIT(EFIT_file_name):
    EFITdict = {}
    f = open(EFIT_file_name,'r')
    eqdsk=f.readlines()
    line1=eqdsk[0].split()
    if len(line1) == 2:
        nwnh = eqdsk[0].split()[1]
        nw = int(nwnh[0:3])
        nh = int(nwnh[3:7])
    else:
        nw = int(line1[-2])
        nh = int(line1[-1])

    EFITdict['nw'] = nw    #number of grid for Rgrid
    EFITdict['nh'] = nh    #number of grid for Zgrid

    entrylength=16
    #note: here rmin is rleft from EFIT
    #print(len(eqdsk[1])/entrylength) is not integer
    try:
        rdim,zdim,rcentr,rleft,zmid = \
            [float(eqdsk[1][j*entrylength:(j+1)*entrylength]) \
            for j in range(len(eqdsk[1])//entrylength)]
    except:
        entrylength=15
        try:
            rdim,zdim,rcentr,rleft,zmid = \
                [float(eqdsk[1][j*entrylength:(j+1)*entrylength]) \
                for j in range(len(eqdsk[1])//entrylength)]
        except:
            exit('Error reading EQDSK file, please check format!')

    rmaxis,zmaxis,simag,sibry,bcentr = \
        [float(eqdsk[2][j*entrylength:(j+1)*entrylength]) \
        for j in range(len(eqdsk[2])//entrylength)]
    current,simag2,xdum,rmaxis2,xdum = \
        [float(eqdsk[3][j*entrylength:(j+1)*entrylength]) \
        for j in range(len(eqdsk[3])//entrylength)]
    zmaxis2,xdum,sibry2,xdum,xdum = \
        [float(eqdsk[4][j*entrylength:(j+1)*entrylength]) \
        for j in range(len(eqdsk[4])//entrylength)]

    EFITdict['rdim'] = rdim
    EFITdict['zdim'] = zdim
    EFITdict['rcentr'] = rcentr
    EFITdict['rleft'] = rleft
    EFITdict['zmid'] = zmid
    EFITdict['rmaxis'] = rmaxis    # R of magnetic axis (m)
    EFITdict['zmaxis'] = zmaxis    # Z of magnetic axis (m)
    EFITdict['simag'] = simag    # poloidal flux at magnetic axis
    EFITdict['sibry'] = sibry    # poloidal flux at plasma boundary
    EFITdict['bcentr'] = bcentr    # vacuum toroidal magnetic field in Telsa
    EFITdict['current'] = current    # plasma current in Ampere

    print('EFIT file Resolution: %d x %d' %(EFITdict['nw'],EFITdict['nh']))
    print('Horizontal dimension(m): %10.4f' %EFITdict['rdim'])
    print('Vertical dimension(m): %10.4f' %EFITdict['zdim'])
    print('Minimum R of rectangular grid: %10.4f' %EFITdict['rleft'])
    print('(R, Z) of magnetic axis: (%10.4f, %10.4f)' %(EFITdict['rmaxis'],EFITdict['zmaxis']))
    print('poloidal flux at magnetic axis in Weber/rad: %10.4f' %EFITdict['simag'])
    print('poloidal flux at the plasma boundary in Weber/rad: %10.4f' %EFITdict['sibry'])
    print('Vacuum toroidal magnetic field at R = %10.4f: %10.4f Tesla' %(EFITdict['rcentr'],EFITdict['bcentr']))
    print('Z of center of rectangular grid: %10.4f' %EFITdict['zmid'])
    print('Plasma current: %10.4f Ampere' %EFITdict['current'])

    Rgrid = np.linspace(0, 1, nw, endpoint = True) * rdim + rleft
    Zgrid = np.linspace(0, 1, nh, endpoint = True) * zdim + (zmid - zdim/2.)
    EFITdict['Rgrid'] = Rgrid    # Rgrid of psi(Z, R)
    EFITdict['Zgrid'] = Zgrid    # Zgrid of psi(Z, R)

    Fpol = empty(nw, dtype = float)
    Pres = empty(nw, dtype = float)
    FFprime = empty(nw, dtype = float)
    Pprime = empty(nw, dtype = float)
    qpsi = empty(nw, dtype = float)
    jtor = empty(nw, dtype = float)

    # psi_pol is written on uniform (R,Z) grid (res=nw(R)*nh(Z))
    psirz_1d = empty(nw * nh, dtype = float)
    
    start_line = 5
    wordsInLine = 5
    lines=range(nw//wordsInLine)
    if nw%wordsInLine!=0: lines=range(nw//wordsInLine+1)
    for i in lines:
        n_entries=len(eqdsk[i+start_line])//entrylength
        Fpol[i*wordsInLine:i*wordsInLine+n_entries] = \
            [float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) \
            for j in range(n_entries)]
    start_line=i+start_line+1
    EFITdict['Fpol'] = Fpol    # poloidal current function F = R * Btor on psipn grid

    for i in lines:
        n_entries=len(eqdsk[i+start_line])//entrylength
        Pres[i*wordsInLine:i*wordsInLine+n_entries] = \
            [float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) \
            for j in range(n_entries)]
    start_line=i+start_line+1
    EFITdict['Pres'] = Pres    # plasma pressure in N / m^2 on psipn grid

    for i in lines:
        n_entries=len(eqdsk[i+start_line])//entrylength
        FFprime[i*wordsInLine:i*wordsInLine+n_entries] = \
            [float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) \
            for j in range(n_entries)]
    start_line=i+start_line+1
    EFITdict['FFprime'] = FFprime    # FFprime on psipn grid

    for i in lines:
        n_entries=len(eqdsk[i+start_line])//entrylength
        Pprime[i*wordsInLine:i*wordsInLine+n_entries] = \
            [float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) \
            for j in range(n_entries)]
    start_line=i+start_line+1
    EFITdict['Pprime'] = Pprime    # Pprime on psipn grid

    lines_twod=range(nw*nh//wordsInLine)
    if nw*nh%wordsInLine!=0: lines_twod=range(nw*nh//wordsInLine+1)
    for i in lines_twod:
        n_entries=len(eqdsk[i+start_line])//entrylength
        psirz_1d[i*wordsInLine:i*wordsInLine+n_entries] = \
            [float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) \
            for j in range(n_entries)]
    start_line=i+start_line+1

    psirz=psirz_1d.reshape(nh,nw)
    EFITdict['psirz'] = psirz    # poloidal flux on the rectangular grid (Rgrid, Zgrid)

    for i in lines:
        n_entries=len(eqdsk[i+start_line])//entrylength
        qpsi[i*wordsInLine:i*wordsInLine+n_entries] = \
            [float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) \
            for j in range(n_entries)]
    start_line=i+start_line+1
    EFITdict['qpsi'] = qpsi    # safety factor on psipn grid

    # even grid of psi_pol, on which all 1D fields are defined
    psipn = np.linspace(0., 1., nw)
    EFITdict['psipn'] = psipn    # uniform psipn grid

    interpol_order = 3 
    psip = psipn * (sibry - simag) + simag
    q_spl_psi = UnivariateSpline(psip, qpsi, k=interpol_order, s=1e-5)
    psi_pol_fine = linspace(psip[0], psip[-1], nw*10)
    psi_tor_fine = empty((nw*10),dtype=float)
    psi_tor_fine[0] = 0.
    for i in range(1, nw * 10):
        x = psi_pol_fine[:i+1]
        y = q_spl_psi(x)
        psi_tor_fine[i] = np.trapz(y,x)

    rhot_n_fine = np.sqrt(psi_tor_fine/(psi_tor_fine[-1]-psi_tor_fine[0]))
    rho_tor_spl = UnivariateSpline(psi_pol_fine, rhot_n_fine, k=interpol_order, s=1e-5)
    rhotn = rho_tor_spl(psip)
    EFITdict['rhotn'] = rhotn    # square root of toroidal flux on psipn grid

    Z0_ind = np.argmin(abs(Zgrid - zmaxis))
    R0_ind = np.argmin(abs(Rgrid - rmaxis - 0.02))
    R_obmp = Rgrid[R0_ind:]
    psirz_obmp = psirz[Z0_ind, R0_ind:]
    psipn_obmp = (psirz_obmp - simag) / (sibry - simag)

    sepInd = np.argmin(abs(psipn_obmp - 1.))
    psipn_obmp = psipn_obmp[:sepInd + 1]
    R_obmp = list(R_obmp[:sepInd + 1])

    R = interp(psipn_obmp, R_obmp, psipn)
    if 1 == 0:
        plt.plot(psipn_obmp, R_obmp, label = 'before')
        plt.plot(psipn, R, label = 'after')
        plt.xlabel('psipn')
        plt.ylabel('R')
        plt.legend(loc = 2)
        plt.show()

    EFITdict['R'] = R    # major radius (m) on psipn grid

    #jtor = rmaxis * Pprime + FFprime / rmaxis
    jtor = R * Pprime + FFprime / R
    EFITdict['jtor'] = jtor    # toroidal current density on psipn grid

    #psirz_spl = interpolate.RectBivariateSpline(Zgrid, Rgrid, psirz)

    Bp_Z_grid = np.empty(np.shape(psirz))
    for i in range(nh):
        Bp_Z_grid_this = first_derivative(psirz[i,:].flatten(), Rgrid) / Rgrid
        Bp_Z_grid[i,:] = Bp_Z_grid_this.copy()

    Bp_R_grid = np.empty(np.shape(psirz))
    for i in range(nw):
        Bp_R_grid_this = - first_derivative(psirz[:,i].flatten(), Zgrid) / Rgrid[i]
        Bp_R_grid[:,i] = Bp_R_grid_this.copy()

    #Bp_R_spl = interpolate.RectBivariateSpline(Zgrid, Rgrid, Bp_R_grid)
    #Bp_Z_spl = interpolate.RectBivariateSpline(Zgrid, Rgrid, Bp_Z_grid)

    Bp_tot_grid = np.sqrt(Bp_R_grid ** 2 + Bp_Z_grid ** 2)
    Bp_obmp = Bp_tot_grid[Z0_ind, R0_ind : R0_ind + sepInd + 1]

    Bpol = interp(psipn_obmp, Bp_obmp, psipn)
    EFITdict['Bpol'] = Bpol    # B_pol on psipn grid

    F_spl = interpolate.UnivariateSpline(psipn, Fpol)
    Btor = F_spl(psipn) / R
    EFITdict['Btor'] = abs(Btor)    #B_tor on psipn grid

    return EFITdict

def magneticShear(EFITdict, show_plots = False):

    rhotn = EFITdict['rhotn']
    q = EFITdict['qpsi']

    #uni_rhot = np.linspace(rhotn[0], rhotn[-1], len(rhotn) * 10)
    uni_rhot = np.linspace(rhotn[0], rhotn[-1], len(rhotn))
    q_unirhot = interp(rhotn, q, uni_rhot)
    shat_unirhot = uni_rhot / q_unirhot * first_derivative(q_unirhot, uni_rhot)
    shat = interp(uni_rhot, shat_unirhot, rhotn)

    R_unirhot = interp(rhotn, EFITdict['R'], uni_rhot)
    Ls_unirhot = q_unirhot * R_unirhot / shat_unirhot
    Ls = interp(uni_rhot, Ls_unirhot, rhotn)

    if show_plots:
        plt.plot(uni_rhot, shat_unirhot)
        plt.ylabel('shat')
        plt.xlabel('rhot')
        plt.axis([0.8, 1., 0., 10.])
        plt.show()
        plt.plot(uni_rhot, Ls_unirhot)
        plt.ylabel('Ls')
        plt.xlabel('rhot')
        plt.axis([0.8, 1., 0., 2.])
        plt.show()

    return uni_rhot, shat_unirhot, Ls_unirhot
