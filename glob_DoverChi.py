#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from finite_differences import *
from interp import *
import optparse as op
from ParIO import * 
from parIOWrapper import *

def readProfiles(speciesName, suffix, subtract_convection = True):

    profilesFilename = 'profiles_'+speciesName+'_'+ suffix
    fluxprofFilename = 'fluxprof'+speciesName+'_'+ suffix+'.dat'

    profiles = np.genfromtxt(profilesFilename)
    rho = profiles[:,0]
    tempProf = profiles[:,2]
    densProf = profiles[:,3]
    omtProf = profiles[:,4]
    omnProf = profiles[:,5]

    fluxprof = np.genfromtxt(fluxprofFilename,skip_footer=3)
    fluxrho = fluxprof[:,0]
    Gamma_es = fluxprof[:,1]
    Qheat_es = fluxprof[:,2]
    Gamma_em = fluxprof[:,3]
    Qheat_em = fluxprof[:,4]
    Gamma_tot = Gamma_es + Gamma_em
    Qheat_tot = Qheat_es + Qheat_em

    suffix = '_'+suffix
    #pars = init_read_parameters_file(suffix)
    par = Parameters()
    par.Read_Pars('parameters'+suffix)
    pars = par.pardict
    tempProf = tempProf / pars['Tref']
    densProf = densProf / pars['nref']

    temp_fluxrho = interp(rho, tempProf, fluxrho)
    dens_fluxrho = interp(rho, densProf, fluxrho)
    omt_fluxrho = interp(rho, omtProf, fluxrho)
    omn_fluxrho = interp(rho, omnProf, fluxrho)

    if subtract_convection:
        Qheat_tot = Qheat_tot - 3./2.*Gamma_tot*temp_fluxrho

    if 1 == 0:
        plt.plot(fluxrho, Gamma_tot, label='Gamma')
        plt.plot(fluxrho, Qheat_tot, label='Q')
        plt.legend()
        plt.title(speciesName)
        plt.show()

    return rho, temp_fluxrho, dens_fluxrho, omt_fluxrho, omn_fluxrho, fluxrho, Gamma_tot, Qheat_tot

def weightedAvg(fluxrho, Gamma, Qheat, omt, omn, t, n):

    weight = np.abs(Qheat)
    norm = np.sum(weight)

    avgGamma = np.sum(Gamma * weight) / norm
    avgQheat = np.sum(Qheat * weight) / norm

    avgdndrho = np.sum(omn * n * weight) / norm
    avgndtdrho = np.sum(omt * t * n * weight) / norm

    D = avgGamma / avgdndrho
    chi = avgQheat / avgndtdrho

    avgt = np.sum(t * weight) / norm
    avgn = np.sum(n * weight) / norm

    return avgGamma, avgQheat, D, chi, avgt, avgn

def weightedScaleLength(fluxrho, Qheat, omt, omn):

    weight = np.abs(Qheat)
    norm = np.sum(weight)

    if 1 == 0:
        plt.plot(fluxrho, omt, '.', label='omt')
        plt.plot(fluxrho, omn, '.', label='omn')
        plt.legend()
        plt.show()

    avg_omt = np.sum(omt * weight) / norm
    avg_omn = np.sum(omn * weight) / norm

    return avg_omt, avg_omn

parser = op.OptionParser()
options, args = parser.parse_args()
suffix = args[0]

subtract_convection = True

species = ['i','e']
Gamma = np.empty(len(species),dtype='float')
Qheat = np.empty(len(species),dtype='float')
D = np.empty(len(species),dtype='float')
chi = np.empty(len(species),dtype='float')
omt = np.empty(len(species),dtype='float')
omn = np.empty(len(species),dtype='float')
t = np.empty(len(species),dtype='float')
n = np.empty(len(species),dtype='float')

for i in range(len(species)):
    s = species[i]
    rho, t_fluxrho, n_fluxrho, omt_fluxrho, omn_fluxrho, fluxrho, Gamma_fluxrho, Qheat_fluxrho = readProfiles(s, suffix, subtract_convection)
    Gamma[i], Qheat[i], D[i], chi[i], t[i], n[i] = weightedAvg(fluxrho, Gamma_fluxrho, Qheat_fluxrho, omt_fluxrho, omn_fluxrho, t_fluxrho, n_fluxrho)
    outstr = 'D_'+s+'/ chi_'+s+' ='
    print outstr, D[i] / chi[i]
    print 'Gamma =', Gamma[i]
    print 'Qheat =', Qheat[i]
    print 'Gamma / Qheat =', Gamma[i] / Qheat[i]

    omt[i], omn[i] = weightedScaleLength(fluxrho, Qheat_fluxrho, omt_fluxrho, omn_fluxrho)
    outstr = 'omt_'+s+', omn_'+s+' ='
    print outstr, omt[i], omn[i]
    

chi_tot = chi[0] + chi[1]
print 'D_i / chi_tot =', D[0] / chi_tot
print 'D_e / chi_tot =', D[1] / chi_tot
#print 'D_z / chi_tot =', D[2] / chi_tot

print 'chi_i / chi_e =', chi[0]/chi[1]
#print 'D_z / chi_e =', D[2]/chi[1]

Qheat_tot = Qheat[0] + Qheat[1]
print Gamma
print Qheat_tot
print 'Gamma_i / Q_tot =', Gamma[0] / Qheat_tot
print 'Gamma_e / Q_tot =', Gamma[1] / Qheat_tot
#print 'Gamma_z / Q_tot =', Gamma[2] / Qheat_tot * n[1] / n[2]

if 1 == 0:
    s = 'i'
    rho, t, n, fluxrho, Gamma, Qheat = readProfiles(s, suffix)
    dndrho, dtdrho, t_fluxrho, n_fluxrho = gradProf(rho, t, n, fluxrho)
    suffix = '_'+suffix
    pars = init_read_parameters_file(suffix)
    Bref, Tref, nref, Lref, mref = read_ref_values(suffix, pars)
    a = pars['Lref'] * fluxrho
    R = a * pars['major_R']
    Area = 2. * np.pi * a * 2. * np.pi * R
    n0 = n_fluxrho * 1.E19
    t0 = t_fluxrho * 1E3 * 1.6E-19
    D0 = 0.035
    chi0i = 0.15
    Gamma0 = Area * D0 * n0 * dndrho / n_fluxrho / a
    Qheat0i = Area * chi0i * n0 * t0 * dtdrho / t_fluxrho / a
    if 1 == 1:
        plt.plot(fluxrho, Gamma0)
        plt.xlabel('rho')
        plt.ylabel('Gamma (particle # / second)')
        plt.show() 
    if 1 == 1:
        plt.plot(fluxrho, Qheat0i)
        plt.xlabel('rho')
        plt.ylabel('Qheati (W)')
        plt.show() 
    if 1 == 1:
        plt.plot(fluxrho, Gamma0*1.6E-19/1E3,label = 'Gamma (MW/keV)')
        plt.plot(fluxrho, Qheat0i*1.E-6, label = 'Qheat (MW)')
        plt.xlabel('rho')
        plt.legend()
        plt.show() 
