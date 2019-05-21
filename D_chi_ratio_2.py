#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import optparse as op
from parIOWrapper import *
from nrgWrapper import *

parser = op.OptionParser()
options, args = parser.parse_args()
suffix = args[0]
if suffix != '.dat':
    suffix = '_' + suffix

pars = init_read_parameters_file(suffix)

def D_over_chi(time,nrg,omn,omt,T,n):
    Gamma_es, Gamma_em, Q_es, Q_em = \
                        read_Gamma_Q(time,nrg,False)
    Gamma = Gamma_es + Gamma_em
    Qtot = Q_es + Q_em
    #Qtot = Qtot - 5./3.*Gamma
    Qtot = Qtot - 3./2.*Gamma*T
    Qes = Q_es - 3./2.*Gamma_es*T
    Qem = Q_em - 3./2.*Gamma_em*T
    D = Gamma / omn / n
    chi = Qtot / omt / n / T
    Dochi = Gamma/Qtot*omt/omn*T

    return Dochi, Qtot, Qes, Qem, Gamma, D, chi
    
def D_over_chi_tot(Gamma_s,omn_s,n_s,Q_e,omt_e,Q_i,omt_i,Ti,ni):
    Ds_chi_tot = Gamma_s/omn_s/n_s/(Q_e/omt_e+Q_i/omt_i/Ti/ni)
    return Ds_chi_tot

if pars['n_spec'] == 1:
    q_e = -1.
    omn_e, omt_e = read_species_gradients(q_e,pars)
    temp_e, dens_e = read_species_tempdens(q_e,pars)
    time,nrge = read_from_nrg_files(pars,suffix,False)
    Dochi, Q, Qes, Qem, Gamma, D, chi = \
        D_over_chi(time,nrge,omn_e,omt_e,temp_e,dens_e)
    print 'electron: Q_em/Q_es = %12.5f' % float(Qem/Qes)
elif pars['n_spec'] == 2:
    q_i = 1.
    omn_i, omt_i = read_species_gradients(q_i,pars)
    temp_i, dens_i = read_species_tempdens(q_i,pars)
    q_e = -1.
    omn_e, omt_e = read_species_gradients(q_e,pars)
    temp_e, dens_e = read_species_tempdens(q_e,pars)

    time, nrgi, nrge = read_from_nrg_files(pars,suffix,False)
    Dochi_i, Q_i, Qes_i, Qem_i, Gamma_i, D_i, chi_i = \
        D_over_chi(time,nrgi,omn_i,omt_i,temp_i,dens_i)
    Dochi_e, Q_e, Qes_e, Qem_e, Gamma_e, D_e, chi_e = \
        D_over_chi(time,nrge,omn_e,omt_e,temp_e,dens_e)
    print 'electron: Q_em/Q_es = %12.5f' % float(Qem_e/Qes_e)
    print 'Q_i/Q_e = %12.5f' % float(Q_i/Q_e)
    print 'chi_i/chi_e =', chi_i / chi_e
    Dsochi_tot_i = D_over_chi_tot(Gamma_i,omn_i,dens_i,Q_e,omt_e,Q_i,omt_i,temp_i, dens_i)
    print 'D_i/chi_tot = %12.5f' % Dsochi_tot_i
    Dsochi_tot_e = D_over_chi_tot(Gamma_e,omn_e,dens_e,Q_e,omt_e,Q_i,omt_i,temp_i, dens_i)
    print 'D_e/chi_tot = %12.5f' % Dsochi_tot_e
    print 'D_i/chi_e =', D_i / chi_e
    print 'D_e/chi_e =', D_e / chi_e
    print 'D_i/chi_i =', D_i / chi_i
    print 'D_e/chi_i =', D_e / chi_i
elif pars['n_spec'] ==3:
    q_i = 1.
    omn_i, omt_i = read_species_gradients(q_i,pars)
    temp_i, dens_i = read_species_tempdens(q_i,pars)
    q_z = pars['charge3']
    omn_z, omt_z = read_species_gradients(q_z,pars)
    temp_z, dens_z = read_species_tempdens(q_z,pars)
    q_e = -1.
    omn_e, omt_e = read_species_gradients(q_e,pars)
    temp_e, dens_e = read_species_tempdens(q_e,pars)

    time, nrgi, nrgz, nrge = read_from_nrg_files(pars,suffix,False)
    print "Actually....species order: first ion, second impurity, third electron"
    Dochi_i, Q_i, Qes_i, Qem_i, Gamma_i, D_i, chi_i = \
        D_over_chi(time,nrgi,omn_i,omt_i,temp_i,dens_i)
    Dochi_e, Q_e, Qes_e, Qem_e, Gamma_e, D_e, chi_e = \
        D_over_chi(time,nrge,omn_e,omt_e,temp_e,dens_e)
    Dochi_z, Q_z, Qes_z, Qem_z, Gamma_z, D_z, chi_z = \
        D_over_chi(time,nrgz,omn_z,omt_z,temp_z,dens_z)
    print 'electron: Q_em/Q_es = %12.5f' % float(Qem_e/Qes_e)
    print 'Q_i/Q_e = %12.5f' % float(Q_i/Q_e)
    print 'chi_i/chi_e =', chi_i / chi_e
    Dsochi_tot_i = D_over_chi_tot(Gamma_i,omn_i,dens_i,Q_e,omt_e,Q_i,omt_i,temp_i, dens_i)
    print 'D_i/chi_tot = %12.5f' % Dsochi_tot_i
    Dsochi_tot_e = D_over_chi_tot(Gamma_e,omn_e,dens_e,Q_e,omt_e,Q_i,omt_i,temp_i, dens_i)
    print 'D_e/chi_tot = %12.5f' % Dsochi_tot_e
    Dsochi_tot_z = D_over_chi_tot(Gamma_z,omn_z,dens_z,Q_e,omt_e,Q_i,omt_i,temp_i, dens_i)
    print 'D_z/chi_tot = %12.5f' % Dsochi_tot_z

    print 'D_i/chi_e =', D_i / chi_e
    print 'D_e/chi_e =', D_e / chi_e
    print 'D_z/chi_e =', D_z / chi_e
    print 'D_i/chi_i =', D_i / chi_i
    print 'D_e/chi_i =', D_e / chi_i
    print 'D_z/chi_i =', D_z / chi_i

