#!/usr/bin/env python
# -*- coding: utf-8 -*-

import optparse as op
from parIOWrapper import *
from geomWrapper import *
import sys

parser=op.OptionParser(description='')
options,args=parser.parse_args()
suffix = args[0]
if suffix !='.dat':
   suffix = '_'+suffix

sys.stdout = open('ref' + suffix, 'w')
pars = init_read_parameters_file(suffix)
if 1 == 0:
    geom_type, geom_pars, geom_coeff = init_read_geometry_file(suffix,pars)
    print((list(geom_coeff.keys())))
    print('n0_global (check) =', float(pars['kymin']) * \
          float(geom_pars['Cy']) / float(pars['rhostar']))
          #geom_coeff['C_y'][0]/pars['rhostar']
else:
    geom_type, geom_pars, geom_coeff = init_read_geometry_file_glob(suffix,pars)
    print((list(geom_coeff.keys())))
#    print 'n0_global (check) =', float(pars['kymin']) * \
#          float(geom_pars['Cy']) / float(pars['rhostar'])
          #geom_coeff['C_y'][0]/pars['rhostar']
if 1 == 0:
    # convert k_theta in cm^-1 to ky in the unit of 1/rho_ref
    ktheta = 0.19
    ktheta_GENE = kthetaConversion(ktheta, pars)
    print(('ktheta_cm * rho_cm = ', ktheta_GENE))
    ky_midplane = kthetaConversion(ktheta,pars)
    print('At outboard mid plane ky =', ky_midplane)
    # compute the factor of ky at mid plane to input ky
    k2_fac = k2_factor(geom_type, geom_coeff, False)
    print('kymin input in GENE should be ', ky_midplane / k2_fac)
Bref, Tref, nref, Lref, mref = read_ref_values(suffix, pars)
print('Bref = ', Bref)
print('Tref = ', Tref)
print('nref = ', nref)
print('Lref = ', Lref)
print('mref = ', mref)
Apar_norm = Bref * float(pars['rhostar'])**2 * Lref
print('Apar_norm =', Apar_norm, 'Tesla*m')
pref, cref, Omegaref, rhoref, rhorefStar, Aparref, Gamma_gb, Q_gb = \
derivedRef(suffix, pars)
print(('rho_ref =', rhoref))
print('freq_ref_kHz = ', cref / Lref / 2. / np.pi / 1000.)
print('time_ref_s = ', Lref / cref)
v_alfven = np.sqrt(2./float(pars['beta']))
print(('Alfven speed v_A / v_ref = ', v_alfven))
sys.stdout.close()
if 1 == 1:
    f = open('vAlfven' + suffix,'w')
    f.write('#1.x0 2.beta 3.v_alfven\n')
    ky = float(pars['kymin'])
    np.savetxt(f, np.column_stack((float(pars['x0']), float(pars['beta']), v_alfven)), fmt='%4.8f', delimiter = "  ")
    f.close()

if 1 == 0:
    omn_e, omt_e = read_species_gradients(-1.,pars)
    print('omn_e = ', omn_e)
    print('omt_e = ', omt_e)
    print('eta_e = ', omt_e/omn_e)
    temp_e, dens_e = read_species_tempdens(-1., pars)
    print('temp_e =', temp_e)
    print('dens_e =', dens_e)
    if pars['n_spec'] > 1:
        omn_i, omt_i = read_species_gradients(1.,pars)
        print('omn_i = ', omn_i)
        print('omt_i = ', omt_i)
        temp_i, dens_i = read_species_tempdens(1., pars)
        print('temp_i =', temp_i)
        print('dens_i =', dens_i)
        if pars['n_spec'] > 2:
            omn_z, omt_z = read_species_gradients(pars['charge3'],pars)
            print('omn_z = ', omn_z)
            print('omt_z = ', omt_z)
            temp_z, dens_z = read_species_tempdens(pars['charge3'], pars)
            print('temp_z =', temp_z)
            print('dens_z =', dens_z)
    print('amhd = ', pars['q0']**2*pars['major_R']*pars['beta']*(omn_e+omt_e+omn_i+omt_i))
    p_n = dens_i*temp_i*omn_i + dens_e*temp_e*omn_e + dens_z*temp_z*omn_z
    p_t = dens_i*temp_i*omt_i + dens_e*temp_e*omt_e + dens_z*temp_z*omt_z

    print('pressure from density =', dens_i*temp_i*omn_i + dens_e*temp_e*omn_e + dens_z*temp_z*omn_z)
    print('pressure from density =', p_n)
    print('pressure from temperature =', p_t)
    print('total pressure =', p_n + p_t)
    factor_pt = 0.9
    factor = (p_n + p_t - factor_pt*p_t)/p_n
    print('factor =', factor)
    #beta_chosen = 0.003
    #print 'amhd = ', pars['q0']**2*pars['major_R']*beta_chosen*(omn_e+omt_e+omn_i+omt_i)

