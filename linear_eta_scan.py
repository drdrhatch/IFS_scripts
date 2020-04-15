import numpy as np
import os
import sys
from ParIO import *

eta_new = 2.2

par = Parameters()

# read in the entire parameters file
param_filename = 'parameters_01'
f = open(param_filename, 'r')
parfile = f.read()
f.close()

par.Read_Pars(param_filename)
pars = par.pardict

# omn, omt, temp, dens for each species

omn_i = float(pars['omn1'])
omt_i = float(pars['omt1'])
temp_i = float(pars['temp1'])
dens_i = float(pars['dens1'])
gradp_ni = dens_i * temp_i * omn_i
gradp_ti = dens_i * temp_i * omt_i

eta_i = omt_i / omn_i
print(('nominal eta_i = ' + str(np.round(eta_i, 4))))

omt_i_new = eta_new * (omt_i + omn_i) / (1. + eta_new)
omn_i_new = (omt_i + omn_i) / (1. + eta_new)

print(('Check omt + omn', omt_i + omn_i, omt_i_new + omn_i_new))
print(('new eta_i = ' + str(np.round(omt_i_new/omn_i_new, 4))))

if int(pars['n_spec']) > 1:
    omn_e = float(pars['omn2'])
    omt_e = float(pars['omt2'])
    temp_e = float(pars['temp2'])
    dens_e = float(pars['dens2'])
    gradp_ne = dens_e * temp_e * omn_e
    gradp_te = dens_e * temp_e * omt_e
    eta_e = omt_e /omn_e
    print(('nominal eta_e = ' + str(np.round(eta_e, 4))))

    omt_e_new = eta_new * (omt_e + omn_e) / (1. + eta_new)
    omn_e_new = (omt_e + omn_e) / (1. + eta_new)

    print(('Check omt + omn', omt_e + omn_e, omt_e_new + omn_e_new))
    print(('new eta_e = ' + str(np.round(omt_e_new/omn_e_new, 4))))

    if int(pars['n_spec']) > 2:
        omn_z = float(pars['omn3'])
        omt_z = float(pars['omt3'])
        temp_z = float(pars['temp3'])
        dens_z = float(pars['dens3'])
        gradp_nz = dens_z * temp_z * omn_z
        gradp_tz = dens_z * temp_z * omt_z
        eta_z = omt_z /omn_z
        print(('nominal eta_z = ' + str(np.round(eta_z, 4))))

        omt_z_new = eta_new * (omt_z + omn_z) / (1. + eta_new)
        omn_z_new = (omt_z + omn_z) / (1. + eta_new)

        print(('Check omt + omn', omt_z + omn_z, omt_z_new + omn_z_new))
        print(('new eta_z = ' + str(np.round(omt_z_new/omn_z_new, 4))))

    

# relative species line index
s_flag = 0 

# s_ind is species index: 0 for ION
# 1 for electron, 2 for impurities
s_ind = -1 
parfile_split = parfile.split('\n')
for i in range(len(parfile_split)):
    #print(i, s_flag, s_ind)
    if s_flag == 1 and 'name' in parfile_split[i]:
        print((parfile_split[i]))
        s_ind += 1
        s_flag += 1
        continue
    if s_flag == 2:
        if s_ind == 0:
            parfile_split[i] = 'omn = '+ str(omn_i_new)
        elif s_ind == 1:
            parfile_split[i] = 'omn = '+ str(omn_e_new)
        elif s_ind == 2:
            parfile_split[i] = 'omn = '+ str(omn_z_new)
        print((parfile_split[i]))
        s_flag += 1
        continue
    if s_flag == 3:
        if s_ind == 0:
            parfile_split[i] = 'omt = '+ str(omt_i_new)
        elif s_ind == 1:
            parfile_split[i] = 'omt = '+ str(omt_e_new)
        elif s_ind == 2:
            parfile_split[i] = 'omt = '+ str(omt_z_new)
        print((parfile_split[i]))
        s_flag = 0
        continue
    if 'species' in parfile_split[i]:
        s_flag = 1
        continue

parfile_out='\n'.join(parfile_split)

if 1 == 1:
    f=open('parameters_test','w')
    f.write(parfile_out)
    f.close()

