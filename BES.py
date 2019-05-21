#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import math
import numpy as np
from parIOWrapper import init_read_parameters_file
import sys
suffix = '_' + sys.argv[1]
file = h5py.File('mom_e'+suffix+'.h5','r')
keys = list(file.keys())
dset = file[keys[0]]
pars = init_read_parameters_file(suffix)
nz0 = pars['nz0']
print 'nz0 =', nz0
dens_data = file['mom_e/dens/']
fluc_data = []
for subgroup in dens_data:
    sub_dens = file['mom_e/dens/'+subgroup]
    fluc_data.append(sub_dens[nz0/2])
time = file['mom_e/time']
t_f = float(time[-1])

t = []
nrg_e = []
final_dens = list(fluc_data[-1][0])
with open('nrg'+suffix) as nrg:
    for line in nrg.read().split("\n")[::4]:
        t.append(line)
#        print line
t.pop()
nrg.close()
with open('nrg'+suffix) as nrg:
    for line in nrg.read().split("\n")[2::4]:
        nrg_e.append(line)
#        print nrg_e
#print range(len(t))
nrg.close()
for i in range(len(t)):
#    print t[i].lstrip()
    t_c =  float(t[i].lstrip())
#    print 't_c',t_c
#    print t_c - t_f
    if abs(t_c - t_f) < 10**-10:
        i_f = i  
        print 't_f', t_f, t_c
#print nrg_e[i_f].split()
nrg_final_data = nrg_e[i_f]
#print 't_f', t_f

gamma_es = float(nrg_final_data[4])
gamma_em = float(nrg_final_data[5])
q_es = float(nrg_final_data[6])
q_em = float(nrg_final_data[7])

gamma_tot = gamma_es + gamma_em
q_tot = q_es + q_em
Tref = float(pars['Tref'])*10**3
mref = float(pars['mref'])
nref = float(pars['nref'])
Bref = float(pars['Bref'])
Lref = float(pars['Lref'])
rhostar = float(pars['rhostar'])
delta_n = float(final_dens[0][0])
cref = math.sqrt(Tref/mref)
n0 = nref*10**19
pref = n0*Tref
#print 'Tref', Tref
#print 'mref', mref
#print 'n0', n0
#print 'rhostar', rhostar
#print 'pref', nref*Tref
#print 'delta_n', delta_n
#print 'n1_real', final_dens[0][0]
c = 3*10**8
k = 1.38*(10**-23)
e = 1.6*10**-19
mp = 1.67*10**-27
T_si = e*Tref
m_si = mp*mref
#print 'T_si', T_si
#print 'm_si', m_si
#print 'gamma_tot', gamma_tot
#print 'q_tot', q_tot
print 'rhostar', rhostar
print 'GENE delta_n/n_0', delta_n
print 'delta_n/n_0', delta_n*rhostar
#print 'cref', cref
c_si = math.sqrt(T_si/m_si)
Lref = pars['Lref']
R = pars['major_R']*Lref
r = pars['minor_r']*Lref
print 'n0 (m^-3)', n0
print 'Tref (eV)', Tref
print 'c_si', c_si
g_ratio = gamma_tot*cref*nref/(delta_n**2)
q_ratio = q_tot*pref/(delta_n**2)
g_ratio_si = gamma_tot*c_si*nref/(delta_n**2*Lref**2)
q_ratio_si = q_tot*nref*T_si*10**19/(delta_n**2*Lref**2)
q_gb = n0*c_si*(T_si)*(rhostar**2)
g_gb = n0*c_si*(rhostar**2)
Area = (2*math.pi*R)*(2*math.pi*r*pars['x0'])
q = q_gb*q_tot*Area*10**-6/(delta_n*rhostar)**2
#print "t_final Density fluctuation at outboard midplane.", delta_n*rhostar*100, "%"
print "Gamma_gb (kW/(eV*m^2))", g_gb*10**-3
print "Q_gb (kW/m^2)",q_gb*10**-3 
#print "GENE Q_tot", q_tot
#print "Q estimate (kW/m^2)",q_gb*q_tot*10**-3
print "Area (m^2)", Area
print "Lref", Lref  
#print "Gamma in MW: ", g_ratio_si*Area*10**-6
print "Q estimate in MW/(delta_n/n)^2: ", q
print "Estimated Q(MW) using 1%: ",q*(0.01)**2 
print "Estimated Q(MW) using 3%: ",q*(0.03)**2
#print "Density Fluctuation Data", dens_data[[dens_data.len()-1]]
#print dens_data.visititems(print_attrs)
