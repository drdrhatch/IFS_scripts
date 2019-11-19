#!/usr/bin/env python
# -*- coding: utf-8 -*-

#Test comment for git#
import math
import numpy as np
from parIOWrapper import *
from momentsWrapper import *
from finite_differences import *
from nrgWrapper import *
from fieldlib import *
from momlib import *
import sys
suffix = '_' + sys.argv[1]
pars = init_read_parameters_file(suffix)
moms = momfile('mom_e'+suffix,pars)
time = moms.tmom[:]
nz0 = pars['nz0']
print 'nz0 =', nz0
deln_t = {}
for t in time:
    upar,deln,tpar,tperp,qpar,qperp = moments_from_mom_file(pars,suffix,True,False,t)
    deln_t[t] = deln
#    print 'deln at outboard midplane', deln_t[t][nz0/2]
t_f = time[-1]
print 't_f', t_f
t = []
nrg_e = []
nspec = pars['n_spec']
final_dens = deln_t[t_f][nz0/2]
#print final_dens
with open('nrg'+suffix) as nrg:
    for line in nrg.read().split("\n")[::nspec+1]:
        t.append(line)
#        print line
t.pop()
nrg.close()
#Gamma_es, Gamma_em, Q_es, Q_em = \
#                        read_Gamma_Q(t_f,'nrg_'+suffix,False,)
with open('nrg'+suffix) as nrg:
    for line in nrg.read().split("\n")[2::nspec+1]:
        nrg_e.append(line)
#        print nrg_e
#print range(len(t))
nrg.close()
for i in range(len(t)):
#    print repr(float(t[i].lstrip()))
    t_c =  float(t[i].lstrip())
#    print 't_c',t_c
#    print t_c - t_f
    if abs(t_c - t_f) < 10**-9:
        i_f = i  
        print 't_f', t_f, t_c
#print nrg_e[i_f].split()
nrg_final_data = nrg_e[i_f].split()
#print nrg_final_data
#print 't_f', t_f
#print 'q_es string',nrg_final_data[6]
#print 'q_es float',float(nrg_final_data[6])
Gamma_es = float(nrg_final_data[4])
Gamma_em = float(nrg_final_data[5])
Q_es = float(nrg_final_data[6])
Q_em = float(nrg_final_data[7])
print 'Gamma_es', Gamma_es
print 'Gamma_em', Gamma_em
print 'Q_es', Q_es
print 'Q_em', Q_em
Gamma_tot = Gamma_es + Gamma_em
Q_tot = Q_es + Q_em
Tref = float(pars['Tref'])*10**3
mref = float(pars['mref'])
nref = float(pars['nref'])
Bref = float(pars['Bref'])
Lref = float(pars['Lref'])
rhostar = float(pars['rhostar'])
delta_n = float(final_dens.real)
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
print 'gamma_tot', Gamma_tot
print 'q_tot', Q_tot
print 'rhostar', rhostar
#print 'GENE delta_n/n_0', delta_n
print 'delta_n/n_0', abs(delta_n*rhostar)
#print 'cref', cref
c_si = math.sqrt(T_si/m_si)
Lref = pars['Lref']
R = pars['major_R']*Lref
r = pars['minor_r']*Lref
print 'n0 (m^-3)', n0
print 'Tref (eV)', Tref
print 'c_si', c_si
G_ratio = Gamma_tot*cref*nref/(delta_n**2)
Q_ratio = Q_tot*pref/(delta_n**2)
G_ratio_si = Gamma_tot*c_si*nref/(delta_n**2*Lref**2)
Q_ratio_si = Q_tot*nref*T_si*10**19/(delta_n**2*Lref**2)
Q_gb = n0*c_si*(T_si)*(rhostar**2)
G_gb = n0*c_si*(rhostar**2)
Area = (2*math.pi*R)*(2*math.pi*r*pars['x0'])
Q = Q_gb*Q_tot*Area*10**-6/(delta_n*rhostar)**2
G = G_gb*Gamma_tot*Area*10**-6/(delta_n*rhostar)**2
#print "t_final Density fluctuation at outboard midplane.", delta_n*rhostar*100, "%"
print "Gamma_gb (kW/(eV*m^2))", G_gb*10**-3
print "Q_gb (kW/m^2)",Q_gb*10**-3 
#print "GENE Q_tot", q_tot
#print "Q estimate (MW/m^2)",q_gb*q_tot*10**-6
print "Area (m^2)", Area
print "Lref", Lref  
#print "Gamma in MW: ", g_ratio_si*Area*10**-6
print "Q estimate in MW/(delta_n/n)^2: ", Q
print "Estimated Q(MW) using 1%: ",Q*(0.01)**2 
print "Estimated Q(MW) using 3%: ",Q*(0.03)**2
print "Gamma in MW/(delta_n/n)^2", G
print "Estimated G(MW) using 1%: ",G*(0.01)**2 
print "Estimated G(MW) using 3%: ",G*(0.03)**2

#print "Density Fluctuation Data", dens_data[[dens_data.len()-1]]
#print dens_data.visititems(print_attrs)
