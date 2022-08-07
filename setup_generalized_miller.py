#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import optparse as op
from finite_differences import *
import os
from interp import *
from ParIO import *

parser=op.OptionParser(description='Sets up a simulation with generalized Miller geometry.  Reads in a parameters file that has been set up for tracer_efit.  Reads in profiles_* for each species to calculate gradients and reference quantities.  Arguments: gfile, r(Miller).')
#parser.add_option('--T0','-t',type = 'float',action='store',dest="T0",help = 'Core T(keV).',default=1)

options,args=parser.parse_args()
if len(args)!=2:
    exit("""
Please include (1) efit file name, (2) r(Miller).
    \n""")

gfile = args[0]
r0 = float(args[1])

os.system('extract_generalized_miller_from_eqdsk.py '+gfile+' -r '+str(r0)+' > gmout.dat')

f = open('gmout.dat','r')
data = f.read().split('\n')
f.close()
geom_out_str = ''
Lref_str = 'empty'
Bref_str = 'empty'
grad_factor = 0.0
for i in range(len(data)):
    if 'Coordinates' in data[i]:
        rho_tor = data[i].split()[-1]
        psi_N = data[i].split()[-3]
        psi_N = psi_N[:-1]
        print("r0",r0)
        print("rho_tor",rho_tor)
        print("psi_N",psi_N)
    elif '&geometry' in data[i]:
        for j in range(100):
            geom_out_str += data[i+j]+'\n'
            if data[i+j] and '/' in data[i+j].split()[0]:
                break
    elif 'Lref' in data[i] and 'efit' not in data[i] and 'factor' not in data[i]:            
        Lref_str = data[i]
    elif 'Bref' in data[i] and 'efit' not in data[i]:
        Bref_str = data[i]
    elif 'Gradient conversion' in data[i]:
        grad_factor = float(data[i].split()[-1])
        
 
print("geometry namelist:",geom_out_str)
print("\n\n"+Lref_str)
print("\n\n"+Bref_str)
print("\n\ngrad_factor: ",grad_factor)

par = Parameters()
par.Read_Pars('parameters')
pars = par.pardict
print('pars',pars)

spec_list = []
spec_num_list = []
charge_list = []
for i in pars:
    if 'name' in i:
        spec_list.append(pars[i][1:-1])
        spec_num_list.append(i[-1])
        charge_list.append(pars['charge'+spec_num_list[-1]])
        if charge_list[-1] == -1:
            ename = spec_list[-1]
        elif charge_list[-1] == 1:
            iname = spec_list[-1]
        elif charge_list[-1] > 1:
            zname = spec_list[-1]

print("spec_list",spec_list)
print("ename",ename)
print("iname",iname)
nspec = len(spec_list)

profile_call_str = 'calc_eta_from_gene_profiles.py profiles_'+ename+' profiles_'+iname
profile_call_str += ' '+str(rho_tor)
if nspec == 3:
    profile_call_str += '-i profiles_'+zname

os.system(profile_call_str+' > prof_info.dat')
    
f = open('prof_info.dat','r')
pinfo = f.read().split('\n')
f.close()

for i in pinfo:
    if 'te' in i and 'om' not in i:
        Tref = float(i.split()[-1] )
    if 'ne' in i and 'om' not in i:
        nref = float(i.split()[-1] )
    if 'ti' in i and 'om' not in i:
        ti = float(i.split()[-1] )
    if 'ni' in i and 'om' not in i:
        ni = float(i.split()[-1] )
    if nspec == 3:
        if 'nz' in i and 'om' not in i:
            nz = float(i.split()[-1] )
        if 'tz' in i and 'om' not in i:
            tz = float(i.split()[-1] )
        
    if 'omte' in i:
        omte = float(i.split()[-1] )
    if 'omne' in i:
        omne = float(i.split()[-1] )
    if 'omti' in i:
        omti = float(i.split()[-1] )
    if 'omni' in i:
        omni = float(i.split()[-1] )
    if nspec == 3:
        if 'omnz' in i:
            omnz = float(i.split()[-1] )
        if 'omtz' in i:
            omtz = float(i.split()[-1] )
 
f = open('par_entries_gen_miller.dat','w')
f.write("Instructions:\n1. Delete x0\n2. Replace geometry namelist with block below.\n3. Fill in gradients and reference values with values below.\n\n#################\n")
f.write(geom_out_str)
f.write("\n####################\nReference values:\n")
f.write(Lref_str+'\n')
f.write(Bref_str+'\n')
f.write('nref = '+str(nref)+'\n')
f.write('Tref = '+str(Tref)+'\n')
for name in spec_list:
    if name == ename:
        f.write("\n####################\nSpecies "+name+"\n")
        f.write("dens = 1.0\n")
        f.write("temp = 1.0\n")
        f.write("omt = "+str(omte*grad_factor)+'\n')
        f.write("omn = "+str(omne*grad_factor)+'\n')
    elif name == iname:
        f.write("\n####################\nSpecies "+name+"\n")
        f.write("dens = "+str(ni/nref)+'\n')
        f.write("temp = "+str(ti/Tref)+'\n')
        f.write("omt = "+str(omti*grad_factor)+'\n')
        f.write("omn = "+str(omni*grad_factor)+'\n')
    else:
        f.write("\n####################\nSpecies "+name+"\n")
        f.write("dens = "+str(nz/nref)+'\n')
        f.write("temp = "+str(tz/Tref)+'\n')
        f.write("omt = "+str(omtz*grad_factor)+'\n')
        f.write("omn = "+str(omnz*grad_factor)+'\n')


        




