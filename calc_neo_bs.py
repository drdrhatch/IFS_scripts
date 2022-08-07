#!/usr/bin/env python
# -*- coding: utf-8 -*-
#import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import optparse as op
from finite_differences import *
import os
from interp import *
from read_pfile import *

parser=op.OptionParser(description='Calculates the bootstrap current using neo at a selected value of psi_norm.  Inputs: gfile, pfile, psi')
#parser.add_option('--T0','-t',type = 'float',action='store',dest="T0",help = 'Core T(keV).',default=1)

options,args=parser.parse_args()
if len(args)!=3:
    exit("""
Please include (1) efit file name, (2) pfile name, and (3) value of psinorm  as arguments.
    \n""")

gfile = args[0]
pfile = args[1]
psinorm = args[2]

print("EQDSK file:",gfile)
print("p file:",pfile)
print("psinorm:",psinorm)

os.system('profiles_gen -g '+gfile+' -i '+pfile)
os.system('profiles_gen -g '+gfile+' -i '+pfile+' -loc_psi '+psinorm+' > pgout.dat')
f = open('pgout.dat','r')
data = f.read()
f.close()

data = data.split('\n')
Ti_info = ''
Te_info = ''
Bunit = 0.0

for i in range(10):
    if 'Ti' in data[-1-i]:
        Ti_info = data[-1-i]
    elif 'Te' in data[-1-i]:
        Te_info = data[-1-i]
    elif 'Bunit' in data[-1-i]:
        line = data[-1-i].split()
        Bunit = float(line[-1])

print("Ti_info",Ti_info)
print("Te_info",Te_info)
print("Bunit",Bunit)

psinorm = float(psinorm)
def read_pg_pfile(suffix):
    f = open('pfile.'+suffix)
    data = f.read()
    data = data.split('\n')
    nr = int(data[0])
    psinorm = np.empty(nr)
    var = np.empty(nr)
    for i in range(nr):
        psinorm[i] = float(data[i+1].split()[0])
        var[i] = float(data[i+1].split()[1])
    return psinorm,var
    
f.close()
#ne is in 10^20 m^-3
psi_ne,ne = read_pg_pfile('ne')
psi_ni,ni = read_pg_pfile('ni')
psi_te,te = read_pg_pfile('te')
psi_ti,ti = read_pg_pfile('ti')

#print('psi_ne',psi_ne)
#print('ne',ne)

ine0 = np.argmin(abs(psi_ne-psinorm))
ne0 = ne[ine0]

ini0 = np.argmin(abs(psi_ni-psinorm))
ni0 = ni[ini0]

ite0 = np.argmin(abs(psi_te-psinorm))
te0 = te[ite0]

iti0 = np.argmin(abs(psi_ti-psinorm))
ti0 = ti[iti0]

print("ne0(1e20m^-3)",ne0)
print("ni0(1e20m^-3)",ni0)
print("te0(1e20m^-3)",te0)
print("ti0(1e20m^-3)",ti0)

os.system('cp input.neo.locpargen input.neo')
os.system('neo -e ./')

f = open('out.neo.transport','r')
ntrans = f.read().split()
JBS_neo = ntrans[2]
f.close()
print("Normalized bootstrap current from NEO:",JBS_neo)

print("Identifying first species")
f = open('out.neo.run','r')
data = f.read().split('\n')
for i in range(len(data)):
    if 'r/a:' in data[i]:
        r_over_a = data[i]
    if 'indx' in data[i]:
        print(data[i])
        print(data[i+1])
        temp = data[i+1].split()
        z = float(temp[1])
        if z==1:
            spec1 = 'i'
        elif z==-1:
            spec1 = 'e'
        else:
            print("Error! First species is not main ion or electron.")
            stop
        break
print("First species is "+spec1+'.')
e = 1.602e-19
md = 3.3452e-27
if spec1 == 'i':
    vnorm = np.sqrt(ti0*e*1000/md)
    nnorm = ni0*1e20
else:
    vnorm = np.sqrt(te0*e*1000/md)
    nnorm = ne0*1e20

f = open('out.neo.theory','r')
ttrans = f.read().split()
f.close()
JBS_HH = ttrans[4]
JBS_S = ttrans[10]
JBS_K = ttrans[17]
print("Normalized bootstrap current from Hinton Hazeltine:",JBS_HH)
print("Normalized bootstrap current from Sauter:",JBS_S)
print("Normalized bootstrap current from Koh:",JBS_K)
print("vnorm",vnorm)

JBS_neo_SI = float(JBS_neo)*e*nnorm*vnorm
print("JBS_neo_SI",JBS_neo_SI)

###Now calculate a measure of Jtor from the equilibrium
os.system('efit_tools.py '+gfile+' -p -n')
binfo = np.genfromtxt('Binfo_'+gfile)
Bpol_max = np.max(binfo[:,2])
i0 = np.argmin(abs(binfo[:,1]-psinorm))
Bpol0 = binfo[i0,2]
Btor0 = binfo[i0,3]
B0 = np.sqrt(Bpol0**2+Btor0**2)
a0 = binfo[-1,0]-binfo[0,0]
mu0 = 1.2566e-6
print("Maximum poloidal B:",Bpol_max)
print("Bpol at psi0:",Bpol0)
print("Btor at psi0:",Btor0)
print("B0 at psi0:",B0)
print("a outboard (m):",a0)
Javg = 2*Bpol_max/mu0/a0
print("Javg (A/m^2)",Javg)
print("JBS_neo_SI*Bunit/Bloc",JBS_neo_SI*Bunit/B0)
Javg_Bnorm = Javg*B0/Bunit
print("Javg Bloc/Bunit (A/m^2)",Javg_Bnorm)

###Now put results in a folder and write a summary to file
os.system('mkdir results_'+str(psinorm))
os.system('cp out.neo.* results_'+str(psinorm))
f = open('./results_'+str(psinorm)+'/bootstrap_summary','w')
f.write("nnorm(m^-3): "+str(nnorm)+'\n')
f.write("mnorm(kg): "+str(md)+'\n')
f.write("vnorm(m/s): "+str(vnorm)+'\n')
f.write("Bunit(T): "+str(Bunit)+'\n')
f.write("B_local(T) efit: "+str(B0)+'\n')
f.write("a outboard (m) efit: "+str(a0)+'\n')
f.write("psi: "+str(psinorm)+'\n')
f.write(str(r_over_a)+'\n')
f.write("##################"+'\n')
f.write("JBS neo: "+str(JBS_neo)+'\n')
f.write("JBS HH: "+str(JBS_HH)+'\n')
f.write("JBS Sauter: "+str(JBS_S)+'\n')
#f.write("JBS Koh: "+str(JBS_K)+'\n')
f.write("####SI units\n")
f.write("JBS neo (<Jpar B>/Bunit) (A/m^2): "+str(JBS_neo_SI)+'\n')
f.write("Javg efit (<Javg B_local>/Bunit) (A/m^2): "+str(Javg_Bnorm)+'\n')
f.write("####Without Bunit\n")
f.write("JBS neo (<Jpar B>/B_local) (A/m^2): "+str(JBS_neo_SI*Bunit/B0)+'\n')
f.write("Javg efit (2 Bpolmax/(mu0 a_outboard)) (A/m^2): "+str(Javg)+'\n')
f.close()






