#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy
import matplotlib.pyplot as plt
import re
import optparse as op
from ParIO import * 
from interp import *

parser=op.OptionParser(description='Plots quasilinear estimates of heat flux and compares with nonlinear fluxes.')

options,args=parser.parse_args()
if len(args)!=2:
    exit("""
Please include run number as argument (e.g., 0001) and scanfiles suffix as second argument."
    \n""")
suffix = args[0]
sfsuffix = args[1]

if 'dat' in suffix:
   suffix = '.dat'
elif '_' not in suffix:
   suffix = '_' + suffix

par = Parameters()
par.Read_Pars('parameters'+suffix)
pars = par.pardict
#print pars

if pars['n_spec'] == 1:
   omt = pars['omt1']
else:
   print "Error: Not ready for multiple species!"
   stop

N=pars['nx0']/2+1
pwd=sys.path[0]
f=numpy.loadtxt('fluxspectrae'+suffix+'.dat')
ky1=[]
Qes=[]
for i in range(N,len(f)):
    ky1.append(f[i][0])
    Qes.append(f[i][2])

#g=open('parameters'+suffix,'r')
#while True:
#    h=g.readline()
#    if 'omn' in h:
#        omn=float(re.findall(r"\d+\.?\d*",h)[0])
#    elif 'omt' in h:
#        omt=float(re.findall(r"\d+\.?\d*",h)[0])
#        break
#g.close()

k=numpy.loadtxt('scanfiles'+sfsuffix+'/gamma_kperp2_ratio_kxcenter0')
ky2=[]
ratio=[]
for j in range(len(k)):
    ky2.append(k[j][0])
    ratio.append(k[j][3])

Qql = np.array(ratio)*omt
kygrid = np.linspace(pars['kymin'],(pars['nky0']-1)*pars['kymin'],num = pars['nky0']-1)
print "kygrid",kygrid
#Cubic spline interpolation
Qql_interp1 = interp(ky2,Qql,kygrid)
#Linear interpolation
Qql_interp2 = np.interp(kygrid,ky2,Qql)

Qql_tot1 = np.sum(Qql_interp1)
Qql_tot2 = np.sum(Qql_interp2)
Qnl_tot = np.sum(Qes)
print "Sum Qql interp1:",np.sum(Qql_interp1)
print "Sum Qql linear interp:",np.sum(Qql_interp2)
print "Sum Qnl:",np.sum(Qes)

fig,ax1=plt.subplots()
ax1.plot(ky1,Qes,c='blue',label='Q_es')
plt.legend(loc=2)
ax2=ax1.twinx()
ax2.plot(ky2,Qql,c='red',label='Q_ql, tot: '+str(Qnl_tot)[0:4])
ax2.plot(kygrid,Qql_interp1,'--',c='black',label='interp Qql: '+str(Qql_tot1)[0:4])
ax2.plot(kygrid,Qql_interp2,'--',c='green',label='interp Qql lin: '+str(Qql_tot2)[0:4])
plt.title('C0 = '+str(Qql_tot2/Qnl_tot)[0:4])
plt.xlabel('ky')
plt.legend(loc=1)
plt.gcf().autofmt_xdate()
plt.show()


