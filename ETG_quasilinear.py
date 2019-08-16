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
parser.add_option('--include_qn2','-q', action='store_true',dest = 'include_qn2', help = 'Includes the weight factor Q/n^2', default=False)
parser.add_option('--sat_rule','-s', action='store',type = int,dest = 'sat_rule', help = 'Selects saturation rule: 1=gamma/kperp2, 2=gamma/(kx^2+ky^2).', default=1)

options,args=parser.parse_args()
if len(args)!=2:
    exit("""
Please include run number as argument (e.g., 0001) and scanfiles suffix as second argument."
    \n""")
suffix = args[0]
sfsuffix = args[1]
include_qn2 = options.include_qn2
sat_rule = options.sat_rule
print "include_qn2",include_qn2

plot_title = 'QL spectrum, sat_rule: '+str(sat_rule)

if include_qn2:
    plot_title += ' with Q/n^2 '

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

ldata = np.genfromtxt('scanfiles'+sfsuffix+'/mode_info_all')

if sat_rule == 1:
    Qql0 = ldata[:,6]*omt
elif sat_rule == 2:
    Qql0 = ldata[:,4]/(ldata[:,0]**2+ldata[:,10]**2)
elif sat_rule == 3:
    Qql0 = np.empty(len(ldata[:,0]))
    for i in range(len(ldata[:,0])):
        denom = max(ldata[i,0]*ldata[i,10],ldata[i,0]**2)
        Qql0[i] = ldata[i,4]/denom
elif sat_rule == 4:
    Qql0 = ldata[:,6]*omt 
    for i in range(len(ldata[:,0])):
        if ldata[i,10] < 0.5*ldata[i,0]:
            Qql0[i] = ldata[i,4]/(ldata[i,0]*ldata[i,2])





Qql = np.empty(0)
ky2 = np.empty(0)
Qn2 = np.empty(0)
nky = 0
for i in range(len(ldata[:,0])):
    if not ldata[i,0] in ky2: 
        ky2 = np.append(ky2,ldata[i,0])
        Qql = np.append(Qql,Qql0[i])
        Qn2 = np.append(Qn2,ldata[i,7])
        nky += 1
    else:
        if Qql0[i] > Qql[nky-1]:
            Qql[nky-1] = Qql0[i]  #Taking maximum over ballooning angle / eigenvalue
            Qn2[nky-1] = ldata[i,7]

#gkp = np.genfromtxt('scanfiles'+sfsuffix+'/gamma_kperp2_ratio_kxcenter0')
#kperp2 = gkp[:,2]
#Test kperp^3
#Qql2 = Qql*kperp2/kperp2**1.5
if include_qn2:
    Qql = Qql*Qn2

kygrid = np.linspace(pars['kymin'],(pars['nky0']-1)*pars['kymin'],num = pars['nky0']-1)
print "kygrid",kygrid
#Cubic spline interpolation
Qql_interp1 = interp(ky2,Qql,kygrid)
#Linear interpolation
Qql_interp2 = np.interp(kygrid,ky2,Qql)
#Qql2_interp2 = np.interp(kygrid,ky2,Qql2)

Qql_tot1 = np.sum(Qql_interp1)
Qql_tot2 = np.sum(Qql_interp2)
Qnl_tot = np.sum(Qes)
#Qql2_tot = np.sum(Qql2_interp2)
print "Sum Qql interp1:",np.sum(Qql_interp1)
print "Sum Qql linear interp:",np.sum(Qql_interp2)
#print "Sum Qql (kp3) linear interp:",np.sum(Qql2_interp2)
print "Sum Qnl:",np.sum(Qes)
#c3 = Qql_tot2/Qql2_tot

fig,ax1=plt.subplots()
ax1.plot(ky1,Qes,c='blue',label='Q_es')
plt.legend(loc=2)
ax2=ax1.twinx()
ax2.plot(ky2,Qql,c='red',label='Q_ql, tot: '+str(Qnl_tot)[0:4])
#ax2.plot(kygrid,Qql_interp1,'--',c='black',label='interp Qql: '+str(Qql_tot1)[0:4])
ax2.plot(kygrid,Qql_interp2,'--',c='green',label='interp Qql lin: '+str(Qql_tot2)[0:4])
#ax2.plot(kygrid,c3*Qql2_interp2,'--',c='purple',label=str(c3)+' x '+'+Qql kp3: '+str(Qql2_tot)[0:4])
#plt.title('C0 = '+str(Qql_tot2/Qnl_tot)[0:4]+'  c3 = '+str(c3))
#plt.title('C0 = '+str(Qql_tot2/Qnl_tot)[0:4])
plt.xlabel('ky')
plt.legend(loc=1)
plt.gcf().autofmt_xdate()
plt.title(plot_title + ', C0 = '+str(Qql_tot2/Qnl_tot)[0:4])
plt.show()

#plt.plot(ldata[:,0],ldata[:,4],'x')
#plt.xlabel('kymin')
#plt.ylabel('gamma')
#plt.show()
#
#plt.plot(ldata[:,0],ldata[:,2],'x')
#plt.xlabel('kymin')
#plt.ylabel('kperp')
#plt.show()
#
#plt.plot(ldata[:,0],ldata[:,6],'x')
#plt.xlabel('kymin')
#plt.ylabel('gamma/kperp2')
#plt.show()
#
#plt.plot(ldata[:,0],ldata[:,7],'x')
#plt.xlabel('kymin')
#plt.ylabel('Q/n^2')
#plt.show()
