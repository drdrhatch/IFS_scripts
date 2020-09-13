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
parser.add_option('--add_all','-a', action='store_true',dest = 'add_all', help = 'Add all contributions for each ky (kx_center scan or EV run).', default=False)
parser.add_option('--noplot','-n', action='store_true',dest = 'noplot', help = 'No plot.', default=False)
parser.add_option('--most_unstable','-m', action='store_true',dest = 'most_unstable', help = 'Takes contribution only from the most unstable mode.', default=False)
#parser.add_option('--width_factor','-w', action='store_true',dest = 'width_factor', help = 'Includes the width factor in the calculation.', default=False)

options,args=parser.parse_args()
if len(args)!=2:
    exit("""
Please include run number as argument (e.g., 0001) and scanfiles suffix as second argument."
    \n""")
suffix = args[0]
sfsuffix = args[1]
include_qn2 = options.include_qn2
add_all = options.add_all
sat_rule = options.sat_rule
#width_factor = options.width_factor
noplot = options.noplot
most_unstable = options.most_unstable
print ("include_qn2",include_qn2)
print ("add_all",add_all)
print ("sat_rule",sat_rule)
print ("most_unstable",most_unstable)


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
elif pars['n_spec'] == 2:
   print ("Warning: Not ready for multiple species!")
   omt = pars['omt2']

N=int(np.floor(pars['nx0']/2))+1
pwd=sys.path[0]
f=numpy.loadtxt('fluxspectrae'+suffix+'.dat')
ky1=[]
Qes=[]
for i in range(N,len(f)):
    ky1.append(f[i][0])
    Qes.append(f[i][2])

ldata = np.genfromtxt('scanfiles'+sfsuffix+'/mode_info_all')
wf = ldata[:,11]  #width
wf = wf/np.max(wf)
parlin = Parameters()
parlin.Read_Pars('scanfiles'+sfsuffix+'/parameters')
parslin = parlin.pardict
IVEV = parslin['comp_type']
print ("IVEV",IVEV)

if sat_rule == 1: #Conventional gamma / kperp2
    Qql0 = ldata[:,6]*omt  
elif sat_rule == 2: #Using kx**2 + ky**2 instead of kperp2
    Qql0 = ldata[:,4]*omt/(ldata[:,0]**2+ldata[:,10]**2)
elif sat_rule == 3:  #Maximum of kx*ky vs ky**2
    Qql0 = np.empty(len(ldata[:,0]))
    for i in range(len(ldata[:,0])):
        denom = max(ldata[i,0]*ldata[i,10],ldata[i,0]**2)
        Qql0[i] = ldata[i,4]*omt/denom
elif sat_rule == 4:  #Standard rule except gamma / ky**2 if kperp < 0.5 ky
    Qql0 = ldata[:,6]*omt 
    for i in range(len(ldata[:,0])):
        #if ldata[i,10] < 0.5*ldata[i,0]:
        if ldata[i,2] < 0.5*ldata[i,0]:
            #Qql0[i] = ldata[i,4]/(ldata[i,0]*ldata[i,2])*omt
            Qql0[i] = ldata[i,4]/(ldata[i,0]**2)*omt
elif sat_rule == 5: #gamma / ky*kperp
    Qql0 = np.empty(len(ldata[:,0]))
    for i in range(len(ldata[:,0])):
        denom = ldata[i,0]*ldata[i,2]
        Qql0[i] = ldata[i,4]/denom*omt
elif sat_rule == 6: #gamma / 2(ky**2+kperp**2)
    Qql0 = ldata[:,4]/(0.5*ldata[:,0]**2+0.5*ldata[:,2]**2)*omt
elif sat_rule == 7:
    Qql0 = ldata[:,6]*omt*abs(ldata[:,5])/ldata[:,4]
elif sat_rule == 8:  #Standard with width factor
    Qql0 = ldata[:,6]*omt*wf[:]

if include_qn2:
    Qql0 = ldata[:,7]*Qql0


Qql = np.empty(0)
ky2 = np.empty(0)
nky = 0
for i in range(len(ldata[:,0])):
    if not ldata[i,0] in ky2 and ldata[i,9] == 0.0: 
        ky2 = np.append(ky2,ldata[i,0])
        Qql = np.append(Qql,Qql0[i])
        nky += 1
        print("wf[i]",wf[i])
        print("ky[i]",ldata[i,0])
    else:
        if add_all and Qql0[i] > 0 and not most_unstable:
            Qql[nky-1] += Qql0[i]  #Summing at each ky for all positive values
        elif Qql0[i] > Qql[nky-1] and not most_unstable:
            Qql[nky-1] = Qql0[i]  #Taking maximum over ballooning angle / eigenvalue

#gkp = np.genfromtxt('scanfiles'+sfsuffix+'/gamma_kperp2_ratio_kxcenter0')
#kperp2 = gkp[:,2]
#Test kperp^3
#Qql2 = Qql*kperp2/kperp2**1.5
#if include_qn2:
#    Qql = Qql*Qn2

kygrid = np.linspace(pars['kymin'],(pars['nky0']-1)*pars['kymin'],num = pars['nky0']-1)
#print "kygrid",kygrid
#print "ky1",ky1
#Cubic spline interpolation
Qql_interp1 = interp(ky2,Qql,kygrid)
#Linear interpolation
Qql_interp2 = np.interp(kygrid,ky2,Qql)
#Qql2_interp2 = np.interp(kygrid,ky2,Qql2)

Qql_tot1 = np.sum(Qql_interp1)
Qql_tot2 = np.sum(Qql_interp2)
Qnl_tot = np.sum(Qes)
#Qql2_tot = np.sum(Qql2_interp2)
print ("Sum Qql interp1:",np.sum(Qql_interp1))
print ("Sum Qql linear interp:",np.sum(Qql_interp2))
#print "Sum Qql (kp3) linear interp:",np.sum(Qql2_interp2)
print ("Sum Qnl:",np.sum(Qes))
C0 = np.sum(Qes)/np.sum(Qql_interp2)
print ("C0:",C0)
ikpeak = np.argmax(Qql_interp2)
print ("Qql peak ky",kygrid[ikpeak])
C0peak = np.sum(Qes)/np.sum(Qql_interp2[ikpeak])
print ("C0peak",C0peak)

#c3 = Qql_tot2/Qql2_tot
outfile = 'QL_summary'+suffix+'_'+sfsuffix+'_iqn2'+str(include_qn2)[0]+'_sat'+str(sat_rule)+'_add'+str(add_all)[0]+'_mu'+str(most_unstable)[0]

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
plt.title(outfile + ', C0 = '+str(Qnl_tot/Qql_tot2)[0:6])
if not noplot:
    plt.show()
else:
    plt.savefig(outfile+'.ps',format='ps')

f=open(outfile,'w')
f.write('#comp_type: '+IVEV+'\n')
f.write('#suffix: '+str(suffix)+'\n')
f.write('#sfsuffix: '+str(sfsuffix)+'\n')
f.write('#include_qn2: '+str(include_qn2)+'\n')
f.write('#sat_rule: '+str(sat_rule)+'\n')
f.write('#add_all: '+str(add_all)+'\n')
f.write('#most_unstable: '+str(most_unstable)+'\n')
f.write('#C0 = '+str(C0)[0:6]+'\n')
f.write('#C0peak = '+str(C0peak)[0:6]+'\n')
f.write('#1.ky 2.Qnl 3.Qql_interp\n')
np.savetxt(f,np.column_stack((kygrid,Qes[1:],Qql_interp2)))
f.close()

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
