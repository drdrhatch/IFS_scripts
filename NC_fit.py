#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import optparse as op
from ParIO import * 
import scipy.optimize as optimize

parser=op.OptionParser(description='Makes estimates of NC particle flux based on a neoclassical simulation that is not entirely converged.')
parser.add_option('--number','-n',type='int',action='store',dest='number',help = 'Select column 1.Gamma 2.Q 3.Momentum 4.Jboot.',default=1)
parser.add_option('--niterations','-i',type='int',action='store',dest='niterations',help = 'Number of variations for each initial guess.',default=10)
#parser.add_option('--time','-t',type = 'float',action='store',dest="time0",help = 'Time to plot mode structure.',default=-1)
options,args=parser.parse_args()
print("options.number",options.number)
column = int(options.number)-1
niterations = int(options.niterations)

if column == 0:
    print("Calculating particle flux.")
elif column == 1:
    print("Calculating heat flux.")
elif column == 2:
    print("Calculating momentum flux.")
elif column == 3:
    print("Calculating Jboot.")

if len(args)!=1:
    exit("""
Please enter file suffix\n""")

suffix = args[0]

if 'dat' in suffix:
   suffix = '.dat'
elif '_' not in suffix:
   suffix = '_'+suffix

filename = 'neoclass'+suffix

par = Parameters()
par.Read_Pars('parameters'+suffix)
pars = par.pardict

f=open(filename,'r')
data = f.read()

lines = data.split('\n')

keep_going = True
i=2
while keep_going:
    ls = lines[i].split()
    if len(ls) == 1:
        num_spec = i-2
        keep_going = False
    i += 1

print("Number of species:",num_spec )
    
time = np.empty(0)
# delete header
lines = np.delete(lines,0)
temp = lines[0::num_spec+1]
for i in range(len(temp)):
    if temp[i]:
        time = np.append(time,float(temp[i]))
ntime = len(time)
fluxes = np.empty((ntime,4,num_spec))
for i in range(num_spec):
    temp  = lines[i+1::num_spec+1]
    for j in range(len(temp)):
        if temp[j]:
            ls = np.array(temp[j].split())
            ls = ls.astype(np.float)
            fluxes[j,:,i] = ls[:]

ambipolarity = 0
for i in range(num_spec):
    print("species "+pars['name'+str(i+1)])
    print( "Charge of species "+str(i+1)+": "+str(pars['charge'+str(i+1)]))
    print("Flux/FluxGB",fluxes[-1,column,i])
    ambipolarity += fluxes[-1,column,i]*pars['charge'+str(i+1)] 
   
print("Ambipolarity:",ambipolarity)
print("Normalized ambipolarity:",ambipolarity/(abs(fluxes[-1,column,0])+abs(fluxes[-1,column,1])+abs(fluxes[-1,column,2])))

def fit_func(t,G0,c0,gam):
    return G0+c0*np.e**(-gam*t)

bvals = np.empty((3,num_spec))
for s in range(num_spec):

    G0 = fluxes[-1,column,s]
    c0 = fluxes[0,column,s] - G0
    #gam = 0.001
    gam = (fluxes[-2,column,s]-fluxes[-1,column,s])/(time[-1]-time[-2])/c0
    if s > 0:
        gam = bvals[2,s-1]
    print("gam0",gam)
    plt.plot(time,fluxes[:,column,s])
    plt.plot(time,G0+c0*np.e**(-gam*(time-time[0]))  )
    plt.show()

    f0 = np.linspace(0.1*G0,10*G0,num=niterations)
    f1 = np.linspace(-10*c0,10*c0,num=niterations)
    f2 = np.linspace(0.1*gam,10*gam,num=niterations)

    minerr = 1.0
    prec = 1e-6
    for i in range(len(f0)):
        if minerr < prec:
            break
        for j in range(len(f1)):
            if minerr < prec:
                break
            for k in range(len(f2)):
                x0 = np.array([f0[i],f1[j],f2[k]])
                print("i,j,k",i,j,k)
                fit,pcov= optimize.curve_fit(fit_func,time-time[0],fluxes[:,column,s],p0=x0,maxfev = 10000) 
                thiserr = np.sum(np.diag(pcov))
                if abs(thiserr) < minerr:
                    minerr = thiserr
                    bestvals = np.array([fit[0],fit[1],fit[2]])
                    print("##############")
                    print(i,j,k)
                    print("fit",fit[0],fit[1],fit[2])
                    print("minerr",minerr)
                    if minerr < prec:
                        break

    plt.plot(time-time[0],fluxes[:,column,s],label='data')
    plt.plot(time-time[0],G0+c0*np.e**(-gam*(time-time[0])),label='guess' )
    plt.plot(time-time[0],bestvals[0]+bestvals[1]*np.e**(-bestvals[2]*(time-time[0])),'--',label='fit' )
    plt.legend()
    plt.show()
    bvals[:,s] = bestvals 



print("bvals",bvals)

ambipolarity = 0
for s in range(num_spec):
    print("species "+pars['name'+str(s+1)])
    print("Charge of species "+str(s+1)+": "+str(pars['charge'+str(s+1)]))
    print("Flux/FluxGB",bvals[0,s])
    ambipolarity += bvals[0,s]*pars['charge'+str(s+1)] 
   
if column==0:
    print("Ambipolarity:",ambipolarity)
    print("Normalized ambipolarity:",ambipolarity/np.sum(abs(bvals[0,:])))


#proceed = input("Proceed with fit (0 = no)?:")
#if proceed != '0':
#    gam0 = (fluxes[-1,0,i] - fluxes[-2,0,i]) / (time[-1] - time[-2])
#    start_time = float(input("Enter start time:"))
#    istart = np.argmin(abs(time - start_time))    
#    x0 = np.array([1.0,1.0,0.01])
#    for i in range(num_spec):
#        fit,dummy= optimize.curve_fit(fit_func,time[istart:],fluxes[istart:,0,i],x0)      
#        print("fit",fit)
#        plt.plot(time,fluxes[:,0,i],label='data')
#        plt.plot(time,fit[0]+fit[1]*np.e**(-fit[2]*time),label='fit')
#        plt.legend()
#        plt.show()
    
if column == 0:
    dGdt = np.empty(num_spec)
    q = np.empty(num_spec)
    num = 0
    denom = 0
    for s in range(num_spec):
        dGdt[s] = (fluxes[-1,column,s] - fluxes[0,column,s])/(time[-1] - time[0])
        q[s] = pars['charge'+str(s+1)] 
        denom += q[s]*dGdt[s] 
        num += q[s]*fluxes[0,column,s]

    t_ambi = -num/denom

    ambi_test = 0
    Gs = np.empty(num_spec)
    print("Estimates based on linear extrapolation to ambipolar condition.")
    for s in range(num_spec):
        temp = fluxes[0,column,s]*q[s] + q[s]*dGdt[s]*t_ambi
        Gs[s] = fluxes[0,column,s] + dGdt[s]*t_ambi
        #print(pars['charge'+str(s+1)],"Gs",Gs[s]) 
        print(pars['name'+str(s+1)],"Gs",Gs[s]) 
        ambi_test += temp

    print("Test of ambipolarity:",ambi_test)

tnew = np.linspace(0,10000,num=1000)
for s in range(num_spec):
    plt.plot(tnew,bvals[0,s]+bvals[1,s]*np.e**(-bvals[2,s]*(tnew)),'x-',label='Exponential fit '+pars['name'+str(s+1)] )
    plt.plot(time-time[0],bvals[0,s]+bvals[1,s]*np.e**(-bvals[2,s]*(time-time[0]))  )
    #plt.hlines(Gs[s],0,10000,label='Ambipolar estimate '+pars['name'+str(s+1)])
    if column == 0:
       plt.plot(tnew,tnew/tnew*Gs[s],label='Ambipolar estimate '+pars['name'+str(s+1)])
plt.legend()
plt.title("Column: "+str(column))
plt.show()



