import numpy as np
import matplotlib.pyplot as plt
from write_iterdb import *
from read_iterdb_file import *
from interp import *

fileName = "13left015NegOmTorCMOD1120815027.iterdb"
profilesName = 'profiles_i'
rhot, te, ti, ne, ni, nz, omega_tor = read_iterdb_file(fileName)
data = np.genfromtxt(profilesName)
rhot0 = data[:,0]
rhop0 = data[:,1]
rhop = interp(rhot0,rhop0,rhot)
Ptot = te * ne + ti * (ni + nz)

# collisionality scan
if 1 == 0:
   factor = 1.4
   newTe = 1./factor * te
   newTi = 1./factor * ti
   newNe = factor * ne
   newNi = factor * ni
   newNz = factor * nz
   newPtot = newTe * newNe + newTi * (newNi + newNz)
   
# alpha times omte & omti
if 1 == 1:
    rhotMidPed = 0.97
    midPedIndex = np.argmin(abs(rhot - rhotMidPed))
    teMidPed = te[midPedIndex]
    tiMidPed = ti[midPedIndex]
    print('rhot =', rhotMidPed)
    print('te =', teMidPed)
    print('ti =', tiMidPed)

    alpha = 1.1
    newTe = teMidPed*np.power(te/teMidPed,alpha)
    newTi = tiMidPed*np.power(ti/tiMidPed,alpha)
    
    #new total pressure with new te & ti and old ni, ne & nz
    newPtot = ne * newTe + newTi * (ni + nz)

    #new density profile to keep total pressure the same
    newNe = ne*Ptot/newPtot
    newNi = ni*Ptot/newPtot
    newNz = nz*Ptot/newPtot
    
    #total pressure with new density and tempeturate profiles
    p = newNe * newTe + newTi * (newNi + newNz)

# alpha times omn
if 1 == 0:
    rhotTopPed = 0.97
    topPedIndex = np.argmin(abs(rhot - rhotTopPed))
    neTopPed = ne[topPedIndex]
    niTopPed = ni[topPedIndex]
    nzTopPed = nz[topPedIndex]
    alpha = 0.4
    newNe = neTopPed * np.power(ne / neTopPed, alpha)
    newNi = niTopPed * np.power(ni / niTopPed, alpha)
    newNz = nzTopPed * np.power(nz / nzTopPed, alpha)
    newPtot = newNe * te + (newNi + newNz) * ti
    newTi = ti * Ptot / newPtot
    newTe = te * Ptot / newPtot
    p = newNe * newTe + (newNi + newNz) * newTi

if 1 == 1:
    plt.plot(rhot,ne,label='ne')
    plt.plot(rhot,newNe,label='new ne')
    plt.legend()
    plt.show()

    plt.plot(rhot,te,label='te')
    plt.plot(rhot,newTe,label='new te')
    plt.legend()
    plt.show()

    plt.plot(rhot,Ptot,label='total P')
    plt.plot(rhot,p,label='new total P')
    plt.legend()
    plt.show()

if 1 == 1:
    file_out_base = '04omn'+'13left015NegOmTor'+'CMOD' 
    base_number = '1120815027'
    time_str = '01075'
    output_iterdb(rhot,rhop,newNe*1.E-19,newTe*1.E-3,newNi*1.E-19,newTi*1.E-3,file_out_base+base_number,base_number,time_str,vrot=omega_tor,nimp=newNz*1.E-19)
