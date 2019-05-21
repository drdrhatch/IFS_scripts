import numpy as np
import matplotlib.pyplot as plt
from write_iterdb import *
from read_iterdb_file import *
from interp import *
from finite_differences import *

#modes:
#mode = 'eta': increase omt(e,i) by alpha, decrease omn to keep pressure fixed
#mode = 'etaTe': increase omte by alpha, decrease omti 
#mode = 'coll': vary collisionality at fixed pressure by decreasing/increasing temperature/density by alpha
#mode = 'TiTe': increase Te by alpha, decrease Ti by alpha
#mode = 'omnz': increase omnz (i.e. impurity density gradient) by alpha and decrease omni 
#mode = 'omte': increase omte without compensation

mode = 'omte'
alpha = 0.5
target_factor = 0.5

#===========================================================
#===========================================================
#===========================================================

fileName = 'jet78697.51005_hager_Z6.0Zeff2.35__negom.iterdb'
profilesName = 'gene_profiles_e_jet78697'
file_out_base = 'jet78697.51005_hager_Z6.0Zeff2.35__negom'
base_number = '78697'
rhotMidPed = 0.97



#fileName = 'geqdsk_89454_Z4.0_negom.iterdb'
#profilesName = 'gene_profiles_e_89454'
#file_out_base = 'geqdsk_89454_Z4.0_negom'
#base_number = '89454'
#rhotMidPed = 0.96


#fileName = 'jet92432.49152_hager_Be_Z6.0Zeff1.03__negom.iterdb'
#profilesName = 'gene_profiles_e_jet92432'
#file_out_base = 'jet92432.49152_hager_Be_Z6.0Zeff1.03__negom'
#base_number = '92432'
#rhotMidPed = 0.975

#fileName = 'jet92432.49152_hager_Be_Z4.0Zeff1.8__negom.iterdb'
#profilesName = 'gene_profiles_e_jet92432.49152_hager_Be'
#file_out_base = 'jet92432.49152_hager_Be_Z4.0Zeff1.8__negom'
#base_number = '92432'
#rhotMidPed = 0.97

#fileName = "jet92300.51580_hager.eq_Z4.0Zeff1.4_.iterdb"
#profilesName = 'gene_profiles_e_jet92300'
#file_out_base = 'jet92300.51580_hager.eq_Z4.0Zeff1.4'
#rhotMidPed = 0.97
#alpha = 0.95

rhot, te, ti, ne, ni, nz, omega_tor = read_iterdb_file(fileName)
data = np.genfromtxt(profilesName)
rhot0 = data[:,0]
rhop0 = data[:,1]
rhop = interp(rhot0,rhop0,rhot)
Ptot = te * ne + ti * (ni + nz)

# collisionality scan
if mode == 'coll':
   newTe = 1./alpha * te
   newTi = 1./alpha * ti
   newNe = alpha * ne
   newNi = alpha * ni
   newNz = alpha * nz
   newPtot = newTe * newNe + newTi * (newNi + newNz)
   plt.title('alpha = '+str(alpha))
   plt.plot(rhot,ne*te**-1.5,label='ne*Te**-1.5 old')
   plt.plot(rhot,newNe*newTe**-1.5,label='ne*Te**-1.5 new')
   plt.plot(rhot,target_factor*ne*te**-1.5,'--',color='black',label='target')
   ax = plt.axis()
   plt.axis([0.9,1.0,0.0,ax[3]])
   plt.legend()
   plt.show()
   
# alpha times omte & omti
if mode == 'eta':
    midPedIndex = np.argmin(abs(rhot - rhotMidPed))
    teMidPed = te[midPedIndex]
    tiMidPed = ti[midPedIndex]
    print 'rhot =', rhotMidPed
    print 'te =', teMidPed
    print 'ti =', tiMidPed

    newTe = teMidPed*np.power(te/teMidPed,alpha)
    newTi = tiMidPed*np.power(ti/tiMidPed,alpha)

    Ptemp = ne * newTe + newTi * (ni + nz)
    #new density profile to keep total pressure the same
    newNe = ne*Ptot/Ptemp
    newNi = ni*Ptot/Ptemp
    newNz = nz*Ptot/Ptemp
    
    #total pressure with new density and tempeturate profiles
    newPtot = newNe * newTe + newTi * (newNi + newNz)

    etae0 = ne/te*fd_d1_o4_uneven(te,rhot)/fd_d1_o4_uneven(ne,rhot)
    newetae = newNe/newTe*fd_d1_o4_uneven(newTe,rhot)/fd_d1_o4_uneven(newNe,rhot)
    plt.title('alpha='+str(alpha))
    plt.plot(rhot,etae0,label='etae old')
    plt.plot(rhot,target_factor*etae0,'--',color='black',label='target')
    plt.plot(rhot,newetae,label='etae new')
    ax = plt.axis()
    plt.axis([0.9,1.0,0.0,6])
    plt.legend()
    plt.show()

if mode == 'omte':
    midPedIndex = np.argmin(abs(rhot - rhotMidPed))
    teMidPed = te[midPedIndex]
    print 'rhot =', rhotMidPed
    print 'te =', teMidPed

    newTe = teMidPed*np.power(te/teMidPed,alpha)
    newTi = ti

    newNe = ne
    newNi = ni
    newNz = nz
    
    #total pressure with new density and tempeturate profiles
    newPtot = newNe * newTe + newTi * (newNi + newNz)

    newomte = -1.0/newTe*fd_d1_o4_uneven(newTe,rhot)
    omte0 = -1.0/te*fd_d1_o4_uneven(te,rhot)
    #plt.plot(omte0)
    #plt.show()

    plt.title('alpha='+str(alpha))
    plt.plot(rhot,omte0,label='omte old')
    plt.plot(rhot,target_factor*omte0,'--',color='black',label='target')
    plt.plot(rhot,newomte,label='omte new')
    ax = plt.axis()
    print "ax",ax
    plt.axis([0.9,1.0,0.0,ax[3]])
    plt.legend()
    plt.show()

if mode == 'etaTe':
    midPedIndex = np.argmin(abs(rhot - rhotMidPed))
    teMidPed = te[midPedIndex]
    tiMidPed = ti[midPedIndex]
    print 'rhot =', rhotMidPed
    print 'te =', teMidPed
    print 'ti =', tiMidPed

    newTe = teMidPed*np.power(te/teMidPed,alpha)
    #new Ti profile to keep total pressure the same
    newTi = (te*ne+ti*(nz+ni) - newTe*ne)/(ni+nz)

    newNe = ne
    newNi = ni
    newNz = nz
    
    #total pressure with new density and tempeturate profiles
    newPtot = newNe * newTe + newTi * (newNi + newNz)

    etae0 = ne/te*fd_d1_o4_uneven(te,rhot)/fd_d1_o4_uneven(ne,rhot)
    newetae = newNe/newTe*fd_d1_o4_uneven(newTe,rhot)/fd_d1_o4_uneven(newNe,rhot)
    plt.title('alpha='+str(alpha))
    plt.plot(rhot,etae0,label='etae old')
    plt.plot(rhot,target_factor*etae0,'--',color='black',label='target')
    plt.plot(rhot,newetae,label='etae new')
    ax = plt.axis()
    plt.axis([0.9,1.0,0.0,6])
    plt.legend()
    plt.show()

if mode == 'omnz':
    Z = float(raw_input('Enter Z of impurity:\n'))
    print "Using Z= ",Z
    midPedIndex = np.argmin(abs(rhot - rhotMidPed))
    nzMidPed = nz[midPedIndex]

    newNz = nzMidPed*np.power(nz/nzMidPed,alpha)
    #Must satisfy quasineutrality and constant pressure
    newNe = (Ptot - ti*newNz*(1-Z))/(ti+te)
    newNi = newNe-Z*newNz

    newTe = te
    newTi = ti
    
    #total pressure with new density and tempeturate profiles
    newPtot = newNe * newTe + newTi * (newNi + newNz)

    omnz0 = fd_d1_o4_uneven(nz,rhot)/nz
    newomnz = fd_d1_o4_uneven(newNz,rhot)/newNz
    plt.title('alpha='+str(alpha))
    plt.plot(rhot,omnz0,label='omnz old')
    plt.plot(rhot,target_factor*omnz0,'--',color='black',label='target')
    plt.plot(rhot,newomnz,label='omnz new')
    ax = plt.axis()
    plt.axis([0.9,1.0,ax[2],ax[3]])
    plt.legend()
    plt.show()


if mode == 'TiTe':

    newTe = te*alpha
    #new Ti profile to keep total pressure the same
    newTi = (te*ne+ti*(nz+ni) - newTe*ne)/(ni+nz)

    newNe = ne
    newNi = ni
    newNz = nz
    
    #total pressure with new density and tempeturate profiles
    newPtot = newNe * newTe + newTi * (newNi + newNz)

    plt.title('alpha='+str(alpha))
    plt.plot(rhot,te/ti,label='Te/Ti old')
    plt.plot(rhot,target_factor*te/ti,'--',color='black',label='target')
    plt.plot(rhot,newTe/newTi,label='Te/Ti new')
    ax = plt.axis()
    plt.axis([0.9,1.0,0.0,ax[3]])
    plt.legend()
    plt.show()

if 1 == 1:
    plt.plot(rhot,ne,label='ne')
    plt.plot(rhot,newNe,label='new ne')
    plt.legend()
    plt.show()

    plt.plot(rhot,te,label='te')
    plt.plot(rhot,newTe,label='new te')
    plt.legend()
    plt.show()

    plt.plot(rhot,ti,label='ti')
    plt.plot(rhot,newTi,label='new ti')
    plt.legend()
    plt.show()

    plt.plot(rhot,Ptot,label='total P')
    plt.plot(rhot,newPtot,label='new total P')
    plt.legend()
    plt.show()

if 1 == 1:
    time_str = '9999'
    output_iterdb(rhot,rhop,newNe*1.E-19,newTe*1.E-3,newNi*1.E-19,newTi*1.E-3,file_out_base+'_alpha'+str(alpha)+'_'+mode,base_number,time_str,vrot=omega_tor,nimp=newNz*1.E-19)
    f=open('gene_profiles_i'+file_out_base+'_alpha'+str(alpha)+'_'+mode,'w')
    np.savetxt(f,np.column_stack((rhot,rhop,newTi*1.0e-3,newNi*1.0e-19)))
    f.close()
    f=open('gene_profiles_e'+file_out_base+'_alpha'+str(alpha)+'_'+mode,'w')
    np.savetxt(f,np.column_stack((rhot,rhop,newTe*1.0e-3,newNe*1.0e-19)))
    f.close()
    f=open('gene_profiles_z'+file_out_base+'_alpha'+str(alpha)+'_'+mode,'w')
    np.savetxt(f,np.column_stack((rhot,rhop,newTi*1.0e-3,newNz*1.0e-19)))
    f.close()
