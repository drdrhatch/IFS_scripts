import numpy as np
from write_iterdb import *
from finite_differences import *
from interp import *

####### Modify ######
base_name = 'JET1_1_t14_w118_Zeff1.6Be_from_base'
base_num = ''
file_in_name = 'rbsProfs'
show_plots = False
test_plots = True
####### Modify ######

dummy = (raw_input('Warning! Modifying electron and ion density profiles to match pressure from rbsProfs (press any key): \n'))
zimp = float(raw_input('Enter Z for impurity species: \n'))
zeff = float(raw_input('Enter Zeff: \n'))

def calc_nfracs(zeff,zimp):
    nc_frac = (zeff-1)/(zimp**2 - zimp)
    ni_frac = 1-zimp*nc_frac

    print "Calculating impurity and main ion (z=1) fractions for:"
    print "Zeff = ",zeff
    print "Zimp = ",zimp
    print "nc/ne:",nc_frac
    print "ni/ne:",ni_frac

    zeff_test = ni_frac + zimp**2*nc_frac
    qn_test = ni_frac+zimp*nc_frac

    print "Zeff test:", zeff_test
    print "Quasineutrality test (should be 1.0):", qn_test
    return ni_frac,nc_frac

def calc_a():
   f=open('rbsProfs','r')
   rbs = f.read()
   f.close()
   rbs = rbs.split('\n')  
   a_factor = float(rbs[1].split()[3])
   print "a_factor",a_factor
   rbs = np.genfromtxt('rbsProfs')
   isep = np.argmin(abs(rbs[:,0]-1.0))
   print "rhotor[isep]",rbs[isep,0]   
   a = rbs[isep,22]
   print "a[isep]",a
   print "a_min_eff",a_factor*a
   return a_factor*a

data = np.genfromtxt(file_in_name)
a = calc_a()
#Making gene profiles:
#0:rho_tor
#1:psi
#2:pres(MKS)
#3:dens(10^20/m^3)
#4:Ti(kev)
#5:Te(kev)
#9: gamE_N   (Neoclassical)
#10: gamE0_N
#16: Er_MKS (V/m)
#17: Er0_MKS (V/m)
#21: R[m]  (major radius of center of flux surface?)
#23: q
#24: R_out [m] (Local major radius at outboard midplane)
#25: B_pol[T]
#26: B_tor[T]

rhot = data[:,0]
rhop = np.sqrt(data[:,1])
#convert to 10^19
dens = data[:,3]*10.0

ti = data[:,4]
te = data[:,5]

alpha1 = ti*((zimp**2-zimp)*zeff/(zeff-1) - zimp**2) + te*(zimp**2-zimp)/(zeff-1) + ti 
ptot = dens*(ti+te)
ne = ptot/alpha1*(zimp**2-zimp)/(zeff-1)  
nz = ptot/alpha1
ni = ptot/alpha1*((zimp**2-zimp)*zeff/(zeff-1) - zimp**2)

print "ion dilution:\n",ni[-1]/ne[-1]
plt.plot(rhot,ni/ne,label='ion dilution ni/ne')
plt.legend()
plt.show()

if test_plots:

    plt.plot(rhot,ne,label='ne')
    plt.plot(rhot,ni,label='ni')
    plt.plot(rhot,nz,label='nz')
    plt.plot(rhot,dens,label='n0',color = 'black',linewidth=2)
    plt.legend()
    plt.show()

    print "Pressure test:"
    plt.plot(ptot,label='ptot')
    plt.plot(ti*ni+ti*nz+te*ne,label='ptot test')
    plt.legend()
    plt.show()

    print "zeff test:"
    plt.plot(rhot,(ni+zimp**2*nz)/ne,label='zeff')
    plt.legend()
    plt.show()

    print "quasineutrality test:"
    plt.plot(rhot,ni+zimp*nz-ne,label='qn test (should be zero)')
    plt.legend()
    plt.show()


f=open('gene_profiles'+file_in_name[10:]+'_i','w')
f.write('# 1.rho_tor 2.rho_pol 3.Ti(kev) 4.ni(10^19m^-3)\n#\n')
np.savetxt(f,np.column_stack((rhot,rhop,ti,ni)))
f.close()

f=open('gene_profiles'+file_in_name[10:]+'_e','w')
f.write('# 1.rho_tor 2.rho_pol 3.Te(kev) 4.ne(10^19m^-3)\n#\n')
np.savetxt(f,np.column_stack((rhot,rhop,te,ne)))
f.close()

f=open('gene_profiles'+file_in_name[10:]+'_imp','w')
f.write('# 1.rho_tor 2.rho_pol 3.Ti(kev) 4.nimp(10^19m^-3)\n#\n')
np.savetxt(f,np.column_stack((rhot,rhop,te,nz)))
f.close()

omega_tor = data[:,16]/(data[:,25]*data[:,24])
if show_plots:
    plt.plot(rhot,omega_tor)
    plt.plot(data[:,0],data[:,27])
    plt.xlabel('rho_tor')
    plt.ylabel('omega_tor (rad/s)')
    plt.show()


#Calculte Hahm-Burrell
psi_new = np.arange(1000)/999.0
#omega_tor_new = interp(data[:,1],omega_tor,psi_new)
#### Hahm-Burrell ####
#### Hahm-Burrell ####
#### Hahm-Burrell ####
R_new = interp(data[:,1],data[:,24],psi_new)
Bpol_new = interp(data[:,1],data[:,25],psi_new)
Btor_new = interp(data[:,1],data[:,26],psi_new)
rhot_new = interp(data[:,1],data[:,0],psi_new)
Er_new = interp(data[:,1],data[:,16],psi_new)
omega_tor_new = Er_new/(Bpol_new*R_new)
omega_tor_prime = fd_d1_o4(omega_tor_new,psi_new)
gamma_HB = omega_tor_prime*(Bpol_new*R_new)**2/((Bpol_new**2+Btor_new**2)**0.5)
cs = np.sqrt(te*1000*1.6e-19/(2.0*1.67e-27))
cs_new = interp(data[:,1],cs,psi_new)
gamma_HB_norm = gamma_HB/(cs_new/a)
#np.savetxt('gamma_HB_'+base_name,np.column_stack((rhot_new,gamma_HB_norm,omega_tor_new)))
#plt.plot(rhot_new,gamma_HB_norm,label='gamma_HB')
#plt.plot(rhot,-data[:,9],label='-gamma_HB(Prashant)')
#plt.title('gamma_HB/(cs/a)')
#plt.legend(loc = 'upper left')
#plt.show()

output_iterdb(rhot[:],rhop[:],ne[:],te[:],ni[:],ti[:],base_name,'9999','0.0',vrot=omega_tor,nimp = nz)

#### Hahm-Burrell estimate ####
#### Hahm-Burrell estimate ####
#### Hahm-Burrell estimate ####
R_maj = data[:,24]
R_new = np.arange(1000)/999.0*(R_maj[-1]-R_maj[0])+R_maj[0]
#plt.plot(R_maj)
#plt.plot(R_new)
#plt.title('R old vs new')
#plt.show()
Btor_new = interp(R_maj,data[:,26],R_new)
Bpol_new = interp(R_maj,data[:,25],R_new)
psi_new = interp(R_maj,data[:,1],R_new)
Er_new = interp(R_maj,data[:,16],R_new)
rhot_new = interp(R_maj,data[:,0],R_new)
gamma_HB = fd_d1_o4(Er_new,R_new)/Btor_new
cs = np.sqrt(te*1000*1.6e-19/(2.0*1.67e-27))
cs_new = interp(R_maj,cs,R_new)
gamma_HB_norm = gamma_HB/(cs_new/a)
#np.savetxt('gamma_HB2_'+base_name,np.column_stack((rhot_new,gamma_HB_norm,omega_tor_new)))

#plt.plot(R_new,gamma_HB_norm,label='gamma_HB2')
#plt.plot(R_maj,data[:,9],label='gamma_HB(Prashant)')
#plt.title('gamma_HB/(cs/a)')
#plt.legend(loc = 'upper left')
#plt.show()

#### Hahm-Burrell using R ####
#### Hahm-Burrell using R ####
#### Hahm-Burrell using R ####

gamma_HB = fd_d1_o4(Er_new/Bpol_new/R_new,R_new)
gamma_HB = gamma_HB*Bpol_new*R_new/(Bpol_new**2+Btor_new**2)**0.5
gamma_HB_norm = gamma_HB/(cs_new/a)
#np.savetxt('gamma_HB3_'+base_name,np.column_stack((rhot_new,gamma_HB_norm,omega_tor_new)))
if show_plots:
    plt.plot(R_new,gamma_HB_norm,label='gamma_HB3')
    plt.plot(R_maj,-data[:,9],label='gamma_HB(Prashant)')
    plt.title('gamma_HB/(cs/a)')
    plt.legend(loc = 'upper left')
    plt.show()

### Calculate GENE shear rates ###
rhot_new = np.arange(1000)/999.0
q_new = interp(data[:,0],data[:,23],rhot_new)
Er_new = interp(data[:,0],data[:,16],rhot_new)
Bpol_new = interp(data[:,0],data[:,25],rhot_new)
R_new = interp(data[:,0],data[:,24],rhot_new)
omega_tor_new = Er_new/(R_new*Bpol_new)
omega_tor_new_prime = fd_d1_o4(omega_tor_new,rhot_new)
gamma_GENE = rhot_new/q_new*omega_tor_new_prime
cs = np.sqrt(data[:,5]*1000*1.6e-19/(2.0*1.67e-27))
cs_new = interp(data[:,0],cs,rhot_new)
gamma_GENE_norm = gamma_GENE/(cs_new/a)
f=open('gamma_GENE_'+base_name,'w')
f.write('# 1.rhot	2.gamma_GENE_norm	3.omega_tor\n')	
np.savetxt(f,np.column_stack((rhot_new,gamma_GENE_norm,omega_tor_new)))
f.close()








