from read_pfile import *
from read_EFIT import *
from finite_differences_x import *
import matplotlib.pyplot as plt
from write_iterdb import *
import sys

def comparePlot(x1, y1, x2, y2, label1, label2, xl, yl, lc = 1):
    plt.plot(x1, y1, label = label1)
    plt.plot(x2, y2, label = label2)
    plt.xlabel(xl)
    plt.ylabel(yl)
    plt.legend(loc = lc)
    plt.show()

####################################### inputs ###########################################

impurityCharge = 6.
efit_file_name = sys.argv[1]
p_file_name = sys.argv[2]

###########################################################################################

psi0, ne0, te0, ni0, ti0, nz0, er0 = read_pfile(p_file_name,impurityCharge,add_impurity=True)

zeff = (ni0 + nz0 * impurityCharge**2) / ne0 

EFITdict = read_EFIT(efit_file_name)
print(list(EFITdict.keys()))

sepInd = np.argmin(abs(EFITdict['psipn'] - 1.))
print('index at psipn = 1 is ', sepInd)
Rsep = EFITdict['R'][sepInd]
print('major R(m) at psipn = 1 is ', Rsep)
print('major R(m) at index = 1 is ', EFITdict['R'][0])

# construct R grid with uniform spacing 
# uniform spacing because first_derivative requires so
# find pressure, temperature, density values on uniform R grid

uni_R = np.linspace(EFITdict['R'][0],Rsep,EFITdict['nw']*10)
psip_uniR = interp(EFITdict['R'], EFITdict['psipn'], uni_R)
rhot_uniR = interp(EFITdict['R'], EFITdict['rhotn'], uni_R)

rhot0 = interp(EFITdict['psipn'], EFITdict['rhotn'], psi0)
pi0 = ni0 * ti0
pi_uniR = interp(rhot0,pi0,rhot_uniR)
ni_uniR = interp(rhot0,ni0,rhot_uniR)
ti_uniR = interp(rhot0,ti0,rhot_uniR)
pe0 = ne0 * te0
pe_uniR = interp(rhot0,pe0,rhot_uniR)
ne_uniR = interp(rhot0,ne0,rhot_uniR)
te_uniR = interp(rhot0,te0,rhot_uniR)

# compute grad P_i / n_i / e, grad P_e / n_e / e 
gradPioverNe = first_derivative(pi_uniR,uni_R)/ni_uniR 
gradPeoverNe = first_derivative(pe_uniR,uni_R)/ne_uniR 

# plot profiles in pedestal
if 1 == 1:
    plt.plot(rhot0,ne0,label='ne (10^20 m^-3)')
    plt.plot(rhot0,te0,label='te (KeV)')
    plt.plot(rhot0,ni0,label='ni (10^20 m^-3)')
    plt.plot(rhot0,nz0,label='nz (10^20 m^-3)')
    plt.plot(rhot0,ti0,label='ti (KeV)')
    plt.axis([0.85,1.,0.,3.])
    plt.xlabel('rhot')
    plt.legend()
    plt.show()

# construct rho_tor grid with uniform spacing 
# find pressure, temperature, density values on uniform rho_tor grid

uni_rhot = np.linspace(min(rhot0),max(rhot0),len(rhot0)*10.)
ti_u = interp(rhot0,ti0,uni_rhot)
te_u = interp(rhot0,te0,uni_rhot)
ne_u = interp(rhot0,ne0,uni_rhot)
ni_u = interp(rhot0,ni0,uni_rhot)
nz_u = interp(rhot0,nz0,uni_rhot)
p_u = (ni_u + nz_u) * ti_u + ne_u * te_u
if 1 == 1:
    comparePlot(uni_rhot, p_u * 1.E03 * 1.6E-19 * 1.E20, EFITdict['rhotn'], EFITdict['Pres'], 
                'Pressure ITERDB', 'Pressure EFIT', 'rhot', '', 1)

tprime_i = -first_derivative(ti_u,uni_rhot)/ti_u
tprime_e = -first_derivative(te_u,uni_rhot)/te_u
nprime_e = -first_derivative(ne_u,uni_rhot)/ne_u
nprime_i = -first_derivative(ni_u,uni_rhot)/ni_u
nprime_z = -first_derivative(nz_u,uni_rhot)/nz_u
eta_i = tprime_i / nprime_i
eta_e = tprime_e / nprime_e
eta_z = tprime_i / nprime_z

if 1 == 1:
    plt.plot(uni_rhot,nprime_e,label='nprime_e')
    plt.plot(uni_rhot,nprime_i,label='nprime_i')
    plt.plot(uni_rhot,nprime_z,label='nprime_z')
    plt.xlabel('rhot')
    plt.legend()
    plt.show()

if 1 == 1:
    comparePlot(uni_rhot, tprime_i, uni_rhot, tprime_e, 'tprime_i', \
                'tprime_e', 'rhot', '', 2)
    comparePlot(uni_rhot, eta_e, uni_rhot, eta_i, 'eta_e', \
                'eta_i', 'rhot', '', 1)
    
# convert from kV/m to V/m
Er_Vm = interp(rhot0,er0,uni_rhot)*1E3

# compare experimental Er with grad P_i / n_i / e
if 1 == 1:
#    comparePlot(rhot_uniR,gradPioverNe,rhot_uniR,gradPeoverNe,'grad Pe / ne /e',\
#                'grad Pi / ni / e', 'rhot', 'kV/m', 3)
    comparePlot(uni_rhot,Er_Vm/1E3, rhot_uniR, gradPioverNe, 'experimental Er',\
                'grad Pi / ni / e', 'rhot', 'kV/m', 3)

R_u = interp(EFITdict['rhotn'],EFITdict['R'],uni_rhot)
Bpol_u = interp(EFITdict['rhotn'],EFITdict['Bpol'],uni_rhot)

# add minus sign for consistency
omega_tor = - Er_Vm / (R_u * Bpol_u)

####################################### outputs ###########################################

# 'profiles_e/i/z' files are profile files GENE could run with
if 1 == 1:
    psi_u = interp(rhot0,psi0,uni_rhot)
    rhop_u = np.sqrt(psi_u)
    f = open('profiles_e','w')
    f.write('# 1.rho_tor 2.rho_pol 3.Te(kev) 4.ne(10^19m^-3)\n#\n')
    np.savetxt(f,np.column_stack((uni_rhot,rhop_u,te_u,ne_u*10.)))
    f.close()

    f = open('profiles_i','w')
    f.write('# 1.rho_tor 2.rho_pol 3.Ti(kev) 4.ni(10^19m^-3)\n#\n')
    np.savetxt(f,np.column_stack((uni_rhot,rhop_u,ti_u,ni_u*10.)))
    f.close()

    f = open('profiles_z','w')
    f.write('# 1.rho_tor 2.rho_pol 3.Ti(kev) 4.ni(10^19m^-3)\n#\n')
    np.savetxt(f,np.column_stack((uni_rhot,rhop_u,ti_u,nz_u*10.)))
    f.close()

# ITERDB file has profiles and vrot for ExB angular velocity
if 1 == 1:
    ############################ modify ############################
    file_out_base = 'DIIID' 
    base_number = '154406'
    time_str = '05170'
    ################################################################
    rhop=np.sqrt(psi0)
    psi_u = interp(rhot0,psi0,uni_rhot)
    rhop_u = np.sqrt(psi_u)
    # densities are multiplied by 10 here 
    # because output_iterdb() expects density in 10^19 m^-3
    output_iterdb(uni_rhot,rhop_u,ne_u*10.,te_u,ni_u*10.,ti_u,
                  file_out_base+base_number,base_number,time_str,
                  vrot=omega_tor,nimp=nz_u*10.)

