import numpy as np
import matplotlib.pyplot as plt
from fastran_data import *
from os.path import exists
from finite_differences import *
from interp import *
import os
import scipy.integrate

#efit_name = 'g101473.05700'
#efit_name = 'g081481.05700'
rhot_omt = 0.3
sparse_factor = 1
plot_num = 500
plot_final = True
nrho_q = 1000
nrho_transport = 50
output_gene = True
output_gene_time = 100
output_restart = True
output_restart_time = 100
show_plots = True

rho_transport = np.linspace(0,1,nrho_transport+1)
e = 1.60218e-19
mu0 = 1.25664e-6
isd = read_instate('input/instate')
aminor = isd['aminor']

if_dict = read_infastran('input/infastran')
dt = if_dict['dt']

trnum = -1
for i in range(30):
    if exists('simulation_results/0/components/fastran_tr_fastran_'+str(i)+'/fastran.nc'):
        trnum = i
        break
if trnum == -1:
    print("ERROR: fastran_tr* not found")

path_end_tr = '/components/fastran_tr_fastran_'+str(trnum)+'/fastran.nc'
#path_end_eq = '/components/fastran_eq_efit_'+str(eqnum)+'/'+efit_name
fdict = get_profiles_ffile('simulation_results/'+'0'+path_end_tr)
nrho = len(fdict['rhot'][:])
print("nrho",nrho)
i=0
file_exists = exists('simulation_results/'+str(i)+path_end_tr)
while file_exists:
    i+=1
    file_exists = exists('simulation_results/'+str(i)+path_end_tr)
nt = i

#file_exists = exists('simulation_results/'+str(0)+path_end_eq)
#if file_exists:
#    os.system('calc_shat_from_efit.py '+efit_name)

ne = np.empty((nt,nrho))
omne = np.empty((nt,nrho))
ni = np.empty((nt,nrho))
omni = np.empty((nt,nrho))
nz = np.empty((nt,nrho))
omnz = np.empty((nt,nrho))
te = np.empty((nt,nrho))
omte = np.empty((nt,nrho))
ti = np.empty((nt,nrho))
omti = np.empty((nt,nrho))
zeff = np.empty((nt,nrho))
fluxe_exp = np.empty(nrho)
fluxi_exp = np.empty(nrho)
fluxe = np.empty((nt,nrho))
fluxi = np.empty((nt,nrho))
fluxe_hr = np.empty((nt,nrho))
fluxi_hr = np.empty((nt,nrho))
sion = np.empty((nt,nrho))
chie = np.empty((nt,nrho))
chii = np.empty((nt,nrho))
chie_exp = np.empty((nt,nrho))
chii_exp = np.empty((nt,nrho))
chie0 = np.empty((nt,nrho_transport+1))
chii0 = np.empty((nt,nrho_transport+1))
chin0 = np.empty((nt,nrho_transport+1))
De = np.empty((nt,nrho))
Di = np.empty((nt,nrho))
qprof = np.empty((nt,nrho))
shatprof = np.empty((nt,nrho))
R = np.empty((nt))
Ip = np.empty((nt))
sn = np.empty((nt))
time = np.empty((nt))
B0 = np.empty((nt))
omti0 = np.empty((nt))
shat0 = np.empty((nt))
beta = np.empty((nt,nrho))
#betap = np.empty((nt,nrho))
#betan_loc = np.empty((nt,nrho))
betan = np.empty((nt))
alpha = np.empty((nt,nrho))
j_tot = np.empty((nt,nrho))
j_bs = np.empty((nt,nrho))
j_oh = np.empty((nt,nrho))
j_nb = np.empty((nt,nrho))
j_rf = np.empty((nt,nrho))
pe_nb = np.empty((nt,nrho))
pe_rf = np.empty((nt,nrho))
p_rad = np.empty((nt,nrho))
p_ohm = np.empty((nt,nrho))
p_ei = np.empty((nt,nrho))
pe_fus = np.empty((nt,nrho))
pe_ionization = np.empty((nt,nrho))
pe_nb = np.empty((nt,nrho))
pi_nb = np.empty((nt,nrho))
pi_rf = np.empty((nt,nrho))
pi_fus = np.empty((nt,nrho))
pi_cx = np.empty((nt,nrho))
pi_ionization = np.empty((nt,nrho))
pe = np.empty((nt,nrho))
pi = np.empty((nt,nrho))
taue = np.empty((nt))
taui = np.empty((nt))
tauth = np.empty((nt))
tautot = np.empty((nt))
we = np.empty((nt))
wi = np.empty((nt))
tau98 = np.empty((nt))

print("shape ne",np.shape(ne))
#def get_flux_transport(flux_tot,nr,nrt):
nth = int( nrho/ nrho_transport )
print("nth",nth)

for i in range(nt):
    fdict = get_profiles_ffile('simulation_results/'+str(i)+path_end_tr)
    time[i] = dt*i
    if i==0:
        rhot = fdict['rhot']
        irhot0 = np.argmin(abs(rhot-rhot_omt))
        zimp = fdict['zimp']
        ne0 = fdict['ne'][0,:]
        te0 = fdict['te'][0,:]
        ti0 = fdict['ti'][0,:]
        ni0 = fdict['ni'][0,:]
        fluxe_exp = fdict['fluxe_exp']
        fluxi_exp = fdict['fluxi_exp']
        chie_exp = fdict['chie_exp']
        chii_exp = fdict['chii_exp']

    taue[i] = fdict['taue']
    taui[i] = fdict['taui']
    tau98[i] = fdict['tau98']
    tauth[i] = fdict['tauth']
    tautot[i] = fdict['tautot']
    we[i] = fdict['we']
    wi[i] = fdict['wi']
    ne[i,:] = fdict['ne'][1,:]
    omne[i,:] = -fd_d1_o4_smoothend(ne[i,:],rhot)/ne[i,:]
    ni[i,:] = fdict['ni2'][1,:]
    omni[i,:] = -fd_d1_o4_smoothend(ni[i,:],rhot)/ni[i,:]
    nz[i,:] = fdict['nz'][1,:]
    omnz[i,:] = -fd_d1_o4_smoothend(nz[i,:],rhot)/nz[i,:]
    te[i,:] = fdict['te'][1,:]
    omte[i,:] = -fd_d1_o4_smoothend(te[i,:],rhot)/te[i,:]
    ti[i,:] = fdict['ti'][1,:]
    omti[i,:] = -fd_d1_o4_smoothend(ti[i,:],rhot)/ti[i,:]
    omti0[i] = omti[i,irhot0]
    zeff[i,:] = fdict['zeff'][1,:]
    qprof[i,:] = fdict['q']
    sion[i,:] = fdict['sion'][1,:]
    j_tot[i,:] = fdict['j_tot'][1,:]
    j_oh[i,:] = fdict['j_oh'][1,:]
    j_rf[i,:] = fdict['j_rf'][1,:]
    j_nb[i,:] = fdict['j_nb'][1,:]
    j_bs[i,:] = fdict['j_bs'][1,:]
    pe_nb[i,:] = fdict['pe_nb']
    pe_rf[i,:] = fdict['pe_rf']
    p_rad[i,:] = fdict['p_rad']
    p_ohm[i,:] = fdict['p_ohm']
    p_ei[i,:] = fdict['p_ei']
    pe_ionization[i,:] = fdict['pe_ionization']
    pi_nb[i,:] = fdict['pi_nb']
    pi_rf[i,:] = fdict['pi_rf']
    pi_fus[i,:] = fdict['pi_fus']
    pi_cx[i,:] = fdict['pi_cx']
    pi_ionization[i,:] = fdict['pi_ionization']
    pe[i,:] = fdict['pe']
    pi[i,:] = fdict['pi']

    print("r0",fdict['r0'])
    R[i] = fdict['r0']
    Ip[i] = fdict['ip']
    sn[i] = fdict['sn']
    #if i == nt-1:
    #    time = fdict['time']
    #    print("time",time)
    B0[i] = fdict['b0']
    #betan_loc[i,:] = fdict['betan_loc']
    betan[i] = fdict['betan']
    print('betan',betan[i])
    #beta[i,:] = betan_loc[i,:]*Ip[i]/aminor/B0[i]
    beta[i,:] = (ni[i,:]*1e19*ti[i,:]*e*1000 + ne[i,:]*1e19*te[i,:]*e*1000)*2*mu0/B0[i]**2
    alpha[i,:] = -fd_d1_o4_smoothend(beta[i,:],rhot)*qprof[i,:]**2*R[i]/aminor
    #betap[i,:] = (ni[i,:]*1e19*ti[i,:]*e*1000 + ne[i,:]*1e19*te[i,:]*e*1000)*2/(mu0*Ip[i]**2*1e6**2)*8*np.pi**2*aminor**2


    file_exists = exists('simulation_results/'+str(i)+path_end_tr)
    #eqpath = 'simulation_results/'+str(i)
    #os.system('cp '+eqpath+path_end_eq+' ./')
    #os.system('calc_shat_from_efit.py '+efit_name)
    #shat_data = np.genfromtxt('shat.dat')
    #if i==0:
    #    rhot0 = shat_data[:,0]
    #    q0 = np.empty((nt,len(rhot0)))
    #    shat0 = np.empty((nt,len(rhot0)))
    #shat0[i,:] = shat_data[:,2]
    #q0[i,:] = shat_data[:,1]
    shatprof[i,:] = rhot/qprof[i,:]*fd_d1_o4_smoothend(qprof[i,:],rhot)
    shat0[i] = shatprof[i,irhot0]
    fluxe[i,:] = fdict['fluxe']
    fluxi[i,:] = fdict['fluxi']
    chie0[i,:] = fdict['chie'][0::nth]
    chii0[i,:] = fdict['chii'][0::nth]
    chin0[i,:] = fdict['chin'][0::nth]
    #chie0[i,:] = interp(rho_transport,fdict['chie'][0::nth],rhot)
    #chii0[i,:] = interp(rho_transport,fdict['chii'][0::nth],rhot)
    #chin0[i,:] = interp(rho_transport,fdict['chin'][0::nth],rhot)
    i+=1

#plt.plot(rhot,fluxe_exp)
#plt.show()

if output_gene:
    itime = np.argmin(abs(time - output_gene_time))
    output_gene_time = time[itime]
    zeros = np.zeros(len(ne[itime,:]))
    f=open('gene_profiles_e_t'+str(output_gene_time),'w')
    f.write('# 1.rho_tor 2.zeros 3.Te(kev) 4.ne(10^19m^-3)\n#\n')
    np.savetxt(f,np.column_stack((rhot,zeros,te[itime,:],ne[itime,:])))
    f.close()
    f=open('gene_profiles_i_t'+str(output_gene_time),'w')
    f.write('# 1.rho_tor 2.zeros 3.Ti(kev) 4.ni(10^19m^-3)\n#\n')
    np.savetxt(f,np.column_stack((rhot,zeros,ti[itime,:],ni[itime,:])))
    f.close()
    f=open('gene_profiles_c_t'+str(output_gene_time),'w')
    f.write('# 1.rho_tor 2.zeros 3.TC(kev) 4.nC(10^19m^-3)\n#\n')
    np.savetxt(f,np.column_stack((rhot,zeros,ti[itime,:],nz[itime,:])))
    f.close()
    #plt.plot(rhot,ne[-1,:])
    #plt.plot(rhot,ni[-1,:] + 6*nz[-1,:])
    #plt.show()

if output_restart:
    itime = np.argmin(abs(time - output_restart_time))
    output_four_col(te[itime,:],'te_restart_'+str(output_restart_time)+'.dat')
    output_four_col(ne[itime,:],'ne_restart_'+str(output_restart_time)+'.dat')
    output_four_col(ti[itime,:],'ti_restart_'+str(output_restart_time)+'.dat')
    output_four_col(qprof[itime,:],'q_restart_'+str(output_restart_time)+'.dat')
    output_four_col(j_tot[itime,:],'j_tot_restart_'+str(output_restart_time)+'.dat')
    p_eq = te[itime,:]*e*1000*ne[itime,:]*1e19+ti[itime,:]*e*1000*ni[itime,:]*1e19+ti[itime,:]*e*1000*nz[itime,:]*1e19
    if show_plots:
        plt.plot(rhot,p_eq)
        plt.xlabel('rhot')
        plt.ylabel('p_eq')
        plt.title('p_eq for restart')
        plt.show()
    output_four_col(p_eq,'p_eq_restart.dat')

if show_plots:
    plt.plot(time,tautot,label = 'tot')
    plt.plot(time,taue,label='e')
    plt.plot(time,taui,label='i')
    plt.plot(time,tauth,label='Therm.')
    plt.plot(time,tau98,label='H98')
    plt.xlabel('t(s)')
    plt.ylabel('Confinement time (s)')
    plt.legend(loc='upper left')
    plt.title('ttot: '+str(tautot[-1])[0:8]+'(s), tH98: '+str(tau98[-1])[0:8]+'(s), ti: '+str(taui[-1])[0:8]+'(s)')
    plt.show()
f=open('tau.dat','w')
f.write('#1.time 2.tautot 3.tauH98 4.taue 5.taui 6.tauth\n')
np.savetxt(f,np.column_stack((time,tautot,tau98,taue,taui,tauth)))
f.close()

if show_plots:
   plt.plot(time,we,label = 'e')
   plt.plot(time,wi,label='i')
   plt.xlabel('t(s)')
   plt.ylabel('Stored energy (MJ)')
   plt.legend(loc='upper left')
   plt.title('elec.: '+str(we[-1])[0:8]+'(MJ), ions: '+str(wi[-1])[0:8]+'(MJ)')
   plt.show()

f=open('stored_energy.dat','w')
f.write('#1.time 2.we 3.wi \n')
np.savetxt(f,np.column_stack((time,we,wi)))
f.close()




print("time",time)
print("B0",B0)
print("Ip",Ip)
if show_plots:
    plt.plot(time,Ip)
    plt.xlabel('t(s)')
    plt.ylabel('Ip(MA)')
    plt.savefig('Ip.pdf')
    plt.show()

    plt.plot(time,omti0)
    plt.xlabel('t(s)')
    plt.ylabel('omt(rhot0)')
    plt.show()
    np.savetxt('omt0.dat',np.column_stack((time,omti0)))


    plt.plot(time,shat0)
    plt.xlabel('t(s)')
    plt.ylabel('shat(rhot0)')
    plt.show()

    plt.plot(time,sn)
    plt.xlabel('t(s)')
    plt.ylabel('sn(10^19/s?)')
    plt.show()

    plt.plot(rhot,pe[-1,:],label='Pe total')
    plt.plot(rhot,pi[-1,:],label='Pe total')
    plt.ylabel('P(MW/m^3)')
    plt.xlabel('rhot')
    plt.show()

np.savetxt('shat0.dat',np.column_stack((time,shat0)))
pi_nb_tot = scipy.integrate.simpson(pi_nb[-1,:]*4*np.pi**2*rhot*fdict['a0']*fdict['r0'],fdict['a0']*rhot)
pe_rf_tot = scipy.integrate.simpson(pe_rf[-1,:]*4*np.pi**2*rhot*fdict['a0']*fdict['r0'],fdict['a0']*rhot)

if show_plots:
    plt.figure(figsize=(8,4))
    plt.suptitle('Power')
    plt.subplot2grid((1,2),(0,0))
    plt.title('Electrons (final)')
    plt.plot(rhot,pe[-1,:],color = 'black',label='total')
    plt.plot(rhot,pe_nb[-1,:],color = 'red',label='NB')
    plt.plot(rhot,pe_rf[-1,:],'--',color = 'blue',label='RF ('+str(pe_rf_tot)[0:6]+'MW)?')
    plt.plot(rhot,p_ohm[-1,:],'--',color = 'magenta',label='Ohmic')
    plt.plot(rhot,p_rad[-1,:],'-.',color = 'green',label='Rad')
    plt.ylabel('P(MW/m^3)')
    plt.xlabel('rhot')
    plt.legend()
    plt.subplot2grid((1,2),(0,1))
    plt.title('Ions (final)')
    plt.plot(rhot,pi[-1,:],color = 'black',label='total')
    plt.plot(rhot,pi_nb[-1,:],color = 'red',label='NB ('+str(pi_nb_tot)[0:6]+'MW)?')
    plt.plot(rhot,pi_rf[-1,:],'--',color = 'blue',label='RF')
    plt.plot(rhot,p_ohm[-1,:],'--',color = 'magenta',label='Ohmic')
    plt.plot(rhot,p_rad[-1,:],'-.',color = 'green',label='Rad')
    plt.legend()
    plt.ylabel('P(MW/m^3)')
    plt.xlabel('rhot')
    plt.savefig('./power.pdf')
    plt.show()



    plt.figure(figsize=(8,4))
    plt.suptitle('Current')
    plt.subplot2grid((1,2),(0,0))
    plt.plot(rhot,j_tot[0,:],color = 'black',label='total')
    jtot_tot0 = scipy.integrate.simpson(j_tot[0,:]*2*np.pi*rhot*fdict['a0'],fdict['a0']*rhot)
    plt.title('First time step: '+str(jtot_tot0)[0:6])
    plt.plot(rhot,j_bs[0,:],color = 'red',label='BS')
    plt.plot(rhot,j_oh[0,:],'--',color = 'blue',label='Ohmic')
    plt.plot(rhot,j_rf[0,:],'--',color = 'magenta',label='RF')
    plt.plot(rhot,j_nb[0,:],'-.',color = 'green',label='NB')
    plt.ylabel('J(MA/m^2)')
    plt.xlabel('rhot')
    plt.legend()
    plt.subplot2grid((1,2),(0,1))
    jtot_tot_end = scipy.integrate.simpson(j_tot[-1,:]*2*np.pi*rhot*fdict['a0'],fdict['a0']*rhot)
    plt.title('Last time step: '+str(jtot_tot_end)[0:6])
    plt.plot(rhot,j_tot[-1,:],color = 'black',label='total')
    plt.plot(rhot,j_bs[-1,:],color = 'red',label='BS')
    plt.plot(rhot,j_oh[-1,:],'--',color = 'blue',label='Ohmic')
    plt.plot(rhot,j_rf[-1,:],'--',color = 'magenta',label='RF')
    plt.plot(rhot,j_nb[-1,:],'-.',color = 'green',label='NB')
    plt.ylabel('J(MA/m^2)')
    plt.xlabel('rhot')
    plt.savefig('./current.pdf')
    plt.show()

f = open('current.dat','w')
f.write('#Contributions to current density at final time step.\n')
f.write('#1.rhot 2.j_tot 3.j_bs 4.j_oh 5.j_rf 6.j_nb\n')
np.savetxt(f,np.column_stack((rhot,j_tot[-1,:],j_bs[-1,:],j_oh[-1,:],j_rf[-1,:],j_nb[-1,:])))
f.close()

if show_plots:
    for i in range(nt):
        siontot = scipy.integrate.simpson(sion[-1,:]*4*np.pi**2*rhot*fdict['a0']*fdict['r0'],fdict['a0']*rhot)
        plt.plot(rhot,sion[i,:])
        #plt.title('sion(10^19 m^-3/s) Total(estimate)(#/s): '+str(siontot)[0:6])
        plt.title('PS(10^19 m^-3/s): '+str(siontot)[0:6])
        plt.xlabel('rhot')
    plt.show()

    for i in range(nt):
        plt.plot(rhot,beta[i,:])
        #plt.plot(rhot,betan_loc[i,:])
        plt.title(r'$\beta$')
        plt.xlabel('rhot')
    plt.plot(rhot,beta[0,:],'--',linewidth=3,color = 'black',label='Initial')
    plt.plot(rhot,beta[-1,:],'--',linewidth=3,color = 'blue',label='Final')
    plt.legend()
    plt.show()

    for i in range(nt):
        plt.plot(rhot,alpha[i,:])
        plt.title(r'$\alpha$')
        plt.xlabel('rhot')
    plt.plot(rhot,alpha[0,:],'--',linewidth=3,color = 'black',label='Initial')
    plt.plot(rhot,alpha[-1,:],'--',linewidth=3,color = 'blue',label='Final')
    plt.legend()
    plt.show()





#for i in range(nt):
#    plt.plot(rhot,betap[i,:])
#    #plt.plot(rhot,betan_loc[i,:])
#    plt.title(r'$\beta_p$')
#    plt.xlabel('rhot')
#plt.show()

    for i in range(nt):
        #plt.plot(rhot0,q0[i,:])
        plt.plot(rhot,qprof[i,:])
        plt.title('q')
        plt.xlabel('rhot')
    plt.plot(rhot,qprof[0,:],'--',linewidth=3,color = 'black',label='Initial')
    plt.plot(rhot,qprof[-1,:],'--',linewidth=3,color = 'blue',label='Final')
    plt.legend()
    plt.savefig('./safety_factor.pdf')
    plt.show()

    for i in range(nt):
        #plt.plot(rhot0,shat0[i,:])
        plt.plot(rhot,shatprof[i,:])
        plt.title('shat')
        plt.xlabel('rhot')
    plt.plot(rhot,shatprof[0,:],'--',linewidth=3,color = 'black',label='Initial')
    plt.plot(rhot,shatprof[-1,:],'--',linewidth=3,color = 'blue',label='Final')
    plt.legend()
    plt.savefig('./shear.pdf')
    plt.show()

f = open('q_shat.dat','w')
f.write('#q and shat at initial and final time steps\n')
f.write('#1.rhot 2.q(t=0) 3.shat(t=0) 4.q(t=final) 5.shat(t=final)\n')
np.savetxt(f,np.column_stack((rhot,qprof[0,:],shatprof[0,:],qprof[-1,:],shatprof[-1,:])))
f.close()
           
########Start here
#Interpolate onto rho grid
#translate n, T into SI units
#calculate chi, D
#Plot Q,chi, Gam, D on same plot
for i in range(nt):
    fluxe_hr[i,:] = interp(rho_transport,fluxe[i,:][0::nth],rhot)
    fluxi_hr[i,:] = interp(rho_transport,fluxi[i,:][0::nth],rhot)
    netemp = ne[i,:]*1e19
    Tetemp = te[i,:]*1000*e
    nitemp = ni[i,:]*1e19
    Titemp = ti[i,:]*1000*e
    chie[i,:] = 1000*fluxe_hr[i,:]/netemp/Tetemp/(omte[i,:]/aminor)
    chii[i,:] = 1000*fluxi_hr[i,:]/nitemp/Titemp/(omti[i,:]/aminor)



#plt.plot(rhot,fluxe_exp,color = 'black',label = 'Exp.')
#for i in range(nt):
#    if i%sparse_factor == 0 and plot_num > i:
#        plt.plot(rho_transport,fluxe[i,:][0::nth])
#        plt.plot(rhot,fluxe_hr[i,:])
#    if plot_final and i == nt-1:
#        plt.plot(rhot,fluxe_hr[i,:],color = 'orange', linewidth=3,label = 'Final')
#        plt.plot(rho_transport,fluxe[i,:][0::nth],'--',color = 'blue', linewidth=3,label = 'Final')
#plt.ylabel(r'$Q_e(kW/m^2)$',size=14)
#plt.legend()
#plt.xlabel(r'$\rho_{tor}$',size=14)
#plt.show()
#
#for i in range(nt):
#    if i%sparse_factor == 0 and plot_num > i:
#        plt.plot(rhot,chie[i,:])
#        plt.plot(rho_transport,chie0[i,:],color = 'red')
#        plt.plot(rho_transport,chin0[i,:],color = 'green')
#    if plot_final and i == nt-1:
#        plt.plot(rhot,chie[i,:],'--',color = 'blue', linewidth=3,label = 'Final')
#plt.ylabel(r'$\chi_e(m^2/s)$',size=14)
#plt.legend()
#plt.xlabel(r'$\rho_{tor}$',size=14)
#plt.show()
#
#
#plt.plot(rhot,fluxi_exp,color = 'black',label = 'Exp.')
#for i in range(nt):
#    if i%sparse_factor == 0 and plot_num > i:
#        plt.plot(rho_transport,fluxi[i,:][0::nth])
#    if plot_final and i == nt-1:
#        plt.plot(rho_transport,fluxi[i,:][0::nth],'--',color = 'blue', linewidth=3,label = 'Final')
#plt.ylabel(r'$Q_i(kW/m^2)$',size=14)
#plt.legend()
#plt.xlabel(r'$\rho_{tor}$',size=14)
#plt.show()
#
#
#for i in range(nt):
#    if i%sparse_factor == 0 and plot_num > i:
#        plt.plot(rhot,chii[i,:])
#        plt.plot(rhot,chii0[i,:],color = 'red')
#        plt.plot(rhot,chin0[i,:],color = 'green')
#    if plot_final and i == nt-1:
#        plt.plot(rhot,chii[i,:],'--',color = 'blue', linewidth=3,label = 'Final')
#plt.ylabel(r'$\chi_i(m^2/s)$',size=14)
#plt.legend()
#plt.xlabel(r'$\rho_{tor}$',size=14)
#plt.show()
#


#plt.plot(rhot,fluxi_exp,color = 'black',label = 'Exp.')
#for i in range(nt):
#    if i%sparse_factor == 0 and plot_num > i:
#        plt.plot(rhot,fluxi[i,:],'x')
#    if plot_final and i == nt-1:
#        plt.plot(rhot,fluxi[i,:],'+',color = 'blue', linewidth=3,label = 'Final')
#plt.ylabel('fluxi',size=14)
#plt.legend()
#plt.xlabel(r'$\rho_{tor}$',size=14)
#plt.show()

if show_plots:
    plt.figure(figsize=(8,7))
    plt.suptitle('Transport')

    plt.subplot2grid((2,2),(0,0))
    for i in range(nt):
        if i%sparse_factor == 0 and plot_num > i:
            plt.plot(rhot,fluxe_hr[i,:])
        if plot_final and i == nt-1:
            plt.plot(rhot,fluxe_hr[i,:],'--',color = 'blue', linewidth=3,label = 'Final')
    plt.plot(rhot,fluxe_exp,color = 'black', linewidth=3,label = 'Exp.')
    plt.ylabel(r'$Q_e(kWm^{-2})$',size=14)
    plt.legend()
    plt.xlabel(r'$\rho_{tor}$',size=14)

    plt.subplot2grid((2,2),(0,1))
    for i in range(nt):
        if i%sparse_factor == 0 and plot_num > i:
            plt.plot(rhot,fluxi_hr[i,:])
        if plot_final and i == nt-1:
            plt.plot(rhot,fluxi_hr[i,:],'--',color = 'blue', linewidth=3,label = 'Final')
    plt.plot(rhot,fluxi_exp,color = 'black', linewidth=3,label = 'Exp.')
    plt.ylabel(r'$Q_i(kWm^{-2})$',size=14)
    plt.legend()
    plt.xlabel(r'$\rho_{tor}$',size=14)

    plt.subplot2grid((2,2),(1,0))
    for i in range(nt):
        if i%sparse_factor == 0 and plot_num > i:
            if i ==0:
                plt.plot(rho_transport,chie0[i,:],label = 'Start')
            else:
                plt.plot(rho_transport,chie0[i,:])
    plt.plot(rho_transport,chie0[i,:],'--',color = 'blue',label = 'Final',linewidth=3)
    plt.plot(rhot,chie_exp,color = 'black',label = 'Exp.')
    plt.ylabel(r'$\chi_e(m^2/s)$',size=14)
    plt.xlabel(r'$\rho_{tor}$',size=14)
    plt.savefig('./transporte.pdf')
    plt.legend()

    plt.subplot2grid((2,2),(1,1))
    for i in range(nt):
        if i%sparse_factor == 0 and plot_num > i:
            if i ==0:
                plt.plot(rho_transport,chii0[i,:],label = 'Start')
            else:
                plt.plot(rho_transport,chii0[i,:])
    plt.plot(rho_transport,chii0[i,:],'--',color = 'blue',label = 'Final',linewidth=3)
    plt.plot(rhot,chii_exp,color = 'black',label = 'Exp.')
    plt.ylabel(r'$\chi_i(m^2/s)$',size=14)
    plt.xlabel(r'$\rho_{tor}$',size=14)
    plt.legend()
    plt.tight_layout()
    plt.savefig('./transporti.pdf')
    plt.show()


if 1==1:
    plt.figure(figsize=(8,4))
    plt.suptitle('Particle Transport')

    plt.subplot2grid((1,2),(0,0))
    for i in range(nt):
        if i%sparse_factor == 0 and plot_num > i:
            plt.plot(rho_transport,chin0[i,:])
        if plot_final and i == nt-1:
            plt.plot(rho_transport,chin0[i,:],'--',color = 'blue', linewidth=3,label = 'Final')
    #plt.plot(rho_transport,chin0[i,:],'--',color = 'blue', linewidth=3,label = 'Final')
    plt.ylabel(r'$D_e(m^2/s)$',size=14)
    plt.legend()
    plt.xlabel(r'$\rho_{tor}$',size=14)

    plt.subplot2grid((1,2),(0,1))
    for i in range(nt):
        if i%sparse_factor == 0 and plot_num > i:
            plt.plot(rho_transport,chin0[i,:]/chie0[i,:])
        if plot_final and i == nt-1:
            plt.plot(rho_transport,chin0[i,:]/chie0[i,:],'--',color = 'blue', linewidth=3,label = 'Final')
    #plt.plot(rhot,fluxi_exp,color = 'black', linewidth=3,label = 'Final')
    plt.ylabel(r'$D_e/\chi_e$',size=14)
    plt.legend()
    plt.xlabel(r'$\rho_{tor}$',size=14)
    plt.tight_layout()
    plt.savefig('./particle_transport.pdf')
    plt.show()


    plt.figure(figsize=(8,7))
    plt.suptitle('Electrons')
    plt.subplot2grid((2,2),(0,0))
    for i in range(nt):
        if i%sparse_factor == 0 and plot_num > i:
            plt.plot(rhot,ne[i,:])
        if plot_final and i == nt-1:
            plt.plot(rhot,ne[i,:],'--',color = 'blue', linewidth=3,label = 'Final')

    plt.plot(rhot,ne0,color = 'red',linewidth=2)
    plt.plot(isd['rhot'],isd['ne'],'--',color = 'black',linewidth=3,label='Initial')
    plt.ylabel(r'$n_e(10^{19}m^{-3})$',size=14)
    plt.legend()
    plt.xlabel(r'$\rho_{tor}$',size=14)

    plt.subplot2grid((2,2),(0,1))
    for i in range(nt):
        if i%sparse_factor == 0 and plot_num > i:
            plt.plot(rhot,te[i,:])
        if plot_final and i == nt-1:
            plt.plot(rhot,te[i,:],'--',color = 'blue', linewidth=2,label = 'Final')

    plt.plot(rhot,te0,color = 'red',linewidth=2)
    plt.plot(isd['rhot'],isd['te'],'--',color = 'black',linewidth=2,label='Initial')
    plt.ylabel(r'$T_e(keV)$',size=14)
    plt.legend()
    plt.xlabel(r'$\rho_{tor}$',size=14)

    plt.subplot2grid((2,2),(1,0))
    for i in range(nt):
        if i%sparse_factor == 0 and plot_num > i:
            if i ==0:
                plt.plot(rhot,omne[i,:],label='n')
                plt.plot(rhot,omte[i,:],'--',label='T')
            else:
                plt.plot(rhot,omne[i,:])
                plt.plot(rhot,omte[i,:],'--')
    plt.plot(rhot,omne[0,:],linewidth = 2,color = 'black',label='Initial')
    plt.plot(rhot,omte[0,:],'--',linewidth=2,color = 'black')
    plt.plot(rhot,omne[-1,:],linewidth = 2,color = 'blue',label='Final')
    plt.plot(rhot,omte[-1,:],'--',linewidth=2,color = 'blue')
    plt.ylabel(r'$a/L_{Te},a/L_{ne}$',size=14)
    plt.xlabel(r'$\rho_{tor}$',size=14)
    plt.legend()

    plt.subplot2grid((2,2),(1,1))
    for i in range(nt):
        if i%sparse_factor == 0 and plot_num > i:
            plt.plot(rhot,omte[i,:]/omne[i,:])
    plt.ylabel(r'$\eta_e$',size=14)
    plt.xlabel(r'$\rho_{tor}$',size=14)
    plt.axis((0,1,-10,10))
    plt.tight_layout()
    plt.savefig('./profiles_e.pdf')
    plt.show()




    plt.figure(figsize=(8,7))
    plt.suptitle('Main Ions')

    plt.subplot2grid((2,2),(0,0))
    for i in range(nt):
        if i%sparse_factor == 0 and plot_num > i:
            plt.plot(rhot,ni[i,:])
        if plot_final and i == nt-1:
            plt.plot(rhot,ni[i,:],'--',color = 'blue', linewidth=2,label = 'Final')
    plt.ylabel(r'$n_i(10^{19}m^{-3})$',size=14)
    plt.plot(rhot,ni[0,:],'--',color = 'black',linewidth=2,label='Initial')
    #plt.plot(isd['rhot'],isd['ni'],color = 'black',linewidth=2)
    plt.xlabel(r'$\rho_{tor}$',size=14)

    plt.subplot2grid((2,2),(0,1))
    for i in range(nt):
        if i%sparse_factor == 0 and plot_num > i:
            plt.plot(rhot,ti[i,:])
        if plot_final and i == nt-1:
            plt.plot(rhot,ti[i,:],'--',color = 'blue', linewidth=2,label = 'Final')
    plt.plot(rhot,ti0,color = 'red',linewidth=2)
    plt.plot(isd['rhot'],isd['ti'],'--',color = 'black',linewidth=2,label='Initial')
    plt.legend()
    plt.ylabel(r'$T_i(keV)$',size=14)
    plt.xlabel(r'$\rho_{tor}$',size=14)

    plt.subplot2grid((2,2),(1,0))
    for i in range(nt):
        if i%sparse_factor == 0 and plot_num > i:
            if i ==0:
                plt.plot(rhot,omti[i,:],'--',label='T')
                plt.plot(rhot,omni[i,:],label='n')
            else:
                plt.plot(rhot,omti[i,:],'--')
                plt.plot(rhot,omni[i,:])
    plt.plot(rhot,omni[0,:],'-',linewidth = 2,color = 'black')
    plt.plot(rhot,omti[0,:],'--',linewidth=2,color = 'black',label='Initial')
    plt.plot(rhot,omni[-1,:],'-',linewidth = 2,color = 'blue')
    plt.plot(rhot,omti[-1,:],'--',linewidth=2,color = 'blue',label='Final')
    plt.ylabel(r'$a/L_{Ti},a/L_{ni}$',size=14)
    plt.xlabel(r'$\rho_{tor}$',size=14)
    plt.legend()

    plt.subplot2grid((2,2),(1,1))
    for i in range(nt):
        if i%sparse_factor == 0 and plot_num > i:
            plt.plot(rhot,omti[i,:]/omni[i,:])
    plt.ylabel(r'$\eta_i$',size=14)
    plt.xlabel(r'$\rho_{tor}$',size=14)
    plt.axis((0,1,-10,10))
    plt.tight_layout()
    plt.savefig('./profiles_i.pdf')
    plt.show()




    plt.figure(figsize=(8,7))
    plt.suptitle('Impurities (Z ='+str(zimp)+')')

    plt.subplot2grid((2,2),(0,0))
    for i in range(nt):
        if i%sparse_factor == 0 and plot_num > i:
            plt.plot(rhot,nz[i,:])
    plt.ylabel(r'$n_z(10^{19}m^{-3})$',size=14)
    plt.xlabel(r'$\rho_{tor}$',size=14)

    plt.subplot2grid((2,2),(0,1))
    for i in range(nt):
        if i%sparse_factor == 0 and plot_num > i:
            plt.plot(rhot,ti[i,:])
    plt.ylabel(r'$T_i(keV)=T_z$',size=14)
    plt.xlabel(r'$\rho_{tor}$',size=14)

    plt.subplot2grid((2,2),(1,0))
    for i in range(nt):
        if i%sparse_factor == 0 and plot_num > i:
            if i ==0:
                plt.plot(rhot,omnz[i,:],label='n')
                plt.plot(rhot,omti[i,:],'--',label='Ti')
            else:
                plt.plot(rhot,omnz[i,:])
                plt.plot(rhot,omti[i,:],'--')
    plt.ylabel(r'$a/L_{Ti},a/L_{nz}$',size=14)
    plt.xlabel(r'$\rho_{tor}$',size=14)
    plt.legend()

    plt.subplot2grid((2,2),(1,1))

    for i in range(nt):
        if i%sparse_factor == 0 and plot_num > i:
            plt.plot(rhot,zeff[i,:])
    plt.ylabel(r'$Z_{eff}$',size=14)
    plt.xlabel(r'$\rho_{tor}$',size=14)

    plt.tight_layout()
    plt.show()






