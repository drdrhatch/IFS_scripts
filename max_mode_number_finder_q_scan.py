from finite_differences import *
import matplotlib.pyplot as plt
from interp import *
import math
import csv
from read_EFIT import *
from read_iterdb_file import *
#Last edited by Max Curie 02/27/2020
#Supported by David R Hatch's script mtmDopplerFreqs.py


#**************Block for user******************************************
#**************Setting up*********************************************
n0_min=1         #minmum mode number (include) that finder will cover
n0_max=100       #maximum mode number (include) that finder will cover
omega_percent=3 #choose the omega within the top that percent defined in(0,100)
f_max=[200,450]        #upper bound of the frequency experimentally observed 
f_min=[100,350]        #lower bound of the frequency experimentally observed 
plot = 0         #set to 1 if you want to plot the result
report=1         #set to 1 if you want to export a csv report
q_max=1.1        #the maximum of q_scale for q scan, wher q=q_original * q_scale
q_min=0.9        #the minimum of q_scale for q scan
q_step=0.0001     #step of q_scale for q scan
#iterdb_file_name = 'DIIID175823.iterdb' #name of the iterdb file
#geomfile = 'g175823.04108_257x257'      #name of the magnetic geometry file
iterdb_file_name = 'DIIID174082.iterdb' #name of the iterdb file
geomfile = 'g174082.3000'      #name of the magnetic geometry file
x0_center=0.96                          #radial location where ky will be calculated.
Lref = 7.3893888417726761E-01           #minor radius of the device in meter
Bref_Gauss = 19774.                     #Magnetic field in Gauss at the center of the simulation
mref = 2.                               # mass of ion in proton mass
#**************End of Setting up*********************************************
#**************End of Block for user******************************************


#*************Loading the data******************************************
rhot0, te0, ti0, ne0, ni0, nz0, vrot0 = read_iterdb_file(iterdb_file_name)
EFITdict = read_EFIT(geomfile)
xgrid = EFITdict['psipn']
q = EFITdict['qpsi']

uni_rhot = np.linspace(min(rhot0),max(rhot0),len(rhot0)*10.)

te_u = interp(rhot0,te0,uni_rhot)
ne_u = interp(rhot0,ne0,uni_rhot)
vrot_u = interp(rhot0,vrot0,uni_rhot)
q      = interp(xgrid,q,uni_rhot)
tprime_e = -fd_d1_o4(te_u,uni_rhot)/te_u
nprime_e = -fd_d1_o4(ne_u,uni_rhot)/ne_u

#*************End of loading the data******************************************


#****************Start setting up ******************
center_index = np.argmin(abs(x0_center-uni_rhot)) 
#******************End setting up ****************

#****************Start scanning mode number*************
# this is for the q that has the mode one need 
ky_range=[]
n0_range=[]
m0_range=[]
f_range=[]
f_GENE_range=[]
x_range=[]
drive_range=[]
q_range=[]


#2 is selected q that has the band that is require
ky_range2=[]
n0_range2=[]
m0_range2=[]
f_range2=[]
f_GENE_range2=[]
x_range2=[]
drive_range2=[]
q_range2=[]

q_temp=q

for q_scale in np.arange(q_min,q_max,q_step):
    ky_range3=[]
    n0_range3=[]
    m0_range3=[]
    f_range3=[]
    f_GENE_range3=[]
    x_range3=[]
    drive_range3=[]
    q_range3=[]

    q=q_temp*q_scale
    q0      = q[center_index]
    ne = ne_u[center_index]
    te = te_u[center_index] #it is in eV
    Bref=float(Bref_Gauss)/10000
    m_SI = mref *1.6726*10**(-27)
    c  = 1
    qref = 1.6*10**(-19)
    nref = ne * 1.E19
    Tref = te * qref
    Cy0 = x0_center/q0
    cref = np.sqrt(Tref / m_SI)
    Omegaref = qref * Bref / m_SI / c
    rhoref = cref / Omegaref 
    for n0 in range(n0_min,n0_max+1):
        n0_TEMP=0

    #***********Calculating the ky********************************
    #From SI_Gauss_GENE_unit.py
        m0 = n0*q0
        ky=n0*q0*rhoref/(Lref*x0_center)
        kymin = ky
        n0_global = n0
        te_mid = te_u[center_index]
        kyGENE =kymin * (q/q0) * np.sqrt(te_u/te_mid) #Add the effect of the q varying
    #***Calculate omeage star********************************
    #from mtm_doppler
        omMTM = kyGENE*(tprime_e+nprime_e)
        gyroFreq = 9.79E3/np.sqrt(mref)*np.sqrt(te_u)/Lref
        mtmFreq = omMTM*gyroFreq/2./np.pi/1000.
        omegaDoppler = vrot_u*n0_global/2./np.pi/1E3
        omega=mtmFreq + omegaDoppler
        if plot==1:
            plt.clf()
            plt.title('mode number finder')
            plt.xlabel('r/a')
            plt.ylabel('frequency (kHz)') 
            plt.plot(uni_rhot,mtmFreq + omegaDoppler,label='Diamagnetic plus Doppler (MTM in lab frame)')
    #***End of Calculate omeage star**************************
    #***Find the range where the omega is within the top omega_percent%
        omega_range=[]
        range_ind=[]
        omega_max=np.max(omega)
        omega_min=omega_max*(100-omega_percent)/100
        for i in range(len(uni_rhot)):
            if omega[i] >= omega_min:
                omega_range.append(uni_rhot[i])
                range_ind.append(i)
        range_min=min(omega_range)
        ind_min  =min(range_ind)
        range_max=max(omega_range)
        ind_max  =max(range_ind)

    #***End of Find the range where the omega is within the top omega_percent%
    
    #*******Find the range within the frequency range************************
        range_ind2=[]
        for i in range(len(uni_rhot)):
            if omega[i] >= min(f_min)and omega[i] <= max(f_max):
                range_ind2.append(i)

    #*******Find the range within the frequency range************************
    
    #find the rational surfaces
    
    #From the plot_mode_structures.py
        qmin = np.min(q[ind_min:ind_max])
        qmax = np.max(q[ind_min:ind_max])
        mmin = math.ceil(qmin*n0)
        mmax = math.floor(qmax*n0)
        mnums = np.arange(mmin,mmax+1)
        qrats = mnums/float(n0)
    
        if ky <= 5:
    #***********End of Calculating the ky************************
            for i in range(len(qrats)):
                n = int(n0)
                m = int(mmin + i)
                ix = np.argmin(abs(q-qrats[i])) 
                if uni_rhot[ix] >= range_min and uni_rhot[ix] <= range_max:
                    if ix in range_ind2:
                        print('ky='+str(ky))
                        print('(m,n)='+str((m,n)))
                        temp_str=str((n,m))
                        if plot==1:
                            plt.axvline(uni_rhot[ix],color='red', label= temp_str)
                        n0_TEMP=n0_TEMP+1
                        ky_range.append(ky)
                        n0_range.append(n)
                        m0_range.append(m)
                        f_range.append(omega[ix])
                        f_GENE_range.append(omega[ix]*2*np.pi*1000/gyroFreq[ix])
                        x_range.append(uni_rhot[ix])
                        drive_range.append(omega[ix]/omega_max)
                        q_range.append(q_scale)

                        ky_range3.append(ky)
                        n0_range3.append(n)
                        m0_range3.append(m)
                        f_range3.append(omega[ix])
                        f_GENE_range3.append(omega[ix]*2*np.pi*1000/gyroFreq[ix])
                        x_range3.append(uni_rhot[ix])
                        drive_range3.append(omega[ix]/omega_max)
                        q_range3.append(q_scale)

    
                    else:
                        if plot==1:
                            plt.axvline(uni_rhot[ix],color='green')
                else:
                    if plot==1:
                        plt.axvline(uni_rhot[ix],color='green')
        else:
           for i in range(len(qrats)):
                n = int(n0)
                m = int(mmin + i)
                ix = np.argmin(abs(q-qrats[i]))
                if plot==1:
                    plt.axvline(uni_rhot[ix],color='green')

        if n0_TEMP > 0 and plot==1:
            plt.xlim(0.9,1)
            plt.ylim(0,max(mtmFreq + omegaDoppler))
            plt.legend()
            plt.savefig('mode_number_finder_n0='+str(n0)+'.png')
    print('**********************************************')
    print('**************Start of report*****************')
    if len(ky_range)==0:
        print('There is no unstabel MTM')
    else:
        print('ky range from '+str(min(ky_range))+' to '+str(max(ky_range)))
        print('n0 range from '+str(min(n0_range))+' to '+str(max(n0_range)))
    print('**********************************************')
    print('***************End of report******************')

#*************Select a band******************

    temp_1=1 #1 if this q is "good", 0 if this q is bad
    
    #if not all band has at least one "line", then it is "bad"
    num_in_band=np.zeros(len(f_max))
    for i in range(len(f_max)):
        for j in range(len(f_range3)):
            if f_range3[j]>= f_min[i] and f_range3[j]<= f_max[i]:
                num_in_band[i]=num_in_band[i]+1
        if num_in_band[i]==0:
            temp_1=0

    #if the f is outside the band then it is "bad"
    for i in range(len(f_range3)):
        for j in range(len(f_max)-1):
            if f_range3[i] > f_max[j] and f_range3[i] < f_min[j+1]:  #this q has a frequency outside of the band
                temp_1=0
                break

    if temp_1 == 1:
        for i in range(len(ky_range3)):
            ky_range2.append(ky_range3[i])
            n0_range2.append(n0_range3[i])
            m0_range2.append(m0_range3[i])
            f_range2.append(f_range3[i])
            f_GENE_range2.append(f_GENE_range3[i])
            x_range2.append(x_range3[i])
            drive_range2.append(drive_range3[i])
            q_range2.append(q_range3[i])
#*************End of select a band******************
            
    
if report==1:
    with open('mode_number_finder_report_all.csv','w') as csvfile:
        data = csv.writer(csvfile, delimiter=',')
        data.writerow(['q_scale        ','n ',       'm ',       'kymin          ','frequency(kHz)           ','location(r/a)            ','omega_GENE    ','Drive'])
        for i in range(len(n0_range)):
            data.writerow([q_range[i]   ,n0_range[i],m0_range[i],ky_range[i],      f_range[i],                  x_range[i],                f_GENE_range[i],  drive_range[i]])
    csvfile.close()
    



    with open('mode_number_finder_report_band.csv','w') as csvfile2:
        data2 = csv.writer(csvfile2, delimiter=',')
        data2.writerow(['q_scale        ', 'n ',       'm ',        'kymin          ','frequency(kHz)           ','location(r/a)            ','omega_GENE    ','Drive'])
        for i in range(len(n0_range2)):
            data2.writerow([q_range2[i]   ,n0_range2[i],m0_range2[i],ky_range2[i],   f_range2[i],                x_range2[i],                f_GENE_range2[i],drive_range2[i]])
    csvfile2.close()

print(q_range2)
    


