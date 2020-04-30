from finite_differences import *
import matplotlib.pyplot as plt
from interp import *
import numpy as np
import math
import csv
from genetools import *
from omega_tool_max import omega_calc

#Last edited by Max Curie 04/27/2020
#Criteria: D_chi, typical frequency, Epar

#**********Start of input setup**********************************************************************
def species_order(suffix): 
    #*****Start of calcuation of diamagnetic frequency******
    paramfpath="parameters_"+str(suffix)
    geneparam=read_parameters(paramfpath)
    #Determine the order of the species
    
    if 'species3' in geneparam: 
        species=range(3)
        for i in range(3):
            if geneparam['species'+str(i+1)]['charge']==-1:
                species[i]='e' #it is eletron
            elif geneparam['species'+str(i+1)]['charge']==1:
                species[i]='i' #it is ion
            elif geneparam['species'+str(i+1)]['charge']>1:
                species[i]='z' #it is impurity
    elif 'species3' not in geneparam: 
        species=range(2)
        for i in range(2):
            if geneparam['species'+str(i+1)]['charge']==-1:
                species[i]='e' #it is eletron
            elif geneparam['species'+str(i+1)]['charge']==1:
                species[i]='i' #it is ion
    return species

def Epar(suffix):
    fieldfpath="field_"+str(suffix)
    fielddata=read_field(fieldfpath)
    epar=field_info(fielddata)['Epar_Cancellation']
    return epar

def f(suffix): #calculate the frequency
    #from genetools.py
    #momfpath_e=
    #momdata_e=read_mom(momfpath)
#******start of reading omega*******
    omegafpath="omega_"+str(suffix)
    if not os.path.isfile(omegafpath):
        print('No omega file, using omega scan')
        #from omega_tool_max2
        omega_calc("_"+str(suffix))
        #from genetools.py
        omegadata = read_omega(omegafpath)
    else:
        omegadata = read_omega(omegafpath)
    
    omega=omegadata['omega'][0]
    gamma=omegadata['gamma'][0]
    kymin=omegadata['kymin'][0]
#******End of reading omega*******

#*****Start of calcuation of diamagnetic frequency******
    paramfpath="parameters_"+str(suffix)
    geneparam=read_parameters(paramfpath)

    species=species_order(suffix)

    omn_e=geneparam['species'+str(species.index('e')+1)]['omn']
    omn_i=geneparam['species'+str(species.index('i')+1)]['omn']
    omt_e=geneparam['species'+str(species.index('e')+1)]['omt']
    omt_i=geneparam['species'+str(species.index('i')+1)]['omt']
    #print(kymin)
    fe = kymin*(omn_e+omt_e) #omega*e in cs/a
    fi = kymin*(omn_i+omt_i) #omega*i in cs/a
    #print("fe is "+str(fe))

    #nu_ei = geneparam['info']['nu_ei']
    nu_ei = 0

    if 'kx_center' in geneparam['box']:
        kx=geneparam['box']['kx_center']
    else:
        kx=0
    
    Cs = npy.sqrt(geneparam['units']['Tref']*1000.0*1.602e-19/geneparam['units']['mref']/1.6726e-27)
    omegaref = Cs/geneparam['units']['Lref']

    f_kHz=omega*omegaref/(2.*3.1415926*1000.)

    #doppler = vrot_u*n0_global/2./np.pi/1E3 #Doppler shift in kHz

    #*****End of calcuation of diamagnetic frequency******
    return fe, fi, omega, gamma, kx, kymin, f_kHz, nu_ei

def D_chi_3(suffix):

    #from genetools.py
    paramfpath="parameters_"+str(suffix)
    geneparam=read_parameters(paramfpath)

    Tref=geneparam['units']['Tref']
    nref=geneparam['units']['nref']

    species=species_order(suffix)
    Te=geneparam['species'+str(species.index('e')+1)]['temp']
    Ti=geneparam['species'+str(species.index('i')+1)]['temp']
    Tz=geneparam['species'+str(species.index('z')+1)]['temp']
    ne=geneparam['species'+str(species.index('e')+1)]['dens']
    ni=geneparam['species'+str(species.index('i')+1)]['dens']
    nz=geneparam['species'+str(species.index('z')+1)]['dens']
    omn_e=geneparam['species'+str(species.index('e')+1)]['omn']
    omn_i=geneparam['species'+str(species.index('i')+1)]['omn']
    omn_z=geneparam['species'+str(species.index('z')+1)]['omn']
    omt_e=geneparam['species'+str(species.index('e')+1)]['omt']
    omt_i=geneparam['species'+str(species.index('i')+1)]['omt']
    omt_z=geneparam['species'+str(species.index('z')+1)]['omt']
    
    #from genetools.py
    nrgfpath="nrg_"+str(suffix)
    nrgdata = read_nrg(nrgfpath)

    D_e=nrgdata[nrgfpath]['e']['PFluxes'][-1]+nrgdata[nrgfpath]['e']['PFluxem'][-1]
    D_i=nrgdata[nrgfpath]['i']['PFluxes'][-1]+nrgdata[nrgfpath]['i']['PFluxem'][-1]
    D_z=nrgdata[nrgfpath]['z']['PFluxes'][-1]+nrgdata[nrgfpath]['z']['PFluxem'][-1]
    D_e=D_e/ omn_e / ne
    D_i=D_i/ omn_i / ni
    D_z=D_z/ omn_z / nz    
    Qes_e=nrgdata[nrgfpath]['e']['HFluxes'][-1]-3./2.*Te*nrgdata[nrgfpath]['e']['PFluxes'][-1]
    Qes_i=nrgdata[nrgfpath]['i']['HFluxes'][-1]-3./2.*Ti*nrgdata[nrgfpath]['i']['PFluxes'][-1]
    Qes_z=nrgdata[nrgfpath]['z']['HFluxes'][-1]-3./2.*Tz*nrgdata[nrgfpath]['z']['PFluxes'][-1]
    Qem_e=nrgdata[nrgfpath]['e']['HFluxem'][-1]-3./2.*Te*nrgdata[nrgfpath]['e']['PFluxem'][-1]
    Qem_i=nrgdata[nrgfpath]['i']['HFluxem'][-1]-3./2.*Ti*nrgdata[nrgfpath]['i']['PFluxem'][-1]
    Qem_z=nrgdata[nrgfpath]['z']['HFluxem'][-1]-3./2.*Tz*nrgdata[nrgfpath]['z']['PFluxem'][-1]
    Q_e = (Qes_e+Qem_e)
    Q_i = (Qes_i+Qem_i)
    Q_z = (Qes_z+Qem_z)
    chi_e = (Qes_e+Qem_e) / omt_e / ne / Te
    chi_i = (Qes_i+Qem_i) / omt_i / ni / Ti
    chi_z = (Qes_z+Qem_z) / omt_z / nz / Tz

    return Qes_e,Qes_i,Qes_z,Qem_e,Qem_i,Qem_z,Q_e,Q_i,Q_z,D_e,D_i,D_z,chi_e,chi_i,chi_z

def D_chi_2(suffix):

    #from genetools.py
    paramfpath="parameters_"+str(suffix)
    geneparam=read_parameters(paramfpath)

    Tref=geneparam['units']['Tref']
    nref=geneparam['units']['nref']

    species=species_order(suffix)
    Te=geneparam['species'+str(species.index('e')+1)]['temp']
    Ti=geneparam['species'+str(species.index('i')+1)]['temp']
    ne=geneparam['species'+str(species.index('e')+1)]['dens']
    ni=geneparam['species'+str(species.index('i')+1)]['dens']
    omn_e=geneparam['species'+str(species.index('e')+1)]['omn']
    omn_i=geneparam['species'+str(species.index('i')+1)]['omn']
    omt_e=geneparam['species'+str(species.index('e')+1)]['omt']
    omt_i=geneparam['species'+str(species.index('i')+1)]['omt']
    
    #from genetools.py
    nrgfpath="nrg_"+str(suffix)
    nrgdata = read_nrg(nrgfpath)

    D_e=nrgdata[nrgfpath]['e']['PFluxes'][-1]+nrgdata[nrgfpath]['e']['PFluxem'][-1]
    D_i=nrgdata[nrgfpath]['i']['PFluxes'][-1]+nrgdata[nrgfpath]['i']['PFluxem'][-1]
    D_e=D_e/ omn_e / ne
    D_i=D_i/ omn_i / ni
    Qes_e=nrgdata[nrgfpath]['e']['HFluxes'][-1]-3./2.*Te*nrgdata[nrgfpath]['e']['PFluxes'][-1]
    Qes_i=nrgdata[nrgfpath]['i']['HFluxes'][-1]-3./2.*Ti*nrgdata[nrgfpath]['i']['PFluxes'][-1]
    Qem_e=nrgdata[nrgfpath]['e']['HFluxem'][-1]-3./2.*Te*nrgdata[nrgfpath]['e']['PFluxem'][-1]
    Qem_i=nrgdata[nrgfpath]['i']['HFluxem'][-1]-3./2.*Ti*nrgdata[nrgfpath]['i']['PFluxem'][-1]
    Q_e = (Qes_e+Qem_e)
    Q_i = (Qes_i+Qem_i)
    chi_e = (Qes_e+Qem_e) / omt_e / ne / Te
    chi_i = (Qes_i+Qem_i) / omt_i / ni / Ti

    return Qes_e,Qes_i,Qem_e,Qem_i,Q_e,Q_i,D_e,D_i,chi_e,chi_i
#**********End of input setup**********************************************************************

#**********Start of judge**************************************************************************
def Epar_judge(epar):
    #from line 78 plot_scan_info_efit.py
    #from Hatch 2016 Microtearing turbulence limiting the JET-ILW pedestal
    epar=abs(epar)
    if epar<0.4:
        mode="KBM"
    elif epar<0.7 and epar>0.3:
        mode="MTM"
    elif epar>0.9:
        mode="ITG/TEM"
    else:
        mode="other"
    return mode



def f_judge(fe, fi, omega):  #judge from the typical frquency
    ion_direction=0
    electron_direction=0
    if omega > 0: 
        ion_direction=1
    elif omega < 0: 
        eletron_direction=1
    else:
    	print("frequency = 0kHz")
    omega=abs(omega)
    ITG=0 #likely hood to be ITG in percentage
    KBM=0 #likely hood to be ITG in percentage
    MTM=0 #likely hood to be MTM in percentage
    ETG=0
    #****************ETG typical frequency*********
    f_ITG=(fe+fi/2)
    f_KBM=fi/2
    f_MTM=fe
    
    
    #ITG=(omega-f_ITG)/np.max((omega,f_ITG))
    #KBM=(omega-f_KBM)/np.max((omega,f_KBM))
    #MTM=(omega-f_MTM)/np.max((omega,f_MTM))
    ITG=abs((omega-f_ITG)/f_ITG) #likelyhood to be ITG in percentage
    ETG=0                        #likelyhood to be ETG in percentage
    KBM=abs((omega-f_KBM)/f_KBM) #likelyhood to be KBM in percentage
    MTM=abs((omega-f_MTM)/f_MTM) #likelyhood to be MTM in percentage

    mode_percent=[ITG,ETG,KBM,MTM]
    if np.argmax(mode_percent) == 0 and ion_direction==1:
        mode="ITG/TEM/KBM"
    elif np.argmax(mode_percent) == 1 and electron_direction==1:
        mode="ETG"
    elif np.argmax(mode_percent) == 2 and ion_direction==1:
        mode="ITG/TEM/KBM"
    elif np.argmax(mode_percent) == 3 and electron_direction==1:
        mode="MTM"
    else:
        mode="other"
    return mode,ITG,ETG,KBM,MTM

def D_chi_2_judge(Qes_e,Qes_i,Qem_e,Qem_i,Q_e,Q_i,D_e,D_i,chi_e,chi_i):
    chi_tot=chi_i+chi_e
    if abs(chi_i/chi_e) < 0.8 and abs(Q_i/Q_e) < 0.8 and abs(Qem_e/Qes_e) >= 1 and abs(D_e/chi_e)<0.6:
        mode="MTM"
    elif abs(chi_i/chi_e) > 0.8 and abs(Q_i/Q_e) > 0.8 and abs(Qem_e/Qes_e) >= 1 and abs(D_e/chi_e)>0.6:
        mode="KBM"
    elif abs(chi_i/chi_e) < 0.8 and abs(Q_i/Q_e) < 0.8 and abs(Qem_e/Qes_e) < 1 and abs(D_e/chi_e)>0.6:
        mode="ETG"
    elif abs(chi_i/chi_e) > 0.8 and abs(Q_i/Q_e) > 0.8 and abs(Qem_e/Qes_e) < 1 and abs(D_e/chi_tot)<0.44:
        mode="ITG/TEM"
    else:
        mode="other"
    return mode

def D_chi_3_judge(Qes_e,Qes_i,Qes_z,Qem_e,Qem_i,Qem_z,Q_e,Q_i,Q_z,D_e,D_i,D_z,chi_e,chi_i,chi_z):
    chi_tot=chi_i+chi_e+chi_z
    if abs(chi_i/chi_e) < 0.8 and abs(Q_i/Q_e) < 0.8 and abs(Qem_e/Qes_e) >= 1 and abs(D_e/chi_e)<0.6:
        mode="MTM"
    elif abs(chi_i/chi_e) > 0.8 and abs(Q_i/Q_e) > 0.8 and abs(Qem_e/Qes_e) >= 1 and abs(D_e/chi_e)>0.6:
        mode="KBM"
    elif abs(chi_i/chi_e) < 0.8 and abs(Q_i/Q_e) < 0.8 and abs(Qem_e/Qes_e) < 1 and abs(D_e/chi_e)>0.6:
        mode="ETG"
    elif abs(chi_i/chi_e) > 0.8 and abs(Q_i/Q_e) > 0.8 and abs(Qem_e/Qes_e) < 1 and abs(D_e/chi_tot)<0.44:
        mode="ITG/TEM"
    else:
        mode="other"
    return mode

#**********End of judge**************************************************************************

def federal_court(suffix):
    print("************"+str(suffix)+"*************")
#************Start of transport ratio case**************
    if len(species_order(suffix))==3:
        Qes_e,Qes_i,Qes_z,Qem_e,Qem_i,Qem_z,Q_e,Q_i,Q_z,D_e,D_i,D_z,chi_e,chi_i,chi_z=D_chi_3(suffix)
        D_chi_mode=D_chi_3_judge(Qes_e,Qes_i,Qes_z,Qem_e,Qem_i,Qem_z,Q_e,Q_i,Q_z,D_e,D_i,D_z,chi_e,chi_i,chi_z)
        chi_tot=chi_i+chi_e+chi_z
    elif len(species_order(suffix))==2:
        Qes_e,Qes_i,Qem_e,Qem_i,Q_e,Q_i,D_e,D_i,chi_e,chi_i=D_chi_2(suffix)
        D_chi_mode=D_chi_2_judge(Qes_e,Qes_i,Qem_e,Qem_i,Q_e,Q_i,D_e,D_i,chi_e,chi_i)
        chi_tot=chi_i+chi_e
    print("Based on Transport ratio, the mode is: "+str(D_chi_mode))
#************End of transport ratio case**************

#************Start of typical frequency case**************
    fe, fi, omega, gamma, kx, kymin, f_kHz, nu_ei=f(suffix)
    f_mode,f_ITG,f_ETG,f_KBM,f_MTM=f_judge(fe, fi, omega)
    print("Based on Typical frequency, the mode is: "+str(f_mode))
    print("The likelyhood of ITG is "+str(f_ITG))
    print("The likelyhood of ETG is "+str(f_ETG))
    print("The likelyhood of KBM is "+str(f_KBM))
    print("The likelyhood of MTM is "+str(f_MTM))
#************End of typical frequency case**************

#***********Start of E paraelle case*******************
    epar=Epar(suffix)
    epar_mode=Epar_judge(epar)
    print("Based on E_par Cancellation, the mode is: "+str(epar_mode))
    print("The e parelle is "+str(epar))
#***********End of E paraelle case*******************

#***********Start of the summary***********************
    if D_chi_mode=="KBM" and epar_mode=="KBM":
        tot_mode="KBM"
    elif D_chi_mode=="MTM" and f_mode=="MTM":
        tot_mode="MTM"
    elif D_chi_mode=="ITG/TEM" and f_mode=="ITG/TEM/KBM":
        tot_mode="ITG/TEM"
    elif D_chi_mode=="ETG" and f_mode=="ETG":
        tot_mode="ETG"
    else:
        tot_mode="other"
    print("Based on all the factors, the mode is: "+str(tot_mode))
    if tot_mode=="other":
        print("Please check the result manually")
#***********End of the summary***********************

#**********Start of report***************************
    with open('Court_results.csv', 'a') as csvfile:
    	csv_data = csv.writer(csvfile, delimiter=',')
        csv_data.writerow([suffix,tot_mode,\
                            kymin,kx,nu_ei,\
                            f_mode,omega,f,fi,fe,gamma,\
                            epar_mode,epar,\
                            D_chi_mode,Qem_e/Qes_e,Qem_i/Qes_i,Qem_z/Qes_z,Q_i/Q_e,chi_i/chi_e,\
                            D_i/chi_tot,D_e/chi_tot,D_z/chi_tot,D_i/chi_e,D_e/chi_e,D_z/chi_e,\
                            D_i/chi_i,D_e/chi_i,D_z/chi_i])
    csvfile.close()
#**********End of report***************************

#if report==1:
#    with open('mode_number_finder_report.csv','w') as csvfile:
#        data = csv.writer(csvfile, delimiter=',')
#        data.writerow(['n ','m ','kymin          ','frequency(kHz)           ','location(r/a)            ','omega_GENE    ','Drive'])
#       for i in range(len(n0_range)):
#            data.writerow([n0_range[i],m0_range[i],ky_range[i],f_range[i],x_range[i],f_GENE_range[i],drive_range[i]])
#    csvfile.close()

def scan_cases():
    cwd = os.getcwd()
    filelist = []
    for filename in os.listdir(cwd):
        if filename.startswith("field"):
            filelist.append(filename[-5:])
    filelist.sort()
    

    #*****Start save the data********* 
    with open('Court_results.csv', 'w') as csvfile:
        csv_data = csv.writer(csvfile, delimiter=',')
        csv_data.writerow(["Suffix","Total mode",\
                            "ky*rhoi","kx*rhoi","nu_ei(cs/a)",\
                            "f mode","omega(cs/a)","f(kHz)","omega*i(cs/a)","omega*e(cs/a)","gamma(cs/a)",\
                            "E_par mode","E_par",\
                            "D_chi mode","Qem_e/Qes_e","Qem_i/Qes_i","Qem_z/Qes_z","Q_i/Q_e","chi_i/chi_e",\
                            "D_i/chi_tot","D_e/chi_tot","D_z/chi_tot","D_i/chi_e","D_e/chi_e","D_z/chi_e",\
                            "D_i/chi_i","D_e/chi_i","D_z/chi_i"])
    csvfile.close()
    #************End save the data**************
    for suffix in filelist:
        if len(suffix)==5:
            suffix=suffix[1:5]
        federal_court(suffix)


#************Operatre***********************
scan_cases()



