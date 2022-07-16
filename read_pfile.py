import numpy as np
from interp import *
import matplotlib.pyplot as plt
import sys
import math
from finite_differences_x import *
from read_EFIT import read_EFIT

#Created by Max T. Curie  11/02/2020
#Last edited by Max Curie 11/02/2020
#Supported by scripts in IFS

def read_pfile(p_file_name,Z,add_impurity=False):

    impurity_charge = float(Z)

    f=open(p_file_name,'r')
    data = f.read()
    f.close()

    sdata = data.split('\n')
    nr = int(sdata[0].split()[0]) 

    print("p-file resolution: nr = "+str(nr) )

    if (nr+1)*int(len(sdata)/(nr+1))!=sdata:
        print('data need to be managed manually, please look into the code output carefully. ')

    name_list=[]
    
    #for i in range(int(len(sdata)/(nr+1))):   #scan all of the quantitites in the p file
        #print(sdata[i*nr+i].split()[2])
        #name_list.append(sdata[i*nr+i].split()[2])

    need_list=['ne','ni','te','ti','er','vtor']   #List of quantities that need for to be read
    need_list_number=[-1]*len(need_list)          #Number of the namelist that matches with the list of quantities that need for to be read
    
    print('All the quantities listed in the p file: ')
    #print(name_list)
    for i in range(len(need_list)):
        for j in range(len(sdata)):
            if len(sdata[j].split()) >= 3:
                if str(need_list[i]) in str(sdata[j].split()[2]):
                    need_list_number[i]=j

    case=0

    print('list of quantities name needed')
    print(need_list)
    print('number coorepsonds to the name list, -1 means missing')
    print(need_list_number)
    #if len(need_list) > len(need_list_number): 
    if -1 in need_list_number:
        print('Error, missing needed quantities in p file, check line 31 need_list in read_qfile.py')
        if need_list_number[4]==-1:
            need_list_number[4]=0
            temp0=input("The array er is missing, force the array to zero, countiune: 1:Yes, 2.No         ")
            if temp0==2:
                sys.exit()
            elif temp0==1:
                case=1    #er missing 
                
        elif need_list_number[5]==-1:
            temp0=input("The array vtor is missing, force the array to zero, countiune: 1:Yes, 2.No        ")
            need_list_number[5]=0
            if temp0==2:
                sys.exit()
            elif temp0==1 and case==1:
                print("Both array vtor and er are missing, cannot calculate doppler shift")
                case=3   #both 
            elif temp0==1 and case==0:
                case=2	 #vtor missing
                

        
        #for j in range(len(need_list)):
        #    if need_list_number need_list[j]
        #        print('Error, missing'+str(need_list[j]))
    print(str(name_list))
    print(str(need_list))
    print(str(need_list_number))

    psi=[]
    f=[]
    df=[]
    
    #print(len(need_list_number))
    for i in range(len(need_list_number)):
        n_temp=int(need_list_number[i])   #take the i_th data set
        for j in range(nr):
            temp = sdata[n_temp+j+1].split()
            psi.append(float(temp[0]))
            f.append(float(temp[1]))
            df.append(float(temp[2]))
            
        #plt.clf()
        #plt.plot(psi,f)
        #plt.show()

# ne is electron density profile, dne is gradient of ne
# psipne is the grid of psi_pol on which ne&dne above is recorded
    temp_i=0
    psipne = psi[temp_i*nr:(temp_i+1)*nr]
    temp_i=temp_i+1
    psipni = psi[temp_i*nr:(temp_i+1)*nr]
    temp_i=temp_i+1
    psipte = psi[temp_i*nr:(temp_i+1)*nr]
    temp_i=temp_i+1
    psipti = psi[temp_i*nr:(temp_i+1)*nr]
    temp_i=temp_i+1
    psiper = psi[temp_i*nr:(temp_i+1)*nr]
    temp_i=temp_i+1
    psipvtor = psi[temp_i*nr:(temp_i+1)*nr]
    
    temp_i=0
    ne = f[temp_i*nr:(temp_i+1)*nr]
    temp_i=temp_i+1
    ni = f[temp_i*nr:(temp_i+1)*nr]
    temp_i=temp_i+1
    te = f[temp_i*nr:(temp_i+1)*nr]
    temp_i=temp_i+1
    ti = f[temp_i*nr:(temp_i+1)*nr]
    temp_i=temp_i+1
    er = f[temp_i*nr:(temp_i+1)*nr]
    temp_i=temp_i+1
    vtor = f[temp_i*nr:(temp_i+1)*nr]

    temp_i=0
    dne = df[temp_i*nr:(temp_i+1)*nr]
    temp_i=temp_i+1
    dni = df[temp_i*nr:(temp_i+1)*nr]
    temp_i=temp_i+1
    dte = df[temp_i*nr:(temp_i+1)*nr]
    temp_i=temp_i+1
    dti = df[temp_i*nr:(temp_i+1)*nr]
    temp_i=temp_i+1
    der = df[temp_i*nr:(temp_i+1)*nr]
    temp_i=temp_i+1
    dvtor = df[temp_i*nr:(temp_i+1)*nr]

    if case==1: #er missing 
        er=np.zeros(len(er))
        der=np.zeros(len(der))
    elif case==2: #vtor missing
        vtor=np.zeros(len(vtor))
        dvtor=np.zeros(len(dvtor))
    elif case==3: #both are missing
        vtor=np.zeros(len(vtor))
        dvtor=np.zeros(len(dvtor))
        er=np.zeros(len(er))
        der=np.zeros(len(der))

    # psipne is the grid of psi_pol on which ne&dne above is recorded
    

    #plt.plot(psipne,'x')
    #plt.plot(psipni,'r.')
    #plt.plot(psipti,'bo')
    #plt.plot(psiper,'g.')
    #plt.plot(psipvtor,'rx')
    #plt.show()

    #psi0 is normalized psi_pol in [0,1] with 1000 points
    #quantities with _out are interpolated on psi0 grid 
    #print "length of arrays",len(ne),len(ni),len(te),len(ti)
    #psi0 = np.arange(1000)/999.0
    #psi0 = np.linspace(0.,1.,len(ne),endpoint=False)
    psi0 = np.linspace(0.,1.,400)
    ne_out = interp(psipne,ne,psi0)
    te_out = interp(psipte,te,psi0)
    ni_out = interp(psipni,ni,psi0)
    ti_out = interp(psipti,ti,psi0)
    er_out = interp(psiper,er,psi0)
    vtor_out = interp(psipvtor,vtor,psi0)

    dne_out = interp(psipne,dne,psi0)
    dte_out = interp(psipte,dte,psi0)
    dni_out = interp(psipni,dni,psi0)
    dti_out = interp(psipti,dti,psi0)
    der_out = interp(psiper,der,psi0)
    dvtor_out= interp(psipvtor,dvtor,psi0)

    
    if add_impurity:
        nz_out = np.empty(len(ne_out))
        for i in range(len(ne_out)):
            nz_out[i] = (ne_out[i]-ni_out[i])/impurity_charge
    else:
        nz_out = np.zeros(len(ne_out))


    # output quantities are psi_pol, ne, te, ni, ti
    # grid: even psi_pol_norm
    # resolution: 1000
    return psi0, ne_out, te_out, ni_out, ti_out, nz_out, er_out, vtor_out




    
def p_to_iterdb_format(p_file_name,geomfile_name):
    impurityCharge=float(input("Impurity Charge:"))
    psi0, ne0, te0, ni0, ti0, nz0, er0, vtor_out = read_pfile(p_file_name,impurityCharge,add_impurity=True)
    case=0
    if sum(er0)==0:
        print('Er is empty, using vtor to calculate Shear')
        case=1
    elif sum(vtor_out)==0:
        print('vtor is empty, using Er to calculate Shear')
        case=2
    elif sum(er0)!=0 and sum(er0)!=0:
        print('Neither Er nor vtor is empty, using both to calculate')
        case=3
    elif sum(er0)==0 and sum(er0)==0:
        print('Both Er and vtor are empty, cannot calculate Shear')
        case=4

    zeff = (ni0 + nz0 * impurityCharge**2.) / ne0 

    EFITdict = read_EFIT(geomfile_name)
    print(str(list(EFITdict.keys())))

    sepInd = np.argmin(abs(EFITdict['psipn'] - 1.))
    print('index at psipn = 1 is '+str(sepInd) )
    Rsep = EFITdict['R'][sepInd]
    print('major R(m) at psipn = 1 is '+str(Rsep))
    print('major R(m) at index = 1 is '+str(EFITdict['R'][0]))
    
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


    uni_rhot = np.linspace(min(rhot0),max(rhot0),len(rhot0)*10)
    ti_u = interp(rhot0,ti0,uni_rhot)
    te_u = interp(rhot0,te0,uni_rhot)
    ne_u = interp(rhot0,ne0,uni_rhot)
    ni_u = interp(rhot0,ni0,uni_rhot)
    nz_u = interp(rhot0,nz0,uni_rhot)
    p_u = (ni_u + nz_u) * ti_u + ne_u * te_u


    tprime_i = -first_derivative(ti_u,uni_rhot)/ti_u
    tprime_e = -first_derivative(te_u,uni_rhot)/te_u
    nprime_e = -first_derivative(ne_u,uni_rhot)/ne_u
    nprime_i = -first_derivative(ni_u,uni_rhot)/ni_u
    nprime_z = -first_derivative(nz_u,uni_rhot)/nz_u
    eta_i = tprime_i / nprime_i
    eta_e = tprime_e / nprime_e
    eta_z = tprime_i / nprime_z


    # convert from kV/m to V/m
    Er_Vm = interp(rhot0,er0,uni_rhot)*1E3


    R_u = interp(EFITdict['rhotn'],EFITdict['R'],uni_rhot)
    Bpol_u = interp(EFITdict['rhotn'],EFITdict['Bpol'],uni_rhot)
    vtor_out_u = interp(psi0,vtor_out,uni_rhot)
    
    # add minus sign for consistency
    omega_tor_Er = - Er_Vm / (R_u * Bpol_u)
    #print(R_u[0])
    omega_tor_Vor = vtor_out_u*1000. / (R_u)

    if case==1:
        omega_tor=omega_tor_Vor
    if case==2:
        omega_tor=omega_tor_Er
    if case==3:
        if sum(abs((omega_tor_Er-omega_tor_Vor)/omega_tor_Vor))>0.05*float(len(omega_tor_Vor)):
            print("Too much difference between omega_tor calculated from Er and vtor")
            plt.clf()
            plt.plot(uni_rhot,omega_tor_Er,label='omega_tor_Er')
            plt.plot(uni_rhot,omega_tor_Vor,label='omega_tor_Vor')
            plt.xlabel('rhot')
            plt.legend()
            plt.show()

            decide=int(input("omega_tor_Er or omega_tor_Vor, 1. omega_tor_Er, 2. omega_tor_Vor:      "))

            if decide==1:
                omega_tor=omega_tor_Er
            elif decide==2:
                omega_tor=omega_tor_Vor
            else:
                print("please input 1 or 2")
        else:
            omega_tor=omega_tor_Vor
    if case==4:
        print('Both Er and vtor are empty, cannot calculate Shear')
        omega_tor=omega_tor_Er
    # densities are multiplied by 10 here 
    # because output_iterdb() expects density in 10^19 m^-3
    
    psi_u = interp(rhot0,psi0,uni_rhot)
    rhop_u = np.sqrt(np.array(psi_u))

    return uni_rhot, rhop_u, 1000.0*te_u, 1000.0*ti_u, 1.0e19*ne_u*10., 1.0e19*ni_u*10., omega_tor

def read_pfile_direct(file_name):
    '''
    This function reads in the pfile with no modification and takes names directly from the file
    '''
    f = open(file_name,'r')
    data = f.read()
    f.close()
    lines = data.split('\n')
    entry_length = int(float(lines[0].split()[0]))
    print('Entry length',entry_length)

    pfile_dict = {}
    count = 0
    for i in range(len(lines)):
        this_line = lines[i]
        ltemp = this_line.split()
        if int(float(ltemp[0])) == entry_length:
            print(count)
            count += 1
            this_psinorm = ltemp[1]+'_'+ltemp[2]
            this_name1 = ltemp[2]
            this_name2 = ltemp[3]
            pfile_dict[this_psinorm] = np.empty(entry_length) 
            pfile_dict[this_name1] = np.empty(entry_length)
            pfile_dict[this_name2] = np.empty(entry_length)
            if count == 1:
                psiname = ltemp[1]
                pfile_dict[psiname] = np.empty(entry_length)
        elif 'SPECIES' in this_line:
            nspec = int(float(this_line.split()[0]))
            pfile_dict['species info'] = lines[i]
            for j in range(nspec):
                pfile_dict['species info'] += '\n'+lines[i+j+1]
            break
        else:
            if count == 1:
                np.append(pfile_dict[psiname],ltemp[0])
            np.append(pfile_dict[this_psinorm],ltemp[0])
            np.append(pfile_dict[this_name1],ltemp[1])
            np.append(pfile_dict[this_name2],ltemp[2])

    return pfile_dict



