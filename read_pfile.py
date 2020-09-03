import numpy as np
from interp import *
import matplotlib.pyplot as plt


def read_pfile(p_file_name,Z,add_impurity=False):

    impurity_charge = float(Z)

    f=open(p_file_name,'r')
    data = f.read()
    f.close()

    sdata = data.split('\n')
    nr = int(sdata[0].split()[0]) 

    print("p-file resolution: nr = "+str(nr) )

    name_list=[]
    
    for i in range(len(sdata)/(nr+1)):   #scan all of the quantitites in the p file
        #print(sdata[i*nr+i].split()[2])
        name_list.append(sdata[i*nr+i].split()[2])

    need_list=['ne','ni','te','ti','er','vtor']   #List of quantities that need for to be read
    need_list_number=[]                           #Number of the namelist that matches with the list of quantities that need for to be read
    
    for i in range(len(name_list)):
        for j in range(len(need_list)):
            #print(name_list[i])
            if str(need_list[j]) in str(name_list[i]):
                need_list_number.append(i)
  # more work on the error check
    if len(need_list) > len(need_list_number): 
        print('Error, missing needed quantities in q file, check line 25 in read_qfile.py')
        #for j in range(len(need_list)):
        #    if need_list_number need_list[j]
        #        print('Error, missing'+str(need_list[j]))
    
    psi=[]
    f=[]
    df=[]
    
    #print(len(need_list_number))
    for i in range(len(need_list_number)):
        n_temp=int(need_list_number[i])   #take the i_th data set
        for j in range(nr):
            temp = sdata[n_temp*nr+j+n_temp+1].split()
            psi.append(float(temp[0]))
            f.append(float(temp[1]))
            df.append(float(temp[2]))

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


    
def read_pfile_raw(p_file_name):

    f=open(p_file_name,'r')
    data = f.read()
    f.close()

    sdata = data.split('\n')
    nr = int(sdata[0].split()[0])
    print("p file resolution: nr = ", nr)

    # ne is electron density profile, dne is gradient of ne
    ne = np.empty(0)
    dne = np.empty(0)
    ni = np.empty(0)
    dni = np.empty(0)
    te = np.empty(0)
    dte = np.empty(0)
    ti = np.empty(0)
    dti = np.empty(0)
    er = np.empty(0)
    der = np.empty(0)
    # psipne is the grid of psi_pol on which ne&dne above is recorded
    psipne = np.empty(0)
    psipni = np.empty(0)
    psipte = np.empty(0)
    psipti = np.empty(0)
    psiper = np.empty(0)

    for i in np.array(range(nr)):
        temp = sdata[i+1].split()
        psipne = np.append(psipne,float(temp[0]))
        ne = np.append(ne,float(temp[1]))
        dne = np.append(dne,float(temp[2]))
        temp = sdata[nr+i+2].split()
        psipte = np.append(psipte,float(temp[0]))
        te = np.append(te,float(temp[1]))
        dte = np.append(dte,float(temp[2]))
        temp = sdata[2*nr+i+3].split()
        psipni = np.append(psipni,float(temp[0]))
        ni = np.append(ni,float(temp[1]))
        dni = np.append(dni,float(temp[2]))
        temp = sdata[3*nr+i+4].split()
        psipti = np.append(psipti,float(temp[0]))
        ti = np.append(ti,float(temp[1]))
        dti = np.append(dti,float(temp[2]))
        temp = sdata[16*nr+i+17].split()
        psiper = np.append(psiper,float(temp[0]))
        er = np.append(er,float(temp[1]))
        der = np.append(der,float(temp[2]))

    return psipne,ne,psipte,te,psipni,ni,psipti,ti,psiper,er
