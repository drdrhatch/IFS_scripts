import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

def get_nc_file(filename):
    ds = nc.Dataset(filename)
    #print metadata:  print(ds.__dict__)
    #print dimensions:for dim in ds.dimensions.values():
    #    print(dim)
    #for var in ds.variables.values():
    #    print(var)
    #Accessing data values: prcp = ds['prcp'][:]
    return ds

def get_profiles_ffile(filename):
    ds = get_nc_file(filename)
    fdict = {}
    rhot = ds['rho'][:]
    te = ds['te'][:]  #KeV
    ti = ds['ti'][:]
    ne = ds['ne'][:]  #10^19
    ni = ds['ni'][:]
    fdict['taue'] = ds['taue'][1]
    fdict['taui'] = ds['taui'][1]
    fdict['tauth'] = ds['tauth'][1]
    fdict['tautot'] = ds['tautot'][1]
    fdict['tau98'] = ds['tau98'][1]
    fdict['we'] = ds['we'][1]
    fdict['wi'] = ds['wi'][1]
    fdict['sion'] = ds['sion'][:]
    fdict['ne'] = ne
    fdict['ni'] = ni
    fdict['te'] = te
    fdict['ti'] = ti
    #####Current
    fdict['j_tot'] = ds['j_tot'][:] 
    fdict['j_bs'] = ds['j_bs'][:] 
    fdict['j_nb'] = ds['j_nb'][:] 
    fdict['j_rf'] = ds['j_rf'][:] 
    fdict['j_oh'] = ds['j_oh'][:] 
    fdict['rhot'] = rhot
    fdict['zeff'] = ds['zeff'][:] 
    fdict['zimp'] = ds['zimp'][:]
    fdict['fluxe_exp'] = ds['fluxe_exp'][1,:]
    fdict['fluxi_exp'] = ds['fluxi_exp'][1,:]
    fdict['fluxe'] = ds['fluxe'][1,:]
    fdict['fluxi'] = ds['fluxi'][1,:]
    fdict['chii'] = ds['chii'][1,:]
    fdict['chie'] = ds['chie'][1,:]
    fdict['chii_exp'] = ds['chii_exp'][1,:]
    fdict['chie_exp'] = ds['chie_exp'][1,:]
    fdict['chin'] = ds['chin'][1,:]
    fdict['q'] = ds['q'][1,:]
    fdict['betan_loc'] = ds['betan_loc'][1,:]  #local betaN
    fdict['betan'] = ds['betan'][1]  
    fdict['r0'] = ds['r0'][1]  #Major R
    fdict['a0'] = ds['a0'][1]  #Major R
    fdict['ip'] = ds['ip'][1]  #Plasma current
    fdict['b0'] = ds['b0'][1]  #B0
    fdict['time'] = ds['time'][:]  
    fdict['nz'] = fdict['ne']*(fdict['zeff']-1)/(fdict['zimp']**2-fdict['zimp'])
    fdict['ni2'] = fdict['ne']-fdict['zimp']*fdict['nz']
    fdict['sn'] = ds['sn'][1]
    fdict['pe_nb'] = ds['pe_nb'][1,:]
    fdict['pe_rf'] = ds['pe_rf'][1,:]
    fdict['p_rad'] = ds['p_rad'][1,:]
    fdict['p_ohm'] = ds['p_ohm'][1,:]
    fdict['p_ei'] = ds['p_ei'][1,:]
    fdict['pe_fus'] = ds['pe_fus'][1,:]
    fdict['pe_ionization'] = ds['pe_ionization'][1,:]
    fdict['pi_nb'] = ds['pi_nb'][1,:]
    fdict['pi_rf'] = ds['pi_rf'][1,:]
    fdict['pi_fus'] = ds['pi_fus'][1,:]
    fdict['pi_cx'] = ds['pi_cx'][1,:]
    fdict['pi_ionization'] = ds['pi_ionization'][1,:]
    fdict['pe'] = ds['pe'][1,:]
    fdict['pi'] = ds['pi'][1,:]
    #for var in ds.variables.values():
    #    print(var)
    #stop
    return fdict

def get_profiles_statefile(filename):
    ds = nc.Dataset(filename)

def read_instate(filename):
    is_data = {}
    f = open(filename,'r')
    ds = f.read().split('\n')
    rhot = np.empty(0)
    ne = np.empty(0)
    te = np.empty(0)
    ti = np.empty(0)
    zeff = np.empty(0)
    se_nb = np.empty(0)
    j_oh = np.empty(0)
    j_tot = np.empty(0)
    p_eq = np.empty(0)
    q0 = np.empty(0)
    rbdry = np.empty(0)
    zbdry = np.empty(0)
    rlim = np.empty(0)
    zlim = np.empty(0)
    for i in range(len(ds)):
        j = 0
        if 'NBDRY' in ds[i]:
            nbdry = float(ds[i].split()[2])
        if 'NLIM' in ds[i]:
            nlim = float(ds[i].split()[2])
        if 'AMINOR' in ds[i]:
            aminor = float(ds[i].split()[2])
        if 'RMAJOR' in ds[i]:
            rmajor = float(ds[i].split()[2])
        if 'Z_IMP' in ds[i]:
            zimp = float(ds[i].split()[2])
        if 'RHO' in ds[i] and 'NRHO' not in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    rhot = np.append(rhot,float(temp[k]))
            while True:
                j += 1 
                if 'NE' in ds[i+j]:
                    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    rhot = np.append(rhot,float(temp[k]))
            #print('rhot',rhot)     
        if 'NE' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    ne = np.append(ne,float(temp[k]))
            while True:
                j += 1 
                if 'TE' in ds[i+j]:
                    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    ne = np.append(ne,float(temp[k]))
            #print('ne',ne)     
        if 'TE' in ds[i] and 'INSTATE' not in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    te = np.append(te,float(temp[k]))
            while True:
                j += 1 
                if 'TI' in ds[i+j]:
                    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    te = np.append(te,float(temp[k]))
            #print('te',te)     
        if 'TI' in ds[i] and 'IONIZATION' not in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    ti = np.append(ti,float(temp[k]))
            while True:
                j += 1 
                if 'ZEFF' in ds[i+j]:
                    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    ti = np.append(ti,float(temp[k]))
        if 'ZEFF' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    zeff = np.append(zeff,float(temp[k]))
            while True:
                j += 1 
                if 'OMEGA' in ds[i+j]:
                    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    zeff = np.append(zeff,float(temp[k]))
        if 'P_EQ' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    p_eq = np.append(p_eq,float(temp[k]))
            while True:
                j += 1 
                if 'NBDRY' in ds[i+j]:
                    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    p_eq = np.append(p_eq,float(temp[k]))

        if 'SE_NB' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    se_nb = np.append(se_nb,float(temp[k]))
            while True:
                j += 1 
                if 'SE_IONIZATION' in ds[i+j]:
                    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    se_nb = np.append(se_nb,float(temp[k]))
        if 'J_OH' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    j_oh = np.append(j_oh,float(temp[k]))
            while True:
                j += 1 
                if 'J_BS' in ds[i+j]:
                    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    j_oh = np.append(j_oh,float(temp[k]))
        if 'J_TOT' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    j_tot = np.append(j_tot,float(temp[k]))
            while True:
                j += 1 
                if 'P_EQ' in ds[i+j]:
                    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    j_tot = np.append(j_tot,float(temp[k]))
        if ' Q ' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    q0 = np.append(q0,float(temp[k]))
            while True:
                j += 1 
                if 'PSIPOL' in ds[i+j]:
                    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    q0 = np.append(q0,float(temp[k]))
        if 'RBDRY' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    rbdry = np.append(rbdry,float(temp[k]))
            while True:
                j += 1 
                if 'ZBDRY' in ds[i+j]:
                    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    rbdry = np.append(rbdry,float(temp[k]))
        if 'ZBDRY' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    zbdry = np.append(zbdry,float(temp[k]))
            while True:
                j += 1 
                if 'NLIM' in ds[i+j]:
                    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    zbdry = np.append(zbdry,float(temp[k]))
        if 'RLIM' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    rlim = np.append(rlim,float(temp[k]))
            while True:
                j += 1 
                if 'ZLIM' in ds[i+j]:
                    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    rlim = np.append(rlim,float(temp[k]))
        if 'ZLIM' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    zlim = np.append(zlim,float(temp[k]))
            while True:
                j += 1 
                if '/' in ds[i+j]:
                    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    zlim = np.append(zlim,float(temp[k]))

    is_data['rhot'] = rhot
    is_data['ne'] = ne
    is_data['te'] = te
    is_data['ti'] = ti
    is_data['zeff'] = zeff
    is_data['nz'] = ne*(zeff-1)/(zimp**2-zimp)
    is_data['ni'] = ne - zimp*is_data['nz']
    is_data['se_nb'] = se_nb
    is_data['j_tot'] = j_tot
    is_data['j_oh'] = j_oh
    is_data['p_eq'] = p_eq
    is_data['q0'] = q0
    is_data['zimp'] = zimp
    is_data['aminor'] = aminor
    is_data['rmajor'] = rmajor
    is_data['zlim'] = zlim
    is_data['rlim'] = rlim
    is_data['zbdry'] = zbdry
    is_data['rbdry'] = rbdry
    is_data['nbdry'] = nbdry
    is_data['nlim'] = nlim
    return is_data
    #plt.plot(rhot,te)
    #plt.plot(rhot,ti)
    #plt.show()
    #plt.plot(rhot,ne)
    #plt.show()

def output_four_col(arr,filename):
    f = open(filename,'w')
    for i in range(len(arr)):
        f.write(str(arr[i])+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
    f.close()

def read_inprof(filename):
    is_data = {}
    f = open(filename,'r')
    ds = f.read().split('\n')
    rhot = np.empty(0)
    ne = np.empty(0)
    te = np.empty(0)
    ti = np.empty(0)
    zeff = np.empty(0)
    for i in range(len(ds)):
        j = 0
        if 'RHO' in ds[i] and 'NRHO' not in ds[i]:
            temp = ds[i].split()
            nrho = int(float(temp[2]))
            while True:
                j += 1 
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    rhot = np.append(rhot,float(temp[k]))
                if len(rhot) == nrho:
                    break
        if 'NE' in ds[i] and 'QIONE' not in ds[i]:
            temp = ds[i].split()
            nne = int(float(temp[2]))
            while True:
                j += 1 
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    ne = np.append(ne,float(temp[k]))
                if len(ne) >= nne:
                    break
            #print('ne',ne)     
        if 'TE' in ds[i]:
            temp = ds[i].split()
            nte = int(float(temp[2]))
            while True:
                j += 1 
                if 'TI' in ds[i+j]:
                    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    te = np.append(te,float(temp[k]))
                if len(te) >= nte:
                    break
            #print('te',te)     
        if 'TI' in ds[i] and 'IONIZATION' not in ds[i]:
            temp = ds[i].split()
            nti = int(float(temp[2]))
            while True:
                j += 1 
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    ti = np.append(ti,float(temp[k]))
                if len(ti) >= nti:
                    break
        if 'ZEFF' in ds[i]:
            temp = ds[i].split()
            nzeff = int(float(temp[2]))
            while True:
                j += 1 
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    zeff = np.append(zeff,float(temp[k]))
                if len(zeff) >= nzeff:
                    break

    is_data['rhot'] = rhot
    is_data['ne'] = ne
    is_data['te'] = te
    is_data['ti'] = ti
    is_data['zeff'] = zeff
    #plt.plot(rhot,te)
    #plt.plot(rhot,ti)
    #plt.show()
    #plt.plot(rhot,ne)
    #plt.show()
    return is_data

def read_infastran(filename):
    f = open(filename,'r')
    data = f.read()
    f.close()
    if_dict = {}
    lines = data.split('\n')
    for i in lines:
        if 'dt' in i and 'relax' not in i:
            dt = float(i.split()[2])
    if_dict['dt'] = dt
    return if_dict





               
            
            
            
    






