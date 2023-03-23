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

def write_instate(is_dict,filename):
    f = open(filename,'w')
    f.write('&INSTATE\n')
    f.write('  TOKAMAK_ID = '+is_dict['tokamak_id']+'\n')
    f.write('  DENSITY_MODEL = '+str(int(is_dict['density_model']))+'\n')
    f.write('  MODEL_SHAPE = '+str(int(is_dict['model_shape']))+'\n')
    f.write('  R0 = '+str((is_dict['R0']))+'\n')
    f.write('  B0 = '+str((is_dict['B0']))+'\n')
    f.write('  IP = '+str((is_dict['IP']))+'\n')
    f.write('  RMAJOR = '+str((is_dict['rmajor']))+'\n')
    f.write('  AMINOR = '+str((is_dict['aminor']))+'\n')
    f.write('  KAPPA = '+str((is_dict['kappa']))+'\n')
    f.write('  DELTA = '+str((is_dict['delta']))+'\n')
    f.write('  N_ION = '+str(int(is_dict['n_ion']))+'\n')
    f.write('  Z_ION = '+str(int(is_dict['z_ion']))+'\n')
    f.write('  A_ION = '+str(int(is_dict['a_ion']))+'\n')
    f.write('  F_ION = '+str((is_dict['f_ion']))+'\n')
    f.write('  N_IMP = '+str(int(is_dict['n_imp']))+'\n')
    f.write('  Z_IMP = '+str(int(is_dict['z_imp']))+'\n')
    f.write('  A_IMP = '+str(int(is_dict['a_imp']))+'\n')
    f.write('  F_IMP = '+str((is_dict['f_imp']))+'\n')
    f.write('  N_MIN = '+str(int(is_dict['n_min']))+'\n')
    f.write('  Z_MIN = '+str(int(is_dict['z_min']))+'\n')
    f.write('  A_MIN = '+str(int(is_dict['a_min']))+'\n')
    f.write('  N_BEAM = '+str(int(is_dict['n_beam']))+'\n')
    f.write('  Z_BEAM = '+str(int(is_dict['z_beam']))+'\n')
    f.write('  A_BEAM = '+str(int(is_dict['a_beam']))+'\n')
    f.write('  N_FUSION = '+str(int(is_dict['n_fusion']))+'\n')
    f.write('  Z_FUSION = '+str(int(is_dict['z_fusion']))+'\n')
    f.write('  A_FUSION = '+str(int(is_dict['a_fusion']))+'\n')
    f.write('  NRHO = '+str(int(is_dict['nrho']))+'\n')
    f.write('  RHO = ')
    for i in range(len(is_dict['rhot'])):
        f.write(str(is_dict['rhot'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['rhot']) - 1:
            f.write('\n')
    f.write('  NE = ')
    for i in range(len(is_dict['ne'])):
        f.write(str(is_dict['ne'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['ne']) - 1:
            f.write('\n')
    f.write('  TE = ')
    for i in range(len(is_dict['te'])):
        f.write(str(is_dict['te'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['te']) - 1:
            f.write('\n')
    f.write('  TI = ')
    for i in range(len(is_dict['ti'])):
        f.write(str(is_dict['ti'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['ti']) - 1:
            f.write('\n')
    f.write('  ZEFF = ')
    for i in range(len(is_dict['zeff'])):
        f.write(str(is_dict['zeff'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['zeff']) - 1:
            f.write('\n')
    f.write('  OMEGA = ')
    for i in range(len(is_dict['omega'])):
        f.write(str(is_dict['omega'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['omega']) - 1:
            f.write('\n')
    f.write('  J_OH = ')
    for i in range(len(is_dict['j_oh'])):
        f.write(str(is_dict['j_oh'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['j_oh']) - 1:
            f.write('\n')
    f.write('  J_BS = ')
    for i in range(len(is_dict['j_bs'])):
        f.write(str(is_dict['j_bs'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['j_bs']) - 1:
            f.write('\n')
    f.write('  SE_NB = ')
    for i in range(len(is_dict['se_nb'])):
        f.write(str(is_dict['se_nb'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['se_nb']) - 1:
            f.write('\n')
    f.write('  Q = ')
    for i in range(len(is_dict['q0'])):
        f.write(str(is_dict['q0'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['q0']) - 1:
            f.write('\n')
    f.write('  PSIPOL = ')
    for i in range(len(is_dict['psipol'])):
        f.write(str(is_dict['psipol'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['psipol']) - 1:
            f.write('\n')
    f.write('  J_TOT = ')
    for i in range(len(is_dict['j_tot'])):
        f.write(str(is_dict['j_tot'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['j_tot']) - 1:
            f.write('\n')
    f.write('  P_EQ = ')
    for i in range(len(is_dict['p_eq'])):
        f.write(str(is_dict['p_eq'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['p_eq']) - 1:
            f.write('\n')
    f.write('  P_RAD = ')
    for i in range(len(is_dict['p_rad'])):
        f.write(str(is_dict['p_rad'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['p_rad']) - 1:
            f.write('\n')
    f.write('  J_NB = ')
    for i in range(len(is_dict['j_nb'])):
        f.write(str(is_dict['j_nb'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['j_nb']) - 1:
            f.write('\n')
    f.write('  J_EC = ')
    for i in range(len(is_dict['j_ec'])):
        f.write(str(is_dict['j_ec'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['j_ec']) - 1:
            f.write('\n')
    f.write('  J_IC = ')
    for i in range(len(is_dict['j_ic'])):
        f.write(str(is_dict['j_ic'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['j_ic']) - 1:
            f.write('\n')
    f.write('  PE_NB = ')
    for i in range(len(is_dict['pe_nb'])):
        f.write(str(is_dict['pe_nb'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['pe_nb']) - 1:
            f.write('\n')
    f.write('  PE_EC = ')
    for i in range(len(is_dict['pe_ec'])):
        f.write(str(is_dict['pe_ec'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['pe_ec']) - 1:
            f.write('\n')
    f.write('  PE_IC = ')
    for i in range(len(is_dict['pe_ic'])):
        f.write(str(is_dict['pe_ic'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['pe_ic']) - 1:
            f.write('\n')
    f.write('  PE_FUS = ')
    for i in range(len(is_dict['pe_fus'])):
        f.write(str(is_dict['pe_fus'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['pe_fus']) - 1:
            f.write('\n')
    f.write('  PE_IONIZATION = ')
    for i in range(len(is_dict['pe_ionization'])):
        f.write(str(is_dict['pe_ionization'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['pe_ionization']) - 1:
            f.write('\n')
    f.write('  PI_NB = ')
    for i in range(len(is_dict['pi_nb'])):
        f.write(str(is_dict['pi_nb'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['pi_nb']) - 1:
            f.write('\n')
    f.write('  PI_EC = ')
    for i in range(len(is_dict['pi_ec'])):
        f.write(str(is_dict['pi_ec'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['pi_ec']) - 1:
            f.write('\n')
    f.write('  PI_FUS = ')
    for i in range(len(is_dict['pi_fus'])):
        f.write(str(is_dict['pi_fus'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['pi_fus']) - 1:
            f.write('\n')
    f.write('  PI_IONIZATION = ')
    for i in range(len(is_dict['pi_ionization'])):
        f.write(str(is_dict['pi_ionization'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['pi_ionization']) - 1:
            f.write('\n')
    f.write('  PI_NB = ')
    for i in range(len(is_dict['pi_nb'])):
        f.write(str(is_dict['pi_nb'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['pi_nb']) - 1:
            f.write('\n')
    f.write('  PI_EC = ')
    for i in range(len(is_dict['pi_ec'])):
        f.write(str(is_dict['pi_ec'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['pi_ec']) - 1:
            f.write('\n')
    f.write('  PI_FUS = ')
    for i in range(len(is_dict['pi_fus'])):
        f.write(str(is_dict['pi_fus'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['pi_fus']) - 1:
            f.write('\n')
    f.write('  PI_IONIZATION = ')
    for i in range(len(is_dict['pi_ionization'])):
        f.write(str(is_dict['pi_ionization'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['pi_ionization']) - 1:
            f.write('\n')
    f.write('  PI_CX = ')
    for i in range(len(is_dict['pi_cx'])):
        f.write(str(is_dict['pi_cx'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['pi_cx']) - 1:
            f.write('\n')
    f.write('  P_OHM = ')
    for i in range(len(is_dict['p_ohm'])):
        f.write(str(is_dict['p_ohm'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['p_ohm']) - 1:
            f.write('\n')
    f.write('  P_EI = ')
    for i in range(len(is_dict['p_ei'])):
        f.write(str(is_dict['p_ei'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['p_ei']) - 1:
            f.write('\n')
    f.write('  TORQUE_NB = ')
    for i in range(len(is_dict['torque_nb'])):
        f.write(str(is_dict['torque_nb'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['torque_nb']) - 1:
            f.write('\n')
    f.write('  TORQUE_IN = ')
    for i in range(len(is_dict['torque_in'])):
        f.write(str(is_dict['torque_in'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['torque_in']) - 1:
            f.write('\n')
    f.write('  SE_IONIZATION = ')
    for i in range(len(is_dict['se_ionization'])):
        f.write(str(is_dict['se_ionization'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['se_ionization']) - 1:
            f.write('\n')
    f.write('  SI_NB = ')
    for i in range(len(is_dict['si_nb'])):
        f.write(str(is_dict['si_nb'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['si_nb']) - 1:
            f.write('\n')
    f.write('  SI_IONIZATION = ')
    for i in range(len(is_dict['si_ionization'])):
        f.write(str(is_dict['si_ionization'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['si_ionization']) - 1:
            f.write('\n')
    f.write('  DENSITY_BEAM = ')
    for i in range(len(is_dict['density_beam'])):
        f.write(str(is_dict['density_beam'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['density_beam']) - 1:
            f.write('\n')
    f.write('  WBEAM = ')
    for i in range(len(is_dict['wbeam'])):
        f.write(str(is_dict['wbeam'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['wbeam']) - 1:
            f.write('\n')
    f.write('  DENSITY_ALPHA = ')
    for i in range(len(is_dict['density_alpha'])):
        f.write(str(is_dict['density_alpha'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['density_alpha']) - 1:
            f.write('\n')
    f.write('  WALPHA = ')
    for i in range(len(is_dict['walpha'])):
        f.write(str(is_dict['walpha'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['walpha']) - 1:
            f.write('\n')
    f.write('  CHIE = ')
    for i in range(len(is_dict['chie'])):
        f.write(str(is_dict['chie'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['chie']) - 1:
            f.write('\n')
    f.write('  CHII = ')
    for i in range(len(is_dict['chii'])):
        f.write(str(is_dict['chii'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['chii']) - 1:
            f.write('\n')
    f.write('  NBDRY = ')
    f.write(str(int(is_dict['nbdry']))+'\n')
    f.write('  RBDRY = ')
    for i in range(len(is_dict['rbdry'])):
        f.write(str(is_dict['rbdry'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['rbdry']) - 1:
            f.write('\n')
    f.write('  ZBDRY = ')
    for i in range(len(is_dict['zbdry'])):
        f.write(str(is_dict['zbdry'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['zbdry']) - 1:
            f.write('\n')
    f.write('  NLIM = ')
    f.write(str(int(is_dict['nlim']))+'\n')
    f.write('  RLIM = ')
    for i in range(len(is_dict['rlim'])):
        f.write(str(is_dict['rlim'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['rlim']) - 1:
            f.write('\n')
    f.write('  ZLIM = ')
    for i in range(len(is_dict['zlim'])):
        f.write(str(is_dict['zlim'][i])[0:10]+' ')
        if (i+1)%4 == 0:
            f.write('\n    ')
        elif i == len(is_dict['zlim']) - 1:
            f.write('\n')
    f.write('/')













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
    j_bs = np.empty(0)
    psipol = np.empty(0)
    j_tot = np.empty(0)
    p_eq = np.empty(0)
    q0 = np.empty(0)
    rbdry = np.empty(0)
    zbdry = np.empty(0)
    rlim = np.empty(0)
    zlim = np.empty(0)
    omega = np.empty(0)
    ######
    p_rad = np.empty(0)
    j_nb = np.empty(0)
    j_ec = np.empty(0)
    j_ic = np.empty(0)
    pe_nb = np.empty(0)
    pe_ec = np.empty(0)
    pe_ic = np.empty(0)
    pe_fus = np.empty(0)
    pe_ionization = np.empty(0)
    pi_nb = np.empty(0)
    pi_ec = np.empty(0)
    pi_fus = np.empty(0)
    pi_ionization = np.empty(0)
    pi_cx = np.empty(0)
    p_ohm = np.empty(0)
    p_ei = np.empty(0)
    torque_nb = np.empty(0)
    torque_in = np.empty(0)
    se_nb = np.empty(0)
    se_ionization = np.empty(0)
    si_nb = np.empty(0)
    si_ionization = np.empty(0)
    density_beam = np.empty(0)
    wbeam = np.empty(0)
    density_alpha = np.empty(0)
    walpha = np.empty(0)
    chie = np.empty(0)
    chii = np.empty(0)

    for i in range(len(ds)):
        j = 0
        if 'NBDRY' in ds[i]:
            nbdry = float(ds[i].split()[2])
        if 'DENSITY_MODEL' in ds[i]:
            density_model = float(ds[i].split()[2])
        if 'MODEL_SHAPE' in ds[i]:
            model_shape = float(ds[i].split()[2])
        if 'TOKAMAK_ID' in ds[i]:
            tokamak_id = (ds[i].split()[2])
        if 'NLIM' in ds[i]:
            nlim = float(ds[i].split()[2])
        if 'AMINOR' in ds[i]:
            aminor = float(ds[i].split()[2])
        if 'RMAJOR' in ds[i]:
            rmajor = float(ds[i].split()[2])
        if 'Z_IMP' in ds[i]:
            zimp = float(ds[i].split()[2])
        if 'R0' in ds[i]:
            R0 = float(ds[i].split()[2])
        if 'NRHO' in ds[i]:
            nrho = float(ds[i].split()[2])
        if 'KAPPA' in ds[i]:
            kappa = float(ds[i].split()[2])
        if 'DELTA' in ds[i]:
            delta = float(ds[i].split()[2])
        if 'N_ION' in ds[i]:
            n_ion = float(ds[i].split()[2])
        if 'Z_ION' in ds[i]:
            z_ion = float(ds[i].split()[2])
        if 'A_ION' in ds[i]:
            a_ion = float(ds[i].split()[2])
        if 'F_ION' in ds[i]:
            f_ion = float(ds[i].split()[2])
        if 'N_IMP' in ds[i]:
            n_imp = float(ds[i].split()[2])
        if 'A_IMP' in ds[i]:
            a_imp = float(ds[i].split()[2])
        if 'F_IMP' in ds[i]:
            f_imp = float(ds[i].split()[2])
        if 'N_MIN' in ds[i]:
            n_min = float(ds[i].split()[2])
        if 'Z_MIN' in ds[i]:
            z_min = float(ds[i].split()[2])
        if 'A_MIN' in ds[i]:
            a_min = float(ds[i].split()[2])
        if 'N_BEAM' in ds[i]:
            n_beam = float(ds[i].split()[2])
        if 'Z_BEAM' in ds[i]:
            z_beam = float(ds[i].split()[2])
        if 'A_BEAM' in ds[i]:
            a_beam = float(ds[i].split()[2])
        if 'N_FUSION' in ds[i]:
            n_fusion = float(ds[i].split()[2])
        if 'Z_FUSION' in ds[i]:
            z_fusion = float(ds[i].split()[2])
        if 'A_FUSION' in ds[i]:
            a_fusion = float(ds[i].split()[2])
        if 'B0' in ds[i]:
            B0 = float(ds[i].split()[2])
        if 'IP' in ds[i] and 'PSIPOL' not in ds[i]:
            IP = float(ds[i].split()[2])
        if 'RHO' in ds[i] and 'NRHO' not in ds[i] and 'RHO_TYPE' not in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    rhot = np.append(rhot,float(temp[k]))
            while True:
                j += 1 
                #if j == nrho-1:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    rhot = np.append(rhot,float(temp[k]))
                if len(rhot) >= nrho:
                    break
            #print('rhot',rhot)     
        if 'NE' in ds[i] and 'TOKAMAK' not in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    ne = np.append(ne,float(temp[k]))
            while True:
                j += 1 
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    ne = np.append(ne,float(temp[k]))
                if len(ne) >= nrho:
                    break
            #print('ne',ne)     
        if 'TE' in ds[i] and 'INSTATE' not in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    te = np.append(te,float(temp[k]))
            while True:
                j += 1 
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    te = np.append(te,float(temp[k]))
                if len(te) >= nrho:
                    break
            #print('te',te)     
        if 'TI' in ds[i] and 'IONIZATION' not in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    ti = np.append(ti,float(temp[k]))
            while True:
                j += 1 
                #if 'ZEFF' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    ti = np.append(ti,float(temp[k]))
                if len(ti) >= nrho:
                    break
        if 'ZEFF' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    zeff = np.append(zeff,float(temp[k]))
            while True:
                j += 1 
                #if 'OMEGA' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    zeff = np.append(zeff,float(temp[k]))
                if len(zeff) >= nrho:
                    break
        if 'OMEGA' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    omega = np.append(omega,float(temp[k]))
            while True:
                j += 1 
                #if 'OMEGA' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    omega = np.append(omega,float(temp[k]))
                if len(omega) >= nrho:
                    break

        if 'P_EQ' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    p_eq = np.append(p_eq,float(temp[k]))
            while True:
                j += 1 
                #if 'NBDRY' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    p_eq = np.append(p_eq,float(temp[k]))
                if len(p_eq) >= nrho:
                    break

        if 'SE_NB' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    se_nb = np.append(se_nb,float(temp[k]))
            while True:
                j += 1 
                #if 'SE_IONIZATION' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    se_nb = np.append(se_nb,float(temp[k]))
                if len(se_nb) >= nrho:
                    break
        if 'J_OH' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    j_oh = np.append(j_oh,float(temp[k]))
            while True:
                j += 1 
                #if 'J_BS' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    j_oh = np.append(j_oh,float(temp[k]))
                if len(j_oh) >= nrho:
                    break
        if 'J_BS' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    j_bs = np.append(j_bs,float(temp[k]))
            while True:
                j += 1 
                #if 'J_BS' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    j_bs = np.append(j_bs,float(temp[k]))
                if len(j_bs) >= nrho:
                    break
        if 'J_TOT' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    j_tot = np.append(j_tot,float(temp[k]))
            while True:
                j += 1 
                #if 'P_EQ' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    #print("len(j_tot)",len(j_tot))
                    #print("(j_tot[-1])",j_tot[-1])
                    j_tot = np.append(j_tot,float(temp[k]))
                if len(j_tot) >= nrho:
                    break
        if 'PSIPOL' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    psipol = np.append(psipol,float(temp[k]))
            while True:
                j += 1 
                #if 'P_EQ' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    #print("len(j_tot)",len(j_tot))
                    #print("(j_tot[-1])",j_tot[-1])
                    psipol = np.append(psipol,float(temp[k]))
                if len(psipol) >= nrho:
                    break
        if ' Q ' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    #print(temp[k])
                    q0 = np.append(q0,float(temp[k]))
            while True:
                j += 1 
                #if 'PSIPOL' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    q0 = np.append(q0,float(temp[k]))
                if len(q0) >= nrho:
                    break
        if 'RBDRY' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    rbdry = np.append(rbdry,float(temp[k]))
            while True:
                j += 1 
                #if 'ZBDRY' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    rbdry = np.append(rbdry,float(temp[k]))
                if len(rbdry) >= nbdry:
                    break
        if 'ZBDRY' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    zbdry = np.append(zbdry,float(temp[k]))
            while True:
                j += 1 
                #if 'NLIM' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    zbdry = np.append(zbdry,float(temp[k]))
                if len(zbdry) >= nbdry:
                    break
        if 'RLIM' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    rlim = np.append(rlim,float(temp[k]))
            while True:
                j += 1 
                #if 'ZLIM' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    rlim = np.append(rlim,float(temp[k]))
                if len(rlim) >= nlim:
                    break
        if 'ZLIM' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    zlim = np.append(zlim,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    zlim = np.append(zlim,float(temp[k]))
                if len(zlim) >= nlim:
                    break
        if 'P_RAD' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    p_rad = np.append(p_rad,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    p_rad = np.append(p_rad,float(temp[k]))
                if len(p_rad) >= nrho:
                    break
        if 'J_NB' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    j_nb = np.append(j_nb,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    j_nb = np.append(j_nb,float(temp[k]))
                if len(j_nb) >= nrho:
                    break
        if 'J_EC' in ds[i] and ds[i][0] != '!':
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    j_ec = np.append(j_ec,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    j_ec = np.append(j_ec,float(temp[k]))
                if len(j_ec) >= nrho:
                    break
        if 'J_IC' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    j_ic = np.append(j_ic,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    j_ic = np.append(j_ic,float(temp[k]))
                if len(j_ic) >= nrho:
                    break
        if 'PE_NB' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    pe_nb = np.append(pe_nb,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    pe_nb = np.append(pe_nb,float(temp[k]))
                if len(pe_nb) >= nrho:
                    break
        if 'PE_EC' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    pe_ec = np.append(pe_ec,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    pe_ec = np.append(pe_ec,float(temp[k]))
                if len(pe_ec) >= nrho:
                    break
        if 'PE_IC' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    pe_ic = np.append(pe_ic,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    pe_ic = np.append(pe_ic,float(temp[k]))
                if len(pe_ic) >= nrho:
                    break
        if 'PE_FUS' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    pe_fus = np.append(pe_fus,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    pe_fus = np.append(pe_fus,float(temp[k]))
                if len(pe_fus) >= nrho:
                    break
        if 'PE_IONIZATION' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    pe_ionization = np.append(pe_ionization,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    pe_ionization = np.append(pe_ionization,float(temp[k]))
                if len(pe_ionization) >= nrho:
                    break
        if 'PI_NB' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    pi_nb = np.append(pi_nb,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    pi_nb = np.append(pi_nb,float(temp[k]))
                if len(pi_nb) >= nrho:
                    break
        if 'PI_EC' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    pi_ec = np.append(pi_ec,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    pi_ec = np.append(pi_ec,float(temp[k]))
                if len(pi_ec) >= nrho:
                    break
        if 'PI_FUS' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    pi_fus = np.append(pi_fus,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    pi_fus = np.append(pi_fus,float(temp[k]))
                if len(pi_fus) >= nrho:
                    break
        if 'PI_IONIZATION' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    pi_ionization = np.append(pi_ionization,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    pi_ionization = np.append(pi_ionization,float(temp[k]))
                if len(pi_ionization) >= nrho:
                    break
        if 'PI_CX' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    pi_cx = np.append(pi_cx,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    pi_cx = np.append(pi_cx,float(temp[k]))
                if len(pi_cx) >= nrho:
                    break
        if 'P_OHM' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    p_ohm = np.append(p_ohm,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    p_ohm = np.append(p_ohm,float(temp[k]))
                if len(p_ohm) >= nrho:
                    break
        if 'P_EI' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    p_ei = np.append(p_ei,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    p_ei = np.append(p_ei,float(temp[k]))
                if len(p_ei) >= nrho:
                    break
        if 'TORQUE_NB' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    torque_nb = np.append(torque_nb,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    torque_nb = np.append(torque_nb,float(temp[k]))
                if len(torque_nb) >= nrho:
                    break
        if 'TORQUE_IN' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    torque_in = np.append(torque_in,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    torque_in = np.append(torque_in,float(temp[k]))
                if len(torque_in) >= nrho:
                    break
#        if 'SE_NB' in ds[i]:
#            temp = ds[i].split()
#            for k in range(len(temp)):
#                if k>1:
#                    se_nb = np.append(se_nb,float(temp[k]))
#            while True:
#                j += 1 
#                #if '/' in ds[i+j]:
#                #    break
#                temp = ds[i+j].split()
#                for k in range(len(temp)):
#                    #print(temp[k])
#                    se_nb = np.append(se_nb,float(temp[k]))
#                if len(se_nb) >= nrho:
#                    break
        if 'SE_IONIZATION' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    se_ionization = np.append(se_ionization,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    se_ionization = np.append(se_ionization,float(temp[k]))
                if len(se_ionization) >= nrho:
                    break
        if 'SI_NB' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    si_nb = np.append(si_nb,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    si_nb = np.append(si_nb,float(temp[k]))
                if len(si_nb) >= nrho:
                    break
        if 'SI_IONIZATION' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    si_ionization = np.append(si_ionization,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    si_ionization = np.append(si_ionization,float(temp[k]))
                if len(si_ionization) >= nrho:
                    break
        if 'DENSITY_BEAM' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    density_beam = np.append(density_beam,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    density_beam = np.append(density_beam,float(temp[k]))
                if len(density_beam) >= nrho:
                    break
        if 'WBEAM' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    wbeam = np.append(wbeam,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    wbeam = np.append(wbeam,float(temp[k]))
                if len(wbeam) >= nrho:
                    break
        if 'DENSITY_ALPHA' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    density_alpha = np.append(density_alpha,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    density_alpha = np.append(density_alpha,float(temp[k]))
                if len(density_alpha) >= nrho:
                    break
        if 'WALPHA' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    walpha = np.append(walpha,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    walpha = np.append(walpha,float(temp[k]))
                if len(walpha) >= nrho:
                    break
        if 'CHIE' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    chie = np.append(chie,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    chie = np.append(chie,float(temp[k]))
                if len(chie) >= nrho:
                    break
        if 'CHII' in ds[i]:
            temp = ds[i].split()
            for k in range(len(temp)):
                if k>1:
                    chii = np.append(chii,float(temp[k]))
            while True:
                j += 1 
                #if '/' in ds[i+j]:
                #    break
                temp = ds[i+j].split()
                for k in range(len(temp)):
                    #print(temp[k])
                    chii = np.append(chii,float(temp[k]))
                if len(chii) >= nrho:
                    break

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
    is_data['j_bs'] = j_bs
    is_data['psipol'] = psipol
    is_data['p_eq'] = p_eq
    is_data['q0'] = q0
    is_data['zimp'] = zimp
    is_data['R0'] = R0
    is_data['aminor'] = aminor
    is_data['rmajor'] = rmajor
    is_data['zlim'] = zlim
    is_data['rlim'] = rlim
    is_data['zbdry'] = zbdry
    is_data['rbdry'] = rbdry
    is_data['nbdry'] = nbdry
    is_data['nlim'] = nlim
    is_data['omega'] = omega
    is_data['tokamak_id'] = tokamak_id
    is_data['density_model'] = density_model
    is_data['model_shape'] = model_shape
    is_data['B0'] = B0
    is_data['IP'] = IP
    is_data['kappa'] = kappa
    is_data['delta'] = delta
    is_data['n_ion'] = n_ion
    is_data['z_ion'] = z_ion
    is_data['a_ion'] = a_ion
    is_data['f_ion'] = f_ion
    is_data['n_ion'] = n_ion
    is_data['a_imp'] = a_imp
    is_data['n_min'] = n_min
    is_data['z_min'] = z_min
    is_data['a_min'] = a_min
    is_data['n_beam'] = n_beam
    is_data['z_beam'] = z_beam
    is_data['a_beam'] = a_beam
    is_data['n_fusion'] = n_fusion
    is_data['z_fusion'] = z_fusion
    is_data['a_fusion'] = a_fusion
    is_data['n_imp'] = n_imp
    is_data['z_imp'] = zimp
    is_data['f_imp'] = f_imp
    is_data['nrho'] = nrho
    is_data['p_rad'] = p_rad
    is_data['j_nb'] =j_nb
    is_data['j_ec'] =j_ec
    is_data['j_ic'] =j_ic
    is_data['pe_nb'] =pe_nb
    is_data['pe_ec'] =pe_ec
    is_data['pe_ic'] =pe_ic
    is_data['pe_fus'] =pe_fus
    is_data['pe_ionization'] =pe_ionization
    is_data['pi_nb'] =pi_nb
    is_data['pi_ec'] =pi_ec
    is_data['pi_fus'] =pi_fus
    is_data['pi_ionization'] =pi_ionization
    is_data['pi_cx'] =pi_cx
    is_data['p_ohm'] =p_ohm
    is_data['p_ei'] =p_ei
    is_data['torque_nb'] =torque_nb
    is_data['torque_in'] =torque_in
    is_data['se_nb'] =se_nb
    is_data['se_ionization'] =se_ionization
    is_data['si_nb'] =si_nb
    is_data['si_ionization'] =si_ionization
    is_data['density_beam'] =density_beam
    is_data['wbeam'] =wbeam
    is_data['density_alpha'] =density_alpha
    is_data['walpha'] =walpha
    is_data['chie'] =chie
    is_data['chii'] =chii




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





               
            
            
            
    






