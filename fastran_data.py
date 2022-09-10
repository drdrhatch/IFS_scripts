import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

def get_nc_file(filename):
    ds = nc.Dataset(filename)
    return ds

def get_profiles_ffile(filename):
    ds = get_nc_file(filename)
    rhot = ds['rho'][:]
    te = ds['te'][:]  #KeV
    ti = ds['ti'][:]
    ne = ds['ne'][:]  #10^19
    ni = ds['ni'][:]
    return rhot,te,ne,ti,ni

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
    for i in range(len(ds)):
        j = 0
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

    is_data['rhot'] = rhot
    is_data['ne'] = ne
    is_data['te'] = te
    is_data['ti'] = ti
    is_data['zeff'] = zeff
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



               
            
            
            
    






