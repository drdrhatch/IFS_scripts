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







