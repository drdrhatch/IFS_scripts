import numpy as np
import matplotlib.pyplot as plt
import re

######################modify######################3
######################modify######################3
######################modify######################3
#file_name='profiles_nshift0.02_t3.25.iterdb'
#file_name='profiles_3.25.iterdb'
#file_name = 'profiles_t3.035_nshift0.02.iterdb'
#file_name = 'geqdsk_89453_negom.iterdb'
#file_name = 'efit_Dial_Nmod1_Zp2_48_new.iterdb'

def profile_e_info(suffix):
    gene_e = 'profiles_e'+suffix
    gene_i = 'profiles_i'+suffix
    gene_z = 'profiles_z'+suffix
    suffix=gene_e
    f = open(suffix, 'r')
    #prof=f.read()#the read from the profile
    prof = np.genfromtxt(suffix, dtype=float, skip_header=2)
    x_a=prof[:,0]
    x_rho_ref=prof[:,1]
    T=prof[:,2]
    n0=prof[:,3]
    omt=prof[:,4]
    omn=prof[:,5]

    return x_a,x_rho_ref,T,n0,omt,omn 

def profile_i_info(suffix):
    gene_e = 'profiles_e'+suffix
    gene_i = 'profiles_i'+suffix
    gene_z = 'profiles_z'+suffix
    suffix=gene_i
    f = open(suffix, 'r')
    #prof=f.read()#the read from the profile
    prof = np.genfromtxt(suffix, dtype=float, skip_header=2)
    x_a=prof[:,0]
    x_rho_ref=prof[:,1]
    T=prof[:,2]
    n0=prof[:,3]
    omt=prof[:,4]
    omn=prof[:,5]

    return x_a,x_rho_ref,T,n0,omt,omn 

def profile_z_info(suffix):
    gene_e = 'profiles_e'+suffix
    gene_i = 'profiles_i'+suffix
    gene_z = 'profiles_z'+suffix
    suffix=gene_z
    f = open(suffix, 'r')
    #prof=f.read()#the read from the profile
    prof = np.genfromtxt(suffix, dtype=float, skip_header=2)
    x_a=prof[:,0]
    x_rho_ref=prof[:,1]
    T=prof[:,2]
    n0=prof[:,3]
    omt=prof[:,4]
    omn=prof[:,5]

    return x_a,x_rho_ref,T,n0,omt,omn 



#profile_e_info()