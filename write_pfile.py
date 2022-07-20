import numpy as np
from interp import *
import matplotlib.pyplot as pl
from finite_differences import *

#ne -- Electron density
#te -- Electron temperature
#ni -- Ion density
#ti -- Ion temperature
#nb -- Fast ion density
#pb -- Fast ion pressure
#ptot -- Total pressure
#omeg -- Toroidal rotation: VTOR/R
#omegp -- Poloidal rotation: Bt * VPOL / (RBp)
#omgvb -- VxB rotation term in the ExB rotation frequency: OMEG + OMEGP
#omgpp -- Diamagnetic term in the ExB rotation frequency: (P_Carbon)/dpsi / (6*n_Carbon)
#omgeb -- ExB rotation frequency: OMGPP + OMGVB = Er/(RBp)
#er -- Radial electric field from force balance: OMGEB * RBp
#ommvb -- Main ion VXB term of Er/RBp, considered a flux function
#ommpp -- Main ion pressure term of Er/RBp, considered a flux function
#omevb -- Electron VXB term of Er/RBp, considered a flux function
#omepp -- Electron pressure term of Er/RBp, considered a flux function
#kpol -- KPOL=VPOL/Bp : V_vector = KPOL*B_vector + OMGEB * PHI_Vector
#omghb -- Han-Burrell form for the ExB velocity shearing rate: OMGHB = (RBp)**2/Bt * d (Er/RBp)/dpsi
#nz1 -- Density of the 1st impurity species
#vtor1 -- Toroidal velocity of the 1st impurity species
#vpol1 -- Poloidal velocity of the 1st impurity species
#N Z A -- N Z A of ION SPECIES



def get_lists():
    quants = ['ne(10^20/m^3)']
    grads = ['dne/dpsiN']
    quants.append('te(KeV)')
    grads.append('dte/dpsiN')
    quants.append('ni(10^20/m^3)')
    grads.append('dni/dpsiN')
    quants.append('ti(KeV)')
    grads.append('dti/dpsiN')
    quants.append('nb(10^20/m^3)')
    grads.append('dnb/dpsiN')
    quants.append('pb(KPa)')
    grads.append('dpb/dpsiN')
    quants.append('ptot(KPa)') 
    grads.append('dptot/dpsiN')
    quants.append('omeg(kRad/s)') 
    grads.append('domeg/dpsiN')
    quants.append('omgvb(kRad/s)') 
    grads.append('domgvb/dpsiN')
    quants.append('omgpp(kRad/s)') 
    grads.append('domgpp/dpsiN')
    quants.append('er(kV/m)') 
    grads.append('der/dpsiN')
    quants.append('ommvb()') 
    grads.append('dommvb/dpsiN')
    quants.append('ommpp()') 
    grads.append('dommpp/dpsiN')
    quants.append('omevb()') 
    grads.append('domevb/dpsiN')
    quants.append('omepp()') 
    grads.append('domepp/dpsiN')
    quants.append('kpol(km/s/T)') 
    grads.append('dkpol/dpsiN')
    quants.append('omghb()') 
    grads.append('domghb/dpsiN')
    quants.append('nz1(10^20/m^3)') 
    grads.append('dnz1/dpsiN')
    quants.append('vtor1(km/s)') 
    grads.append('dvtor1/dpsiN')
    quants.append('vpol1(km/s)') 
    grads.append('dvpol1/dpsiN')
    return quants, grads

def add_to_pdict(pdict,psiN,field,field_str):
    quants,grads = get_lists()
    if field_str in quants:
        pdict[field_str] = field
        pdict['psinorm_'+field_str] = psiN
        dfield = fd_d1_o4_smoothend(field,psiN)
        qind = quants.index(field_str)  
        pdict[grads[qind]] = dfield
    else:
        print("Error! field_str is not a pfile quantity.")
        stop
    return pdict

def format_psinorm(psi_out):
    temp = str(psi_out)
    ldiff = 8-len(temp)
    if ldiff < 0:
        temp = temp[0:8]
    elif ldiff == 0:
        pass
    else:
        for i in range(ldiff):
            temp += '0'
    return temp

def write_pfile(pdict):

    if 'entry_length' not in pdict:
        if 'ne(10^20/m^3)' in pdict:
            pdict['entry_length'] = len(pdict['ne(10^20/m^3)'])
        else:
            print("Error! pdict does not have sufficient info")
            stop

    if 'file_name' in pdict:
        file_name = pdict['file_name']
    else:
        file_name = 'pfile_default'
    f = open(file_name+'new','w')
    quants,grads = get_lists()
    for i in range(len(quants)):
        if quants[i] in pdict and grads[i] in pdict:
            #f.write(str(pdict['entry_length']) + ' ' + quants[i] + ' ' + grads[i] + '\n') 
            f.write(str(pdict['entry_length']) + ' psinorm '+ quants[i] + ' ' + grads[i] + '\n') 
            for j in range(pdict['entry_length']):
                psi_string = format_psinorm(pdict['psinorm_'+quants[i]][j])
                f.write(psi_string+'   '\
                        +'{:.7}'.format(pdict[quants[i]][j])+'   '\
                        +'{:.7}'.format(pdict[grads[i]][j])+'\n')
    f.write(pdict['species info'])
    f.close()




