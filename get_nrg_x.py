import numpy as np
import sys

def get_nrg0(suffix,nspec=2,ncols=10):

    print( 'Getting data from nrg file')
    time=np.empty(0,dtype='float')
    nrgi=np.empty((0,ncols),dtype='float')
    nrge=np.empty((0,ncols),dtype='float')

    if nspec == 1:
        print( "nspec =", nspec)
        print( "Assume the kinetic species to be electrons.")
    if nspec == 2:
        print( "nspec =", nspec)
        print( "Species order: first: ions, second: electrons.")
    if nspec==3:
        nrgz=np.empty((0,ncols),dtype='float')
        print( "nspec =", nspec)
        print( "Species order: first: ions, second: electrons, third: impurity.")
    if nspec>=4:
        print( "nspec=",nspec)
        sys.exit("nspec must be 1, 2 or 3")

    f=open('nrg'+suffix,'r')
    nrg_in = f.read()
    nrg0 = np.empty((1,ncols))
    nrg_in_lines = nrg_in.split('\n')
    for j in range(len(nrg_in_lines)):
        if nrg_in_lines[j] and j % (nspec+1) == 0:
            time = np.append(time,float(nrg_in_lines[j]))
        elif nrg_in_lines[j] and j % (nspec+1) == 1:
            nline = nrg_in_lines[j].split()
            for i in range(ncols):
                nrg0[0,i] = nline[i]
            if nspec == 1:
               nrge=np.append(nrge,nrg0,axis=0)
            else:
               nrgi=np.append(nrgi,nrg0,axis=0)
        elif nspec>=2 and nrg_in_lines[j] and j % (nspec+1) ==2:
            nline=nrg_in_lines[j].split()
            for i in range(ncols):
                nrg0[0,i]=nline[i]
            nrge=np.append(nrge,nrg0,axis=0)
        elif nspec==3 and nrg_in_lines[j] and j % (nspec+1) ==3:
            nline=nrg_in_lines[j].split()
            for i in range(ncols):
                nrg0[0,i]=nline[i]
            nrgz = np.append(nrgz,nrg0,axis=0)

    if nspec==1:
        return time,nrge
    elif nspec==2:
        return time,nrgi,nrge
    else:
        return time,nrgi,nrge,nrgz
