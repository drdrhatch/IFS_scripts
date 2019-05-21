#!/usr/bin/env python
# -*- coding: utf-8 -*-

from read_write_geometry import * 
from ParIO import * 
import matplotlib.pyplot as plt
import sys
import optparse as op

def smoothWdiff(geomInput,s = 0.3,nmax = 3):

    #plt.plot(geomInput,label='input')
    n = 0
    newArray = np.empty(len(geomInput),dtype='complex128')
    while n < nmax:
        for i in range(len(geomInput)):
            if i == 0 or i == len(geomInput) - 1:
                newArray[i] = geomInput[i]
            else:
                newArray[i] = geomInput[i] + \
                    s * (geomInput[i + 1] + geomInput[i - 1]\
                    - 2. * geomInput[i])
        #plt.plot(newArray,label='after '+str(n+1)+' iterations')
        geomInput = newArray
        n = n + 1
    #plt.legend()
    #plt.show()
    return geomInput
               
def smoothWhypdiff(geomInput,s = 0.003,nmax = 3):

    #plt.plot(geomInput,label='input')
    n = 0
    newArray = np.empty(len(geomInput),dtype='complex128')
    while n < nmax:
        for i in range(len(geomInput)):
            if i == 0 or i == len(geomInput)-1:
                newArray[i] = geomInput[i]
            elif i == 1 or i == len(geomInput)-2:
                newArray[i] = geomInput[i] + \
                    s*(geomInput[i+1]+geomInput[i-1]\
                    - 2.*geomInput[i])
            else:
                newArray[i] = geomInput[i] + \
                    s*(geomInput[i+2]-4.*geomInput[i+1] \
                    + 6.*geomInput[i] - \
                    - 4.*geomInput[i-1] + geomInput[i-2])
        #plt.plot(newArray,label='after '+str(n+1)+' iterations')
        geomInput = newArray
        n = n + 1
    #plt.legend()
    #plt.show()
    return geomInput


parser=op.OptionParser(description='')
options,args=parser.parse_args()
suffix = args[0]
if not suffix =='.dat':
   suffix = '_'+suffix

par = Parameters()
par.Read_Pars('parameters'+suffix)
pars = par.pardict

if 'x_local' in pars:
    if pars['x_local']:
        x_local = True
    else:
        x_local = False
else:
    x_local = True

if 'lilo' in pars:
    if pars['lilo']:
        lilo = True
    else:
        lilo = False
else:
    lilo = False

if x_local:
    geom_type = pars['magn_geometry'][1:-1]
    geom_file = geom_type + suffix
    geom_pars, geom_coeff = read_geometry_local(geom_file)
    gdBdx = geom_coeff['gdBdx']
    gdBdz = geom_coeff['gdBdz']

    if 1 == 1:
        plt.plot(gdBdx)
        plt.title('dBdx (before)')
        plt.show()

        smooth = False
        dBdx_smooth = gdBdx.copy()
        while not smooth:
            selection = int(raw_input\
                ('Smoothed?\n1.True\n2.False\n'))
            if selection == 1:
                smooth = True
                break
            elif selection == 2:
                smoothRate = float(raw_input\
                    ('Enter smoothing rate (< 0.5): '))
                numIteration = float(raw_input\
                    ('Enter number of iteration (>= 1): '))
    
            dBdx_smooth = smoothWdiff(gdBdx,smoothRate,numIteration)
            dBdx_smooth = np.real(dBdx_smooth)

            plt.plot(gdBdx,label='before')
            plt.plot(dBdx_smooth,label='after')
            plt.title('dBdx')
            plt.legend()
            plt.show()

        geom_coeff['gdBdx'] = dBdx_smooth

    if 1 == 1:
        plt.plot(gdBdz)
        plt.title('dBdz (before)')
        plt.show()

        smooth = False
        dBdz_smooth = gdBdz.copy()
        while not smooth:
            selection = int(raw_input\
                ('Smoothed?\n1.True\n2.False\n'))
            if selection == 1:
                smooth = True
                break
            elif selection == 2:
                smoothRate = float(raw_input\
                    ('Enter smoothing rate (< 0.5): '))
                numIteration = float(raw_input\
                    ('Enter number of iteration (>= 1): '))

            dBdz_smooth = smoothWdiff(gdBdz,smoothRate,numIteration)
            dBdz_smooth = np.real(dBdz_smooth)

            plt.plot(gdBdz,label='before')
            plt.plot(dBdz_smooth,label='after')
            plt.title('dBdz')
            plt.legend()
            plt.show()
    
        geom_coeff['gdBdz'] = dBdz_smooth

    file_name = 'smooth_'+geom_type+suffix
    write_tracer_efit_file_local(geom_pars,geom_coeff,file_name)

    if 1 == 1:
        plt.plot(geom_coeff['gdBdx'])
        plt.title('dBdx (after)')
        plt.show()
        plt.plot(geom_coeff['gdBdz'])
        plt.title('dBdz (after)')
        plt.show()

elif lilo:
    geom_type = pars['magn_geometry'][1:-1]
    geom_file = geom_type + suffix
    geom_pars, geom_coeff = read_geometry_global(geom_file)

    if 1 == 1:
        plt.plot(geom_coeff['dBdx'])
        plt.title('dBdx (before)')
        plt.show()

        dBdx = geom_coeff['dBdx']
        smooth = False
        dBdx_smooth = dBdx[:,0].copy()
        while not smooth:
            selection = int(raw_input\
                ('Smoothed?\n1.True\n2.False\n'))
            if selection == 1:
                smooth = True
                break
            elif selection == 2:
                smoothRate = float(raw_input\
                    ('Enter smoothing rate (< 0.5): '))
                numIteration = float(raw_input\
                    ('Enter number of iteration (>= 1): '))

            dBdx_smooth = smoothWdiff\
                (dBdx[:,0],smoothRate,numIteration)
            dBdx_smooth = np.real(dBdx_smooth)

            plt.plot(dBdx[:,0],label='before')
            plt.plot(dBdx_smooth,label='after')
            plt.title('dBdx')
            plt.legend()
            plt.show()

        for i in range(pars['nx0']):
            geom_coeff['dBdx'][:,i] = dBdx_smooth

    if 1 == 1:
        plt.plot(geom_coeff['dBdz'])
        plt.title('dBdz (before)')
        plt.show()

        dBdz = geom_coeff['dBdz']
        smooth = False
        dBdz_smooth = dBdz[:,0].copy()
        while not smooth:
            selection = int(raw_input\
                ('Smoothed?\n1.True\n2.False\n'))
            if selection == 1:
                smooth = True
                break
            elif selection == 2:
                smoothRate = float(raw_input\
                    ('Enter smoothing rate (< 0.5): '))
                numIteration = float(raw_input\
                    ('Enter number of iteration (>= 1): '))

            dBdz_smooth = smoothWdiff(dBdz[:,0],smoothRate,numIteration)
            dBdz_smooth = np.real(dBdz_smooth)

            plt.plot(dBdz[:,0],label='before')
            plt.plot(dBdz_smooth,label='after')
            plt.title('dBdz')
            plt.legend()
            plt.show()

        for i in range(pars['nx0']):
            geom_coeff['dBdz'][:,i] = dBdz_smooth


    file_name = 'smooth_'+geom_type+suffix
    write_tracer_efit_file(geom_pars,geom_coeff,file_name)

    if 1 == 1:
        plt.plot(geom_coeff['dBdx'])
        plt.title('dBdx (after)')
        plt.show()
        plt.plot(geom_coeff['dBdz'])
        plt.title('dBdz (after)')
        plt.show()

else:
    geom_type = pars['magn_geometry'][1:-1]
    geom_file = geom_type + suffix
    geom_pars, geom_coeff = read_geometry_global(geom_file)

    if 1 == 1:
        plt.plot(geom_coeff['dBdx'])
        plt.title('dBdx (before)')
        plt.show()

        dBdx_smooth = np.zeros((pars['nz0'],pars['nx0']),\
            dtype='float')
        for i in range(pars['nx0']):
            dBdx_smooth[:,i] = geom_coeff['dBdx'][:,i].copy()

        smooth = False
        while not smooth:
            selection = int(raw_input\
               ('Smoothed?\n1.True\n2.False\n'))
            if selection == 1:
                smooth = True
                break
            elif selection == 2:
                smoothRate = float(raw_input\
                    ('Enter smoothing rate (< 0.5): '))
                numIteration = float(raw_input\
                    ('Enter number of iteration (>= 1): '))
            for i in range(pars['nx0']):
                dBdx = geom_coeff['dBdx'][:,i].copy()
                dBdx_smooth[:,i] = np.real(smoothWdiff(dBdx,smoothRate,numIteration))
                if i == pars['nx0']*4/5:
                    plt.plot(dBdx,label='before')
                    plt.plot(dBdx_smooth[:,i],label='after')
                    plt.title('dBdx')
                    plt.legend()
                    plt.show()

        geom_coeff['dBdx'] = dBdx_smooth

    if 1 == 1:
        plt.plot(geom_coeff['dBdz'])
        plt.title('dBdz (before)')
        plt.show()

        dBdz_smooth = np.zeros((pars['nz0'],pars['nx0']),\
            dtype='float')
        for i in range(pars['nx0']):
            dBdz_smooth[:,i] = geom_coeff['dBdz'][:,i].copy()

        smooth = False
        while not smooth:
            selection = int(raw_input\
                ('Smoothed?\n1.True\n2.False\n'))
            if selection == 1:
                smooth = True
                break
            elif selection == 2:
                smoothRate = float(raw_input\
                    ('Enter smoothing rate (< 0.5): '))
                numIteration = float(raw_input\
                    ('Enter number of iteration (>= 1): '))

            for i in range(pars['nx0']):
                dBdz = geom_coeff['dBdz'][:,i].copy()
                dBdz_smooth[:,i] = np.real(smoothWdiff(dBdz,smoothRate,numIteration))
                if i == pars['nx0']*4/5:
                    plt.plot(dBdz,label='before')
                    plt.plot(dBdz_smooth[:,i],label='after')
                    plt.title('dBdz')
                    plt.legend()
                    plt.show()

        geom_coeff['dBdz'] = dBdz_smooth

    file_name = 'smooth_'+geom_type+suffix
    write_tracer_efit_file(geom_pars,geom_coeff,file_name)

    if 1 == 1:
        plt.plot(geom_coeff['dBdx'])
        plt.title('dBdx (after)')
        plt.show()
        plt.plot(geom_coeff['dBdz'])
        plt.title('dBdz (after)')
        plt.show()

