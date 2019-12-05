#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import traceback

import numpy as npy
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdfh

from scipy.integrate import simps
from scipy.interpolate import interp1d,interp2d
from scipy.interpolate import CubicSpline,RectBivariateSpline

import read_iterdb


def bisection(fx,xmin,xmax,root=0.0,Nmax=100,eps=1.0e-16):
    for it in range(Nmax):
        x = (xmax+xmin)/2.0
        if abs(fx(x)-root) < eps: break
        if   (fx(xmax)-root)*(fx(x)-root) < 0.0: xmin = x
        else:                                    xmax = x
    return x

def rk4(f1,f2,czt,slns,stp):
    k01 = f1(slns[0],slns[1])
    k11 = f2(slns[0],slns[1])
    k02 = f1(slns[0]+stp*k01/2.0,slns[1]+stp*k11/2.0)
    k12 = f2(slns[0]+stp*k01/2.0,slns[1]+stp*k11/2.0)
    k03 = f1(slns[0]+stp*k02/2.0,slns[1]+stp*k12/2.0)
    k13 = f2(slns[0]+stp*k02/2.0,slns[1]+stp*k12/2.0)
    k04 = f1(slns[0]+stp*k03,slns[1]+stp*k13)
    k14 = f2(slns[0]+stp*k03,slns[1]+stp*k13)
    slns[0] = slns[0]+stp*(k01+2*k02+2*k03+k04)/6.0
    slns[1] = slns[1]+stp*(k11+2*k12+2*k13+k14)/6.0
    return slns

def rk5(f1,f2,czt,slns,stp):
    k01 = f1(slns[0],slns[1])
    k11 = f2(slns[0],slns[1])

    k02 = f1(slns[0]+stp*k01/4.0,slns[1]+stp*k11/4.0)
    k12 = f2(slns[0]+stp*k01/4.0,slns[1]+stp*k11/4.0)

    k03 = f1(slns[0]+stp*(3.0*k01+9.0*k02)/32.0,slns[1]+stp*(3.0*k11+9.0*k12)/32.0)
    k13 = f2(slns[0]+stp*(3.0*k01+9.0*k02)/32.0,slns[1]+stp*(3.0*k11+9.0*k12)/32.0)

    k04 = f1(slns[0]+stp*(1932*k01-7200.0*k02+7296.0*k03)/2197.0,slns[1]+stp*(1932*k11-7200.0*k12+7296.0*k13)/2197.0)
    k14 = f2(slns[0]+stp*(1932*k01-7200.0*k02+7296.0*k03)/2197.0,slns[1]+stp*(1932*k11-7200.0*k12+7296.0*k13)/2197.0)

    k05 = f1(slns[0]+stp*(439*k01/216.0-8.0*k02+3680.0*k03/513.0-845.0*k04/4104.0),slns[1]+stp*(439*k11/216.0-8.0*k12+3680.0*k13/513.0-845.0*k14/4104.0))
    k15 = f2(slns[0]+stp*(439*k01/216.0-8.0*k02+3680.0*k03/513.0-845.0*k04/4104.0),slns[1]+stp*(439*k11/216.0-8.0*k12+3680.0*k13/513.0-845.0*k14/4104.0))

    k06 = f1(slns[0]+stp*(-8.0*k01/27.0+2.0*k02-3544.0*k03/2565.0+1859.0*k04/4104.0-11.0*k05/40.0),slns[1]+stp*(-8.0*k11/27.0+2.0*k12-3544.0*k13/2565.0+1859.0*k14/4104.0-11.0*k15/40.0))
    k16 = f2(slns[0]+stp*(-8.0*k01/27.0+2.0*k02-3544.0*k03/2565.0+1859.0*k04/4104.0-11.0*k05/40.0),slns[1]+stp*(-8.0*k11/27.0+2.0*k12-3544.0*k13/2565.0+1859.0*k14/4104.0-11.0*k15/40.0))

    slns[0] = slns[0]+stp*(16.0*k01/135.0+6656.0*k03/12825.0+28561.0*k04/56430.0-9.0*k05/50.0+2*k06/55.0)
    slns[1] = slns[1]+stp*(16.0*k11/135.0+6656.0*k13/12825.0+28561.0*k14/56430.0-9.0*k15/50.0+2*k16/55.0)
    return slns

def read_iterdb_file(filename):
    '''
    This Code is Written by: David R. Hatch
    It reads the iterdb file and returns three dictionaries,
    each diectionary has five quantities:
    electron density (NE) and temperature (TE),
    ion density (NM1) and temperatures (TI),
    impurity density (NM2), if any,
    rotational velocity (VROT).
    The three dictionaries provide the toroidal coordinate (rhotor),
    profiles, and units for each quantity.
    '''
    f=open(filename,'r')
    data_in=f.read()
    data_linesplit=data_in.split('\n')

    keep_going=1
    i=0
    while keep_going:
        test=re.search(';-# OF X PTS',data_linesplit[i])
        if test:
            num=data_linesplit[i].split()[0]
            num=float(num)
            num=int(num)
            keep_going=(1==2)
        if i == len(data_linesplit):
            keep_going=(1==2)
        i=i+1

    lnum=0
    try_again=1
    prof_out = {}
    rhot_out = {}
    units_out = {}
    while try_again:
        lnum,try_again,quantity,units,rhot,arr=get_next(data_linesplit,lnum,num)
        prof_out[quantity]=arr
        units_out[quantity]=units
        rhot_out[quantity]=rhot
    return rhot_out,prof_out,units_out

def get_next(data_linesplit,lnum,num):
    sec_num_lines = num/6
    if num % 6 != 0:
        sec_num_lines += 1
    keep_going=1
    while keep_going:
        test=re.search('-DEPENDENT VARIABLE LABEL',data_linesplit[lnum])
        if test :
            quantity=data_linesplit[lnum].split()[0]
            units=data_linesplit[lnum].split()[1]
        test=re.search('DATA FOLLOW',data_linesplit[lnum])
        if test:
            keep_going=(1==2)
        lnum=lnum+1

    rhot=npy.empty(0)
    lnum0 = lnum
    for j in range(lnum0,lnum0+int(sec_num_lines)):
        for k in range(6):
            str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
            if(str_temp):
                temp=npy.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                rhot=npy.append(rhot,temp)
        lnum=lnum+1
    lnum=lnum+1

    arr=npy.empty(0)
    lnum0 = lnum
    for j in range(lnum0,lnum0+int(sec_num_lines)):
        for k in range(6):
            str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
            if(str_temp):
                temp=npy.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                arr=npy.append(arr,temp)
        lnum=lnum+1

    lnum_out=lnum
    try_again=1
    if len(data_linesplit)-lnum < 10:
        try_again=False
    return lnum_out, try_again,quantity,units,rhot,arr

def plot_iterdb(iterdbdata,reportpath=''):
    if len(reportpath)==0:
       reportpath = './'
    elif 'report' not in reportpath:
       if    reportpath[-1] != "/": reportpath += "/report/"
       else:                        reportpath += "report/"
    rhotor,profiles,units = iterdbdata

    figure = plt.figure('Density')
    axhand = figure.add_subplot(1,1,1)
    axhand.plot(rhotor['NE'],profiles['NE'],label='$n_e$')
    axhand.plot(rhotor['NM1'],profiles['NM1'],label='$n_i$')
    if len(rhotor.keys()) >= 6:
       axhand.plot(rhotor['NM2'],profiles['NM2'],label='$n_z$')
    axhand.set_title('Density Profiles')
    axhand.set_xlabel('$\\rho_{tor}$')
    axhand.set_ylabel('Density('+units['NE']+')')
    axhand.legend()
    figure.savefig(reportpath+'iterdb_density.png') 
    plt.close(figure)

    figure = plt.figure('Temperature')
    axhand = figure.add_subplot(1,1,1)
    axhand.plot(rhotor['TE'],profiles['TE'],label='$T_e$')
    axhand.plot(rhotor['TI'],profiles['TI'],label='$T_i$')
    axhand.set_title('Temperature Profiles')
    axhand.set_xlabel('$\\rho_{tor}$')
    axhand.set_ylabel('Temperature('+units['TE']+')')
    axhand.legend()
    figure.savefig(reportpath+'iterdb_temperature.png') 
    plt.close(figure)

    return 1

def read_profiles_file(pfpath,setParam={}):
   #Developed by Ehab Hassan on 2019-02-17
    from scipy.interpolate import CubicSpline
    if not os.path.isfile(pfpath):
       print('Fatal: file %s not found.' % pfpath)
       sys.exit()

    ofh = open(pfpath,'r')

    profiles = {}
    units    = {}
    while True:
          recs = ofh.readline().split()
          if   len(recs)>4:
               nrec = int(recs[0])
               ary0=npy.zeros(nrec)
               ary1=npy.zeros(nrec)
               ary2=npy.zeros(nrec)
               var0 = str(recs[1]).lower()
               var1 = str(recs[2]).lower()
               var2 = str(recs[3]).lower()
               for i in range(nrec):
                   recs = ofh.readline().split()
                   ary0[i] = float(recs[0])
                   ary1[i] = float(recs[1])
                   ary2[i] = float(recs[2])
               profiles[var0]=ary0
               profiles[var1]=ary1
               profiles[var2]=ary2
          elif len(recs)==4:
               nrec = int(recs[0])
               ary0=npy.zeros(nrec)
               ary1=npy.zeros(nrec)
               ary2=npy.zeros(nrec)
               var0 = str(recs[1]).lower()
               temp = str(recs[2]).lower()
               var1 = temp[:temp.index("(")]
               if   var1.strip() in ['ne','ni','nb','nz1']:
                    powr = int(temp[temp.index("^")+1:temp.index("/")])
                    unit = 'm^{-3}'
               elif var1.strip() in ['te','ti']:
                    unit = 'eV'
               elif var1.strip() in ['pb','ptot']:
                    unit = 'Pa'
               elif var1.strip() in ['vtor1','vpol1']:
                    unit = 'm/s'
               else:
                    unit = temp[temp.index("(")+1:temp.index(")")]
               var2 = str(recs[3]).lower()
               for i in range(nrec):
                   recs    = ofh.readline().split()
                   ary0[i] = float(recs[0])
                   ary1[i] = float(recs[1])
                   ary2[i] = float(recs[2])
               if var0 in profiles.keys():
                  CS = CubicSpline(ary0,ary1)
                  units[var1]    = unit
                  profiles[var1] = CS(profiles[var0])
                  CS = CubicSpline(ary0,ary2)
                  profiles[var2] = CS(profiles[var0])
               else:
                  profiles[var0] = ary0
                  units[var1]    = unit
                  profiles[var1] = ary1
                  profiles[var2] = ary2
               if   var1.strip() in ['ne','ni','nb','nz1']:
                    profiles[var1] *= 10.0**powr
               elif var1.strip() in ['te','ti','pb','ptot','vtor1','vpol1']:
                    profiles[var1] *= 1.0e3
          else:
               break
    ofh.close()

    if 'rhotor' in setParam:
       if 'eqdskfpath' not in setParam:
          eqdskfpath = raw_input('Path to EQDSK file: ')
       elif 'eqdskfpath' in setParam:
          eqdskfpath = setParam['eqdskfpath']

       eqdskdata = read_efit_file(eqdskfpath)

       qpsifn  = interp1d(eqdskdata['PSIN'],eqdskdata['qpsi'])
       qpsi    = qpsifn(profiles['psinorm'])
       qpsifn  = interp1d(profiles['psinorm'],qpsi)
       psinorm = npy.linspace(profiles['psinorm'][0],profiles['psinorm'][-1],setParam['rhotor'])
       qpsi    = qpsifn(psinorm)
       phinorm = npy.zeros_like(psinorm)
       for i in range(1,npy.size(qpsi)):
           x = psinorm[:i+1]
           y = qpsi[:i+1]
           phinorm[i] = npy.trapz(y,x)
       profiles['rhotor'] = npy.sqrt(phinorm)
       
    return profiles,units

def plot_profiles(pfpath):
    from matplotlib.backends.backend_pdf import PdfPages
    profiles = read_profiles(pfpath)
    nvars = sorted(profiles.keys())
    print('Plotting profiles in: %s ...' % pfpath)
    figs  = PdfPages('profiles.pdf')
    for i in nvars:
        if i in ['n','z','a','psinorm']: continue
        figname = plt.figure(i)
        plt.plot(profiles['psinorm'],profiles[i])
        plt.xlabel('$\Psi$')
        plt.ylabel(i)
        figs.savefig(figname)
    figs.close()
    return profiles

def read_efit_file(eqdskfpath,setParam={}):
   #Developed by Ehab Hassan on 2019-02-27
    if os.path.isfile(eqdskfpath) == False:
       errorFunc = traceback.extract_stack(limit=2)[-2][3]
       errorLine = traceback.extract_stack(limit=2)[-2][1]
       errorFile = traceback.extract_stack(limit=2)[-2][2]
       errMSG    = 'Call %s line %5d in file %s Failed.\n'
       errMSG   += 'Fatal: file %s not found.'
       raise IOError(errMSG %(errorFunc,errorLine,errorFile,eqdskfpath))
    
    ofh = open(eqdskfpath,'r')
    eqdskdata = {}
    cline = ofh.readline()
    eqdskdata['idum']   = int(cline[48:52])
    eqdskdata['RDIM']   = int(cline[52:56])
    eqdskdata['ZDIM']   = int(cline[56:61])
    cline = ofh.readline()
    eqdskdata['RLEN']   = float(cline[0:16])
    eqdskdata['ZLEN']   = float(cline[16:32])
    eqdskdata['RCTR']   = float(cline[32:48])
    eqdskdata['RLFT']   = float(cline[48:64])
    eqdskdata['ZMID']   = float(cline[64:80])
    cline = ofh.readline()
    eqdskdata['RMAX']   = float(cline[0:16])
    eqdskdata['ZMAX']   = float(cline[16:32])
    eqdskdata['PSIMAX'] = float(cline[32:48])
    eqdskdata['PSIBND'] = float(cline[48:64])
    eqdskdata['BCTR']   = float(cline[64:80])
    cline = ofh.readline()
    eqdskdata['CURNT']  = float(cline[0:16])
    eqdskdata['PSIMAX'] = float(cline[16:32])
    eqdskdata['XDUM']   = float(cline[32:48])
    eqdskdata['RMAX']   = float(cline[48:64])
    eqdskdata['XDUM']   = float(cline[64:80])
    cline = ofh.readline()
    eqdskdata['ZMAX']   = float(cline[0:16])
    eqdskdata['XDUM']   = float(cline[16:32])
    eqdskdata['PSIBND'] = float(cline[32:48])
    eqdskdata['XDUM']   = float(cline[48:64])
    eqdskdata['XDUM']   = float(cline[64:80])

    nlines1D = int(npy.ceil(eqdskdata['RDIM']/5.0))

    eqdskdata['fpol'] = npy.zeros(eqdskdata['RDIM'])
    for iline in range(nlines1D):
        cline = ofh.readline()
        try:
            eqdskdata['fpol'][iline*5+0] = float(cline[0:16])
            eqdskdata['fpol'][iline*5+1] = float(cline[16:32])
            eqdskdata['fpol'][iline*5+2] = float(cline[32:48])
            eqdskdata['fpol'][iline*5+3] = float(cline[48:64])
            eqdskdata['fpol'][iline*5+4] = float(cline[64:80])
        except:
            error = 'empty records'

    eqdskdata['pressure'] = npy.zeros(eqdskdata['RDIM'])
    for iline in range(nlines1D):
        cline = ofh.readline()
        try:
            eqdskdata['pressure'][iline*5+0] = float(cline[0:16])
            eqdskdata['pressure'][iline*5+1] = float(cline[16:32])
            eqdskdata['pressure'][iline*5+2] = float(cline[32:48])
            eqdskdata['pressure'][iline*5+3] = float(cline[48:64])
            eqdskdata['pressure'][iline*5+4] = float(cline[64:80])
        except:
            error = 'empty records'

    eqdskdata['ffprime'] = npy.zeros(eqdskdata['RDIM'])
    for iline in range(nlines1D):
        cline = ofh.readline()
        try:
            eqdskdata['ffprime'][iline*5+0] = float(cline[0:16])
            eqdskdata['ffprime'][iline*5+1] = float(cline[16:32])
            eqdskdata['ffprime'][iline*5+2] = float(cline[32:48])
            eqdskdata['ffprime'][iline*5+3] = float(cline[48:64])
            eqdskdata['ffprime'][iline*5+4] = float(cline[64:80])
        except:
            error = 'empty records'

    eqdskdata['pprime'] = npy.zeros(eqdskdata['RDIM'])
    for iline in range(nlines1D):
        cline = ofh.readline()
        try:
            eqdskdata['pprime'][iline*5+0] = float(cline[0:16])
            eqdskdata['pprime'][iline*5+1] = float(cline[16:32])
            eqdskdata['pprime'][iline*5+2] = float(cline[32:48])
            eqdskdata['pprime'][iline*5+3] = float(cline[48:64])
            eqdskdata['pprime'][iline*5+4] = float(cline[64:80])
        except:
            error = 'empty records'

    nlines2D = int(npy.ceil(eqdskdata['RDIM']*eqdskdata['ZDIM']/5.0))

    eqdskdata['psiRZ'] = npy.zeros(eqdskdata['RDIM']*eqdskdata['ZDIM'])
    for iline in range(nlines2D):
        cline = ofh.readline()
        try:
            eqdskdata['psiRZ'][iline*5+0] = float(cline[0:16])
            eqdskdata['psiRZ'][iline*5+1] = float(cline[16:32])
            eqdskdata['psiRZ'][iline*5+2] = float(cline[32:48])
            eqdskdata['psiRZ'][iline*5+3] = float(cline[48:64])
            eqdskdata['psiRZ'][iline*5+4] = float(cline[64:80])
        except:
            error = 'empty records'
    eqdskdata['psiRZ'] = npy.reshape(eqdskdata['psiRZ'],(eqdskdata['ZDIM'],eqdskdata['RDIM']))

    eqdskdata['qpsi'] = npy.zeros(eqdskdata['RDIM'])
    for iline in range(nlines1D):
        cline = ofh.readline()
        try:
            eqdskdata['qpsi'][iline*5+0] = float(cline[0:16])
            eqdskdata['qpsi'][iline*5+1] = float(cline[16:32])
            eqdskdata['qpsi'][iline*5+2] = float(cline[32:48])
            eqdskdata['qpsi'][iline*5+3] = float(cline[48:64])
            eqdskdata['qpsi'][iline*5+4] = float(cline[64:80])
        except:
            error = 'empty records'
 
    cline = ofh.readline()
    eqdskdata['nbound'] = int(cline[0:5])
    eqdskdata['nlimit'] = int(cline[5:10])

    if eqdskdata['nbound'] > 0:
       nlines1D = int(npy.ceil(2*eqdskdata['nbound']/5.0))

       Ary1D = npy.zeros(2*eqdskdata['nbound'])
       for iline in range(nlines1D):
           cline = ofh.readline()
           try:
               Ary1D[iline*5+0] = float(cline[0:16])
               Ary1D[iline*5+1] = float(cline[16:32])
               Ary1D[iline*5+2] = float(cline[32:48])
               Ary1D[iline*5+3] = float(cline[48:64])
               Ary1D[iline*5+4] = float(cline[64:80])
           except:
               error = 'empty records'

       eqdskdata['rbound'] = Ary1D[0::2]
       eqdskdata['zbound'] = Ary1D[1::2]


    if eqdskdata['nlimit'] > 0:
       nlines1D = int(npy.ceil(2*eqdskdata['nlimit']/5.0))

       Ary1D = npy.zeros(2*eqdskdata['nlimit'])
       for iline in range(nlines1D):
           cline = ofh.readline()
           try:
               Ary1D[iline*5+0] = float(cline[0:16])
               Ary1D[iline*5+1] = float(cline[16:32])
               Ary1D[iline*5+2] = float(cline[32:48])
               Ary1D[iline*5+3] = float(cline[48:64])
               Ary1D[iline*5+4] = float(cline[64:80])
           except:
               error = 'empty records'

       eqdskdata['rlimit'] = Ary1D[0::2]
       eqdskdata['zlimit'] = Ary1D[1::2]

 
    eqdskdata['ZR1D']  = npy.arange(eqdskdata['ZDIM'],dtype=float)*eqdskdata['ZLEN']/(eqdskdata['ZDIM']-1.0)
    eqdskdata['ZR1D'] += eqdskdata['ZMID']-eqdskdata['ZMID']/2.0

    eqdskdata['RR1D']  = npy.arange(eqdskdata['RDIM'],dtype=float)*eqdskdata['RLEN']/(eqdskdata['RDIM']-1.0)
    eqdskdata['RR1D'] += eqdskdata['RLFT']

    eqdskdata['psiRZ'] = (eqdskdata['psiRZ']-eqdskdata['PSIMAX'])/(eqdskdata['PSIBND']-eqdskdata['PSIMAX'])

    eqdskdata['PSI']    = (eqdskdata['PSIBND']-eqdskdata['PSIMAX'])*npy.arange(eqdskdata['RDIM'])/(eqdskdata['RDIM']-1.0)
    eqdskdata['PSIN']   = (eqdskdata['PSI']-eqdskdata['PSI'][0])/(eqdskdata['PSI'][-1]-eqdskdata['PSI'][0])
    eqdskdata['rhopsi'] = npy.sqrt(eqdskdata['PSIN'])

    extendPSI    = npy.linspace(eqdskdata['PSI'][0],eqdskdata['PSI'][-1],10*npy.size(eqdskdata['PSI']))
    extendPHI    = npy.empty_like(extendPSI)
    extendPHI[0] = 0.0
    qfunc        = CubicSpline(eqdskdata['PSI'],eqdskdata['qpsi'])
    for i in range(1,npy.size(extendPSI)):
        x           = extendPSI[:i+1]
        y           = qfunc(x)
        extendPHI[i]= npy.trapz(y,x)

    eqdskdata['PHI'] = npy.empty_like(eqdskdata['PSI'])
    phifunc          = CubicSpline(extendPSI,extendPHI)
    for i in range(npy.size(eqdskdata['PSI'])):
        eqdskdata['PHI'][i] = phifunc(eqdskdata['PSI'][i])

    eqdskdata['PHIN']   = (eqdskdata['PHI']-eqdskdata['PHI'][0])/(eqdskdata['PHI'][-1]-eqdskdata['PHI'][0])
    eqdskdata['rhotor'] = npy.sqrt(eqdskdata['PHIN'])

    return eqdskdata

def plot_eqdsk(eqdskdata,reportpath=''):
   #Developed by Ehab Hassan on 2019-03-11
    if len(reportpath)==0:
       reportpath = './'
    else:
       if    reportpath[-1] != "/": reportpath += "/report/"
       else:                        reportpath += "report/"

    figure = plt.figure('Magnetic Surface Boundary')
    axhand = figure.add_subplot(1,1,1)
    axhand.plot(eqdskdata['rbound'],eqdskdata['zbound'])
    axhand.set_title('Magnetic Surface Boundary')
    axhand.set_xlabel('r')
    axhand.set_ylabel('z')
    figure.savefig(reportpath+'eqdsk_magsurfbound.png')
    plt.close(figure)

    figure = plt.figure('Safety Factor')
    axhand = figure.add_subplot(1,1,1)
    axhand.plot(eqdskdata['rhopsi'],eqdskdata['qpsi'])
    axhand.set_title('Safety Factor')
    axhand.set_xlabel('$\\rho_{\\psi}$')
    axhand.set_ylabel('q')
    figure.savefig(reportpath+'eqdsk_safetyfactor.png')
    plt.close(figure)

    figure = plt.figure('Pressure Gradient')
    axhand = figure.add_subplot(1,1,1)
    axhand.plot(eqdskdata['rhopsi'],eqdskdata['pprime'])
    axhand.set_title('Pressure Gradient')
    axhand.set_xlabel('$\\rho_{\\psi}$')
    axhand.set_ylabel("P'")
    figure.savefig(reportpath+'eqdsk_pprime.png')
    plt.close(figure)

    figure = plt.figure('Current Flux')
    axhand = figure.add_subplot(1,1,1)
    axhand.plot(eqdskdata['rhopsi'],eqdskdata['ffprime'])
    axhand.set_title('Cuurent Flux')
    axhand.set_xlabel('$\\rho_{\\psi}$')
    axhand.set_ylabel("FF'")
    figure.savefig(reportpath+'eqdsk_ffprime.png')
    plt.close(figure)

    return 1

def profile_mapping(psi,phi,profpsi=[],profphi=[]):
    if   any(profpsi):
         profpsifn=interp1d(psi,profpsi,kind='linear')
         phipsifn =interp1d(psi,phi,kind='linear')
         profphi  =profpsifn(phipsifn(psi))
    elif any(profphi):
         profphifn=interp1d(phi,profphi,kind='linear')
         psiphifn =interp1d(phi,psi,kind='linear')
         profpsi  =profphifn(psiphifn(phi))

def psi2phi(q,psi):
    qpsifn  = interp1d(psi,q)
    psinorm = npy.linspace(psi[0],psi[-1],10*len(psi))
    qpsi    = qpsifn(psinorm)
    phinorm = npy.zeros_like(psinorm)
    for i in range(1,npy.size(qpsi)):
        x = psinorm[1:i+1]
        y = qpsi[1:i+1]
        phinorm[i] = npy.trapz(y,x)
    phinorm  = (phinorm-phinorm[0])/(phinorm[-1]-phinorm[0])
    phipsifn = interp1d(psinorm,phinorm)
    phi      = phipsifn(psi)
    return phi

def phi2psi(q,phi):
    qphifn  = interp1d(phi,q)
    phinorm = npy.linspace(phi[0],phi[-1],10*len(phi))
    qphi    = qphifn(phinorm)
    psinorm = npy.zeros_like(phinorm)
    for i in range(1,npy.size(qphi)):
        x = phinorm[1:i+1]
        y = 1./qphi[1:i+1]
        psinorm[i] = npy.trapz(y,x)
    psinorm  = (psinorm-psinorm[0])/(psinorm[-1]-psinorm[0])
    psiphifn = interp1d(phinorm,psinorm)
    psi      = psiphifn(phi)
    return psi


def findmonotonic(A,kind="increasing"):
    if kind.lower()=="increasing":
       bgnloc=0
       endloc=len(A)-1
       for i in range(len(A)):
           if A[i+1]>A[i]:
              bgnloc=i
              break
       for i in range(bgnloc,len(A)-1):
           if A[i+1]<A[i]:
              endloc=i
              break
    return bgnloc,endloc


def magsurf_contours(eqdskfpath):
    eqdskdata = read_efit_file(eqdskfpath.strip())

    minZ = min(eqdskdata['zbound'])+eqdskdata['ZLEN']/2
    maxZ = max(eqdskdata['zbound'])+eqdskdata['ZLEN']/2
    minR = min(eqdskdata['rbound'])
    maxR = max(eqdskdata['rbound'])

    R2D,Z2D = npy.meshgrid(eqdskdata['RR1D'],eqdskdata['ZR1D'])

    fig,ax = plt.subplots(1)
    cs = plt.contour(R2D,Z2D,eqdskdata['psiRZ'],1000)
    plt.close(fig)
    lines = []
    for j in range(0,len(cs.levels)):
      if cs.levels[j]>1.0:
         break
      else:
         for line in cs.collections[j].get_paths():
           lines.append(line.vertices)
           if any(lines[-1][:,0] <= minR):  continue
           if any(lines[-1][:,0] >= maxR):  continue
           if any(lines[-1][:,1] <= minZ):  continue
           if any(lines[-1][:,1] >= maxZ):  continue
           rcval = lines[-1][:,0]
           zcval = lines[-1][:,1]
    return rcval,zcval-eqdskdata['ZLEN']/2

def magsurf_solvflines(eqdskfpath='',eqdskdata={},psi=1.0,eps=1.0e-16):
    if not eqdskdata.keys():
       if eqdskfpath:
          eqdskdata = read_efit_file(eqdskfpath.strip())
       else:
          print('FATAL: EQDSK FILE NOT FOUND. EXIT!')
          sys.exit()

    qpsifn = interp1d(eqdskdata['PSIN'],eqdskdata['qpsi'],kind="linear")

    minZ = min(eqdskdata['zbound'])+eqdskdata['ZLEN']/2
    maxZ = max(eqdskdata['zbound'])+eqdskdata['ZLEN']/2
    minR = min(eqdskdata['rbound'])
    maxR = max(eqdskdata['rbound'])

    psiRZ    = interp2d(eqdskdata['RR1D'],eqdskdata['ZR1D'],eqdskdata['psiRZ'],kind="cubic")
    dpsiRZdR = psiRZ.__call__(x=eqdskdata['RR1D'],y=eqdskdata['ZR1D'],dx=1)
    dpsiRZdZ = psiRZ.__call__(x=eqdskdata['RR1D'],y=eqdskdata['ZR1D'],dy=1)
    BR       =-dpsiRZdZ/eqdskdata['BCTR']/eqdskdata['RCTR']
    BZ       = dpsiRZdR/eqdskdata['BCTR']/eqdskdata['RCTR']
    f1       = interp2d(eqdskdata['RR1D'],eqdskdata['ZR1D'],BR,kind="cubic")
    f2       = interp2d(eqdskdata['RR1D'],eqdskdata['ZR1D'],BZ,kind="cubic")

    cs_psi = 1.0
    for k in range(npy.size(eqdskdata['ZR1D'])):
        z    = eqdskdata['ZR1D'][npy.size(eqdskdata['ZR1D'])-k-1]
        psir = interp1d(eqdskdata['RR1D'],psiRZ(eqdskdata['RR1D'],z),kind="cubic")
        r    = bisection(psir,eqdskdata['RR1D'][0],eqdskdata['RR1D'][-1],root=psi,eps=eps)
        if r <= minR or r >= maxR:  continue
        if z <= minZ or z >= maxZ:  continue
        if abs(psiRZ(r,z)-psi)<eps or abs(psiRZ(r,z)-psi)<cs_psi:
           cs_psi = abs(psiRZ(r,z)-psi)
           cs_r   = r
           cs_z   = z

    npts = 1024

    rsln = []
    zsln = []
    slns = npy.array([cs_r,cs_z])

    if type(psiRZ(cs_r,cs_z))==float:
       maxzta = 4.0*npy.pi*eqdskdata['RLEN']*eqdskdata['ZLEN']*psiRZ(cs_r,cs_z)
    else:
       maxzta = 4.0*npy.pi*eqdskdata['RLEN']*eqdskdata['ZLEN']*max(psiRZ(cs_r,cs_z))
    zta = npy.linspace(0.0,maxzta,npts)
    stp = npy.diff(zta)[0]
    rflipped = 0; rflag = False
    zflipped = 0; zflag = False
    nzta = npy.size(zta)
    for j in range(nzta):
       slns = rk5(f1,f2,zta[j],slns,stp)
       if rflipped >= 2: rflag = True
       if zflipped >= 2: zflag = True
       if (rflag and zflag) and (abs(rref-slns[0])<0.05 and abs(zref-slns[1])<0.05):
          break
       rsln.append(slns[0])
       zsln.append(slns[1])
       if j==0:
          rref = eqdskdata['rbound'][npy.argmin(abs(rsln[0]-eqdskdata['rbound']))]
          zref = eqdskdata['zbound'][npy.argmin(abs(zsln[0]-eqdskdata['zbound']))]
       elif j==1:
          rrefsign = npy.sign(rsln[-1]-rsln[-2])
          zrefsign = npy.sign(zsln[-1]-zsln[-2])
       elif j>2:
          rcrntsign = npy.sign(rsln[-1]-rsln[-2])
          zcrntsign = npy.sign(zsln[-1]-zsln[-2])
          if rrefsign ==-rcrntsign: rflipped += 1; rrefsign = rcrntsign 
          if zrefsign ==-zcrntsign: zflipped += 1; zrefsign = zcrntsign 

    return npy.array(rsln),npy.array(zsln)-eqdskdata['ZLEN']/2

def magsurf_searching(eqdskfpath='',eqdskdata={},psi=1.0,eps=1.0e-6):
    if not eqdskdata.keys():
       if eqdskfpath:
          eqdskdata = read_efit_file(eqdskfpath.strip())
       else:
          print('FATAL: EQDSK FILE NOT FOUND. EXIT!')
          sys.exit()

    minZ = min(eqdskdata['zbound'])+eqdskdata['ZLEN']/2
    maxZ = max(eqdskdata['zbound'])+eqdskdata['ZLEN']/2
    minR = min(eqdskdata['rbound'])
    maxR = max(eqdskdata['rbound'])

    psiRZ    = interp2d(eqdskdata['RR1D'],eqdskdata['ZR1D'],eqdskdata['psiRZ'],kind="cubic")

    rbound = []
    zbound = []

#  #Z1D = npy.linspace(eqdskdata['ZR1D'][0],eqdskdata['ZR1D'][-1],1024)
#  #R1D = npy.linspace(minR,maxR,1024)
#  #R1D = npy.linspace(eqdskdata['RR1D'][0],eqdskdata['RR1D'][-1],1024)

#   Z1D = npy.linspace(minZ,maxZ,1024)
#  #medR= (minR+maxR)/2.0
#  #R1D = npy.linspace(medR,maxR,1024)
#   R1D = npy.linspace(eqdskdata['RR1D'][0],eqdskdata['RR1D'][npy.size(eqdskdata['RR1D'])/2] ,1024)
#   for k in range(npy.size(Z1D)):
#       z    = Z1D[k]
#       psir = interp1d(R1D,psiRZ(R1D,z),kind="cubic")
#       r    = bisection(psir,R1D[0],R1D[-1],root=psi,eps=eps)
#       if r <= minR or r >= maxR:  continue
#       if z <= minZ or z >= maxZ:  continue
#       if abs(psiRZ(r,z)-psi)<eps:
#          rbound.append(r)
#          zbound.append(z)

#  #R1D = npy.linspace(minR,medR,1024)
#   R1D = npy.linspace(eqdskdata['RR1D'][npy.size(eqdskdata['RR1D'])/2],eqdskdata['RR1D'][-1],1024)
#   for k in range(npy.size(Z1D)):
#       z    = Z1D[k]
#       psir = interp1d(R1D,psiRZ(R1D,z),kind="cubic")
#       r    = bisection(psir,R1D[0],R1D[-1],root=psi,eps=eps)
#       if r <= minR or r >= maxR:  continue
#       if z <= minZ or z >= maxZ:  continue
#       if abs(psiRZ(r,z)-psi)<eps:
#          rbound.append(r)
#          zbound.append(z)

    ZMAX = eqdskdata['ZMAX']
    RMAX = eqdskdata['RMAX']
    dens = 512
    for i in range(npy.size(eqdskdata['rbound'])):
        rho = npy.sqrt((eqdskdata['rbound'][i]-RMAX)**2+(eqdskdata['zbound'][i]-ZMAX)**2)
        ang = npy.arctan2(eqdskdata['zbound'][i]-ZMAX,eqdskdata['rbound'][i]-RMAX)
        rho1D = npy.linspace(0.0,rho,dens)
        psi1D = npy.zeros(dens)
        for irho in range(dens):
            r = RMAX+rho1D[irho]*npy.cos(ang)
            z = ZMAX+rho1D[irho]*npy.sin(ang)
            psi1D[i] = psiRZ(r,z)
        psirho = interp1d(rho1D,psi1D,kind="cubic")
        rhopsi = bisection(psirho,rho1D[0],rho1D[-1],root=psi,eps=eps)
        if abs(psirho(rhopsi)-psi) <= eps:
           rbound.append(RMAX+rhopsi*npy.cos(ang))
           zbound.append(ZMAX+rhopsi*npy.sin(ang))

    return npy.array(rbound),npy.array(zbound)

def magsurf_interp(eqdskfpath):
    eqdskdata = read_efit_file(eqdskfpath.strip())

   #Finding the magnetic surface using interpolation technique
    R1D = eqdskdata['RR1D']-eqdskdata['RLFT']-(eqdskdata['RLEN']/2.0)
    Z1D = eqdskdata['ZR1D']-(eqdskdata['ZLEN']/2.0)
    psiRZfn  = interp2d(R1D,Z1D,eqdskdata['psiRZ'],kind="linear")

    rho1D = npy.sqrt(R1D[npy.size(R1D)/2:]**2+Z1D[npy.size(Z1D)/2:]**2)
    ang1D = npy.linspace(0.0,2.0*npy.pi,npy.size(rho1D))

    pct      = 1.00000
    psi      = (eqdskdata['PSIBND']-eqdskdata['PSIMAX'])*npy.arange(npy.size(rho1D))/(npy.size(rho1D)-1.0)
    psiN     = ((psi-psi[0])/(psi[-1]-psi[0]))*pct

    rho2D    = npy.zeros(npy.size(psiN))
    nwR2D    = npy.zeros(npy.size(psiN))
    nwZ2D    = npy.zeros(npy.size(psiN))

    for j in range(npy.size(ang1D)):
        nwpsi = []
        nwrho = []
        for i in range(npy.size(rho1D)):
            r = rho1D[i]*npy.cos(ang1D[j])
            z = rho1D[i]*npy.sin(ang1D[j])
            psival = psiRZfn(r,z)
            if psival[0] <= 1.0:
               nwpsi.append(psival[0])
               nwrho.append(rho1D[i])
            else:
               continue
        nwpsi = npy.array(nwpsi)
        nwrho = npy.array(nwrho)
        bgnloc,endloc=findmonotonic(nwpsi)
        psifn = interp1d(nwpsi[bgnloc:endloc+1],nwrho[bgnloc:endloc+1],kind="linear",fill_value="extrapolate")
        rho2D = npy.vstack((rho2D,psifn(psiN)))
        nwR2D = npy.vstack((nwR2D,rho2D[-1,:]*npy.cos(ang1D[j])))
        nwZ2D = npy.vstack((nwZ2D,rho2D[-1,:]*npy.sin(ang1D[j])))

    rho2D = npy.delete(rho2D,0,0)
    nwR2D = npy.delete(nwR2D,0,0)
    nwZ2D = npy.delete(nwZ2D,0,0)

    return nwR2D[:,-1].nwZ2D[-1,:]

