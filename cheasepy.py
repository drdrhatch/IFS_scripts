#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import time
import h5py
import traceback
import subprocess

import numpy             as npy
import matplotlib.pyplot as plt

from glob import glob
from efittools import read_efit_file
from efittools import psi2phi,phi2psi
from efittools import magsurf_solvflines
from efittools import read_profiles_file
from scipy.integrate   import trapz,simps,quad
from scipy.interpolate import splrep,splev
from scipy.interpolate import interp1d,interp2d
from scipy.interpolate import CubicSpline,RectBivariateSpline
from matplotlib.backends.backend_pdf import PdfPages

mu0 = 4.0e-7*npy.pi

def current_correction(chease,expeq={},setParam={}):
    if   'nsttp' in setParam:
         if   type(setParam['nsttp'])==list:
              current_type = setParam['nsttp'][0]
              current_src  = setParam['nsttp'][1]
         elif type(setParam['nsttp']) in [int,str]:
              current_type = setParam['nsttp']
              current_src  = 'chease'
    else:
              current_type = 'iprln'
              current_src  = 'chease'

    if type(current_type) == str: current_type.lower()
    if type(current_src)  == str: current_src.lower()

    if   expeq: expeqavail = True
    else:       expeqavail = False

    if current_src in [2,'expeq'] and not expeqavail:
       print('FATAL: Data from EXPEQ file is required but not provided')
       sys.exit()

    if 'ITEXP' in setParam: ITEXP = setParam['ITEXP']
    error = (chease['ITOR']/ITEXP)-1.0
    
    if   current_type in [3,'iprl','iprln','iparallel','iparalleln']:
         if   current_src in [0,'chease']:
              IOHMIC  = (1.0-error)*(chease['IPRLN']-chease['IBSN'])
         elif current_src in [2,'expeq']:
              IOHMIC  = (1.0-error)*(expeq['iprlN']-chease['IBSN'])
         IPRLN   = chease['IBSN'] + IOHMIC
         correction = IPRLN[:]
    elif current_type in [4,'jprl','jprln','jparallel','jparalleln']:
         if   current_src in [0,'chease']:
              JOHMIC  = (1.0-error)*(chease['JPRLN']-chease['JBSN'])
         elif current_src in [2,'expeq']:
              JOHMIC  = (1.0-error)*(expeq['jprlN']-chease['JBSN'])
         JPRLN   = chease['JBSN'] + JOHMIC
         correction = JPRLN[:]

    if 'correction' not in locals():
       print('WARNING: return correction is not calculated!')
       correction = npy.zeros(npy.size(chease['IPRLN']))

    return correction

def interp(xin,fxin,xnew,ynew=[],yout=[]):
    #splrep returns knots and coefficients for cubic spline
    #Use these knots and coefficients with splev to get a new y

    x_knots = splrep(xin,fxin)
    fxnew   = splev(xnew,x_knots,der=0)

    if len(ynew)>0 and len(yout)>0:
       y_knots = splrep(ynew,fxnew)
       fxout   = splev(yout,y_knots,der=0)

       x_knots = splrep(ynew,xnew)
       xout    = splev(yout,x_knots,der=0)
    else:
       fxout = fxnew[:]

    return fxout


def findall(inlist,item):
   #Developed by Ehab Hassan on 2019-03-28
    inds = []
    for i,character in enumerate(inlist):
        if character==item: inds.append(i)
    return inds

def namelistcreate(csvfn,rec,setParam={}):
    infid = open(csvfn, "r")
    import csv
    table = {}
    for row in csv.reader(infid):
        table[row[0]] = row[1:]
    infid.close()

    namelistkeys = table.keys()
    setParamKeys = setParam.keys()

    nrec = len(table[list(namelistkeys)[0]])
    NBSEXPQ = []
    for irec in range(nrec):
        NBSEXPQ.append(1111)

    wfh = open('chease_namelist','w')
    wfh.write('*** for EQDSK file copied onto EXPEQ file \n')
    wfh.write('*** cp this file to "chease_namelist" and run chease \n')
    wfh.write('***  \n')
    wfh.write('***  \n')
    wfh.write(' &EQDATA \n')
    if   'NTNOVA'    in namelistkeys: wfh.write(' NTNOVA=%4d,     \n' % (int(table['NTNOVA'][rec])))
    if   'RELAX'     in setParamKeys: wfh.write(' RELAX=%2.2f, '      % (float(setParam['RELAX'])))
    elif 'RELAX'     in namelistkeys: wfh.write(' RELAX=%2.2f, '      % (float(table['RELAX'][rec])))
    else:                             wfh.write(' RELAX=%2.2f, '      % (float(0.0)))
    if   'NBSEXPQ'   in setParamKeys: wfh.write(' NBSEXPQ=%04d, '     % (int(setParam['NBSEXPQ'])))
    elif 'NBSEXPQ'   in namelistkeys: wfh.write(' NBSEXPQ=%04d, '     % (int(table['NBSEXPQ'][rec])))
    else:   table['NBSEXPQ']=NBSEXPQ; wfh.write(' NBSEXPQ=1111, ')
    if   'NEQDSK'    in setParamKeys: wfh.write(' NEQDSK=%1d,     \n' % (int(setParam['NEQDSK'])))
    elif 'NEQDSK'    in namelistkeys: wfh.write(' NEQDSK=%1d,     \n' % (int(table['NEQDSK'][rec])))
    else:                             wfh.write(' NEQDSK=%1d,     \n' % (int(1)))

    if   'NS'        in setParamKeys: wfh.write(' NS=%4d, '           % (int(setParam['NS'])))
    elif 'NS'        in namelistkeys: wfh.write(' NS=%4d, '           % (int(table['NS'][rec])))
    if   'NT'        in setParamKeys: wfh.write(' NT=%4d,         \n' % (int(setParam['NT'])))
    elif 'NT'        in namelistkeys: wfh.write(' NT=%4d,         \n' % (int(table['NT'][rec])))
    if   'NPSI'      in setParamKeys: wfh.write(' NPSI=%4d, '         % (int(setParam['NPSI'])))
    elif 'NPSI'      in namelistkeys: wfh.write(' NPSI=%4d, '         % (int(table['NPSI'][rec])))
    if   'NCHI'      in setParamKeys: wfh.write(' NCHI=%4d, '         % (int(setParam['NCHI'])))
    elif 'NCHI'      in namelistkeys: wfh.write(' NCHI=%4d, '         % (int(table['NCHI'][rec])))
    if   'NISO'      in setParamKeys: wfh.write(' NISO=%4d,       \n' % (int(setParam['NISO'])))
    elif 'NISO'      in namelistkeys: wfh.write(' NISO=%4d,       \n' % (int(table['NISO'][rec])))
    if   'NRBOX'     in setParamKeys: wfh.write(' NRBOX=%4d, '        % (int(setParamKeys['NRBOX'])))
    elif 'NRBOX'     in namelistkeys: wfh.write(' NRBOX=%4d, '        % (int(table['NRBOX'][rec])))
    if   'NZBOX'     in setParamKeys: wfh.write(' NZBOX=%4d,      \n' % (int(setParamKeys['NZBOX'])))
    elif 'NZBOX'     in namelistkeys: wfh.write(' NZBOX=%4d,      \n' % (int(table['NZBOX'][rec])))

    if   'NCSCAL'    in setParamKeys: wfh.write(' NCSCAL=%1d, '       % (float(setParam['NCSCAL'])))
    elif 'NCSCAL'    in namelistkeys: wfh.write(' NCSCAL=%1d, '       % (int(table['NCSCAL'][rec])))
    else:                             wfh.write(' NCSCAL=4, ')
    if   'NOPT'      in setParamKeys: wfh.write(' NOPT=%1d,       \n' % (int(setParam['NOPT'])))
    elif 'NOPT'      in namelistkeys: wfh.write(' NOPT=%1d,       \n' % (int(table['NOPT'][rec])))
    else:                             wfh.write(' NOPT=0,         \n')
    if   'NSURF'     in setParamKeys: wfh.write(' NSURF=%1d, '        % (int(setParam['NSURF'])))
    elif 'NSURF'     in namelistkeys: wfh.write(' NSURF=%1d, '        % (int(table['NSURF'][rec])))
    else:                             wfh.write(' NSURF=6, ')
    if   'NFUNC'     in setParamKeys: wfh.write(' NFUNC=%1d,      \n' % (int(setParam['NFUNC'])))
    elif 'NFUNC'     in namelistkeys: wfh.write(' NFUNC=%1d,      \n' % (int(table['NFUNC'][rec])))
    else:                             wfh.write(' NFUNC=4,        \n')
    if   'NPPFUN'    in setParamKeys: 
                                      wfh.write(' NPPFUN=%1d,     \n' % (int(setParam['NPPFUN'])))
                                      nppfun = int(setParam['NPPFUN'])
    elif 'NPPFUN'    in namelistkeys:
                                      wfh.write(' NPPFUN=%1d,     \n' % (int(table['NPPFUN'][rec])))
                                      nppfun = int(table['NPPFUN'][rec])
    else:                             
                                      wfh.write(' NPPFUN=4,       \n')
                                      nppfun = 4
    if   'NFUNRHO'   in setParamKeys: wfh.write(' NFUNRHO=%1d, '      % (int(setParam['NFUNRHO'])))
    elif 'NFUNRHO'   in namelistkeys: wfh.write(' NFUNRHO=%1d, '      % (int(table['NFUNRHO'][rec])))
    else:                             wfh.write(' NFUNRHO=0, ')
    if   'NRHOMESH'  in setParamKeys: wfh.write(' NRHOMESH=%1d,   \n' % (int(setParam['NRHOMESH'])))
    elif 'NRHOMESH'  in namelistkeys: wfh.write(' NRHOMESH=%1d,   \n' % (int(table['NRHOMESH'][rec])))
    else:                             wfh.write(' NRHOMESH=0,     \n')
    if   'NSTTP'     in setParamKeys: wfh.write(' NSTTP=%1d, '        % (int(setParam['NSTTP'])))
    elif 'NSTTP'     in namelistkeys: wfh.write(' NSTTP=%1d, '        % (int(table['NSTTP'][rec])))
    else:                             wfh.write(' NSTTP=%1d, '        % (int(1)))
    if   'NPROPT'    in setParamKeys:
                     if   nppfun== 4: wfh.write(' NPROPT=%1d,     \n' % ( int(setParam['NPROPT'])))
                     elif nppfun== 8: wfh.write(' NPROPT=%1d,     \n' % (-int(setParam['NPROPT'])))
    elif 'NPROPT'    in namelistkeys:
                     if   nppfun== 4: wfh.write(' NPROPT=%1d,     \n' % ( int(table['NPROPT'][rec])))
                     elif nppfun== 8: wfh.write(' NPROPT=%1d,     \n' % (-int(table['NPROPT'][rec])))
    else:                             
        if   'NSTTP' in setParamKeys:
                     if   nppfun== 4: wfh.write(' NPROPT=%1d,     \n' % ( int(setParam['NSTTP'])))
                     elif nppfun== 8: wfh.write(' NPROPT=%1d,     \n' % (-int(setParam['NSTTP'])))
        elif 'NSTTP' in namelistkeys:
                     if   nppfun== 4: wfh.write(' NPROPT=%1d,     \n' % ( int(table['NSTTP'][rec])))
                     elif nppfun== 8: wfh.write(' NPROPT=%1d,     \n' % (-int(table['NSTTP'][rec])))

    if   'NVERBOSE'  in setParamKeys: wfh.write(' NVERBOSE=%1d,   \n' % (int(setParam['NVERBOSE'])))
    elif 'NVERBOSE'  in namelistkeys: wfh.write(' NVERBOSE=%1d,   \n' % (int(table['NVERBOSE'][rec])))
    else:                             wfh.write(' NVERBOSE=%1d,   \n' % (int(4)))

    if   'QSPEC'     in setParamKeys: wfh.write(' QSPEC=%3.3f, '      % (float(setParam['QSPEC'])))
    elif 'QSPEC'     in namelistkeys: wfh.write(' QSPEC=%3.3f, '      % (float(table['QSPEC'][rec])))
    if   'CSSPEC'    in setParamKeys: wfh.write(' CSSPEC=%3.3f,   \n' % (float(setParam['CSSPEC'])))
    elif 'CSSPEC'    in namelistkeys: wfh.write(' CSSPEC=%3.3f,   \n' % (float(table['CSSPEC'][rec])))
    if   'TRIANG'    in setParamKeys: wfh.write(' TRIANG=%3.3f,   \n' % (float(setParam['TRIANG'])))
    elif 'TRIANG'    in namelistkeys: wfh.write(' TRIANG=%3.3f,   \n' % (float(table['TRIANG'][rec])))

    if   'R0'        in setParamKeys: wfh.write(' R0=%10.8f, '        % (float(setParam['R0'])))
    elif 'R0'        in namelistkeys: wfh.write(' R0=%10.8f, '        % (float(table['R0'][rec])))
    if   'RZ0'       in setParamKeys: wfh.write(' RZ0=%10.8f,     \n' % (float(setParam['RZ0'])))
    elif 'RZ0'       in namelistkeys: wfh.write(' RZ0=%10.8f,     \n' % (float(table['RZ0'][rec])))

    if   'RBOXLEN'   in setParamKeys: wfh.write(' RBOXLEN=%3.3f, '    % (float(setParam['RBOXLEN'])))
    elif 'RBOXLEN'   in namelistkeys: wfh.write(' RBOXLEN=%3.3f, '    % (float(table['RBOXLEN'][rec])))
    if   'ZBOXLEN'   in setParamKeys: wfh.write(' ZBOXLEN=%3.3f, '    % (float(setParam['ZBOXLEN'])))
    elif 'ZBOXLEN'   in namelistkeys: wfh.write(' ZBOXLEN=%3.3f, '    % (float(table['ZBOXLEN'][rec])))
    if   'RBOXLFT'   in setParamKeys: wfh.write(' RBOXLFT=%3.3f,  \n' % (float(setParam['RBOXLFT'])))
    elif 'RBOXLFT'   in namelistkeys: wfh.write(' RBOXLFT=%3.3f,  \n' % (float(table['RBOXLFT'][rec])))

    if   'R0EXP'     in setParamKeys: wfh.write(' R0EXP=%3.3f, '      % (float(setParam['R0EXP'])))
    elif 'R0EXP'     in namelistkeys: wfh.write(' R0EXP=%3.3f, '      % (float(table['R0EXP'][rec])))
    if   'B0EXP'     in setParamKeys: wfh.write(' B0EXP=%3.3f,    \n' % (float(setParam['B0EXP'])))
    elif 'B0EXP'     in namelistkeys: wfh.write(' B0EXP=%3.3f,    \n' % (float(table['B0EXP'][rec])))

    if   'NDIAGOP'   in setParamKeys: wfh.write(' NDIAGOP=%1d, '      % (int(setParam['NDIAGOP'])))
    elif 'NDIAGOP'   in namelistkeys: wfh.write(' NDIAGOP=%1d, '      % (int(table['NDIAGOP'][rec])))
    else:                             wfh.write(' NDIAGOP=%1d, '      % (int(1)))
    if   'NIDEAL'    in setParamKeys: wfh.write(' NIDEAL=%1d,     \n' % (int(setParam['NIDEAL'])))
    elif 'NIDEAL'    in namelistkeys: wfh.write(' NIDEAL=%1d,     \n' % (int(table['NIDEAL'][rec])))
    else:                             wfh.write(' NIDEAL=9,       \n')
    if   'NDIFPS'    in setParamKeys: wfh.write(' NDIFPS=%1d, '       % (int(setParam['NDIFPS'])))
    elif 'NDIFPS'    in namelistkeys: wfh.write(' NDIFPS=%1d, '       % (int(table['NDIFPS'][rec])))
    else:                             wfh.write(' NDIFPS=%1d,     \n' % (int(0)))
    if   'NDIFT'     in setParamKeys: wfh.write(' NDIFT=%1d,      \n' % (int(setParam['NDIFT'])))
    elif 'NDIFT'     in namelistkeys: wfh.write(' NDIFT=%1d,      \n' % (int(table['NDIFT'][rec])))
    else:                             wfh.write(' NDIFT=%1d,      \n' % (int(1)))

    if   'NMESHC'    in namelistkeys: wfh.write(' NMESHC=%1d, '       % (int(table['NMESHC'][rec])))
    else:                             wfh.write(' NMESHC=1,')
    if   'NPOIDC'    in namelistkeys: wfh.write(' NPOIDC=%1d, '       % (int(table['NPOIDC'][rec])))
    else:                             wfh.write(' NPOIDC=2,')
    if   'SOLPDC'    in namelistkeys: wfh.write(' SOLPDC=%2.2f,   \n' % (float(table['SOLPDC'][rec])))
    else:                             wfh.write(' SOLPDC=0.7, \n')
    if   'CPLACE'    in namelistkeys: 
       CPLACE = [float(f) for f in table['CPLACE'][rec].split(',')]
       wfh.write(' CPLACE=')
       for cplace in CPLACE:
                                      wfh.write('%4.4f,'              % float(cplace))
       wfh.write('\n')
    else:                             wfh.write(' CPLACE=0.9,0.95,0.99,1.0, \n')
    if   'CWIDTH'    in namelistkeys: 
       CWIDTH = [float(f) for f in table['CWIDTH'][rec].split(',')]
       wfh.write(' CWIDTH=')
       for cwidth    in CWIDTH:
                                      wfh.write('%4.4f,'              % float(cwidth))
       wfh.write('\n')
    else:                             wfh.write(' CWIDTH=0.1,0.05,0.01,0.025, \n')

    if   'NMESHD'    in namelistkeys: wfh.write(' NMESHD=%1d, '       % (int(table['NMESHD'][rec])))
    if   'NPOIDD'    in namelistkeys: wfh.write(' NPOIDD=%1d, '       % (int(table['NPOIDD'][rec])))
    if   'SOLPDD'    in namelistkeys: wfh.write(' SOLPDD=%2.2f,   \n' % (float(table['SOLPDD'][rec])))
    if   'DPLACE'    in namelistkeys: 
       DPLACE = [float(f) for f in table['DPLACE'][rec].split(',')]
       wfh.write(' DPLACE=')
       for dplace    in DPLACE:
                                    wfh.write('%4.4f,'              % float(dplace))
       wfh.write('\n')
    if 'DWIDTH'      in namelistkeys: 
       DWIDTH = [float(f) for f in table['DWIDTH'][rec].split(',')]
       wfh.write(' DWIDTH=')
       for dwidth in DWIDTH:
                                    wfh.write('%4.4f,'              % float(dwidth))
       wfh.write('\n')
    if   'NMESHPOL'  in setParamKeys: wfh.write(' NMESHPOL=%4d, '     % (int(setParam['NMESHPOL'])))
    elif 'NMESHPOL'  in namelistkeys: wfh.write(' NMESHPOL=%4d, '     % (int(table['NMESHPOL'][rec])))
    else:                             wfh.write(' NMESHPOL=%4d, '     % (int(1)))
    if   'SOLPDPOL'  in setParamKeys: wfh.write(' SOLPDPOL=%2.2f, \n' % (float(setParam['SOLPDPOL'])))
    elif 'SOLPDPOL'  in namelistkeys: wfh.write(' SOLPDPOL=%2.2f, \n' % (float(table['SOLPDPOL'][rec])))
    else:                             wfh.write(' SOLPDPOL=%2.2f, \n' % (float(0.1)))
    if   'NTURN'     in setParamKeys: wfh.write(' NTURN=%2d, '        % (int(setParam['NTURN'])))
    elif 'NTURN'     in namelistkeys: wfh.write(' NTURN=%2d, '        % (int(table['NTURN'][rec])))
    else:                             wfh.write(' NTURN=%2d, '        % (int(20)))
    if   'NBLC0'     in setParamKeys: wfh.write(' NBLC0=%2d, '        % (int(setParam['NBLC0'])))
    elif 'NBLC0'     in namelistkeys: wfh.write(' NBLC0=%2d, '        % (int(table['NBLC0'][rec])))
    else:                             wfh.write(' NBLC0=%2d, '        % (int(16)))
    if   'NPPR'      in setParamKeys: wfh.write(' NPPR=%2d,       \n' % (int(setParam['NPPR'])))
    elif 'NPPR'      in namelistkeys: wfh.write(' NPPR=%2d,       \n' % (int(table['NPPR'][rec])))
    else:                             wfh.write(' NPPR=%2d,       \n' % (int(24)))
    if   'NINMAP'    in setParamKeys: wfh.write(' NINMAP=%2d, '       % (int(setParam['NINMAP'])))
    elif 'NINMAP'    in namelistkeys: wfh.write(' NINMAP=%2d, '       % (int(table['NINMAP'][rec])))
    else:                             wfh.write(' NINMAP=%2d, '       % (int(40)))
    if   'NINSCA'    in setParamKeys: wfh.write(' NINSCA=%2d,     \n' % (int(setParam['NINSCA'])))
    elif 'NINSCA'    in namelistkeys: wfh.write(' NINSCA=%2d,     \n' % (int(table['NINSCA'][rec])))
    else:                             wfh.write(' NINSCA=%2d,     \n' % (int(40)))
    if   'NSYM'      in setParamKeys: wfh.write(' NSYM=%1d, '         % (int(setParam['NSYM'])))
    elif 'NSYM'      in namelistkeys: wfh.write(' NSYM=%1d, '         % (int(table['NSYM'][rec])))
    else:                             wfh.write(' NSYM=%1d, '         % (int(0)))
    if   'NEGP'      in setParamKeys: wfh.write(' NEGP=%1d, '         % (int(setParam['NEGP'])))
    elif 'NEGP'      in namelistkeys: wfh.write(' NEGP=%1d, '         % (int(table['NEGP'][rec])))
    else:                             wfh.write(' NEGP=%1d, '         % (int(0)))
    if   'NER'       in setParamKeys: wfh.write(' NER=%1d,        \n' % (int(setParam['NER'])))
    elif 'NER'       in namelistkeys: wfh.write(' NER=%1d,        \n' % (int(table['NER'][rec])))
    else:                             wfh.write(' NER=%1d,        \n' % (int(2)))
    if   'EPSLON'    in setParamKeys: wfh.write(' EPSLON=%6.2E,   \n' % (float(setParam['EPSLON'])))
    elif 'EPSLON'    in namelistkeys: wfh.write(' EPSLON=%6.2E,   \n' % (float(table['EPSLON'][rec])))
    else:                             wfh.write(' EPSLON=%6.2E,   \n' % (float(1.0E-10)))
    if   'ETAEI'     in setParamKeys: wfh.write(' ETAEI=%2.1f, '      % (float(setParam['ETAEI'])))
    elif 'ETAEI'     in namelistkeys: wfh.write(' ETAEI=%2.1f, '      % (float(table['ETAEI'][rec])))
    else:                             wfh.write(' ETAEI=%2.1f, '      % (float(3.0)))
    if   'RPEOP'     in setParamKeys: wfh.write(' RPEOP=%2.1f, '      % (float(setParam['RPEOP'])))
    elif 'RPEOP'     in namelistkeys: wfh.write(' RPEOP=%2.1f, '      % (float(table['RPEOP'][rec])))
    else:                             wfh.write(' RPEOP=%2.1f, '      % (float(0.5)))
    if   'RZION'     in setParamKeys: wfh.write(' RZION=%2.1f, '      % (float(setParam['RZION'])))
    elif 'RZION'     in namelistkeys: wfh.write(' RZION=%2.1f, '      % (float(table['RZION'][rec])))
    else:                             wfh.write(' RZION=%2.1f, '      % (float(1.0)))
    if   'GAMMA'     in setParamKeys: wfh.write(' GAMMA=%12.11f,  \n' % (float(setParam['GAMMA'])))
    elif 'GAMMA'     in namelistkeys: wfh.write(' GAMMA=%12.11f,  \n' % (float(table['GAMMA'][rec])))
    else:                             wfh.write(' GAMMA=%12.11f,  \n' % (float(1.6666666667)))
    if   'AT3(1)'    in setParamKeys: wfh.write(' AT3(1)=%2.2f,   \n' % (float(setParam['AT3(1)'])))
    elif 'AT3(1)'    in namelistkeys: wfh.write(' AT3(1)=%2.2f,   \n' % (float(table['AT3(1)'][rec])))
    else:                             wfh.write(' AT3(1)=%2.2f,   \n' % (float(-0.69)))
    if   'TENSPROF'  in setParamKeys: wfh.write(' TENSPROF=%2.2f, \n' % (float(setParam['TENSPROF'])))
    elif 'TENSPROF'  in namelistkeys: wfh.write(' TENSPROF=%2.2f, \n' % (float(table['TENSPROF'][rec])))
    if   'TENSBND'   in setParamKeys: wfh.write(' TENSBND=%2.2f,  \n' % (float(setParam['TENSBND'])))
    elif 'TENSBND'   in namelistkeys: wfh.write(' TENSBND=%2.2f,  \n' % (float(table['TENSBND'][rec])))
    if   'cocos_in'  in setParamKeys: wfh.write(' cocos_in=%1d,   \n' % (int(setParam['cocos_in'])))
    elif 'cocos_in'  in namelistkeys: wfh.write(' cocos_in=%1d,   \n' % (int(table['cocos_in'][rec])))
    if   'cocos_out' in setParamKeys: wfh.write(' cocos_out=%2d   \n' % (int(setParam['cocos_out'])))
    elif 'cocos_out' in namelistkeys: wfh.write(' cocos_out=%2d   \n' % (int(table['cocos_out'][rec])))
    wfh.write(' &END \n')
    wfh.write('\n')
    wfh.close()

    return table

def derivative(x,fx,axis=0,dorder=1,method='gradient'):
    fxShape = npy.shape(fx)
    nDim    = len(fxShape)
    if   method=='gradient':
         if   nDim == 1:
              dfdx = npy.gradient(fx,x)
    elif method=='CubicSpline':
         if   nDim == 1:
              CS = CubicSpline(x,fx)
              dfdx = CS(x,dorder)
         elif nDim == 2:
              m = fxShape[0]
              n = fxShape[1]
              dfdx=npy.zeros((m,n))
              if   axis == 0:
                   for j in range(n):
                      CS = CubicSpline(x,fx[:,j])
                      dfdx[:,j] = CS(x,dorder)
              elif axis == 1:
                   for i in range(m):
                      CS = CubicSpline(x,fx[i,:])
                      dfdx[i,:] = CS(x,dorder)
    return dfdx

def integrate(x,fx,axis=0,method='trapz'):
    fxShape = npy.shape(fx)
    nDim    = len(fxShape)
    if   method in ['trapz','trapzoid']:
         if   nDim == 1:
              intf = trapz(y=fx,x=x,dx=npy.argmin(diff(x)))
    elif method in ['simps','simpson']:
         if   nDim == 1:
              intf = simps(y=fx,x=x,dx=npy.argmin(diff(x)))
    elif method=='CubicSpline':
         m = fxShape[0]
         n = fxShape[1]
         if   nDim == 1:
                   intf = CubicSpline(x,fx).integrate(x[0],x[-1])
         elif nDim == 2:
              if   axis == 0:
                   try:
                      CS   = CubicSpline(x,fx,axis=0,bc_type='periodic',extrapolate='periodic')
                      intf = CS.integrate(x[0],x[-1],extrapolate='periodic')
                   except ValueError:
                      intf = npy.zeros(m)
                      for j in range(m):
                          func  = interp1d(x,fx[:,j],kind='linear')
                          fyfnc = lambda z: func(z)
                          intf[j] = quad(fyfnc,x[0],x[-1])
              elif axis == 1:
                   try:
                      CS   = CubicSpline(x,fx,axis=1)
                      intf = CS.integrate(x[0],x[-1])
                   except ValueError:
                      intf = npy.zeros(n)
                      for j in range(n):
                          func  = interp1d(x,fx[:,j],kind='linear')
                          fyfnc = lambda z: func(z)
                          intf[j] = quad(fyfnc,x[0],x[-1])
    return intf

def read_chease(cheasefpath,setParam={},**kwargs):
    if os.path.isfile(cheasefpath) == False:
       errorFunc = traceback.extract_stack(limit=2)[-2][3]
       errorLine = traceback.extract_stack(limit=2)[-2][1]
       errorFile = traceback.extract_stack(limit=2)[-2][2]
       errMSG    = 'Call %s line %5d in file %s Failed.\n'
       errMSG   += 'Fatal: file %s not found.'
       raise IOError(errMSG %(errorFunc,errorLine,errorFile,cheasefpath))

    extendCHI  = False
    reverseCHI = False

    hdffh = h5py.File(cheasefpath,'r')
    hdffhKeys = list(hdffh.keys())[0]
    datag = hdffh[hdffhKeys]

    CHEASEdata                = {}

    CHEASEdata['NPSI']        = datag.attrs["NPSI"]
    CHEASEdata['NCHI']        = datag.attrs["NCHI"]
    CHEASEdata['NRBOX']       = datag.attrs["NRBOX"]
    CHEASEdata['NZBOX']       = datag.attrs["NZBOX"]
    CHEASEdata['NBOUND']      = datag.attrs["NBOUND"]
    CHEASEdata['B0EXP']       = datag.attrs["B0EXP"]
    CHEASEdata['R0EXP']       = datag.attrs["R0EXP"]

    CHEASEdata['PSI']         = npy.array(datag["grid"]["PSI"])
    CHEASEdata['rhopsi']      = npy.sqrt(CHEASEdata['PSI']/CHEASEdata['PSI'][-1])
    CHEASEdata['rhotor']      = npy.array(datag["var1d"]["rho_tor"])

    CHEASEdata['CHI']         = npy.array(datag["grid"]["CHI"])
    if CHEASEdata['CHI'][0]>CHEASEdata['CHI'][1]:
       CHEASEdata['CHI']      = CHEASEdata["CHI"][::-1]
       reverseCHI             = True
    if min(abs(CHEASEdata['CHI']-2.0*npy.pi))>=min(npy.diff(CHEASEdata['CHI'])):
       CHEASEdata['CHI']      = npy.append(CHEASEdata['CHI'],2.0*npy.pi)
       extendCHI              = True

    CHEASEdata['JBS']         = npy.array(datag["var1d"]["jbsBav"])
    CHEASEdata['Zeff']        = npy.array(datag["var1d"]["zeff"])
    CHEASEdata['kappa']       = npy.array(datag["var1d"]["kappa"])
    CHEASEdata['shear']       = npy.array(datag["var1d"]["shear"])
    CHEASEdata['signeo']      = npy.array(datag["var1d"]["signeo"])
    CHEASEdata['PPrime']      = 2.0*npy.pi*npy.array(datag["var1d"]["dpdpsi"])

    CHEASEdata['P']           = npy.array(datag["var1d"]["p"])
    CHEASEdata['q']           = npy.array(datag["var1d"]["q"])
    CHEASEdata['R_av']        = npy.array(datag["var1d"]["R_av"])
    CHEASEdata['ageom']       = npy.array(datag["var1d"]["ageom"])
    CHEASEdata['Rgeom']       = npy.array(datag["var1d"]["Rgeom"])
    CHEASEdata['Volume']      = npy.array(datag["var1d"]["Volume"])
    CHEASEdata['GDPSI_av']    = npy.array(datag["var1d"]["GDPSI_av"])
    CHEASEdata['radius_av']   = npy.array(datag["var1d"]["radius_av"])

    CHEASEdata['Te']          = npy.array(datag["var1d"]["Te"])
    CHEASEdata['Ti']          = npy.array(datag["var1d"]["Ti"])
    CHEASEdata['ne']          = npy.array(datag["var1d"]["ne"])
    CHEASEdata['ni']          = npy.array(datag["var1d"]["ni"])

    CHEASEdata['Zi']          = 1.0
    CHEASEdata['Zz']          = 6.0
    CHEASEdata['nz']          = CHEASEdata['Zeff']*CHEASEdata['ne']
    CHEASEdata['nz']         -= CHEASEdata['ni']*CHEASEdata['Zi']**2
    CHEASEdata['nz']         /= CHEASEdata['Zz']**2

    CHEASEdata['TePrime']     = npy.array(datag["var1d"]["dTedpsi"])
    CHEASEdata['TiPrime']     = npy.array(datag["var1d"]["dTidpsi"])
    CHEASEdata['nePrime']     = npy.array(datag["var1d"]["dnedpsi"])
    CHEASEdata['niPrime']     = npy.array(datag["var1d"]["dnidpsi"])

    CHEASEdata['T']           = npy.array(datag["var1d"]["f"])
    CHEASEdata['TTPrime']     = 2.0*npy.pi*npy.array(datag["var1d"]["fdfdpsi"])
    CHEASEdata['TPrime']      = CHEASEdata['TTPrime']/CHEASEdata['T']

    CHEASEdata['rmesh']       = npy.array(datag["var1d"]["rmesh"])
    CHEASEdata['zmesh']       = npy.array(datag["var1d"]["zmesh"])
    CHEASEdata['rbound']      = npy.array(datag["var1d"]["rboundplasma"])
    CHEASEdata['zbound']      = npy.array(datag["var1d"]["zboundplasma"])
    CHEASEdata['delta_upper'] = npy.array(datag["var1d"]["delta_upper"])
    CHEASEdata['delta_lower'] = npy.array(datag["var1d"]["delta_lower"])

    CHEASEdata['dpdpsi']      = npy.array(datag["var1d"]["dpdpsi"])
    CHEASEdata['dqdpsi']      = npy.array(datag["var1d"]["dqdpsi"])
    CHEASEdata['dVdpsi']      = npy.array(datag["var1d"]["dVdpsi"])
    CHEASEdata['d2qdpsi2']    = npy.array(datag["var1d"]["d2qdpsi2"])
    CHEASEdata['dsheardpsi']  = npy.array(datag["var1d"]["dsheardpsi"])
    CHEASEdata['dpsidrhotor'] = npy.array(datag["var1d"]["dpsidrhotor"])


   #THE DIMENSION OF ALL THE FOLLOWING QUNATITIES ARE (NCHI,NPSI)
    CHEASEdata['R']           = npy.array(datag["var2d"]["R"])
    CHEASEdata['Z']           = npy.array(datag["var2d"]["Z"])
    CHEASEdata['B']           = npy.array(datag["var2d"]["B"])
    CHEASEdata['J']           = npy.array(datag["var2d"]["Jacobian"])
    CHEASEdata['g11']         = npy.array(datag["var2d"]["g11"])
    CHEASEdata['g22']         = npy.array(datag["var2d"]["g22"])
    CHEASEdata['g33']         = npy.array(datag["var2d"]["g33"])
    CHEASEdata['dBdpsi']      = npy.array(datag["var2d"]["dBdpsi"])
    CHEASEdata['dBdchi']      = npy.array(datag["var2d"]["dBdchi"])
    CHEASEdata['dChidZ']      = npy.array(datag["var2d"]["dChidZ"])
    CHEASEdata['dPsidZ']      = npy.array(datag["var2d"]["dPsidZ"])
    CHEASEdata['dChidR']      = npy.array(datag["var2d"]["dChidR"])
    CHEASEdata['dPsidR']      = npy.array(datag["var2d"]["dPsidR"])

    if extendCHI:
       CHEASEdata['R']        = npy.vstack((CHEASEdata['R'],CHEASEdata['R'][0,:]))
       CHEASEdata['Z']        = npy.vstack((CHEASEdata['Z'],CHEASEdata['Z'][0,:]))
       CHEASEdata['B']        = npy.vstack((CHEASEdata['B'],CHEASEdata['B'][0,:]))
       CHEASEdata['J']        = npy.vstack((CHEASEdata['J'],CHEASEdata['J'][0,:]))
       CHEASEdata['g11']      = npy.vstack((CHEASEdata['g11'],CHEASEdata['g11'][0,:]))
       CHEASEdata['g22']      = npy.vstack((CHEASEdata['g22'],CHEASEdata['g22'][0,:]))
       CHEASEdata['g33']      = npy.vstack((CHEASEdata['g33'],CHEASEdata['g33'][0,:]))
       CHEASEdata['dBdpsi']   = npy.vstack((CHEASEdata['dBdpsi'],CHEASEdata['dBdpsi'][0,:]))
       CHEASEdata['dBdchi']   = npy.vstack((CHEASEdata['dBdchi'],CHEASEdata['dBdchi'][0,:]))
       CHEASEdata['dChidZ']   = npy.vstack((CHEASEdata['dChidZ'],CHEASEdata['dChidZ'][0,:]))
       CHEASEdata['dPsidZ']   = npy.vstack((CHEASEdata['dPsidZ'],CHEASEdata['dPsidZ'][0,:]))
       CHEASEdata['dChidR']   = npy.vstack((CHEASEdata['dChidR'],CHEASEdata['dChidR'][0,:]))
       CHEASEdata['dPsidR']   = npy.vstack((CHEASEdata['dPsidR'],CHEASEdata['dPsidR'][0,:]))

   #THE DIMENSION OF ALL THE FOLLOWING QUNATITIES ARE (NRBOX,NZBOX)
    CHEASEdata['psiRZ']       = npy.array(datag["var2d"]["psiRZ"])
    CHEASEdata['chiRZ']       = npy.array(datag["var2d"]["chiRZ"])

    CHEASEdata['C0']          = npy.trapz(y=CHEASEdata['J']/CHEASEdata['R'],                    x=CHEASEdata['CHI'],axis=0)
    CHEASEdata['C1']          = npy.trapz(y=CHEASEdata['J'],                                    x=CHEASEdata['CHI'],axis=0)
    CHEASEdata['C2']          = npy.trapz(y=CHEASEdata['J']/CHEASEdata['R']**2,                 x=CHEASEdata['CHI'],axis=0)
    CHEASEdata['C3']          = npy.trapz(y=CHEASEdata['J']*CHEASEdata['g11']*CHEASEdata['g33'],x=CHEASEdata['CHI'],axis=0)

    CHEASEdata['y1']          = 1.0+CHEASEdata['C3']/CHEASEdata['C2']/CHEASEdata['T']**2/4.0/npy.pi**2

    CHEASEdata['<B2>']        = npy.trapz(y=CHEASEdata['J']*CHEASEdata['B']**2,x=CHEASEdata['CHI'],axis=0)/CHEASEdata['C1']
    CHEASEdata['<JdotB>']     =-CHEASEdata['T']*CHEASEdata['PPrime']-CHEASEdata['TPrime']*CHEASEdata['<B2>']/mu0
    CHEASEdata['JPRL']        = CHEASEdata['<JdotB>']/CHEASEdata['B0EXP']

    CHEASEdata['<T/R2>']      = npy.trapz(y=CHEASEdata['J']*CHEASEdata['T']*CHEASEdata['g33'],x=CHEASEdata['CHI'],axis=0)/CHEASEdata['C1']
    CHEASEdata['IPRL']        = CHEASEdata['R0EXP']*CHEASEdata['<JdotB>']/CHEASEdata['<T/R2>']

    CHEASEdata['ISTR']        =-((CHEASEdata['C2']/CHEASEdata['C0'])*(CHEASEdata['TTPrime']/mu0))
    CHEASEdata['ISTR']       +=-((CHEASEdata['C1']/CHEASEdata['C0'])*CHEASEdata['PPrime'])
    CHEASEdata['ISTR']       *= CHEASEdata['R0EXP']**2

    CHEASEdata['JPHI']        =-(CHEASEdata['R']*CHEASEdata['PPrime'])-(CHEASEdata['TTPrime']/(mu0*CHEASEdata['R']))

    CHEASEdata['JTOR']        = npy.trapz(y=CHEASEdata['JPHI']*CHEASEdata['J']/CHEASEdata['R'],x=CHEASEdata['CHI'],axis=0)
    CHEASEdata['JTORN']       = CHEASEdata['JTOR']*mu0*CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
    CHEASEdata['ITOR']        = npy.trapz(y=CHEASEdata['JTOR'],x=CHEASEdata['PSI'],axis=0)


    CHEASEdata['IBS']         = CHEASEdata['R0EXP']*CHEASEdata['JBS']/CHEASEdata['<T/R2>']
    CHEASEdata['IOHMIC']      = CHEASEdata['IPRL']-CHEASEdata['IBS']
    CHEASEdata['JOHMIC']      = CHEASEdata['JPRL']-CHEASEdata['JBS']

    extendPSI    = npy.linspace(CHEASEdata['PSI'][0],CHEASEdata['PSI'][-1],10*npy.size(CHEASEdata['PSI']))
    extendPHI    = npy.empty_like(extendPSI)
    extendPHI[0] = 0.0
    qfunc        = CubicSpline(CHEASEdata['PSI'],CHEASEdata['q'])
    for i in range(1,npy.size(extendPSI)):
        x           = extendPSI[:i+1]
        y           = qfunc(x)
        extendPHI[i]= npy.trapz(y,x)

    CHEASEdata['PHI'] = npy.empty_like(CHEASEdata['PSI'])
    phifunc           = CubicSpline(extendPSI,extendPHI)
    for i in range(npy.size(CHEASEdata['PSI'])):
        CHEASEdata['PHI'][i] = phifunc(CHEASEdata['PSI'][i])

    CHEASEdata['PHIN']     = (CHEASEdata['PHI']-CHEASEdata['PHI'][0])/(CHEASEdata['PHI'][-1]-CHEASEdata['PHI'][0])
    CHEASEdata['PSIN']     = (CHEASEdata['PSI']-CHEASEdata['PSI'][0])/(CHEASEdata['PSI'][-1]-CHEASEdata['PSI'][0])
    CHEASEdata['CHIN']     = (CHEASEdata['CHI']-CHEASEdata['CHI'][0])/(CHEASEdata['CHI'][-1]-CHEASEdata['CHI'][0])

    CHEASEdata['rmeshN']   = CHEASEdata['rmesh']/CHEASEdata['R0EXP'] 
    CHEASEdata['zmeshN']   = CHEASEdata['zmesh']/CHEASEdata['R0EXP'] 
    CHEASEdata['rhopsiN']  = (CHEASEdata['rhopsi']-CHEASEdata['rhopsi'][0])/(CHEASEdata['rhopsi'][-1]-CHEASEdata['rhopsi'][0])
    CHEASEdata['rhotorN']  = (CHEASEdata['rhotor']-CHEASEdata['rhotor'][0])/(CHEASEdata['rhotor'][-1]-CHEASEdata['rhotor'][0])
    CHEASEdata['rhotorN']  = npy.sqrt(CHEASEdata['PHIN'])
    CHEASEdata['rboundN']  = CHEASEdata['rbound']/CHEASEdata['R0EXP'] 
    CHEASEdata['zboundN']  = CHEASEdata['zbound']/CHEASEdata['R0EXP']

    #Implementing Interpolation to EQDSK Grid
    rhopsiflag = False; rhotorflag = False
    if 'nrhomesh' in setParam:
       if   setParam['nrhomesh'] in [0,'rhopsi']: rhopsiflag = True
       elif setParam['nrhomesh'] in [1,'rhotor']: rhotorflag = True

    eqdskflag    = False
    expeqflag    = False
    interpflag   = False
    for key,value in kwargs.items():
        if   key in ['eqdsk','eqdskdata','eqdskfpath']:
             if   type(value)==str and os.path.isfile(value.strip()):
                                       eqdskdata = read_eqdsk(eqdskfpath=value.strip())
             if    'rhopsi' in eqdskdata and rhopsiflag:
                   interprho = eqdskdata['PSIN'][:]; interpflag = True
             elif  'rhotor' in eqdskdata and rhotorflag:
                   interprho = eqdskdata['PHIN'][:]; interpflag = True
             psi    = eqdskdata['PSI'];   phi    = eqdskdata['PHI']
             rhopsi = eqdskdata['rhopsi'];rhotor = eqdskdata['rhotor']
             eqdskflag = True

    if interpflag:
       if   rhopsiflag:
            CHEASEdata['q']       = interp(CHEASEdata['PSI'],CHEASEdata['q'],psi)
            CHEASEdata['T']       = interp(CHEASEdata['PSI'],CHEASEdata['T'],psi)
            CHEASEdata['P']       = interp(CHEASEdata['PSI'],CHEASEdata['P'],psi)
            CHEASEdata['Te']      = interp(CHEASEdata['PSI'],CHEASEdata['Te'],psi)
            CHEASEdata['Ti']      = interp(CHEASEdata['PSI'],CHEASEdata['Ti'],psi)
            CHEASEdata['ne']      = interp(CHEASEdata['PSI'],CHEASEdata['ne'],psi)
            CHEASEdata['ni']      = interp(CHEASEdata['PSI'],CHEASEdata['ni'],psi)
            CHEASEdata['nz']      = interp(CHEASEdata['PSI'],CHEASEdata['nz'],psi)
            CHEASEdata['IBS']     = interp(CHEASEdata['PSI'],CHEASEdata['IBS'],psi)
            CHEASEdata['JBS']     = interp(CHEASEdata['PSI'],CHEASEdata['JBS'],psi)
            CHEASEdata['Zeff']    = interp(CHEASEdata['PSI'],CHEASEdata['Zeff'],psi)
            CHEASEdata['ISTR']    = interp(CHEASEdata['PSI'],CHEASEdata['ISTR'],psi)
            CHEASEdata['IPRL']    = interp(CHEASEdata['PSI'],CHEASEdata['IPRL'],psi)
            CHEASEdata['JPRL']    = interp(CHEASEdata['PSI'],CHEASEdata['JPRL'],psi)
            CHEASEdata['kappa']   = interp(CHEASEdata['PSI'],CHEASEdata['kappa'],psi)
            CHEASEdata['shear']   = interp(CHEASEdata['PSI'],CHEASEdata['shear'],psi)
            CHEASEdata['signeo']  = interp(CHEASEdata['PSI'],CHEASEdata['signeo'],psi)
            CHEASEdata['TPrime']  = interp(CHEASEdata['PSI'],CHEASEdata['TPrime'],psi)
            CHEASEdata['IOHMIC']  = interp(CHEASEdata['PSI'],CHEASEdata['IOHMIC'],psi)
            CHEASEdata['JOHMIC']  = interp(CHEASEdata['PSI'],CHEASEdata['JOHMIC'],psi)
            CHEASEdata['PPrime']  = interp(CHEASEdata['PSI'],CHEASEdata['PPrime'],psi)
            CHEASEdata['TTPrime'] = interp(CHEASEdata['PSI'],CHEASEdata['TTPrime'],psi)
            CHEASEdata['rhopsiN'] = rhopsi[:]
            CHEASEdata['rhotorN'] = rhotor[:]
       elif rhotorflag:
            CHEASEdata['q']       = interp(CHEASEdata['PSI'],CHEASEdata['q'],psi,phi,phi)
            CHEASEdata['T']       = interp(CHEASEdata['PSI'],CHEASEdata['T'],psi,phi,phi)
            CHEASEdata['P']       = interp(CHEASEdata['PSI'],CHEASEdata['P'],psi,phi,phi)
            CHEASEdata['Te']      = interp(CHEASEdata['PSI'],CHEASEdata['Te'],psi,phi,phi)
            CHEASEdata['Ti']      = interp(CHEASEdata['PSI'],CHEASEdata['Ti'],psi,phi,phi)
            CHEASEdata['ne']      = interp(CHEASEdata['PSI'],CHEASEdata['ne'],psi,phi,phi)
            CHEASEdata['ni']      = interp(CHEASEdata['PSI'],CHEASEdata['ni'],psi,phi,phi)
            CHEASEdata['nz']      = interp(CHEASEdata['PSI'],CHEASEdata['nz'],psi,phi,phi)
            CHEASEdata['IBS']     = interp(CHEASEdata['PSI'],CHEASEdata['IBS'],psi,phi,phi)
            CHEASEdata['JBS']     = interp(CHEASEdata['PSI'],CHEASEdata['JBS'],psi,phi,phi)
            CHEASEdata['Zeff']    = interp(CHEASEdata['PSI'],CHEASEdata['Zeff'],psi,phi,phi)
            CHEASEdata['ISTR']    = interp(CHEASEdata['PSI'],CHEASEdata['ISTR'],psi,phi,phi)
            CHEASEdata['IPRL']    = interp(CHEASEdata['PSI'],CHEASEdata['IPRL'],psi,phi,phi)
            CHEASEdata['JPRL']    = interp(CHEASEdata['PSI'],CHEASEdata['JPRL'],psi,phi,phi)
            CHEASEdata['kappa']   = interp(CHEASEdata['PSI'],CHEASEdata['kappa'],psi,phi,phi)
            CHEASEdata['shear']   = interp(CHEASEdata['PSI'],CHEASEdata['shear'],psi,phi,phi)
            CHEASEdata['signeo']  = interp(CHEASEdata['PSI'],CHEASEdata['signeo'],psi,phi,phi)
            CHEASEdata['TPrime']  = interp(CHEASEdata['PSI'],CHEASEdata['TPrime'],psi,phi,phi)
            CHEASEdata['IOHMIC']  = interp(CHEASEdata['PSI'],CHEASEdata['IOHMIC'],psi,phi,phi)
            CHEASEdata['JOHMIC']  = interp(CHEASEdata['PSI'],CHEASEdata['JOHMIC'],psi,phi,phi)
            CHEASEdata['PPrime']  = interp(CHEASEdata['PSI'],CHEASEdata['PPrime'],psi,phi,phi)
            CHEASEdata['TTPrime'] = interp(CHEASEdata['PSI'],CHEASEdata['TTPrime'],psi,phi,phi)
            CHEASEdata['rhotorN'] = rhopsi[:]
            CHEASEdata['rhopsiN'] = rhotor[:]

    CHEASEdata['RN']       = CHEASEdata['R']/CHEASEdata['R0EXP']
    CHEASEdata['ZN']       = CHEASEdata['Z']/CHEASEdata['R0EXP']
    CHEASEdata['BN']       = CHEASEdata['B']/CHEASEdata['B0EXP']
    CHEASEdata['PN']       = CHEASEdata['P']*mu0/CHEASEdata['B0EXP']**2
    CHEASEdata['IBSN']     = CHEASEdata['IBS']*mu0/CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
    CHEASEdata['JBSN']     = CHEASEdata['JBS'] *mu0*CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
    CHEASEdata['IPRLN']    = CHEASEdata['IPRL']*mu0/CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
    CHEASEdata['JPRLN']    = CHEASEdata['JPRL']*mu0*CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
    CHEASEdata['ISTRN']    = CHEASEdata['ISTR']*mu0/CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
    CHEASEdata['JPHIN']    = CHEASEdata['JPHI']*mu0*CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
    CHEASEdata['ITORN']    = CHEASEdata['ITOR']*mu0/CHEASEdata['R0EXP']/CHEASEdata['B0EXP']
    CHEASEdata['IOHMICN']  = CHEASEdata['IPRLN']-CHEASEdata['IBSN']
    CHEASEdata['JOHMICN']  = CHEASEdata['JPRLN']-CHEASEdata['JBSN']
    CHEASEdata['PPrimeN']  = CHEASEdata['PPrime']*mu0*CHEASEdata['R0EXP']**2/CHEASEdata['B0EXP']
    CHEASEdata['TTPrimeN'] = CHEASEdata['TTPrime']/CHEASEdata['B0EXP']
 
    return CHEASEdata

def read_eqdsk(eqdskfpath,setParam={},**kwargs):
    if os.path.isfile(eqdskfpath) == False:
       errorFunc = traceback.extract_stack(limit=2)[-2][3]
       errorLine = traceback.extract_stack(limit=2)[-2][1]
       errorFile = traceback.extract_stack(limit=2)[-2][2]
       errMSG    = 'Call %s line %5d in file %s Failed.\n'
       errMSG   += 'Fatal: file %s not found.'
       raise IOError(errMSG %(errorFunc,errorLine,errorFile,eqdskfpath))
    
    rhopsiflag = False; rhotorflag = False
    if 'nrhomesh' in setParam:
        if   setParam['nrhomesh'] in [0,'rhopsi']: rhopsiflag = True
        elif setParam['nrhomesh'] in [1,'rhotor']: rhotorflag = True

    interpflag   = False
    cheaseflag   = False
    for key,value in kwargs.items():
        if   key in ['chease','cheasedata','cheasefpath']:
             if   type(value)==str and os.path.isfile(value.strip()):
                                       cheasedata = read_chease(cheasefpath=value.strip())
             elif type(value)==dict:   cheasedata = value.copy()
             if    'rhopsiN' in cheasedata and rhopsiflag:
                   interprho = cheasedata['PSIN'][:]; interpflag = True
             elif  'rhotorN' in cheasedata and rhotorflag:
                   interprho = cheasedata['PHIN'][:]; interpflag = True
             psi     = cheasedata['PSI'][:];    phi     = cheasedata['PHI'][:]
             rhopsiN = cheasedata['rhopsiN'][:];rhotorN = cheasedata['rhotorN'][:]
             cheaseflag = True

    EQDSKdata = read_efit_file(eqdskfpath)

    if interpflag:
       if   rhopsiflag:
            EQDSKdata['qpsi']     = interp(EQDSKdata['PSI'],EQDSKdata['qpsi'],psi) 
            EQDSKdata['fpol']     = interp(EQDSKdata['PSI'],EQDSKdata['fpol'],psi)
            EQDSKdata['pprime']   = interp(EQDSKdata['PSI'],EQDSKdata['pprime'],psi)
            EQDSKdata['ffprime']  = interp(EQDSKdata['PSI'],EQDSKdata['ffprime'],psi)
            EQDSKdata['pressure'] = interp(EQDSKdata['PSI'],EQDSKdata['pressure'],psi) 
            EQDSKdata['rhopsi']   = rhopsiN[:]
            EQDSKdata['rhotor']   = rhotorN[:]
       elif rhotorflag:
            EQDSKdata['qpsi']     = interp(EQDSKdata['PSI'],EQDSKdata['qpsi'],psi,phi,phi) 
            EQDSKdata['fpol']     = interp(EQDSKdata['PSI'],EQDSKdata['fpol'],psi,phi,phi)
            EQDSKdata['pprime']   = interp(EQDSKdata['PSI'],EQDSKdata['pprime'],psi,phi,phi)
            EQDSKdata['ffprime']  = interp(EQDSKdata['PSI'],EQDSKdata['ffprime'],psi,phi,phi)
            EQDSKdata['pressure'] = interp(EQDSKdata['PSI'],EQDSKdata['pressure'],psi,phi,phi)
            EQDSKdata['rhopsi']   = rhopsiN[:]
            EQDSKdata['rhotor']   = rhotorN[:]

    return EQDSKdata

def read_expeq(expeqfpath,setParam={},**kwargs):
    if os.path.isfile(expeqfpath) == False:
       errorFunc = traceback.extract_stack(limit=2)[-2][3]
       errorLine = traceback.extract_stack(limit=2)[-2][1]
       errorFile = traceback.extract_stack(limit=2)[-2][2]
       errMSG    = 'Call %s line %5d in file %s Failed.\n'
       errMSG   += 'Fatal: file %s not found.'
       raise IOError(errMSG %(errorFunc,errorLine,errorFile,expeqfpath))


    rhopsiflag = False; rhotorflag = False
    if 'nrhomesh' in setParam:
        if   setParam['nrhomesh'] in [0,'rhopsi']: rhopsiflag = True
        elif setParam['nrhomesh'] in [1,'rhotor']: rhotorflag = True

    eqdskflag    = False
    cheaseflag   = False
    interpflag   = False
    for key,value in kwargs.items():
        if   key in ['chease','cheasedata','cheasefpath']:
             if   type(value)==str and os.path.isfile(value.strip()):
                                       cheasedata = read_chease(cheasefpath=value.strip())
             elif type(value)==dict:   cheasedata = value.copy()
             if    'rhopsiN' in cheasedata and rhopsiflag:
                 interprho = cheasedata['rhopsiN'][:]; interpflag = True
             elif  'rhotorN' in cheasedata and rhotorflag:
                 interprho = cheasedata['rhotorN'][:]; interpflag = True
             rhopsi = cheasedata['rhopsiN']; rhotor = cheasedata['rhotorN'] 
             cheaseflag = True
        elif key in ['eqdsk','eqdskdata','eqdskfpath']:
             if   type(value)==str and os.path.isfile(value.strip()):
                                       eqdskdata = read_eqdsk(eqdskfpath=value.strip())
             elif  type(value)==dict:  eqdskdata = value.copy()
             if    'rhopsi' in eqdskdata and rhopsiflag:
                   interprho = eqdskdata['rhopsi'][:]; interpflag = True
             elif  'rhotor' in eqdskdata and rhotorflag:
                   interprho = eqdskdata['rhotor'][:]; interpflag = True
             rhopsi = eqdskdata['rhopsi']; rhotor = eqdskdata['rhotor'] 
             eqdskflag = True

    ofh = open(expeqfpath,'r')
    EXPEQOUT = ofh.readlines()
    ofh.close()

    EXPEQdata                  = {} 
    EXPEQdata['aspect']        = float(EXPEQOUT[0])
    EXPEQdata['zgeom']         = float(EXPEQOUT[1])
    EXPEQdata['pedge']         = float(EXPEQOUT[2])
    nRZmesh                    =   int(EXPEQOUT[3])
    EXPEQdata['rbound']        = npy.array([irec.split()[0] for irec in EXPEQOUT[4:nRZmesh+4]],dtype=float)
    EXPEQdata['zbound']        = npy.array([irec.split()[1] for irec in EXPEQOUT[4:nRZmesh+4]],dtype=float)
    
    nrhomesh                   = int(EXPEQOUT[nRZmesh+4].split()[0])
    EXPEQdata['nppfun']        = int(EXPEQOUT[nRZmesh+4].split()[1])
    EXPEQdata['nsttp']         = int(EXPEQOUT[nRZmesh+5].split()[0])
    EXPEQdata['nrhotype']      = int(EXPEQOUT[nRZmesh+5].split()[1])

    if   EXPEQdata['nrhotype'] == 0:
         EXPEQdata['rhopsiN']  = npy.array(EXPEQOUT[nRZmesh+6+0*nrhomesh:nRZmesh+6+1*nrhomesh],dtype=float)
    elif EXPEQdata['nrhotype'] == 1:
         EXPEQdata['rhotorN']  = npy.array(EXPEQOUT[nRZmesh+6+0*nrhomesh:nRZmesh+6+1*nrhomesh],dtype=float)
    
    if   EXPEQdata['nppfun']   == 4:
         EXPEQdata['pprimeN']  = npy.array(EXPEQOUT[nRZmesh+6+1*nrhomesh:nRZmesh+6+2*nrhomesh],dtype=float)
    elif EXPEQdata['nppfun']   == 8:
         EXPEQdata['pN']       = npy.array(EXPEQOUT[nRZmesh+6+1*nrhomesh:nRZmesh+6+2*nrhomesh],dtype=float)
  
    if   EXPEQdata['nsttp']    == 1:
         EXPEQdata['ttprimeN'] = npy.array(EXPEQOUT[nRZmesh+6+2*nrhomesh:nRZmesh+6+3*nrhomesh],dtype=float)
    elif EXPEQdata['nsttp']    == 2:
         EXPEQdata['istrN']    = npy.array(EXPEQOUT[nRZmesh+6+2*nrhomesh:nRZmesh+6+3*nrhomesh],dtype=float)
    elif EXPEQdata['nsttp']    == 3:
         EXPEQdata['iprlN']    = npy.array(EXPEQOUT[nRZmesh+6+2*nrhomesh:nRZmesh+6+3*nrhomesh],dtype=float)
    elif EXPEQdata['nsttp']    == 4:
         EXPEQdata['jprlN']    = npy.array(EXPEQOUT[nRZmesh+6+2*nrhomesh:nRZmesh+6+3*nrhomesh],dtype=float)
    elif EXPEQdata['nsttp']    == 5:
         EXPEQdata['q']        = npy.array(EXPEQOUT[nRZmesh+6+2*nrhomesh:nRZmesh+6+3*nrhomesh],dtype=float)

    if interpflag:
       if   rhopsiflag and EXPEQdata['nrhotype']==0:
            if   EXPEQdata['nppfun']   == 4:
                 EXPEQdata['pprimeN']  = interp(EXPEQdata['rhopsiN'],EXPEQdata['pprimeN'],rhopsi)
            elif EXPEQdata['nppfun']   == 8:
                 EXPEQdata['pN']       = interp(EXPEQdata['rhopsiN'],EXPEQdata['pN'],rhopsi)
            if   EXPEQdata['nsttp']    == 1:
                 EXPEQdata['ttprimeN'] = interp(EXPEQdata['rhopsiN'],EXPEQdata['ttprimeN'],rhopsi)
            elif EXPEQdata['nsttp']    == 2:
                 EXPEQdata['istrN']    = interp(EXPEQdata['rhopsiN'],EXPEQdata['istrN'],rhopsi)
            elif EXPEQdata['nsttp']    == 3:
                 EXPEQdata['iprlN']    = interp(EXPEQdata['rhopsiN'],EXPEQdata['iprlN'],rhopsi)
            elif EXPEQdata['nsttp']    == 4:
                 EXPEQdata['jprlN']    = interp(EXPEQdata['rhopsiN'],EXPEQdata['jprlN'],rhopsi)
            elif EXPEQdata['nsttp']    == 5:
                 EXPEQdata['q']        = interp(EXPEQdata['rhopsiN'],EXPEQdata['q'],rhopsi)
            EXPEQdata['rhopsiN']       = rhopsi[:]
            EXPEQdata['rhotorN']       = rhotor[:]
       elif rhotorflag and EXPEQdata['nrhotype']==1:
            if   EXPEQdata['nppfun']   == 4:
                 EXPEQdata['pprimeN']  = interp(EXPEQdata['rhotorN'],EXPEQdata['pprimeN'],rhotor)
            elif EXPEQdata['nppfun']   == 8:
                 EXPEQdata['pN']       = interp(EXPEQdata['rhotorN'],EXPEQdata['pN'],rhotor)
            if   EXPEQdata['nsttp']    == 1:
                 EXPEQdata['ttprimeN'] = interp(EXPEQdata['rhotorN'],EXPEQdata['ttprimeN'],rhotor)
            elif EXPEQdata['nsttp']    == 2:
                 EXPEQdata['istrN']    = interp(EXPEQdata['rhotorN'],EXPEQdata['istrN'],rhotor)
            elif EXPEQdata['nsttp']    == 3:
                 EXPEQdata['iprlN']    = interp(EXPEQdata['rhotorN'],EXPEQdata['iprlN'],rhotor)
            elif EXPEQdata['nsttp']    == 4:
                 EXPEQdata['jprlN']    = interp(EXPEQdata['rhotorN'],EXPEQdata['jprlN'],rhotor)
            elif EXPEQdata['nsttp']    == 5:
                 EXPEQdata['q']        = interp(EXPEQdata['rhotorN'],EXPEQdata['q'],rhotor)
            EXPEQdata['rhopsiN']       = rhopsi[:]
            EXPEQdata['rhotorN']       = rhotor[:]
       elif rhopsiflag and EXPEQdata['nrhotype']==1:
            if   EXPEQdata['nppfun']   == 4:
                 EXPEQdata['pprimeN']  = interp(EXPEQdata['rhotorN'],EXPEQdata['pprimeN'],rhotor,rhopsi,rhopsi)
            elif EXPEQdata['nppfun']   == 8:
                 EXPEQdata['pN']       = interp(EXPEQdata['rhotorN'],EXPEQdata['pN'],rhotor,rhopsi,rhopsi)
            if   EXPEQdata['nsttp']    == 1:
                 EXPEQdata['ttprimeN'] = interp(EXPEQdata['rhotorN'],EXPEQdata['ttprimeN'],rhotor,rhopsi,rhopsi)
            elif EXPEQdata['nsttp']    == 2:
                 EXPEQdata['istrN']    = interp(EXPEQdata['rhotorN'],EXPEQdata['istrN'],rhotor,rhopsi,rhopsi)
            elif EXPEQdata['nsttp']    == 3:
                 EXPEQdata['iprlN']    = interp(EXPEQdata['rhotorN'],EXPEQdata['iprlN'],rhotor,rhopsi,rhopsi)
            elif EXPEQdata['nsttp']    == 4:
                 EXPEQdata['jprlN']    = interp(EXPEQdata['rhotorN'],EXPEQdata['jprlN'],rhotor,rhopsi,rhopsi)
            elif EXPEQdata['nsttp']    == 5:
                 EXPEQdata['q']        = interp(EXPEQdata['rhotorN'],EXPEQdata['q'],rhotor,rhopsi,rhopsi)
            EXPEQdata['rhopsiN']       = rhopsi[:]
            EXPEQdata['rhotorN']       = rhotor[:]
       elif rhotorflag and EXPEQdata['nrhotype']==0:
            if   EXPEQdata['nppfun']   == 4:
                 EXPEQdata['pprimeN']  = interp(EXPEQdata['rhopsiN'],EXPEQdata['pprimeN'],rhopsi,rhotor,rhotor)
            elif EXPEQdata['nppfun']   == 8:
                 EXPEQdata['pN']       = interp(EXPEQdata['rhopsiN'],EXPEQdata['pN'],rhopsi,rhotor,rhotor)
            if   EXPEQdata['nsttp']    == 1:
                 EXPEQdata['ttprimeN'] = interp(EXPEQdata['rhopsiN'],EXPEQdata['ttprimeN'],rhopsi,rhotor,rhotor)
            elif EXPEQdata['nsttp']    == 2:
                 EXPEQdata['istrN']    = interp(EXPEQdata['rhopsiN'],EXPEQdata['istrN'],rhopsi,rhotor,rhotor)
            elif EXPEQdata['nsttp']    == 3:
                 EXPEQdata['iprlN']    = interp(EXPEQdata['rhopsiN'],EXPEQdata['iprlN'],rhopsi,rhotor,rhotor)
            elif EXPEQdata['nsttp']    == 4:
                 EXPEQdata['jprlN']    = interp(EXPEQdata['rhopsiN'],EXPEQdata['jprlN'],rhopsi,rhotor,rhotor)
            elif EXPEQdata['nsttp']    == 5:
                 EXPEQdata['q']        = interp(EXPEQdata['rhopsiN'],EXPEQdata['q'],rhopsi,rhotor,rhotor)
            EXPEQdata['rhopsiN']       = rhopsi[:]
            EXPEQdata['rhotorN']       = rhotor[:]

    return EXPEQdata

def write_expeq(setParam={},counter=0,outfile=True,**kwargs):
    '''
    nrhomesh=[rho_type(0:rhopsi,1:rhotor),rho_src(0:chease,1:eqdsk)]
    nppfun  =[pressure_type(8:pressure,4:pprime),pressure_src(0:chease,1:eqdsk,2:expeq,3:exptnz,4:iterdb,5:profiles)]
    nsttp   =[current_type(1:ttprime,2:istar,3:iprl,4:jprl,5:q),current_src(0:chease,1:eqdsk,2:expeq)]
    boundary= boundary_type(0:asis,1,interpolate)
    '''

    counter += 1

    eqdskflag    = False
    expeqflag    = False
    interpflag   = False
    cheaseflag   = False
    iterdbflag   = False
    exptnzflag   = False
    profilesflag = False
    for key,value in kwargs.items():
        if key in ['chease','cheasedata','cheasefpath']:
           if os.path.isfile(value.strip()):
              cheasepath = value.strip()
              cheasedata = read_chease(cheasefpath=cheasepath)
              cheaseflag = True

        if key in ['eqdsk','eqdskdata','eqdskfpath']:
           if os.path.isfile(value.strip()):
              eqdskpath = value.strip()
              eqdskdata = read_eqdsk(eqdskfpath=eqdskpath)
              eqdskflag = True

        if key in ['expeq','expeqdata','expeqfpath']:
           if os.path.isfile(value.strip()):
              expeqpath = value.strip()
              expeqdata = read_expeq(expeqfpath=expeqpath)
              expeqflag = True

        if key in ['exptnz','exptnzdata','exptnzfpath']:
           if os.path.isfile(value.strip()):
              exptnzpath = value.strip()
              exptnzflag = True

        if key in ['profiles','profilesdata','profilesfpath']:
           if os.path.isfile(value.strip()):
              profilespath = value.strip()
              profilesflag = True

        if key in ['iterdb','iterdbdata','iterdbfpath']:
           if os.path.isfile(value.strip()):
              iterdbpath = value.strip()
              iterdbflag = True


    if not (cheaseflag or expeqflag  or eqdskflag or exptnzflag or profilesflag or iterdbflag):
       raise IOError('FATAL: NO VALID INPUT PROFILES AVAILABLE. EXIT!')


    if 'nrhomesh' in setParam.keys():
       if   type(setParam['nrhomesh'])==list:
            if   type(setParam['nrhomesh'][0])==float: setParam['nrhomesh'][0] = int(setParam['nrhomesh'][0])
            elif type(setParam['nrhomesh'][0])==str:   setParam['nrhomesh'][0] = setParam['nrhomesh'][0].lower()
            if   type(setParam['nrhomesh'][1])==float: setParam['nrhomesh'][1] = int(setParam['nrhomesh'][1])
            elif type(setParam['nrhomesh'][1])==str:   setParam['nrhomesh'][1] = setParam['nrhomesh'][1].lower()
            nrhotype = setParam['nrhomesh'][:]
       elif type(setParam['nrhomesh']) in [int,str,float]:
            if   type(setParam['nrhomesh'])==float: setParam['nrhomesh'] = int(setParam['nrhomesh'])
            elif type(setParam['nrhomesh'])==str:   setParam['nrhomesh'] = setParam['nrhomesh'].lower()
            if   cheaseflag:   nrhotype=[setParam['nrhomesh'],0]
            elif eqdskflag:    nrhotype=[setParam['nrhomesh'],1]
    else:
            if   cheaseflag:   nrhotype=[0,0]
            elif eqdskflag:    nrhotype=[0,1]


    if 'nppfun' in setParam.keys():
       if   type(setParam['nppfun'])==list:
            if   type(setParam['nppfun'][0])==float: setParam['nppfun'][0] = int(setParam['nppfun'][0])
            elif type(setParam['nppfun'][0])==str:   setParam['nppfun'][0] = setParam['nppfun'][0].lower()
            if   type(setParam['nppfun'][1])==float: setParam['nppfun'][1] = int(setParam['nppfun'][1])
            elif type(setParam['nppfun'][1])==str:   setParam['nppfun'][1] = setParam['nppfun'][1].lower()
            nppfun=setParam['nppfun'][:]
       elif type(setParam['nppfun']) in [int,float,str]:
            if   type(setParam['nppfun'])==float: setParam['nppfun'] = int(setParam['nppfun'])
            elif type(setParam['nppfun'])==str:   setParam['nppfun'] = setParam['nppfun'].lower()
            if   cheaseflag:   nppfun=[setParam['nppfun'],0]
            elif eqdskflag:    nppfun=[setParam['nppfun'],1]
            elif expeqflag:    nppfun=[setParam['nppfun'],2]
            elif exptnzflag:   nppfun=[setParam['nppfun'],3]
            elif profilesflag: nppfun=[setParam['nppfun'],4]
            elif iterdbflag:   nppfun=[setParam['nppfun'],5]
    else:
            if   cheaseflag:   nppfun=[4,0]
            elif eqdskflag:    nppfun=[4,1]
            elif expeqflag:    nppfun=[4,2]
            elif exptnzflag:   nppfun=[4,3]
            elif profilesflag: nppfun=[4,4]
            elif iterdbflag:   nppfun=[4,5]

    if 'nsttp' in setParam.keys():
       if   type(setParam['nsttp'])==list:
            if   type(setParam['nsttp'][0])==float: setParam['nsttp'][0] = int(setParam['nsttp'][0])
            elif type(setParam['nsttp'][0])==str:   setParam['nsttp'][0] = setParam['nsttp'][0].lower()
            if   type(setParam['nsttp'][1])==float: setParam['nsttp'][1] = int(setParam['nsttp'][1])
            elif type(setParam['nsttp'][1])==str:   setParam['nsttp'][1] = setParam['nsttp'][1].lower()
            nsttp=setParam['nsttp'][:]
       elif type(setParam['nsttp']) in [int,float,str]:
            if   type(setParam['nsttp'])==float: setParam['nsttp'] = int(setParam['nsttp'])
            elif type(setParam['nsttp'])==str:   setParam['nsttp'] = setParam['nsttp'].lower()
            if   cheaseflag:   nsttp=[setParam['nsttp'],0]
            elif eqdskflag:    nsttp=[setParam['nsttp'],1]
            elif expeqflag:    nsttp=[setParam['nsttp'],2]
    else:
            if   cheaseflag:   nsttp=[1,0]
            elif eqdskflag:    nsttp=[1,1]
            elif expeqflag:    nsttp=[1,2]

    if 'boundary' in setParam.keys():
         if   type(setParam['boundary'])==float:
              boundary = int(setParam['boundary'])
         elif type(setParam['boundary'])==str:
              boundary = setParam['boundary'].lower()
    else:
              boundary = 0

    if   'cheasemode' in setParam.keys():
         cheasemode = int(setParam['cheasemode'])
    else:
         cheasemode = 1

    if   nrhotype[1] in [0,'chease'] and not cheaseflag:
         raise IOError('chease.h5 FILE IS NOT PROVIDED. EXIT!')
    elif nrhotype[1] in [1,'eqdsk'] and not eqdskflag:
         raise IOError('EQDSK FILE IS NOT PROVIDED. EXIT!')

    if nrhotype[0] in [1,'rhotor'] and nsttp[0] in [1,'ttprime','ttprimen']:
       raise NameError('FATAL: nrhotype (must) = 0 or rhopsi. Exit!')

    if   nsttp[1] in [0,'chease'] and not cheaseflag:
         raise IOError('chease.h5 FILE IS NOT PROVIDED. EXIT!')
    elif nsttp[1] in [1,'eqdsk'] and not eqdskflag:
         raise IOError('EQSSK FILE IS NOT PROVIDED. EXIT!')
    elif nsttp[1] in [2,'expeq'] and not expeqflag:
         raise IOError('EXPEQ FILE IS NOT PROVIDED. EXIT!')

    if   nppfun[1] in [0,'chease'] and not cheaseflag:
         raise IOError('chease.h5 FILE IS NOT PROVIDED. EXIT!')
    elif nppfun[1] in [1,'eqdsk'] and not eqdskflag:
         raise IOError('EQDSK FILE IS NOT PROVIDED. EXIT!')
    elif nppfun[1] in [2,'expeq'] and not expeqflag:
         raise IOError('EXPEQ FILE IS NOT PROVIDED. EXIT!')
    elif nppfun[1] in [3,'exptnz'] and not exptnzflag:
         raise IOError('EXPTNZ FILE IS NOT PROVIDED. EXIT!')
    elif nppfun[1] in [4,'profiles'] and not profilesflag:
         raise IOError('PROFILES FILE IS NOT PROVIDED. EXIT!')
    elif nppfun[1] in [5,'iterdb'] and not iterdbflag:
         raise IOError('ITERDB FILE IS NOT PROVIDED. EXIT!')

    currents = []
    currents.append([2,'istr','istrn','istar','istarn'])
    currents.append([3,'iprl','iprln','iparallel','iparalleln'])
    currents.append([4,'jprl','jprln','jparallel','jparalleln'])
    if (nsttp[1] in [1,'eqdsk']) and (nsttp[0] in currents):
       raise IOError("FATAL: eqdsk option is not accepted with nsttp = 2, 3 or 4")

    expeq = {}

    if   'cheasedata' in locals():
         expeq['R0EXP'] = cheasedata['R0EXP']
         expeq['B0EXP'] = cheasedata['B0EXP']
    elif 'eqdskdata' in locals():
         expeq['R0EXP'] = abs(eqdskdata['RCTR'])
         expeq['B0EXP'] = abs(eqdskdata['BCTR'])
    else:
         if   'R0EXP' in setParam.keys():
              expeq['R0EXP'] = setParam['R0EXP']
         else:
              expeq['R0EXP'] = 1.0
         if   'B0EXP' in setParam.keys():
              expeq['B0EXP'] = setParam['B0EXP']
         else:
              expeq['B0EXP'] = 1.0

    if   nrhotype[1] in [0,'chease'] and cheaseflag:
         expeq['nRZmesh'] = npy.size(cheasedata['rboundN'])
         expeq['rbound']  = cheasedata['rboundN'][:]
         expeq['zbound']  = cheasedata['zboundN'][:]

         expeq['aspect']  = (max(cheasedata['rbound'])-min(cheasedata['rbound']))
         expeq['aspect'] /= (max(cheasedata['rbound'])+min(cheasedata['rbound']))
         expeq['zgeom']   = npy.mean(cheasedata['zmesh'])/expeq['R0EXP']

         if   nrhotype[0] in [0,'rhopsi','rhopsin']:
              expeq['rhopsiN'] = cheasedata['rhopsiN'][:]
              expeq['nrhomesh']= npy.size(expeq['rhopsiN'])
              rhopsiflag = True
              rhotorflag = False
         elif nrhotype[0] in [1,'rhotor','rhotorn']:
              expeq['rhotorN'] = cheasedata['rhotorN'][:]
              expeq['nrhomesh']= npy.size(expeq['rhotorN'])
              rhopsiflag = False
              rhotorflag = True
    elif nrhotype[1] in [1,'eqdsk'] and eqdskflag:
         eqdskParam = {'boundary_type':boundary}
         rbound,zbound    = find_boundary(eqdskdata,setParam=eqdskParam)

         expeq['nRZmesh'] = npy.size(rbound)
         expeq['rbound']  = rbound/expeq['R0EXP']
         expeq['zbound']  = zbound/expeq['R0EXP']

         expeq['aspect']  = (max(rbound)-min(rbound))
         expeq['aspect'] /= (max(rbound)+min(rbound))
         expeq['zgeom']   = eqdskdata['ZMAX']/expeq['R0EXP']

         if   nrhotype[0] in [0,'rhopsi','rhopsin']:
              expeq['rhopsiN'] = eqdskdata['rhopsi'][:]
              expeq['nrhomesh']= npy.size(expeq['rhopsiN'])
              rhopsiflag = True
              rhotorflag = False
         elif nrhotype[0] in [1,'rhotor','rhotorn']:
              expeq['rhotorN'] = eqdskdata['rhotor'][:]
              expeq['nrhomesh']= npy.size(expeq['rhotorN'])
              rhopsiflag = False
              rhotorflag = True
    elif expeqflag:
         expeq['nRZmesh'] = npy.size(expeqdata['rbound'])
         expeq['rbound']  = expeqdata['rbound'][:]
         expeq['zbound']  = expeqdata['zbound'][:]

         expeq['aspect']  = expeqdata['aspect']
         expeq['zgeom']   = expeqdata['zgeom']

         if   eqdskflag:
              if   nrhotype[0] in [0,'rhopsi','rhopsin']:
                   expeq['rhopsiN'] = eqdskdata['rhopsi'][:]
                   expeq['nrhomesh']= npy.size(expeq['rhopsiN'])
                   rhopsiflag = True
                   rhotorflag = False
              elif nrhotype[0] in [1,'rhotor','rhotorn']:
                   expeq['rhotorN'] = eqdskdata['rhotor'][:]
                   expeq['nrhomesh']= npy.size(expeq['rhotorN'])
                   rhopsiflag = False
                   rhotorflag = True
         elif cheaseflag:
              if   nrhotype[0] in [0,'rhopsi','rhopsin']:
                   expeq['rhopsiN'] = cheasedata['rhopsiN'][:]
                   expeq['nrhomesh']= npy.size(expeq['rhopsiN'])
                   rhopsiflag = True
                   rhotorflag = False
              elif nrhotype[0] in [1,'rhotor','rhotorn']:
                   expeq['rhotorN'] = cheasedata['rhotorN'][:]
                   expeq['nrhomesh']= npy.size(expeq['rhotorN'])
                   rhopsiflag = False
                   rhotorflag = True

    SetParam = {'nrhomesh':nrhotype[0]}


    if   nppfun[1] in [0,'chease'] and cheaseflag:
         if   nrhotype[1] in [0,'chease']:
              cheasedata = read_chease(cheasefpath=cheasepath,setParam=SetParam)
         elif nrhotype[1] in [1,'eqdsk'] and eqdskflag:
              cheasedata = read_chease(cheasefpath=cheasepath,setParam=SetParam,eqdsk=eqdskpath)
         else:
              raise IOError('No adequate source of rhomesh available!')

         expeq['pedge']   = cheasedata['PN'][-1]

         if   nppfun[0] in [4,'pprime','pprimen']:
              expeq['PPrimeN'] = cheasedata['PPrimeN'][:]
         elif nppfun[0] in [8,'pressure','pressuren','p','pn']:
              expeq['pressureN'] = cheasedata['PN'][:]

    elif nppfun[1] in [1,'eqdsk'] and eqdskflag:
         if   nrhotype[1] in [0,'chease'] and cheaseflag:
              eqdskdata = read_eqdsk(eqdskfpath=eqdskpath,setParam=SetParam,chease=cheasepath)
         elif nrhotype[1] in [1,'eqdsk']:
              eqdskdata = read_eqdsk(eqdskfpath=eqdskpath,setParam=SetParam)
         else:
              raise IOError('No adequate source of rhomesh available!')

         expeq['pedge']   = mu0*eqdskdata['pressure'][-1]/expeq['B0EXP']**2

         if   nppfun[0] in [4,'pprime','pprimen']:
              expeq['PPrimeN'] = mu0*eqdskdata['pprime']*expeq['R0EXP']**2/expeq['B0EXP']
         elif nppfun[0] in [8,'pressure','pressuren','p','pn']:
              expeq['pressureN'] = mu0*eqdskdata['pressure']/expeq['B0EXP']**2

    elif nppfun[1] in [2,'expeq'] and expeqflag:
         if   nrhotype[1] in [0,'chease'] and cheaseflag:
              expeqdata = read_expeq(expeqfpath=expeqpath,setParam=SetParam,chease=cheasepath)
         elif nrhotype[1] in [1,'eqdsk']  and eqdskflag:
              expeqdata = read_expeq(expeqfpath=expeqpath,setParam=SetParam,eqdsk=eqdskpath)
         else:
              raise IOError('No adequate source of rhomesh available!')

         if   nrhotype[0] in [0,'rhopsi','rhopsin']:
              expeq['rhopsiN'] = expeqdata['rhopsiN'][:]
              expeq['nrhomesh']= npy.size(expeq['rhopsiN'])
         elif nrhotype[0] in [1,'rhotor','rhotorn']:
              expeq['rhotorN'] = expeqdata['rhotorN'][:]
              expeq['nrhomesh']= npy.size(expeq['rhotorN'])

         expeq['pedge']   = expeqdata['pedge']

         if   nppfun[0] in [4,'pprime','pprimen']:
              expeq['PPrimeN'] = expeqdata['pprimeN'][:]
         elif nppfun[0] in [8,'pressure','pressuren','p','pn']:
              expeq['pressureN'] = expeqdata['pN'][:]

    elif nppfun[1] in [3,'exptnz'] and exptnzflag:
         if   nrhotype[1] in [0,'chease'] and cheaseflag:
              exptnzdata = read_exptnz(exptnzfpath=exptnzpath,setParam=SetParam,chease=cheasepath)
         elif nrhotype[1] in [1,'eqdsk']  and eqdskflag:
              exptnzdata = read_exptnz(exptnzfpath=exptnzpath,setParam=SetParam,eqdsk=eqdskpath)
         else:
              raise IOError('No adequate source of rhomesh available!')

         if   nrhotype[0] in [0,'rhopsi','rhopsin']:
              expeq['rhopsiN'] = exptnzdata['rhopsiN'][:]
              expeq['nrhomesh']= npy.size(expeq['rhopsiN'])
         elif nrhotype[0] in [1,'rhotor','rhotorn']:
              expeq['rhotorN'] = exptnzdata['rhotorN'][:]
              expeq['nrhomesh']= npy.size(expeq['rhotorN'])

         expeq['pedge']   = mu0*exptnzdata['pressure'][-1]/expeq['B0EXP']**2

         if   nppfun[0] in [4,'pprime','pprimen']:
              expeq['PPrimeN'] = mu0*exptnzdata['pprime']*expeq['R0EXP']**2/expeq['B0EXP']
         elif nppfun[0] in [8,'pressure','pressuren','p','pn']:
              expeq['pressureN'] = mu0*exptnzdata['pressure']/expeq['B0EXP']**2

    elif nppfun[1] in [4,'profiles'] and profilesflag:
         if   nrhotype[1] in [0,'chease'] and cheaseflag:
              profilesdata = read_profiles(profilesfpath=profilespath,setParam=SetParam,chease=cheasepath)
         elif nrhotype[1] in [1,'eqdsk']  and eqdskflag:
              profilesdata = read_profiles(profilesfpath=profilespath,setParam=SetParam,eqdsk=eqdskpath)
         else:
              raise IOError('No adequate source of rhomesh available!')

         if   nrhotype[0] in [0,'rhopsi','rhopsin']:
              expeq['rhopsiN'] = profilesdata['rhopsiN'][:]
              expeq['nrhomesh']= npy.size(expeq['rhopsiN'])
         elif nrhotype[0] in [1,'rhotor','rhotorn']:
              expeq['rhotorN'] = profilesdata['rhotorN'][:]
              expeq['nrhomesh']= npy.size(expeq['rhotorN'])

         expeq['pedge']   = mu0*profilesdata['pressure'][-1]/expeq['B0EXP']**2

         if   nppfun[0] in [4,'pprime','pprimen']:
              expeq['PPrimeN'] = mu0*profilesdata['pprime']*expeq['R0EXP']**2/expeq['B0EXP']
         elif nppfun[0] in [8,'pressure','pressuren','p','pn']:
              expeq['pressureN'] = mu0*profilesdata['pressure']/expeq['B0EXP']**2

    elif nppfun[1] in [5,'iterdb'] and iterdbflag:
         if   nrhotype[1] in [0,'chease'] and cheaseflag:
              iterdbdata = read_iterdb(iterdbfpath=iterdbpath,setParam=SetParam,chease=cheasepath)
         elif nrhotype[1] in [1,'eqdsk']  and eqdskflag:
              iterdbdata = read_iterdb(iterdbfpath=iterdbpath,setParam=SetParam,eqdsk=eqdskpath)

         if   nrhotype[0] in [0,'rhopsi','rhopsin']:
              expeq['rhopsiN'] = iterdbdata['rhopsiN'][:]
              expeq['nrhomesh']= npy.size(expeq['rhopsiN'])
         elif nrhotype[0] in [1,'rhotor','rhotorn']:
              expeq['rhotorN'] = iterdbdata['rhotorN'][:]
              expeq['nrhomesh']= npy.size(expeq['rhotorN'])

         expeq['pedge']   = mu0*iterdbdata['pressure'][-1]/expeq['B0EXP']**2

         if   nppfun[0] in [4,'pprime','pprimen']:
              expeq['PPrimeN'] = mu0*iterdbdata['pprime']*expeq['R0EXP']**2/expeq['B0EXP']
         elif nppfun[0] in [8,'pressure','pressuren','p','pn']:
              expeq['pressureN'] = mu0*iterdbdata['pressure']/expeq['B0EXP']**2
         else:
              raise IOError('No adequate source of rhomesh available!')

    if   nsttp[1] in [0,'chease'] and cheaseflag:
         if   nrhotype[1] in [0,'chease'] and cheaseflag:
              cheasedata = read_chease(cheasefpath=cheasepath,setParam=SetParam)
         elif nrhotype[1] in [1,'eqdsk'] and eqdskflag:
              cheasedata = read_chease(cheasefpath=cheasepath,setParam=SetParam,eqdsk=eqdskpath)

         if   cheasemode == 1:
              if   nsttp[0] in [1,'ttprime','ttprimen']:
                   expeq['TTPrimeN'] = cheasedata['TTPrimeN'][:]
              elif nsttp[0] in [2,'istr','istrn','istar','istarn']:
                   expeq['ISTRN'] = cheasedata['ISTRN'][:]
              elif nsttp[0] in [3,'iprl','iprln','iparallel','iparalleln']:
                   expeq['IPRLN'] = cheasedata['IPRLN'][:]
              elif nsttp[0] in [4,'jprl','jprln','jparallel','jparalleln']:
                   expeq['JPRLN'] = cheasedata['JPRLN'][:]
              elif nsttp[0] in [5,'q','qpsi','qtor']:
                   expeq['q'] = cheasedata['q'][:]
         elif cheasemode == 2:
              correctParam = {'nsttp':nsttp,'ITEXP':setParam['ITEXP']}
              correction   = current_correction(chease=cheasedata,setParam=correctParam)
              if   nsttp[0] in [1,'ttprime','ttprimen']:
                   expeq['TTPrimeN'] = correction[:]
              elif nsttp[0] in [2,'istr','istrn','istar','istarn']:
                   expeq['ISTRN'] = correction[:]
              elif nsttp[0] in [3,'iprl','iprln','iparallel','iparalleln']:
                   expeq['IPRLN'] = correction[:]
              elif nsttp[0] in [4,'jprl','jprln','jparallel','jparalleln']:
                   expeq['JPRLN'] = correction[:]
    elif nsttp[1] in [2,'expeq'] and expeqflag:
         if   nrhotype[1] in [0,'chease'] and cheaseflag:
              expeqdata = read_expeq(expeqfpath=expeqpath,setParam=SetParam,chease=cheasepath)
         elif nrhotype[1] in [1,'eqdsk'] and eqdskflag:
              expeqdata = read_expeq(expeqfpath=expeqpath,setParam=SetParam,eqdsk=eqdskpath)

         if   cheasemode == 1:
              if   nsttp[0] in [1,'ttprime','ttprimen']:
                   expeq['TTPrimeN'] = expeqdata['ttprimeN'][:]
              elif nsttp[0] in [2,'istr','istrn','istar','istarn']:
                   expeq['ISTRN'] = expeqdata['istrN'][:]
              elif nsttp[0] in [3,'iprl','iprln','iparallel','iparalleln']:
                   expeq['IPRLN'] = expeqdata['iprlN'][:]
              elif nsttp[0] in [4,'jprl','jprln','jparallel','jparalleln']:
                   expeq['JPRLN'] = expeqdata['jprlN'][:]
              elif nsttp[0] in [5,'q','qpsi','qtor']:
                   expeq['q'] = expeqdata['q'][:]
         elif cheasemode == 2:
              correctParam = {'nsttp':nsttp,'ITEXP':setParam['ITEXP']}
              correction   = current_correction(chease=cheasedata,expeq=expeqdata,setParam=correctParam)
              if   nsttp[0] in [1,'ttprime','ttprimen']:
                   expeq['TTPrimeN'] = correction[:]
              elif nsttp[0] in [2,'istr','istrn','istar','istarn']:
                   expeq['ISTRN'] = correction[:]
              elif nsttp[0] in [3,'iprl','iprln','iparallel','iparalleln']:
                   expeq['IPRLN'] = correction[:]
              elif nsttp[0] in [4,'jprl','jprln','jparallel','jparalleln']:
                   expeq['JPRLN'] = correction[:]
    elif nsttp[1] in [1,'eqdsk'] and eqdskflag:
         if   nrhotype[1] in [0,'chease'] and cheaseflag:
              eqdskdata = read_eqdsk(eqdskfpath=eqdskpath,setParam=SetParam,chease=cheasepath)
         elif nrhotype[1] in [1,'eqdsk']:
              eqdskdata = read_eqdsk(eqdskfpath=eqdskpath,setParam=SetParam)

         if   cheasemode==1:
              if   nsttp[0] in [1,'ttprime','ttprimen']:
                   expeq['TTPrimeN'] = eqdskdata['ffprime']/expeq['B0EXP']
              elif nsttp[0] in [5,'q','qpsi','qtor']:
                   expeq['q'] = eqdskdata['qpsi']
              else:
                   raise ValueError("FATAL: nsttp[0] (must) = 1 or 5.")
         elif cheasemode in [2,3]:
              raise ValueError("FATAL: eqdsk option is not accepted with cheasemode = 2 or 3")
            
    if outfile:
       ofh = open("EXPEQ",'w')
       ofh.write('%18.8E\n'           % expeq['aspect'])
       ofh.write('%18.8E\n'           % expeq['zgeom'])
       ofh.write('%18.8E\n'           % expeq['pedge'])
       ofh.write('%5d\n'              % expeq['nRZmesh'])
       for i in range(expeq['nRZmesh']):
           ofh.write('%18.8E%18.8E\n' % (expeq['rbound'][i],expeq['zbound'][i]))
       ofh.write('%5d%5d\n'           % (expeq['nrhomesh'],nppfun[0]))
       ofh.write('%5d%5d\n'           % (nsttp[0],nrhotype[0]))


       if   nrhotype[0] in [0,'rhopsi']:
            for i in range(expeq['nrhomesh']):
                ofh.write('%18.8E\n'       % expeq['rhopsiN'][i])
       elif nrhotype[0] in [1,'rhotor']:
            for i in range(expeq['nrhomesh']):
                ofh.write('%18.8E\n'       % expeq['rhotorN'][i])

       if   nppfun[0] in [4,'pprime','pprimen']:
            for i in range(expeq['nrhomesh']):
                ofh.write('%18.8E\n'       % expeq['PPrimeN'][i])
       elif nppfun[0] in [8,'pressure','pressuren','p','pn']:
            for i in range(expeq['nrhomesh']):
                ofh.write('%18.8E\n'       % expeq['pressureN'][i])

       if   nsttp[0] in [1,'ttprime','ttprimen']:
            for i in range(expeq['nrhomesh']):
                ofh.write('%18.8E\n'       % expeq['TTPrimeN'][i])
       elif nsttp[0] in [2,'istar','istarn','istr','istrn']:
            for i in range(expeq['nrhomesh']):
                ofh.write('%18.8E\n'       % expeq['ISTRN'][i])
       elif nsttp[0] in [3,'iprl','iprln','iparallel','iparalleln']:
            for i in range(expeq['nrhomesh']):
                ofh.write('%18.8E\n'       % expeq['IPRLN'][i])
       elif nsttp[0] in [4,'jprl','jprln','jparallel','jparalleln']:
            for i in range(expeq['nrhomesh']):
                ofh.write('%18.8E\n'       % expeq['JPRLN'][i])
       elif nsttp[0] in [5,'q']:
            for i in range(expeq['nrhomesh']):
                ofh.write('%18.8E\n'       % expeq['q'][i])
       ofh.close()

    return expeq

def find_boundary(eqdsk='',setParam={}):
    if eqdsk:
       eqdskflag  = True
       if   type(eqdsk)==str and os.path.isfile(eqdsk.strip()):
                               eqdskdata = read_eqdsk(eqdskfpath=eqdsk.strip())
       elif type(eqdsk)==dict: eqdskdata = eqdsk.copy()
       else:
            eqdskflag = False
    else:
       eqdskflag = False

    if not eqdskflag: raise IOError('FATAL: EQDSK FILE IS NOT PROVIDED. EXIT!')

    asisflag = False; interpflag = False
    if 'boundary_type' in setParam:
       if   setParam['boundary_type'] in [0,'asis']:   asisflag   = True
       elif setParam['boundary_type'] in [1,'interp']: interpflag = True
       else:                                           asisflag   = True

    if   asisflag:
         rbound = eqdskdata['rbound']
         zbound = eqdskdata['zbound']
    elif interpflag:
       rbndtst = int(eqdskdata['RLEN']/(max(eqdskdata['rbound'])-abs(min(eqdskdata['rbound']))))
       zbndtst = int(eqdskdata['ZLEN']/(max(eqdskdata['zbound'])+abs(min(eqdskdata['zbound']))))
       if  rbndtst==1 and zbndtst==1:
           rbound,zbound = magsurf_solvflines(eqdskdata=eqdskdata,psi=0.9999,eps=1.0e-16)
       else:
           rbound=npy.zeros(2*len(eqdskdata['rbound'])-1)
           zbound=npy.zeros(2*len(eqdskdata['zbound'])-1)
           rbound[0] = eqdskdata['rbound'][0]
           zbound[0] = eqdskdata['zbound'][0]
           for i in range(1,len(eqdskdata['rbound'])):
               rbound[i]  = eqdskdata['rbound'][i]
               rbound[-i] = eqdskdata['rbound'][i]
               zbound[i]  = eqdskdata['zbound'][i]
               zbound[-i] =-eqdskdata['zbound'][i]

    return rbound,zbound

def read_exptnz(exptnzfpath,setParam={},**kwargs):
    if os.path.isfile(exptnzfpath) == False:
       errorFunc = traceback.extract_stack(limit=2)[-2][3]
       errorLine = traceback.extract_stack(limit=2)[-2][1]
       errorFile = traceback.extract_stack(limit=2)[-2][2]
       errMSG    = 'Call %s line %5d in file %s Failed.\n'
       errMSG   += 'Fatal: file %s not found.'
       raise IOError(errMSG %(errorFunc,errorLine,errorFile,exptnzfpath))

    rhopsiflag = False; rhotorflag = False
    if 'nrhomesh' in setParam:
        if   setParam['nrhomesh'] in [0,'rhopsi']: rhopsiflag = True
        elif setParam['nrhomesh'] in [1,'rhotor']: rhotorflag = True

    eqdskflag    = False
    cheaseflag   = False
    interpflag   = False
    for key,value in kwargs.items():
        if   key in ['chease','cheasedata','cheasefpath']:
             if   type(value)==str and os.path.isfile(value.strip()):
                                       cheasedata = read_chease(cheasefpath=value.strip())
             elif type(value)==dict:   cheasedata = value.copy()
             if    'rhopsiN' in cheasedata and rhopsiflag:
                   interprho = cheasedata['rhopsiN'][:]; interpflag = True
             elif  'rhotorN' in cheasedata and rhotorflag:
                   interprho = cheasedata['rhotorN'][:]; interpflag = True
             rhopsi = cheasedata['rhopsiN']; rhotor = cheasedata['rhotorN'] 
             cheaseflag = True
        elif key in ['eqdsk','eqdskdata','eqdskfpath']:
             if   type(value)==str and os.path.isfile(value.strip()):
                                       eqdskdata = read_eqdsk(eqdskfpath=value.strip())
             elif  type(value)==dict:  eqdskdata = value.copy()
             if    'rhopsi' in eqdskdata and rhopsiflag:
                   interprho = eqdskdata['rhopsi'][:]; interpflag = True
             elif  'rhotor' in eqdskdata and rhotorflag:
                   interprho = eqdskdata['rhotor'][:]; interpflag = True
             rhopsi = eqdskdata['rhopsi']; rhotor = eqdskdata['rhotor']
             eqdskflag = True

    ofh = open(exptnzfpath,'r')
    EXPTNZOUT = ofh.readlines()
    ofh.close()

    n_rho   = int(EXPTNZOUT[0].split()[0])
    rhotype =     EXPTNZOUT[0].split()[1].strip()[0:6]

    EXPTNZdata               = {}

    EXPTNZdata['Zi'] = 1.0
    EXPTNZdata['Zz'] = 6.0

    if   rhotype=='rhopsi':
         EXPTNZdata['rhopsiN'] = npy.array(EXPTNZOUT[0*n_rho+1:1*n_rho+1],dtype=float)
    elif rhotype=='rhotor':
         EXPTNZdata['rhotorN'] = npy.array(EXPTNZOUT[0*n_rho+1:1*n_rho+1],dtype=float)
    EXPTNZdata['Te']           = npy.array(EXPTNZOUT[1*n_rho+1:2*n_rho+1],dtype=float)
    EXPTNZdata['ne']           = npy.array(EXPTNZOUT[2*n_rho+1:3*n_rho+1],dtype=float)
    EXPTNZdata['Zeff']         = npy.array(EXPTNZOUT[3*n_rho+1:4*n_rho+1],dtype=float)
    EXPTNZdata['Ti']           = npy.array(EXPTNZOUT[4*n_rho+1:5*n_rho+1],dtype=float)
    EXPTNZdata['ni']           = npy.array(EXPTNZOUT[5*n_rho+1:6*n_rho+1],dtype=float)

    if interpflag:
       if   rhopsiflag and rhotype=='rhopsi':
            EXPTNZdata['Te']      = interp(EXPTNZdata['rhopsiN'],EXPTNZdata['Te'],rhopsi)
            EXPTNZdata['Ti']      = interp(EXPTNZdata['rhopsiN'],EXPTNZdata['Ti'],rhopsi)
            EXPTNZdata['ne']      = interp(EXPTNZdata['rhopsiN'],EXPTNZdata['ne'],rhopsi) 
            EXPTNZdata['ni']      = interp(EXPTNZdata['rhopsiN'],EXPTNZdata['ni'],rhopsi)
            EXPTNZdata['Zeff']    = interp(EXPTNZdata['rhopsiN'],EXPTNZdata['Zeff'],rhopsi)
            EXPTNZdata['rhopsiN'] = rhopsi[:]
            EXPTNZdata['rhotorN'] = rhotor[:]
       elif rhotorflag and rhotype=='rhotor':
            EXPTNZdata['Te']      = interp(EXPTNZdata['rhotorN'],EXPTNZdata['Te'],rhotor)
            EXPTNZdata['Ti']      = interp(EXPTNZdata['rhotorN'],EXPTNZdata['Ti'],rhotor)
            EXPTNZdata['ne']      = interp(EXPTNZdata['rhotorN'],EXPTNZdata['ne'],rhotor) 
            EXPTNZdata['ni']      = interp(EXPTNZdata['rhotorN'],EXPTNZdata['ni'],rhotor)
            EXPTNZdata['Zeff']    = interp(EXPTNZdata['rhotorN'],EXPTNZdata['Zeff'],rhotor)
            EXPTNZdata['rhopsiN'] = rhopsi[:]
            EXPTNZdata['rhotorN'] = rhotor[:]
       elif rhotorflag and rhotype=='rhopsi':
            EXPTNZdata['Te']      = interp(EXPTNZdata['rhopsiN'],EXPTNZdata['Te'],rhopsi,rhotor,rhotor)
            EXPTNZdata['Ti']      = interp(EXPTNZdata['rhopsiN'],EXPTNZdata['Ti'],rhopsi,rhotor,rhotor)
            EXPTNZdata['ne']      = interp(EXPTNZdata['rhopsiN'],EXPTNZdata['ne'],rhopsi,rhotor,rhotor) 
            EXPTNZdata['ni']      = interp(EXPTNZdata['rhopsiN'],EXPTNZdata['ni'],rhopsi,rhotor,rhotor)
            EXPTNZdata['Zeff']    = interp(EXPTNZdata['rhopsiN'],EXPTNZdata['Zeff'],rhopsi,rhotor,rhotor)
            EXPTNZdata['rhopsiN'] = rhopsi[:]
            EXPTNZdata['rhotorN'] = rhotor[:]
       elif rhopsiflag and rhotype=='rhotor':
            EXPTNZdata['Te']      = interp(EXPTNZdata['rhotorN'],EXPTNZdata['Te'],rhotor,rhopsi,rhopsi)
            EXPTNZdata['Ti']      = interp(EXPTNZdata['rhotorN'],EXPTNZdata['Ti'],rhotor,rhopsi,rhopsi)
            EXPTNZdata['ne']      = interp(EXPTNZdata['rhotorN'],EXPTNZdata['ne'],rhotor,rhopsi,rhopsi) 
            EXPTNZdata['ni']      = interp(EXPTNZdata['rhotorN'],EXPTNZdata['ni'],rhotor,rhopsi,rhopsi)
            EXPTNZdata['Zeff']    = interp(EXPTNZdata['rhotorN'],EXPTNZdata['Zeff'],rhotor,rhopsi,rhopsi)
            EXPTNZdata['rhopsiN'] = rhopsi[:]
            EXPTNZdata['rhotorN'] = rhotor[:]
    elif rhopsiflag and rhotype=='rhotor':
         print("WARNING: setParam['nrhomesh'] = 0 or rhopsi, but the path to a target rhopsi is not provided.")
         print("         Converting the profiles to poloidal (psi) coordinates could not be done, and")
         print("         all profiles are provided in the toroidal (phi) coordinates.")
    elif rhotorflag and rhotype=='rhopsi':
         print("WARNING: setParam['nrhomesh'] = 1 or rhotor, but the path to a target rhotor is not provided.")
         print("         Converting the profiles to toroidal (phi) coordinates could not be done, and")
         print("         all profiles are provided in the poloidal (psi) coordinates.")

    EXPTNZdata['nz']         = EXPTNZdata['Zeff']*EXPTNZdata['ne']
    EXPTNZdata['nz']        -= EXPTNZdata['ni']*EXPTNZdata['Zi']**2
    EXPTNZdata['nz']        /= EXPTNZdata['Zz']**2

    EXPTNZdata['pressure']   = EXPTNZdata['Te']*EXPTNZdata['ne']
    EXPTNZdata['pressure']  += EXPTNZdata['Ti']*EXPTNZdata['ni']
    EXPTNZdata['pressure']  += EXPTNZdata['Ti']*EXPTNZdata['nz']
    EXPTNZdata['pressure']  *= 1.602e-19

    if   rhopsiflag:
         EXPTNZdata['pprime']= derivative(x=EXPTNZdata['rhopsiN'],fx=EXPTNZdata['pressure'],method='CubicSpline')
    elif rhotorflag:
         EXPTNZdata['pprime']= derivative(x=EXPTNZdata['rhopsiN'],fx=EXPTNZdata['pressure'],method='CubicSpline')
    elif rhotype=='rhopsi':
         EXPTNZdata['pprime']= derivative(x=EXPTNZdata['rhopsiN'],fx=EXPTNZdata['pressure'],method='CubicSpline')
    elif rhotype=='rhotor':
         EXPTNZdata['pprime']= derivative(x=EXPTNZdata['rhotorN'],fx=EXPTNZdata['pressure'],method='CubicSpline')

    return EXPTNZdata

def write_exptnz(setParam={},outfile=True,**kwargs):
    '''
    nrhomesh=[rho_type(0:rhopsi,1:rhotor),rho_src(0:chease,1:eqdsk)]
    eprofile=[eprofile_src(0:chease,3:exptnz,4:profiles,5:iterdb)]
    iprofile=[iprofile_src(0:chease,3:exptnz,4:profiles,5:iterdb)]
    '''
    eqdskflag    = False
    interpflag   = False
    cheaseflag   = False
    iterdbflag   = False
    exptnzflag   = False
    profilesflag = False
    for key,value in kwargs.items():
        if key in ['chease','cheasedata','cheasefpath']:
           if os.path.isfile(value.strip()):
              cheasepath = value.strip()
              cheaseflag = True

        if key in ['eqdsk','eqdskdata','eqdskfpath']:
           if os.path.isfile(value.strip()):
              eqdskpath = value.strip()
              eqdskflag = True

        if key in ['expeq','expeqdata','expeqfpath']:
           if os.path.isfile(value.strip()):
              expeqpath = value.strip()
              expeqflag = True

        if key in ['exptnz','exptnzdata','exptnzfpath']:
           if os.path.isfile(value.strip()):
              exptnzpath = value.strip()
              exptnzflag = True

        if key in ['profiles','profilesdata','profilesfpath']:
           if os.path.isfile(value.strip()):
              profilespath = value.strip()
              profilesflag = True

        if key in ['iterdb','iterdbdata','iterdbfpath']:
           if os.path.isfile(value.strip()):
              iterdbpath = value.strip()
              iterdbflag = True

    if not (cheaseflag or exptnzflag or profilesflag or iterdbflag):
       raise IOError('FATAL: NO VALID INPUT PROFILES AVAILABLE. EXIT!')

    if 'nrhomesh' in setParam.keys():
       if   type(setParam['nrhomesh'])==list: 
            if   type(setParam['nrhomesh'][0])==float: setParam['nrhomesh'][0] = int(setParam['nrhomesh'][0])
            elif type(setParam['nrhomesh'][0])==str:   setParam['nrhomesh'][0] = setParam['nrhomesh'][0].lower()
            if   type(setParam['nrhomesh'][1])==float: setParam['nrhomesh'][1] = int(setParam['nrhomesh'][1])
            elif type(setParam['nrhomesh'][1])==str:   setParam['nrhomesh'][1] = setParam['nrhomesh'][1].lower()
            nrhotype = setParam['nrhomesh'][:]
       elif type(setParam['nrhomesh'])==int:
            if   cheaseflag:   nrhotype=[setParam['nrhomesh'],0]
            elif eqdskflag:    nrhotype=[setParam['nrhomesh'],1]
    else:
            if   cheaseflag:   nrhotype=[0,0]
            elif eqdskflag:    nrhotype=[0,1]

    if 'eprofile' in setParam.keys(): 
            if   type(setParam['eprofile'])==float: setParam['eprofile'] = int(setParam['eprofile'])
            elif type(setParam['eprofile'])==str:   setParam['eprofile'] = setParam['eprofile'].lower()
            eprofile= setParam['eprofile']
    else:
            if   cheaseflag:   eprofile=0
            elif exptnzflag:   eprofile=3
            elif profilesflag: eprofile=4
            elif iterdbflag:   eprofile=5

    if 'iprofile' in setParam.keys():
            if   type(setParam['iprofile'])==float: setParam['iprofile'] = int(setParam['iprofile'])
            elif type(setParam['iprofile'])==str:   setParam['iprofile'] = setParam['iprofile'].lower()
            iprofile= setParam['iprofile']
    else:
             if   cheaseflag:   iprofile=0
             elif exptnzflag:   iprofile=3
             elif profilesflag: iprofile=4
             elif iterdbflag:   iprofile=5

    if   nrhotype[1] in [0,'chease'] and not cheaseflag:
         raise IOError("FATAL: nrhotype=chease and chease.h5 file is not available. EXIT!")
    elif nrhotype[1] in [1,'eqdsk'] and not eqdskflag:
         raise IOError("FATAL: nrhotype=eqdsk and eqdsk file is not available. EXIT!")

    if   eprofile in [0,'chease'] and not cheaseflag:
         raise IOError("FATAL: eprofile=chease and chease.h5 is not available. EXIT!")
    elif eprofile in [3,'exptnz'] and not exptnzflag:
         raise IOError("FATAL: eprofile=exptnz and exptnz is not available. EXIT!")
    elif eprofile in [4,'profiles'] and not profilesflag:
         raise IOError("FATAL: eprofile=profiples and profiles is not available. EXIT!")
    elif eprofile in [5,'iterdb'] and not iterdbflag:
         raise IOError("FATAL: eprofile=iterdb and iterdb is not available. EXIT!")

    if   iprofile in [0,'chease'] and not cheaseflag:
         raise IOError("FATAL: iprofile=chease and chease.h5 is not available. EXIT!")
    elif iprofile in [3,'exptnz'] and not exptnzflag:
         raise IOError("FATAL: iprofile=exptnz and exptnz is not available. EXIT!")
    elif iprofile in [4,'profiles'] and not profilesflag:
         raise IOError("FATAL: iprofile=profiples and profiles is not available. EXIT!")
    elif iprofile in [5,'iterdb'] and not iterdbflag:
         raise IOError("FATAL: iprofile=iterdb and iterdb is not available. EXIT!")

    if   nrhotype[0] in [0,'rhopsi']:
         rhopsiflag = True
         rhotorflag = False
    elif nrhotype[0] in [1,'rhotor']:
         rhopsiflag = False
         rhotorflag = True

    SetParam={'nrhomesh':nrhotype[0]}

    exptnz = {}

    if   eprofile in [0,'chease'] and cheaseflag:
         if   nrhotype[1] in [0,'chease']:
              cheasedata = read_chease(cheasefpath=cheasepath,setParam=SetParam)
         elif nrhotype[1] in [1,'eqdsk'] and eqdskflag:
              cheasedata = read_chease(cheasefpath=cheasepath,setParam=SetParam,eqdsk=eqdskpath)
         if   rhopsiflag: exptnz['rhopsiN'] = cheasedata['rhopsiN']
         elif rhotorflag: exptnz['rhotorN'] = cheasedata['rhotorN']
         exptnz['Te'] = cheasedata['Te']
         exptnz['ne'] = cheasedata['ne']

    elif eprofile in [3,'exptnz'] and exptnzflag:
         if   nrhotype[1] in [0,'chease'] and cheaseflag:
              exptnzdata = read_exptnz(exptnzfpath=exptnzpath,setParam=SetParam,chease=cheasepath)
         elif nrhotype[1] in [1,'eqdsk'] and eqdskflag:
              exptnzdata = read_exptnz(exptnzfpath=exptnzpath,setParam=SetParam,eqdsk=eqdskpath)
         if   rhopsiflag: exptnz['rhopsiN'] = exptnzdata['rhopsiN']
         elif rhotorflag: exptnz['rhotorN'] = exptnzdata['rhotorN']
         exptnz['Te'] = exptnzdata['Te']
         exptnz['ne'] = exptnzdata['ne']

    elif eprofile in [4,'profiles'] and profilesflag:
         if   nrhotype[1] in [0,'chease'] and cheaseflag:
              profilesdata = read_profiles(profilesfpath=profilespath,setParam=SetParam,chease=cheasepath)
         elif nrhotype[1] in [1,'eqdsk'] and eqdskflag:
              profilesdata = read_profiles(profilesfpath=profilespath,setParam=SetParam,eqdsk=eqdskpath)
         if   rhopsiflag: exptnz['rhopsiN'] = profilesdata['rhopsiN']
         elif rhotorflag: exptnz['rhotorN'] = profilesdata['rhotorN']
         exptnz['Te'] = profilesdata['Te']
         exptnz['ne'] = profilesdata['ne']

    elif eprofile in [5,'iterdb'] and iterdbflag:
         if   nrhotype[1] in [0,'chease'] and cheaseflag:
              iterdbdata = read_iterdb(iterdbfpath=iterdbpath,setParam=SetParam,chease=cheasepath)
         elif nrhotype[1] in [1,'eqdsk'] and eqdskflag:
              iterdbdata = read_iterdb(iterdbfpath=iterdbpath,setParam=SetParam,eqdsk=eqdskpath)
         if   rhopsiflag: exptnz['rhopsiN'] = iterdbdata['rhopsiN']
         elif rhotorflag: exptnz['rhotorN'] = iterdbdata['rhotorN']
         exptnz['Te'] = iterdbdata['Te']
         exptnz['ne'] = iterdbdata['ne']

    if   iprofile in [0,'chease'] and cheaseflag:
         if   nrhotype[1] in [0,'chease']:
              cheasedata = read_chease(cheasefpath=cheasepath,setParam=SetParam)
         elif nrhotype[1] in [1,'eqdsk'] and eqdskflag:
              cheasedata = read_chease(cheasefpath=cheasepath,setParam=SetParam,eqdsk=eqdskpath)
         exptnz['Ti']   = cheasedata['Ti']
         exptnz['ni']   = cheasedata['ni']
         exptnz['nz']   = cheasedata['nz']
         exptnz['Zeff'] = cheasedata['Zeff']

    elif iprofile in [3,'exptnz'] and exptnzflag:
         if   nrhotype[1] in [0,'chease'] and cheaseflag:
              exptnzdata = read_exptnz(exptnzfpath=exptnzpath,setParam=SetParam,chease=cheasepath)
         elif nrhotype[1] in [1,'eqdsk'] and eqdskflag:
              exptnzdata = read_exptnz(exptnzfpath=exptnzpath,setParam=SetParam,eqdsk=eqdskpath)
         exptnz['Ti']   = exptnzdata['Ti']
         exptnz['ni']   = exptnzdata['ni']
         exptnz['nz']   = exptnzdata['nz']
         exptnz['Zeff'] = exptnzdata['Zeff']

    elif iprofile in [4,'profiles'] and profilesflag:
         if   nrhotype[1] in [0,'chease'] and cheaseflag:
              profilesdata = read_profiles(profilesfpath=profilespath,setParam=SetParam,chease=cheasepath)
         elif nrhotype[1] in [1,'eqdsk'] and eqdskflag:
              profilesdata = read_profiles(profilesfpath=profilespath,setParam=SetParam,eqdsk=eqdskpath)
         exptnz['Ti']   = profilesdata['Ti']
         exptnz['ni']   = profilesdata['ni']
         exptnz['nz']   = profilesdata['nz']
         exptnz['Zeff'] = profilesdata['Zeff']

    elif iprofile in [5,'iterdb'] and iterdbflag:
         if   nrhotype[1] in [0,'chease'] and cheaseflag:
              iterdbdata = read_iterdb(iterdbfpath=iterdbpath,setParam=SetParam,chease=cheasepath)
         elif nrhotype[1] in [1,'eqdsk'] and eqdskflag:
              iterdbdata = read_iterdb(iterdbfpath=iterdbpath,setParam=SetParam,eqdsk=eqdskpath)
         exptnz['Ti']   = iterdbdata['Ti']
         exptnz['ni']   = iterdbdata['ni']
         exptnz['nz']   = iterdbdata['nz']
         exptnz['Zeff'] = iterdbdata['Zeff']

    if 'Zi' not in exptnz: exptnz['Zi'] = 1.0
    if 'Zz' not in exptnz: exptnz['Zz'] = 6.0

    if 'nz' not in exptnz:
       exptnz['nz'] = exptnz['Zeff']*exptnz['ne']
       exptnz['nz']-= exptnz['ni']*exptnz['Zi']**2
       exptnz['nz']/= exptnz['Zz']**2

    if   rhopsiflag: rhosize = npy.size(exptnz['rhopsiN'])
    elif rhotorflag: rhosize = npy.size(exptnz['rhotorN'])

    if outfile:
       ofh = open("EXPTNZ",'w')
       if   nrhotype[0] in [0,'rhopsi']:
            rhosize = npy.size(exptnz['rhopsiN'])
            ofh.write('%5d rhopsi,  Te,   ne,   Zeff,   Ti,   ni  profiles\n' % rhosize)
            for i in range(rhosize): ofh.write('%16.6E\n' % exptnz['rhopsiN'][i])
       elif nrhotype[0] in [1,'rhotor']:
            rhosize = npy.size(exptnz['rhotorN'])
            ofh.write('%5d rhotor,  Te,   ne,   Zeff,   Ti,   ni  profiles\n' % rhosize)
            for i in range(rhosize): ofh.write('%16.6E\n' % exptnz['rhotorN'][i])
       for i in range(rhosize): ofh.write('%16.6E\n' % exptnz['Te'][i])
       for i in range(rhosize): ofh.write('%16.6E\n' % exptnz['ne'][i])
       for i in range(rhosize): ofh.write('%16.6E\n' % exptnz['Zeff'][i])
       for i in range(rhosize): ofh.write('%16.6E\n' % exptnz['Ti'][i])
       for i in range(rhosize): ofh.write('%16.6E\n' % exptnz['ni'][i])
       ofh.close()

    return exptnz

def read_profiles(profilesfpath,setParam={},Zeffprofile=True,**kwargs):
    if os.path.isfile(profilesfpath) == False:
       errorFunc = traceback.extract_stack(limit=2)[-2][3]
       errorLine = traceback.extract_stack(limit=2)[-2][1]
       errorFile = traceback.extract_stack(limit=2)[-2][2]
       errMSG    = 'Call %s line %5d in file %s Failed.\n'
       errMSG   += 'Fatal: file %s not found.'
       raise IOError(errMSG %(errorFunc,errorLine,errorFile,profilesfpath))

    rhopsiflag = False; rhotorflag = False
    if 'nrhomesh' in setParam:
        if   setParam['nrhomesh'] in [0,'rhopsi']: rhopsiflag = True
        elif setParam['nrhomesh'] in [1,'rhotor']: rhotorflag = True

    eqdskflag    = False
    cheaseflag   = False
    interpflag   = False
    for key,value in kwargs.items():
        if   key in ['chease','cheasedata','cheasefpath']:
             if   type(value)==str and os.path.isfile(value.strip()):
                                       cheasedata = read_chease(cheasefpath=value.strip())
             elif type(value)==dict:   cheasedata = value.copy()
             if   'rhopsiN' in cheasedata and rhopsiflag:
                  interprho = cheasedata['PSIN'][:]; interpflag = True
             elif 'rhotorN' in cheasedata and rhotorflag:
                  interprho = cheasedata['PHIN'][:]; interpflag = True
             else:                                      interpflag = False
             psi     = cheasedata['PSIN'][:];   phi     = cheasedata['PHIN'][:]
             rhopsiN = cheasedata['rhopsiN'][:];rhotorN = cheasedata['rhotorN'][:]
             cheaseflag = True
        elif key in ['eqdsk','eqdskdata','eqdskfpath']:
             if   type(value)==str and os.path.isfile(value.strip()):
                                       eqdskdata = read_eqdsk(eqdskfpath=value.strip())
             elif type(value)==dict:   eqdskdata = value.copy()
             if   'rhopsi' in eqdskdata and rhopsiflag:
                  interprho = eqdskdata['PSIN'][:]; interpflag = True
             elif 'rhotor' in eqdskdata and rhotorflag:
                  interprho = eqdskdata['PHIN'][:]; interpflag = True
             else:                                    interpflag = False
             psi     = eqdskdata['PSIN'][:];  phi     = eqdskdata['PHIN'][:]
             rhopsiN = eqdskdata['rhopsi'][:];rhotorN = eqdskdata['rhotor'][:]
             eqdskflag = True

    profiles,units = read_profiles_file(profilesfpath.strip())

    PROFILESdata            = {}
    PROFILESdata['rhopsiN'] = npy.sqrt(profiles['psinorm'])
    PROFILESdata['PSIN']    = profiles['psinorm'][:]
    PROFILESdata['Te']      = profiles['te'][:]
    PROFILESdata['ne']      = profiles['ne'][:]
    PROFILESdata['Ti']      = profiles['ti'][:]
    PROFILESdata['ni']      = profiles['ni'][:]
    PROFILESdata['nb']      = profiles['nb'][:]
    PROFILESdata['nz']      = profiles['nz1'][:]

    PROFILESdata['Pb']      = profiles['pb'][:]
    PROFILESdata['VPOL']    = profiles['vpol1'][:]
    PROFILESdata['VTOR']    = profiles['vtor1'][:]

    PROFILESdata['pressure']= profiles['ptot'][:]

    nrhopsi                 = npy.size(PROFILESdata['rhopsiN'])

    PROFILESdata['Zi']      = profiles['z'][1]
    PROFILESdata['Zz']      = profiles['z'][0]

    PROFILESdata['Zeff']    = PROFILESdata['ni']*PROFILESdata['Zi']**2
    PROFILESdata['Zeff']   += PROFILESdata['nz']*PROFILESdata['Zz']**2
    PROFILESdata['Zeff']   /= PROFILESdata['ne']

    if not Zeffprofile:
       Zeff_array = npy.array([1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0])
       Zeff_mean  = npy.mean(PROFILESdata['Zeff'])
       Zeff_diff  = abs(Zeff_array-Zeff_mean)
       Zeff_value = Zeff_array[Zeff_diff==min(Zeff_diff)][0]
       PROFILESdata['Zeff']    = npy.ones(nrhopsi)*Zeff_value

    PROFILESdata['pprime']  = derivative(x=PROFILESdata['rhopsiN']**2,fx=PROFILESdata['pressure'],method='CubicSpline')

    if   interpflag:
         if   rhopsiflag:
              PROFILESdata['Te']      = interp(PROFILESdata['PSIN'],PROFILESdata['Te'],psi)
              PROFILESdata['Ti']      = interp(PROFILESdata['PSIN'],PROFILESdata['Ti'],psi)
              PROFILESdata['ne']      = interp(PROFILESdata['PSIN'],PROFILESdata['ne'],psi)
              PROFILESdata['ni']      = interp(PROFILESdata['PSIN'],PROFILESdata['ni'],psi)
              PROFILESdata['nz']      = interp(PROFILESdata['PSIN'],PROFILESdata['nz'],psi)
              PROFILESdata['Pb']      = interp(PROFILESdata['PSIN'],PROFILESdata['Pb'],psi)
              PROFILESdata['VTOR']    = interp(PROFILESdata['PSIN'],PROFILESdata['VTOR'],psi)
              PROFILESdata['VPOL']    = interp(PROFILESdata['PSIN'],PROFILESdata['VPOL'],psi) 
              PROFILESdata['Zeff']    = interp(PROFILESdata['PSIN'],PROFILESdata['Zeff'],psi)
              PROFILESdata['pprime']  = interp(PROFILESdata['PSIN'],PROFILESdata['pprime'],psi)
              PROFILESdata['pressure']= interp(PROFILESdata['PSIN'],PROFILESdata['pressure'],psi)
              PROFILESdata['rhopsiN'] = rhopsiN[:]
              PROFILESdata['rhotorN'] = rhotorN[:]
         elif rhotorflag:
              PROFILESdata['Te']       = interp(PROFILESdata['PSIN'],PROFILESdata['Te'],psi,phi,phi)
              PROFILESdata['Ti']       = interp(PROFILESdata['PSIN'],PROFILESdata['Ti'],psi,phi,phi)
              PROFILESdata['ne']       = interp(PROFILESdata['PSIN'],PROFILESdata['ne'],psi,phi,phi)
              PROFILESdata['ni']       = interp(PROFILESdata['PSIN'],PROFILESdata['ni'],psi,phi,phi)
              PROFILESdata['nz']       = interp(PROFILESdata['PSIN'],PROFILESdata['nz'],psi,phi,phi)
              PROFILESdata['Pb']       = interp(PROFILESdata['PSIN'],PROFILESdata['Pb'],psi,phi,phi)
              PROFILESdata['VTOR']     = interp(PROFILESdata['PSIN'],PROFILESdata['VTOR'],psi,phi,phi)
              PROFILESdata['VPOL']     = interp(PROFILESdata['PSIN'],PROFILESdata['VPOL'],psi,phi,phi)
              PROFILESdata['Zeff']     = interp(PROFILESdata['PSIN'],PROFILESdata['Zeff'],psi,phi,phi)
              PROFILESdata['pprime']   = interp(PROFILESdata['PSIN'],PROFILESdata['pprime'],psi,phi,phi)
              PROFILESdata['pressure'] = interp(PROFILESdata['PSIN'],PROFILESdata['pressure'],psi,phi,phi)
              PROFILESdata['rhopsiN']  = rhopsiN[:]
              PROFILESdata['rhotorN']  = rhotorN[:]
    elif rhotorflag and not interpflag:
         print("WARNING: setParam['nrhomesh'] = 1 or rhotor, but the path to a target rhotor is not provided.")
         print("         Converting the profiles to toroidal (phi) coordinates could not be done, and")
         print("         all profiles are provided in the poloidal (psi) coordinates.")

    return PROFILESdata

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

def read_iterdb(iterdbfpath,setParam={},**kwargs):
    if os.path.isfile(iterdbfpath) == False:
       errorFunc = traceback.extract_stack(limit=2)[-2][3]
       errorLine = traceback.extract_stack(limit=2)[-2][1]
       errorFile = traceback.extract_stack(limit=2)[-2][2]
       errMSG    = 'Call %s line %5d in file %s Failed.\n'
       errMSG   += 'Fatal: file %s not found.'
       raise IOError(errMSG %(errorFunc,errorLine,errorFile,iterdbfpath))

    rhopsiflag = False; rhotorflag = False
    if   'nrhomesh' in setParam:
         if   setParam['nrhomesh'] in [0,'rhopsi']: rhopsiflag = True
         elif setParam['nrhomesh'] in [1,'rhotor']: rhotorflag = True

    rhotors,profiles,units = read_iterdb_file(iterdbfpath)

    '''
    Normalizing the rhotor vectors before using them to interpolate the physical quantities
    '''
    rhotorNE    = (rhotors['NE']-rhotors['NE'][0])/(rhotors['NE'][-1]-rhotors['NE'][0])
    rhotorTE    = (rhotors['TE']-rhotors['TE'][0])/(rhotors['TE'][-1]-rhotors['TE'][0])
    rhotorNM1   = (rhotors['NM1']-rhotors['NM1'][0])/(rhotors['NM1'][-1]-rhotors['NM1'][0])
    rhotorTI    = (rhotors['TI']-rhotors['TI'][0])/(rhotors['TI'][-1]-rhotors['TI'][0])
    if 'NM2' in profiles:
     rhotorNM2  = (rhotors['NM2']-rhotors['NM2'][0])/(rhotors['NM2'][-1]-rhotors['NM2'][0])
    if 'VROT' in profiles:
     rhotorVROT = (rhotors['VROT']-rhotors['VROT'][0])/(rhotors['VROT'][-1]-rhotors['VROT'][0])

    eqdskflag    = False
    cheaseflag   = False
    interpflag   = False
    for key,value in kwargs.items():
        if   key in ['chease','cheasedata','cheasefpath']:
             if   type(value)==str and os.path.isfile(value.strip()):
                                       cheasedata = read_chease(cheasefpath=value.strip())
             elif type(value)==dict:   cheasedata = value.copy()
             if   'rhopsiN' in cheasedata and rhopsiflag:
                  interprho = cheasedata['rhopsiN'][:]; interpflag = True
             elif 'rhotorN' in cheasedata and rhotorflag:
                  interprho = cheasedata['rhotorN'][:]; interpflag = True
             rhopsi = cheasedata['rhopsiN'][:];rhotor = cheasedata['rhotorN'][:]
             cheaseflag = True
        elif key in ['eqdsk','eqdskdata','eqdskfpath']:
             if   type(value)==str and os.path.isfile(value.strip()):
                                       eqdskdata = read_eqdsk(eqdskfpath=value.strip())
             elif type(value)==dict:   eqdskdata = value.copy()
             if   'rhopsi' in eqdskdata and rhopsiflag:
                  interprho = eqdskdata['rhopsi'][:]; interpflag = True
             elif 'rhotor' in eqdskdata and rhotorflag:
                  interprho = eqdskdata['rhotor'][:]; interpflag = True
             rhopsi = eqdskdata['rhopsi'][:]; rhotor = eqdskdata['rhotor'][:]
             eqdskflag = True
        else:
             interprho  = rhotorNE[:]

    ITERDBdata       = {}
    ITERDBdata['Zi'] = 1.0
    ITERDBdata['Zz'] = 6.0

    if   rhopsiflag and interpflag:
         ITERDBdata['rhopsiN'] = rhopsi
         ITERDBdata['rhotorN'] = rhotor
    elif rhotorflag and interpflag:
         ITERDBdata['rhopsiN'] = rhopsi
         ITERDBdata['rhotorN'] = rhotor
    elif rhopsiflag and not interpflag:
         print("WARNING: setParam['nrhomesh'] = 0 or rhopsi, but the path to a target rhopsi is not provided.")
         print("         Converting the profiles to poloidal (psi) coordinates could not be done, and")
         print("         all profiles are provided in the toroidal (phi) coordinates.")
    else:
         ITERDBdata['rhopsiN'] = interp(rhotor,rhopsi,rhotorNE)
         ITERDBdata['rhotorN'] = rhotorNE

    nrhosize = npy.size(interprho)

    if   rhopsiflag and interpflag:
         ITERDBdata['Te']      = interp(rhotors['TE'],profiles['TE'],rhotor,rhopsi,rhopsi)
         ITERDBdata['Ti']      = interp(rhotors['TI'],profiles['TI'],rhotor,rhopsi,rhopsi)
         ITERDBdata['ne']      = interp(rhotors['NE'],profiles['NE'],rhotor,rhopsi,rhopsi)
         ITERDBdata['ni']      = interp(rhotors['NM1'],profiles['NM1'],rhotor,rhopsi,rhopsi)
         if 'NM2' in profiles: 
            ITERDBdata['nz']   = interp(rhotors['NM2'],profiles['NM2'],rhotor,rhopsi,rhopsi)
         else:
            ITERDBdata['nz']   = npy.zeros(nrhosize)
         if 'VROT' in profiles :
            ITERDBdata['VROT']   = interp(rhotors['VROT'],profiles['VROT'],rhotor,rhopsi,rhopsi)
         else:
            ITERDBdata['VROT'] = npy.zeros(nrhosize)
    elif rhotorflag and interpflag:
         ITERDBdata['ne']      = interp(rhotors['NE'],profiles['NE'],rhotor)
         ITERDBdata['Te']      = interp(rhotors['TE'],profiles['TE'],rhotor)
         ITERDBdata['ni']      = interp(rhotors['NM1'],profiles['NM1'],rhotor)
         ITERDBdata['Ti']      = interp(rhotors['TI'],profiles['TI'],rhotor)
         if 'NM2' in profiles: 
            ITERDBdata['nz']   = interp(rhotors['NM2'],profiles['NM2'],rhotor)
         else:
            ITERDBdata['nz']   = npy.zeros(nrhosize)
         if 'VROT' in profiles :
            ITERDBdata['VROT'] = interp(rhotors['VROT'],profiles['VROT'],rhotor)
         else:
            ITERDBdata['VROT'] = npy.zeros(nrhosize)
    else:
         ITERDBdata['ne']      = interp(rhotors['NE'],profiles['NE'],rhotorNE)
         ITERDBdata['Te']      = interp(rhotors['TE'],profiles['TE'],rhotorNE)
         ITERDBdata['ni']      = interp(rhotors['NM1'],profiles['NM1'],rhotorNE)
         ITERDBdata['Ti']      = interp(rhotors['TI'],profiles['TI'],rhotorNE)
         if 'NM2' in profiles: 
            ITERDBdata['nz']   = interp(rhotors['NM2'],profiles['NM2'],rhotorNE)
         else:
            ITERDBdata['nz']   = npy.zeros(nrhosize)
         if 'VROT' in profiles :
            ITERDBdata['VROT'] = interp(rhotors['VROT'],profiles['VROT'],rhotorNE)
         else:
            ITERDBdata['VROT'] = npy.zeros(nrhosize)

    ITERDBdata['Zeff']    = ITERDBdata['ni']*ITERDBdata['Zi']**2
    ITERDBdata['Zeff']   += ITERDBdata['nz']*ITERDBdata['Zz']**2
    ITERDBdata['Zeff']   /= ITERDBdata['ne']

    ITERDBdata['pressure']   = ITERDBdata['Te']*ITERDBdata['ne']
    ITERDBdata['pressure']  += ITERDBdata['Ti']*ITERDBdata['ni']
    ITERDBdata['pressure']  += ITERDBdata['Ti']*ITERDBdata['nz']
    ITERDBdata['pressure']  *= 1.602e-19

    if   rhopsiflag and interpflag:
         ITERDBdata['pprime']= derivative(x=ITERDBdata['rhopsiN'],fx=ITERDBdata['pressure'],method='CubicSpline')
    else:
         ITERDBdata['pprime']= derivative(x=ITERDBdata['rhotorN'],fx=ITERDBdata['pressure'],method='CubicSpline')

    return ITERDBdata

def plot_chease(OSPATH,reportpath='',skipfigs=1):
    from matplotlib.backends.backend_pdf import PdfPages
    from glob import glob

    if reportpath == '':
       report = False
    else:
       report = True

    if not os.path.isfile(OSPATH):
       srhpath    = os.path.join(OSPATH,'*.h5')
       h5list     = sorted(glob(srhpath))
    else:
       h5list     = [OSPATH]

    if not report:
       srhpath    = os.path.join(OSPATH,'*_EQDSK')
       eqdsklist  = sorted(glob(srhpath))

       profpath   = os.path.join(OSPATH,'*_PROFILES')
       proflist   = sorted(glob(profpath))

       srhpath    = os.path.join(OSPATH,'EXPEQ_iter*.OUT')
       expeqlist  = sorted(glob(srhpath))

       srhpath    = os.path.join(OSPATH,'EXPTNZ_iter*.OUT')
       exptnzlist = sorted(glob(srhpath))

       icounter = 1
       for h5fid in h5list[0::skipfigs+1]:
           print('Plotting CHEASE data in: %s ...' % h5fid)
           if   h5fid[13:17] == 'iter':
                caselabel  = h5fid[13:21]
           elif h5fid[8:14] in ['KEFITD','MEFITD']:
                caselabel  = h5fid[8:21]
           else:
                caselabel  = h5fid

           CHEASEdata = read_chease(cheasefpath=h5fid)
           CHEASEdataKeys = CHEASEdata.keys()

           EXPTNZdata   = read_exptnz(exptnzlist[h5list.index(h5fid)],eqdsk=eqdsklist[0])
           EXPEQdata    = read_expeq(expeqlist[h5list.index(h5fid)])
           EQDSKdata    = read_eqdsk(eqdsklist[0])
           if proflist:
              PROFILESdata = read_profiles(proflist[0])

           EDENfig = plt.figure("Electron Density")
           plt.plot(CHEASEdata['rhopsiN'],CHEASEdata['ne'],linestyle='-',label='CHESAE-'+caselabel[-6:-3])
           plt.plot(EXPTNZdata['rhopsiN'],EXPTNZdata['ne'],linestyle=':',label='EXPTNZ-'+caselabel[-6:-3])
           plt.title('Electron Density Profiles')
           plt.xlabel('$\\rho_{\psi_N}$')
           plt.ylabel('$n_e$')
           plt.legend()
       
           GDNEfig = plt.figure("Electron Density Gradient")
           plt.plot(CHEASEdata['rhopsiN'],CHEASEdata['nePrime'],label=caselabel)
           plt.title('Electron Density Gradient Profiles')
           plt.xlabel('$\\rho_{\psi_N}$')
           plt.ylabel('$\\nabla{n_e}$')
           plt.legend()
       
           ETMPfig = plt.figure("Electron Temperature")
           plt.plot(CHEASEdata['rhopsiN'],CHEASEdata['Te'],linestyle='-',label='CHESAE-'+caselabel[-6:-3])
           plt.plot(EXPTNZdata['rhopsiN'],EXPTNZdata['Te'],linestyle='-',label='EXPTNZ-'+caselabel[-6:-3])
           plt.title('Electron Temperature Profiles')
           plt.xlabel('$\\rho_{\psi_N}$')
           plt.ylabel('$T_e$')
           plt.legend()

           GDTEfig = plt.figure("Electron Temperature Gradient")
           plt.plot(CHEASEdata['rhopsiN'],CHEASEdata['TePrime'],label=caselabel)
           plt.title('Electron Temperature Gradient Profiles')
           plt.xlabel('$\\rho_{\psi_N}$')
           plt.ylabel('$\\nabla{T_e}$')
           plt.legend()
       
           SFTYfig = plt.figure("Safety Factor (q)")
           plt.plot(CHEASEdata['rhopsiN'],CHEASEdata['q'],linestyle='-',label='CHESAE-'+caselabel[-6:-3])
           plt.plot(EQDSKdata['rhopsi'], EQDSKdata['qpsi'], linestyle=':',label='EQDSK')
           if 'q' in EXPEQdata:
              plt.plot(EXPEQdata['rhopsiN'],EXPEQdata['q'],linestyle='--',label='EXPEQ-'+caselabel[-6:-3])
           plt.title("Safety Factor Profiles")
           plt.xlabel('$\\rho_{\psi_N}$')
           plt.ylabel("q")
           plt.legend()

           EQDSKPressure = mu0*EQDSKdata['pressure']/EQDSKdata['BCTR']**2
           EXPTNZPressure = mu0*EXPTNZdata['pressure']/EQDSKdata['BCTR']**2
           if proflist:
              PROFILESPressure = mu0*PROFILESdata['pressure']/EQDSKdata['BCTR']**2
           TPPRfig = plt.figure("Plasma Pressure")
           plt.plot(CHEASEdata['rhopsiN'],CHEASEdata['PN'],  linestyle='solid', label='CHESAE-'+caselabel[-6:-3])
           plt.plot(EQDSKdata['rhopsi'], EQDSKPressure,      linestyle='dotted', label='EQDSK')
           plt.plot(EXPTNZdata['rhopsiN'],EXPTNZPressure,    linestyle='dashdot',label='EXPTNZ-'+caselabel[-6:-3])
           if proflist:
              plt.plot(PROFILESdata['rhopsiN'],PROFILESPressure,linestyle=(0,(5,1)),label='PROFILES')
           if 'pN' in EXPEQdata:
              plt.plot(EXPEQdata['rhopsiN'],EXPEQdata['pN'],linestyle='dashed',label='EXPEQ-'+caselabel[-6:-3])
           plt.title('Plasma Pressure Profiles')
           plt.xlabel('$\\rho_{psi_N}$')
           plt.ylabel('$P$')
           plt.legend()
       
           EQDSKPPrime = mu0*EQDSKdata['pprime']*EQDSKdata['RCTR']**2/EQDSKdata['BCTR']
           EXPTNZPPrime = mu0*EXPTNZdata['pprime']*EQDSKdata['RCTR']**2/EQDSKdata['BCTR']
           if proflist:
              PROFILESPPrime = mu0*PROFILESdata['pprime']*EQDSKdata['RCTR']**2/EQDSKdata['BCTR']
           PPRMfig = plt.figure("P'")
           plt.plot(CHEASEdata['rhopsiN'],CHEASEdata['PPrimeN'], linestyle='solid', label='CHESAE-'+caselabel[-6:-3])
           plt.plot(EQDSKdata['rhopsi'], EQDSKPPrime,            linestyle='dotted', label='EQDSK')
          #plt.plot(EXPTNZdata['rhopsiN'],EXPTNZPPrime,          linestyle='dashdot',label='EXPTNZ-'+caselabel[-6:-3])
           if proflist:
              plt.plot(PROFILESdata['rhopsiN'],PROFILESPPrime,      linestyle=(0,(5,1)),label='PROFILES')
           if 'pprimeN' in EXPEQdata:
              plt.plot(EXPEQdata['rhopsiN'],EXPEQdata['pprimeN'],linestyle='dashed',label='EXPEQ-'+caselabel[-6:-3])
           plt.title("P' Profiles")
           plt.xlabel('$\\rho_{\psi_N}$')
           plt.ylabel("P'")
           plt.legend()
       
           EQDSKTTPrime = EQDSKdata['ffprime']/EQDSKdata['BCTR']
           TTPMfig = plt.figure("TT'")
           plt.plot(CHEASEdata['rhopsiN'],CHEASEdata['TTPrimeN'],linestyle='-',label='CHESAE-'+caselabel[-6:-3])
           plt.plot(EQDSKdata['rhopsi'], EQDSKTTPrime,          linestyle=':',label='EQDSK')
           if 'ttprimeN' in EXPEQdata:
              plt.plot(EXPEQdata['rhopsiN'],EXPEQdata['ttprimeN'],linestyle='-',label='EXPEQ-'+caselabel[-6:-3])
           plt.title("TT' Profiles")
           plt.xlabel('$\\rho_{\psi_N}$')
           plt.ylabel("TT'")
           plt.legend()
       
           ISTRfig = plt.figure("I*")
           plt.plot(CHEASEdata['rhopsiN'],CHEASEdata['ISTRN'],linestyle='-',label='CHESAE-'+caselabel[-6:-3])
           if 'istrN' in EXPEQdata:
              plt.plot(EXPEQdata['rhopsiN'],EXPEQdata['istrN'],linestyle='--',label='EXPEQ-'+caselabel[-6:-3])
           plt.title("I* Profiles")
           plt.xlabel('$\\rho_{\psi_N}$')
           plt.ylabel("I*")
           plt.legend()
       
           ICRTfig = plt.figure("Parallel Current")
           plt.plot(CHEASEdata['rhopsiN'],CHEASEdata['IPRLN'],linestyle='-',label='CHESAE-'+caselabel[-6:-3])
           if 'iprlN' in EXPEQdata:
              plt.plot(EXPEQdata['rhopsiN'],EXPEQdata['iprlN'],linestyle='--',label='EXPEQ-'+caselabel[-6:-3])
           plt.title('Parallel Current Profiles')
           plt.xlabel('$\\rho_{\psi_N}$')
           plt.ylabel('$I_{||}$')
           plt.legend()
       
           JCRTfig = plt.figure("Parallel Current Density")
           plt.plot(CHEASEdata['rhopsiN'],CHEASEdata['JPRLN'],linestyle='-',label='CHESAE-'+caselabel[-6:-3])
           if 'jprlN' in EXPEQdata:
              plt.plot(EXPEQdata['rhopsiN'],EXPEQdata['jprlN'],linestyle='--',label='EXPEQ-'+caselabel[-6:-3])
           plt.title('Parallel Current Density Profiles')
           plt.xlabel('$\\rho_{\psi_N}$')
           plt.ylabel('$J_{||}$')
           plt.legend()
       
           SCRTfig = plt.figure("Bootstrap Currents")
           plt.plot(CHEASEdata['rhopsiN'],CHEASEdata['JBSN'],label='$J_{BS}$-'+caselabel)
           plt.title('Bootstrap Current Density Profiles')
           plt.xlabel('$\\rho_{\psi_N}$')
           plt.ylabel('$J_{BS}$')
           plt.legend()
       
           (CHEASEdata['PSIN2D'],CHEASEdata['CHIN2D']) = npy.meshgrid(CHEASEdata['PSIN'],CHEASEdata['CHIN'])
           BF2Dfig = plt.figure("Magnetic Field, B($\psi$,$\chi$)")
           plt.contour(CHEASEdata['CHIN2D'],CHEASEdata['PSIN2D'],CHEASEdata['BN'])
           plt.title('Magnetic Field Profiles')
           plt.xlabel('$\chi$')
           plt.ylabel('$\psi$')
       
           JPHIfig = plt.figure("Toroidal Current")
           plt.contour(CHEASEdata['CHIN2D'],CHEASEdata['PSIN2D'],CHEASEdata['JPHIN'],cmap=plt.cm.hot)
           plt.title('Toroidal Current Profiles')
           plt.xlabel('$\\rho_{\psi_N}$')
           plt.ylabel('$J_{\phi}$')
       
           BFRZfig = plt.figure("Magnetic Field, B(R,Z}")
           plt.contour(CHEASEdata['RN'],CHEASEdata['ZN'],CHEASEdata['BN'])
           plt.title('Magnetic Field Profiles, B(R,Z}')
           plt.xlabel('$R$')
           plt.ylabel('$Z$')

           PSRZfig = plt.figure("Magnetic Poloidal Boundary Surface")
           plt.plot(CHEASEdata['rbound'],CHEASEdata['zbound'],label=caselabel)
           plt.title('Magnetic Poloidal Boundary Surface')
           plt.xlabel('$R$')
           plt.ylabel('$Z$')
           plt.legend()

           del(CHEASEdata)

           chsfigs = PdfPages('cheaseresults.pdf')
           chsfigs.savefig(EDENfig)
           chsfigs.savefig(GDNEfig)
           chsfigs.savefig(ETMPfig)
           chsfigs.savefig(GDTEfig)
           chsfigs.savefig(SFTYfig)
           chsfigs.savefig(TPPRfig)
           chsfigs.savefig(PPRMfig)
           chsfigs.savefig(TTPMfig)
           chsfigs.savefig(ISTRfig)
           chsfigs.savefig(ICRTfig)
           chsfigs.savefig(JPHIfig)
           chsfigs.savefig(JCRTfig)
           chsfigs.savefig(SCRTfig)
           chsfigs.savefig(BF2Dfig)
           chsfigs.savefig(BFRZfig)
           chsfigs.savefig(PSRZfig)
           chsfigs.close()

    elif report:
       for h5fid in h5list[0::skipfigs+1]:
           print('Plotting CHEASE data in: %s ...' % h5fid)

           if reportpath[-1] != "/": reportpath += "/"
           if 'report' not in reportpath:
              reportpath += "report/"
           if not os.path.isdir(reportpath):
              os.system('mkdir %s' % reportpath)

           CHEASEdata = read_chease(cheasefpath=h5fid)
           CHEASEdataKeys = CHEASEdata.keys()

           ETMPfig = plt.figure("Electron Temperature")
           plt.plot(CHEASEdata['PSIN'],CHEASEdata['Te'],label='Te')
           plt.plot(CHEASEdata['PSIN'],CHEASEdata['Ti'],label='Ti')
           plt.title('Temperature Profiles')
           plt.xlabel('$\psi$')
           plt.ylabel('$T_e$')
           plt.legend()
           ETMPfig.savefig(reportpath+"chease_temperature.png")
           plt.close(ETMPfig)

           EDENfig = plt.figure("Electron Density")
           plt.plot(CHEASEdata['PSIN'],CHEASEdata['ne'],label='ne')
           plt.plot(CHEASEdata['PSIN'],CHEASEdata['ni'],label='ni')
           plt.title('Density Profiles')
           plt.xlabel('$\psi$')
           plt.ylabel('$n_e$')
           plt.legend()
           EDENfig.savefig(reportpath+"chease_density.png")
           plt.close(EDENfig)
           
           PPRMfig = plt.figure("P'")
           plt.plot(CHEASEdata['PSIN'],CHEASEdata['PPrimeN'])
           plt.title("P' Profiles")
           plt.xlabel('$\psi$')
           plt.ylabel("P'")
           PPRMfig.savefig(reportpath+"chease_pprime.png")
           plt.close(PPRMfig)
       
           TTPMfig = plt.figure("TT'")
           plt.plot(CHEASEdata['PSIN'],CHEASEdata['TTPrimeN'])
           plt.title("TT' Profiles")
           plt.xlabel('$\psi$')
           plt.ylabel("TT'")
           TTPMfig.savefig(reportpath+"chease_ttprime.png")
           plt.close(TTPMfig)
       
           SFTYfig = plt.figure("Safety Factor (q)")
           plt.plot(CHEASEdata['PSIN'],CHEASEdata['q'])
           plt.title("Safety Factor Profiles")
           plt.xlabel('$\psi$')
           plt.ylabel("q")
           SFTYfig.savefig(reportpath+"chease_safetyfactor.png")
           plt.close(SFTYfig)

           PSRZfig = plt.figure("Magnetic Poloidal Boundary Surface")
           plt.plot(CHEASEdata['rbound'],CHEASEdata['zbound'])
           plt.title('Magnetic Poloidal Boundary Surface')
           plt.xlabel('$R$')
           plt.ylabel('$Z$')
           PSRZfig.savefig(reportpath+"chease_magsurfbound.png")
           plt.close(PSRZfig)

           del(CHEASEdata)

    return 1

def cheasepy(srcVals={},namelistVals={},pltVals={},cheaseVals={}):
    CGRN = '\x1b[32m'
    CYLW = '\x1b[33m'
    CBLU = '\x1b[34m'
    CRED = '\x1b[91m'
    CEND = '\x1b[0m'

    while True:
          print(CGRN+'Select on of the following options:'+CEND)
          print(CGRN+'(1) Run the Code and Plot the Outputs.'+CEND)
          print(CGRN+'(2) Plot the Available Outputs.'+CEND)
          print(CGRN+'(3) Remove Input/Output Files.'+CEND)
          print(CGRN+'(4) Exit.'+CEND)
          try:
              if 'runmode' in cheaseVals:
                 selection = int(cheaseVals['runmode'])
                 print('Selected Option: %d' %selection)
              else:
                 selection = int(input('Selected Option: '))
              if    selection in [1,2,3,4]:
                    if    selection == 3:
                       if glob('./NGA'):                   os.system('rm NGA')
                       if glob('./NDES'):                  os.system('rm NDES')
                       if glob('./*OUT*'):                 os.system('rm *OUT*')
                       if glob('./*.pdf'):                 os.system('rm *.pdf')
                       if glob('./EXPEQ'):                 os.system('rm EXPEQ')
                       if glob('./EXPTNZ'):                os.system('rm EXPTNZ')
                       if glob('./*_EQDSK'):               os.system('rm *_EQDSK')
                       if glob('./*_EXPTNZ'):              os.system('rm *_EXPTNZ')
                       if glob('./*CHEASE'):               os.system('rm *CHEASE')
                       if glob('./*ITERDB'):               os.system('rm *ITERDB')
                       if glob('./*PROFILES'):             os.system('rm *PROFILES')
                       if glob('./ogyropsi*'):             os.system('rm ogyropsi*')
                       if glob('./EXPEQ_EQDSK*'):          os.system('rm EXPEQ_EQDSK*')
                       if glob('./EXPEQ_EXPEQ*'):          os.system('rm EXPEQ_EXPEQ*')
                       if glob('./chease_namelist*'):      os.system('rm chease_namelist*')
                       if glob('./chease_parameters.csv'): os.system('rm chease_parameters.csv')
                       continue
                    elif  selection == 4: sys.exit()
                    else: break

              else: raise(NameError)
          except NameError:
             print(CRED+'Select between 1,2,3, and 4 options.'+CEND)
             continue

    '''
    "iterTotal     = n" where n is the number of chease iterations
                              n = 0, 1, 2, ...
                              n = 0 just reproduce the equilibrium
                                    given in the input sources based
                                    on the input set of parameters.
    "boundary_type  = n" where n in geometry options
    Boundary Surface options:
    n = 0 or 'asis' --- Take the boundary surface 'as is'
    n = 1 or 'interp' - Interpolate the given boundary surface

    "current_src   = n" where n in geometry options
    Geometry Source options:
    n = 0 or 'chease'
    n = 1 or 'eqdsk'
    n = 2 or 'expeq'

    "rhomesh_src   = n" where n in meshgrid options
    Meshgrid Source options:
    n = 0 or 'chease'
    n = 1 or 'eqdsk'

    "pressure_src  = n" where n in profiles options
    "eprofiles_src = n" where n in profiles options
    "iprofiles_src = n" where n in profiles options
    Profiles Source options:
    n = 0 or 'chease'
    n = 1 or 'eqdsk'
    n = 2 or 'expeq'
    n = 3 or 'exptnz'
    n = 4 or 'profiles'
    n = 5 or 'iterdb'
    '''

    srcValsKeys     = srcVals.keys()
    if   'iterTotal'   in srcVals: iterTotal     = srcVals['iterTotal']
    else:                          iterTotal     = 0  

    if 'boundary_type' in srcVals: boundary_type = srcVals['boundary_type']   
    else:                          boundary_type = 0

    rhomesh_update = False
    if   'rhomesh_src' in srcVals:
         if   type(srcVals['rhomesh_src']) == list:
              rhomesh_src    = srcVals['rhomesh_src'][0]
              rhomesh_update = True
         else:
              rhomesh_src    = srcVals['rhomesh_src']
    else:
              rhomesh_src    = 1

    current_update = False
    if   'current_src' in srcVals:
         if   type(srcVals['current_src']) == list:
              current_src    = srcVals['current_src'][0]
              current_update = True
         else:
              current_src    = srcVals['current_src']
    else:
              current_src    = 1

    pressure_update = False
    if   'pressure_src' in srcVals:
         if   type(srcVals['pressure_src']) == list:
              pressure_src    = srcVals['pressure_src'][0]
              pressure_update = True
         else:
              pressure_src    = srcVals['pressure_src']
    else:
              pressure_src    = 1

    eprofiles_update = False
    if   'eprofiles_src' in srcVals:
         if   type(srcVals['eprofiles_src']) == list:
              eprofiles_src    = srcVals['eprofiles_src'][0]
              eprofiles_update = True
         else:
              eprofiles_src    = srcVals['eprofiles_src']
    else:
              eprofiles_src    = 3

    iprofiles_update = False
    if   'iprofiles_src' in srcVals:
         if   type(srcVals['iprofiles_src']) == list:
              iprofiles_src    = srcVals['iprofiles_src'][0]
              iprofiles_update = True
         else:
              iprofiles_src    = srcVals['iprofiles_src']
    else:
              iprofiles_src    = 3

    eqdskrequired    = False
    expeqrequired    = False
    cheaserequired   = False
    exptnzrequired   = False
    iterdbrequired   = False
    profilesrequired = False

    if   rhomesh_src  in [0,'chease']:   cheaserequired   = True
    elif rhomesh_src  in [1,'eqdsk']:    eqdskrequired    = True

    if   current_src  in [0,'chease']:   cheaserequired   = True
    elif current_src  in [1,'eqdsk']:    eqdskrequired    = True
    elif current_src  in [2,'expeq']:    expeqrequired    = True

    if not (cheaserequired or eqdskrequired):
       raise IOError('FATAL: No Geometry File is provided!')

    if   pressure_src in [0,'chease']:   cheaserequired   = True
    elif pressure_src in [1,'eqdsk']:    eqdskrequired    = True
    elif pressure_src in [2,'expeq']:    expeqrequired    = True
    elif pressure_src in [3,'exptnz']:   exptnzrequired   = True
    elif pressure_src in [4,'profiles']: profilesrequired = True
    elif pressure_src in [5,'iterdb']:   iterdbrequired   = True

    if   eprofiles_src in [0,'chease']:   cheaserequired   = True
    elif eprofiles_src in [3,'exptnz']:   exptnzrequired   = True
    elif eprofiles_src in [4,'profiles']: profilesrequired = True
    elif eprofiles_src in [5,'iterdb']:   iterdbrequired   = True

    if   iprofiles_src in [0,'chease']:   cheaserequired   = True
    elif iprofiles_src in [3,'exptnz']:   exptnzrequired   = True
    elif iprofiles_src in [4,'profiles']: profilesrequired = True
    elif iprofiles_src in [5,'iterdb']:   iterdbrequired   = True

    if not (exptnzrequired or profilesrequired or iterdbrequired):
       raise IOError('FATAL: No Profiles File is provided!')

    if selection == 1:
       if glob('./EQSDK')!=[]:     EQDSKexist  = True
       else:                       EQDSKexist  = False

       if glob('./EXPEQ')!=[]:     EXPEQexist  = True
       else:                       EXPEQexist  = False

       if glob('./EXPTNZ')!=[]:    EXPTNZexist = True
       else:                       EXPTNZexist = False

       if glob('./*CHEASE')!=[]:   CHEASEexist = True
       else:                       CHEASEexist = False

       if glob('./*ITERDB')!=[]:   ITERDBexist = True
       else:                       ITERDBexist = False

       if glob('./*PROFILES')!=[]: PROFILESexist = True
       else:                       PROFILESexist = False

       if 'removeinputs' in cheaseVals:
          if type(cheaseVals['removeinputs'])==str:
             removeinputs = cheaseVals['removeinputs'].lower()
          else:
             removeinputs = cheaseVals['removeinputs']
       else:
          if EXPEQexist or EXPTNZexist or EQDSKexist or ITERDBexist or PROFILESexist:
             os.system('ls')
             if   sys.version_info[0] > 3:
                removeinputs = str(input(CRED+'Remove input files from a previous run (yes/no)? '+CEND)).lower()
             elif sys.version_info[0] < 3:
                removeinputs = raw_input(CRED+'Remove input files from a previous run (yes/no)? '+CEND).lower()
          else:
             removeinputs = True
 
       while True:
             print(CYLW+'Select CHEASE running mode:'+CEND)
             print(CYLW+'(1) Check Equilibrium Preservation Over Multiple Iterations.'+CEND)
             print(CYLW+'(2) Converge to Total Current by correcting Current.'+CEND)
             print(CYLW+'(3) Converge to Total Current by correcting Pressure.'+CEND)
             try:
                if 'cheasemode' in cheaseVals:
                   cheasemode = cheaseVals['cheasemode']
                   print('Selected Option: %d' %cheasemode)
                else:
                   cheasemode = int(input('Selected Option: '))
                if    cheasemode in [1,2,3]: break
                else: raise(NameError)
             except NameError:
                print(CRED+'Select between 1, 2, and 3 options.'+CEND)
                continue
    
       namelistParam           = {}
       if 'removeoutputs' in cheaseVals:
          if type(cheaseVals['removeoutputs'])==str:
             removeoutputs = cheaseVals['removeoutputs'].lower()
          else:
             removeoutputs = cheaseVals['removeoutputs']
          print(CBLU+'Remove output files of previous a run (yes/no)? '+str(removeoutputs)+CEND)
       elif removeinputs in ['yes','y',1,True]:
          removeoutputs = True
       else:
          os.system('ls')
          if   sys.version_info[0] >=3:
             removeoutputs = str(input(CBLU+'Remove output files from previous a run (yes/no)? '+CEND)).lower()
          elif sys.version_info[0] < 3:
             removeoutputs = raw_input(CBLU+'Remove output files from previous a run (yes/no)? '+CEND).lower()

       if removeoutputs in ['yes','y',1,True]:
          if glob('./NGA'):                   os.system('rm NGA')
          if glob('./NDES'):                  os.system('rm NDES')
          if glob('./*OUT*'):                 os.system('rm *OUT*')
          if glob('./*.pdf'):                 os.system('rm *.pdf')
          if glob('./EXPEQ'):                 os.system('rm EXPEQ')
          if glob('./EXPTNZ'):                os.system('rm EXPTNZ')
          if glob('./*CHEASE'):               os.system('rm *CHEASE')
          if glob('./*ITERDB'):               os.system('rm *ITERDB')
          if glob('./*PROFILES'):             os.system('rm *PROFILES')
          if glob('./ogyropsi*'):             os.system('rm ogyropsi*')
          if glob('./EXPEQ_EQDSK*'):          os.system('rm EXPEQ_EQDSK*')
          if glob('./EXPEQ_EXPEQ*'):          os.system('rm EXPEQ_EXPEQ*')
          if glob('./chease_namelist*'):      os.system('rm chease_namelist*')
          print(CRED+'List of Available CHEASE Files:'+CEND)
          os.system('ls')
    
       if len(glob('./*_EQDSK'))==0 or len(glob('./chease_parameters.csv'))==0:
          shotlist = sorted(glob('./shots/*'))
          if 'shotpath' in cheaseVals:
             shotrec = cheaseVals['shotpath']
             if os.path.isdir(shotrec):
                 shotpath = shotrec
                 if shotpath[-1]!='/': shotpath+='/'
                 slashpos = findall(shotpath,'/')
                 if len(slashpos)>=2:
                     shotfile = shotpath[slashpos[-2]+1:-1]
                 else:
                     shotfile = shotpath[0:]
             elif type(shotrec)==int:
                 print('Select Shot Number: %d' %shotrec)
                 shotpath = shotlist[shotrec-1]
                 shotfile = shotlist[shotrec-1][8:]
             else:
                 raise IOError(shotrec+' IS NOT A DIRECTORY OR SHOT NUMBER. Exit!')
          else:
             while True:
                print(CYLW+'List of the available shots:'+CEND)
                for ishot in range(len(shotlist)):
                    print(CYLW+'(%02d) %s' % (ishot+1,shotlist[ishot][8:])+CEND)
                try:
                   shotrec = int(input('Select Shot Number: '))
                   if shotrec-1 in range(len(shotlist)):
                      print(CGRN+'Chease runs the %s shot.' % shotlist[shotrec-1][8:]+CEND)
                      shotpath = shotlist[shotrec-1]
                      shotfile = shotlist[shotrec-1][8:]
                      break
                   else:
                      raise(NameError)
                except NameError:
                   print(CRED+'Choose ONLY from the %0d available shots.' % len(shotlist)+CEND)
                   continue
          if os.path.isfile('%s/chease_parameters.csv' % (shotpath)):
             os.system('cp %s/chease_parameters.csv .' % shotpath)
             namelist = namelistcreate('chease_parameters.csv',0,namelistParam)
          else:
             raise IOError('chease namelist file NOT FOUND in the given path!')

          if   'NFUNRHO'  in namelist: rhomesh_type = int(namelist['NFUNRHO'][0])
          elif 'NRHOMESH' in namelist: rhomesh_type = int(namelist['NRHOMESH'][0])
          else:                        rhomesh_type = 0

          if   'NPPFUN'   in namelist: pressure_type = int(namelist['NPPFUN'][0])
          else:                        pressure_type = 4

          if   'NSTTP'    in namelist: current_type  = int(namelist['NSTTP'][0])
          else:                        current_type  = 1

          if   eqdskrequired and os.path.isfile('%s/%s_EQDSK' % (shotpath,shotfile)):
               os.system('cp   %s/%s_EQDSK .'          % (shotpath,shotfile))
               eqdskfpath =      '%s_EQDSK'            %          (shotfile)
          elif eqdskrequired:
               raise IOError('EQDSK file NOT FOUND in the given path!')

          if   cheaserequired and os.path.isfile('%s/%s_CHEASE' % (shotpath,shotfile)):
               os.system('cp   %s/%s_CHEASE .'         % (shotpath,shotfile))
               cheasefpath =     '%s_CHEASE'           %          (shotfile)
          elif cheaserequired:
               raise IOError('CHEASE file NOT FOUND in the given path!')

          if   expeqrequired and os.path.isfile(  '%s/%s_EXPEQ' % (shotpath,shotfile)):
               os.system('cp     %s/%s_EXPEQ ./EXPEQ'    % (shotpath,shotfile))
               expeqfpath = 'EXPEQ'
          elif expeqrequired:
               raise IOError('EXPEQ file NOT FOUND in the given path!')

          exptnzParam              = {}
          exptnzParam['nrhomesh']  = [rhomesh_type,rhomesh_src]
          exptnzParam['eprofiles'] = eprofiles_src
          exptnzParam['iprofiles'] = iprofiles_src

          if int(namelist['NBSEXPQ'][0]) != 0:
             if   exptnzrequired   and os.path.isfile('%s/%s_EXPTNZ'  %(shotpath,shotfile)):
                  os.system('cp   %s/%s_EXPTNZ .'    % (shotpath,shotfile))
                  exptnzfpath =     '%s_EXPTNZ'      %          (shotfile)
                  if   rhomesh_src in [0,'chease']:
                       write_exptnz(setParam=exptnzParam,exptnz=exptnzfpath,chease=cheasefpath)
                  elif rhomesh_src in [1,'eqdsk']:
                       write_exptnz(setParam=exptnzParam,exptnz=exptnzfpath,eqdsk=eqdskfpath)
             elif exptnzrequired:
                  raise IOError('EXPTNZ file NOT FOUND in the given path!')

             if   profilesrequired and os.path.isfile('%s/%s_PROFILES'%(shotpath,shotfile)):
                  os.system('cp   %s/%s_PROFILES .'  % (shotpath,shotfile))
                  profilesfpath =   '%s_PROFILES'    %          (shotfile)
                  if   rhomesh_src in [0,'chease']:
                       write_exptnz(setParam=exptnzParam,profiles=profilesfpath,chease=cheasefpath)
                  elif rhomesh_src in [1,'eqdsk']:
                       write_exptnz(setParam=exptnzParam,profiles=profilesfpath,eqdsk=eqdskfpath)
             elif profilesrequired:
                  raise IOError('Profiles file NOT FOUND in the given path!')

             if   iterdbrequired   and os.path.isfile('%s/%s_ITERDB'  %(shotpath,shotfile)):
                  os.system('cp   %s/%s_ITERDB .'    % (shotpath,shotfile))
                  iterdbfpath =     '%s_ITERDB'      %          (shotfile)
                  if   rhomesh_src in [0,'chease']:
                       write_exptnz(setParam=exptnzParam,iterdb=iterdbfpath,chease=cheasefpath)
                  elif rhomesh_src in [1,'eqdsk']:
                       write_exptnz(setParam=exptnzParam,iterdb=iterdbfpath,eqdsk=eqdskfpath)
             elif iterdbrequired:
                  raise IOError('ITERDB file NOT FOUND in the given path!')

       elif removeinputs in ['yes','y',1,True]:
          if glob('./*_EQDSK'):               os.system('rm *_EQDSK')
          if glob('./*_EXPTNZ'):              os.system('rm *_EXPTNZ')
          if glob('./*_PROFILES'):            os.system('rm *_PROFILES')
          if glob('./chease_parameters.csv'): os.system('rm chease_parameters.csv')
          shotlist = sorted(glob('./shots/*'))
          if 'shotpath' in cheaseVals:
             shotrec = cheaseVals['shotpath']
             if os.path.isdir(shotrec):
                 shotpath = shotrec
                 if shotpath[-1]!='/': shotpath+='/'
                 slashpos = findall(shotpath,'/')
                 if len(slashpos)>=2:
                     shotfile = shotpath[slashpos[-2]+1:-1]
                 else:
                     shotfile = shotpath[0:]
             elif type(shotrec)==int:
                 print('Select Shot Number: %d' %shotrec)
                 shotpath = shotlist[shotrec-1]
                 shotfile = shotlist[shotrec-1][8:]
             else:
                 raise IOError(shotrec+' IS NOT A DIRECTORY OR SHOT NUMBER. Exit!')
          else:
             while True:
                print(CYLW+'List of the available shots:'+CEND)
                for ishot in range(len(shotlist)):
                    print(CYLW+'(%02d) %s' % (ishot+1,shotlist[ishot][8:])+CEND)
                try:
                      shotrec = int(input('Select Shot Number: '))
                      if shotrec-1 in range(len(shotlist)):
                         print(CGRN+'Chease runs the %s shot.' % shotlist[shotrec-1][8:]+CEND)
                         shotpath = shotlist[shotrec-1]
                         shotfile = shotlist[shotrec-1][8:]
                         break
                      else:
                         raise(NameError)
                except NameError:
                   print(CRED+'Choose ONLY from the %0d available shots.' % len(shotlist)+CEND)
                   continue
          if os.path.isfile('%s/chease_parameters.csv' % (shotpath)):
             os.system('cp %s/chease_parameters.csv .' % shotpath)
             namelist = namelistcreate('chease_parameters.csv',0,namelistParam)
          else:
             raise IOError('chease namelist file NOT FOUND in the given path!')

          if   'NFUNRHO'  in namelist: rhomesh_type = int(namelist['NFUNRHO'][0])
          elif 'NRHOMESH' in namelist: rhomesh_type = int(namelist['NRHOMESH'][0])
          else:                        rhomesh_type = 0

          if   'NPPFUN'   in namelist: pressure_type = int(namelist['NPPFUN'][0])
          else:                        pressure_type = 4

          if   'NSTTP'    in namelist: current_type  = int(namelist['NSTTP'][0])
          else:                        current_type  = 1

          if   eqdskrequired and os.path.isfile('%s/%s_EQDSK' % (shotpath,shotfile)):
               os.system('cp   %s/%s_EQDSK .'          % (shotpath,shotfile))
               eqdskfpath =      '%s_EQDSK'            %          (shotfile)
          elif eqdskrequired:
             raise IOError('EQDSK file NOT FOUND in the given path!')

          if   cheaserequired and os.path.isfile('%s/%s_CHEASE' % (shotpath,shotfile)):
               os.system('cp   %s/%s_CHEASE .'         % (shotpath,shotfile))
               cheasefpath =     '%s_CHEASE'           %          (shotfile)
          elif cheaserequired:
               raise IOError('CHEASE file NOT FOUND in the given path!')

          if   expeqrequired and os.path.isfile('%s/%s_EXPEQ' % (shotpath,shotfile)):
               os.system('cp     %s/%s_EXPEQ ./EXPEQ'    % (shotpath,shotfile))
               expeqfpath = 'EXPEQ'
          elif expeqrequired:
               raise IOError('EXPEQ file NOT FOUND in the given path!')

          exptnzParam              = {}
          exptnzParam['nrhomesh']  = [rhomesh_type,rhomesh_src]
          exptnzParam['eprofiles'] = eprofiles_src
          exptnzParam['iprofiles'] = iprofiles_src

          if int(namelist['NBSEXPQ'][0]) != 0:
             if   exptnzrequired   and os.path.isfile('%s/%s_EXPTNZ'  %(shotpath,shotfile)):
                  os.system('cp   %s/%s_EXPTNZ .'    % (shotpath,shotfile))
                  exptnzfpath =     '%s_EXPTNZ'      %          (shotfile)
                  if   rhomesh_src in [0,'chease']:
                       write_exptnz(setParam=exptnzParam,exptnz=exptnzfpath,chease=cheasefpath)
                  elif rhomesh_src in [1,'eqdsk']:
                       write_exptnz(setParam=exptnzParam,exptnz=exptnzfpath,eqdsk=eqdskfpath)
             elif exptnzrequired:
                  raise IOError('EXPTNZ file NOT FOUND in the given path!')

             if   profilesrequired and os.path.isfile('%s/%s_PROFILES' %(shotpath,shotfile)):
                  os.system('cp   %s/%s_PROFILES .'  % (shotpath,shotfile))
                  profilesfpath =   '%s_PROFILES'    %          (shotfile)
                  if   rhomesh_src in [0,'chease']:
                       write_exptnz(setParam=exptnzParam,profiles=profilesfpath,chease=cheasefpath)
                  elif rhomesh_src in [1,'eqdsk']:
                       write_exptnz(setParam=exptnzParam,profiles=profilesfpath,eqdsk=eqdskfpath)
             elif profilesrequired:
                  raise IOError('Profiles file NOT FOUND in the given path!')

             if   iterdbrequired   and os.path.isfile('%s/%s_ITERDB'   %(shotpath,shotfile)):
                  os.system('cp   %s/%s_ITERDB .'    % (shotpath,shotfile))
                  iterdbfpath =     '%s_ITERDB'      %          (shotfile)
                  if   rhomesh_src in [0,'chease']:
                       write_exptnz(setParam=exptnzParam,iterdb=iterdbfpath,chease=cheasefpath)
                  elif rhomesh_src in [1,'eqdsk']:
                       write_exptnz(setParam=exptnzParam,iterdb=iterdbfpath,eqdsk=eqdskfpath)
             elif iterdbrequired:
                  raise IOError('ITERDB file NOT FOUND in the given path!')
       else:
          namelist = namelistcreate('chease_parameters.csv',0,namelistParam)
          if   'NFUNRHO'  in namelist: rhomesh_type = int(namelist['NFUNRHO'][0])
          elif 'NRHOMESH' in namelist: rhomesh_type = int(namelist['NRHOMESH'][0])
          else:                        rhomesh_type = 0

          if   'NPPFUN'   in namelist: pressure_type = int(namelist['NPPFUN'][0])
          else:                        pressure_type = 4

          if   'NSTTP'    in namelist: current_type  = int(namelist['NSTTP'][0])
          else:                        current_type  = 1

          if   cheaserequired   and os.path.isfile('%s/%s_CHEASE'      %(shotpath,shotfile)):
               os.system('cp   %s/%s_CHEASE .'    % (shotpath,shotfile))
               cheasefpath =     '%s_CHEASE'      %          (shotfile)
          elif cheaserequired:
               raise IOError('CHEASE file NOT FOUND in the given path!')

          if   eqdskrequired    and os.path.isfile('%s/%s_EQDSK'   %(shotpath,shotfile)):
               os.system('cp   %s/%s_EQDSK .'     % (shotpath,shotfile))
               eqdskfpath =      '%s_EQDSK'       %          (shotfile)
          elif eqdskrequired:
               raise IOError('EQDSK file NOT FOUND in the given path!')

          if   expeqrequired    and os.path.isfile('%s/%s_EXPEQ'   %(shotpath,shotfile)):
               os.system('cp   %s/%s_EXPEQ .'     % (shotpath,shotfile))
               expeqfpath =      '%s_EXPEQ'       %          (shotfile)
          elif expeqrequired:
               raise IOError('EXPEQ file NOT FOUND in the given path!')

          if   exptnzrequired   and os.path.isfile('%s/%s_EXPTNZ'   %(shotpath,shotfile)):
               os.system('cp   %s/%s_EXPTNZ .'    % (shotpath,shotfile))
               exptnzfpath =     '%s_EXPTNZ'      %          (shotfile)
          elif exptnzrequired:
               raise IOError('EXPTNZ file NOT FOUND in the given path!')

          if   profilesrequired and os.path.isfile('%s/%s_PROFILES' %(shotpath,shotfile)):
               os.system('cp   %s/%s_PROFILES .'  % (shotpath,shotfile))
               profilesfpath =   '%s_PROFILES'    %          (shotfile)
          elif profilesrequired:
               raise IOError('Profiles file NOT FOUND in the given path!')

          if   iterdbrequired   and os.path.isfile('%s/%s_ITERDB'   %(shotpath,shotfile)):
               os.system('cp   %s/%s_ITERDB .'    % (shotpath,shotfile))
               iterdbfpath =     '%s_ITERDB'      %          (shotfile)
          elif iterdbrequired:
               raise IOError('ITERDB file NOT FOUND in the given path!')

       if cheaseVals!={}:
          runchease = 'yes'
          print(CBLU+'Do you want to continue? (yes/no)? '+str(runchease)+CEND)
       else:
          os.system('ls')
          if   sys.version_info[0] >=3:
             runchease = str(input(CBLU+'Do you want to continue? (yes/no)? '+CEND)).lower()
          elif sys.version_info[0] < 3:
             runchease = raw_input(CBLU+'Do you want to continue? (yes/no)? '+CEND).lower()
       if runchease not in ['yes','y',1,True]: sys.exit()

       eqdskdata = read_eqdsk(eqdskfpath=eqdskfpath)
       R0EXP = abs(eqdskdata['RCTR'])
       B0EXP = abs(eqdskdata['BCTR'])
       ITEXP = abs(eqdskdata['CURNT'])

       expeqParam               = {}
       expeqParam['nrhomesh']   = [rhomesh_type,rhomesh_src]
       expeqParam['nppfun']     = [pressure_type,pressure_src]
       expeqParam['nsttp']      = [current_type,current_src]
       expeqParam['boundary']   =  boundary_type
       expeqParam['cheasemode'] =  1
       expeqParam['ITEXP']      =  ITEXP
       expeqParam['R0EXP']      =  R0EXP
       expeqParam['B0EXP']      =  B0EXP

       if   int(namelist['NEQDSK'][0]) == 1:
            print('Reading from EQDSK file.')
            os.system('cp *_EQDSK  EXPEQ')
       elif int(namelist['NEQDSK'][0]) == 0:
            print('Reading from EXPEQ file.')
            if not os.path.isfile('./EXPEQ'):
               if   pressure_src in [1,'eqdsk']:
                    expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath)
               elif pressure_src in [3,'exptnz']:
                    expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,exptnz=exptnzfpath)
               elif pressure_src in [4,'profiles']:
                    expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,profiles=profilesfpath)
               elif pressure_src in [5,'iterdb']:
                    expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,iterdb=iterdbfpath)

       namelistParam['R0EXP'] = R0EXP
       namelistParam['B0EXP'] = B0EXP

       namelistParam['CSSPEC'] = 0.0
       if namelist['NSTTP'] == 5: 
           if   current_src in [0,'chease']:
                if   rhomesh_src in [0,'chease']:
                     cheasedata = read_chease(cheasefpath=cheasefpath)
                elif rhomesh_src in [1,'eqdsk']:
                     cheasedata = read_chease(cheasefpath=cheasefpath,eqdsk=eqdskfpath)
                namelistParam['QSPEC']  = cheasedata['q'][0]
           elif current_src in [1,'eqdsk']:
                if   rhomesh_src in [0,'chease']:
                     eqdskdata = read_eqdsk(eqdskfpath=eqdskfpath)
                elif rhomesh_src in [1,'eqdsk']:
                     eqdskdata = read_eqdsk(eqdskfpath=eqdskfpath,chease=cheasefpath)
                namelistParam['QSPEC']  = eqdskdata['qpsi'][0]
           elif current_src in [2,'expeq']:
                if   rhomesh_src in [0,'chease']:
                     expeqdata = read_expeq(expeqfpath=expeqfpath,chease=cheasefpath)
                elif rhomesh_src in [1,'eqdsk']:
                     expeqdata = read_expeq(expeqfpath=expeqfpath,eqdsk=eqdskfpath)
                namelistParam['QSPEC']  = expeqdata['q'][0]
       namelistParam['NCSCAL'] = 4
       namelist = namelistcreate('chease_parameters.csv',0,namelistParam)
       
       cheasefname = sorted(glob('./ogyropsi_*.h5'))
       if   len(cheasefname)==0:
            it=0
            exit_status = os.system('./chease_hdf5 chease_namelist > iter%03d.OUT' % it)
           #exit_status = os.system('./chease_hdf5 chease_namelist')
           #exit_status = subprocess.call(['./chease_hdf5','chease_namelist'])
            if abs(exit_status) > 0: sys.exit()
            if os.path.isfile('./chease_namelist'): os.system('mv ./chease_namelist ./chease_namelist_iter%03d' % it)
            if os.path.isfile('./ogyropsi.dat'): os.system('mv ./ogyropsi.dat ogyropsi_iter%03d.dat' % it)
            if os.path.isfile('./ogyropsi.h5'): os.system('mv ./ogyropsi.h5 ogyropsi_iter%03d.h5' % it)
            if os.path.isfile('./EXPEQ.OUT'): os.system('mv ./EXPEQ.OUT EXPEQ_iter%03d.OUT' % it)
            if os.path.isfile('./EXPTNZ.OUT'): os.system('mv ./EXPTNZ.OUT EXPTNZ_iter%03d.OUT' % it)
            if os.path.isfile('./EXPEQ_EXPEQ.IN'): os.system('mv ./EXPEQ_EXPEQ.IN EXPEQ_EXPEQ_iter%03d.IN' % it)
       else:
            it=int(cheasefname[-1][-6:-3])+1

       eqdskrequired    = False
       expeqrequired    = False
       cheaserequired   = False
       exptnzrequired   = False
       iterdbrequired   = False
       profilesrequired = False

       if   rhomesh_update:
            rhomesh_src = srcVals['rhomesh_src'][min(len(srcVals['rhomesh_src'])-1,it+1)]
       if   rhomesh_src  in [0,'chease']:   cheaserequired   = True
       elif rhomesh_src  in [1,'eqdsk']:    eqdskrequired    = True
       elif rhomesh_src  in [2,'expeq']:    expeqrequired    = True

       if   current_update:
            current_src = srcVals['current_src'][min(len(srcVals['current_src'])-1,it+1)]
       if   current_src  in [0,'chease']:   cheaserequired   = True
       elif current_src  in [1,'eqdsk']:    eqdskrequired    = True
       elif current_src  in [2,'expeq']:    expeqrequired    = True

       if not (cheaserequired or eqdskrequired or expeqrequired):
          raise IOError('FATAL: No Geometry File is provided!')

       if   pressure_update:
            pressure_src = srcVals['pressure_src'][min(len(srcVals['pressure_src'])-1,it+1)]
       if   pressure_src in [0,'chease']:   cheaserequired   = True
       elif pressure_src in [1,'eqdsk']:    eqdskrequired    = True
       elif pressure_src in [2,'expeq']:    expeqrequired    = True
       elif pressure_src in [3,'exptnz']:   exptnzrequired   = True
       elif pressure_src in [4,'profiles']: profilesrequired = True
       elif pressure_src in [5,'iterdb']:   iterdbrequired   = True

       if   eprofiles_update:
            eprofiles_src = srcVals['eprofiles_src'][min(len(srcVals['eprofiles_src'])-1,it+1)]
       if   eprofiles_src in [0,'chease']:   cheaserequired   = True
       elif eprofiles_src in [3,'exptnz']:   exptnzrequired   = True
       elif eprofiles_src in [4,'profiles']: profilesrequired = True
       elif eprofiles_src in [5,'iterdb']:   iterdbrequired   = True

       if   iprofiles_update:
            iprofiles_src = srcVals['iprofiles_src'][min(len(srcVals['iprofiles_src'])-1,it+1)]
       if   iprofiles_src in [0,'chease']:   cheaserequired   = True
       elif iprofiles_src in [3,'exptnz']:   exptnzrequired   = True
       elif iprofiles_src in [4,'profiles']: profilesrequired = True
       elif iprofiles_src in [5,'iterdb']:   iterdbrequired   = True

       if not (exptnzrequired or profilesrequired or iterdbrequired or cheaserequired):
          raise IOError('FATAL: No Profiles File is provided!')

       if eqdskrequired    and os.path.isfile('%s_EQDSK'   %(shotfile)):
          if os.path.isfile('%s/%s_EQDSK'   %(shotpath,shotfile)):
             os.system('cp   %s/%s_EQDSK .'     % (shotpath,shotfile))
             eqdskfpath =      '%s_EQDSK'       %          (shotfile)
          elif eqdskrequired:
             raise IOError('EQDSK file NOT FOUND in the given path!')

       if profilesrequired and not os.path.isfile('%s_PROFILES' %(shotfile)):
          if os.path.isfile('%s/%s_PROFILES' %(shotpath,shotfile)):
             os.system('cp   %s/%s_PROFILES .'  % (shotpath,shotfile))
             profilesfpath =   '%s_PROFILES'    %          (shotfile)
          elif profilesrequired:
               raise IOError('Profiles file NOT FOUND in the given path!')

       if iterdbrequired   and not os.path.isfile('%s_ITERDB'   %(shotfile)):
          if os.path.isfile('%s/%s_ITERDB'   %(shotpath,shotfile)):
             os.system('cp   %s/%s_ITERDB .'    % (shotpath,shotfile))
             iterdbfpath =     '%s_ITERDB'      %          (shotfile)
          elif iterdbrequired:
             raise IOError('ITERDB file NOT FOUND in the given path!')

       if 'NPROPT' in namelist:
          current_type = int(namelist['NPROPT'][min(len(namelist['NPROPT'])-1,it+1)])
       else:
          current_type = int(namelist['NSTTP'][min(len(namelist['NSTTP'])-1,it+1)])
       pressure_type   = int(namelist['NPPFUN'][min(len(namelist['NPPFUN'])-1,it+1)])
       rhomesh_type    = int(namelist['NRHOMESH'][min(len(namelist['NRHOMESH'])-1,it+1)])

       expeqParam['nsttp']      = [current_type,current_src]
       expeqParam['nppfun']     = [pressure_type,pressure_src]
       expeqParam['nrhomesh']   = [rhomesh_type,rhomesh_src]
       expeqParam['boundary']   =  boundary_type
       expeqParam['cheasemode'] =  cheasemode

       exptnzParam['nrhomesh']  = [rhomesh_type,rhomesh_src]
       exptnzParam['eprofiles'] = eprofiles_src
       exptnzParam['iprofiles'] = iprofiles_src

       cheasefpath = 'ogyropsi_iter%03d.h5' % it
       cheasedata = read_chease(cheasefpath=cheasefpath)

       ITErr = (cheasedata['ITOR']-ITEXP)/ITEXP
       print('Iter  = ', it)
       print('ITOR  = ', cheasedata['ITOR'])
       print('ITEXP = ', ITEXP)
       print('ITErr = ', abs(ITErr))
    
       while (abs(ITErr) > 1.0e-6):
           if (cheasemode == 1) and (it >= iterTotal): break

           expeqfpath  = 'EXPEQ_iter%03d.OUT'   % it
           exptnzfpath = 'EXPTNZ_iter%03d.OUT'  % it
           cheasefpath = 'ogyropsi_iter%03d.h5' % it

           if   rhomesh_src in [0,'chease']:
                if   pressure_src in [0,'chease']   and current_src in [0,'chease']:
                     expeqdata = write_expeq(setParam=expeqParam,chease=cheasefpath)
                elif pressure_src in [0,'chease']   and current_src in [1,'eqdsk']:
                     expeqdata = write_expeq(setParam=expeqParam,chease=cheasefpath,eqdsk=eqdskfpath)
                elif pressure_src in [0,'chease']   and current_src in [2,'expeq']:
                     expeqdata = write_expeq(setParam=expeqParam,chease=cheasefpath,expeq=expeqfpath)
                elif pressure_src in [1,'eqdsk']    and current_src in [0,'chease']:
                     expeqdata = write_expeq(setParam=expeqParam,chease=cheasefpath,eqdsk=eqdskfpath)
                elif pressure_src in [1,'eqdsk']    and current_src in [1,'eqdsk']:
                     expeqdata = write_expeq(setParam=expeqParam,chease=cheasefpath,eqdsk=eqdskfpath)
                elif pressure_src in [1,'eqdsk']    and current_src in [2,'expeq']:
                     expeqdata = write_expeq(setParam=expeqParam,chease=cheasefpath,expeq=expeqfpath,eqdsk=eqdskfpath)
                elif pressure_src in [2,'expeq']    and current_src in [0,'chease']:
                     expeqdata = write_expeq(setParam=expeqParam,chease=cheasefpath,expeq=expeqfpath)
                elif pressure_src in [2,'expeq']    and current_src in [1,'eqdsk']:
                     expeqdata = write_expeq(setParam=expeqParam,chease=cheasefpath,expeq=expeqfpath,eqdsk=eqdskfpath)
                elif pressure_src in [2,'expeq']    and current_src in [2,'expeq']:
                     expeqdata = write_expeq(setParam=expeqParam,chease=cheasefpath,expeq=expeqfpath)
                elif pressure_src in [3,'exptnz']   and current_src in [0,'chease']:
                     expeqdata = write_expeq(setParam=expeqParam,chease=cheasefpath,exptnz=exptnzfpath)
                elif pressure_src in [3,'exptnz']   and current_src in [1,'eqdsk']:
                     expeqdata = write_expeq(setParam=expeqParam,chease=cheasefpath,eqdsk=eqdskfpath,exptnz=exptnzfpath)
                elif pressure_src in [3,'exptnz']   and current_src in [2,'expeq']:
                     expeqdata = write_expeq(setParam=expeqParam,chease=cheasefpath,expeq=expeqfpath,exptnz=exptnzfpath)
                elif pressure_src in [4,'profiles'] and current_src in [0,'chease']:
                     expeqdata = write_expeq(setParam=expeqParam,chease=cheasefpath,profiles=profilesfpath)
                elif pressure_src in [4,'profiles'] and current_src in [1,'eqdsk']:
                     expeqdata = write_expeq(setParam=expeqParam,chease=cheasefpath,eqdsk=eqdskfpath,profiles=profilesfpath)
                elif pressure_src in [4,'profiles'] and current_src in [2,'expeq']:
                     expeqdata = write_expeq(setParam=expeqParam,chease=cheasefpath,expeq=expeqfpath,profiles=profilesfpath)
                elif pressure_src in [5,'iterdb']   and current_src in [0,'chease']:
                     expeqdata = write_expeq(setParam=expeqParam,chease=cheasefpath,iterdb=iterdbfpath)
                elif pressure_src in [5,'iterdb']   and current_src in [1,'eqdsk']:
                     expeqdata = write_expeq(setParam=expeqParam,chease=cheasefpath,eqdsk=eqdskfpath,iterdb=iterdbfpath)
                elif pressure_src in [5,'iterdb']   and current_src in [2,'expeq']:
                     expeqdata = write_expeq(setParam=expeqParam,chease=cheasefpath,expeq=expeqfpath,iterdb=iterdbfpath)

           elif rhomesh_src in [1,'eqdsk']:
                if   pressure_src in [0,'chease']   and current_src in [0,'chease']:
                     expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,chease=cheasefpath)
                elif pressure_src in [0,'chease']   and current_src in [1,'eqdsk']:
                     expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,chease=cheasefpath)
                elif pressure_src in [0,'chease']   and current_src in [2,'expeq']:
                     expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,expeq=expeqfpath,chease=cheasefpath)
                elif pressure_src in [1,'eqdsk']    and current_src in [0,'chease']:
                     expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,chease=cheasefpath)
                elif pressure_src in [1,'eqdsk']    and current_src in [1,'eqdsk']:
                     expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath)
                elif pressure_src in [1,'eqdsk']    and current_src in [2,'expeq']:
                     expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,expeq=expeqfpath)
                elif pressure_src in [2,'expeq']    and current_src in [0,'chease']:
                     expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,expeq=expeqfpath,chease=cheasefpath)
                elif pressure_src in [2,'expeq']    and current_src in [1,'eqdsk']:
                     expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,expeq=expeqfpath)
                elif pressure_src in [2,'expeq']    and current_src in [2,'expeq']:
                     expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,expeq=expeqfpath)
                elif pressure_src in [3,'exptnz']   and current_src in [0,'chease']:
                     expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,chease=cheasefpath,exptnz=exptnzfpath)
                elif pressure_src in [3,'exptnz']   and current_src in [1,'eqdsk']:
                     expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,exptnz=exptnzfpath)
                elif pressure_src in [3,'exptnz']   and current_src in [2,'expeq']:
                     expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,expeq=expeqfpath,exptnz=exptnzfpath)
                elif pressure_src in [4,'profiles'] and current_src in [0,'chease']:
                     expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,chease=cheasefpath,profiles=profilesfpath)
                elif pressure_src in [4,'profiles'] and current_src in [1,'eqdsk']:
                     expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,profiles=profilesfpath)
                elif pressure_src in [4,'profiles'] and current_src in [2,'expeq']:
                     expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,expeq=expeqfpath,profiles=profilesfpath)
                elif pressure_src in [5,'iterdb']   and current_src in [0,'chease']:
                     expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,chease=cheasefpath,iterdb=iterdbfpath)
                elif pressure_src in [5,'iterdb']   and current_src in [1,'eqdsk']:
                     expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,iterdb=iterdbfpath)
                elif pressure_src in [5,'iterdb']   and current_src in [2,'expeq']:
                     expeqdata = write_expeq(setParam=expeqParam,eqdsk=eqdskfpath,expeq=expeqfpath,iterdb=iterdbfpath)

           if   rhomesh_src in [0,'chease']: 
                if   eprofiles_src in [3,'exptnz']   and iprofiles_src in [3,'exptnz']:
                     exptnzdata = write_exptnz(setParam=exptnzParam,chease=cheasefpath,exptnz=exptnzfpath)
                elif eprofiles_src in [3,'exptnz']   and iprofiles_src in [4,'profiles']:
                     exptnzdata = write_exptnz(setParam=exptnzParam,chease=cheasefpath,exptnz=exptnzfpath,profiles=profilesfpath)
                elif eprofiles_src in [3,'exptnz']   and iprofiles_src in [5,'iterdb']:
                     exptnzdata = write_exptnz(setParam=exptnzParam,chease=cheasefpath,exptnz=exptnzfpath,iterdb=iterdbfpath)
                elif eprofiles_src in [4,'profiles'] and iprofiles_src in [3,'exptnz']:
                     exptnzdata = write_exptnz(setParam=exptnzParam,chease=cheasefpath,profiles=profilesfpath,exptnz=exptnzfpath)
                elif eprofiles_src in [4,'profiles'] and iprofiles_src in [4,'profiles']:
                     exptnzdata = write_exptnz(setParam=exptnzParam,chease=cheasefpath,profiles=profilesfpath)
                elif eprofiles_src in [4,'profiles'] and iprofiles_src in [5,'iterdb']:
                     exptnzdata = write_exptnz(setParam=exptnzParam,chease=cheasefpath,profiles=profilesfpath,iterdb=iterdbfpath)
                elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [3,'exptnz']:
                     exptnzdata = write_exptnz(setParam=exptnzParam,chease=cheasefpath,iterdb=iterdbfpath,exptnz=exptnzfpath)
                elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [4,'profiles']:
                     exptnzdata = write_exptnz(setParam=exptnzParam,chease=cheasefpath,iterdb=iterdbfpath,profiles=profilesfpath)
                elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [5,'iterdb']:
                     exptnzdata = write_exptnz(setParam=exptnzParam,chease=cheasefpath,iterdb=iterdbfpath)
           elif rhomesh_src in [1,'eqdsk']: 
                if   eprofiles_src in [3,'exptnz']   and iprofiles_src in [3,'exptnz']:
                     exptnzdata = write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,exptnz=exptnzfpath)
                elif eprofiles_src in [3,'exptnz']   and iprofiles_src in [4,'profiles']:
                     exptnzdata = write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,exptnz=exptnzfpath,profiles=profilesfpath)
                elif eprofiles_src in [3,'exptnz']   and iprofiles_src in [5,'iterdb']:
                     exptnzdata = write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,exptnz=exptnzfpath,iterdb=iterdbfpath)
                elif eprofiles_src in [4,'profiles'] and iprofiles_src in [3,'exptnz']:
                     exptnzdata = write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,profiles=profilesfpath,exptnz=exptnzfpath)
                elif eprofiles_src in [4,'profiles'] and iprofiles_src in [4,'profiles']:
                     exptnzdata = write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,profiles=profilesfpath)
                elif eprofiles_src in [4,'profiles'] and iprofiles_src in [5,'iterdb']:
                     exptnzdata = write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,profiles=profilesfpath,iterdb=iterdbfpath)
                elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [3,'exptnz']:
                     exptnzdata = write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,iterdb=iterdbfpath,exptnz=exptnzfpath)
                elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [4,'profiles']:
                     exptnzdata = write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,iterdb=iterdbfpath,profiles=profilesfpath)
                elif eprofiles_src in [5,'iterdb']   and iprofiles_src in [5,'iterdb']:
                     exptnzdata = write_exptnz(setParam=exptnzParam,eqdsk=eqdskfpath,iterdb=iterdbfpath)

           namelistParam['QSPEC'] = cheasedata['q'][0]


           namelist = namelistcreate('chease_parameters.csv',min(len(namelist['fname'])-1,it+1),namelistParam)
           os.system('cp ./chease_namelist chease_namelist_iter%3d' % (min(len(namelist['fname'])-1,it+1)))
           exit_status = os.system('./chease_hdf5 chease_namelist > iter%03d.OUT' % (it+1))
          #exit_status = os.system('./chease_hdf5 chease_namelist')
          #exit_status = subprocess.call(['./chease_hdf5','chease_namelist'])
           if abs(exit_status) > 0: sys.exit()
           if os.path.isfile('./ogyropsi.dat'):   os.system('mv ./ogyropsi.dat ogyropsi_iter%03d.dat'     % (it+1))
           if os.path.isfile('./ogyropsi.h5'):    os.system('mv ./ogyropsi.h5 ogyropsi_iter%03d.h5'       % (it+1))
           if os.path.isfile('./EXPEQ.OUT'):      os.system('mv ./EXPEQ.OUT EXPEQ_iter%03d.OUT'           % (it+1))
           if os.path.isfile('./EXPTNZ.OUT'):     os.system('mv ./EXPTNZ.OUT EXPTNZ_iter%03d.OUT'         % (it+1))
           if os.path.isfile('./EXPEQ_EXPEQ.IN'): os.system('mv ./EXPEQ_EXPEQ.IN EXPEQ_EXPEQ_iter%03d.IN' % (it+1))

           cheasepath = 'ogyropsi_iter%03d.h5' % (it+1)
           cheasedata = read_chease(cheasefpath=cheasepath)

           ITErr = (cheasedata['ITOR']-ITEXP)/ITEXP
           print('Iter  = ', it+1)
           print('ITOR = ', cheasedata['ITOR'])
           print('ITEXP = ', ITEXP)
           print('ITErr = ', abs(ITErr))
    
           it+=1
    
       #REMOVING CHEASE FILES NOT NEEDED IN CHEASEPY
       if glob('./NGA'):            os.system('rm NGA')
       if glob('./NDES'):           os.system('rm NDES')
       if glob('./NOUT'):           os.system('rm NOUT')
       if glob('./EQDSK_COCOS_*'):  os.system('rm EQDSK_COCOS_*')
       if glob('./EXPEQ.OUT.TOR'):  os.system('rm EXPEQ.OUT.TOR')
       if glob('./EXPEQ_EQDSK.IN'): os.system('rm EXPEQ_EQDSK.IN')
    
       pltValsKeys = pltVals.keys()
       if 'skipfigs' in pltValsKeys: skipfigs = pltVals['skipfigs']     
       else:                         skipfigs = 0  
       plot_chease('./',skipfigs=0)
    elif selection == 2:
       plot_chease('./',skipfigs=0)
    
    return 1
