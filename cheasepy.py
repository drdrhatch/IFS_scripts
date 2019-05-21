#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import time
import numpy as np  
import efittools
import read_iterdb

import matplotlib.pyplot as plt
plt.switch_backend('agg')

from glob import glob
from scipy.interpolate import interp1d,interp2d
from matplotlib.backends.backend_pdf import PdfPages

def namelistcreate(csvfn,rec,setParam={}):
    infid = open(csvfn, "r")
    import csv
    table = {}
    for row in csv.reader(infid):
        table[row[0]] = row[1:]
    infid.close()

    namelistkeys = table.keys()
    setParamKeys = setParam.keys()

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
    else:                             wfh.write(' NBSEXPQ=1111, ')
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
    if   'NPPFUN'    in setParamKeys: wfh.write(' NPPFUN=%1d,     \n' % (int(setParam['NPPFUN'])))
    elif 'NPPFUN'    in namelistkeys: wfh.write(' NPPFUN=%1d,     \n' % (int(table['NPPFUN'][rec])))
    else:                             wfh.write(' NPPFUN=4,       \n')
    if   'NFUNRHO'   in setParamKeys: wfh.write(' NFUNRHO=%1d, '      % (int(setParam['NFUNRHO'])))
    elif 'NFUNRHO'   in namelistkeys: wfh.write(' NFUNRHO=%1d, '      % (int(table['NFUNRHO'][rec])))
    else:                             wfh.write(' NFUNRHO=0, ')
    if   'NRHOMESH'  in setParamKeys: wfh.write(' NRHOMESH=%1d,   \n' % (int(setParam['NRHOMESH'])))
    elif 'NRHOMESH'  in namelistkeys: wfh.write(' NRHOMESH=%1d,   \n' % (int(table['NRHOMESH'][rec])))
    else:                             wfh.write(' NRHOMESH=0,     \n')
    if   'NSTTP'     in setParamKeys: wfh.write(' NSTTP=%1d, '        % (int(setParam['NSTTP'])))
    elif 'NSTTP'     in namelistkeys: wfh.write(' NSTTP=%1d, '        % (int(table['NSTTP'][rec])))
    else:                             wfh.write(' NSTTP=%1d, '        % (int(3)))
    if   'NPROPT'    in setParamKeys: wfh.write(' NPROPT=%1d,     \n' % (int(setParam['NSTTP'])))
    elif 'NPROPT'    in namelistkeys: wfh.write(' NPROPT=%1d,     \n' % (int(table['NSTTP'][rec])))
    else:                             wfh.write(' NPROPT=%1d,     \n' % (int(1)))
    if   'NVERBOSE'  in setParamKeys: wfh.write(' NVERBOSE=%1d,   \n' % (int(setParam['NVERBOSE'])))
    elif 'NVERBOSE'  in namelistkeys: wfh.write(' NVERBOSE=%1d,   \n' % (int(table['NVERBOSE'][rec])))
    else:                             wfh.write(' NVERBOSE=%1d,   \n' % (int(4)))

    if   'QSPEC'     in setParamKeys: wfh.write(' QSPEC=%3.3f, '      % (float(setParam['QSPEC'])))
    elif 'QSPEC'     in namelistkeys: wfh.write(' QSPEC=%3.3f, '      % (float(table['QSPEC'][rec])))
    if   'CSSPEC'    in setParamKeys: wfh.write(' CSSPEC=%3.3f,   \n' % (float(setParam['CSSPEC'])))
    elif 'CSSPEC'    in namelistkeys: wfh.write(' CSSPEC=%3.3f,   \n' % (float(table['CSSPEC'][rec])))

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
    if   'NPOIDC'    in namelistkeys: wfh.write(' NPOIDC=%1d, '       % (int(table['NPOIDC'][rec])))
    if   'SOLPDC'    in namelistkeys: wfh.write(' SOLPDC=%2.2f,   \n' % (float(table['SOLPDC'][rec])))
    if   'CPLACE'    in namelistkeys: 
       CPLACE = [float(f) for f in table['CPLACE'][rec].split(',')]
       wfh.write(' CPLACE=')
       for cplace in CPLACE:
                                    wfh.write('%4.4f,'              % float(cplace))
       wfh.write('\n')
    if   'CWIDTH'    in namelistkeys: 
       CWIDTH = [float(f) for f in table['CWIDTH'][rec].split(',')]
       wfh.write(' CWIDTH=')
       for cwidth    in CWIDTH:
                                    wfh.write('%4.4f,'              % float(cwidth))
       wfh.write('\n')
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


def CubicSplineDerivative1D(x,fx,dorder=1):
    from scipy.interpolate import CubicSpline
   #CS = CubicSpline(x,fx,extrapolate=bool)
    CS = CubicSpline(x,fx)
    return CS(x,dorder)

def CubicSplineDerivative2D(x,fx,axis=0,dorder=1):
    from scipy.interpolate import CubicSpline
    m,n = np.shape(fx)
    dfx=np.zeros((m,n))
    if   axis == 0:
         for j in range(n):
           #CS = CubicSpline(x,fx[:,j],extrapolate=bool)
            CS = CubicSpline(x,fx[:,j])
            dfx[:,j] = CS(x,dorder)
    elif axis == 1:
         for i in range(m):
           #CS = CubicSpline(x,fx[i,:],extrapolate=bool)
            CS = CubicSpline(x,fx[i,:])
            dfx[i,:] = CS(x,dorder)
    return dfx

def CubicSplineIntegral(x,fxy,intaxis=0):
    from scipy.interpolate import CubicSpline
    if   len(np.shape(fxy)) == 1:
              fy = CubicSpline(x,fxy).integrate(x[0],x[-1])
    elif len(np.shape(fxy)) == 2:
         if   intaxis == 0:
              fy = CubicSpline(x,fxy,axis=0,bc_type='periodic',extrapolate='periodic').integrate(x[0],x[-1],extrapolate='periodic')
         elif intaxis == 1:
              fy = CubicSpline(x,fxy,axis=1).integrate(x[0],x[-1])
    return fy

#def read_profiles(pfpath):
#    from scipy.interpolate import CubicSpline
#    if os.path.isfile(pfpath) == False:
#       print('Fatal: file %s not found.' % pfpath)
#       sys.exit()
#
#    ofh = open(DIIIDpfpath,'r')
#    print('Reading %s file ...' % DIIIDpfpath)
#
#    profiles = {}
#    while True:
#          recs = ofh.readline().split()
#          if len(recs)>=4:
#             nrec = int(recs[0]);
#             var0 = str(recs[1]).lower(); ary0=np.zeros(nrec)
#             var1 = str(recs[2]).lower(); ary1=np.zeros(nrec)
#             var2 = str(recs[3]).lower(); ary2=np.zeros(nrec)
#             for i in range(nrec):
#                 recs    = ofh.readline().split()
#                 ary0[i] = float(recs[0])
#                 ary1[i] = float(recs[1])
#                 ary2[i] = float(recs[2])
#             if var0 in profiles.keys():
#                CS = CubicSpline(ary0,ary1,extrapolate=bool)
#                profiles[var1] = CS(profiles[var0])
#                CS = CubicSpline(ary0,ary2,extrapolate=bool)
#                profiles[var2] = CS(profiles[var0])
#             else:
#                profiles[var0] = ary0
#                profiles[var1] = ary1
#                profiles[var2] = ary2
#          else:
#             break
#
#    ofh.close()
#    return profiles
 

#def profilesPlotter(pfpath):
#    from matplotlib.backends.backend_pdf import PdfPages
#    DIIIDprofiles = read_profiles(pfpath)
#    nvars = sorted(DIIIDprofiles.keys())
#    print('Plotting DIIID profiles in: %s ...' % DIIIDpfpath)
#    DIIIDfigs  = PdfPages('DIIIDprofiles.pdf')
#    for i in nvars:
#        if i in ['n','z','a','psinorm']: continue
#        figname = plt.figure(i)
#        plt.plot(DIIIDprofiles['psinorm'],DIIIDprofiles[i])
#        plt.xlabel('$\Psi$')
#        plt.ylabel(i)
#        DIIIDfigs.savefig(figname)
#    DIIIDfigs.close()
#    return 1
    

def read_cheaseh5(h5fpath):
    import h5py
    from scipy.interpolate import CubicSpline,RectBivariateSpline

    mu0 = 4.0e-7*np.pi

    CHEASEdata                = {}

    if os.path.isfile(h5fpath) == False:
       print('Fatal: file %s not found.' % h5fpath)
       sys.exit()

    hdffh = h5py.File(h5fpath,'r')
    datag = hdffh[hdffh.keys()[0]]


    CHEASEdata['NPSI']        = datag.attrs["NPSI"]
    CHEASEdata['NCHI']        = datag.attrs["NCHI"]
    CHEASEdata['NRBOX']       = datag.attrs["NRBOX"]
    CHEASEdata['NZBOX']       = datag.attrs["NZBOX"]
    CHEASEdata['NBOUND']      = datag.attrs["NBOUND"]
    CHEASEdata['B0EXP']       = datag.attrs["B0EXP"]
    CHEASEdata['R0EXP']       = datag.attrs["R0EXP"]

    CHEASEdata['PSI']         = np.array(datag["grid"]["PSI"])
    CHEASEdata['rhopsi']      = np.sqrt(CHEASEdata['PSI']/CHEASEdata['PSI'][-1])
    CHEASEdata['rhotor']      = np.array(datag["var1d"]["rho_tor"])
    CHEASEdata['drhodpsi']    = CubicSplineDerivative1D(x=CHEASEdata['PSI'],fx=CHEASEdata['rhopsi'])

    CHEASEdata['CHI']         = np.array(datag["grid"]["CHI"])
    CHEASEdata['CHI']         = np.append(CHEASEdata['CHI'],2.0*np.pi)

    CHEASEdata['PSIN']        = (CHEASEdata['PSI']-CHEASEdata['PSI'][0])/(CHEASEdata['PSI'][-1]-CHEASEdata['PSI'][0])
    CHEASEdata['CHIN']        = (CHEASEdata['CHI']-CHEASEdata['CHI'][0])/(CHEASEdata['CHI'][-1]-CHEASEdata['CHI'][0])
    CHEASEdata['rhopsiN']     = (CHEASEdata['rhopsi']-CHEASEdata['rhopsi'][0])/(CHEASEdata['rhopsi'][-1]-CHEASEdata['rhopsi'][0])
    CHEASEdata['rhotorN']     = (CHEASEdata['rhotor']-CHEASEdata['rhotor'][0])/(CHEASEdata['rhotor'][-1]-CHEASEdata['rhotor'][0])

    CHEASEdata['JBS']         = np.array(datag["var1d"]["jbsBav"])
    CHEASEdata['Zeff']        = np.array(datag["var1d"]["zeff"])
    CHEASEdata['kappa']       = np.array(datag["var1d"]["kappa"])
    CHEASEdata['shear']       = np.array(datag["var1d"]["shear"])
    CHEASEdata['signeo']      = np.array(datag["var1d"]["signeo"])
    CHEASEdata['PPrime']      = 2.0*np.pi*np.array(datag["var1d"]["dpdpsi"])

    CHEASEdata['P']           = np.array(datag["var1d"]["p"])
    CHEASEdata['q']           = np.array(datag["var1d"]["q"])
    CHEASEdata['R_av']        = np.array(datag["var1d"]["R_av"])
    CHEASEdata['ageom']       = np.array(datag["var1d"]["ageom"])
    CHEASEdata['Rgeom']       = np.array(datag["var1d"]["Rgeom"])
    CHEASEdata['Volume']      = np.array(datag["var1d"]["Volume"])
    CHEASEdata['GDPSI_av']    = np.array(datag["var1d"]["GDPSI_av"])
    CHEASEdata['radius_av']   = np.array(datag["var1d"]["radius_av"])

    CHEASEdata['T']           = np.array(datag["var1d"]["f"])
    CHEASEdata['TTPrime']     = 2.0*np.pi*np.array(datag["var1d"]["fdfdpsi"])
    CHEASEdata['TPrime']      = CHEASEdata['TTPrime']/CHEASEdata['T']

    CHEASEdata['Te']          = np.array(datag["var1d"]["Te"])
    CHEASEdata['Ti']          = np.array(datag["var1d"]["Ti"])
    CHEASEdata['ne']          = np.array(datag["var1d"]["ne"])
    CHEASEdata['ni']          = np.array(datag["var1d"]["ni"])
    CHEASEdata['TePrime']     = np.array(datag["var1d"]["dTedpsi"])
    CHEASEdata['TiPrime']     = np.array(datag["var1d"]["dTidpsi"])
    CHEASEdata['nePrime']     = np.array(datag["var1d"]["dnedpsi"])
    CHEASEdata['niPrime']     = np.array(datag["var1d"]["dnidpsi"])

    CHEASEdata['rmesh']       = np.array(datag["var1d"]["rmesh"])
    CHEASEdata['zmesh']       = np.array(datag["var1d"]["zmesh"])
    CHEASEdata['rbound']      = np.array(datag["var1d"]["rboundplasma"])
    CHEASEdata['zbound']      = np.array(datag["var1d"]["zboundplasma"])
    CHEASEdata['delta_upper'] = np.array(datag["var1d"]["delta_upper"])
    CHEASEdata['delta_lower'] = np.array(datag["var1d"]["delta_lower"])

    CHEASEdata['rmeshN']      = CHEASEdata['rmesh']/CHEASEdata['R0EXP'] 
    CHEASEdata['zmeshN']      = CHEASEdata['zmesh']/CHEASEdata['R0EXP'] 
    CHEASEdata['rboundN']     = CHEASEdata['rbound']/CHEASEdata['R0EXP'] 
    CHEASEdata['zboundN']     = CHEASEdata['zbound']/CHEASEdata['R0EXP']
 
    CHEASEdata['dpdpsi']      = np.array(datag["var1d"]["dpdpsi"])
    CHEASEdata['dqdpsi']      = np.array(datag["var1d"]["dqdpsi"])
    CHEASEdata['dVdpsi']      = np.array(datag["var1d"]["dVdpsi"])
    CHEASEdata['d2qdpsi2']    = np.array(datag["var1d"]["d2qdpsi2"])
    CHEASEdata['dsheardpsi']  = np.array(datag["var1d"]["dsheardpsi"])
    CHEASEdata['dpsidrhotor'] = np.array(datag["var1d"]["dpsidrhotor"])


   #THE DIMENSION OF ALL THE FOLLOWING QUNATITIES ARE (NCHI,NPSI)
    CHEASEdata['R']           = np.array(datag["var2d"]["R"])
    CHEASEdata['Z']           = np.array(datag["var2d"]["Z"])
    CHEASEdata['B']           = np.array(datag["var2d"]["B"])
    CHEASEdata['J']           = np.array(datag["var2d"]["Jacobian"])
    CHEASEdata['g11']         = np.array(datag["var2d"]["g11"])
    CHEASEdata['g22']         = np.array(datag["var2d"]["g22"])
    CHEASEdata['g33']         = np.array(datag["var2d"]["g33"])
    CHEASEdata['dBdpsi']      = np.array(datag["var2d"]["dBdpsi"])
    CHEASEdata['dBdchi']      = np.array(datag["var2d"]["dBdchi"])
    CHEASEdata['dChidZ']      = np.array(datag["var2d"]["dChidZ"])
    CHEASEdata['dPsidZ']      = np.array(datag["var2d"]["dPsidZ"])
    CHEASEdata['dChidR']      = np.array(datag["var2d"]["dChidR"])
    CHEASEdata['dPsidR']      = np.array(datag["var2d"]["dPsidR"])

    CHEASEdata['R']           = np.vstack((CHEASEdata['R'],CHEASEdata['R'][0,:]))
    CHEASEdata['Z']           = np.vstack((CHEASEdata['Z'],CHEASEdata['Z'][0,:]))
    CHEASEdata['B']           = np.vstack((CHEASEdata['B'],CHEASEdata['B'][0,:]))
    CHEASEdata['J']           = np.vstack((CHEASEdata['J'],CHEASEdata['J'][0,:]))
    CHEASEdata['g11']         = np.vstack((CHEASEdata['g11'],CHEASEdata['g11'][0,:]))
    CHEASEdata['g22']         = np.vstack((CHEASEdata['g22'],CHEASEdata['g22'][0,:]))
    CHEASEdata['g33']         = np.vstack((CHEASEdata['g33'],CHEASEdata['g33'][0,:]))
    CHEASEdata['dBdpsi']      = np.vstack((CHEASEdata['dBdpsi'],CHEASEdata['dBdpsi'][0,:]))
    CHEASEdata['dBdchi']      = np.vstack((CHEASEdata['dBdchi'],CHEASEdata['dBdchi'][0,:]))
    CHEASEdata['dChidZ']      = np.vstack((CHEASEdata['dChidZ'],CHEASEdata['dChidZ'][0,:]))
    CHEASEdata['dPsidZ']      = np.vstack((CHEASEdata['dPsidZ'],CHEASEdata['dPsidZ'][0,:]))
    CHEASEdata['dChidR']      = np.vstack((CHEASEdata['dChidR'],CHEASEdata['dChidR'][0,:]))
    CHEASEdata['dPsidR']      = np.vstack((CHEASEdata['dPsidR'],CHEASEdata['dPsidR'][0,:]))

   #THE DIMENSION OF ALL THE FOLLOWING QUNATITIES ARE (NRBOX,NZBOX)
    CHEASEdata['psiRZ']       = np.array(datag["var2d"]["psiRZ"])
    CHEASEdata['chiRZ']       = np.array(datag["var2d"]["chiRZ"])

    CHEASEdata['RN']          = CHEASEdata['R']/CHEASEdata['R0EXP']
    CHEASEdata['ZN']          = CHEASEdata['Z']/CHEASEdata['R0EXP']
    CHEASEdata['BN']          = CHEASEdata['B']/CHEASEdata['B0EXP']

    CHEASEdata['C0']          = CubicSplineIntegral(x=CHEASEdata['CHI'],fxy=CHEASEdata['J']/CHEASEdata['R'],                     intaxis=0)
    CHEASEdata['C1']          = CubicSplineIntegral(x=CHEASEdata['CHI'],fxy=CHEASEdata['J'],                                     intaxis=0)
    CHEASEdata['C2']          = CubicSplineIntegral(x=CHEASEdata['CHI'],fxy=CHEASEdata['J']/CHEASEdata['R']**2,                  intaxis=0)
    CHEASEdata['C3']          = CubicSplineIntegral(x=CHEASEdata['CHI'],fxy=CHEASEdata['J']*CHEASEdata['g11']*CHEASEdata['g33'], intaxis=0)

    CHEASEdata['y1']          = 1.0+CHEASEdata['C3']/CHEASEdata['C2']/CHEASEdata['T']**2/4.0/np.pi**2

    CHEASEdata['PN']          = mu0*CHEASEdata['P']/CHEASEdata['B0EXP']**2

    CHEASEdata['<B2>']        = CubicSplineIntegral(x=CHEASEdata['CHI'],fxy=CHEASEdata['J']*CHEASEdata['B']**2,intaxis=0)/CHEASEdata['C1']
    CHEASEdata['<JdotB>']     =-CHEASEdata['T']*CHEASEdata['PPrime']-CHEASEdata['TPrime']*CHEASEdata['<B2>']/mu0
    CHEASEdata['JPRL']        = CHEASEdata['<JdotB>']/CHEASEdata['B0EXP']
    CHEASEdata['JPRLN']       = CHEASEdata['JPRL']*mu0*CHEASEdata['R0EXP']/CHEASEdata['B0EXP']

    CHEASEdata['<T/R2>']      = CubicSplineIntegral(x=CHEASEdata['CHI'],fxy=CHEASEdata['J']*CHEASEdata['T']*CHEASEdata['g33'],intaxis=0)/CHEASEdata['C1']
    CHEASEdata['IPRL']        = CHEASEdata['R0EXP']*CHEASEdata['<JdotB>']/CHEASEdata['<T/R2>']
    CHEASEdata['IPRLN']       = CHEASEdata['IPRL']*mu0/CHEASEdata['R0EXP']/CHEASEdata['B0EXP']

    CHEASEdata['ISTR']        =-((CHEASEdata['C2']/CHEASEdata['C0'])*(CHEASEdata['TTPrime']/mu0))
    CHEASEdata['ISTR']       +=-((CHEASEdata['C1']/CHEASEdata['C0'])*CHEASEdata['PPrime'])
    CHEASEdata['ISTR']       *= CHEASEdata['R0EXP']**2
    CHEASEdata['ISTRN']       = CHEASEdata['ISTR']*mu0/CHEASEdata['R0EXP']/CHEASEdata['B0EXP']

    CHEASEdata['JPHI']        =-(CHEASEdata['R']*CHEASEdata['PPrime'])-(CHEASEdata['TTPrime']/(mu0*CHEASEdata['R']))
    CHEASEdata['JPHIN']       = CHEASEdata['JPHI']*mu0*CHEASEdata['R0EXP']/CHEASEdata['B0EXP']

    PSImin                    = CHEASEdata['PSI'][0]
    PSImax                    = CHEASEdata['PSI'][-1]
    CHImin                    = CHEASEdata['CHI'][0]
    CHImax                    = CHEASEdata['CHI'][-1]
    fchipsi                   = CHEASEdata['JPHI']*CHEASEdata['J']/CHEASEdata['R']
    CHEASEdata['ITOR']        = RectBivariateSpline(CHEASEdata['CHI'],CHEASEdata['PSI'],fchipsi).integral(CHImin,CHImax,PSImin,PSImax)
   #CHEASEdata['ITOR']        = abs(RectBivariateSpline(CHEASEdata['CHI'],CHEASEdata['PSI'],fchipsi).integral(CHImin,CHImax,PSImin,PSImax))
    CHEASEdata['ITORN']       = CHEASEdata['ITOR']*mu0/CHEASEdata['R0EXP']/CHEASEdata['B0EXP']


    CHEASEdata['IBS']         = CHEASEdata['R0EXP']*CHEASEdata['JBS']/CHEASEdata['<T/R2>']
    CHEASEdata['IBSN']        = CHEASEdata['IBS']*mu0/CHEASEdata['R0EXP']/CHEASEdata['B0EXP']

    CHEASEdata['JBSN']        = CHEASEdata['JBS'] *mu0*CHEASEdata['R0EXP']/CHEASEdata['B0EXP']

    CHEASEdata['PPrimeN']     = CHEASEdata['PPrime']*mu0*CHEASEdata['R0EXP']**2/CHEASEdata['B0EXP']
    CHEASEdata['TTPrimeN']    = CHEASEdata['TTPrime']/CHEASEdata['B0EXP']

    return CHEASEdata


def read_expeq(expeqfpath):
    if os.path.isfile(expeqfpath) == False:
       print('Fatal: file %s not found.' % expeqfpath)
       sys.exit()

    ofh = open(expeqfpath,'r')
    EXPEQOUT = ofh.readlines()
    ofh.close()

    EXPEQdata                  = {} 
    EXPEQdata['aspect']        = float(EXPEQOUT[0])
    EXPEQdata['zgeom']         = float(EXPEQOUT[1])
    EXPEQdata['pedge']         = float(EXPEQOUT[2])
    nRZmesh                    =   int(EXPEQOUT[3])
    EXPEQdata['rbound']         = np.array([irec.split()[0] for irec in EXPEQOUT[4:nRZmesh+4]],dtype=float)
    EXPEQdata['zbound']         = np.array([irec.split()[1] for irec in EXPEQOUT[4:nRZmesh+4]],dtype=float)
    
    nrhomesh                   = int(EXPEQOUT[nRZmesh+4].split()[0])
    EXPEQdata['nppfun']        = int(EXPEQOUT[nRZmesh+4].split()[1])
    EXPEQdata['nsttp']         = int(EXPEQOUT[nRZmesh+5].split()[0])
    EXPEQdata['nrhotype']      = int(EXPEQOUT[nRZmesh+5].split()[1])

    if   EXPEQdata['nrhotype'] == 0:
         EXPEQdata['rhopsiN']  = np.array(EXPEQOUT[nRZmesh+6+0*nrhomesh:nRZmesh+6+1*nrhomesh],dtype=float)
    elif EXPEQdata['nrhotype'] == 1:
         EXPEQdata['rhotorN']  = np.array(EXPEQOUT[nRZmesh+6+0*nrhomesh:nRZmesh+6+1*nrhomesh],dtype=float)
    
    if   EXPEQdata['nppfun']   == 4:
         EXPEQdata['pprimeN']  = np.array(EXPEQOUT[nRZmesh+6+1*nrhomesh:nRZmesh+6+2*nrhomesh],dtype=float)
    elif EXPEQdata['nppfun']   == 8:
         EXPEQdata['pN']       = np.array(EXPEQOUT[nRZmesh+6+1*nrhomesh:nRZmesh+6+2*nrhomesh],dtype=float)
  
    if   EXPEQdata['nsttp']    == 1:
         EXPEQdata['ttprimeN'] = np.array(EXPEQOUT[nRZmesh+6+2*nrhomesh:nRZmesh+6+3*nrhomesh],dtype=float)
    elif EXPEQdata['nsttp']    == 2:
         EXPEQdata['istrN']    = np.array(EXPEQOUT[nRZmesh+6+2*nrhomesh:nRZmesh+6+3*nrhomesh],dtype=float)
    elif EXPEQdata['nsttp']    == 3:
         EXPEQdata['iprlN']    = np.array(EXPEQOUT[nRZmesh+6+2*nrhomesh:nRZmesh+6+3*nrhomesh],dtype=float)
    elif EXPEQdata['nsttp']    == 4:
         EXPEQdata['jprlN']    = np.array(EXPEQOUT[nRZmesh+6+2*nrhomesh:nRZmesh+6+3*nrhomesh],dtype=float)

    return EXPEQdata


def write_expeq(h5fpath="",expeqfpath="",exptnzfpath="",setParam={}):
    if h5fpath:
       CHEASEdata = read_cheaseh5(h5fpath)
       CHEASEflag = True
    if exptnzfpath:
       EXPTNZdata = read_exptnz(exptnzfpath)
       EXPTNZflag = True
    if expeqfpath:
       EXPEQdata  = read_expeq(expeqfpath)
       EXPEQflag  = True

    expeq            = {}
    expeq['aspect']  = (max(CHEASEdata['rbound'])-min(CHEASEdata['rbound']))
    expeq['aspect'] /= (max(CHEASEdata['rbound'])+min(CHEASEdata['rbound']))
    expeq['zgeom']   = np.mean(CHEASEdata['zmesh'])/CHEASEdata['R0EXP']
    expeq['pedge']   = CHEASEdata['PN'][-1]

    if 'msbound' in setParam.keys():
       expeq['msbound'] = setParam['msbound']
    else:
       expeq['msbound'] = 0
    if 'nsttp' in setParam.keys():
       expeq['nsttp'] = setParam['nsttp']
    else:
       expeq['nsttp']    = [1,0]
    if 'nppfun' in setParam.keys():
       expeq['nppfun'] = setParam['nppfun']
    else:
       expeq['nppfun']   = [4,0]
    if 'nrhomesh' in setParam.keys():
       expeq['nrhotype'] = setParam['nrhotype']
    else:
       expeq['nrhotype'] = [0,0]
    expeq['nrhomesh']    = np.size(CHEASEdata['rhopsiN'])

    if   expeq['msbound'] == 0 or expeq['msbound'] =='h5':
         expeq['nRZmesh'] = np.size(CHEASEdata['rboundN'])
         expeq['rbound']  = CHEASEdata['rboundN']
         expeq['zbound']  = CHEASEdata['zboundN']
    elif expeq['msbound'] == 1 or expeq['msbound'] =='expeq':
         expeq['nRZmesh'] = np.size(EXPEQdata['rbound'])
         expeq['rbound']  = EXPEQdata['rbound']
         expeq['zbound']  = EXPEQdata['zbound']

    ofh = open("EXPEQ",'w')
    ofh.write('%18.8E\n'           % expeq['aspect'])
    ofh.write('%18.8E\n'           % expeq['zgeom'])
    ofh.write('%18.8E\n'           % expeq['pedge'])
    ofh.write('%5d\n'              % expeq['nRZmesh'])
    for i in range(expeq['nRZmesh']):
        ofh.write('%18.8E%18.8E\n' % (expeq['rbound'][i],expeq['zbound'][i]))
    ofh.write('%5d%5d\n'           % (expeq['nrhomesh'],expeq['nppfun'][0]))
    ofh.write('%5d%5d\n'           % (expeq['nsttp'][0],expeq['nrhotype'][0]))

    if expeq['nrhotype'][0] == 0 or expeq['nrhotype'][0] == 'rhopsi':
       if   expeq['nrhotype'][1] == 0 or expeq['nrhotype'][1] == 'h5':
          for i in range(np.size(CHEASEdata['rhopsiN'])):
              ofh.write('%18.8E\n'       % CHEASEdata['rhopsiN'][i])
       elif expeq['nrhotype'][1] == 1 or expeq['nrhotype'][1] == 'expeq':
          for i in range(np.size(EXPEQdata['rhopsiN'])):
              ofh.write('%18.8E\n'       % EXPEQdata['rhopsiN'][i])
       elif expeq['nrhotype'][1] == 2 or expeq['nrhotype'][1] == 'exptnz':
          for i in range(np.size(EXPTNZdata['rhopsiN'])):
              ofh.write('%18.8E\n'       % EXPTNZdata['rhopsiN'][i])
    elif expeq['nrhotype'][0] == 1 or expeq['nrhotype'][0] == 'rhotor':
       if   expeq['nrhotype'][1] == 0 or expeq['nrhotype'][1] == 'h5':
          for i in range(np.size(CHEASEdata['rhotorN'])):
              ofh.write('%18.8E\n'       % CHEASEdata['rhotorN'][i])
       elif expeq['nrhotype'][1] == 1 or expeq['nrhotype'][1] == 'expeq':
          for i in range(np.size(EXPEQdata['rhotorN'])):
              ofh.write('%18.8E\n'       % EXPEQdata['rhotorN'][i])
       elif expeq['nrhotype'][1] == 2 or expeq['nrhotype'][1] == 'exptnz':
          for i in range(np.size(EXPTNZdata['rhotorN'])):
              ofh.write('%18.8E\n'       % EXPTNZdata['rhotorN'][i])

    if expeq['nppfun'][0] == 4 or expeq['nppfun'][0] == 'p':
       if   expeq['nppfun'][1] == 0 or expeq['nppfun'][1] == 'h5':
            expeq['PPrime']=CHEASEdata['PPrimeN']
            for i in range(np.size(CHEASEdata['PPrimeN'])):
                ofh.write('%18.8E\n'       % CHEASEdata['PPrimeN'][i])
       elif expeq['nppfun'][1] == 1 or expeq['nppfun'][1] == 'expeq':
            for i in range(np.size(EXPEQdata['pprimeN'])):
                ofh.write('%18.8E\n'       % EXPEQdata['pprimeN'][i])
       elif expeq['nppfun'][1] == 2 or expeq['nppfun'][1] == 'exptnz':
            PtPrime  = CHEASEdata['TePrime']*CHEASEdata['ne']+CHEASEdata['Te']*CHEASEdata['nePrime']
            PtPrime += CHEASEdata['TiPrime']*CHEASEdata['ni']+CHEASEdata['Ti']*CHEASEdata['niPrime']
            PtPrime *= 1.602e-19
            PtPrimeN = 2.0*np.pi*PtPrime*(4.0e-7*np.pi)*CHEASEdata['R0EXP']**2/CHEASEdata['B0EXP']
            expeq['PPrime']=PtPrimeN
            for i in range(np.size(CHEASEdata['PtPrimeN'])):
                ofh.write('%18.8E\n'       % PtPrimeN[i])
    elif expeq['nppfun'][0] == 8 or expeq['nppfun'][0] == 'pprime':
       if   expeq['nppfun'][1] == 0 or expeq['nppfun'][1] == 'h5':
            for i in range(np.size(CHEASEdata['PN'])):
                ofh.write('%18.8E\n'       % CHEASEdata['PN'][i])
       elif expeq['nppfun'][1] == 1 or expeq['nppfun'][1] == 'expeq':
            for i in range(np.size(EXPEQdata['pN'])):
                ofh.write('%18.8E\n'       % EXPEQdata['pN'][i])
       elif expeq['nppfun'][1] == 2 or expeq['nppfun'][1] == 'exptnz':
            pN  = CHEASEdata['exptnz_Te']*CHEASEdata['exptnz_ne']
            pN += CHEASEdata['exptnz_Ti']*CHEASEdata['exptnz_ni']
            pN += CHEASEdata['Ti']*(CHEASEdata['ne'] - CHEASEdata['ni'])
            pN *= 1.602e-19*(4.0e-7*np.pi)/CHEASEdata['B0EXP']**2
            for i in range(np.size(pN)):
                ofh.write('%18.8E\n'       % pN[i])

    if expeq['nsttp'][0] == 1 or expeq['nsttp'][0] == 'ttprime':
         if expeq['nrhotype'][0] == 1 or expeq['nrhotype'] == 'rhotor':
            print 'Runtime FATAL Error: nrhotype (must) = 0 or rhopsi'
            print 'CheasePy stopped running!'
            sys.exit()
         if setParam['cheasemode'] in [2,3]:
            ITErr = (CHEASEdata['ITOR']-setParam['ITEXP'])/setParam['ITEXP']
            if   expeq['nsttp'][1] == 0 or expeq['nsttp'][1] == 'h5':
                 if   setParam['cheasemode'] == 2:
                      IOHMIC  = (1.0-ITErr)*(CHEASEdata['IPRLN']-CHEASEdata['IBSN'])
                      IPRLN   = CHEASEdata['IBSN'] + IOHMIC
            elif expeq['nsttp'][1] == 1 or expeq['nsttp'][1] == 'expeq':
                 if   setParam['cheasemode'] == 2:
                      IOHMIC  = (1.0-ITErr)*(EXPEQdata['iprlN']-CHEASEdata['IBSN'])
                      IPRLN   = CHEASEdata['IBSN'] + IOHMIC
            y0       = (1.0+CHEASEdata['C3']/CHEASEdata['C2']/CHEASEdata['T']**2/4.0/np.pi**2)*CHEASEdata['R']
            y1       = (CHEASEdata['C1']/CHEASEdata['C2'])-(y0*CHEASEdata['R']**2)
            JPHIN    = (IPRLN+CHEASEdata['PPrime']*y1)/y0
            TTPM     =-JPHIN*CHEASEdata['R']-CHEASEdata['PPrimeN']*CHEASEdata['R']**2
            TTPMN    = CubicSplineIntegral(x=CHEASEdata['CHI'],fxy=TTPM*CHEASEdata['J'], intaxis=0)
            for i in range(np.size(TTPMN)):
                ofh.write('%18.8E\n'       % TTPMN[i])
         else:
            if   expeq['nsttp'][1] == 0 or expeq['nsttp'][1] == 'h5':
                 for i in range(np.size(CHEASEdata['TTPrimeN'])):
                     ofh.write('%18.8E\n'       % CHEASEdata['TTPrimeN'][i])
            elif expeq['nsttp'][1] == 1 or expeq['nsttp'][1] == 'expeq':
                 for i in range(np.size(EXPEQdata['ttprimeN'])):
                     ofh.write('%18.8E\n'       % EXPEQdata['ttprimeN'][i])
    elif expeq['nsttp'][0] == 2 or expeq['nsttp'][0] == 'istar':
         if setParam['cheasemode'] in [2,3]:
            ITErr = (CHEASEdata['ITOR']-setParam['ITEXP'])/setParam['ITEXP']
            if   expeq['nsttp'][1] == 0 or expeq['nsttp'][1] == 'h5':
                 if   setParam['cheasemode'] == 2:
                      IOHMIC  = (1.0-ITErr)*(CHEASEdata['IPRLN']-CHEASEdata['IBSN'])
                      IPRLN   = CHEASEdata['IBSN'] + IOHMIC
            elif expeq['nsttp'][1] == 1 or expeq['nsttp'][1] == 'expeq':
                 if   setParam['cheasemode'] == 2:
                      IOHMIC  = (1.0-ITErr)*(EXPEQdata['iprlN']-CHEASEdata['IBSN'])
                      IPRLN   = CHEASEdata['IBSN'] + IOHMIC
            y0       = (1.0+CHEASEdata['C3']/CHEASEdata['C2']/CHEASEdata['T']**2/4.0/np.pi**2)*CHEASEdata['R']
            y1       = (CHEASEdata['C1']/CHEASEdata['C2'])-(y0*CHEASEdata['R']**2)
            JPHIN    = (IPRLN+CHEASEdata['PPrimeN']*y1)/y0
            INTGR1  = CubicSplineIntegral(x=CHEASEdata['CHI'],fxy=JPHIN*CHEASEdata['J']/CHEASEdata['R'], intaxis=0)
            INTGR2  = CubicSplineIntegral(x=CHEASEdata['CHI'],fxy=      CHEASEdata['J']/CHEASEdata['R'], intaxis=0)
            ISTRN    = INTGR1/INTGR2
            for i in range(np.size(ISTRN)):
                ofh.write('%18.8E\n'       % ISTRN[i])
         else:
            if   expeq['nsttp'][1] == 0 or expeq['nsttp'][1] == 'h5':
                 for i in range(np.size(CHEASEdata['ISTRN'])):
                     ofh.write('%18.8E\n'       % CHEASEdata['ISTRN'][i])
            elif expeq['nsttp'][1] == 1 or expeq['nsttp'][1] == 'expeq':
                 for i in range(np.size(EXPEQdata['istrN'])):
                     ofh.write('%18.8E\n'       % EXPEQdata['istrN'][i])
    elif expeq['nsttp'][0] == 3 or expeq['nsttp'][0] == 'iparallel':
         if setParam['cheasemode'] in [2,3]:
            ITErr = (CHEASEdata['ITOR']-setParam['ITEXP'])/setParam['ITEXP']
           #if    ITErr/write_expeq.errorval > 1.0: ITErr = 2.0*write_expeq.errorval
           #elif  ITErr/write_expeq.errorval < 0.5: ITErr = 2.0*ITErr
           #write_expeq.errorval = ITErr
            if   expeq['nsttp'][1] == 0 or expeq['nsttp'][1] == 'h5':
                 if   setParam['cheasemode'] == 2:
                      IOHMIC  = (1.0-ITErr)*(CHEASEdata['IPRLN']-CHEASEdata['IBSN'])
                      IPRLN   = CHEASEdata['IBSN'] + IOHMIC
            elif expeq['nsttp'][1] == 1 or expeq['nsttp'][1] == 'expeq':
                 if   setParam['cheasemode'] == 2:
                      IOHMIC  = (1.0-ITErr)*(EXPEQdata['iprlN']-CHEASEdata['IBSN'])
                      IPRLN   = CHEASEdata['IBSN'] + IOHMIC
            for i in range(np.size(IPRLN)):
                ofh.write('%18.8E\n'       % IPRLN[i])
         else:
            if   expeq['nsttp'][1] == 0 or expeq['nsttp'][1] == 'h5':
                 for i in range(np.size(CHEASEdata['IPRLN'])):
                     ofh.write('%18.8E\n'       % CHEASEdata['IPRLN'][i])
            elif expeq['nsttp'][1] == 1 or expeq['nsttp'][1] == 'expeq':
                 for i in range(np.size(EXPEQdata['iprlN'])):
                     ofh.write('%18.8E\n'       % EXPEQdata['iprlN'][i])
    elif expeq['nsttp'][0] == 4 or expeq['nsttp'][0] == 'jparallel':
         if setParam['cheasemode'] in [2,3]:
            ITErr = (CHEASEdata['ITOR']-setParam['ITEXP'])/setParam['ITEXP']
            if   expeq['nsttp'][1] == 0 or expeq['nsttp'][1] == 'h5':
                 if   setParam['cheasemode'] == 2:
                      JOHMIC  = (1.0-ITErr)*(CHEASEdata['JPRLN']-CHEASEdata['JBSN'])
                      JPRLN   = CHEASEdata['JBSN'] + JOHMIC
            elif expeq['nsttp'][1] == 1 or expeq['nsttp'][1] == 'expeq':
                 if   setParam['cheasemode'] == 2:
                      JOHMIC  = (1.0-ITErr)*(EXPEQdata['jprlN']-CHEASEdata['JBSN'])
                      JPRLN   = CHEASEdata['JBSN'] + JOHMIC
            for i in range(np.size(JPRLN)):
                 ofh.write('%18.8E\n'       % JPRLN[i])
         else:
            if   expeq['nsttp'][1] == 0 or expeq['nsttp'][1] == 'h5':
                 for i in range(np.size(CHEASEdata['JPRLN'])):
                     ofh.write('%18.8E\n'       % CHEASEdata['JPRLN'][i])
            elif expeq['nsttp'][1] == 1 or expeq['nsttp'][1] == 'expeq':
                 for i in range(np.size(EXPEQdata['jprlN'])):
                     ofh.write('%18.8E\n'       % EXPEQdata['jprlN'][i])

    ofh.close()
    return expeq


def read_exptnz(exptnzfpath):
    if os.path.isfile(exptnzfpath) == False:
       print('Fatal: file %s not found.' % exptnzfpath)
       sys.exit()

    ofh = open(exptnzfpath,'r')
    EXPTNZOUT = ofh.readlines()
    ofh.close()

    n_rho   = int(EXPTNZOUT[0].split()[0])
    rhotype =     EXPTNZOUT[0].split()[1].strip()[0:6]
    EXPTNZdata             = {}
    if rhotype=='rhopsi':
       EXPTNZdata['rhopsiN']  = np.array(EXPTNZOUT[0*n_rho+1:1*n_rho+1],dtype=float)
    else:
       EXPTNZdata['rhotorN']  = np.array(EXPTNZOUT[0*n_rho+1:1*n_rho+1],dtype=float)
    EXPTNZdata['Te']          = np.array(EXPTNZOUT[1*n_rho+1:2*n_rho+1],dtype=float)
    EXPTNZdata['ne']          = np.array(EXPTNZOUT[2*n_rho+1:3*n_rho+1],dtype=float)
    EXPTNZdata['Zeff']        = np.array(EXPTNZOUT[3*n_rho+1:4*n_rho+1],dtype=float)
    EXPTNZdata['Ti']          = np.array(EXPTNZOUT[4*n_rho+1:5*n_rho+1],dtype=float)
    EXPTNZdata['ni']          = np.array(EXPTNZOUT[5*n_rho+1:6*n_rho+1],dtype=float)

    return EXPTNZdata

def write_exptnz(h5fpath='',exptnzfpath='',setParam={}):
    #nrhotype=[rhotype(0:rhopsi,1:rhotor),rhosrc(0:chease,1:exptnz)]
    #profiles=[eprofilesrc(0:chease,1:exptnz),iprofilesrc(0:chease,1:exptnz)]
    if not h5fpath and exptnzfpath:
       EXPTNZdata = read_exptnz(exptnzfpath)
       if 'nrhomesh' in setParam.keys():
           nrhotype = [setParam['nrhomesh'][0],1]
       else:
           nrhotype = [0,1]
       profiles   = [1,1]
    elif h5fpath and not exptnzfpath:
       CHEASEdata = read_cheaseh5(h5fpath)
       if 'nrhomesh' in setParam.keys():
           nrhotype = [setParam['nrhomesh'][0],0]
       else:
           nrhotype = [0,0]
       profiles   = [0,0]
    else:
       CHEASEdata = read_cheaseh5(h5fpath)
       EXPTNZdata = read_exptnz(exptnzfpath)
       if 'nrhomesh' in setParam.keys():
           nrhotype = setParam['nrhomesh']
       else:
           nrhotype = [0,1]
       if 'profiles' in setParam.keys():
           profiles = setParam['profiles']
       else:
           profiles   = [1,1]

    ofh = open("EXPTNZ",'w')
    if   nrhotype[1] == 0 or nrhotype[1] == 'h5':
         if   nrhotype[0] == 0 or nrhotype[0] == 'rhopsi':
              n_rho = np.size(CHEASEdata['rhopsiN'])
              ofh.write('%5d rhopsi,  Te,   ne,   Zeff,   Ti,   ni  profiles\n' % n_rho)
              for i in range(n_rho): ofh.write('%16.6E\n' % CHEASEdata['rhopsiN'][i])
         elif nrhotype[0] == 1 or nrhotype[0] == 'rhotor':
              n_rho = np.size(CHEASEdata['rhotorN'])
              ofh.write('%5d rhotor,  Te,   ne,   Zeff,   Ti,   ni  profiles\n' % n_rho)
              for i in range(n_rho): ofh.write('%16.6E\n' % CHEASEdata['rhotorN'][i])
    elif nrhotype[1] == 1 or nrhotype[1] == 'exptnz':
         if   nrhotype[0] == 0 or nrhotype[0] == 'rhopsi':
              n_rho = np.size(EXPTNZdata['rhopsiN'])
              ofh.write('%5d rhopsi,  Te,   ne,   Zeff,   Ti,   ni  profiles\n' % n_rho)
              for i in range(n_rho): ofh.write('%16.6E\n' % EXPTNZdata['rhopsiN'][i])
         elif nrhotype[0] == 1 or nrhotype[0] == 'rhotor':
              n_rho = np.size(EXPTNZdata['rhotorN'])
              ofh.write('%5d rhotor,  Te,   ne,   Zeff,   Ti,   ni  profiles\n' % n_rho)
              for i in range(n_rho): ofh.write('%16.6E\n' % EXPTNZdata['rhotorN'][i])

    if   profiles[0] == 0 or profiles[0] == 'h5':
         for i in range(n_rho): ofh.write('%16.6E\n' % CHEASEdata['Te'][i])
         for i in range(n_rho): ofh.write('%16.6E\n' % CHEASEdata['ne'][i])
    elif profiles[0] == 1 or profiles[0] == 'exptnz':
         for i in range(n_rho): ofh.write('%16.6E\n' % EXPTNZdata['Te'][i])
         for i in range(n_rho): ofh.write('%16.6E\n' % EXPTNZdata['ne'][i])

    if   profiles[1] == 0 or profiles[1] == 'h5':
         for i in range(n_rho): ofh.write('%16.6E\n' % CHEASEdata['Zeff'][i])
         for i in range(n_rho): ofh.write('%16.6E\n' % CHEASEdata['Ti'][i])
         for i in range(n_rho): ofh.write('%16.6E\n' % CHEASEdata['ni'][i])
    elif profiles[1] == 1 or profiles[1] == 'exptnz':
         for i in range(n_rho): ofh.write('%16.6E\n' % EXPTNZdata['Zeff'][i])
         for i in range(n_rho): ofh.write('%16.6E\n' % EXPTNZdata['Ti'][i])
         for i in range(n_rho): ofh.write('%16.6E\n' % EXPTNZdata['ni'][i])

    ofh.close()

    return 1



def plot_cheasedata(OSPATH,reportpath='',skipfigs=1):
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
       shotnam    = sorted(glob(srhpath))

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

           CHEASEdata = read_cheaseh5(h5fid)
           CHEASEdataKeys = CHEASEdata.keys()

           EXPTNZdata = read_exptnz(exptnzlist[h5list.index(h5fid)])
           EXPEQdata  = read_expeq(expeqlist[h5list.index(h5fid)])
           EXPTNZdataKeys = EXPTNZdata.keys()
           EXPEQdataKeys  = EXPEQdata.keys()

           PSINfig = plt.figure("PSIN")
           plt.plot(CHEASEdata['PSIN'],CHEASEdata['rhopsiN'],label=caselabel)
           plt.suptitle(shotnam[0][2:-6])
           plt.title('$\psi vs \\rho(\psi)$ ')
           plt.xlabel('$\psi$')
           plt.ylabel('$\\rho(\psi)$')
           plt.legend()
       
           EDENfig = plt.figure("Electron Density")
           plt.plot(CHEASEdata['PSIN'],CHEASEdata['ne'],label=caselabel)
           plt.suptitle(shotnam[0][2:-6])
           plt.title('Electron Density Profiles')
           plt.xlabel('$\psi$')
           plt.ylabel('$n_e$')
           plt.legend()
       
           dnedpsi = CubicSplineDerivative1D(x=CHEASEdata['PSI'],fx=CHEASEdata['ne'])

           GDNEfig = plt.figure("Electron Density Gradient")
           plt.plot(CHEASEdata['PSIN'],CHEASEdata['nePrime'],label=caselabel)
           plt.suptitle(shotnam[0][2:-6])
           plt.title('Electron Density Gradient Profiles')
           plt.xlabel('$\psi$')
           plt.ylabel('$\\nabla{n_e}$')
           plt.legend()
       
           ETMPfig = plt.figure("Electron Temperature")
           plt.plot(CHEASEdata['PSIN'],CHEASEdata['Te'],label=caselabel)
           plt.suptitle(shotnam[0][2:-6])
           plt.title('Electron Temperature Profiles')
           plt.xlabel('$\psi$')
           plt.ylabel('$T_e$')
           plt.legend()

           dTedpsi = CubicSplineDerivative1D(x=CHEASEdata['PSI'],fx=CHEASEdata['Te'])

           GDTEfig = plt.figure("Electron Temperature Gradient")
           plt.plot(CHEASEdata['PSIN'],CHEASEdata['TePrime'],label=caselabel)
           plt.suptitle(shotnam[0][2:-6])
           plt.title('Electron Temperature Gradient Profiles')
           plt.xlabel('$\psi$')
           plt.ylabel('$\\nabla{T_e}$')
           plt.legend()
       
           SFTYfig = plt.figure("Safety Factor (q)")
           plt.plot(CHEASEdata['PSIN'],CHEASEdata['q'],label=caselabel)
           plt.suptitle(shotnam[0][2:-6])
           plt.title("Safety Factor Profiles")
           plt.xlabel('$\psi$')
           plt.ylabel("q")
           plt.legend()

           TPPRfig = plt.figure("Plasma Pressure")
           plt.plot(CHEASEdata['PSIN'],CHEASEdata['PN'],label=caselabel)
           if 'pN' in EXPEQdataKeys:
              plt.plot(EXPEQdata['rhopsiN']**2,EXPEQdata['pN'],'--',label=caselabel)
           plt.suptitle(shotnam[0][2:-6])
           plt.title('Plasma Pressure Profiles')
           plt.xlabel('$\psi$')
           plt.ylabel('$P(\psi)$')
           plt.legend()
       
           PPRMfig = plt.figure("P'")
           plt.plot(CHEASEdata['PSIN'],CHEASEdata['PPrimeN'],label=caselabel)
           plt.suptitle(shotnam[0][2:-6])
           plt.title("P' Profiles")
           plt.xlabel('$\psi$')
           plt.ylabel("P'")
           plt.legend()
       
           TTPMfig = plt.figure("TT'")
           plt.plot(CHEASEdata['PSIN'],CHEASEdata['TTPrimeN'],label=caselabel)
           if 'ttprimeN' in EXPEQdataKeys:
              plt.plot(EXPEQdata['rhopsiN']**2,EXPEQdata['ttprimeN'],'--',label=caselabel)
           plt.suptitle(shotnam[0][2:-6])
           plt.title("TT' Profiles")
           plt.xlabel('$\psi$')
           plt.ylabel("TT'")
           plt.legend()
       
           ISTRfig = plt.figure("I*")
           plt.plot(CHEASEdata['PSIN'],CHEASEdata['ISTRN'],label=caselabel)
           if 'istrN' in EXPEQdataKeys:
              plt.plot(EXPEQdata['rhopsiN']**2,EXPEQdata['istrN'],'--',label=caselabel)
           plt.suptitle(shotnam[0][2:-6])
           plt.title("I* Profiles")
           plt.xlabel('$\psi$')
           plt.ylabel("I*")
           plt.legend()
       
           ICRTfig = plt.figure("Parallel Current")
           plt.plot(CHEASEdata['PSIN'],CHEASEdata['IPRLN'],label=caselabel)
           plt.suptitle(shotnam[0][2:-6])
           plt.title('Parallel Current Profiles')
           plt.xlabel('$\psi$')
           plt.ylabel('$I_{||}$')
           plt.legend()
       
           JCRTfig = plt.figure("Parallel Current Density")
           plt.plot(CHEASEdata['PSIN'],CHEASEdata['JPRLN'],label=caselabel)
           plt.suptitle(shotnam[0][2:-6])
           plt.title('Parallel Current Density Profiles')
           plt.xlabel('$\psi$')
           plt.ylabel('$J_{||}$')
           plt.legend()
       
           SCRTfig = plt.figure("Bootstrap Currents")
           plt.plot(CHEASEdata['PSIN'],CHEASEdata['JBSN'],label='$J_{BS}$-'+caselabel)
           plt.suptitle(shotnam[0][2:-6])
           plt.title('Bootstrap Current Density Profiles')
           plt.xlabel('$\psi$')
           plt.ylabel('$J_{BS}$')
           plt.legend()
       
           (CHEASEdata['PSIN2D'],CHEASEdata['CHIN2D']) = np.meshgrid(CHEASEdata['PSIN'],CHEASEdata['CHIN'])
           BF2Dfig = plt.figure("Magnetic Field, B($\psi$,$\chi$)")
           plt.contour(CHEASEdata['CHIN2D'],CHEASEdata['PSIN2D'],CHEASEdata['BN'],label=caselabel)
           plt.suptitle(shotnam[0][2:-6])
           plt.title('Magnetic Field Profiles')
           plt.xlabel('$\chi$')
           plt.ylabel('$\psi$')
       
           JPHIfig = plt.figure("Toroidal Current")
           plt.contour(CHEASEdata['CHIN2D'],CHEASEdata['PSIN2D'],CHEASEdata['JPHIN'],cmap=plt.cm.hot)
           plt.suptitle(shotnam[0][2:-6])
           plt.title('Toroidal Current Profiles')
           plt.xlabel('$\psi$')
           plt.ylabel('$J_{\phi}$')
       
           BFRZfig = plt.figure("Magnetic Field, B(R,Z}")
           plt.contour(CHEASEdata['RN'],CHEASEdata['ZN'],CHEASEdata['BN'],label=caselabel)
           plt.suptitle(shotnam[0][2:-6])
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
           chsfigs.savefig(PSINfig)
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

           CHEASEdata = read_cheaseh5(h5fid)
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

def mapping(f1,f2):
    from scipy.interpolate import CubicSpline
    fobj01 = CubicSpline(f1[:,0],f2,axis=0,bc_type='not-a-knot',extrapolate=None)
    plt.plot(fobj01(f1[:,0]),fobj01(f1[:,1]))
    plt.show()
    f3 = 1
    return f3

def profiles2exptnz(pfpath):
    profiles,units = efittools.read_profiles(pfpath.strip())
    EXPTNZdata = {}
    EXPTNZdata['rhopsiN'] = np.sqrt(profiles['psinorm'])
    EXPTNZdata['Te']      = profiles['te']
    EXPTNZdata['ne']      = profiles['ne']
    EXPTNZdata['Ti']      = profiles['ti']
    EXPTNZdata['ni']      = profiles['ni']

    nrhopsi               = np.size(EXPTNZdata['rhopsiN'])

    EXPTNZdata['Zeff']    = profiles['z'][1]**2*profiles['ni']
    EXPTNZdata['Zeff']   += profiles['z'][0]**2*profiles['nz1']
    EXPTNZdata['Zeff']   /= profiles['ne']
   #Zeff_array = np.array([1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0])
   #Zeff_mean  = np.mean(EXPTNZdata['Zeff'])
   #Zeff_diff  = abs(Zeff_array-Zeff_mean)
   #Zeff_value = Zeff_array[Zeff_diff==min(Zeff_diff)][0]
   #EXPTNZdata['Zeff']    = np.ones(nrhopsi)*Zeff_value

    ofh = open("%sEXPTNZ" % pfpath[-28:-8:1],'w')
    ofh.write('%5d rhopsi,  Te,   ne,   Zeff,   Ti,   ni  profiles\n' % nrhopsi)
    for i in range(nrhopsi): ofh.write('%16.6E\n' % EXPTNZdata['rhopsiN'][i])
    for i in range(nrhopsi): ofh.write('%16.6E\n' % EXPTNZdata['Te'][i])
    for i in range(nrhopsi): ofh.write('%16.6E\n' % EXPTNZdata['ne'][i])
    for i in range(nrhopsi): ofh.write('%16.6E\n' % EXPTNZdata['Zeff'][i])
    for i in range(nrhopsi): ofh.write('%16.6E\n' % EXPTNZdata['Ti'][i])
    for i in range(nrhopsi): ofh.write('%16.6E\n' % EXPTNZdata['ni'][i])
    ofh.close()

    return EXPTNZdata

def iterdb2exptnz(iterdbfpath,eqdskfpath):
    rhotor,profiles,units = read_iterdb.read_iterdb(iterdbfpath)
    eqdskdata = efittools.read_eqdsk(eqdskfpath)
    EXPTNZdata = {}

    nerhotorfn            = interp1d(rhotor['NE'], profiles['NE'], kind='linear')
    terhotorfn            = interp1d(rhotor['TE'], profiles['TE'], kind='linear')
    nirhotorfn            = interp1d(rhotor['NM1'],profiles['NM1'],kind='linear')
    tirhotorfn            = interp1d(rhotor['TI'], profiles['TI'], kind='linear')
    nzrhotorfn            = interp1d(rhotor['NM2'],profiles['NM2'],kind='linear')
    EXPTNZdata['ne']      = nerhotorfn(eqdskdata['rhotor'])
    EXPTNZdata['Te']      = terhotorfn(eqdskdata['rhotor'])
    EXPTNZdata['ni']      = nirhotorfn(eqdskdata['rhotor'])
    EXPTNZdata['Ti']      = tirhotorfn(eqdskdata['rhotor'])
    EXPTNZdata['nz']      = nzrhotorfn(eqdskdata['rhotor'])

    EXPTNZdata['rhopsiN'] = eqdskdata['rhopsi']
    nrhopsi               = np.size(EXPTNZdata['rhopsiN'])

    EXPTNZdata['Zeff']    = (2**2)*EXPTNZdata['ni']
    EXPTNZdata['Zeff']   += (6**2)*EXPTNZdata['nz']
    EXPTNZdata['Zeff']   /= EXPTNZdata['ne']

    ofh = open("%s_EXPTNZ" % iterdbfpath[:-7],'w')
    ofh.write('%5d rhopsi,  Te,   ne,   Zeff,   Ti,   ni  profiles\n' % nrhopsi)
    for i in range(nrhopsi): ofh.write('%16.6E\n' % EXPTNZdata['rhopsiN'][i])
    for i in range(nrhopsi): ofh.write('%16.6E\n' % EXPTNZdata['Te'][i])
    for i in range(nrhopsi): ofh.write('%16.6E\n' % EXPTNZdata['ne'][i])
    for i in range(nrhopsi): ofh.write('%16.6E\n' % EXPTNZdata['Zeff'][i])
    for i in range(nrhopsi): ofh.write('%16.6E\n' % EXPTNZdata['Ti'][i])
    for i in range(nrhopsi): ofh.write('%16.6E\n' % EXPTNZdata['ni'][i])
    ofh.close()

    return EXPTNZdata


def eqdsk2expeq(eqdskdata={},eqdskfpath='',setParam={}):
    if not eqdskdata:
       if eqdskfpath:
          eqdskdata = efittools.read_eqdsk(eqdskfpath)
       else:
          print('FATAL: EQDSK FILE NOT FOUND. EXIT!')
          sys.exit()

    mu0 = 4.0e-7*np.pi

    print('Converting EQDSK to EXPEQ ...')
    expeqdata  = {}

    if 'nsttp' in setParam: expeqdata['nsttp'] = setParam['nsttp']
    else:                   expeqdata['nsttp'] = 1

    if 'nrhotype' in setParam: expeqdata['nrhotype'] = setParam['nrhotype']
    else:                      expeqdata['nrhotype'] = 0

    if 'nppfun' in setParam: expeqdata['nppfun'] = setParam['nppfun']
    else:                    expeqdata['nppfun'] = 4

    rbound,zbound = efittools.magsurf_solvflines(eqdskdata=eqdskdata,psi=0.9999,eps=1.0e-16)

    expeqdata['R0EXP']   = eqdskdata['RCTR']
    expeqdata['B0EXP']   = eqdskdata['BCTR']
    expeqdata['zgeom']   = eqdskdata['ZMAX']/eqdskdata['RCTR']
    expeqdata['pedge']   = eqdskdata['pressure'][-1]
    expeqdata['aspect']  = (max(eqdskdata['rbound'])-min(eqdskdata['rbound']))
    expeqdata['aspect'] /= (max(eqdskdata['rbound'])+min(eqdskdata['rbound']))
    expeqdata['nbound']  = np.size(rbound)

    ofh = open("EXPEQ",'w')
    ofh.write('%18.8E\n'           % expeqdata['aspect'])
    ofh.write('%18.8E\n'           % expeqdata['zgeom'])
    ofh.write('%18.8E\n'           % expeqdata['pedge'])
    ofh.write('%5d\n'              % expeqdata['nbound'])
    expeqdata['rbound'] = rbound/eqdskdata['RCTR']
    expeqdata['zbound'] = zbound/eqdskdata['RCTR']
    for i in range(expeqdata['nbound']):
        ofh.write('%18.8E%18.8E\n' % (expeqdata['rbound'][i],expeqdata['zbound'][i]))
    ofh.write('%5d%5d\n'           % (np.size(eqdskdata['rhopsi']),expeqdata['nppfun']))
    ofh.write('%5d%5d\n'           % (expeqdata['nsttp'],expeqdata['nrhotype']))

    if expeqdata['nrhotype'] == 0 or expeqdata['nrhotype'] == 'rhopsi':
       expeqdata['rhopsi'] = eqdskdata['rhopsi']
       for i in range(np.size(eqdskdata['rhopsi'])):
           ofh.write('%18.8E\n'       % expeqdata['rhopsi'][i])
    elif expeqdata['nrhotype'] == 1 or expeqdata['nrhotype'] == 'rhotor':
       expeqdata['rhotor'] = eqdskdata['rhotor']
       for i in range(np.size(eqdskdata['rhotor'])):
           ofh.write('%18.8E\n'       % eqdskdata['rhotor'][i])
    else:
	    raise ValueError('FATAL: nrhotype takes only 0 and 1. Exit!')

    if expeqdata['nppfun'] == 4 or expeqdata['nppfun'][0] == 'pprime':
       expeqdata['PPrime'] = mu0*eqdskdata['pprime']*eqdskdata['RCTR']**2/eqdskdata['BCTR']
       for i in range(np.size(expeqdata['PPrime'])):
           ofh.write('%18.8E\n'       % expeqdata['PPrime'][i])
    elif expeqdata['nppfun'] == 8 or expeqdata['nppfun'] == 'p':
       expeqdata['p'] = mu0*eqdskdata['pressure']/eqdskdata['BCTR']**2
       for i in range(np.size(expeqdata['p'])):
           ofh.write('%18.8E\n'       % expeqdata['p'][i])
    else:
       raise ValueError('FATAL: nppfun takes only 4 and 8. Exit!')

    if expeqdata['nsttp'] == 1 or expeqdata['nsttp'] == 'ttprime':
       expeqdata['TTPrime'] = eqdskdata['ffprime']/eqdskdata['BCTR']
       if expeqdata['nrhotype'] == 1 or expeqdata['nrhotype'] == 'rhotor':
          raise ValueError('FATAL: with ttprime or ffprime, nrhotype (must) = 0 or rhopsi. EXIT!')
       for i in range(np.size(expeqdata['TTPrime'])):
           ofh.write('%18.8E\n'       % expeqdata['TTPrime'][i])
    elif expeqdata['nsttp'] == 5 or expeqdata['nsttp'] == 'q':
       expeqdata['q'] = eqdskdata['qpsi']
       for i in range(np.size(eqdskdata['qpsi'])):
           ofh.write('%18.8E\n'       % expeqdata['qpsi'][i])
    else:
       raise ValueError('FATAL: nsttp takes only 1 and 5. Exit!')

    ofh.close()

    return expeqdata


def cheasepy(srcVals={},namelistVals={},pltVals={}):
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
              selection = input('Selected Option: ')
              if    selection in [1,2,3,4]: break
              else: raise(NameError)
          except NameError:
             print(CRED+'Select between 1,2,3, and 4 options.'+CEND)
             continue

    # "iterTotal     = n" where n = 0, 1, 2, ...
    # "rhopsi_src    = n" where n = 0 or 'h5', n = 1 or 'expeq', n = 1 or 'exptnz'
    # "current_src   = n" where n = 0 or 'h5', n = 1 or 'expeq'
    # "pressure_src  = n" where n = 0 or 'h5', n = 1 or 'expeq', n = 1 or 'exptnz'
    # "eprofiles_src = n" where n = 0 or 'h5', n = 1 or 'exptnz'
    # "iprofiles_src = n" where n = 0 or 'h5', n = 1 or 'exptnz'
    srcValsKeys     = srcVals.keys()
    if 'iterTotal'     in srcValsKeys: iterTotal     = srcVals['iterTotal']     
    else:                              iterTotal     = 0  
    if 'rhopsi_src'    in srcValsKeys: rhopsi_src    = srcVals['rhopsi_src']   
    else:                              rhopsi_src    = 0  
    if 'current_src'   in srcValsKeys: current_src   = srcVals['current_src']   
    else:                              current_src   = 0  
    if 'pressure_src'  in srcValsKeys: pressure_src  = srcVals['pressure_src']   
    else:                              pressure_src  = 0  
    if 'eprofiles_src' in srcValsKeys: eprofiles_src = srcVals['eprofiles_src']   
    else:                              eprofiles_src = 0  
    if 'iprofiles_src' in srcValsKeys: iprofiles_src = srcVals['iprofiles_src']   
    else:                              iprofiles_src = 0  
    
    if selection == 1:
       while True:
             print(CYLW+'Select CHEASE running mode:'+CEND)
             print(CYLW+'(1) Check Equilibrium Preservation Over Multiple Iterations.'+CEND)
             print(CYLW+'(2) Converge to Total Current by correcting Ohmic current.'+CEND)
             print(CYLW+'(3) Converge to Total Current by correcting Bootstrap current.'+CEND)
             try:
                cheasemode = input('Selected Option: ')
                if    cheasemode in [1,2,3]: break
                else: raise(NameError)
             except NameError:
                print(CRED+'Select between 1, 2, and 3 options.'+CEND)
                continue
    
       namelistParam           = {}
       namelistParam['CSSPEC'] = 0
       namelistParam['NCSCAL'] = 4
       print(CRED+'List of files from previous run(s):'+CEND)
       os.system('ls')
       if raw_input(CBLU+'Remove output (h5, pdf, dat, OUT) files of previous runs (yes/no)? '+CEND).lower() in ['yes','y']:
          if glob('./NGA'):                   os.system('rm NGA')
          if glob('./NDES'):                  os.system('rm NDES')
          if glob('./*OUT*'):                 os.system('rm *OUT*')
          if glob('./*.pdf'):                 os.system('rm *.pdf')
          if glob('./EXPEQ'):                 os.system('rm EXPEQ')
          if glob('./EXPTNZ'):                os.system('rm EXPTNZ')
          if glob('./*.iterdb'):              os.system('rm *.iterdb')
          if glob('./ogyropsi*'):             os.system('rm ogyropsi*')
          if glob('./EXPEQ_EQDSK*'):          os.system('rm EXPEQ_EQDSK*')
          if glob('./EXPEQ_EXPEQ*'):          os.system('rm EXPEQ_EXPEQ*')
          if glob('./chease_namelist*'):      os.system('rm chease_namelist*')
          print(CRED+'List of Available CHEASE Files:'+CEND)
          os.system('ls')
    
       if len(glob('./*_EQDSK'))==0 or len(glob('./chease_parameters.csv'))==0:
          shotlist = sorted(glob('./shots/*'))
          while True:
                print(CYLW+'List of the available shots:'+CEND)
                for ishot in range(len(shotlist)):
                    print(CYLW+'(%02d) %s' % (ishot+1,shotlist[ishot][8:])+CEND)
                try:
                   shotrec = input('Select Shot Number: ')
                   if shotrec-1 in range(len(shotlist)):
                      print(CGRN+'Chease runs the %s shot.' % shotlist[shotrec-1][8:]+CEND)
                      break
                   else:
                      raise(NameError)
                except NameError:
                   print(CRED+'Choose ONLY from the %0d available shots.' % len(shotlist)+CEND)
                   continue
          if os.path.isfile('%s/chease_parameters.csv' % (shotlist[shotrec-1])):
             os.system('cp %s/chease_parameters.csv .' % shotlist[shotrec-1])
             namelist = namelistcreate('chease_parameters.csv',0,namelistParam)
          else:
             raise IOError('chease namelist file NOT FOUND in the given path!')
          print('%s/%s_EQDSK'   % (shotlist[shotrec-1],shotlist[shotrec-1][8:]))
          if os.path.isfile('%s/%s_EQDSK'   % (shotlist[shotrec-1],shotlist[shotrec-1][8:])):
             os.system('cp %s/%s_EQDSK .'   % (shotlist[shotrec-1],shotlist[shotrec-1][8:]))
          else:
             raise IOError('EQDSK file NOT FOUND in the given path!')
          if int(namelist['NBSEXPQ'][0]) != 0:
             if   os.path.isfile('%s/%s_EXPTNZ'      % (shotlist[shotrec-1],shotlist[shotrec-1][8:])):
                  os.system('cp %s/%s_EXPTNZ .'      % (shotlist[shotrec-1],shotlist[shotrec-1][8:]))
             elif os.path.isfile('%s/%s.iterdb'      % (shotlist[shotrec-1],shotlist[shotrec-1][8:])):
                  os.system('cp %s/%s.iterdb .'      % (shotlist[shotrec-1],shotlist[shotrec-1][8:]))
                  iterdb2exptnz('%s.iterdb' % shotlist[shotrec-1][8:],'%s_EQDSK' % shotlist[shotrec-1][8:])
             elif os.path.isfile('%s/%s_Profiles'    % (shotlist[shotrec-1],shotlist[shotrec-1][8:])):
                  os.system('cp %s/%s_Profiles .'    % (shotlist[shotrec-1],shotlist[shotrec-1][8:]))
                  profiles2exptnz('%s_Profiles'      %                     (shotlist[shotrec-1][8:]))
             else:
                  raise IOError('Profiles (EXPTNZ or Profiles) files NOT FOUND in the given path!')
          os.system('ls')
          if raw_input(CBLU+'Do you want to continue? (yes/no)? '+CEND).lower() not in ['yes','y']: sys.exit()
       elif raw_input(CRED+'Do you want to use the avaiable shot EXPEQ and EXPTNZ? (yes/no)? '+CEND).lower() in ['no','n']:
          if glob('./*_EQDSK'):               os.system('rm *_EQDSK')
          if glob('./*_EXPTNZ'):              os.system('rm *_EXPTNZ')
          if glob('./*_Profiles'):            os.system('rm *_Profiles')
          if glob('./chease_parameters.csv'): os.system('rm chease_parameters.csv')
          shotlist = sorted(glob('./shots/*'))
          while True:
                print(CYLW+'List of the available shots:'+CEND)
                for ishot in range(len(shotlist)):
                    print(CYLW+'(%02d) %s' % (ishot+1,shotlist[ishot][8:])+CEND)
                try:
                   shotrec = input('Select Shot Number: ')
                   if shotrec-1 in range(len(shotlist)):
                      print(CGRN+'Chease runs the %s shot.' % shotlist[shotrec-1][8:]+CEND)
                      break
                   else:
                      raise(NameError)
                except NameError:
                   print(CRED+'Choose ONLY from the %0d available shots.' % len(shotlist)+CEND)
                   continue
          if   os.path.isfile('%s/chease_parameters.csv'   % shotlist[shotrec-1]):
               os.system('cp   %s/chease_parameters.csv .' % shotlist[shotrec-1])
               namelist = namelistcreate('chease_parameters.csv',0)
          else:
               raise IOError('chease_parameters.csv file NOT FOUND in the given path!')
          print('%s/%s_EQDSK'   % (shotlist[shotrec-1],shotlist[shotrec-1][8:]))
          if os.path.isfile('%s/%s_EQDSK'       %(shotlist[shotrec-1],shotlist[shotrec-1][8:])):
             os.system('cp   %s/%s_EQDSK .'     %(shotlist[shotrec-1],shotlist[shotrec-1][8:]))
          else:
             raise IOError('EQDSK (or KEFIT) file NOT FOUND in the given path!')
          if int(namelist['NBSEXPQ'][0]) != 0:
             if   os.path.isfile('%s/%s_EXPTNZ'      % (shotlist[shotrec-1],shotlist[shotrec-1][8:])):
                  os.system('cp   %s/%s_EXPTNZ .'    % (shotlist[shotrec-1],shotlist[shotrec-1][8:]))
             elif os.path.isfile('%s/%s.iterdb'      % (shotlist[shotrec-1],shotlist[shotrec-1][8:])):
                  os.system('cp %s/%s.iterdb .'      % (shotlist[shotrec-1],shotlist[shotrec-1][8:]))
                  iterdb2exptnz('%s.iterdb' % shotlist[shotrec-1][8:],'%s_EQDSK' % shotlist[shotrec-1][8:])
             elif os.path.isfile('%s/%s_Profiles'    % (shotlist[shotrec-1],shotlist[shotrec-1][8:])):
                  os.system('cp   %s/%s_Profiles .'  % (shotlist[shotrec-1],shotlist[shotrec-1][8:]))
                  profiles2exptnz('  %s_Profiles'    % (                    shotlist[shotrec-1][8:]))
             else:
                 raise IOError('Profiles (EXPTNZ or Profiles) files NOT FOUND in the given path!')
          os.system('ls')
          if raw_input(CBLU+'Do you want to continue? (yes/no)? '+CEND).lower() not in ['yes','y']: sys.exit()
       else:
          namelist = namelistcreate('chease_parameters.csv',0,namelistParam)

       EQDSKfname = glob('./*_EQDSK')
       if   int(namelist['NEQDSK'][0]) == 1:
            print('Reading from EQDSK file.')
            os.system('cp *_EQDSK  EXPEQ')
       elif int(namelist['NEQDSK'][0]) == 0:
            print('Reading from EXPEQ file.')
            if len(glob('./EXPEQ'))==0:
               expeqdata = eqdsk2expeq(eqdskfpath=EQDSKfname[0])
               namelistParam['R0EXP'] = expeqdata['R0EXP']
               namelistParam['B0EXP'] = expeqdata['B0EXP']
            else:
               eqdskdata=efittools.read_eqdsk(EQDSKfname[0])
               namelistParam['R0EXP'] = eqdskdata['RCTR']
               namelistParam['B0EXP'] = eqdskdata['BCTR']
            namelist = namelistcreate('chease_parameters.csv',0,namelistParam)

       if len(glob('./EXPTNZ'))==0:
          os.system('cp *_EXPTNZ EXPTNZ')

       ofh = open(EQDSKfname[0],'r')
       ITEXP = float(ofh.readlines()[3][0:16])
       ofh.close()

       expeqParam               = {}
       exptnzParam              = {}
    
       expeqParam['ITEXP']      = ITEXP
       expeqParam['cheasemode'] = cheasemode
       
       #Define an attribute for write_expeq function to improve the convergence to the total current
       #by tracking its error value.
       write_expeq.errorval = float('inf')

       OGYROPSIfname = sorted(glob('./ogyropsi_*.h5'))
       if   len(OGYROPSIfname)==0:
            it=0
       else:
            it=int(OGYROPSIfname[-1][-6:-3])+1
       exit_status = os.system('./chease_hdf5 chease_namelist > iter%03d.OUT' % it)
       if abs(exit_status) > 0: sys.exit()
       if os.path.isfile('./chease_namelist'): os.system('mv ./chease_namelist ./chease_namelist_iter%03d' % it)
       if os.path.isfile('./ogyropsi.dat'): os.system('mv ./ogyropsi.dat ogyropsi_iter%03d.dat' % it)
       if os.path.isfile('./ogyropsi.h5'): os.system('mv ./ogyropsi.h5 ogyropsi_iter%03d.h5' % it)
       if os.path.isfile('./EXPEQ.OUT'): os.system('mv ./EXPEQ.OUT EXPEQ_iter%03d.OUT' % it)
       if os.path.isfile('./EXPTNZ.OUT'): os.system('mv ./EXPTNZ.OUT EXPTNZ_iter%03d.OUT' % it)
       if os.path.isfile('./EXPEQ_EXPEQ.IN'): os.system('mv ./EXPEQ_EXPEQ.IN EXPEQ_EXPEQ_iter%03d.IN' % it)
    
       cheasedata = read_cheaseh5('ogyropsi_iter%03d.h5' % it)
       #expeqParam['nsttp']     = [NSTTP_Namelist,   Data_Source]  
       #expeqParam['nppfun']    = [NPPFUN_Namelist,  Data_Source]
       #expeqParam['nrhomesh']  = [NRHOMESH_Namelist,Data_Source]
       #Data_Source = 0 (ogyropsi_iterxxx.h5), 1 (EXPEQ_iterxxx.OUT), 2 (EXPTNZ_iterxxx.OUT)
       expeqParam['nsttp']      = [int(namelist['NSTTP'][0]),current_src]
       expeqParam['nppfun']     = [int(namelist['NPPFUN'][0]),pressure_src]
       expeqParam['nrhotype']   = [int(namelist['NRHOMESH'][0]),rhopsi_src]
       expeqParam['profiles']   = [eprofiles_src,iprofiles_src]
       expeqflag                = write_expeq('ogyropsi_iter%03d.h5' % it,'EXPEQ_iter%03d.OUT' % it,'EXPTNZ_iter%03d.OUT' % it,setParam=expeqParam)
       #exptnzParam['nrhomesh'] = [NRHOMESH_Namelist,Data_Source]
       #Data_Source = 0 (ogyropsi_iterxxx.h5), 1 (EXPTNZ_iterxxx.OUT)
       exptnzParam['nrhotype']  = [int(namelist['NRHOMESH'][0]),rhopsi_src]
       exptnzParam['profiles']  = [eprofiles_src,iprofiles_src]
       exptnzflag               = write_exptnz('ogyropsi_iter%03d.h5' % it,'EXPTNZ_iter%03d.OUT' % it,setParam=exptnzParam)
    
       #namelistParam['RBOXLFT'] = cheasedata['rmesh'][0]
       #namelistParam['RBOXLEN'] = cheasedata['rmesh'][-1]-cheasedata['rmesh'][0]
       #namelistParam['ZBOXLEN'] = cheasedata['zmesh'][-1]*2.0
       namelistParam['R0EXP']   = cheasedata['R0EXP']
       namelistParam['B0EXP']   = cheasedata['B0EXP']
       namelistParam['QSPEC']   = cheasedata['q'][0]
    
       ITErr = (cheasedata['ITOR']-ITEXP)/ITEXP
       print 'Iter  = ', 0
       print 'ITOR  = ', cheasedata['ITOR']
       print 'ITEXP = ', ITEXP
       print 'ITErr = ', abs(ITErr)
    
       it+=1
       ITErr = (cheasedata['ITOR']-ITEXP)/ITEXP
       while (abs(ITErr) > 1.0e-6):
           if (cheasemode != 2) and (it >= iterTotal+1): break
           namelist = namelistcreate('chease_parameters.csv',min(len(namelist['fname'])-1,it),namelistParam)
           os.system('cp chease_namelist chease_namelist_iter%3d' % (min(len(namelist['fname'])-1,it)))
           exit_status = os.system('./chease_hdf5 chease_namelist > iter%03d.OUT' % it)
           if abs(exit_status) > 0: sys.exit()
           if os.path.isfile('./ogyropsi.dat'): os.system('mv ./ogyropsi.dat ogyropsi_iter%03d.dat' % it)
           if os.path.isfile('./ogyropsi.h5'): os.system('mv ./ogyropsi.h5 ogyropsi_iter%03d.h5' % it)
           if os.path.isfile('./EXPEQ.OUT'): os.system('mv ./EXPEQ.OUT EXPEQ_iter%03d.OUT' % it)
           if os.path.isfile('./EXPTNZ.OUT'): os.system('mv ./EXPTNZ.OUT EXPTNZ_iter%03d.OUT' % it)
           if os.path.isfile('./EXPEQ_EXPEQ.IN'): os.system('mv ./EXPEQ_EXPEQ.IN EXPEQ_EXPEQ_iter%03d.IN' % it)
    
           cheasedata = read_cheaseh5('ogyropsi_iter%03d.h5' % it)
           expeqParam['nsttp']      = [int(namelist['NSTTP'][min(len(namelist['fname'])-1,it)]),current_src]
           expeqParam['nppfun']     = [int(namelist['NPPFUN'][min(len(namelist['fname'])-1,it)]),pressure_src]
           expeqParam['nrhotype']   = [int(namelist['NRHOMESH'][min(len(namelist['fname'])-1,it)]),rhopsi_src]
           expeqflag                = write_expeq('ogyropsi_iter%03d.h5' % it,'EXPEQ_iter%03d.OUT' % it,'EXPTNZ_iter%03d.OUT' % it,setParam=expeqParam)
           exptnzParam['nrhotype']  = [int(namelist['NRHOMESH'][0]),rhopsi_src]
           exptnzParam['profiles']  = [eprofiles_src,iprofiles_src]
           exptnzflag               = write_exptnz('ogyropsi_iter%03d.h5' % it,'EXPTNZ_iter%03d.OUT' % it,setParam=exptnzParam)
    
           ITErr = (cheasedata['ITOR']-ITEXP)/ITEXP
           print 'Iter  = ', it
           print 'ITOR = ', cheasedata['ITOR']
           print 'ITEXP = ', ITEXP
           print 'ITErr = ', abs(ITErr)
    
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
       plot_cheasedata('./',skipfigs=0)
    elif selection == 2:
       plot_cheasedata('./',skipfigs=0)
    elif selection == 3:
       if glob('./NGA'):                   os.system('rm NGA')
       if glob('./NDES'):                  os.system('rm NDES')
       if glob('./*OUT*'):                 os.system('rm *OUT*')
       if glob('./*.pdf'):                 os.system('rm *.pdf')
       if glob('./EXPEQ'):                 os.system('rm EXPEQ')
       if glob('./EXPTNZ'):                os.system('rm EXPTNZ')
       if glob('./*_EQDSK'):               os.system('rm *_EQDSK')
       if glob('./*_EXPTNZ'):              os.system('rm *_EXPTNZ')
       if glob('./*.iterdb'):              os.system('rm *.iterdb')
       if glob('./ogyropsi*'):             os.system('rm ogyropsi*')
       if glob('./*_Profiles'):            os.system('rm *_Profiles')
       if glob('./EXPEQ_EQDSK*'):          os.system('rm EXPEQ_EQDSK*')
       if glob('./EXPEQ_EXPEQ*'):          os.system('rm EXPEQ_EXPEQ*')
       if glob('./chease_namelist*'):      os.system('rm chease_namelist*')
       if glob('./chease_parameters.csv'): os.system('rm chease_parameters.csv')
    
    return 1
