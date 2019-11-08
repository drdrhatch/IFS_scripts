#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy
import cheasepy


importedVals = {}

eqdskpath  = 'shots/JET_KEFITD_78697/PROFILES/jet78697.51005_hager.eqdsk'
eqdskdata  = cheasepy.read_eqdsk(eqdskpath)
qin   = eqdskdata['q']
rhoin = eqdskdata['rhopsi']

eqdskpath  = 'shots/DIIID_KEFITD_162940/DIIID_KEFITD_162940_EQDSK'
eqdskdata  = cheasepy.read_eqdsk(eqdskpath)

iterdbpath1p0 = 'shots/JET_KEFITD_78697/PROFILES/jet78697.51005_hager_Z6.0Zeff2.35_negom_alpha0.7_omti_x0_0.96.iterdb'
iterdbdata1p0 = cheasepy.read_iterdb(iterdbfpath=iterdbpath1p0,eqdsk=eqdskpath,setParam={'nrhomesh':0})

Pdiff = eqdskdata['pressure']-iterdbdata1p0['pressure']

iterdbpath1p3 = 'shots/JET_KEFITD_78697/PROFILES/jet78697.51005_hager_Z6.0Zeff2.35_negom_alpha0.7_omti_x0_0.96_alpha1.3_omte_x0_0.97.iterdb'
iterdbdata1p3 = cheasepy.read_iterdb(iterdbfpath=iterdbpath1p3,eqdsk=eqdskpath,setParam={'nrhomesh':0})

importedVals['rhotor']   = eqdskdata['rhotor']
importedVals['rhopsi']   = eqdskdata['rhopsi']
importedVals['q']        = numpy.interp(eqdskdata['rhopsi'],rhoin,qin)
importedVals['pressure'] = iterdbdata1p3['pressure']+Pdiff


'''
Source Options:
n = 0 and chease
n = 1 and eqdsk
n = 2 and expeq
n = 3 and exptnz
n = 4 and profiles
n = 5 and iterdb
n = 7 and imported
'''
srcVals = {}
#Total Number of Iterations of the CHEASE Code
srcVals['iterTotal']     = 0
#Source of Electron Profile:
if srcVals['iterTotal']==0: srcVals['eprofiles_src'] =  'iterdb'
else:                       srcVals['eprofiles_src'] = ['iterdb','chease']
#Source of Ion Profile:
if srcVals['iterTotal']==0: srcVals['iprofiles_src'] =  'iterdb'
else:                       srcVals['iprofiles_src'] = ['iterdb','chease']
#Source of Pressure or Pressure Gradient:
if srcVals['iterTotal']==0: srcVals['pressure_src']  =  'iterdb'
else:                       srcVals['pressure_src']  = ['iterdb','chease']
#Source of Poloidal or Toroidal Magnetic Flux Surface
if srcVals['iterTotal']==0: srcVals['rhomesh_src']   =  'eqdsk'
else:                       srcVals['rhomesh_src']   = ['eqdsk','chease']
#Source of Current Flux Density
if srcVals['iterTotal']==0: srcVals['current_src']   =  'expeq'
else:                       srcVals['current_src']   = ['expeq','chease']

'''
Source of Boundary Source
n = 0 or 'asis'
n = 1 or 'interp'
'''
srcVals['boundary_type'] = 'asis'
#srcVals['boundary_type'] = 'interp'


namelistVals = {}
if srcVals['iterTotal']==0:
   namelistVals['NS']        = 256
   namelistVals['NT']        = 256
   namelistVals['NISO']      = 256
   namelistVals['NPSI']      = 1024
   namelistVals['NCHI']      = 1024
   namelistVals['NRBOX']     = 60
   namelistVals['NZBOX']     = 60
   namelistVals['RELAX']     = 0.0
   namelistVals['NSTTP']     = 1
   namelistVals['NPROPT']    = 3
   namelistVals['NPPFUN']    = 8
   namelistVals['NEQDSK']    = 0
   namelistVals['NFUNRHO']   = 0
   namelistVals['TENSBND']   = 0
   namelistVals['TENSPROF']  = 0
   namelistVals['NRHOMESH']  = 0
   namelistVals['cocos_in']  = 2
   namelistVals['cocos_out'] = 12
else:
   namelistVals['NS']        = [256,256]
   namelistVals['NT']        = [256,256]
   namelistVals['NISO']      = [256,256]
   namelistVals['NPSI']      = [1024,1024]
   namelistVals['NCHI']      = [1024,1024]
   namelistVals['NRBOX']     = [60,60]
   namelistVals['NZBOX']     = [60,60]
   namelistVals['RELAX']     = [0.0,0.0]
   namelistVals['NSTTP']     = [1,3]
   namelistVals['NPROPT']    = [3,3]
   namelistVals['NPPFUN']    = [8,8]
   namelistVals['NEQDSK']    = [0,0]
   namelistVals['NFUNRHO']   = [0,0]
   namelistVals['TENSBND']   = [0,0]
   namelistVals['TENSPROF']  = [0,0]
   namelistVals['NRHOMESH']  = [0,0]
   namelistVals['cocos_in']  = [2,12]
   namelistVals['cocos_out'] = [12,12]

cheaseVals = {}
#Select on of the following options:
#   (1) Run the Code and Plot the Outputs.
#   (2) Plot the Available Outputs.
#   (3) Remove Input/Output Files.
cheaseVals['runmode'] = 1
#Select CHEASE running mode:
#   (1) Check Equilibrium Preservation Over Multiple Iterations.
#   (2) Converge to Total Current by correcting Current.
#   (3) Converge to Total Current by correcting Pressure.
cheaseVals['cheasemode'] = 1
#Remove input (EQDSK, EXPEQ, and EXPTNZ) files from previous run (yes/no)?
cheaseVals['removeinputs'] = 'yes'
#Remove output (h5, pdf, dat, OUT) files from previous run (yes/no)?
cheaseVals['removeoutputs'] = 'yes'
#Path to the selected shot.
#cheaseVals['shotpath'] = './shots/DIIID_KEFITD_162940'
cheaseVals['shotpath'] = './shots/JET_KEFITD_78697'
#cheaseVals['shotpath'] = './shots/CASE3_20190905'
#cheaseVals['shotpath'] = './shots/CASE2_20190509'
 
cheasepy.cheasepy(srcVals=srcVals,namelistVals=namelistVals,cheaseVals=cheaseVals,importedVals=importedVals)

sys.exit()
