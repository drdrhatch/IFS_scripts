#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import cheasepy
import argparse

parser = argparse.ArgumentParser(description='Running CHEASEPY for MHD Equlibrium Code.')
parser.add_argument('--manual','-m',action='store_const',const=1,help='Run CHEASEPY and Select Input Parameters Manually')
parser.add_argument('--submit','-s',action='store_const',const=1,help='Run CHEASEPY with Predefined Input Parameters')

if   parser.parse_args():
     args = parser.parse_args()
     manual_chease = args.manual 
     submit_chease = args.submit 
else:
     manual_chease = 1 
     submit_chease = 0

'''
Source Options:
n = 0 and chease
n = 1 and eqdsk
n = 2 and expeq
n = 3 and exptnz
n = 4 and profiles
n = 5 and iterdb
'''
srcVals = {}
#Total Number of Iterations of the CHEASE Code
srcVals['iterTotal']     = 3
#Source of Electron Profile:
srcVals['eprofiles_src'] = ['profiles','chease']
#Source of Ion Profile:
srcVals['iprofiles_src'] = ['profiles','chease']
#Source of Pressure or Pressure Gradient:
srcVals['pressure_src']  = ['eqdsk','chease']
#Source of Poloidal or Toroidal Magnetic Flux Surface
srcVals['rhomesh_src']   = ['eqdsk','chease']
#Source of Current Flux Density
srcVals['current_src']   = ['eqdsk','chease']

'''
Source of Boundary Source
n = 0 or 'asis'
n = 1 or 'interp'
'''
srcVals['boundary_type'] = 'asis'

if   submit_chease:
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
     cheaseVals['shotpath'] = './shots/DIIID_KEFITD_162940'
     
     cheasepy.cheasepy(srcVals=srcVals,cheaseVals=cheaseVals)
else:
     cheasepy.cheasepy(srcVals=srcVals)

sys.exit()
