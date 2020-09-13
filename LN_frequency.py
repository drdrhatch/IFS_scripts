#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Created by Max T. Curie 09/13/2020

import optparse as op
import sys

from ParIO import Parameters 
from LN_tools import start_end_time
from LN_tools import LN_apar_frequency_nz




pic_path='pic'        #path one want to store the picture and video in
csv_path='csv'        #path one want to store the picture and video in


parser=op.OptionParser(description='Some infrastructure for reading in, manipulating, and plotting nonlinear field data.')
#parser.add_option('--plot_theta','-g',action='store_const',const=False,help = 'Plot global mode structures decomposed in poloidal m number.',default=True)
options,args=parser.parse_args()
print("options",options)
print("args",args)
if len(args)!=1:
    exit("""
Please include run number as argument (e.g., 0001)."
    \n""")
suffix = args[0]

if suffix in ['dat','.dat']:
     suffix = '.dat'
else:
     suffix = '_'+suffix

par = Parameters()
par.Read_Pars('parameters'+suffix)
pars = par.pardict

inz=int(pars['nz0']/2.)

time_start,time_end=start_end_time(suffix,pars)
frequency_kHZ,amplitude_frequency_sum,amplitude_growth_sum=LN_apar_frequency_nz(suffix,inz,time_start,time_end,plot=True,pic_path='pic',csv_path='csv',output_csv=True,show=True)