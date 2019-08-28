#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from omega_tool import omega_calc
import os
import optparse as op
from subprocess import call

parser=op.OptionParser(description='Calculates various quasilinear estimates of heat flux and outputs plots and summary files.')
parser.add_option('--archive','-a', action='store_true',dest = 'archive', help = 'Move to home directory.', default=True)
options,args=parser.parse_args()
if len(args)!=2:
    exit("""
Please include nonlinear run number (i.e. 1) and scanfiles number (i.e., 0001)
    \n""")
suffix = args[0]
sfsuffix = args[1]
archive = options.archive

arch_dir = '/marconi/home/userexternal/dhatch00/ETGQL/'

if 'dat' in suffix:
   suffix = '.dat'

thisdir=os.getcwd().split('/')[-1]
outdir = 'ETGQL_'+thisdir
if not os.path.isfile(outdir):
    call(['mkdir',outdir])

#s1 add
add_all = True
include_qn2 = False
sat_rule = 1
most_unstable = False 
file_name = 'QL_summary_'+suffix+'_'+sfsuffix+'_iqn2'+str(include_qn2)[0]+'_sat'+str(sat_rule)+'_add'+str(add_all)[0]+'_mu'+str(most_unstable)[0]
call(['ETG_quasilinear.py',suffix,sfsuffix,'-n','-a'])
call(['mv',file_name,outdir])
call(['mv',file_name+'.ps',outdir])

#s1 add Qn2
add_all = True
include_qn2 = True
sat_rule = 1
most_unstable = False 
file_name = 'QL_summary_'+suffix+'_'+sfsuffix+'_iqn2'+str(include_qn2)[0]+'_sat'+str(sat_rule)+'_add'+str(add_all)[0]+'_mu'+str(most_unstable)[0]
call(['ETG_quasilinear.py',suffix,sfsuffix,'-n','-a','-q'])
call(['mv',file_name,outdir])
call(['mv',file_name+'.ps',outdir])

#s1 no add
add_all = False
include_qn2 = False
sat_rule = 1
most_unstable = False 
file_name = 'QL_summary_'+suffix+'_'+sfsuffix+'_iqn2'+str(include_qn2)[0]+'_sat'+str(sat_rule)+'_add'+str(add_all)[0]+'_mu'+str(most_unstable)[0]
call(['ETG_quasilinear.py',suffix,sfsuffix,'-n'])
call(['mv',file_name,outdir])
call(['mv',file_name+'.ps',outdir])

#s1 no add Qn2
add_all = False
include_qn2 = True
sat_rule = 1
most_unstable = False 
file_name = 'QL_summary_'+suffix+'_'+sfsuffix+'_iqn2'+str(include_qn2)[0]+'_sat'+str(sat_rule)+'_add'+str(add_all)[0]+'_mu'+str(most_unstable)[0]
call(['ETG_quasilinear.py',suffix,sfsuffix,'-n','-q'])
call(['mv',file_name,outdir])
call(['mv',file_name+'.ps',outdir])

#s4 add 
add_all = True
include_qn2 = False
sat_rule = 4
most_unstable = False 
file_name = 'QL_summary_'+suffix+'_'+sfsuffix+'_iqn2'+str(include_qn2)[0]+'_sat'+str(sat_rule)+'_add'+str(add_all)[0]+'_mu'+str(most_unstable)[0]
call(['ETG_quasilinear.py',suffix,sfsuffix,'-n','-a','-s 4'])
call(['mv',file_name,outdir])
call(['mv',file_name+'.ps',outdir])

#s4 add  qn2
add_all = True
include_qn2 = True
sat_rule = 4
most_unstable = False 
file_name = 'QL_summary_'+suffix+'_'+sfsuffix+'_iqn2'+str(include_qn2)[0]+'_sat'+str(sat_rule)+'_add'+str(add_all)[0]+'_mu'+str(most_unstable)[0]
call(['ETG_quasilinear.py',suffix,sfsuffix,'-n','-a','-s 4','-q'])
call(['mv',file_name,outdir])
call(['mv',file_name+'.ps',outdir])

#s4 no add  qn2
add_all = False
include_qn2 = True
sat_rule = 4
most_unstable = False 
file_name = 'QL_summary_'+suffix+'_'+sfsuffix+'_iqn2'+str(include_qn2)[0]+'_sat'+str(sat_rule)+'_add'+str(add_all)[0]+'_mu'+str(most_unstable)[0]
call(['ETG_quasilinear.py',suffix,sfsuffix,'-n','-s 4','-q'])
call(['mv',file_name,outdir])
call(['mv',file_name+'.ps',outdir])

#s4 no add  
add_all = False
include_qn2 = False
sat_rule = 4
most_unstable = False 
file_name = 'QL_summary_'+suffix+'_'+sfsuffix+'_iqn2'+str(include_qn2)[0]+'_sat'+str(sat_rule)+'_add'+str(add_all)[0]+'_mu'+str(most_unstable)[0]
call(['ETG_quasilinear.py',suffix,sfsuffix,'-n','-s 4'])
call(['mv',file_name,outdir])
call(['mv',file_name+'.ps',outdir])

#s4 no add most unstable only 
add_all = False
include_qn2 = False
sat_rule = 4
most_unstable = True
file_name = 'QL_summary_'+suffix+'_'+sfsuffix+'_iqn2'+str(include_qn2)[0]+'_sat'+str(sat_rule)+'_add'+str(add_all)[0]+'_mu'+str(most_unstable)[0]
call(['ETG_quasilinear.py',suffix,sfsuffix,'-n','-s 4','-m'])
call(['mv',file_name,outdir])
call(['mv',file_name+'.ps',outdir])


call(['cp','summary_'+suffix+'.csv',outdir])
call(['cp','parameters_'+suffix,outdir])
call(['cp','scanfiles'+sfsuffix+'/mode_info_all',outdir+'/mode_info_all_'+sfsuffix])

if archive:
    call(['cp','-r',outdir,arch_dir])

