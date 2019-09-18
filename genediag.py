#!/usr/bin/env bash
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import numpy as npy

import genetools
import geneplots

import argparse

import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='GENE Diagnostic Tools.')
parser.add_argument('--quick',       '-quick',       action='store_const',const=1,help='Trace the Fastest Growing Mode')
parser.add_argument('--findtau',     '-findtau',     action='store_const',const=1,help='Find tau = Zeff*Te/Ti')
parser.add_argument('--siunits',     '-siunits',     action='store_const',const=1,help='Covert to the SI Units')
parser.add_argument('--plotnrg',     '-plotnrg',     action='store_const',const=1,help='Plot the profiles in nrg_xxxx file')
parser.add_argument('--display',     '-display',     action='store_const',const=1,help='Display the plots')
parser.add_argument('--modeinfo',    '-modeinfo',    action='store_const',const=1,help='Retrieve information about a mode')
parser.add_argument('--findarea',    '-findarea',    action='store_const',const=1,help='Find the surface area of magnetic surface')
parser.add_argument('--logscale',    '-logscale',    action='store_const',const=1,help='Plot in log scale')
parser.add_argument('--plotgeom',    '-plotgeom',    action='store_const',const=1,help='Plot the Geometry')
parser.add_argument('--plotmodes',   '-plotmodes',   action='store_const',const=1,help='Plot the mode structures')
parser.add_argument('--findomega',   '-findomega',   action='store_const',const=1,help='Find omega and growth-rate of a mode')
parser.add_argument('--plotneoclass','-plotneoclass',action='store_const',const=1,help='Plot the profiles in neoclass_xxxx file')
parser.add_argument('inputs',nargs='*')

if parser.parse_args():
   args = parser.parse_args()
   quick     =    args.quick
   inputs    =    args.inputs
   findtau   =    args.findtau
   plotnrg   =    args.plotnrg
   display   =    args.display
   siunits   =    args.siunits
   modeinfo  =    args.modeinfo
   findarea  =    args.findarea
   logscale  =    args.logscale
   plotgeom  =    args.plotgeom
   plotmodes =    args.plotmodes 
   findomega =    args.findomega
   plotneoclass = args.plotneoclass
else:
   print('You need to select a function to implement'); sys.exit()

if display or logscale:
   if not (plotnrg or plotneoclass or plotgeom or plotmodes):
      print('WARNING: logscale can not be used independently.')
      print('         You need to accompany logscale with a plot.')

modeorder = []

if findtau:
   psiloc = inputs[0]
   profilefpath = inputs[1]
   tau = genetools.calc_tau(psiloc,profilepath=profilefpath)
   print('tau = Zeff*Te/Ti = %7.5f' % tau)
else:
   orderlist = [str('%04d') % (i+1) for i in range(9999)]
   orderlist.append('dat')
   orderlist.append('.dat')
   for item in inputs:
      if item in orderlist:
         if item == 'dat': item = '.dat'
         modeorder.append(item)

if not modeorder: sys.exit()

for mode in modeorder:
    if   mode.isdigit():
         modenumber = 'mode_%04d'       % int(mode)
         paramfname = 'parameters_%04d' % int(mode)
    else:
         modenumber = 'mode.dat'
         paramfname = 'parameters.dat'
    paramfpath = os.path.abspath(paramfname)
    if siunits or findomega:
       conv_units = genetools.units_conversion(paramfpath=paramfpath)

    if   findomega:
         if   mode.isdigit():
              fieldfname = 'field_%04d'     % int(mode)
         else:
              fieldfname = 'field.dat'
         fieldfpath = os.path.abspath(fieldfname)
         if not os.path.isfile(fieldfpath):
            print('File: %s is not in the given path.' % fieldfname); sys.exit()
         if quick:
            modefreq = genetools.find_mode_frequency(fieldfpath,method='fast-mode')
         else:
            modefreq = genetools.find_mode_frequency(fieldfpath,method='thorough')

         ky        = float(modefreq[modenumber]['ky'])
         omega_phi = float(modefreq[modenumber]['omega_phi'])
         omega_apr = float(modefreq[modenumber]['omega_apr'])
         gamma_phi = float(modefreq[modenumber]['gamma_phi'])
         gamma_apr = float(modefreq[modenumber]['gamma_apr'])

         omegaref = conv_units['cref']/conv_units['Lref']
         print('From Electric Potential Field:')
         print('ky    = %7.4f (Normalized), %7.4f (1/m)' % (ky,ky/conv_units['Lref']))
         print('Omega = %7.4f (Normalized), %7.4f (Hz)'  % (omega_phi,omega_phi*omegaref/(2.0*npy.pi)))
         print('Gamma = %7.4f (Normalized), %7.4f (Hz)'  % (gamma_phi,gamma_phi*omegaref/(2.0*npy.pi)))
         print('From Magnetic Potential Field:')
         print('ky    = %7.4f (Normalized), %7.4f (1/m)' % (ky,ky/conv_units['Lref']))
         print('Omega = %7.4f (Normalized), %7.4f (Hz)'  % (omega_apr,omega_apr*omegaref/(2.0*npy.pi)))
         print('Gamma = %7.4f (Normalized), %7.4f (Hz)'  % (gamma_apr,gamma_apr*omegaref/(2.0*npy.pi)))

         while True:
               if   sys.version_info[0] >=3:
                    saveomega = str(input('Do you want to update omega file?(Yes/No) ')).lower()
               elif sys.version_info[0] < 3:
                    saveomega = raw_input('Do you want to update omega file?(Yes/No) ').lower()
               if saveomega in ['yes','y']:
                  if   sys.version_info[0] >=3:
                       omegatype = input('Source?\n(1)Electric Potential,\n(2)Magnetic Potential.\nSelection:  ')
                  elif sys.version_info[0] < 3:
                       omegatype = input('Source?\n(1)Electric Potential,\n(2)Magnetic Potential.\nSelection:  ')
               if saveomega in ['yes','y','no','n']: break
         if saveomega in ['yes','y']:
            if 'dat' in fieldfpath:
               omegafpath = fieldfpath[:-9]+'omega'+fieldfpath[-4:]
            else:
               omegafpath = fieldfpath[:-10]+'omega'+fieldfpath[-5:]
            ofhand = open(omegafpath,'w')
            if omegatype == 1:
               ofhand.write('%5.3f\t%5.3f\t%5.3f\n' % (ky,gamma_phi,omega_phi))
            else:
               ofhand.write('%5.3f\t%5.3f\t%5.3f\n' % (ky,gamma_apr,omega_apr))
            ofhand.close()

    elif plotmodes:
         if   mode.isdigit():
              fieldfname = 'field_%04d'     % int(mode)
         else:
              fieldfname = 'field.dat'
         fieldfpath = os.path.abspath(fieldfname)
         if not os.path.isfile(fieldfpath):
            print('File: %s is not in the given path.' % fieldfname); sys.exit()
         fielddata = genetools.read_field(fieldfpath=fieldfpath)
         plotParam = {}
         if display: plotParam['display'] = True
         fieldplot = geneplots.plot_field(field=fielddata,setParam=plotParam)
    elif plotnrg:
         if   mode.isdigit():
              nrgfname   = 'nrg_%04d'       % int(mode)
         else:
              nrgfname   = 'nrg.dat'
         nrgfpath   = os.path.abspath(nrgfname)
         if not os.path.isfile(nrgfpath):
            print('File: %s is not in the given path.' % nrgfname); sys.exit()
         mergeplots = True
         plotParam = {}
         if display:    plotParam['display']    = True
         if logscale:   plotParam['logplots']   = True
         if mergeplots: plotParam['mergeplots'] = True
         if siunits:
            plotParam['siunits'] = True
            nrgdata    = genetools.read_nrg(nrgfpath,normalized=False)
         else:
            nrgdata    = genetools.read_nrg(nrgfpath,normalized=True)
         nrgplot    = geneplots.plot_nrg(nrgdata,setParam=plotParam)
         if display:
            for fig in nrgplot: fig.show()

    elif plotneoclass:
         if   mode.isdigit():
              neoclassfname   = 'neoclass_%04d'       % int(mode)
         else:
              neoclassfname   = 'neoclass.dat'
         neoclassfpath   = os.path.abspath(neoclassfname)
         if not os.path.isfile(neoclassfpath):
            print('File: %s is not in the given path.' % neoclassfname); sys.exit()
         plotParam = {}
         if logscale: plotParam['logplots'] = True
         if display:  plotParam['display']  = True
         if siunits:
            plotParam['siunits'] = True
            neoclassdata    = genetools.read_neoclass(neoclassfpath,normalized=False)
         else:
            neoclassdata    = genetools.read_neoclass(neoclassfpath,normalized=True)
         neoclassplot    = geneplots.plot_neoclass(neoclassdata,setParam=plotParam)

    elif plotgeom:
         if   mode.isdigit():
              geomfname  = 'tracer_efit_%04d' % int(mode)
              geomfpath  = os.path.abspath(geomfname)
              if not os.path.isfile(geomfpath):
                 geomfname  = 's_alpha_%04d' % int(mode)
                 geomfpath  = os.path.abspath(geomfname)
                 if not os.path.isfile(geomfpath):
                    geomfname  = 'chease_%04d' % int(mode)
                    geomfpath  = os.path.abspath(geomfname)
         else:
              geomfname  = 'tracer_efit.dat'
              geomfpath  = os.path.abspath(geomfname)
              if not os.path.isfile(geomfpath):
                 geomfname  = 's_alpha.dat'
                 geomfpath  = os.path.abspath(geomfname)
                 if not os.path.isfile(geomfpath):
                    geomfname  = 'chease.dat'
                    geomfpath  = os.path.abspath(geomfname)

         if not os.path.isfile(geomfpath):
            print('File: %s is not in the given path.' % geomfname); sys.exit()
         plotParam = {}
         if display: plotParam['display'] = True
         geomplot   = geneplots.plot_geometry(geometryfpath=geomfpath,setParam=plotParam)

    elif findarea:
         if   mode.isdigit():
              geomfname  = 'tracer_efit_%04d' % int(mode)
              geomfpath  = os.path.abspath(geomfname)
              if not os.path.isfile(geomfpath):
                 geomfname  = 's_alpha_%04d' % int(mode)
                 geomfpath  = os.path.abspath(geomfname)
                 if not os.path.isfile(geomfpath):
                    geomfname  = 'chease_%04d' % int(mode)
                    geomfpath  = os.path.abspath(geomfname)
         else:
              geomfname  = 'tracer_efit.dat'
              geomfpath  = os.path.abspath(geomfname)
              if not os.path.isfile(geomfpath):
                 geomfname  = 's_alpha.dat'
                 geomfpath  = os.path.abspath(geomfname)
                 if not os.path.isfile(geomfpath):
                    geomfname  = 'chease.dat'
                    geomfpath  = os.path.abspath(geomfname)

         geomfpath  = os.path.abspath(geomfname)
         area  = genetools.calculate_surface_area(geomfpath,paramfpath)
         print('Magnetic Surface Area = %7.5f' % area)

    elif modeinfo:
         genepath = os.path.abspath('./')
         if genepath[-1]!='/': genepath+='/'

         fluxinfo = genetools.flux_info(genepath)
         kymin = fluxinfo.keys()[0]
         print('Mode Flux Info:')
         print("ky = %5.3f" % (kymin))
         print("Xi/Xe = %5.3f" % (fluxinfo[kymin]['i']['Chi']/fluxinfo[kymin]['e']['Chi']))
         print("De/Xe = %5.3f" % (fluxinfo[kymin]['e']['Dee']/fluxinfo[kymin]['e']['Chi']))
         print("Dz/Xe = %5.3f" % (fluxinfo[kymin]['z']['Dee']/fluxinfo[kymin]['e']['Chi']))
         print("De/(Xe+Xi) = %5.3f" % (fluxinfo[kymin]['e']['Dee']/(fluxinfo[kymin]['e']['Chi']+fluxinfo[kymin]['i']['Chi'])))
         print("Dz/(Xe+Xi) = %5.3f" % (fluxinfo[kymin]['z']['Dee']/(fluxinfo[kymin]['e']['Chi']+fluxinfo[kymin]['i']['Chi'])))
         print("Instability Type: %s" % (fluxinfo[kymin]['Type']))

         modeinfo = genetools.mode_info(genepath)
         kymin = modeinfo.keys()[0]
         print('Mode General Info:')
         print("ky = %5.3f" % (kymin))
         print("Qem/Qes = %5.3f" % (modeinfo[kymin]['Qem/Qes']))
         print("Correlation Length = %5.3f" % (modeinfo[kymin]['corr_len']))
         print("Parity Factor (A||) = %5.3f" % (modeinfo[kymin]['parity_factor_apar']))
         print("Parity Factor (Phi) = %5.3f" % (modeinfo[kymin]['parity_factor_phi']))
         print("E|| Cancelllation = %5.3f" % (modeinfo[kymin]['Epar_Cancellation']))
         print("Zavg = %5.3f" % (modeinfo[kymin]['zavg']))
         print("Gamma = %5.3f" % (modeinfo[kymin]['gamma']))
         print("Omega = %5.3f" % (modeinfo[kymin]['omega']))
         print("Sqrt(Gamma**2+Omega**2) = %5.3f" % (npy.sqrt(modeinfo[kymin]['gamma']**2+modeinfo[kymin]['omega']**2)))
         print("Instability Type: %s" % (modeinfo[kymin]['Type']))

         if   sys.version_info[0] >=3:
              saveinfo = str(input('Do you want to save info to file?(Yes/No) ')).lower()
         elif sys.version_info[0] < 3:
              saveinfo = raw_input('Do you want to save info to file?(Yes/No) ').lower()

         if saveinfo in ['yes','y']:
            if not os.path.isdir(genepath+"report"):
               os.system('mkdir '+genepath+"report")
               reportpath = genepath+"report/"
            else:
               reportpath = genepath+"report/"

            infofpath = reportpath+'gene.info'
            ofhand = open(infofpath,'w')

            ofhand.write("Mode Flux Info:\n")
            ofhand.write("ky = %5.3f\n" % kymin)
            ofhand.write("Xi/Xe = %5.3f\n" % (fluxinfo[kymin]['i']['Chi']/fluxinfo[kymin]['e']['Chi']))
            ofhand.write("De/Xe = %5.3f\n" % (fluxinfo[kymin]['e']['Dee']/fluxinfo[kymin]['e']['Chi']))
            ofhand.write("Dz/Xe = %5.3f\n" % (fluxinfo[kymin]['z']['Dee']/fluxinfo[kymin]['e']['Chi']))
            ofhand.write("De/(Xe+Xi) = %5.3f\n" % (fluxinfo[kymin]['e']['Dee']/(fluxinfo[kymin]['e']['Chi']+fluxinfo[kymin]['i']['Chi'])))
            ofhand.write("Dz/(Xe+Xi) = %5.3f\n" % (fluxinfo[kymin]['z']['Dee']/(fluxinfo[kymin]['e']['Chi']+fluxinfo[kymin]['i']['Chi'])))
            ofhand.write("Instability Type: %s\n" % (fluxinfo[kymin]['Type']))
            ofhand.write("\n")

            ofhand.write("Mode General Info:\n")
            ofhand.write("Qem/Qes = %5.3f\n" % (modeinfo[kymin]['Qem/Qes']))
            ofhand.write("Correlation Length = %5.3f\n" % (modeinfo[kymin]['corr_len']))
            ofhand.write("Parity Factor (A||) = %5.3f\n" % (modeinfo[kymin]['parity_factor_apar']))
            ofhand.write("Parity Factor (Phi) = %5.3f\n" % (modeinfo[kymin]['parity_factor_phi']))
            ofhand.write("E|| Cancelllation = %5.3f\n" % (modeinfo[kymin]['Epar_Cancellation']))
            ofhand.write("Zavg = %5.3f\n" % (modeinfo[kymin]['zavg']))
            ofhand.write("Gamma = %5.3f\n" % (modeinfo[kymin]['gamma']))
            ofhand.write("Omega = %5.3f\n" % (modeinfo[kymin]['omega']))
            ofhand.write("Sqrt(Gamma**2+Omega**2) = %5.3f\n" % (npy.sqrt(modeinfo[kymin]['gamma']**2+modeinfo[kymin]['omega']**2)))
            ofhand.write("Instability Type : %s\n" % (modeinfo[kymin]['Type']))

            ofhand.close()





