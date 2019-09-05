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
parser.add_argument('--siunits',     '-siunits',     action='store_const',const=1,help='Covert to the SI Units')
parser.add_argument('--plotnrg',     '-plotnrg',     action='store_const',const=1,help='Plot the profiles in nrg_xxxx file')
parser.add_argument('--display',     '-display',     action='store_const',const=1,help='Display the plots')
parser.add_argument('--findarea',    '-findarea',    action='store_const',const=1,help='Find the surface area of magnetic surface')
parser.add_argument('--logscale',    '-logscale',    action='store_const',const=1,help='Plot in log scale')
parser.add_argument('--plotgeom',    '-plotgeom',    action='store_const',const=1,help='Plot the Geometry')
parser.add_argument('--plotmodes',   '-plotmodes',   action='store_const',const=1,help='Plot the mode structures')
parser.add_argument('--findomega',   '-findomega',   action='store_const',const=1,help='Find omega and growth-rate of a mode')
parser.add_argument('--plotneoclass','-plotneoclass',action='store_const',const=1,help='Plot the profiles in neoclass_xxxx file')
parser.add_argument('modeorder',nargs='+')

if parser.parse_args():
   args = parser.parse_args()
   quick     =    args.quick
   plotnrg   =    args.plotnrg
   display   =    args.display
   siunits   =    args.siunits
   findarea  =    args.findarea
   logscale  =    args.logscale
   plotgeom  =    args.plotgeom
   plotmodes =    args.plotmodes 
   findomega =    args.findomega
   modeorder =    args.modeorder
   plotneoclass = args.plotneoclass
else:
   print('You need to select a function to implement'); sys.exit()

if display or logscale:
   if not (plotnrg or plotneoclass or plotgeom or plotmodes):
      print('WARNING: logscale can not be used independently.')
      print('         You need to accompany logscale with a plot.')

if modeorder[0]=='dat': modeorder[0] = '.dat'

for mode in modeorder:
    if   mode.isdigit():
         modenumber = 'mode_%04d'       % int(mode)
         paramfname = 'parameters_%04d' % int(mode)
    else:
         modenumber = 'mode.dat'
         paramfname = 'parameters.dat'
    paramfpath = os.path.abspath(paramfname)
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
         omegaref = conv_units['cref']/conv_units['Lref']
         print('From Electric Potential Field:')
         print('Omega = %7.4f (Normalized), %7.4f (Hz)' % (modefreq[modenumber]['omega_phi'],modefreq[modenumber]['omega_phi']*omegaref/(2.0*npy.pi)))
         print('Gamma = %7.4f (Normalized), %7.4f (Hz)' % (modefreq[modenumber]['gamma_phi'],modefreq[modenumber]['gamma_phi']*omegaref/(2.0*npy.pi)))
         print('From Magnetic Potential Field:')
         print('Omega = %7.4f (Normalized), %7.4f (Hz)' % (modefreq[modenumber]['omega_apr'],modefreq[modenumber]['omega_apr']*omegaref/(2.0*npy.pi)))
         print('Gamma = %7.4f (Normalized), %7.4f (Hz)' % (modefreq[modenumber]['gamma_apr'],modefreq[modenumber]['gamma_apr']*omegaref/(2.0*npy.pi)))
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
         else:
              geomfname  = 'tracer_efit.dat'
         geomfpath  = os.path.abspath(geomfname)
         if not os.path.isfile(geomfpath):
            print('File: %s is not in the given path.' % geomfname); sys.exit()
         geomplot   = geneplots.plot_geometry(geometryfpath=geomfpath)
    elif findarea:
         if   mode.isdigit():
              geomfname  = 'tracer_efit_%04d' % int(mode)
         else:
              geomfname  = 'tracer_efit.dat'
         geomfpath  = os.path.abspath(geomfname)
         area  = genetools.calculate_surface_area(geomfpath,paramfpath)
         print area

