#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import math
import warnings
import genetools
import efittools
import read_iterdb
import matplotlib.cbook
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdfh

import numpy as npy

from read_write_geometry import *
from scipy.interpolate import interp1d

warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

def create_k_grid(x_grid):
    nx  = npy.size(x_grid)
    dx  = min(npy.diff(x_grid))
    dKx = 2.0*npy.pi/nx/dx
    iKx = npy.zeros(nx)
    for ix in range(nx):
        if   ix<=nx/2-1:
             iKx[ix] = ix*dx
        elif ix>=nx/2:
             iKx[ix] =-(nx-ix)*dx
    return npy.fft.fftshift(iKx)


def findall(inlist,item):
    inds = []
    for i,character in enumerate(inlist):
        if character==item: inds.append(i)
    return inds

def str2bool(vin):
    return (vin.strip()).lower() in ('t','.t.','true','.true.')


def plot_nrg(nrgdata,reportpath='',bgn_t=None,end_t=None,fraction=0.9,setParam={}):
   #Developed by Ehab Hassan on 2019-07-22
    inrgf=nrgdata.keys()[0]
    if    inrgf[-3:] == 'dat': inrgfpath = inrgf[:-7]
    else:                      inrgfpath = inrgf[:-8]
    if 'report' not in reportpath:
       if not os.path.isdir(inrgfpath+"report"):
          os.system('mkdir '+inrgfpath+"report")
          reportpath = inrgfpath+"report/"
       else:
          reportpath = inrgfpath+"report/"
    elif reportpath[-1] != "/":
       reportpath += "/"

    if  'mergeplots' in setParam:
         mergeplots = setParam['mergeplots']
    else:
         mergeplots = False

    if  'logplots' in setParam:
         logplots = setParam['logplots']
    else:
         logplots = False

    if  'display' in setParam:
         display = setParam['display']
    else:
         display = False

    if  'siunits' in setParam:
         siunits = setParam['siunits']
    else:
         siunits = False

    time      = nrgdata[inrgf]['time']
    specstype = nrgdata[inrgf].keys()
    specstype.remove('time')

    for inrgf in nrgdata:
        if 'dat' in inrgf:
           isimfpath = inrgf[:-7]
           inrgfpath = inrgf[-7:]
           inrgfext  = inrgf[-4:]
        else:
           isimfpath = inrgf[:-8]
           inrgfpath = inrgf[-8:]
           inrgfext  = inrgf[-5:]

        parameters = genetools.read_parameters("%sparameters%s" %(isimfpath,inrgfext))
        if 'n0_global' in parameters['box']:
           titletxt = '(k_y=%6.3f,n_0=%5d)' % (parameters['box']['kymin'],parameters['box']['n0_global'])
        else:
           titletxt = '(k_y=%6.3f)' % (parameters['box']['kymin'])

        if siunits:
           geomfpath  = "%stracer_efit%s" %(isimfpath,inrgfext)
           area = genetools.calculate_surface_area(parameters=parameters,geometry=geomfpath)

        time = npy.array(nrgdata[inrgf]['time'])
        if bgn_t == None: bgn_t = time[0]
        if end_t == None: end_t = time[-1]
        nonlinear = False
        if 'nonlinear' in parameters['general']:
           if parameters['general']['nonlinear']:
              nonlinear = True
              avg_bgn_t_ind = int(npy.size(time)*fraction)
             #avg_bgn_t = time[avg_bgn_t_ind] 
        bgn_t_ind = npy.argmin(abs(npy.array(time)-bgn_t))
        end_t_ind = npy.argmin(abs(npy.array(time)-end_t))
        ntimes    = end_t_ind-bgn_t_ind+1

        for ispecs in specstype:
            nfig     = plt.figure('n-'+inrgfpath)
            axhand01 = nfig.add_subplot(1,1,1)
            dens     = nrgdata[inrgf][ispecs]['n']
            axhand01.plot(time[bgn_t_ind:end_t_ind+1],dens[bgn_t_ind:end_t_ind+1],label=ispecs)
            axhand01.set_title('n%s' %(titletxt))
            axhand01.set_xlabel('Time')
            axhand01.set_ylabel('Density')
            if logplots: axhand01.set_yscale('symlog')
            axhand01.legend()

            uparafig = plt.figure('upara-'+inrgfpath)
            axhand02 = uparafig.add_subplot(1,1,1)
            upara    = nrgdata[inrgf][ispecs]['upara']
            axhand02.plot(time[bgn_t_ind:end_t_ind+1],upara[bgn_t_ind:end_t_ind+1],label=ispecs)
            axhand02.set_title('$U_{||}%s$' %(titletxt))
            axhand02.set_xlabel('Time')
            axhand02.set_ylabel('Parallel Velocity')
            if logplots: axhand02.set_yscale('symlog')
            axhand02.legend()

            if mergeplots:
               Tfig     = plt.figure('T-'+inrgfpath)
               axhand03 = Tfig.add_subplot(1,1,1)
               Tpara = nrgdata[inrgf][ispecs]['Tpara']
               axhand03.plot(time[bgn_t_ind:end_t_ind+1],Tpara[bgn_t_ind:end_t_ind+1],linestyle='-', label='$T_{\\parallel,%s}$' % ispecs)
               Tperp = nrgdata[inrgf][ispecs]['Tperp']
               axhand03.plot(time[bgn_t_ind:end_t_ind+1],Tperp[bgn_t_ind:end_t_ind+1],linestyle='--',label='$T_{\\perp,%s}$' % ispecs)
               axhand03.set_title('$T_{\\parallel,\\perp}%s$' %(titletxt))
               axhand03.set_xlabel('Time')
               axhand03.set_ylabel('Temperature')
               if logplots: axhand03.set_yscale('symlog')
               axhand03.legend()
            else:
               Tparafig = plt.figure('Tpara-'+inrgfpath)
               axhand03 = Tparafig.add_subplot(1,1,1)
               Tpara = nrgdata[inrgf][ispecs]['Tpara']
               axhand03.plot(time[bgn_t_ind:end_t_ind+1],T[parabgn_t_ind:end_t_ind+1],label=ispecs)
               axhand03.set_title('$T_{\\parallel}%s$' %(titletxt))
               axhand03.set_xlabel('Time')
               axhand03.set_ylabel('Parallel Temperature')
               if logplots: axhand03.set_yscale('symlog')
               axhand03.legend()

               Tperpfig = plt.figure('Tperp-'+inrgfpath)
               axhand03 = Tperpfig.add_subplot(1,1,1)
               Tperp = nrgdata[inrgf][ispecs]['Tperp']
               axhand03.plot(time[bgn_t_ind:end_t_ind+1],Tperp[bgn_t_ind:end_t_ind+1],label=ispecs)
               axhand03.set_title('$T_{\\perp}%s$' %(titletxt))
               axhand03.set_xlabel('Time')
               axhand03.set_ylabel('Transverse Temperature')
               if logplots: axhand03.set_yscale('symlog')
               axhand03.legend()

            if mergeplots:
               PFluxfig = plt.figure('PFlux-'+inrgfpath)
               axhand04 = PFluxfig.add_subplot(1,1,1)
               if siunits:
                  PFluxes = nrgdata[inrgf][ispecs]['PFluxes']*area
                  PFluxem = nrgdata[inrgf][ispecs]['PFluxem']*area
                 #PFluxes*= 1.0e6/(1.0e3*1.0e-19)
                 #PFluxem*= 1.0e6/(1.0e3*1.0e-19)
                  axhand04.set_ylabel('Particle Flux (Particles/s)')
                 #axhand04.set_ylabel('Particle Flux (MW/keV)')
               else:
                  PFluxes = nrgdata[inrgf][ispecs]['PFluxes']
                  PFluxem = nrgdata[inrgf][ispecs]['PFluxem']
                  axhand04.set_ylabel('Particle Flux)')
               if nonlinear:
                  axhand04.plot(time[bgn_t_ind:end_t_ind+1],PFluxes[bgn_t_ind:end_t_ind+1],linestyle='-',label='$\\Gamma_{es,%s}$=%7.5e' % (ispecs,npy.mean(PFluxes[avg_bgn_t_ind:])))
                  axhand04.plot(time[bgn_t_ind:end_t_ind+1],PFluxem[bgn_t_ind:end_t_ind+1],linestyle=':',label='$\\Gamma_{em,%s}$=%7.5e' % (ispecs,npy.mean(PFluxes[avg_bgn_t_ind:])))
               else:
                  axhand04.plot(time[bgn_t_ind:end_t_ind+1],PFluxes[bgn_t_ind:end_t_ind+1],linestyle='-',label='$\\Gamma_{es,%s}$=%7.5e' % (ispecs,PFluxes[-1]))
                  axhand04.plot(time[bgn_t_ind:end_t_ind+1],PFluxem[bgn_t_ind:end_t_ind+1],linestyle=':',label='$\\Gamma_{em,%s}$=%7.5e' % (ispecs,PFluxem[-1]))
               axhand04.set_title('$\\Gamma_{es,em}%s$' %(titletxt))
               axhand04.set_xlabel('Time')
               if logplots: axhand04.set_yscale('symlog')
               axhand04.legend()
            else:
               PFluxesfig = plt.figure('PFluxes-'+inrgfpath)
               axhand04 = PFluxesfig.add_subplot(1,1,1)
               if siunits:
                  PFluxes = nrgdata[inrgf][ispecs]['PFluxes']*area
                 #PFluxes*= 1.0e6/(1.0e3*1.0e-19)
                  axhand04.set_ylabel('Electrostatic Particle Flux (Particles/s)')
                 #axhand04.set_ylabel('Electrostatic Particle Flux (MW/keV)')
               else:
                  PFluxes = nrgdata[inrgf][ispecs]['PFluxes']
                  axhand04.set_ylabel('Electrostatic Particle Flux')
               if nonlinear:
                  axhand04.plot(time[bgn_t_ind:end_t_ind+1],PFluxes[bgn_t_ind:end_t_ind+1],linestyle='-',label='$\\Gamma_{es,%s}$=%7.5e' % (ispecs,npy.mean(PFluxes[avg_bgn_t_ind:])))
               else:
                  axhand04.plot(time[bgn_t_ind:end_t_ind+1],PFluxes[bgn_t_ind:end_t_ind+1],label='$\\Gamma_{es,%s}$=%7.5e' % (ispecs,PFluxes[-1]))
               axhand04.set_title('$\\Gamma_{es}%s$' %(titletxt))
               axhand04.set_xlabel('Time')
               if logplots: axhand04.set_yscale('symlog')
               axhand04.legend()

               PFluxemfig = plt.figure('PFluxem-'+inrgfpath)
               axhand04 = PFluxemfig.add_subplot(1,1,1)
               if siunits:
                  PFluxem = nrgdata[inrgf][ispecs]['PFluxem']*area
                 #PFluxem*= 1.0e6/(1.0e3*1.0e-19)
                  axhand04.set_ylabel('Electromagnetic Particle Flux (Particles/s)')
                 #axhand04.set_ylabel('Electromagnetic Particle Flux (MW/keV)')
               else:
                  PFluxem = nrgdata[inrgf][ispecs]['PFluxem']
                  axhand04.set_ylabel('Electromagnetic Particle Flux')
               if nonlinear:
                  axhand04.plot(time[bgn_t_ind:end_t_ind+1],PFluxem[bgn_t_ind:end_t_ind+1],linestyle=':',label='$\\Gamma_{em,%s}$=%7.5e' % (ispecs,npy.mean(PFluxes[avg_bgn_t_ind:])))
               else:
                  axhand04.plot(time[bgn_t_ind:end_t_ind+1],PFluxem[bgn_t_ind:end_t_ind+1],label='$\\Gamma_{em,%s}$=%7.5e' % (ispecs,PFluxem[-1]))
               axhand04.set_title('$\\Gamma_{em}%s$' %(titletxt))
               axhand04.set_xlabel('Time')
               if logplots: axhand04.set_yscale('symlog')
               axhand04.legend()

            if mergeplots:
               HFluxfig = plt.figure('HFlux-'+inrgfpath)
               axhand05 = HFluxfig.add_subplot(1,1,1)
               if siunits:
                  HFluxes = nrgdata[inrgf][ispecs]['HFluxes']*area
                  HFluxes/= 1.0e6
                  HFluxem = nrgdata[inrgf][ispecs]['HFluxem']*area
                  HFluxem/= 1.0e6
                  axhand05.set_ylabel('Heat Flux (MW)')
               else:
                  HFluxes = nrgdata[inrgf][ispecs]['HFluxes']
                  HFluxem = nrgdata[inrgf][ispecs]['HFluxem']
                  axhand05.set_ylabel('Heat Flux')
               if nonlinear:
                  axhand05.plot(time[bgn_t_ind:end_t_ind+1],HFluxes[bgn_t_ind:end_t_ind+1],linestyle='-',label='$Q_{es,%s}$=%7.5e' % (ispecs,npy.mean(HFluxes[avg_bgn_t_ind:])))
                  axhand05.plot(time[bgn_t_ind:end_t_ind+1],HFluxem[bgn_t_ind:end_t_ind+1],linestyle=':',label='$Q_{em,%s}$=%7.5e' % (ispecs,npy.mean(HFluxem[avg_bgn_t_ind:])))
               else:
                  axhand05.plot(time[bgn_t_ind:end_t_ind+1],HFluxes[bgn_t_ind:end_t_ind+1],linestyle='-',label='$Q_{es,%s}$=%7.5e' % (ispecs,HFluxes[-1]))
                  axhand05.plot(time[bgn_t_ind:end_t_ind+1],HFluxem[bgn_t_ind:end_t_ind+1],linestyle=':',label='$Q_{em,%s}$=%7.5e' % (ispecs,HFluxem[-1]))
               axhand05.set_title('$Q_{es,em}%s$' %(titletxt))
               axhand05.set_xlabel('Time')
               if logplots: axhand05.set_yscale('symlog')
               axhand05.legend()
            else:
               HFluxesfig = plt.figure('HFluxes-'+inrgfpath)
               axhand05 = HFluxesfig.add_subplot(1,1,1)
               if siunits:
                  HFluxes = nrgdata[inrgf][ispecs]['HFluxes']*area
                  HFluxes/= 1.0e6
                  axhand05.set_ylabel('Electrostatic Heat Flux (MW)')
               else:
                  HFluxes = nrgdata[inrgf][ispecs]['HFluxes']
                  axhand05.set_ylabel('Electrostatic Heat Flux')
               if nonlinear:
                  axhand05.plot(time[bgn_t_ind:end_t_ind+1],HFluxes[bgn_t_ind:end_t_ind+1],linestyle='-',label='$Q_{es,%s}$=%7.5e' % (ispecs,npy.mean(HFluxes[avg_bgn_t_ind:])))
               else:
                  axhand05.plot(time[bgn_t_ind:end_t_ind+1],HFluxes[bgn_t_ind:end_t_ind+1],label='$Q_{es,%s}$=%7.5e' % (ispecs,HFluxes[-1]))
               axhand05.set_title('$Q_{es}%s$' %(titletxt))
               axhand05.set_xlabel('Time')
               if logplots: axhand05.set_yscale('symlog')
               axhand05.legend()

               HFluxemfig = plt.figure('HFluxem-'+inrgfpath)
               axhand05 = HFluxemfig.add_subplot(1,1,1)
               if siunits:
                  HFluxem = nrgdata[inrgf][ispecs]['HFluxem']*area
                  HFluxem/= 1.0e6
                  axhand05.set_ylabel('Electromagnetic Heat Flux (MW)')
               else:
                  HFluxem = nrgdata[inrgf][ispecs]['HFluxem']
                  axhand05.set_ylabel('Electromagnetic Heat Flux')
               if nonlinear:
                  axhand04.plot(time[bgn_t_ind:end_t_ind+1],HFluxem[bgn_t_ind:end_t_ind+1],linestyle='-',label='$Q_{em,%s}$=%7.5e' % (ispecs,npy.mean(HFluxem[avg_bgn_t_ind:])))
               else:
                  axhand05.plot(time[bgn_t_ind:end_t_ind+1],HFluxem[bgn_t_ind:end_t_ind+1],label='$Q_{em,%s}$=%7.5e' % (ispecs,HFluxem[-1]))
               axhand05.set_title('$Q_{em}%s$' %(titletxt))
               axhand05.set_xlabel('Time')
               if logplots: axhand05.set_yscale('symlog')
               axhand05.legend()

            if mergeplots:
               Viscosfig = plt.figure('Viscos-'+inrgfpath)
               axhand06 = Viscosfig.add_subplot(1,1,1)
               if siunits:
                  Viscoses = nrgdata[inrgf][ispecs]['Viscoses']*area
                  Viscosem = nrgdata[inrgf][ispecs]['Viscosem']*area
                  axhand06.set_ylabel('Stress Tensor (N.m)')
               else:
                  Viscoses = nrgdata[inrgf][ispecs]['Viscoses']
                  Viscosem = nrgdata[inrgf][ispecs]['Viscosem']
                  axhand06.set_ylabel('Stress Tensor')
               if nonlinear:
                  axhand06.plot(time[bgn_t_ind:end_t_ind+1],Viscoses[bgn_t_ind:end_t_ind+1],linestyle='-',label='$\\Pi_{es,%s}$=%7.5e' % (ispecs,npy.mean(Viscoses[avg_bgn_t_ind:])))
                  axhand06.plot(time[bgn_t_ind:end_t_ind+1],Viscosem[bgn_t_ind:end_t_ind+1],linestyle=':',label='$\\Pi_{em,%s}$=%7.5e' % (ispecs,npy.mean(Viscosem[avg_bgn_t_ind:])))
               else:
                  axhand06.plot(time[bgn_t_ind:end_t_ind+1],Viscoses[bgn_t_ind:end_t_ind+1],linestyle='-',label='$\\Pi_{es,%s}$=%7.5e' % (ispecs,Viscoses[-1]))
                  axhand06.plot(time[bgn_t_ind:end_t_ind+1],Viscosem[bgn_t_ind:end_t_ind+1],linestyle=':',label='$\\Pi_{em,%s}$=%7.5e' % (ispecs,Viscosem[-1]))
               axhand06.set_title('$\\Pi_{es,em}%s$' %(titletxt))
               axhand06.set_xlabel('Time')
               if logplots: axhand06.set_yscale('symlog')
               axhand06.legend()
            else:
               Viscosesfig = plt.figure('Viscoses-'+inrgfpath)
               axhand06 = Viscosesfig.add_subplot(1,1,1)
               if siunits:
                  Viscoses = nrgdata[inrgf][ispecs]['Viscoses']*area
                  axhand06.set_ylabel('Electrostatic Stress Tensor(N.m)')
               else:
                  Viscoses = nrgdata[inrgf][ispecs]['Viscoses']
                  axhand06.set_ylabel('Electrostatic Stress Tensor')
               if nonlinear:
                  axhand06.plot(time[bgn_t_ind:end_t_ind+1],Viscoses[bgn_t_ind:end_t_ind+1],linestyle='-',label='$\\Pi_{es,%s}$=%7.5e' % (ispecs,npy.mean(Viscoses[avg_bgn_t_ind:])))
               else:
                  axhand06.plot(time[bgn_t_ind:end_t_ind+1],Viscoses[bgn_t_ind:end_t_ind+1],label='$\\Pi_{es,%s}$=%7.5e' % (ispecs,Viscoses[-1]))
               axhand06.set_title('$\\Pi_{es}%s$' %(titletxt))
               axhand06.set_xlabel('Time')
               if logplots: axhand06.set_yscale('symlog')
               axhand06.legend()

               Viscosemfig = plt.figure('Viscosem-'+inrgfpath)
               axhand06 = Viscosemfig.add_subplot(1,1,1)
               if siunits:
                  Viscosem = nrgdata[inrgf][ispecs]['Viscosem']*area
                  axhand06.set_ylabel('Electromagnetic Stress Tensor(N.m)')
               else:
                  Viscosem = nrgdata[inrgf][ispecs]['Viscosem']
                  axhand06.set_ylabel('Electromagnetic Stress Tensor')
               if nonlinear:
                  axhand06.plot(time[bgn_t_ind:end_t_ind+1],Viscosem[bgn_t_ind:end_t_ind+1],linestyle='-',label='$\\Pi_{em,%s}$=%7.5e' % (ispecs,npy.mean(Viscosem[avg_bgn_t_ind:])))
               else:
                  axhand06.plot(time[bgn_t_ind:end_t_ind+1],Viscosem[bgn_t_ind:end_t_ind+1],label='$\\Pi_{em,%s}$=%7.5e' % (ispecs,Viscosem[-1]))
               axhand06.set_title('$\\Pi_{em}%s$' %(titletxt))
               axhand06.set_xlabel('Time')
               if logplots: axhand06.set_yscale('symlog')
               axhand06.legend()

        if display: plt.show()

        figlist = []

        figlist.append(nfig)
        nfig.savefig(reportpath+'n_%s.png' % (inrgfpath))
        plt.close(nfig)

        figlist.append(uparafig)
        uparafig.savefig(reportpath+'upara_%s.png' % (inrgfpath))
        plt.close(uparafig)

        if mergeplots:
           figlist.append(Tfig)
           Tfig.savefig(reportpath+'T_%s.png' % (inrgfpath))
           plt.close(Tfig)
        else:
           figlist.append(Tparafig)
           Tparafig.savefig(reportpath+'Tpara_%s.png' % (inrgfpath))
           plt.close(Tparafig)

           figlist.append(Tperpfig)
           Tperpfig.savefig(reportpath+'Tperp_%s.png' % (inrgfpath))
           plt.close(Tperpfig)

        if mergeplots:
           figlist.append(PFluxfig)
           PFluxfig.savefig(reportpath+'PFlux_%s.png' % (inrgfpath))
           plt.close(PFluxfig)
        else:
           figlist.append(PFluxesfig)
           PFluxesfig.savefig(reportpath+'PFluxes_%s.png' % (inrgfpath))
           plt.close(PFluxesfig)

           figlist.append(PFluxemfig)
           PFluxemfig.savefig(reportpath+'PFluxem_%s.png' % (inrgfpath))
           plt.close(PFluxemfig)

        if mergeplots:
           figlist.append(HFluxfig)
           HFluxfig.savefig(reportpath+'HFlux_%s.png' % (inrgfpath))
           plt.close(HFluxfig)
        else:
           figlist.append(PFluxesfig)
           HFluxesfig.savefig(reportpath+'HFluxes_%s.png' % (inrgfpath))
           plt.close(HFluxesfig)

           figlist.append(PFluxemfig)
           HFluxemfig.savefig(reportpath+'HFluxem_%s.png' % (inrgfpath))
           plt.close(HFluxemfig)

        if mergeplots:
           figlist.append(Viscosfig)
           Viscosfig.savefig(reportpath+'Viscos_%s.png' % (inrgfpath))
           plt.close(Viscosfig)
        else:
           figlist.append(Viscosesfig)
           Viscosesfig.savefig(reportpath+'Viscoses_%s.png' % (inrgfpath))
           plt.close(Viscosesfig)

           figlist.append(Viscosemfig)
           Viscosemfig.savefig(reportpath+'Viscosem_%s.png' % (inrgfpath))
           plt.close(Viscosemfig)

    return figlist 


def plot_neoclass(neoclassdata,reportpath='',setParam={}):
   #Developed by Ehab Hassan on 2019-08-20
    ineoclassf=neoclassdata.keys()[0]
    if    ineoclassf[-3:] == 'dat': ineoclassfpath = ineoclassf[:-12]
    else:                           ineoclassfpath = ineoclassf[:-13]
    if 'report' not in reportpath:
       if not os.path.isdir(ineoclassfpath+"report"):
          os.system('mkdir '+ineoclassfpath+"report")
          reportpath = ineoclassfpath+"report/"
       else:
          reportpath = ineoclassfpath+"report/"
    elif reportpath[-1] != "/":
       reportpath += "/"

    if  'logplots' in setParam:
         logplots = setParam['logplots']
    else:
         logplots = False

    if  'display' in setParam:
         display = setParam['display']
    else:
         display = False

    time      = neoclassdata[ineoclassf]['time']
    specstype = neoclassdata[ineoclassf].keys()
    specstype.remove('time')

    for ineoclassf in neoclassdata:
        if 'dat' in ineoclassf:
           isimfpath = ineoclassf[:-12]
           ineoclassfpath = ineoclassf[-12:]
           ineoclassfext  = ineoclassf[-4:]
        else:
           isimfpath = ineoclassf[:-13]
           ineoclassfpath = ineoclassf[-13:]
           ineoclassfext  = ineoclassf[-5:]

        paramfpath = "%sparameters%s"  %(isimfpath,ineoclassfext)
        parameters = genetools.read_parameters(paramfpath)
        if 'x0' in parameters['box']:
           titletxt = '(x_0=%7.5f)' % parameters['box']['x0']

        geomfpath  = "%stracer_efit%s" %(isimfpath,ineoclassfext)
        area = genetools.calculate_surface_area(parameters=parameters,geometry=geomfpath)

        bgnind = 3000
        time = neoclassdata[ineoclassf]['time']
        for ispecs in specstype:
            PFluxfig = plt.figure('PFlux-'+ineoclassfpath)
            axhand = PFluxfig.add_subplot(1,1,1)
            PFlux = neoclassdata[ineoclassf][ispecs]['PFlux']*area
           #axhand.plot(time[bgnind:],PFlux[bgnind:],linestyle='-', label='$\\Gamma_{neo,%s}$=%9.7E (MW/keV)' % (ispecs,PFlux[-1]))
           #PFlux*= 1.0e6/(1.0e3*1.0e-19)
            axhand.plot(time[bgnind:],PFlux[bgnind:],linestyle='-', label='$\\Gamma_{neo,%s}$=%9.7E (particles/s)' % (ispecs,PFlux[-1]))
            axhand.set_title('$\\Gamma_{neo}%s$' %(titletxt))
            axhand.set_xlabel('Time (seconds)')
           #axhand.set_ylabel('Particle Flux (MW/keV)')
            axhand.set_ylabel('Particle Flux (particles/s)')
            if logplots: axhand.set_yscale('symlog')
            axhand.legend()


            HFluxfig = plt.figure('HFlux-'+ineoclassfpath)
            axhand = HFluxfig.add_subplot(1,1,1)
            HFlux = neoclassdata[ineoclassf][ispecs]['HFlux']*area
            HFlux/= 1.0e6
            axhand.plot(time[bgnind:],HFlux[bgnind:],linestyle='-', label='$Q_{neo,%s}$=%9.7E (MW)' % (ispecs,HFlux[-1]))
            axhand.set_title('$Q_{neo}%s$' %(titletxt))
            axhand.set_xlabel('Time (seconds)')
            axhand.set_ylabel('Heat Flux (MW)')
            if logplots: axhand.set_yscale('symlog')
            axhand.legend()

            Viscosfig = plt.figure('Viscos-'+ineoclassfpath)
            axhand = Viscosfig.add_subplot(1,1,1)
            Viscos = neoclassdata[ineoclassf][ispecs]['Viscos']*area
            axhand.plot(time[bgnind:],Viscos[bgnind:],linestyle='-', label='$\\Pi_{neo,%s}$=%9.7E (N.m)' % (ispecs,Viscos[-1]))
            axhand.set_title('$\\Pi_{neo}%s$' %(titletxt))
            axhand.set_xlabel('Time (seconds)')
            axhand.set_ylabel('Stress Tensor (N.m)')
            if logplots: axhand.set_yscale('symlog')
            axhand.legend()

            JBSfig = plt.figure('JBS-'+ineoclassfpath)
            axhand = JBSfig.add_subplot(1,1,1)
            JBS = neoclassdata[ineoclassf][ispecs]['JBS']
            axhand.plot(time[bgnind:],JBS[bgnind:],linestyle='-', label='$J_{bs,neo,%s}$=%9.7E (A)' % (ispecs,JBS[-1]))
            axhand.set_title('$J_{bs,neo}%s$' %(titletxt))
            axhand.set_xlabel('Time (seconds)')
            axhand.set_ylabel('Booststrap Current (Ampere)')
            if logplots: axhand.set_yscale('symlog')
            axhand.legend()

        figlist = []

        if display: plt.show()

        figlist.append(PFluxfig)
        PFluxfig.savefig(reportpath+'PFlux_%s.png' % (ineoclassfpath))
        plt.close(PFluxfig)

        figlist.append(HFluxfig)
        HFluxfig.savefig(reportpath+'HFlux_%s.png' % (ineoclassfpath))
        plt.close(HFluxfig)

        figlist.append(Viscosfig)
        Viscosfig.savefig(reportpath+'Viscos_%s.png' % (ineoclassfpath))
        plt.close(Viscosfig)

        figlist.append(JBSfig)
        JBSfig.savefig(reportpath+'JBS_%s.png' % (ineoclassfpath))
        plt.close(JBSfig)

        return figlist


def plot_scandata(scandata,params={},normalize=True):
   #Developed by Ehab Hassan on 2019-02-05
    #Modified by Ehab Hassan on 2019-03-11
    slashinds=findall(scandata['scanpath'],"/")
    scan_title = scandata['scanpath'][slashinds[-2]+1:slashinds[-1]]
    if not params.keys():
       params = genetools.read_parameters(scandata['scanpath'])
    vlist = list(set(scandata.keys())-{'omega','gamma','scanpath'})
    if   len(vlist) == 1:
         if   'kymin' in vlist:
              if   'x0' in params['box']:
                   plabel = '$x_0$='+str(params['box']['x0'])
              elif 'flux_pos' in params['geometry']:
                   plabel = '$x_0$='+str(params['geometry']['flux_pos'])
              else:
                   plabel = '$x_0$=?'
              if 'kx_center' in params['box'].keys():
                 plabel += ", $k_{x_center}$="+str(params['box']['kx_center'])
         else:
              plabel = '$\\rho_{ref}k_y$ = '+str(rhoref*params['box']['kymin'])
         gammafig = [plt.figure(1)]
         omegafig = [plt.figure(2)]
         axhand01 = gammafig[0].add_subplot(1,1,1)
         axhand02 = omegafig[0].add_subplot(1,1,1)
         if normalize:
              axhand01.plot(scandata[vlist[0]],scandata['gamma'],marker='*',label=plabel)
              axhand02.plot(scandata[vlist[0]],scandata['omega'],marker='*',label=plabel)
         elif not normalize:
              normfactor = params['units']['Lref']/npy.sqrt(params['units']['Tref']/params['units']['mref'])
              axhand01.plot(scandata[vlist[0]],normfactor*scandata['gamma'],marker='*',label=plabel)
              axhand02.plot(scandata[vlist[0]],normfactor*scandata['omega'],marker='*',label=plabel)
         if vlist[0] == 'kymin':
            axhand01.set_title('$k_y$ vs $\gamma$')
            axhand02.set_title('$k_y$ vs $\omega$')
            axhand01.set_xlabel("$\\rho_{ref}k_y$")
            axhand02.set_xlabel("$\\rho_{ref}k_y$")
         elif vlist[0] == 'kx_center':
            axhand01.set_title('$k_{x_{center}}$ vs $\gamma$')
            axhand02.set_title('$k_{x_{center}}$ vs $\omega$')
            axhand01.set_xlabel("$\\rho_{ref}k_{x_{center}}$")
            axhand02.set_xlabel("$\\rho_{ref}k_{x_{center}}$")
         else:
            axhand01.set_title(vlist[0]+' vs $\gamma$')
            axhand02.set_title(vlist[0]+' vs $\omega$')
            axhand01.set_xlabel(vlist[0])
            axhand02.set_xlabel(vlist[0])
         axhand01.set_ylabel('$\gamma/\Omega_{ref}$')
         axhand02.set_ylabel('$\omega/\Omega_{ref}$')
         axhand01.legend()
         axhand02.legend()
    elif len(vlist) == 2:
         if   'kymin' not in vlist:
              iscan = 'kymin = '+str(params['box']['kymin'])+': '
         else:
              iscan = ''
         gammafig = []
         omegafig = []
         maxrepeat = 1
         for var in vlist:
             for irepeat in range(1,len(scandata[var])):
                 if scandata[var][irepeat] == scandata[var][0]: break
             if irepeat >= maxrepeat: maxrepeatvar = var
             maxrepeat = max(maxrepeat,irepeat)
         gammafig.append(plt.figure(1))
         omegafig.append(plt.figure(2))
         axhand01 = gammafig[0].add_subplot(1,1,1)
         axhand02 = omegafig[0].add_subplot(1,1,1)
         for iplot in range(npy.size(scandata[vlist[0]])/maxrepeat):
             sindex = iplot*maxrepeat
             eindex = (iplot+1)*maxrepeat
             for var in list(set(vlist)-{maxrepeatvar}):
                 plabel = iscan+var+'='+str(scandata[var][sindex:eindex][0])+' '
             axhand01.plot(scandata[maxrepeatvar][sindex:eindex],scandata['gamma'][sindex:eindex],marker='*',label=plabel)
             axhand02.plot(scandata[maxrepeatvar][sindex:eindex],scandata['omega'][sindex:eindex],marker='*',label=plabel)
         axhand01.set_title(maxrepeatvar+' vs $\gamma$')
         axhand02.set_title(maxrepeatvar+' vs $\omega$')
         axhand01.set_xlabel(maxrepeatvar)
         axhand02.set_xlabel(maxrepeatvar)
         axhand01.set_ylabel('$\gamma$')
         axhand02.set_ylabel('$\omega$')
         axhand01.legend()
         axhand02.legend()

         gammafig.append(plt.figure(3))
         omegafig.append(plt.figure(4))
         pvarname = list(set(vlist)-{maxrepeatvar})
         axhand01 = gammafig[1].add_subplot(1,1,1)
         axhand02 = omegafig[1].add_subplot(1,1,1)
         for iplot in range(maxrepeat):
             sindex = iplot
             stprng = npy.size(scandata[pvarname[0]])/maxrepeat
             plabel = iscan+maxrepeatvar+'='+str(scandata[maxrepeatvar][sindex:eindex][0])+' '
             axhand01.plot(scandata[pvarname[0]][sindex::stprng],scandata['gamma'][sindex::stprng],marker='*',label=plabel)
             axhand02.plot(scandata[pvarname[0]][sindex::stprng],scandata['omega'][sindex::stprng],marker='*',label=plabel)
         axhand01.set_title(pvarname[0]+' vs $\gamma$')
         axhand02.set_title(pvarname[0]+' vs $\omega$')
         axhand01.set_xlabel(pvarname[0])
         axhand02.set_xlabel(pvarname[0])
         axhand01.set_ylabel('$\gamma$')
         axhand02.set_ylabel('$\omega$')
         axhand01.legend()
         axhand02.legend()

    return gammafig,omegafig

def plot_scans(scanfpaths,geneparams={},normalize=True):
   #Developed by Ehab Hassan on 2019-02-07
    #Modified by Ehab Hassan on 2019-03-07
    if type(scanfpaths)==list:
       for iscan in scanfpaths:
           scanvals = genetools.read_scanfile(iscan)
           gammafig,omegafig = plot_scandata(scandata=scanvals,params=geneparams,normalize=normalize)
       reportpath = "./"
    else:
           scanvals = genetools.read_scanfile(scanfpaths)
           gammafig,omegafig = plot_scandata(scandata=scanvals,params=geneparams,normalize=normalize)
           if scanfpaths[-1] != "/": scanfpaths += "/"
           reportpath = scanfpaths+"report/"
           if not os.path.isdir(reportpath):
              os.system('mkdir '+reportpath)

    pdfpages = pdfh.PdfPages(reportpath+'scanfigs.pdf')
    for item in range(len(gammafig)):
        gammafig[item].savefig(reportpath+'gamma%02d.png' % (item+1))
        omegafig[item].savefig(reportpath+'omega%02d.png' % (item+1))
        pdfpages.savefig(gammafig[item])
        pdfpages.savefig(omegafig[item])
    pdfpages.close()

    return gammafig,omegafig


from ParIO import *
from momlib import *
from fieldlib import *
from finite_differences import *
def my_corr_func_complex(v1,v2,time,show_plot=False,v1eqv2=True):
   #Developed by David Hatch on ????-??-??
    dt=time[1]-time[0]
    N=len(time)
    cfunc=npy.zeros(N,dtype='complex')
    for i in range(N):
        i0=i+1
        cfunc[-i0]=npy.sum(npy.conj(v1[-i0:])*v2[:i0])
    tau=npy.arange(N)
    tau=tau*dt
    if v1eqv2:
        cfunc=npy.real(cfunc)
    max_corr=max(npy.abs(cfunc))
    corr_time=0.0
    i=0
    while corr_time==0.0:
        if (abs(cfunc[i])-max_corr/npy.e) > 0.0 and \
           (abs(cfunc[i+1])-max_corr/npy.e) <= 0.0:
            slope=(cfunc[i+1]-cfunc[i])/(tau[i+1]-tau[i])
            zero=cfunc[i]-slope*tau[i]
            corr_time=(max_corr/npy.e-zero)/slope
        i+=1
    neg_loc = 10000.0
    i=0
    while neg_loc==10000.0 and i < N:
        if cfunc[i] < 0.0:
            neg_loc = tau[i]
        i+=1

    if neg_loc < corr_time:
        print "WARNING: neg_loc < corr_time"
        corr_time = neg_loc

    if show_plot:
        plt.plot(tau,cfunc,'x-')
        ax=plt.axis()
        plt.vlines(corr_time,ax[2],ax[3])
        plt.show()
    return cfunc,tau,corr_time

def plot_mom(mom,param={},reportpath=''):
   #Developed by Ehab Hassan on 2019-03-14
    imomf=mom.keys()[0]
    if reportpath == '':
       if not os.path.isdir(imomf[:-10]+"report"):
          os.system('mkdir '+imomf[:-10]+"report")
          reportpath = imomf[:-10]+"report/"
       else:
          reportpath = imomf[:-10]+"report/"
    elif reportpath[-1] != "/":
       reportpath += "/"

    for imomf in mom:
        if not param.keys():
           param = genetools.read_parameters(imomf[:-10]+"parameters_"+imomf[-4:])

        if 'x_local' in param['general']:
            if param['general']['x_local']:
                x_local = True
            else:
                x_local = False
        else:
            x_local = True

        if x_local:
           dens  = mom[imomf]['dens']
           tpar  = mom[imomf]['tpar'] 
           tperp = mom[imomf]['tperp']
           qpar  = mom[imomf]['qpar'] 
           qperp = mom[imomf]['qperp']
           upar  = mom[imomf]['upar'] 

           nx = param['box']['nx0']
           ny = param['box']['nky0']
           nz = param['box']['nz0']

           if 'lx_a' in param['box']:
                xgrid = param['box']['x0']+npy.arange(nx)/float(nx-1)*param['box']['lx_a']-param['box']['lx_a']/2.0
           elif 'lx' in param['box']:
                xgrid = npy.arange(nx)/float(nx-1)*param['box']['lx']-param['box']['lx']/2.0
           elif 'adapt_lx' in param['box']:
                lx = 1.0/(param['box']['kymin']*param['geometry']['shat'])
                xgrid = npy.arange(nx)/float(nx-1)*lx-lx/2.0
           zgrid = npy.arange(nz)/float(nz-1)*(2.0-(2.0/nz))-1.0

           plt.contourf(xgrid,zgrid,npy.abs(npy.fft.fftshift(dens[:,0,:])))
           plt.show()
           plt.contourf(xgrid,zgrid,npy.abs(npy.fft.fftshift(tpar[:,0,:])))
           plt.show()
           plt.contourf(xgrid,zgrid,npy.abs(npy.fft.fftshift(tperp[:,0,:])))
           plt.show()

    return

def plot_field(field,param={},reportpath='',setParam={}):
   #Developed by Ehab Hassan on 2019-03-14
    ifieldf=field.keys()[0]
    if 'report' not in reportpath:
       if not os.path.isdir(ifieldf[:-10]+"report"):
          os.system('mkdir '+ifieldf[:-10]+"report")
          reportpath = ifieldf[:-10]+"report/"
       else:
          reportpath = ifieldf[:-10]+"report/"
    elif reportpath[-1] != "/":
       reportpath += "/"

    if 'display' in setParam: display = True
    else:                     display = False

    for ifieldf in field:
        if not param.keys():
           if 'dat' in ifieldf:
              paramfpath = ifieldf[:-9]+"parameters"+ifieldf[-4:]
           else:
              paramfpath = ifieldf[:-10]+"parameters"+ifieldf[-5:]
           param = genetools.read_parameters(paramfpath)

        if 'x_local' in param['general']:
            if param['general']['x_local']:
                x_local = True
            else:
                x_local = False
        else:
            x_local = True

        if x_local:
           nx      = field[ifieldf]['nx']
           ny      = field[ifieldf]['ny']
           nz      = field[ifieldf]['nz']
           phi     = field[ifieldf]['phi']
           apar    = field[ifieldf]['apar']
           nfields = field[ifieldf]['nfields']
            
           zgrid = npy.arange(nx*nz)/float(nx*nz-1)*(2.0*nx-(2.0/nz))-nx

           phi1d  = npy.zeros(nx*nz,dtype='complex128')
           apar1d = npy.zeros(nx*nz,dtype='complex128')

           if 'n0_global' in param['box']:
               phase = -npy.e**(-2.0*npy.pi*(0.0+1.0J)*param['box']['n0_global']*param['geometry']['q0'])
           else:
               phase = -npy.e**(-npy.pi*(0.0+1.0J)*param['geometry']['shat']*param['box']['kymin']*param['box']['lx'])

           shatsgn = int(npy.sign(param['geometry']['shat']))
           for i in range(nx/2):
               phi1d[(i+nx/2)*nz:(i+nx/2+1)*nz]=phi[:,0,i*shatsgn]*phase**i
               if i < nx/2:
                   phi1d[(nx/2-i-1)*nz:(nx/2-i)*nz]=phi[:,0,-(i+1)*shatsgn]*phase**(-(i+1))
               if int(nfields)>1:
                  apar1d[(i+nx/2)*nz:(i+nx/2+1)*nz]=apar[:,0,i*shatsgn]*phase**i
                  if i < nx/2:
                       apar1d[(nx/2-i-1)*nz:(nx/2-i)*nz]=apar[:,0,-(i+1)*shatsgn]*phase**(-(i+1))

           phi1d  = phi1d/phi[nz/2,0,0]
           apar1d = apar1d/apar[nz/2,0,0]

           phinds = (abs(phi1d)>=1.0e-4)
           PHIfig = plt.figure('Phi_'+ifieldf[-4:])
           axhand = PHIfig.add_subplot(1,1,1)
           axhand.plot(zgrid[phinds],npy.real(phi1d[phinds]),color='red',label=r'$Re[\phi]$')
           axhand.plot(zgrid[phinds],npy.imag(phi1d[phinds]),color='blue',label=r'$Im[\phi]$')
           axhand.plot(zgrid[phinds],npy.abs(phi1d[phinds]),color='black',label=r'$|\phi|$')
           axhand.set_title(r'$\phi(k_y=%1.3f)$' % float(param['box']['kymin']) )
           axhand.set_xlabel(r'$z/\pi$',size=18)
           axhand.legend()

           aparinds = (abs(apar1d)>=1.0e-4)
           APARfig = plt.figure('Apar_'+ifieldf[-4:])
           axhand = APARfig.add_subplot(1,1,1)
           axhand.plot(zgrid[aparinds],npy.real(apar1d[aparinds]),color='red',label=r'$Re[A_{||}]$')
           axhand.plot(zgrid[aparinds],npy.imag(apar1d[aparinds]),color='blue',label=r'$Im[A_{||}]$')
           axhand.plot(zgrid[aparinds],npy.abs(apar1d[aparinds]),color='black',label=r'$|A_{||}|$')
           axhand.set_title(r'$A_{||}(k_y=%1.3f)$' % float(param['box']['kymin']))
           axhand.set_xlabel(r'$z/\pi$',size=18)
           axhand.legend()

           if 'dat' in ifieldf:
              omegafpath = ifieldf[:-9]+'omega'+ifieldf[-4:]
           else:
              omegafpath = ifieldf[:-10]+'omega'+ifieldf[-5:]

           if os.path.isfile(omegafpath):
               om = np.genfromtxt(omegafpath)
           omega_complex = (om[2]*(0.0+1.0J) + om[1])

          ##Note:  the complex frequency is (gamma + i*omega)
           if 'dat' in ifieldf:
              geomfpath = ifieldf[:-9]+param['geometry']['magn_geometry']+ifieldf[-4:]
           else:
              geomfpath = ifieldf[:-10]+param['geometry']['magn_geometry']+ifieldf[-5:]
           gpars,geometry = read_geometry_local(geomfpath)
           jacxB = geometry['gjacobian']*geometry['gBfield']

           gradphi = fd_d1_o4(phi1d,zgrid)
           for i in range(int(param['box']['nx0'])):
               gradphi[int(param['box']['nz0'])*i:int(param['box']['nz0'])*(i+1)] = gradphi[int(param['box']['nz0'])*i:int(param['box']['nz0'])*(i+1)]/jacxB[:]/npy.pi

           genlist = list(zip(abs(npy.real(gradphi))>=1.0e-6,abs(npy.real(omega_complex*apar1d))>=1.0e-6,abs(npy.imag(gradphi))>=1.0e-6,abs(npy.imag(omega_complex*apar1d))>=1.0e-6))
           geninds = [(i or j or k or l) for (i,j,k,l) in genlist]
           wAPARfig = plt.figure('wApar_dPhi_'+ifieldf[-4:])
           axhand = wAPARfig.add_subplot(1,1,1)
           axhand.plot(zgrid[geninds],npy.real(gradphi)[geninds],'-',color = 'red',label=r'$Re[\nabla \phi]$')
           axhand.plot(zgrid[geninds],npy.imag(gradphi)[geninds],'-.',color = 'red',label=r'$Im[\nabla \phi]$')
           axhand.plot(zgrid[geninds],-npy.real(omega_complex*apar1d)[geninds],'-',color = 'black',label=r'$Re[\omega A_{||}]$')
           axhand.plot(zgrid[geninds],-npy.imag(omega_complex*apar1d)[geninds],'-.',color = 'black',label=r'$Im[\omega A_{||}]$')
           axhand.set_title(r'$\nabla\phi,\partial_tA_{||}(k_y=%1.3f)$' % float(param['box']['kymin']))
           axhand.set_xlabel(r'$z/\pi$',size=18)
           axhand.legend()

           if display: plt.show()

           PHIfig.savefig(reportpath+'phi_mode_%s.png' % (ifieldf[-4:]))
           plt.close(PHIfig)
           APARfig.savefig(reportpath+'apar_mode_%s.png' % (ifieldf[-4:]))
           plt.close(APARfig)
           wAPARfig.savefig(reportpath+'wApar_dPhi_mode_%s.png' % (ifieldf[-4:]))
           plt.close(wAPARfig)

        elif not x_local:
           nx      = field[ifieldf]['nx']
           ny      = field[ifieldf]['ny']
           nz      = field[ifieldf]['nz']
           phi     = field[ifieldf]['phi']
           apar    = field[ifieldf]['apar']

           zgrid = npy.arange(nz+4)/float(nz+4-1)*(2.0+3.0*(2.0/nz))-(1.0+2.0*(2.0/nz))
           xgrid = npy.arange(nx)/float(nx-1)*param['box']['lx_a']+param['box']['x0']-param['box']['lx_a']/2.0

           if 'dat' in ifieldf:
              geomfpath = ifieldf[:-9]+param['geometry']['magn_geometry']+ifieldf[-4:]
           else:
              geomfpath = ifieldf[:-10]+param['geometry']['magn_geometry']+ifieldf[-5:]

           gpars,geometry = read_geometry_global(geomfpath)
           #Find rational q surfaces
           qmin = npy.min(geometry['q'])
           qmax = npy.max(geometry['q'])
           mmin = math.ceil( qmin*param['box']['n0_global'])
           mmax = math.floor(qmax*param['box']['n0_global'])
           mnums = npy.arange(mmin,mmax+1)
           qrats = mnums/float(param['box']['n0_global'])
           zgridm = np.arange(mmax*20)/float(mmax*20)*2.0-1.0
           nm = int(mmax*20)

           phase = (0.0+1.0J)*param['box']['n0_global']*2.0*npy.pi*geometry['q']

           phi_bnd = npy.zeros((nz+4,ny,nx),dtype = 'complex128')
           phi_bnd[2:-2,:,:] = phi
           for j in range(nx):
               phi_bnd[-2,0,j] = phi_bnd[ 2,0,j]*npy.e**(-phase[j])
               phi_bnd[-1,0,j] = phi_bnd[ 3,0,j]*npy.e**(-phase[j])
               phi_bnd[ 0,0,j] = phi_bnd[-4,0,j]*npy.e**( phase[j])
               phi_bnd[ 1,0,j] = phi_bnd[-3,0,j]*npy.e**( phase[j])

           gradphi= npy.zeros((nz+4,ny,nx),dtype = 'complex128')     
           for i in range(nx):
               gradphi[:,0,i] = fd_d1_o4(phi_bnd[:,0,i],zgrid)
               gradphi[2:-2:,0,i] = gradphi[2:-2,0,i]/npy.pi/(geometry['jacobian'][:,i]*geometry['Bfield'][:,i])

           phi_theta = np.zeros((nm,ny,nx),dtype = 'complex128')     
           for i in range(len(xgrid)):
              #phi_theta[:,0,i] = interp(zgrid,phi_bnd[:,0,i],zgridm)
              #Check using the real, imaginary, or absolute value in the interpolation
               phi_theta[:,0,i] = interp(zgrid,npy.real(phi_bnd[:,0,i]),zgridm)
               phi_theta[:,0,i] = phi_theta[:,0,i]*npy.e**(1J*float(param['box']['n0_global'])*geometry['q'][i]*npy.pi*zgridm)

           phi_m = npy.zeros((nm,ny,nx),dtype = 'complex128')     
           for i in range(len(xgrid)):
               phi_m[:,0,i] = npy.fft.fft(phi_theta[:,0,i])

           imax = np.unravel_index(np.argmax(abs(phi)),(nz,nx))
           plot_ballooning = True
           if plot_ballooning:
              PHI2Dfig = plt.figure('Phi2d_'+ifieldf[-4:])
              axhand = PHI2Dfig.add_subplot(3,1,1)
              cp=axhand.contourf(xgrid,zgrid,npy.abs(phi_bnd[:,0,:]),70)
              for i in range(len(qrats)):
                  ix = np.argmin(abs(geometry['q']-qrats[i])) 
                  axhand.axvline(xgrid[ix],color='white')
              axhand.plot(xgrid[imax[1]],zgrid[imax[0]],'x')
              cbar=PHI2Dfig.colorbar(cp)
              cbar.set_label(r'$|\phi|$')
              axhand.set_ylabel(r'$z/\pi$',fontsize=13)
              axhand.set_title(r'Electric Potential ($\phi$)',fontsize=13)
              axhand = PHI2Dfig.add_subplot(3,1,2)
              cp=axhand.contourf(xgrid,zgrid,npy.real(phi_bnd[:,0,:]),70)
              cbar=PHI2Dfig.colorbar(cp)
              cbar.set_label(r'$Re[\phi]$')
              axhand.set_ylabel(r'$z/\pi$',fontsize=13)
              axhand = PHI2Dfig.add_subplot(3,1,3)
              cp=axhand.contourf(xgrid,zgrid,npy.imag(phi_bnd[:,0,:]),70)
              cbar=PHI2Dfig.colorbar(cp)
              cbar.set_label(r'$Im[\phi]$')
              axhand.set_ylabel(r'$z/\pi$',fontsize=13)
              axhand.set_xlabel(r'$\rho_{tor}$',fontsize=13)

#          if plot_ballooning:
#              plt.xlabel(r'$z/\pi$',fontsize=13)
#              plt.ylabel(r'$|\phi|$')
#              plt.plot(zgrid,npy.abs(phi_bnd[:,0,param['box']['nx0']/4]),label='nx0/4')
#              plt.plot(zgrid,npy.abs(phi_bnd[:,0,param['box']['nx0']/2]),label='nx0/2')
#              plt.plot(zgrid,npy.abs(phi_bnd[:,0,3*param['box']['nx0']/4]),label='3/4*nx0')
#              plt.legend()
#              plt.show()

           plot_theta = True
           if plot_theta:
              PHI1Dfig = plt.figure('Phi1d_'+ifieldf[-4:])
              axhand = PHI1Dfig.add_subplot(1,1,1)
              for i in range(int(mmin),int(mmax)+1):
                  axhand.plot(xgrid,npy.abs(phi_m[i,0,:]))
              for i in range(len(qrats)):
                  ix = npy.argmin(abs(geometry['q']-qrats[i])) 
                  axhand.axvline(xgrid[ix],color='black')
              axhand.set_xlabel(r'$\rho_{tor}$',fontsize=13)
              axhand.set_ylabel(r'$\phi_{m}$',fontsize=13)
              axhand.set_title(r'$\phi_m$')

           zgrid = np.arange(nz)/float(nz-1)*(2.0-(2.0/nz))-1.0

           apar_theta = npy.zeros((nm,nx),dtype = 'complex128')     
           for i in range(len(xgrid)):
              #apar_theta[:,i] = interp(zgrid,apar[:,0,i],zgridm)
              #Check using the real, imaginary, or absolute value in the interpolation
               apar_theta[:,i] = interp(zgrid,npy.real(apar[:,0,i]),zgridm)
               apar_theta[:,i] = apar_theta[:,i]*npy.e**(1J*float(param['box']['n0_global'])*geometry['q'][i]*npy.pi*zgridm)

           apar_m     = npy.zeros((nm,nx),dtype = 'complex128')     
           for i in range(len(xgrid)):
               apar_m[:,i] = npy.fft.fft(apar_theta[:,i])

           if plot_theta:
              APAR1Dfig = plt.figure('Apar1d_'+ifieldf[-4:])
              axhand = APAR1Dfig.add_subplot(1,1,1)
              for i in range(int(mmin),int(mmax)+1):
                  axhand.plot(xgrid,npy.abs(apar_m[i,:]))
              for i in range(len(qrats)):
                  ix = npy.argmin(abs(geometry['q']-qrats[i])) 
                  axhand.axvline(xgrid[ix],color='black')
              axhand.set_xlabel(r'$\rho_{tor}$',fontsize=13)
              axhand.set_ylabel(r'$A_{||m}$',fontsize=13)
              axhand.set_title(r'$A_{||m}$')

           if plot_ballooning:
              APAR2Dfig = plt.figure('Apar2d_'+ifieldf[-4:])
              axhand = APAR2Dfig.add_subplot(3,1,1)
              cp=axhand.contourf(xgrid,zgrid,npy.abs(apar[:,0,:]),70)
              for i in range(len(qrats)):
                  ix = npy.argmin(abs(geometry['q']-qrats[i])) 
                  axhand.axvline(xgrid[ix],color='white')
              cbar=APAR2Dfig.colorbar(cp)
              cbar.set_label(r'$|A_{||}|$')
              axhand.set_ylabel(r'$z/\pi$',fontsize=13)
              axhand.set_title(r'Magnetic Potential ($A_{||}$)',fontsize=13)
              axhand = APAR2Dfig.add_subplot(3,1,2)
              cp=axhand.contourf(xgrid,zgrid,npy.real(apar[:,0,:]),70)
              cbar=APAR2Dfig.colorbar(cp)
              cbar.set_label(r'$Re[A_{||}]$')
              axhand.set_ylabel(r'$z/\pi$',fontsize=13)
              axhand = APAR2Dfig.add_subplot(3,1,3)
              cp=axhand.contourf(xgrid,zgrid,npy.imag(apar[:,0,:]),70)
              cbar=APAR2Dfig.colorbar(cp)
              cbar.set_label(r'$Im[A_{||}]$')
              axhand.set_ylabel(r'$z/\pi$',fontsize=13)
              axhand.set_xlabel(r'$\rho_{tor}$',fontsize=13)

           if 'dat' in ifieldf:
              omegafpath = ifieldf[:-9]+'omega'+ifieldf[-4:]
           else:
              omegafpath = ifieldf[:-10]+'omega'+ifieldf[-5:]

           if os.path.isfile(omegafpath):
               om = np.genfromtxt(omegafpath)

          ##Note:  the complex frequency is (gamma + i*omega)
           omega_complex = (om[2]*(0.0+1.0J) + om[1])

           if 'ExBrate' in param['external_contr'] and param['external_contr']['ExBrate'] == -1111: 
              if os.path.isfile(param['in_out']['iterdb_file']):
                 rhot_idb,profs_idb,units_idb = read_iterdb.read_iterdb(param['in_out']['iterdb_file'])
              else:
                 raise IOError('Iterdb file NOT FOUND in the given path!')
              if param['geometry']['geomdir'][-1]=='/':
                 prof_file = param['geometry']['geomdir'][:-6]+'PROFILES/'+param['geometry']['geomfile'][:-5]+'Profiles'
                 geom_file = param['geometry']['geomdir']+param['geometry']['geomfile']
              else:
                 prof_file = param['geometry']['geomdir'][:-5]+'PROFILES/'+param['geometry']['geomfile'][:-5]+'Profiles'
                 geom_file = param['geometry']['geomdir']+'/'+param['geometry']['geomfile']
              if os.path.isfile(prof_file) and os.path.isfile(geom_file):
                 profs,units = efittools.read_profiles_file(prof_file,setParam={'rhotor':param['box']['nx0'],'eqdskfpath':geom_file})
              else:
                 raise IOError('Profile and/or Geometry files NOT FOUND in the given path!')
              omegator0 = interp(rhot_idb['VROT'],profs_idb['VROT'],profs['rhotor'])
              mi = 1.673e-27
              ee = 1.602e-19
              mref = param['units']['mref']*mi
              time_ref = param['units']['Lref']/(param['units']['Tref']*1000.0*ee/mref)**0.5
              apar_cont = 0.0
              diff = 0.0
              apar_cont_2D = npy.empty(npy.shape(apar),dtype='complex128')
              for i in range(param['box']['nx0']):
                  diff += npy.sum(npy.abs(gradphi[2:-2,:,i] + (omega_complex+(0.0+1.0J)*param['box']['n0_global']*omegator0[i]*time_ref)*apar[:,:,i]))
                  apar_cont += npy.sum(npy.abs((omega_complex+(0.0+1.0J)*param['box']['n0_global']*omegator0[i]*time_ref)*apar[:,:,i]))
                  apar_cont_2D[:,:,i] =  (omega_complex+(0.0+1.0J)*param['box']['n0_global']*omegator0[i]*time_ref)*apar[:,:,i]
           else:
              diff = npy.sum(npy.abs(gradphi[2:-2,:,:] + omega_complex*apar[:,:,:]))
              apar_cont = npy.sum(npy.abs(omega_complex*apar[:,:,:]))
           phi_cont = np.sum(npy.abs(gradphi[2:-2,:,:]))

           if 'ExBrate' in param['external_contr'] and param['external_contr']['ExBrate'] == -1111 and plot_ballooning: 
              REPAR2Dfig = plt.figure('real_Epar'+ifieldf[-4:])
              axhand = REPAR2Dfig.add_subplot(3,1,1)
              cp=axhand.contourf(xgrid,zgrid,npy.real(gradphi[2:-2,0,:]),70,vmin = npy.min(np.real(gradphi[2:-2,0,:])),vmax = npy.max(np.real(gradphi[2:-2,0,:])))
              for i in range(len(qrats)):
                  ix = npy.argmin(abs(geometry['q']-qrats[i])) 
                  axhand.axvline(xgrid[ix],color='white')
              axhand.plot(xgrid[imax[1]],zgrid[imax[0]],'x')
              cbar=REPAR2Dfig.colorbar(cp)
              cbar.set_label(r'$Re(grad phi)$')
              axhand.set_ylabel(r'$z/\pi$',fontsize=13)
              axhand.set_title(r'Parallel Electric Field ($Re[E_{||}]$)',fontsize=13)
              axhand = REPAR2Dfig.add_subplot(3,1,2)
              cp=axhand.contourf(xgrid,zgrid,-npy.real(apar_cont_2D[:,0,:]),70,vmin = npy.min(np.real(gradphi[2:-2,0,:])),vmax = npy.max(np.real(gradphi[2:-2,0,:])))
              cbar= REPAR2Dfig.colorbar(cp)
              cbar.set_label(r'$Re[ omega Apar]$')
              axhand.set_ylabel(r'$z/\pi$',fontsize=13)
              axhand=REPAR2Dfig.add_subplot(3,1,3)
              cp=axhand.contourf(xgrid,zgrid,npy.real(gradphi[2:-2,0,:]+apar_cont_2D[:,0,:]),70,vmin = npy.min(npy.real(gradphi[2:-2,0,:])),vmax = npy.max(np.real(gradphi[2:-2,0,:])))
              cbar=REPAR2Dfig.colorbar(cp)
              cbar.set_label(r'$Re[Diff]$')
              axhand.set_ylabel(r'$z/\pi$',fontsize=13)
              axhand.set_xlabel(r'$\rho_{tor}$',fontsize=13)
           
              IEPAR2Dfig = plt.figure('imag_Epar'+ifieldf[-4:])
              axhand = IEPAR2Dfig.add_subplot(3,1,1)
              cp=axhand.contourf(xgrid,zgrid,npy.imag(gradphi[2:-2,0,:]),70,vmin = npy.min(npy.imag(gradphi[2:-2,0,:])),vmax = npy.max(npy.imag(gradphi[2:-2,0,:])))
              for i in range(len(qrats)):
                  ix = npy.argmin(abs(geometry['q']-qrats[i])) 
                  axhand.axvline(xgrid[ix],color='white')
              axhand.plot(xgrid[imax[1]],zgrid[imax[0]],'x')
              cbar=IEPAR2Dfig.colorbar(cp)
              cbar.set_label(r'$Im(grad phi)$')
              axhand.set_ylabel(r'$z/\pi$',fontsize=13)
              axhand.set_title(r'Parallel Electric Field ($Im[E_{||}]$)',fontsize=13)
              axhand = IEPAR2Dfig.add_subplot(3,1,2)
              cp=axhand.contourf(xgrid,zgrid,-npy.imag(apar_cont_2D[:,0,:]),70,vmin = npy.min(npy.imag(gradphi[2:-2,0,:])),vmax = npy.max(npy.imag(gradphi[2:-2,0,:])))
              cbar=IEPAR2Dfig.colorbar(cp)
              cbar.set_label(r'$Im[ omega Apar]$')
              axhand.set_ylabel(r'$z/\pi$',fontsize=13)
              axhand = IEPAR2Dfig.add_subplot(3,1,3)
              cp=axhand.contourf(xgrid,zgrid,npy.imag(gradphi[2:-2,0,:]+apar_cont_2D[:,0,:]),70,vmin = npy.min(npy.imag(gradphi[2:-2,0,:])),vmax = npy.max(npy.imag(gradphi[2:-2,0,:])))
              cbar=IEPAR2Dfig.colorbar(cp)
              cbar.set_label(r'$Im[Diff]$')
              axhand.set_ylabel(r'$z/\pi$',fontsize=13)
              axhand.set_xlabel(r'$\rho_{tor}$',fontsize=13)

              if display: plt.show()

              PHI2Dfig.savefig(reportpath+'phi_mode_%s_2d.png' % (ifieldf[-4:]))
              plt.close(PHI2Dfig)
              PHI1Dfig.savefig(reportpath+'phi_mode_%s.png' % (ifieldf[-4:]))
              plt.close(PHI1Dfig)
              APAR1Dfig.savefig(reportpath+'apar_mode_%s.png' % (ifieldf[-4:]))
              plt.close(APAR1Dfig)
              APAR2Dfig.savefig(reportpath+'apar_mode_%s_2d.png' % (ifieldf[-4:]))
              plt.close(APAR2Dfig)
              REPAR2Dfig.savefig(reportpath+'real_Epar_mode_%s_2d.png' % (ifieldf[-4:]))
              plt.close(REPAR2Dfig)
              IEPAR2Dfig.savefig(reportpath+'imag_Epar_mode_%s_2d.png' % (ifieldf[-4:]))
              plt.close(IEPAR2Dfig)

    return 1

def plot_geometry(geometryfpath,reportpath='',setParam={}):
    if type(geometryfpath)==str:
       gfpathlist=[geometryfpath]
    elif type(geometryfpath)==list:
        gfpathlist=geometryfpath[:]

    if 'dat' in geometryfpath:
       if   'chease' in geometryfpath:
            geomfprefix = geometryfpath[:-10]
       elif 's_alpha' in geometryfpath:
            geomfprefix = geometryfpath[:-11]
       elif 'tracer_efit' in geometryfpath:
            geomfprefix = geometryfpath[:-15]
    else:
       if   'chease' in geometryfpath:
            geomfprefix = geometryfpath[:-11]
       elif 's_alpha' in geometryfpath:
            geomfprefix = geometryfpath[:-12]
       elif 'tracer_efit' in geometryfpath:
            geomfprefix = geometryfpath[:-16]

    if 'report' not in reportpath:
       if not os.path.isdir(geomfprefix+"report"):
          os.system('mkdir '+geomfprefix+"report")
          reportpath = geomfprefix+"report/"
       else:
          reportpath = geomfprefix+"report/"
    elif reportpath[-1] != "/":
       reportpath += "/"

    if  'display' in setParam:
         display = setParam['display']
    else:
         display = False

    geomfigs = pdfh.PdfPages(reportpath+'geometry.pdf')

    for ifile in gfpathlist:
        try:
           parameters,geometry = read_geometry_local(ifile)
           x_local = True
        except ValueError:
           parameters,geometry = read_geometry_global(ifile)
           x_local = False

        Lref = parameters['Lref']
        nz = parameters['gridpoints']
        zgrid = npy.arange(nz)/float(nz-1)*(2.0-(2.0/nz))-1.0

        if ifile[-1]!="/": ifile+="/"
        slashinds=findall(ifile,"/")
        gfname=ifile[slashinds[-2]:slashinds[-1]]

        GGfig = plt.figure("ggCoeff",dpi=500)
        ax  = GGfig.add_subplot(2,3,1)
        if x_local:
           ax.plot(zgrid,geometry['ggxx'],label=gfname)
           ax.set_title('ggxx')
        else:
           ax.plot(zgrid,geometry['gxx'],label=gfname)
           ax.set_title('gxx')
        ax.set_xticks([])

        ax  = GGfig.add_subplot(2,3,2)
        if x_local:
           ax.plot(zgrid,geometry['ggyy'],label=gfname)
           ax.set_title('ggyy')
        else:
           ax.plot(zgrid,geometry['gyy'],label=gfname)
           ax.set_title('gyy')
        ax.set_xticks([])

        ax  = GGfig.add_subplot(2,3,3)
        if x_local:
           ax.plot(zgrid,geometry['ggzz'],label=gfname)
           ax.set_title('ggzz')
        else:
           ax.plot(zgrid,geometry['gzz'],label=gfname)
           ax.set_title('gzz')
        ax.set_xticks([])

        ax  = GGfig.add_subplot(2,3,4)
        if x_local:
           ax.plot(zgrid,geometry['ggxy'],label=gfname)
        else:
           ax.plot(zgrid,geometry['gxy'],label=gfname)
        ax.set_xlabel('$Z/\\pi$')
        ax.set_title('ggxy')

        ax  = GGfig.add_subplot(2,3,5)
        if x_local:
           ax.plot(zgrid,geometry['ggxz'],label=gfname)
        else:
           ax.plot(zgrid,geometry['gxz'],label=gfname)
        ax.set_xlabel('$Z/\\pi$')
        ax.set_title('ggxz')

        ax  = GGfig.add_subplot(2,3,6)
        if x_local:
           ax.plot(zgrid,geometry['ggyz'],label=gfname)
        else:
           ax.plot(zgrid,geometry['gyz'],label=gfname)
        ax.set_xlabel('$Z/\\pi$')
        ax.set_title('ggyz')

        if x_local:
           GLfig = plt.figure("gl_R",dpi=500)
           ax  = GLfig.add_subplot(2,2,1)
           ax.plot(zgrid,geometry['gl_R'],label=gfname)
           ax.set_title('gl_R')
           ax.set_xticks([])

           ax  = GLfig.add_subplot(2,2,2)
           ax.plot(zgrid,geometry['gl_z'],label=gfname)
           ax.set_title('gl_R')
           ax.set_xticks([])

           ax  = GLfig.add_subplot(2,2,3)
           ax.plot(zgrid,geometry['gl_dxdR'],label=gfname)
           ax.set_title('gl_dxdR')
           ax.set_xlabel('$Z/\\pi$')

           ax  = GLfig.add_subplot(2,2,4)
           ax.plot(zgrid,geometry['gl_dxdZ'],label=gfname)
           ax.set_title('gl_dxdZ')
           ax.set_xlabel('$Z/\\pi$')

        gBfig = plt.figure("gBfield",dpi=500)
        ax  = gBfig.add_subplot(2,2,1)
        if x_local:
           ax.plot(zgrid,geometry['gBfield'],label=gfname)
           ax.set_title('gBfield')
        else:
           ax.plot(zgrid,geometry['Bfield'],label=gfname)
           ax.set_title('Bfield')
        ax.set_xticks([])

        ax  = gBfig.add_subplot(2,2,2)
        if x_local:
           ax.plot(zgrid,geometry['gdBdx'],label=gfname)
           ax.set_title('gdBdx')
        else:
           ax.plot(zgrid,geometry['dBdx'],label=gfname)
           ax.set_title('dBdx')
        ax.set_xticks([])

        ax  = gBfig.add_subplot(2,2,3)
        if x_local:
           ax.plot(zgrid,geometry['gdBdy'],label=gfname)
           ax.set_title('gdBdy')
        else:
           ax.plot(zgrid,geometry['dBdy'],label=gfname)
           ax.set_title('dBdy')
        ax.set_xlabel('$Z/\\pi$')
        plt.legend()

        ax12  = gBfig.add_subplot(2,2,4)
        if x_local:
           ax12.plot(zgrid,geometry['gdBdz'],label=gfname)
           ax12.set_title('gdBdz')
        else:
           ax12.plot(zgrid,geometry['dBdz'],label=gfname)
           ax12.set_title('dBdz')
        ax12.set_xlabel('$Z/\\pi$')

    geomfigs.savefig(gBfig)
    gBfig.savefig(reportpath+'geometry_gB.png')
    geomfigs.savefig(GGfig)
    GGfig.savefig(reportpath+'geometry_GG.png')
    if x_local:
       geomfigs.savefig(GLfig)
       GLfig.savefig(reportpath+'geometry_GL.png')
    
    geomfigs.close()


    if display:
       plt.show()

    return 1

