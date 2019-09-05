#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import numpy as npy

from cheasepy  import *
from genetools import *
from geneplots import *
from efittools import *


def create_report(genefpath="./"):
   #Developed by Ehab Hassan on 2019-03-07
    if genefpath[-1] != "/": genefpath += "/"
    slashinds=findall(genefpath,"/")
    if len(slashinds) < 2:
       genefpath = "./"+genefpath
    slashinds=findall(genefpath,"/")

    report_title = genefpath[slashinds[-2]+1:slashinds[-1]]
    reportpath = genefpath+"report/"
    if not os.path.isdir(reportpath):
       os.system('mkdir '+reportpath)


    geneparams = read_parameters(genefpath)

    if 'x_local' in geneparams['general']:
       if geneparams['general']['x_local']:
          x_local = True
       else:
          x_local = False
    else:
       x_local = True
  

    if geneparams['geometry']['magn_geometry']=='tracer_efit':
      iterdb = True
      chease = False
      if 'iterdb_file' in geneparams['in_out'] and geneparams['in_out']['iterdb_file']!='':
         if os.path.isfile(geneparams['in_out']['iterdb_file']):
            iterdbpath = geneparams['in_out']['iterdb_file']
         else:
            if geneparams['in_out']['diagdir'][-1] != '/':
               iterdbpath = geneparams['in_out']['diagdir']+'/'
            else:
               iterdbpath = geneparams['in_out']['diagdir']
            inds = findall(iterdbpath,"/")
            iterdbpath = iterdbpath.replace(iterdbpath[inds[-2]:inds[-1]],'')
            if os.path.isfile(iterdbpath+geneparams['in_out']['iterdb_file']):
               iterdbpath += geneparams['in_out']['iterdb_file']
            else:
               inds = findall(iterdbpath,"/")
               iterdbpath = iterdbpath.replace(iterdbpath[inds[-2]:inds[-1]],'')
               if os.path.isfile(iterdbpath+geneparams['in_out']['iterdb_file']):
                  iterdbpath += geneparams['in_out']['iterdb_file']
               else:
                   iterdbpath = raw_input('iterdb file not found, enter full path for iterdb file: ')
                   if not os.path.isfile(iterdbpath):
                      print('iterdb not found in path provided. EXIT!')
                      sys.exit()
         iterdbdata = read_iterdb.read_iterdb(iterdbpath)
         plot_iterdb(iterdbdata,reportpath)
    elif geneparams['geometry']['magn_geometry']=='chease':
       iterdb = False
       chease = True
       cheasepath = geneparams['geometry']['geomdir']+geneparams['geometry']['geomfile']
       plot_cheasedata(OSPATH=cheasepath,reportpath=reportpath)


    print('Creating GENE Report ...')
    texfname = 'generesults.tex'
    texfhand = open(reportpath+texfname,'w')
    texfhand.write("\\documentclass[]{report} \n")
    texfhand.write("\\usepackage{float} \n")    
    texfhand.write("\\usepackage{placeins} \n")    
    texfhand.write("\\usepackage{graphicx} \n")    
    texfhand.write("\\usepackage{geometry} \n")    
    texfhand.write("\\geometry{legalpaper, landscape, margin=0.7in} \n")    
    texfhand.write("\\usepackage{makecell} \n")    
    texfhand.write("\\usepackage{multirow,tabularx} \n")    
    texfhand.write("\\setcellgapes{4pt} \n")    

    texfhand.write("\\title{GENE Scan Run Results\\\\- "+report_title.replace('_','\_')+" -} \n")
    texfhand.write("\\begin{document}\n")
    texfhand.write("\\maketitle \n")

    print('Reading GENE Parameter(s) files ...')
    texfhand.write("\\section{List of Parameters in %s:} \n" % (genefpath[slashinds[-2]+1:slashinds[-1]].replace('_','\_')))
    texfhand.write("\\begin{tabular}{ l  l  l }\n")
    texfhand.write("\\begin{tabular}{| p{2.5cm} | p{5.5cm} |}\n")
    texfhand.write("\\hline \n")
    texfhand.write("\\multicolumn{2}{|c|}{General}\\\\ \n")
    texfhand.write("\\hline \n")
    for subitem in geneparams['general'].keys():
        texfhand.write(subitem.replace('_','\_')+"&"+str(geneparams['general'][subitem]).replace('_','\_')+"\\\\ \n")
    texfhand.write("\\hline \n")
    if 'species1' in geneparams.keys():
       texfhand.write("\\multicolumn{2}{|c|}{Species-1}\\\\ \n")
       texfhand.write("\\hline \n")
       for subitem in geneparams['species1'].keys():
           texfhand.write(subitem.replace('_','\_')+"&"+str(geneparams['species1'][subitem]).replace('_','\_')+"\\\\ \n")
       texfhand.write("\\hline \n")
    texfhand.write("\\end{tabular} \n")

    texfhand.write("& \n")

    texfhand.write("\\begin{tabular}{| p{2.5cm} | p{5.5cm} |}\n")
    texfhand.write("\\hline \n")
    texfhand.write("\\multicolumn{2}{|c|}{Parallelization}\\\\ \n")
    texfhand.write("\\hline \n")
    for subitem in geneparams['parallelization'].keys():
        texfhand.write(subitem.replace('_','\_')+"&"+str(geneparams['parallelization'][subitem]).replace('_','\_')+"\\\\ \n")
    texfhand.write("\\hline \n")
    texfhand.write("\\multicolumn{2}{|c|}{Box}\\\\ \n")
    texfhand.write("\\hline \n")
    for subitem in geneparams['box'].keys():
        texfhand.write(subitem.replace('_','\_')+"&"+str(geneparams['box'][subitem]).replace('_','\_')+"\\\\ \n")
    texfhand.write("\\hline \n")
    if 'species2' in geneparams.keys():
       texfhand.write("\\multicolumn{2}{|c|}{Species-2}\\\\ \n")
       texfhand.write("\\hline \n")
       for subitem in geneparams['species2'].keys():
           texfhand.write(subitem.replace('_','\_')+"&"+str(geneparams['species2'][subitem]).replace('_','\_')+"\\\\ \n")
       texfhand.write("\\hline \n")
    texfhand.write("\\end{tabular} \n")

    texfhand.write("& \n")

    texfhand.write("\\begin{tabular}{| p{3cm} | p{10cm} |}\n")
    texfhand.write("\\hline \n")
    texfhand.write("\\multicolumn{2}{|c|}{Geometry}\\\\ \n")
    texfhand.write("\\hline \n")
    for subitem in geneparams['geometry'].keys():
        if subitem == "geomdir": continue
        texfhand.write(subitem.replace('_','\_')+"&"+str(geneparams['geometry'][subitem]).replace('_','\_')+"\\\\ \n")
    texfhand.write("\\hline \n")
    texfhand.write("\\multicolumn{2}{|c|}{In-Out}\\\\ \n")
    texfhand.write("\\hline \n")
    for subitem in geneparams['in_out'].keys():
        if subitem == "diagdir": continue
        texfhand.write(subitem.replace('_','\_')+"&"+str(geneparams['in_out'][subitem]).replace('_','\_')+"\\\\ \n")
    texfhand.write("\\hline \n")
    texfhand.write("\\multicolumn{2}{|c|}{Units}\\\\ \n")
    texfhand.write("\\hline \n")
    for subitem in geneparams['units'].keys():
        texfhand.write(subitem.replace('_','\_')+"&"+str(geneparams['units'][subitem]).replace('_','\_')+"\\\\ \n")
    texfhand.write("\\hline \n")
    if 'external_contr' in geneparams.keys():
       texfhand.write("\\multicolumn{2}{|c|}{External Profiles}\\\\ \n")
       texfhand.write("\\hline \n")
       for subitem in geneparams['external_contr'].keys():
           texfhand.write(subitem.replace('_','\_')+"&"+str(geneparams['external_contr'][subitem]).replace('_','\_')+"\\\\ \n")
       texfhand.write("\\hline \n")
    if 'species3' in geneparams.keys():
       texfhand.write("\\multicolumn{2}{|c|}{Species-3}\\\\ \n")
       texfhand.write("\\hline \n")
       for subitem in geneparams['species3'].keys():
           texfhand.write(subitem.replace('_','\_')+"&"+str(geneparams['species3'][subitem]).replace('_','\_')+"\\\\ \n")
       texfhand.write("\\hline \n")
    texfhand.write("\\end{tabular} \n")

    texfhand.write("\\end{tabular} \n")
    texfhand.write("\\\\ \n")

    print('Reading Profiles ...')
    texfhand.write("\\clearpage \n")
    texfhand.write("\\section{Temperature and Density Profiles:} \n")
    Te_x0_gn = npy.nan; Ne_x0_gn = npy.nan; Ze = npy.nan
    Ti_x0_gn = npy.nan; Ni_x0_gn = npy.nan; Zi = npy.nan
    Tz_x0_gn = npy.nan; Nz_x0_gn = npy.nan; Zz = npy.nan
    for ispecid in range(geneparams['box']['n_spec']):
        specname = 'species'+str(ispecid+1)
        if 'e' in geneparams[specname]['name']:
           if type(geneparams[specname]['temp'])==list:
              Te_x0_gn = max(geneparams[specname]['temp'])*max(geneparams['units']['Tref'])*1.0e3
              Ne_x0_gn = max(geneparams[specname]['dens'])*max(geneparams['units']['nref'])*1.0e19
           else:
              Te_x0_gn = geneparams[specname]['temp']*geneparams['units']['Tref']*1.0e3
              Ne_x0_gn = geneparams[specname]['dens']*geneparams['units']['nref']*1.0e19
           Ze       = geneparams[specname]['charge']
        if 'i' in geneparams[specname]['name']:
           if type(geneparams[specname]['temp'])==list:
              Ti_x0_gn = max(geneparams[specname]['temp'])*max(geneparams['units']['Tref'])*1.0e3
              Ni_x0_gn = max(geneparams[specname]['dens'])*max(geneparams['units']['nref'])*1.0e19
           else:
              Ti_x0_gn = geneparams[specname]['temp']*geneparams['units']['Tref']*1.0e3
              Ni_x0_gn = geneparams[specname]['dens']*geneparams['units']['nref']*1.0e19
           Zi       = geneparams[specname]['charge']
        if 'z' in geneparams[specname]['name']:
           if type(geneparams[specname]['temp'])==list:
              Tz_x0_gn = max(geneparams[specname]['temp'])*max(geneparams['units']['Tref'])*1.0e3
              Nz_x0_gn = max(geneparams[specname]['dens'])*max(geneparams['units']['nref'])*1.0e19
           else:
              Tz_x0_gn = geneparams[specname]['temp']*geneparams['units']['Tref']*1.0e3
              Nz_x0_gn = geneparams[specname]['dens']*geneparams['units']['nref']*1.0e19
           Zz       = geneparams[specname]['charge']
    Zeff_gn   = (Zi**2*Ni_x0_gn+Zz**2*Nz_x0_gn)/Ne_x0_gn
    tau_x0_gn = Zeff_gn*(Te_x0_gn/Ti_x0_gn)
    texfhand.write("\\begin{center}\n")
    texfhand.write("\\begin{tabular}{| p{4cm} | p{5cm} |}\n")
    texfhand.write("\\hline \n")
    if 'x0' in geneparams['box']:
       texfhand.write("\\multicolumn{2}{|c|}{Species Properties at $\\rho_{tor}$($x_0$) = "+str(geneparams['box']['x0'])+"}\\\\ \n")
    elif 'flux_pos' in geneparams['geometry']:
       texfhand.write("\\multicolumn{2}{|c|}{Species Properties at $\\rho_{tor}$($x_0$) = "+str(geneparams['geometry']['flux_pos'])+"}\\\\ \n")
    else:
       texfhand.write("\\multicolumn{2}{|c|}{Species Properties at $\\rho_{tor}$($x_0$) = ?}\\\\ \n")
    texfhand.write("\\hline \n")
    texfhand.write("$T_e (eV)$                                & "+str(Te_x0_gn)+"\\\\ \n")
    texfhand.write("$T_i (eV)$                                & "+str(Ti_x0_gn)+"\\\\ \n")
    texfhand.write("$N_e (/m^3)$                              & "+str(Ne_x0_gn)+"\\\\ \n")
    texfhand.write("$N_i (/m^3)$                              & "+str(Ni_x0_gn)+"\\\\ \n")
    texfhand.write("$N_z (/m^3)$                              & "+str(Nz_x0_gn)+"\\\\ \n")
    texfhand.write("$Z_{eff}=\\frac{Z_i^2n_i+Z_z^2*n_z}{n_e}$ & "+str(Zeff_gn)+"\\\\ \n")
    texfhand.write("$\\tau=Z_{eff}\\frac{T_e}{T_i}$           & "+str(tau_x0_gn)+"\\\\ \n")
    texfhand.write("\\hline \n")
    texfhand.write("\\end{tabular} \n")
    texfhand.write("\\end{center} \n")
    texfhand.write("\\begin{figure}[!ht]\n")
    texfhand.write("\\begin{center} \n")
    if iterdb:
       texfhand.write("\\includegraphics[scale=0.8]{"+reportpath+"iterdb_density.png} \n")
       texfhand.write("\\includegraphics[scale=0.8]{"+reportpath+"iterdb_temperature.png} \n")
    elif chease:
       texfhand.write("\\includegraphics[scale=0.8]{"+reportpath+"chease_density.png} \n")
       texfhand.write("\\includegraphics[scale=0.8]{"+reportpath+"chease_temperature.png} \n")
    texfhand.write("\\end{center} \n")
    texfhand.write("\\end{figure} \n")


    print('Reading Geometry ...')
    texfhand.write("\\clearpage \n")
    texfhand.write("\\section{Geometry:} \n")
    if 'geomfile' in geneparams['geometry'].keys():
       if geneparams['geometry']['geomdir'][-1] == '/':
          geomfpath = geneparams['geometry']['geomdir']+geneparams['geometry']['geomfile']
       else:
          geomfpath = geneparams['geometry']['geomdir']+"/"+geneparams['geometry']['geomfile']
    texfhand.write("\\begin{figure}[!ht] \n")
    texfhand.write("\\begin{center} \n")
    texfhand.write("\\begin{tabular}{c c} \n")
    if geneparams['geometry']['magn_geometry']=='tracer_efit':
       eqdskdata = read_eqdsk(geomfpath)
       plot_eqdsk(eqdskdata,genefpath)
       texfhand.write("\\includegraphics[scale=0.6]{"+reportpath+"eqdsk_safetyfactor.png}& \n")
       texfhand.write("\\includegraphics[scale=0.6]{"+reportpath+"eqdsk_magsurfbound.png}\\\\ \n")
       texfhand.write("\\includegraphics[scale=0.6]{"+reportpath+"eqdsk_pprime.png}& \n")
       texfhand.write("\\includegraphics[scale=0.6]{"+reportpath+"eqdsk_ffprime.png}\\ \n")
    elif geneparams['geometry']['magn_geometry']=='chease':
       plot_cheasedata(OSPATH=geomfpath,reportpath=genefpath)
       texfhand.write("\\includegraphics[scale=0.6]{"+reportpath+"chease_safetyfactor.png}& \n")
       texfhand.write("\\includegraphics[scale=0.6]{"+reportpath+"chease_magsurfbound.png}\\\\ \n")
       texfhand.write("\\includegraphics[scale=0.6]{"+reportpath+"chease_pprime.png}& \n")
       texfhand.write("\\includegraphics[scale=0.6]{"+reportpath+"chease_ttprime.png}\\ \n")
    texfhand.write("\\end{tabular} \n")
    texfhand.write("\\end{center} \n")
    texfhand.write("\\end{figure} \n")


    print('Reading Omega files ...')
    texfhand.write("\\clearpage \n")
    scanstatus = plot_scans(genefpath,geneparams)
    texfhand.write("\\section{Growth rate and Oscillation Frequency:} \n")
    texfhand.write("\\begin{figure}[!ht]\n")
    texfhand.write("\\begin{center} \n")
    texfhand.write("\\begin{tabular}{c c c} \n")
    texfhand.write("\\begin{minipage}{0.50\\textwidth}\n")
    texfhand.write("\\includegraphics[scale=0.8]{"+reportpath+"omega01.png} \n")
    texfhand.write("\\end{minipage} \n")
    texfhand.write("& \n")
    texfhand.write("\\bigskip \n")
    texfhand.write("& \n")
    texfhand.write("\\begin{minipage}{0.50\\textwidth}\n")
    texfhand.write("\\includegraphics[scale=0.8]{"+reportpath+"gamma01.png} \n")
    texfhand.write("\\end{minipage}\\\\ \n")
    texfhand.write("\\begin{minipage}{0.50\\textwidth}\n")
    texfhand.write("Oscillation frequency as a function of $k_y$ for the unstable modes excited in the tokamaks. \n")
    texfhand.write("\\end{minipage} \n")
    texfhand.write("& \n")
    texfhand.write("\\bigskip \n")
    texfhand.write("& \n")
    texfhand.write("\\begin{minipage}{0.50\\textwidth}\n")
    texfhand.write("Growth rate as a function of $k_y$ for the unstable modes excited in the tokamaks. \n")
    texfhand.write("\\end{minipage}\n")
    texfhand.write("\\end{tabular} \n")
    texfhand.write("\\end{center} \n")
    texfhand.write("\\end{figure} \n")

    texfhand.write("\\clearpage \n")
    texfhand.write("\\section{Modes Classification} \n")
    texfhand.write("\\subsection{Method I:} \n")
    print('Finding Modes Classification from Fields ...')
    modeinfo = mode_info(genefpath)
    texfhand.write("\\begin{table}[!h] \n")
    if geneparams['box']['n_spec'] >= 3:
       texfhand.write("\\begin{tabular}{| m{1.5cm} | m{1.5cm} | m{1.5cm} | m{1.5cm} | m{1.5cm} | m{2.0cm} | m{1.5cm} | m{1.5cm} | m{1.5cm} | m{2.0cm} | m{1.5cm} |} \n") 
       texfhand.write("\\hline \n")
       texfhand.write("$k_y$ & $\\frac{Q_{EM}}{Q_{ES}}$ & $C_{\\ell}(\\phi,\\phi)$ & $A_{{\\parallel}parity}$ & $\\phi_{parity}$ & ")
       texfhand.write("$\\frac{E_{\\parallel}}{|\\nabla\phi|+|\partial_tA_{\\parallel}|}$ & $Z_{avg}$ & $\\gamma$ & $\\omega$ & $|\\omega+j\\gamma|$ & Mode \\\\ \n")
       texfhand.write("\\hline \n")
       if type(geneparams['box']['kymin'])==list:
          kyminlist = geneparams['box']['kymin']
       else:
          kyminlist = [geneparams['box']['kymin']]
       for item in kyminlist:
           texfhand.write("%5.3f & \n" % (modeinfo[item]['kymin']))
           texfhand.write("%5.3f & \n" % (modeinfo[item]['Qem/Qes']))
           texfhand.write("%5.3f & \n" % (modeinfo[item]['corr_len']))
           texfhand.write("%5.3f & \n" % (modeinfo[item]['parity_factor_apar']))
           texfhand.write("%5.3f & \n" % (modeinfo[item]['parity_factor_phi']))
           texfhand.write("%5.3f & \n" % (modeinfo[item]['Epar_Cancellation']))
           texfhand.write("%5.3f & \n" % (modeinfo[item]['zavg']))
           texfhand.write("%5.3f & \n" % (modeinfo[item]['gamma']))
           texfhand.write("%5.3f & \n" % (modeinfo[item]['omega']))
           texfhand.write("%5.3f & \n" % (npy.sqrt(modeinfo[item]['gamma']**2+modeinfo[item]['omega']**2)))
           texfhand.write("%s \\\\ \n" % (modeinfo[item]['Type']))
           texfhand.write("\\hline \n")
       texfhand.write("\\end{tabular} \n")
    else:
       texfhand.write("\\begin{table}[!h] \n")
       texfhand.write("\\begin{tabular}{| m{1.5cm} | m{1.5cm} |} \n") 
       texfhand.write("\\hline \n")
       texfhand.write("$k_y$ & Mode \\\\ \n")
       texfhand.write("\\hline \n")
       for item in sorted(modetype.keys()):
           texfhand.write("%5.3f & \n" % (item))
           texfhand.write("%s \\\\ \n" % (modetype[item]))
           texfhand.write("\\hline \n")
       texfhand.write("\\end{tabular} \n")
    texfhand.write("\\end{table} \n")

    if x_local:
       texfhand.write("\\clearpage \n")
       texfhand.write("\\section{Modes Classification} \n")
       texfhand.write("\\subsection{Method II:} \n")
       print('Finding Modes Classification from Fluxes ...')
       fluxinfo = flux_info(genefpath)
       texfhand.write("\\begin{table}[!h] \n")
       if geneparams['box']['n_spec'] >= 3:
          texfhand.write("\\begin{tabular}{| m{2.0cm} | m{2.0cm} | m{2.0cm} | m{2.0cm} | m{2.0cm} | m{2.0cm} | m{2.0cm} | m{2.0cm} |} \n")
          texfhand.write("\\hline \n")
          texfhand.write("$k_y$ & $\\chi_i/\\chi_e$ & $\\chi_e/\\chi_i$ & ")
          texfhand.write("$D_e/\\chi_e$ & $D_z/\\chi_e$ & $D_e/(\\chi_e+\\chi_i)$ & $D_z/(\\chi_e+\\chi_i)$ & Mode \\\\ \n")
          texfhand.write("\\hline \n")
          if type(geneparams['box']['kymin'])==list:
             kyminlist = geneparams['box']['kymin']
          else:
             kyminlist = [geneparams['box']['kymin']]
          for item in kyminlist:
              texfhand.write("%5.3f & \n" % (item))
              texfhand.write("%5.3f & \n" % (fluxinfo[item]['i']['Chi']/fluxinfo[item]['e']['Chi']))
              texfhand.write("%5.3f & \n" % (fluxinfo[item]['e']['Chi']/fluxinfo[item]['i']['Chi']))
              texfhand.write("%5.3f & \n" % (fluxinfo[item]['e']['Dee']/fluxinfo[item]['e']['Chi']))
              texfhand.write("%5.3f & \n" % (fluxinfo[item]['z']['Dee']/fluxinfo[item]['e']['Chi']))
              texfhand.write("%5.3f & \n" % (fluxinfo[item]['e']['Dee']/(fluxinfo[item]['e']['Chi']+fluxinfo[item]['i']['Chi'])))
              texfhand.write("%5.3f & \n" % (fluxinfo[item]['z']['Dee']/(fluxinfo[item]['e']['Chi']+fluxinfo[item]['i']['Chi'])))
              texfhand.write("%s \\\\ \n" % (fluxinfo[item]['Type']))
              texfhand.write("\\hline \n")
          texfhand.write("\\end{tabular} \n")
       else:
          texfhand.write("\\begin{table}[!h] \n")
          texfhand.write("\\begin{tabular}{| m{1.5cm} | m{1.5cm} |} \n")
          texfhand.write("\\hline \n")
          texfhand.write("$k_y$ & Mode \\\\ \n")
          texfhand.write("\\hline \n")
          for item in sorted(modetype.keys()):
              texfhand.write("%5.3f & \n" % (item))
              texfhand.write("%s \\\\ \n" % (modetype[item]))
              texfhand.write("\\hline \n")
          texfhand.write("\\end{tabular} \n")
       texfhand.write("\\end{table} \n")


    print('Ploting Mode Structures ...')
    field = read_field(fieldfpath=genefpath)
    plot_field(field=field,reportpath=genefpath)
    if x_local:
       texfhand.write("\\clearpage \n")
       texfhand.write("\\section{Modes $\\phi(k_y)$ Structure} \n")
       if type(geneparams['box']['kymin'])==list:
          kyminlist = geneparams['box']['kymin']
       else:
          kyminlist = [geneparams['box']['kymin']]
       nmodes = len(kyminlist)
       for imode in range(0,nmodes,6):
           if (imode+1) % 6 == 0:
              texfhand.write("\\clearpage \n")
              texfhand.write("\\section*{Modes $\\phi(k_y)$ Structure (continue ...)} \n")
           texfhand.write("\\begin{figure}[!h] \n")
           texfhand.write("\\begin{center} \n")
           texfhand.write("\\begin{tabular}{c c c} \n")
           if imode+1 <= nmodes:
              texfhand.write("\\includegraphics[scale=0.65]{"+reportpath+"phi_mode_%04d.png} & \n"    % (imode+1))
           if imode+2 <= nmodes:
              texfhand.write("\\includegraphics[scale=0.65]{"+reportpath+"phi_mode_%04d.png} & \n"    % (imode+2))
           if imode+3 <= nmodes:
              texfhand.write("\\includegraphics[scale=0.65]{"+reportpath+"phi_mode_%04d.png} \\\\ \n" % (imode+3))
           if imode+4 <= nmodes:
              texfhand.write("\\includegraphics[scale=0.65]{"+reportpath+"phi_mode_%04d.png} & \n"    % (imode+4))
           if imode+5 <= nmodes:
              texfhand.write("\\includegraphics[scale=0.65]{"+reportpath+"phi_mode_%04d.png} & \n"    % (imode+5))
           if imode+6 <= nmodes:
              texfhand.write("\\includegraphics[scale=0.65]{"+reportpath+"phi_mode_%04d.png} \\\\ \n" % (imode+6))
           texfhand.write("\\end{tabular} \n")
           texfhand.write("\\end{center} \n")
           texfhand.write("\\end{figure} \n")

       texfhand.write("\\clearpage \n")
       texfhand.write("\\section{Modes $A_{||}(k_y)$ Structure} \n")
       for imode in range(0,nmodes,6):
           if (imode+1) % 6 == 0:
              texfhand.write("\\clearpage \n")
              texfhand.write("\\section*{Modes $A_{||}(k_y)$ Structure (continue ...)} \n")
           texfhand.write("\\begin{figure}[!h] \n")
           texfhand.write("\\begin{center} \n")
           texfhand.write("\\begin{tabular}{c c c} \n")
           if imode+1 <= nmodes:
              texfhand.write("\\includegraphics[scale=0.65]{"+reportpath+"apar_mode_%04d.png} & \n"    % (imode+1))
           if imode+2 <= nmodes:
              texfhand.write("\\includegraphics[scale=0.65]{"+reportpath+"apar_mode_%04d.png} & \n"    % (imode+2))
           if imode+3 <= nmodes:
              texfhand.write("\\includegraphics[scale=0.65]{"+reportpath+"apar_mode_%04d.png} \\\\ \n" % (imode+3))
           if imode+4 <= nmodes:
              texfhand.write("\\includegraphics[scale=0.65]{"+reportpath+"apar_mode_%04d.png} & \n"    % (imode+4))
           if imode+5 <= nmodes:
              texfhand.write("\\includegraphics[scale=0.65]{"+reportpath+"apar_mode_%04d.png} & \n"    % (imode+5))
           if imode+6 <= nmodes:
              texfhand.write("\\includegraphics[scale=0.65]{"+reportpath+"apar_mode_%04d.png} \\\\ \n" % (imode+6))
           texfhand.write("\\end{tabular} \n")
           texfhand.write("\\end{center} \n")
           texfhand.write("\\end{figure} \n")

       texfhand.write("\\clearpage \n")
       texfhand.write("\\section{Modes $\\nabla\phi(k_y),\partial_tA_{||}(k_y)$ Structure} \n")
       for imode in range(0,nmodes,6):
           if (imode+1) % 6 == 0:
              texfhand.write("\\clearpage \n")
              texfhand.write("\\section*{Modes $\\nabla\phi(k_y),\partial_tA_{||}(k_y)$ Structure (continue ...)} \n")
           texfhand.write("\\begin{figure}[!h] \n")
           texfhand.write("\\begin{center} \n")
           texfhand.write("\\begin{tabular}{c c c} \n")
           if imode+1 <= nmodes:
              texfhand.write("\\includegraphics[scale=0.65]{"+reportpath+"wApar_dPhi_mode_%04d.png} & \n"    % (imode+1))
           if imode+2 <= nmodes:
              texfhand.write("\\includegraphics[scale=0.65]{"+reportpath+"wApar_dPhi_mode_%04d.png} & \n"    % (imode+2))
           if imode+3 <= nmodes:
              texfhand.write("\\includegraphics[scale=0.65]{"+reportpath+"wApar_dPhi_mode_%04d.png} \\\\ \n" % (imode+3))
           if imode+4 <= nmodes:
              texfhand.write("\\includegraphics[scale=0.65]{"+reportpath+"wApar_dPhi_mode_%04d.png} & \n"    % (imode+4))
           if imode+5 <= nmodes:
              texfhand.write("\\includegraphics[scale=0.65]{"+reportpath+"wApar_dPhi_mode_%04d.png} & \n"    % (imode+5))
           if imode+6 <= nmodes:
              texfhand.write("\\includegraphics[scale=0.65]{"+reportpath+"wApar_dPhi_mode_%04d.png} \\\\ \n" % (imode+6))
           texfhand.write("\\end{tabular} \n")
           texfhand.write("\\end{center} \n")
           texfhand.write("\\end{figure} \n")
    if not x_local:
       texfhand.write("\\clearpage \n")
       texfhand.write("\\section{Modes $\\phi(k_y)$ Structure} \n")
       if type(geneparams['box']['kymin'])==list:
          kyminlist = geneparams['box']['kymin']
       else:
          kyminlist = [geneparams['box']['kymin']]
       nmodes = len(kyminlist)
       for imode in range(0,nmodes):
           texfhand.write("\\begin{figure}[!h] \n")
           texfhand.write("\\begin{center} \n")
           texfhand.write("\\begin{tabular}{c} \n")
           texfhand.write("\\includegraphics[scale=0.65]{"+reportpath+"phi_mode_%04d.png} \\\\     \n"    % (imode+1))
           texfhand.write("\\includegraphics[scale=0.65]{"+reportpath+"apar_mode_%04d.png}   \n"                              % (imode+1))
           texfhand.write("\\end{tabular} \n")
           texfhand.write("\\begin{tabular}{c c} \n")
           texfhand.write("\\includegraphics[height=17cm,width=10cm]{"+reportpath+"phi_mode_%04d_2d.png}  & \n" % (imode+1))
           texfhand.write("\\includegraphics[height=17cm,width=10cm]{"+reportpath+"apar_mode_%04d_2d.png}   \n" % (imode+1))
           texfhand.write("\\end{tabular} \n")
           texfhand.write("\\end{center} \n")
           texfhand.write("\\end{figure} \n")
           if imode<nmodes-1: texfhand.write("\\clearpage \n")

    texfhand.write("\\end{document} \n")
    texfhand.close()

    print('Processing TEX to PDF ...')
   #os.system('pdflatex -interaction=batchmode -output-directory='+reportpath+" "+reportpath+texfname)
    os.system('pdflatex -output-directory='+reportpath+" "+reportpath+texfname+" > /dev/null")
   #os.system('pdflatex -output-directory='+reportpath+" "+reportpath+texfname)
    print('GENE Report Created.')


from scipy.interpolate import interp1d
if   len(sys.argv) >= 2:
     create_report(sys.argv[1])
     sys.exit()
    #read_mom(sys.argv[1])
    #sys.exit()
    #field = read_field(sys.argv[1])
    #plot_field(field)
    #plot_scans(sys.argv[1])
    #Find rhotor from psiN
     eqdskdata   = read_eqdsk(sys.argv[1])
     psirhotorfn = interp1d(eqdskdata['PSIN'],eqdskdata['rhotor'],kind="linear")
     psiN = 0.96
     print psiN, psirhotorfn(psiN)
     psiN = 0.97
     print psiN, psirhotorfn(psiN)
     psiN = 0.98
     print psiN, psirhotorfn(psiN)
     psiN = 0.99
     print psiN, psirhotorfn(psiN)
     sys.exit()
else:
    while True:
          genefpath=raw_input('Path to GENE files: ')
          if os.path.isdir(genefpath):
             create_report(genefpath)
          else:
             print('Path for GENE files is not found.\nPlease try again ...')
             continue
   
sys.exit()
