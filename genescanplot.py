import os
import sys
import glob
import numpy as npy
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdfh


def read_parameters(paramfpath):
   #Developed by Ehab Hassan on 2019-02-05
    if "parameters" not in paramfpath:
       if    paramfpath[-1] == "/": paramfpath+="parametrers"
       else:                        paramfpath+="/parameters"
    if os.path.isfile(paramfpath):
         print(paramfpath+' FILE FOUND ...')
    else:
         print(paramfpath+' FILE NOT FOUND. Exit!'); sys.exit()

    ofh = open(paramfpath,'r')
    lines = ofh.readlines()
    ofh.close()

    geneparam = {'filepath':paramfpath}

    nspecs = 0
    for line in lines:
        if   line[0] == '&':
             ckey = line[1:].strip()
             if ckey=="species":
                nspecs += 1
                ckey += "_"+str(nspecs)
             geneparam[ckey] = {}
        elif line[0] == '/' or line[0] == '!':
             continue
        elif line.strip() == '':
             continue
        else:
             items = line.split('=')
             values = items[1].split()
             if "!scanlist" in items[1]:
                cvalues = []
                for item in values[2:]:
                    if item == '!': break
                    cvalues.append(float(item[0:-1]))
                cvalues = npy.array(cvalues)
             else:
               if    values[0].isdigit(): cvalues = float(values[0].strip())
               else:                      cvalues = values[0].strip()
             geneparam[ckey][items[0].strip()] = cvalues

    return geneparam


def read_scanfile(scanfpath):
   #Developed by Ehab Hassan on 2019-02-01
    if "scan.log" not in scanfpath:
       if    scanfpath[-1] == "/": scanfpath+="scan.log"
       else:                       scanfpath+="/scan.log"
    if os.path.isfile(scanfpath):
         print(scanfpath+' FILE FOUND ...')
    else:
         print(scanfpath+' FILE NOT FOUND. Exit!'); sys.exit()

    ofh = open(scanfpath,'r')
    lines = ofh.readlines()
    ofh.close()

    scandata = {'filepath':scanfpath}
    hlist = lines[0].split('|')[1:]
    vlist = []
    for ihead in hlist:
        vlist.append(ihead.split()[0])
    nrecs = npy.size(lines)-1
    vscan = npy.zeros((len(vlist),nrecs))
    gamma = npy.zeros(nrecs)
    omega = npy.zeros(nrecs)
    for irec in range(nrecs):
        line = lines[irec+1].split('|')
        for ivar in range(1,len(vlist)+1):
            vscan[ivar-1,irec] = float(line[ivar])
        gamma[irec] = float(line[-1].split()[0])
        omega[irec] = float(line[-1].split()[1])
    scandata['gamma'] = gamma
    scandata['omega'] = omega
    for ivar in range(len(vlist)):
        scandata[vlist[ivar]] = vscan[ivar,:]
    return scandata

def plot_scandata(scandata):
   #Developed by Ehab Hassan on 2019-02-05
    vlist = list(set(scandata.keys())-{'omega','gamma','filepath'})
    if   len(vlist) == 1:
         if   'kymin' in vlist:
              plabel = scandata['filepath'][0:-9]
         else:
              params = read_parameters(scandata['filepath'][0:-9]+'/parameters')
              plabel = 'kymin = '+params['box']['kymin']
         gammafig = [plt.figure(1)]
         omegafig = [plt.figure(2)]
         axhand01 = gammafig[0].add_subplot(1,1,1)
         axhand02 = omegafig[0].add_subplot(1,1,1)
         axhand01.plot(scandata[vlist[0]],scandata['gamma'],label=plabel)
         axhand02.plot(scandata[vlist[0]],scandata['omega'],label=plabel)
         axhand01.set_title(vlist[0]+' vs $\gamma$')
         axhand02.set_title(vlist[0]+' vs $\omega$')
         axhand01.set_xlabel(vlist[0])
         axhand02.set_xlabel(vlist[0])
         axhand01.set_ylabel('$\gamma$')
         axhand02.set_ylabel('$\omega$')
         axhand01.legend()
         axhand02.legend()
    elif len(vlist) == 2:
         if   'kymin' not in vlist:
              params = read_parameters(scandata['filepath'][0:-9]+'/parameters')
              iscan = 'kymin = '+params['box']['kymin']+': '
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
             axhand01.plot(scandata[maxrepeatvar][sindex:eindex],scandata['gamma'][sindex:eindex],label=plabel)
             axhand02.plot(scandata[maxrepeatvar][sindex:eindex],scandata['omega'][sindex:eindex],label=plabel)
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
             axhand01.plot(scandata[pvarname[0]][sindex::stprng],scandata['gamma'][sindex::stprng],label=plabel)
             axhand02.plot(scandata[pvarname[0]][sindex::stprng],scandata['omega'][sindex::stprng],label=plabel)
         axhand01.set_title(pvarname[0]+' vs $\gamma$')
         axhand02.set_title(pvarname[0]+' vs $\omega$')
         axhand01.set_xlabel(pvarname[0])
         axhand02.set_xlabel(pvarname[0])
         axhand01.set_ylabel('$\gamma$')
         axhand02.set_ylabel('$\omega$')
         axhand01.legend()
         axhand02.legend()

    return gammafig,omegafig


def plot_scans(scanfiles):
   #Developed by Ehab Hassan on 2019-02-07
    if type(scanfiles)==list:
       for iscan in scanfiles:
           scanvals = read_scanfile(iscan)
           gammafig,omegafig = plot_scandata(scanvals)
    else:       
           scanvals = read_scanfile(scanfiles)
           gammafig,omegafig = plot_scandata(scanvals)

    plt.show(gammafig)
    plt.show(omegafig)

    if raw_input('Do you want to save these figures [Yes/No]? ').lower() in ['yes','y']:
       pdfpages = pdfh.PdfPages('scanfigs.pdf')
       for item in range(len(gammafig)):
           gammafig[item].savefig('gamma%02d.png' % (item+1))
           omegafig[item].savefig('omega%02d.png' % (item+1))
           pdfpages.savefig(gammafig[item])
           pdfpages.savefig(omegafig[item])
       pdfpages.close()

    return 1


scanfpath = []
if   len(sys.argv[1:]) >= 1:
     for iarg in sys.argv[1:]:
         scanfpath.append(iarg)
else:
     if raw_input('Do you want to compare several scans in several files? (yes/no)').lower() in ['yes','y']:
        while True:
              infpath = raw_input('Path to scan file: ')
              if len(infpath) > 0: scanfpath.append(infpath)
              else: break
     else:
        scanfpath.append(raw_input('Path to the scan file: '))

plot_scans(scanfpath)


