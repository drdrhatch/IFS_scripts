import numpy as np
from subprocess import call
import os
from interp import *
from finite_differences import *
import sys
import matplotlib.pyplot as plt

kx_center_scan = False
setup_global = False

x0_scan = True
setup_lilo = False

######Modify
include_impurity = False
casedir = 'ion_scale_174082'
case = 'ion_scale_174082'
qfile = 'Binfo_g174082.3000'
efit_file_name = 'g174082.3000'
gene_file_name = 'profiles_e'
iterdb_file = 'DIIID174082.iterdb'

#casedir = 'AUG_27963'
#case = 'AUG_27963_S1'
#qfile = 'Binfo_27963.03250'
#efit_file_name = 'g027963.03250'
#gene_file_name = 'profiles_t3.25_nshift0.02_vpole.gene'
#iterdb_file = 'profiles_t3.25_nshift0.02_vpol.iterdb'

x0_values=[0.94,0.95,0.96,0.97]
ky_scan_string = '0.02, 0.04, 0.06, 0.08, 0.1, 0.125, 0.15, 0.2'
edge_opt = 2.0
#ky_scan_string = '1.0, 5.0, 10.0, 20.0, 40.0, 60.0, 100.0, 140.0'
template_prob_dir_loc  = 'prob_template_local'
template_prob_dir_glob = 'prob_template_global'
submit_runs = False
#x0_glob = '0.975'
x0_glob = '0.9275'
lx_a = '0.135'
nx0_glob = 320
batch_script_job_prefix = '#SBATCH -J '
submit_command = 'sbatch'
ExB_glob = 0.0   #'scan' for both 0.0 and -1111.0
### For kx_center scan
num_kxcenter = 5
gene_dirname = 'knl-refactored/'
######Modify

#####Setup
#####Setup
#####Setup
x0_scan_string = str(x0_values[0])
for i in range(1,len(x0_values)):
    x0_scan_string += ', '+str(x0_values[i])
print "x0_scan_string:",x0_scan_string
basedir='/global/cscratch1/sd/halfmoon/'
diagdir = basedir+case
#homedir = '/home1/01658/drhatch/'
homedir = '/global/homes/h/halfmoon/'
genedir = homedir + gene_dirname
probdirloc = genedir + 'prob_loc_' + case
probdirglob = genedir + 'prob_glob_' + case
probdirlilo = genedir + 'prob_lilo_' + case
probdirkxc = []
for i in range(len(x0_values)):
    probdirkxc.append(genedir+'prob_kxc_'+case+'_'+str(x0_values[i]))
pmvdir = homedir
eqs_dir = homedir+gene_dirname+case+'/'
print "Checking existence of efit file."
print eqs_dir
if os.path.isfile(eqs_dir+efit_file_name):
    print "Efit file exists:",efit_file_name
else:
    sys.exit("Efit file does not exist.  Select different efit file.") 
if os.path.isfile(eqs_dir+iterdb_file):
    print "Iterdb file exists:",case+'.iterdb'
else:
    sys.exit("Iterdb file does not exist.  Select different iterdb file.") 

#####Setup
#####Setup
#####Setup

def calc_shat(q,rhot):
    rhot0 = np.arange(10000.0)/9999.0
    q0 = interp(rhot,q,rhot0)

    qprime = fd_d1_o4(q0,rhot0)
    shat = rhot0/q0*qprime
    return rhot0,q0,shat 

def calc_shat_wpsi(qin,psiin,rhot,rhop,rhot_range=[0.87,1.027]):
    drhot = rhot_range[1]-rhot_range[0]
    rhot0 = np.arange(10000.0)/9999.0*drhot+rhot_range[0]
    ind0 = np.argmin(abs(rhot[:]-rhot_range[0]))
    q0 = full_interp(qin,psiin,rhop[ind0:]**2,rhot[ind0:],rhot0,verify_interp = True)

    qprime = fd_d1_o4(q0,rhot0)
    shat = rhot0/q0*qprime
    plt.plot(rhot0,q0)
    plt.title('q')
    plt.show()
    plt.plot(rhot0,shat)
    plt.title('shat')
    plt.show()
    return rhot0,q0,shat 

#####Execute
#####Execute
#####Execute

if not os.path.exists(diagdir):
    print diagdir," does not exist."  
    os.makedirs(diagdir)
    print "Now it does."

if x0_scan:
    if include_impurity:
       print "Warning: include_impurity not ready for x0_scan!"
       stop
    if os.path.exists(probdirloc):
        call(['rm','-r',probdirloc])
    call(['cp','-r',template_prob_dir_loc,probdirloc])
    call(['cp','-r',homedir+'/scripts/setup_linear_runs.py',probdirloc])

if setup_global:
    if include_impurity:
       print "Warning: include_impurity not ready for global!"
       stop
    if os.path.exists(probdirglob):
        call(['rm','-r',probdirglob])
    call(['cp','-r',template_prob_dir_glob,probdirglob])
    call(['cp','-r',homedir+'/scripts/setup_linear_runs.py',probdirglob])
if setup_lilo:
    if include_impurity:
       print "Warning: include_impurity not ready for global!"
       stop
    if os.path.exists(probdirlilo):
        call(['rm','-r',probdirlilo])
    call(['cp','-r',template_prob_dir_glob,probdirlilo])
    call(['cp','-r',homedir+'/scripts/setup_linear_runs.py',probdirlilo])

if kx_center_scan:
    for dir in probdirkxc:
        if os.path.exists(dir):
            call(['rm','-r',dir])
        call(['cp','-r',template_prob_dir_loc,dir])
        call(['cp','-r',homedir+'/scripts/setup_linear_runs.py',dir])

#### LOCAL ####
if x0_scan:
    os.chdir(probdirloc)
    call(['ln','-s',eqs_dir+efit_file_name])
    call(['ln','-s',eqs_dir+iterdb_file])
    call(['ln','-s',eqs_dir+'gene_profiles_e'])
    call(['ln','-s',eqs_dir+'gene_profiles_i'])

    #change submit.cmd file
    f=open('submit.cmd','r')
    cmd_data=f.read()
    f.close()
    cmd_data_split = cmd_data.split('\n')
    cmd_data_split[1] = batch_script_job_prefix+case+'         # Job Name'
    cmd_file_out='\n'.join(cmd_data_split)
    f=open('submit.cmd','w')
    f.write(cmd_file_out)
    f.close()

    #change diagdir
    f=open('parameters','r')
    parfile=f.read()
    f.close()
    parfile_split = parfile.split('\n')
    for i in range(len(parfile_split)):
        if 'diagdir' in parfile_split[i]:
            parfile_split[i] = 'diagdir = \''+diagdir+'\''
        if 'edge_opt' in parfile_split[i]:
            parfile_split[i] = 'edge_opt = '+str(edge_opt)
        if 'x0' in parfile_split[i] and 'nx0' not in parfile_split[i]:
            parfile_split[i] = 'x0 = 0.95  !scanlist: '+x0_scan_string
        if 'kymin' in parfile_split[i]: 
            parfile_split[i] = 'kymin = 0.05  !scanlist: '+ky_scan_string
        if 'geomfile' in parfile_split[i]: 
            parfile_split[i] = 'geomfile = ' + '\''+efit_file_name + '\''
        if 'iterdb_file' in parfile_split[i]: 
            parfile_split[i] = 'iterdb_file = ' + '\''+ iterdb_file + '\''
    parfile_out='\n'.join(parfile_split)
    f=open('parameters','w')
    f.write(parfile_out)
    f.close()

    if submit_runs:
       call([submit_command,'submit.cmd',])


##### GLOBAL #####
if setup_global:
    os.chdir(probdirglob)
    call(['ln','-s',eqs_dir+efit_file_name])
    call(['ln','-s',eqs_dir+iterdb_file])
    call(['ln','-s','gene_profiles_e'])
    call(['ln','-s','gene_profiles_i'])


    #change submit.cmd file
    f=open('submit.cmd','r')
    cmd_data=f.read()
    f.close()
    cmd_data_split = cmd_data.split('\n')
    cmd_data_split[1] = batch_script_job_prefix+case+'_glob     # Job Name'
    cmd_file_out='\n'.join(cmd_data_split)
    f=open('submit.cmd','w')
    f.write(cmd_file_out)
    f.close()

    #change diagdir
    f=open('parameters','r')
    parfile=f.read()
    f.close()
    parfile_split = parfile.split('\n')
    for i in range(len(parfile_split)):
        if 'diagdir' in parfile_split[i]:
            parfile_split[i] = 'diagdir = \''+diagdir+'\''
        if 'edge_opt' in parfile_split[i]:
            parfile_split[i] = 'edge_opt = '+str(edge_opt)
        if 'x0' in parfile_split[i] and 'nx0' not in parfile_split[i]:
            parfile_split[i] = 'x0 = '+x0_glob
        if 'nx0' in parfile_split[i]:
            parfile_split[i] = 'nx0 = '+str(nx0_glob)
        if 'lx_a' in parfile_split[i]:
            parfile_split[i] = 'lx_a = '+lx_a
        if 'kymin' in parfile_split[i]: 
            parfile_split[i] = 'kymin = 0.05  !scanlist: '+ky_scan_string
        if 'geomfile' in parfile_split[i]: 
            parfile_split[i] = 'geomfile = ' + '\''+efit_file_name + '\''
        if 'iterdb_file' in parfile_split[i]: 
            parfile_split[i] = 'iterdb_file = ' + '\''+ iterdb_file + '\''
        if 'ExBrate' in parfile_split[i] and ExB_glob: 
            if ExB_glob == 'scan':
                parfile_split[i] = 'ExBrate = -1111.0 !scanlist: 0.0, -1111.0' 
            else:
                parfile_split[i] = 'ExBrate = -1111.0' 
    parfile_out='\n'.join(parfile_split)
    f=open('parameters','w')
    f.write(parfile_out)
    f.close()

    if submit_runs:
       call([submit_command,'submit.cmd'])

##### LILO #####
if setup_lilo:
    os.chdir(probdirlilo)
    call(['ln','-s',eqs_dir+efit_file_name])
    call(['ln','-s',eqs_dir+case+'.iterdb'])
    call(['ln','-s','gene_profiles_e'])
    call(['ln','-s','gene_profiles_i'])

    #change submit.cmd file
    f=open('submit.cmd','r')
    cmd_data=f.read()
    f.close()
    cmd_data_split = cmd_data.split('\n')
    cmd_data_split[1] = batch_script_job_prefix+case+'_lilo     # Job Name'
    cmd_file_out='\n'.join(cmd_data_split)
    f=open('submit.cmd','w')
    f.write(cmd_file_out)
    f.close()

    #change diagdir
    f=open('parameters','r')
    parfile=f.read()
    f.close()
    parfile_split = parfile.split('\n')
    for i in range(len(parfile_split)):
        if 'diagdir' in parfile_split[i]:
            parfile_split[i] = 'diagdir = \''+diagdir+'\''
        if 'edge_opt' in parfile_split[i]:
            parfile_split[i] = 'edge_opt = '+str(edge_opt)
        if 'x0' in parfile_split[i] and 'nx0' not in parfile_split[i]:
            parfile_split[i] = 'x0 = '+x0_glob
        if 'lx_a' in parfile_split[i]:
            parfile_split[i] = 'lx_a = '+lx_a
        if 'kymin' in parfile_split[i]: 
            parfile_split[i] = 'kymin = 0.05  !scanlist: '+ky_scan_string
        if 'geomfile' in parfile_split[i]: 
            parfile_split[i] = 'geomfile = ' + '\''+efit_file_name + '\''
        if 'iterdb_file' in parfile_split[i]: 
            parfile_split[i] = 'iterdb_file = ' + '\''+ iterdb_file + '\''
        if 'ExBrate' in parfile_split[i] and ExB_glob: 
            if ExB_glob == 'scan':
                parfile_split[i] = 'ExBrate = -1111.0 !scanlist: 0.0, -1111.0' 
            else:
                parfile_split[i] = 'ExBrate = -1111.0' 
        if 'lilo' in parfile_split[i]: 
            parfile_split[i] = 'lilo = T' 
        if 'rad_bc_type' in parfile_split[i]: 
            parfile_split[i] = 'rad_bc_type = -1' 
        if 'drive_buffer' in parfile_split[i]: 
            parfile_split[i] = 'drive_buffer = F' 
        if 'mag_prof' in parfile_split[i]: 
            parfile_split[i] = 'mag_prof = F' 
    parfile_out='\n'.join(parfile_split)
    f=open('parameters','w')
    f.write(parfile_out)
    f.close()

    if submit_runs:
       call([submit_command,'submit.cmd'])

##### kx_center_scan #####
if kx_center_scan:
    call(['ln','-s',homedir+'scripts/scan_info_efit.py',diagdir])
    call(['ln','-s',homedir+'scripts/plot_scan_info_efit.py',diagdir])
    call(['ln','-s',homedir+'scripts/plot_contour_x0_ky_efit.py',diagdir])
    for j in range(len(x0_values)):
        os.chdir(probdirkxc[j])
        call(['ln','-s',eqs_dir+efit_file_name])
        call(['ln','-s',eqs_dir+iterdb_file])

        #Calculating shat
        qfilein = np.genfromtxt(eqs_dir+qfile)
        q = qfilein[:,4]
        psi = qfilein[:,1]
        genefile = np.genfromtxt(eqs_dir+gene_file_name)
        ind8 = np.argmin(abs(genefile[:,0]-0.8))
        rhot0,q0,shat = calc_shat_wpsi(q,psi,genefile[ind8:,0],genefile[ind8:,1])
        xind = np.argmin(abs(rhot0-x0_values[j]))
        shat_out = shat[xind]
        print "Assuming shat = ",shat_out
        if j==0:
            print "Saving shat.dat"
            np.savetxt('shat.dat',np.column_stack((rhot0,q0,shat)))

        #change submit.cmd file
        f=open('submit.cmd','r')
        cmd_data=f.read()
        f.close()
        cmd_data_split = cmd_data.split('\n')
        cmd_data_split[1] = batch_script_job_prefix+case+'kxc_x0'+str(x0_values[j])+'         # Job Name'
        cmd_file_out='\n'.join(cmd_data_split)
        f=open('submit.cmd','w')
        f.write(cmd_file_out)
        f.close()

        kx_center_scan_string = '  !scanlist: 0.0 '
        for i in range(num_kxcenter-1):
            kx_center_scan_string += ', '+str(0.8*(i+1)/float(num_kxcenter)*2*np.pi)+'*'+str(shat_out)+'*kymin(1)'

        print "kx_center_scan_string",kx_center_scan_string
        #change diagdir
        f=open('parameters','r')
        parfile=f.read()
        f.close()
        parfile_split = parfile.split('\n')
        for i in range(len(parfile_split)):
            if 'diagdir' in parfile_split[i]:
                parfile_split[i] = 'diagdir = \''+diagdir+'\''
            if 'edge_opt' in parfile_split[i]:
                parfile_split[i] = 'edge_opt = '+str(edge_opt)
            if 'x0' in parfile_split[i] and 'nx0' not in parfile_split[i]:
                parfile_split[i] = 'x0 = '+str(x0_values[j])
            if 'kymin' in parfile_split[i]: 
                parfile_split[i] = 'kymin = 0.05  !scanlist: '+ky_scan_string
            if 'geomfile' in parfile_split[i]: 
                parfile_split[i] = 'geomfile = ' + '\''+efit_file_name + '\''
            if 'iterdb_file' in parfile_split[i]: 
                parfile_split[i] = 'iterdb_file = ' + '\''+ iterdb_file + '\''
            if 'kx_center' in parfile_split[i]: 
                parfile_split[i] = 'kx_center = 0.0  ' + kx_center_scan_string
            if include_impurity and 'n_spec' in parfile_split[i]:
                parfile_split[i] = 'n_spec = 3 ' 
            if 'species' in parfile_split[i]:
                spec_ind.append(i)
            if include_impurity and 'n_procs_s' in parfile_split[i]:
                parfile_split[i] = 'n_procs_s = 1 ' 
            if include_impurity and 'n_procs_z' in parfile_split[i]:
                parfile_split[i] = 'n_procs_z = 8 ' 

        if include_impurity:
            imp_charge = raw_input('Enter impurity charge:\n')
            imp_mass = raw_input('Enter impurity mass:\n')
            imp_label = raw_input('Enter impurity label:\n')
        #print 'spec_ind',spec_ind
            parfile_split.insert(spec_ind[1],'\n')
            parfile_split.insert(spec_ind[1],'/')
            parfile_split.insert(spec_ind[1],'prof_type = -2')
            parfile_split.insert(spec_ind[1],'temp = -1')
            parfile_split.insert(spec_ind[1],'dens = -1')
            parfile_split.insert(spec_ind[1],'omn = -1')
            parfile_split.insert(spec_ind[1],'omt = -1')
            parfile_split.insert(spec_ind[1],'charge = '+str(imp_charge))
            parfile_split.insert(spec_ind[1],'mass = '+str(imp_mass))
            parfile_split.insert(spec_ind[1],'name = '+'\''+imp_label+'\'')
            parfile_split.insert(spec_ind[1],'&species')

        parfile_out='\n'.join(parfile_split)
        f=open('parameters','w')
        f.write(parfile_out)
        f.close()

        if submit_runs:
           call([submit_command,'submit.cmd',])

