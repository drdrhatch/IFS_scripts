import numpy as np
import matplotlib.pyplot as plt
import re

######################modify######################3
######################modify######################3
######################modify######################3
#file_name='profiles_nshift0.02_t3.25.iterdb'
#file_name='profiles_3.25.iterdb'
#file_name = 'profiles_t3.035_nshift0.02.iterdb'
file_name = 'JET92174.iterdb'
#file_name = 'efit_Dial_Nmod1_Zp2_48_new.iterdb'
gene_e = 'gene_profiles_92174_e.dat'
gene_i = 'gene_profiles_92174_i.dat'
#file_name='profiles_3.5.iterdb'
#file_name='iterdb.NSTX_129016A03_460'
#If you want to compare with gene profile files:
gene_plots=True    #set to 0 for no gene plots
plot_impurity = False
if plot_impurity:
    gene_imp = 'profiles_z'
#gene_i='profiles_3.25i.gene'
#gene_e='profiles_3.25e.gene'
#gene_i='profiles_nshift0.02_t3.25i.gene'
#gene_e='profiles_nshift0.02_t3.25e.gene'
######################modify######################3
######################modify######################3
######################modify######################3
if gene_plots:
    gprof_i=np.genfromtxt(gene_i)
    gprof_e=np.genfromtxt(gene_e)
    if plot_impurity:
        gprof_imp=np.genfromtxt(gene_imp)

f=open(file_name,'r')
data_in=f.read()
#print data_in
data_linesplit=data_in.split('\n')

keep_going=1
i=0
while keep_going:
    test=re.search(';-# OF X PTS',data_linesplit[i])
    if test:
        num=data_linesplit[i].split()[0]
        num=float(num)
        num=int(num)
        print ("number of points:",num)
        keep_going=(1==2)
    if i == len(data_linesplit):
        keep_going=(1==2)
    i=i+1


def plot_next(data_linesplit,lnum,num):
    sec_num_lines = int(num/6)
    if num % 6 != 0:
        sec_num_lines += 1
    keep_going=1
    while keep_going:
        test=re.search('-DEPENDENT VARIABLE LABEL',data_linesplit[lnum])
        #test2=re.search('INDEPENDENT',data_linesplit[lnum])
        if test :
            quantity=data_linesplit[lnum].split()[0]
            units=data_linesplit[lnum].split()[1]
            print ("Plotting :",quantity)
        test=re.search('DATA FOLLOW',data_linesplit[lnum])
        if test:
            keep_going=(1==2)
        lnum=lnum+1
    
    rhot=np.empty(0)
    lnum0 = lnum
    for j in range(lnum0,lnum0+sec_num_lines):
        for k in range(6):
            str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
            if(str_temp):
                temp=np.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                rhot=np.append(rhot,temp)
        lnum=lnum+1
    #print rhot
    print ("time=",data_linesplit[lnum])
    #print "var=",quantity
    lnum=lnum+1
    
    arr=np.empty(0)
    lnum0 = lnum
    for j in range(lnum0,lnum0+sec_num_lines):
        for k in range(6):
            str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
            if(str_temp):
                temp=np.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                arr=np.append(arr,temp)
        #arr_str=data_linesplit[j].split()
        #arr_flt=np.array(arr_str,dtype='float')
        #arr=np.append(arr,arr_flt)
        lnum=lnum+1
    #print 'rhot',rhot
    #print 'arr',arr
    #print len(rhot),len(arr)
    if quantity=='VROT':
        vout=np.empty((len(arr),2),dtype='float')
        vout[:,0]=rhot
        vout[:,1]=arr
        f=open(quantity+'.dat','w')
        f.write('#'+file_name+'\n'+'#rhot '+quantity+'\n')
        np.savetxt(f,vout)
        f.close()
    
    plt.plot(rhot,arr,label=quantity+'('+units+')')
    plt.xlabel(r'$\rho_{tor}$',size=18)
    plt.legend(loc='lower left')
    plt.show()
    
    if gene_plots and quantity=='TI':
        plt.plot(rhot,arr/1000.0,label=quantity+'('+units+')')
        plt.plot(gprof_i[:,0],gprof_i[:,2],'r.',label='$T_i$ gene')
        plt.legend(loc='lower left')
        plt.show()
    if gene_plots and quantity=='NM1':
        plt.plot(rhot,arr/1.0e19,label=quantity+'('+units+')/1.0e19')
        plt.plot(gprof_i[:,0],gprof_i[:,3],'r.',label='$n_i$ gene')
        plt.legend(loc='lower left')
        plt.show()
    if gene_plots and quantity=='NM2' and plot_impurity:
        plt.plot(rhot,arr/1.0e19,label=quantity+'('+units+')/1.0e19')
        plt.plot(gprof_imp[:,0],gprof_imp[:,3],'r.',label='$n_i$ gene')
        plt.legend(loc='lower left')
        plt.show()
    if gene_plots and quantity=='TE':
        plt.plot(rhot,arr/1000.0,label=quantity+'('+units+')')
        plt.plot(gprof_e[:,0],gprof_e[:,2],'r.',label='$T_e$ gene')
        plt.legend(loc='lower left')
        plt.show()
    if gene_plots and quantity=='NE':
        plt.plot(rhot,arr/1.0e19,label=quantity+'('+units+')/1.0e19')
        plt.plot(gprof_e[:,0],gprof_e[:,3],'r.',label='$n_e$ gene')
        plt.legend(loc='lower left')
        plt.show()



    lnum_out=lnum
    try_again=1
    if len(data_linesplit)-lnum < 10:
        try_again=(1==2)
    return lnum_out,try_again

lnum=0
try_again=1
while try_again:
    lnum,try_again=plot_next(data_linesplit,lnum,num)
    #print lnum,len(data_linesplit)



