import numpy as np
import matplotlib.pyplot as plt
import re

def read_iterdb(filename):
    f=open(filename,'r')
    data_in=f.read()
    data_linesplit=data_in.split('\n')

    keep_going=1
    i=0
    while keep_going:
        test=re.search(';-# OF X PTS',data_linesplit[i])
        if test:
            num=data_linesplit[i].split()[0]
            num=float(num)
            num=int(num)
            print "number of points:",num
            keep_going=(1==2)
        if i == len(data_linesplit):
            keep_going=(1==2)
        i=i+1

    lnum=0
    try_again=1
    prof_out = {}
    rhot_out = {}
    units_out = {}
    while try_again:
        lnum,try_again,quantity,units,rhot,arr=get_next(data_linesplit,lnum,num)
        prof_out[quantity]=arr
        units_out[quantity]=units
        rhot_out[quantity]=rhot
    return rhot_out,prof_out,units_out   
     

def get_next(data_linesplit,lnum,num):
    sec_num_lines = num/6
    if num % 6 != 0:
        sec_num_lines += 1
    keep_going=1
    while keep_going:
        test=re.search('-DEPENDENT VARIABLE LABEL',data_linesplit[lnum])
        #test2=re.search('INDEPENDENT',data_linesplit[lnum])
        if test :
            quantity=data_linesplit[lnum].split()[0]
            units=data_linesplit[lnum].split()[1]
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
    #if quantity=='VROT':
    #    vout=np.empty((len(arr),2),dtype='float')
    #    vout[:,0]=rhot
    #    vout[:,1]=arr
    #    f=open(quantity+'.dat','w')
    #    f.write('#'+file_name+'\n'+'#rhot '+quantity+'\n')
    #    np.savetxt(f,vout)
    #    f.close()
    
    #plt.plot(rhot,arr,label=quantity+'('+units+')')
    #plt.xlabel(r'$\rho_{tor}$',size=18)
    #plt.legend(loc='lower left')
    #plt.show()
    
    #if gene_plots and quantity=='TI':
    #    plt.plot(rhot,arr/1000.0,label=quantity+'('+units+')')
    #    plt.plot(gprof_i[:,0],gprof_i[:,2],'r.',label='$T_i$ gene')
    #    plt.legend(loc='lower left')
    #    plt.show()

    lnum_out=lnum
    try_again=1
    if len(data_linesplit)-lnum < 10:
        try_again=False
    return lnum_out, try_again,quantity,units,rhot,arr


