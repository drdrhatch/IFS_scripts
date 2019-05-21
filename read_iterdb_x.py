import matplotlib.pyplot as plt
import numpy as np
import re

def read_iterdb_x(iterdb_filename):

    ITERDBdict = {}

    iterdb_file=open(iterdb_filename,'r')

    data_in=iterdb_file.read()
    data_linesplit=data_in.split('\n')

    keep_going=1
    i=0
    while keep_going:
        test=re.search(';-# OF X PTS',data_linesplit[i])
        if test:
            num=data_linesplit[i].split()[0]
            num=float(num)
            num=int(num)
            print("number of points in iterdb file:",num)
            keep_going=(1==2)
        if i == len(data_linesplit):
            keep_going=(1==2)
        i=i+1

    lnum=0
    while len(data_linesplit)-lnum > 10:
        sec_num_lines = num//6
        if num % 6 != 0:
            sec_num_lines += 1
        keep_going=1
        while keep_going:
            test=re.search('-DEPENDENT VARIABLE LABEL',data_linesplit[lnum])
            #test2=re.search('INDEPENDENT',data_linesplit[lnum])
            if test :
                quantity=data_linesplit[lnum].split()[0]
                units=data_linesplit[lnum].split()[1]
            test2=re.search('DATA FOLLOW',data_linesplit[lnum])
            if test2:
                keep_going=(1==2)
            lnum=lnum+1

        if quantity=='TE':
           print("Reading :",quantity)
           rhot_te=np.empty(0)
           lnum0 = lnum
           for j in range(lnum0,lnum0+sec_num_lines):
               for k in range(6):
                   str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
                   if(re.search('e',str_temp)):
                       temp=np.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                       rhot_te=np.append(rhot_te,temp)
               lnum=lnum+1
           ITERDBdict['rhot_te'] = rhot_te

           lnum=lnum+1

           te=np.empty(0)
           lnum0 = lnum
           for j in range(lnum0,lnum0+sec_num_lines):
               for k in range(6):
                   str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
                   if(re.search('e',str_temp)):
                      temp=np.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                      te=np.append(te,temp)
               lnum=lnum+1
           ITERDBdict['te'] = te

        if quantity=='TI':
           print("Reading :",quantity)
           rhot_ti=np.empty(0)
           lnum0 = lnum
           for j in range(lnum0,lnum0+sec_num_lines):
               for k in range(6):
                   str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
                   if(re.search('e',str_temp)):
                       temp=np.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                       rhot_ti=np.append(rhot_ti,temp)
               lnum=lnum+1
           ITERDBdict['rhot_ti'] = rhot_ti

           lnum=lnum+1

           ti=np.empty(0)
           lnum0 = lnum
           for j in range(lnum0,lnum0+sec_num_lines):
               for k in range(6):
                   str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
                   if(re.search('e',str_temp)):
                      temp=np.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                      ti=np.append(ti,temp)
               lnum=lnum+1
           ITERDBdict['ti'] = ti
#        else:
#           lnum = lnum + 2*sec_num_lines + 1

        if quantity=='NE':
           print("Reading :",quantity)
           rhot_ne=np.empty(0)
           lnum0 = lnum
           for j in range(lnum0,lnum0+sec_num_lines):
               for k in range(6):
                   str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
                   if(re.search('e',str_temp)):
                       temp=np.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                       rhot_ne=np.append(rhot_ne,temp)
               lnum=lnum+1
           ITERDBdict['rhot_ne'] = rhot_ne

           lnum=lnum+1

           ne=np.empty(0)
           lnum0 = lnum
           for j in range(lnum0,lnum0+sec_num_lines):
               for k in range(6):
                   str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
                   if(re.search('e',str_temp)):
                      temp=np.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                      ne=np.append(ne,temp)
               lnum=lnum+1
           ITERDBdict['ne'] = ne

        if quantity=='NM1':
           print("Reading :",quantity)
           rhot_ni=np.empty(0)
           lnum0 = lnum
           for j in range(lnum0,lnum0+sec_num_lines):
               for k in range(6):
                   str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
                   if(re.search('e',str_temp)):
                       temp=np.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                       rhot_ni=np.append(rhot_ni,temp)
               lnum=lnum+1
           ITERDBdict['rhot_ni'] = rhot_ni

           lnum=lnum+1

           ni=np.empty(0)
           lnum0 = lnum
           for j in range(lnum0,lnum0+sec_num_lines):
               for k in range(6):
                   str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
                   if(re.search('e',str_temp)):
                      temp=np.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                      ni=np.append(ni,temp)
               lnum=lnum+1
           ITERDBdict['ni'] = ni

        if quantity=='NM2':
           print("Reading :",quantity)
           rhot_nz=np.empty(0)
           lnum0 = lnum
           for j in range(lnum0,lnum0+sec_num_lines):
               for k in range(6):
                   str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
                   if(re.search('e',str_temp)):
                       temp=np.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                       rhot_nz=np.append(rhot_nz,temp)
               lnum=lnum+1
           ITERDBdict['rhot_nz'] = rhot_nz

           lnum=lnum+1

           nz=np.empty(0)
           lnum0 = lnum
           for j in range(lnum0,lnum0+sec_num_lines):
               for k in range(6):
                   str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
                   if(re.search('e',str_temp)):
                      temp=np.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                      nz=np.append(nz,temp)
               lnum=lnum+1
           ITERDBdict['nz'] = nz

        if quantity=='VROT':
           print("Reading :",quantity)
           rhot_vrot=np.empty(0)
           lnum0 = lnum
           for j in range(lnum0,lnum0+sec_num_lines):
               for k in range(6):
                   str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
                   if(re.search('e',str_temp)):
                       temp=np.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                       rhot_vrot=np.append(rhot_vrot,temp)
               lnum=lnum+1
           ITERDBdict['rhot_vrot'] = rhot_vrot

           lnum=lnum+1

           vrot=np.empty(0)
           lnum0 = lnum
           for j in range(lnum0,lnum0+sec_num_lines):
               for k in range(6):
                   str_temp=data_linesplit[j][1+k*13:1+(k+1)*13]
                   if(re.search('e',str_temp)):
                      temp=np.array(data_linesplit[j][1+k*13:1+(k+1)*13],dtype='float')
                      vrot=np.append(vrot,temp)
               lnum=lnum+1
           ITERDBdict['vrot'] = vrot
    return ITERDBdict
