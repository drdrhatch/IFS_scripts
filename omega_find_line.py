import optparse as op
import math
import cmath
import sys
import numpy as np
#********************************************************
#***********Define the function**************************
def choose_list(gamma_avg,avg_delta):
  istart=0
  iend=len(gamma_avg)
  #######################################
  #*************Parameters**************#
  #binsize=100*avg_delta
  #binsize = (gamma_avg[iend-1]-gamma_avg[0])/4
  binsize=0.02
  space = 3 #the space between point /density of points is accepted
  length = 5 #minimum length of the line 
  num_of_line=5
  #######################################
  l_bin=int(max(gamma_avg)//binsize) # #of the bins for omega_t
  omega_bin_ID=np.zeros((iend-istart))
  omega_bin=np.zeros(l_bin+1)

  bin_win=np.zeros(l_bin+1) #the bin_number that is has the most count

  #print(gamma_avg)
  l_omega=0
  l_omega_min=l_omega
  for t in range(1,iend-istart+1):
    #print(l_omega)
    bin_num = int(gamma_avg[l_omega]//binsize)  # to put into the right bin
    #print(bin_num)
    omega_bin_ID[l_omega]=bin_num
    #print(omega_bin_ID)
    omega_bin[bin_num]=omega_bin[bin_num]+1
    #print(omega_bin)
    l_omega = l_omega + 1
  #print(max(omega_bin))
  l_omega_max = l_omega
  #print(l_omega_min ,l_omega_max)
  #print(len(omega_bin_ID))
  #print(omega_bin_ID)
  
  line=np.zeros(int(max(omega_bin))) #the line that is plateau
  for t in range(0,l_bin+1):
    if omega_bin[t] == max(omega_bin):
      i = 0
      for s in range(l_omega_min,l_omega_max):
        if int(omega_bin_ID[s]) == t :
          #print(s)
          line[i]= s
          i=i+1
  #print(line)
  #got the line that 
  #########
  #try to exclude the dots that has too much spacing
  lines=np.zeros(num_of_line*len(line)).reshape((num_of_line,len(line)))
  #liens=[line,line,line,line,line]
  l=0 # the # of the line [l, ]
  j=0 # the index of the element of the line
  len_lines=np.zeros(num_of_line) #the length of the line
  temp_l=0 #if line reaches the length then temp_l=1
  m=1 # m is to keep track of the iteration # 
  #print(len(line))
  for n in range(0,num_of_line):
    for i in range(m,len(line)):
      #print(i)
      if line[i]-line[i-1] <= space:
        lines[l,j]=line[i-1]
        j=j+1
        lines[l,j]=line[i]
        if j == length:
          temp_l=1
        temp_j=j
      else:
        j=0
        i=i+1
        break
    m=i
    
    if temp_l==1:
      len_lines[l]=temp_j+1
      l=l+1
      temp_l=0
    if m==len(line)-1:
      break
    
  
  #print("\n number of the lines is: \n", l)
  #print("\n the lines are: \n",lines)
  #print("\n the length of the lines are: \n", len_lines)
  #we got lines, now we need to delet the element of the zero
  
  temp=0  # which line should we choose
  for i in range(0,l):
    if len_lines[i] >= len_lines[temp]:
      temp=i
  #print(lines)
  the_one_temp=lines[temp]
  #int(the_one_temp[0])
  begin=int(the_one_temp[0])
  end=int(the_one_temp[int(len_lines[temp]-1)])
  #print(begin, ",", end)
  the_one=gamma_avg[begin:end]
  begin_end=np.zeros(2)
  begin_end[0]=begin
  begin_end[1]=end+1
  print("*****************************")
  print(begin)
  print(end)
  print((float(begin)/float(iend)))
  print((float(end)/float(iend)))
  #print(gamma_avg)
  #print("\n The chosen segment of the line is: \n",the_one)
  return(begin_end)

#********************************************************
#***********The End of Define the function***************
#gamma_avgest=[1, 10, 10.1, 4, 10, 6, 7, 8, 10.4, 9.9, 7, 11, 10.4, 10.2, 10.2, 10.5, 9.8, 3.4]
#print(choose_list(gamma_avgest))

#line2=np.zeros(int(max(omega_bin)+1)) #the line that is 
#l_omega = (t-1)*(iend-istart)*zmax+z
