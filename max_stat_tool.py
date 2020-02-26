import optparse as op
import math
import cmath
import sys
import numpy as np

#********************************************************
#***********Define the function**************************
#Function average:
#Use the weighted average to "smooth" out the bump and rescaling
def smooth(avg_list,bin_size):
    l_list=len(avg_list)
    list_avg=np.zeros(l_list-bin_size+1)
    list_dev=np.zeros(l_list-bin_size+1)
    avg=0
    for i in range(0,len(list_avg)):
      list_avg[i]=np.mean(avg_list[i:i+bin_size]) 
      list_dev[i]=np.std(avg_list[i:i+bin_size]) 

    return list_avg, list_dev


def avg(avg_list,bin_size):
#weighted average check chapter 7 in an introduction to error analysisby John R. Taylor
  l_list=len(avg_list)
  #print avg_list
  #print 'l_list', l_list
  avg=0
  dev_sum=0
  #print 'Check point', l_list-bin_size+1
  for i in range(0,l_list-bin_size+1):
    #print i
    avg_temp=np.mean(avg_list[i:i+bin_size])
    std_temp=np.std(avg_list[i:i+bin_size])
    #print avg_temp
    #print std_temp
    avg=avg+avg_temp/(std_temp)**2
    dev_sum=dev_sum+1/(std_temp)**2
  avg=avg/dev_sum #average
  dev_sum=np.sqrt(1/dev_sum)
  output=np.zeros(2)
  output[0]=avg
  output[1]=dev_sum
  return(output)

def avg_dev(avg_list,dev_list):
  bin_size=3 #bin size
  l_list=len(avg_list)
  avg=0
  dev=0
  dev_sum=0
  for i in range(0,l_list-bin_size+1):
    avg_temp=np.mean(avg_list[i:i+bin_size])
    dev_temp=np.mean(dev_list[i:i+bin_size])
    std_temp=np.std(avg_list[i:i+bin_size])
    avg=avg+avg_temp/(std_temp)**2
    dev=dev+dev_temp/(std_temp)**2
    dev_sum=dev_sum+1/(std_temp)**2
  avg=avg/dev_sum #average
  dev=dev/dev_sum #standard dev
  output=np.zeros(2)
  output[0]=avg
  output[1]=dev
  return(output)

def norm(a_list): #Normalized to the 1
  return(a_list/np.max(abs(a_list)))


#def step(list,center,): 

def zoom1D(x,y,x_min,x_max):
#this function zoon in the a 1D plot
  x_zoom=[]
  y_zoom=[]
  for i in range(len(x)):
    if x[i]<=x_max and x[i]>=x_min:
      x_zoom.append(x[i])
      y_zoom.append(y[i])
  return x_zoom, y_zoom

def zoom2D(x,y,z,x_min,x_max,y_min,y_max):
#this function zoon in the a 1D plot
  x_zoom=[]
  y_zoom=[]
  z_zoom=[]
  for i in range(len(x)):
    if x[i]<=x_max and x[i]>=x_min:
      if y[i]<=y_max and y[i]>=y_min:
        x_zoom.append(x[i])
        y_zoom.append(y[i])
        z_zoom.append(z[i])
  return x_zoom, y_zoom, z_zoom

#this function is from unutbu https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
def find_nearest_index(array, value): #this function return the index of the a value(or nearest)
  array = np.asarray(array)
  idx = (np.abs(array - value)).argmin()
  return idx

def find_nearest(array, value): #this function return the index of the a value(or nearest)
  array = np.asarray(array)
  idx = (np.abs(array - value)).argmin()
  return array[idx]


def loop(x,f,x_min,x_max):#this function make the function return [sum f(x+n*x0)] where x0 is the period
  x0=x_max-x_min

  x = np.asarray(x)#same as find_nearest_index
  nx_min = (np.abs(x - x_min)).argmin()#same as find_nearest_index
  nx_max = (np.abs(x - x_max)).argmin()#same as find_nearest_index

  x_loop=x[nx_min:nx_max]
  f_loop=np.zeros(len(x_loop))

  x_loop = np.asarray(x_loop)
  for i in range(len(x)):
    xtemp=(x[i]-x_min)%x0+x_min
    nxtemp = np.argmin(np.abs(x_loop - xtemp)) #same as find_nearest_index
    f_loop[nxtemp]=f_loop[nxtemp]+f[i]

  return x_loop, f_loop
