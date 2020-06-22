#This script is intended to find the top and the mid pedestal of the H mod plasma profile for the pre and post processing of the simulation

#Developed by Max Curie on 01/22/2020

import numpy as np
import matplotlib.pyplot as plt
import re
from efittools import *

def find_pedestal(file_name, path_name, plot):
    if len(path_name)==0:
       path_name = './'
    eqdskdata=read_efit_file(path_name+file_name)

    p  = eqdskdata['pressure'] #pressure read from eqdsk
    x  = eqdskdata['rhotor']   #radial location read from eqdsk
    dp = eqdskdata['pprime']   #First order of pressure

    ddp   = np.gradient(dp,x)  #Second order of pressure

    midped = x[np.argmin(dp)]
    topped = x[np.argmin(ddp)]

    print(topped)
    print(midped)
    
    if plot == 1:
        plt.clf()
        plt.title(r'$Compare$')
        plt.xlabel(r'$r/a$',fontsize=10)
        plt.ylabel(r'$ab$',fontsize=13)
        plt.plot(x,p,label="p")
        #plt.plot(x,-norm(dp),label="dp")
        #plt.plot(x,norm(ddp),label="ddp")
        plt.axvline(x=midped, label="Mid-pedestal", color="red")
        plt.axvline(x=topped, label="Top-pedestal", color="green")
        plt.legend()
        #plt.savefig('compare.png')
        plt.show()

    return midped, topped


#file_name = 'g175823.04108_257x257'
#find_pedestal(file_name, path_name='', plot=1)


