import numpy as np
import matplotlib.pyplot as plt
from sys import path
from sys import exit
from get_nrg import *
import os
from finite_differences import *


def calc_gr2(suffix,nspec=2,ncols=10):
    """dens,upar,tpar,tperp,Ges,Gem,Qes,Qem,Pes,Pem"""

    if nspec == 2:
        time,nrgi,nrge=get_nrg0(suffix,nspec=nspec,ncols=ncols)
    elif nspec == 3:
        time,nrgi,nrg2,nrge=get_nrg0(suffix,nspec=nspec,ncols=ncols)
        #sys.exit("Must have n_spec=2")
    else:
        time,nrgi = get_nrg0(suffix,nspec=nspec,ncols=ncols)

    #print  "1. nrgi[-1,0]/nrgi[-2,0]",nrgi[-1,0]/nrgi[-2,0]
    if nrgi[-1,0]/nrgi[-2,0] < 1.0e-10:
        nrgi=np.delete(nrgi,-1,0)
        if nspec > 1:
            nrge=np.delete(nrge,-1,0)
        time=np.delete(time,-1,0)

    if time[-1] > 80:
        start_time=time[-1]-2.0
    else:
        start_time=0.95*time[-1]
    start_index=np.argmin(abs(time-start_time))
    ntime=len(time)-start_index

    dlogdt=np.zeros((ntime,5))
    for i in range(ntime-1):
        i0=i+start_index
        dlogdt[i,0]=0.5*(nrgi[i0+1,0]-nrgi[i0,0])/(time[i0+1]-time[i0])/(0.5*(nrgi[i0+1,0]+nrgi[i0,0]))
        #print dlogdt[i,0]
        dlogdt[i,1]=0.5*(nrgi[i0+1,2]-nrgi[i0,2])/(time[i0+1]-time[i0])/(0.5*(nrgi[i0+1,2]+nrgi[i0,2]))
        dlogdt[i,2]=0.5*(nrgi[i0+1,3]-nrgi[i0,3])/(time[i0+1]-time[i0])/(0.5*(nrgi[i0+1,3]+nrgi[i0,3]))
        dlogdt[i,3]=0.5*(nrgi[i0+1,6]-nrgi[i0,6])/(time[i0+1]-time[i0])/(0.5*(nrgi[i0+1,6]+nrgi[i0,6]))
        dlogdt[i,4]=0.5*(nrgi[i0+1,7]-nrgi[i0,7])/(time[i0+1]-time[i0])/(0.5*(nrgi[i0+1,7]+nrgi[i0,7]))

    avg_gr=np.zeros(10)

    for i in range(5):
        avg_gr[i]=np.sum(dlogdt[:-1,i])/len(dlogdt[:-1,i])

    for i in range(ntime-1):
        if nspec > 1:
            i0=i+start_index
            dlogdt[i,0]=0.5*(nrge[i0+1,0]-nrge[i0,0])/(time[i0+1]-time[i0])/(0.5*(nrge[i0+1,0]+nrge[i0,0]))
            dlogdt[i,1]=0.5*(nrge[i0+1,2]-nrge[i0,2])/(time[i0+1]-time[i0])/(0.5*(nrge[i0+1,2]+nrge[i0,2]))
            dlogdt[i,2]=0.5*(nrge[i0+1,3]-nrge[i0,3])/(time[i0+1]-time[i0])/(0.5*(nrge[i0+1,3]+nrge[i0,3]))
            dlogdt[i,3]=0.5*(nrge[i0+1,6]-nrge[i0,6])/(time[i0+1]-time[i0])/(0.5*(nrge[i0+1,6]+nrge[i0,6]))
            dlogdt[i,4]=0.5*(nrge[i0+1,7]-nrge[i0,7])/(time[i0+1]-time[i0])/(0.5*(nrge[i0+1,7]+nrge[i0,7]))

    for i in range(5,10):
        avg_gr[i]=np.sum(dlogdt[:-1,i-5])/len(dlogdt[:-1,i-5])
 
    
    momname=list()
    momname.append('ni')
    momname.append('Tpar_i  ')
    momname.append('Tperp_i ')
    momname.append('Qes_i   ')
    momname.append('Qem_i   ')
    momname.append('ne      ')
    momname.append('Tpar_e  ')
    momname.append('Tperp_e ')
    momname.append('Qes_e   ')
    momname.append('Qem_e  ')

    #print avg_gr
    #print "Select growth rate to keep:"
    #print "Average Growth Rates:"
    #print "0:ni      ",avg_gr[0]
    #print "1:Tpar_i  ",avg_gr[1]
    #print "2:Tperp_i ",avg_gr[2]
    #print "3:Qes_i   ",avg_gr[3]
    #print "4:Qem_i   ",avg_gr[4]
    #print "5:ne      ",avg_gr[5]
    #print "6:Tpar_e  ",avg_gr[6]
    #print "7:Tperp_e ",avg_gr[7]
    #print "8:Qes_e   ",avg_gr[8]
    #print "9:Qem_e  ",avg_gr[9]
    #print "-1:none"

    fit =  np.e**(2.0*avg_gr[0]*time[start_index:])*(nrgi[-1,0]/np.e**(2.0*avg_gr[0]*time[-1]))
    err = abs(np.sum(nrgi[start_index:,0]-fit[:])/np.sum(nrgi[start_index:,0]))

    print( "Calculated growth rate:",avg_gr[0])
    print( "Error:",err)
    if err > 1.0e-2:
        plt.semilogy(time,nrgi[:,0],'-x')
        plt.semilogy(time[start_index-500:],np.e**(2.0*avg_gr[0]*time[start_index-500:])*(nrgi[-1,0]/np.e**(2.0*avg_gr[0]*time[-1])),'--',color='green')
        plt.show()
        test = input("Accept calculation? (y=yes)")
        if test=='y':
            return avg_gr[0]
        else:
            return calc_gr(suffix,nspec=nspec,ncols=ncols)
    else:
        return avg_gr[0]


def calc_gr(suffix,nspec=2,ncols=10):
    """dens,upar,tpar,tperp,Ges,Gem,Qes,Qem,Pes,Pem"""

    if nspec == 1:
        time,nrgi=get_nrg0(suffix,nspec=nspec,ncols=ncols)
    if nspec == 2:
        time,nrgi,nrge=get_nrg0(suffix,nspec=nspec,ncols=ncols)
    elif nspec == 3:
        time,nrgi,nrge,nrg3=get_nrg0(suffix,nspec=nspec,ncols=ncols)
    plt.semilogy(time[1:],nrgi[1:,0],label='n')
    plt.semilogy(time[1:],nrgi[1:,2],label='tpar')
    plt.semilogy(time[1:],nrgi[1:,3],label='tperp')
    plt.semilogy(time[1:],nrgi[1:,6],label='Qes')
    plt.semilogy(time[1:],nrgi[1:,7],label='Qem')
    plt.title('ions')
    plt.legend(loc='upper left')
    plt.xlabel('t')
    plt.show()
    #plt.savefig('nrgi_plot.ps')
    #plt.close()

    if nspec==2:
        plt.semilogy(time[1:],nrge[1:,0],label='n')
        plt.semilogy(time[1:],nrge[1:,2],label='tpar')
        plt.semilogy(time[1:],nrge[1:,3],label='tperp')
        plt.semilogy(time[1:],nrge[1:,6],label='Qes')
        plt.semilogy(time[1:],nrge[1:,7],label='Qem')
        plt.title('electrons')
        plt.legend(loc='upper left')
        plt.xlabel('t')
        plt.show()
        #plt.savefig('nrge_plot.ps')
        #plt.close()

    start_time=input("Enter start time:")
    start_time=int(float(start_time))
    start_index=np.argmin(abs(time-start_time))
    ntime=len(time)-start_index
    #print "start_time",start_time
    #print "start_index",start_index
    #print "time at start_index",time[start_index]
    #temp = 0.5*fd_d1_o4(np.log(nrgi[start_index:,0]),time[start_index:])
    #print "nrgi[start_index:,0]",nrgi[start_index:,0]
    #print "time[start_index:]",time[start_index:]
    #print "temp",temp
    
    dlogdt = np.zeros((ntime,5))
    dlogdt[:,0] = 0.5*fd_d1_o4(np.log(nrgi[start_index:,0]),time[start_index:])
    dlogdt[:,1] = 0.5*fd_d1_o4(np.log(nrgi[start_index:,2]),time[start_index:])
    dlogdt[:,2] = 0.5*fd_d1_o4(np.log(nrgi[start_index:,3]),time[start_index:])
    dlogdt[:,3] = 0.5*fd_d1_o4(np.log(nrgi[start_index:,6]),time[start_index:])
    dlogdt[:,4] = 0.5*fd_d1_o4(np.log(nrgi[start_index:,7]),time[start_index:])

    #for i in range(ntime-1):
    #    i0=i+start_index
    #    dlogdt[i,0]=0.5*(nrgi[i0+1,0]-nrgi[i0,0])/(time[i0+1]-time[i0])/(0.5*(nrgi[i0+1,0]+nrgi[i0,0]))
    #    dlogdt[i,1]=0.5*(nrgi[i0+1,2]-nrgi[i0,2])/(time[i0+1]-time[i0])/(0.5*(nrgi[i0+1,2]+nrgi[i0,2]))
    #    dlogdt[i,2]=0.5*(nrgi[i0+1,3]-nrgi[i0,3])/(time[i0+1]-time[i0])/(0.5*(nrgi[i0+1,3]+nrgi[i0,3]))
    #    dlogdt[i,3]=0.5*(nrgi[i0+1,6]-nrgi[i0,6])/(time[i0+1]-time[i0])/(0.5*(nrgi[i0+1,6]+nrgi[i0,6]))
    #    dlogdt[i,4]=0.5*(nrgi[i0+1,7]-nrgi[i0,7])/(time[i0+1]-time[i0])/(0.5*(nrgi[i0+1,7]+nrgi[i0,7]))

    avg_gr=np.zeros(10)
    for i in range(5):
        avg_gr[i]=np.sum(dlogdt[2:-3,i])/len(dlogdt[2:-3,i])

    norm = nrgi[start_index+2,1]*np.e**(-2.0*avg_gr[0]*time[start_index+2])
    plt.semilogy(time[start_index+2:-3],norm*np.e**(2.0*avg_gr[0]*time[start_index+2:-3]),'-x',label='n(fit)')
    plt.semilogy(time[start_index:-1],nrgi[start_index:-1,1],'-x',label='n')
    plt.legend(loc='upper left')
    plt.xlabel('t')
    plt.show()


    if nspec==2:
        dlogdt=np.zeros((ntime,5))
        dlogdt[:,0] = 0.5*fd_d1_o4(np.log(nrge[start_index:,0]),time[start_index:])
        dlogdt[:,1] = 0.5*fd_d1_o4(np.log(nrge[start_index:,2]),time[start_index:])
        dlogdt[:,2] = 0.5*fd_d1_o4(np.log(nrge[start_index:,3]),time[start_index:])
        dlogdt[:,3] = 0.5*fd_d1_o4(np.log(nrge[start_index:,6]),time[start_index:])
        dlogdt[:,4] = 0.5*fd_d1_o4(np.log(nrge[start_index:,7]),time[start_index:])
        #for i in range(ntime-1):
        #    i0=i+start_index
        #    dlogdt[i,0]=0.5*(nrge[i0+1,0]-nrge[i0,0])/(time[i0+1]-time[i0])/(0.5*(nrge[i0+1,0]+nrge[i0,0]))
        #    dlogdt[i,1]=0.5*(nrge[i0+1,2]-nrge[i0,2])/(time[i0+1]-time[i0])/(0.5*(nrge[i0+1,2]+nrge[i0,2]))
        #    dlogdt[i,2]=0.5*(nrge[i0+1,3]-nrge[i0,3])/(time[i0+1]-time[i0])/(0.5*(nrge[i0+1,3]+nrge[i0,3]))
        #    dlogdt[i,3]=0.5*(nrge[i0+1,6]-nrge[i0,6])/(time[i0+1]-time[i0])/(0.5*(nrge[i0+1,6]+nrge[i0,6]))
        #    dlogdt[i,4]=0.5*(nrge[i0+1,7]-nrge[i0,7])/(time[i0+1]-time[i0])/(0.5*(nrge[i0+1,7]+nrge[i0,7]))

        #plt.plot(time[start_index:-1],dlogdt[:-1,0],label='n')
        #plt.plot(time[start_index:-1],dlogdt[:-1,1],label='tpar')
        #plt.plot(time[start_index:-1],dlogdt[:-1,2],label='tperp')
        #plt.plot(time[start_index:-1],dlogdt[:-1,3],label='Qes')
        #plt.plot(time[start_index:-1],dlogdt[:-1,4],label='Qem')
        #plt.title('electrons 0.5*logarithmic derivative')
        #plt.legend(loc='upper left')
        #plt.xlabel('t')
        #plt.show()
        for i in range(5,10):
            avg_gr[i]=np.sum(dlogdt[2:-3,i-5])/len(dlogdt[2:-3,i-5])

 
    momname=list()
    momname.append('ni')
    momname.append('Tpar_i  ')
    momname.append('Tperp_i ')
    momname.append('Qes_i   ')
    momname.append('Qem_i   ')
    momname.append('ne      ')
    momname.append('Tpar_e  ')
    momname.append('Tperp_e ')
    momname.append('Qes_e   ')
    momname.append('Qem_e  ')

    #print avg_gr
    print( "Select growth rate to keep:")
    print( "Average Growth Rates:")
    print( "0:ni      ",avg_gr[0])
    print( "1:Tpar_i  ",avg_gr[1])
    print( "2:Tperp_i ",avg_gr[2])
    print( "3:Qes_i   ",avg_gr[3])
    print( "4:Qem_i   ",avg_gr[4])
    print( "5:ne      ",avg_gr[5])
    print( "6:Tpar_e  ",avg_gr[6])
    print( "7:Tperp_e ",avg_gr[7])
    print( "8:Qes_e   ",avg_gr[8])
    print( "9:Qem_e  ",avg_gr[9])
    print( "-1:none")


    selection=input()
    sel=int(float(selection)) 
    if sel<10 and sel > -1:
        print( "Returning growth rate for"+momname[sel]+':',avg_gr[sel])
        return avg_gr[sel]
    elif sel==-1:
        print( "Returning -1: not sufficently converged.")
        return -1
    else:
        print( "Invalid selection.  Returning -1: not sufficiently converged.")
        return -1


