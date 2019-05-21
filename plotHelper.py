import matplotlib.pyplot as plt
from pylab import rcParams

import numpy as np

def singlePlotABS2D(xgrid, zgrid, \
                 field1, \
                 name1, \
                 bigtitle, filename, \
                 xlabel, ylabel, \
                 output = 'display'):

    rcParams['figure.figsize'] = 8., 8.
    plt.figure()

    #plt.ylabel(ylabel,fontsize=13)
    #plt.xlabel(xlabel,fontsize=13)
    #plt.title('abs ' + name1)
    plt.contourf(xgrid,zgrid,np.abs(field1),100, cmap = "coolwarm")
    plt.tick_params(axis='both',which='major',labelsize=20,\
        length=5,width=2,direction='out')
    plt.locator_params(nbins = 7)
    #plt.colorbar()
    plt.tight_layout()
    plt.suptitle(bigtitle)
    if output == 'display':
        plt.show()
    elif output == 'ps':
        fig=plt.gcf()
        fig.savefig(filename, format = 'ps', bbox_inches = 'tight')

def singlePlot2D(xgrid, zgrid, \
                 field1, \
                 name1, \
                 bigtitle, filename, \
                 xlabel, ylabel, \
                 output = 'display'):

    rcParams['figure.figsize'] = 8., 8.
    plt.figure()

    plt.subplot(3,1,1)
    plt.ylabel(ylabel,fontsize=13)
    plt.xlabel(xlabel,fontsize=13)
    plt.title('abs ' + name1)
    plt.contourf(xgrid,zgrid,np.abs(field1),70, cmap = 'Blues')
    plt.colorbar()
    plt.subplot(3,1,2)
    plt.ylabel(ylabel,fontsize=13)
    plt.xlabel(xlabel,fontsize=13)
    plt.title('real ' + name1)
    plt.contourf(xgrid,zgrid,np.real(field1),70)
    plt.colorbar()
    plt.subplot(3,1,3)
    plt.ylabel(ylabel,fontsize=13)
    plt.xlabel(xlabel,fontsize=13)
    plt.title('imag ' + name1)
    plt.contourf(xgrid,zgrid,np.imag(field1),70)
    plt.colorbar()
    plt.tight_layout()
    plt.suptitle(bigtitle)
    if output == 'display':
        plt.show()
    elif output == 'ps':
        fig=plt.gcf()
        fig.savefig(filename, format = 'ps', bbox_inches = 'tight')

def doublePlot2D(xgrid, zgrid, \
                 field1, field2, \
                 name1, name2, \
                 bigtitle, filename, \
                 xlabel, ylabel, \
                 output = 'display'):

    rcParams['figure.figsize'] = 8., 8.
    plt.figure()

    plt.subplot(3,2,1)
    plt.ylabel(ylabel,fontsize=13)
    plt.xlabel(xlabel,fontsize=13)
    plt.title('abs ' + name1)
    plt.contourf(xgrid,zgrid,np.abs(field1),70)
    plt.colorbar()
    plt.subplot(3,2,3)
    plt.ylabel(ylabel,fontsize=13)
    plt.xlabel(xlabel,fontsize=13)
    plt.title('real ' + name1)
    plt.contourf(xgrid,zgrid,np.real(field1),70)
    plt.colorbar()
    plt.subplot(3,2,5)
    plt.ylabel(ylabel,fontsize=13)
    plt.xlabel(xlabel,fontsize=13)
    plt.title('imag ' + name1)
    plt.contourf(xgrid,zgrid,np.imag(field1),70)
    plt.colorbar()
    plt.subplot(3,2,2)
    plt.ylabel(ylabel,fontsize=13)
    plt.xlabel(xlabel,fontsize=13)
    plt.title('abs ' + name2)
    plt.contourf(xgrid,zgrid,np.abs(field2),70)
    plt.colorbar()
    plt.subplot(3,2,4)
    plt.ylabel(ylabel,fontsize=13)
    plt.xlabel(xlabel,fontsize=13)
    plt.title('real ' + name2)
    plt.contourf(xgrid,zgrid,np.real(field2),70)
    plt.colorbar()
    plt.subplot(3,2,6)
    plt.ylabel(ylabel,fontsize=13)
    plt.xlabel(xlabel,fontsize=13)
    plt.title('imag ' + name2)
    plt.contourf(xgrid,zgrid,np.imag(field2),70)
    plt.colorbar()
    plt.tight_layout()
    plt.suptitle(bigtitle)
    if output == 'display':
        plt.show()
    elif output == 'ps':
        fig=plt.gcf()
        fig.savefig(filename, format = 'ps', bbox_inches = 'tight')

def doublePlot1D(xgrid, \
                 field1, field2, \
                 name1, name2, \
                 bigtitle, filename, \
                 output = 'display'):
    rcParams['figure.figsize'] = 8., 8.
    plt.figure()

    plt.subplot(3,2,1)
    plt.xlabel(r'$\rho_{tor}$',fontsize=13)
    plt.title('abs ' + name1)
    plt.plot(xgrid,np.abs(field1))
    plt.subplot(3,2,3)
    plt.xlabel(r'$\rho_{tor}$',fontsize=13)
    plt.title('real ' + name1)
    plt.plot(xgrid,np.real(field1))
    plt.subplot(3,2,5)
    plt.xlabel(r'$\rho_{tor}$',fontsize=13)
    plt.title('imag ' + name1)
    plt.plot(xgrid,np.imag(field1))
    plt.subplot(3,2,2)
    plt.xlabel(r'$\rho_{tor}$',fontsize=13)
    plt.title('abs ' + name2)
    plt.plot(xgrid,np.abs(field2))
    plt.subplot(3,2,4)
    plt.xlabel(r'$\rho_{tor}$',fontsize=13)
    plt.title('real ' + name2)
    plt.plot(xgrid,np.real(field2))
    plt.subplot(3,2,6)
    plt.xlabel(r'$\rho_{tor}$',fontsize=13)
    plt.title('imag ' + name2)
    plt.plot(xgrid,np.imag(field2))
    plt.tight_layout()
    plt.suptitle(bigtitle)
    if output == 'display':
        plt.show()
    elif output == 'ps':
        fig=plt.gcf()
        fig.savefig(filename, format = 'ps', bbox_inches = 'tight')
