from read_iterdb_x import *
from read_EFIT import *
from finite_differences_x import *
from interp import *
import matplotlib.pyplot as plt
import numpy as np
import sys

ITERDBdict1 = read_iterdb_x(sys.argv[1])
ITERDBdict2 = read_iterdb_x(sys.argv[2])
impurityCharge = 6.

if 1 == 1:
    plt.plot(ITERDBdict1['rhot_vrot'], ITERDBdict1['vrot'], label = '1')
    plt.plot(ITERDBdict2['rhot_vrot'], ITERDBdict2['vrot'], '.', label = '2')
    plt.axis((0.9, 1., -1E5, 1E5))
    plt.title('vrot (rad/s)')
    plt.xlabel('rhot')
    plt.legend()
    plt.show()

if 1 == 1:
    plt.plot(ITERDBdict1['rhot_ti'], ITERDBdict1['ti'], label = '1')
    plt.plot(ITERDBdict2['rhot_ti'], ITERDBdict2['ti'], '.', label = '2')
    plt.title('ti (eV)')
    plt.xlabel('rhot')
    plt.legend()
    plt.show()
if 1 == 1:
    plt.plot(ITERDBdict1['rhot_te'], ITERDBdict1['te'], label = '1')
    plt.plot(ITERDBdict2['rhot_te'], ITERDBdict2['te'], '.', label = '2')
    plt.title('te (eV)')
    plt.xlabel('rhot')
    plt.legend()
    plt.show()
if 1 == 1:
    plt.plot(ITERDBdict1['rhot_ni'], ITERDBdict1['ni']*1E-20*10, label = '1')
    plt.plot(ITERDBdict2['rhot_ni'], ITERDBdict2['ni']*1E-20*10, '.', label = '2')
    plt.title('ni (m^-3)')
    plt.xlabel('rhot')
    plt.legend()
    plt.show()
if 1 == 1:
    plt.plot(ITERDBdict1['rhot_ne'], ITERDBdict1['ne']*1E-20*10, label = '1')
    plt.plot(ITERDBdict2['rhot_ne'], ITERDBdict2['ne']*1E-20*10, '.', label = '2')
    plt.title('ne (m^-3)')
    plt.xlabel('rhot')
    plt.legend()
    plt.show()
if 1 == 1:
    plt.plot(ITERDBdict1['rhot_nz'], ITERDBdict1['nz'], label = '1')
    plt.plot(ITERDBdict2['rhot_nz'], ITERDBdict2['nz'], '.', label = '2')
    plt.title('nz (m^-3)')
    plt.xlabel('rhot')
    plt.legend()
    plt.show()
if 1 == 1:
    plt.plot(ITERDBdict1['rhot_ne'], (ITERDBdict1['ni'] + impurityCharge**2 * ITERDBdict1['nz']) / ITERDBdict1['ne'], label = '1')
    plt.plot(ITERDBdict2['rhot_ne'], (ITERDBdict2['ni'] + impurityCharge**2 * ITERDBdict2['nz']) / ITERDBdict2['ne'], '.', label = '2')
    plt.title('Zeff')
    plt.xlabel('rhot')
    plt.legend()
    plt.show()
