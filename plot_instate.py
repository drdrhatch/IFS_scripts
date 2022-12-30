import numpy as np
import matplotlib.pyplot as plt
from fastran_data import *
import scipy.integrate

mod_jtot = False
jtot_factor = 0.8

e = 1.6022e-19

data = read_instate('instate')
plt.plot(data['rhot'],data['ne'])
plt.xlabel('rhot')
plt.title('ne(10^19 m^-3)')
plt.show()

plt.plot(data['rhot'],data['te'])
plt.title('te(kev)')
plt.xlabel('rhot')
plt.show()

plt.plot(data['rhot'],data['ti'])
plt.title('ti(kev)')
plt.xlabel('rhot')
plt.show()

plt.plot(data['rhot'],data['zeff'])
plt.title('ZEFF')
plt.xlabel('rhot')
plt.show()

ptot = data['ne']*1e19*data['te']*1000*e + data['ni']*1e19*data['ti']*1000*e

plt.plot(data['rhot'],data['p_eq'],label='from instate')
plt.plot(data['rhot'],ptot,label='reconstructed')
plt.legend()
plt.title('P_EQ')
plt.show()

Psource_tot = scipy.integrate.simpson(data['se_nb']*4*np.pi**2*data['rhot']*data['aminor']*data['rmajor'],data['aminor']*data['rhot'])
Psource_tot *= 1e19

plt.plot(data['rhot'],data['se_nb'])
plt.title('Total particle source (particles / s): '+str(Psource_tot))
plt.ylabel('Particles/m^3/s')
plt.xlabel('rhot')
plt.show()

plt.scatter(data['rbdry'],data['zbdry'])
plt.scatter(data['rlim'],data['zlim'])
plt.show()
#print('len(rlim)',len(data['rlim']))
#print('len(zlim)',len(data['zlim']))
#print('len(rbdry)',len(data['rbdry']))
#print('len(zbdry)',len(data['zbdry']))
#print('nbdry',data['nbdry'])
#print('nlim',data['nlim'])
#print('rbdry',data['rbdry'])
#print('zbdry',data['zbdry'])
#print('rlim',data['rlim'])
#print('zlim',data['zlim'])

Ip = scipy.integrate.simpson(data['j_tot']*2*np.pi*data['rhot']*data['aminor'],data['aminor']*data['rhot'])

plt.plot(data['rhot'],data['j_tot'],color = 'black',label='j_tot')
plt.plot(data['rhot'],data['j_oh'],color = 'red',label='j_oh')
plt.title('Integrated Ip:'+str(Ip)[0:5])
plt.ylabel('J(MA/m^2)')
plt.xlabel('rhot')
plt.legend()
plt.show()

plt.plot(data['rhot'],data['q0'],label='q')
plt.ylabel('q')
plt.xlabel('rhot')
plt.show()

if mod_jtot:
    plt.plot(data['rhot'],data['j_tot'],label='OG')
    plt.plot(data['rhot'],jtot_factor*data['j_tot'],label='Mod')
    output_four_col(jtot_factor*data['j_tot'],'jtot.dat')
    plt.xlabel('rhot')
    plt.ylabel('jtot')
    plt.legend()
    plt.show()





