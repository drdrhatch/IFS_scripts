import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
from finite_differences import *
from max_stat_tool import gaussian
from max_stat_tool import gaussian_fit
from max_stat_tool import Poly_fit
from MTMDispersion_tools import Parameter_reader

#*********************Start of user block********************************
profile_type="ITERDB"           # "ITERDB" "pfile" 
geomfile_type="gfile"          # "gfile"  "GENE_tracor"

path='/global/u1/m/maxcurie/max/Cases/jet78697/'
profile_name =path+'jet78697.51005_hager_Z6.0Zeff2.35__negom_alpha1.2_TiTe.iterdb' 		#name of the profile file
                                            #DIIID175823.iterdb
                                            #p000000
geomfile_name = 'jet78697.51005_hager.eqdsk'

manual_ped=0
mid_ped0=0.958
q_scale=1.0
#*********************End of user block********************************



uni_rhot,nu,eta,shat,beta,ky,q,mtmFreq,omegaDoppler,\
	omega_n,omni_GENE,omti_GENE,omne_GENE,omte_GENE,omega_n_GENE,\
	xstar,Lref,rhoref=Parameter_reader(profile_type,profile_name,geomfile_type,\
	geomfile_name,q_scale,manual_ped,mid_ped0,plot=False,output_csv=False,suffix='dat')



#********start of q polyfit********************
#https://numpy.org/doc/stable/reference/generated/numpy.polyfit.html
print("start of q polyfit")
order=5

coeff=Poly_fit(uni_rhot, q, order, show=True)


f = open('0q_coeff_parameter','w')
f.write('q_coeffs =')
for i in range(len(coeff)):
    if i==len(coeff)-1:
        f.write(str(coeff[i]))
    else:
        f.write(str(coeff[i])+', ')
f.close()
#********end of of q polyfit********************


xgrid = np.linspace(0.,1.0,500)
#*********start of omne fit***************
print("start of omne fit")
amplitude,mean,stddev=gaussian_fit(uni_rhot,omne_GENE)
print(amplitude,mean,stddev)
omn=gaussian(xgrid, amplitude, mean, stddev)  #analytical omn
#*********end of omne fit***************

#*********start of omte fit***************
print("start of omte fit")
amplitude,mean,stddev=gaussian_fit(uni_rhot,omte_GENE)
print(amplitude,mean,stddev)
omt=gaussian(xgrid, amplitude, mean, stddev)
#*********end of omte fit***************

#*********start of omni fit***************
print("start of omni fit")
amplitude,mean,stddev=gaussian_fit(uni_rhot,omni_GENE)
print(amplitude,mean,stddev)
omni=gaussian(xgrid, amplitude, mean, stddev)  #analytical omn
#*********end of omni fit***************

#*********start of omti fit***************
print("start of omti fit")
amplitude,mean,stddev=gaussian_fit(uni_rhot,omti_GENE)
print(amplitude,mean,stddev)
omti=gaussian(xgrid, amplitude, mean, stddev)
#*********end of omti fit***************



T0 = 2.14
KT = 0.00
KTi = 6.92
oT = 0.25
R0 = 2.8


Tint = np.empty(len(xgrid))
Tiint = np.empty(len(xgrid))
for i in range(1,len(xgrid)):
    Tint[i] = scipy.integrate.simps(omt[0:i]/R0,xgrid[0:i])
    Tiint[i] = scipy.integrate.simps(omti[0:i]/R0,xgrid[0:i])
T = np.e**(-Tint)
Ti = np.e**(-Tiint)
ix0 = np.argmin(abs(xgrid-0.5))
T = T0/T[ix0]*T
Ti = T0/Ti[ix0]*Ti

Ktest = -R0*fd_d1_o4(np.log(T),xgrid)

plt.plot(xgrid,T,label='Te')
plt.plot(xgrid,Ti,label='Ti')
plt.xlabel('x')
plt.ylabel('T')
plt.legend()
plt.show()

plt.plot(xgrid,Ktest,'x',label='Ktest')
plt.plot(xgrid,omt,label='omt')
plt.plot(xgrid,omti,label='omti')
plt.xlabel('x')
plt.ylabel('omt')
plt.legend()
plt.show()

n0 = 6.483e-2
Kn = 2.22
on = 0.25

nint = np.empty(len(xgrid))
for i in range(1,len(xgrid)):
    nint[i] = scipy.integrate.simps(omn[0:i]/R0,xgrid[0:i])
n = np.e**(-nint)
n = n0/n[ix0]*n

Ktest = -R0*fd_d1_o4(np.log(n),xgrid)

plt.plot(xgrid,n)
plt.xlabel('x')
plt.ylabel('n')
plt.show()

plt.plot(xgrid,Ktest,'x',label='Ktest')
plt.plot(xgrid,omn,label='omn')
plt.xlabel('x')
plt.ylabel('omn')
plt.legend()
plt.show()

f = open('profiles_e_fit','w')
f.write('# 1.xgrid 2.xgrid 3.T(kev) 4.n(10^19m^-3)\n')
np.savetxt(f,np.column_stack((xgrid,xgrid,T,n)))
f.close()

f = open('profiles_i_fit','w')
f.write('# 1.xgrid 2.xgrid 3.Ti(kev) 4.n(10^19m^-3)\n')
np.savetxt(f,np.column_stack((xgrid,xgrid,Ti,n)))
f.close()
