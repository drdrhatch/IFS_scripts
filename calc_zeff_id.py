import numpy as np
import matplotlib.pyplot as plt

de = np.genfromtxt('profiles_e')
di = np.genfromtxt('profiles_i')
dz = np.genfromtxt('profiles_z')



ne = de[:,3]
ni = di[:,3]
nz = dz[:,3]

Te = de[:,2]
Ti = di[:,2]
Tz = dz[:,2]


Z = float(raw_input("Enter Z for impurity:\n"))

zeff = (ni+Z**2*nz)/ne

tau = zeff*Te/Ti

id = ni/ne

plt.plot(de[:,0],zeff)
plt.title('zeff')
plt.show()

plt.plot(de[:,0],id)
plt.title('id')
plt.show()

plt.plot(de[:,0],tau)
plt.title('tau')
plt.show()

rho = float(raw_input("Enter location of interest:\n"))

index=np.argmin(abs(rho-de[:,0]))

print('zeff(x/r='+str(rho)+')='+str(zeff[index]))
print('id(x/r='+str(rho)+')='+str(id[index]))
print('tau(x/r='+str(rho)+')='+str(tau[index]))