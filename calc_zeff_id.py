import numpy as np
import matplotlib.pyplot as plt

de = np.genfromtxt('profiles_e')
di = np.genfromtxt('profiles_i')
dz = np.genfromtxt('profiles_z')

ne = de[:,3]
ni = di[:,3]
nz = dz[:,3]

Z = float(raw_input("Enter Z for impurity:\n"))

zeff = (ni+Z**2*nz)/ne
id = ni/ne

plt.plot(de[:,0],zeff)
plt.title('zeff')
plt.show()

plt.plot(de[:,0],id)
plt.title('id')
