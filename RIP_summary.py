import numpy as np
import csv
from max_stat_tool import *
from interp import *
import matplotlib.pyplot as plt

n=5
xaxis=0
yaxis=2
bin_size=20

input0=range(n)
reader=range(n)
z=range(n)
RIP=range(n)
gamma=range(n)

input0[0]="RIP_gauss_n0_19.csv"
input0[1]="RIP_gauss_n0_20.csv"
input0[2]="RIP_gauss_n0_21.csv"
input0[3]="RIP_gauss_n0_22.csv"
input0[4]="RIP_gauss_n0_23.csv"

gamma[0]=0.038647633414777421
gamma[1]=0.023891619743183427
gamma[2]=0.012999999999999999
gamma[3]=0.036999999999999998
gamma[4]=0.0167634273889

RIP_sum = 0
RIP_sum_gamma_weight = 0


for i in range(n):
    
    reader = csv.reader(open(input0[i], "rb"), delimiter=",")
    x = list(reader)
    #print len(x)
    x=x[1:len(x)]
    print((len(x)))
    result = np.array(x).astype("float")
    (nx,ny)=np.shape(result)
    print(nx,ny)
    RIP[i]=result[:,yaxis]
    z[i]=result[:,xaxis]
    rhot0=z[0]
    uni_rhot = np.linspace(min(rhot0),max(rhot0),len(rhot0)*10.)
    RIP[i] = interp(z[i],RIP[i],uni_rhot)

    RIP_sum = RIP_sum+RIP[i]
    RIP_sum_gamma_weight = RIP_sum_gamma_weight+RIP[i]*gamma[i]

RIP_sum = RIP_sum/float(n)
RIP_sum_gamma_weight = RIP_sum_gamma_weight/float(n)

plt.clf()
plt.title('fluctuation Ratio Based On RIP and BES(sum)')
plt.xlabel(r'$\frac{\int n_e \delta B_r dR}{B_0\int n_e dR} \frac{n_e}{\delta n_e}$',fontsize=10)
plt.ylabel('Height (m)',fontsize=13)
plt.plot(smooth(RIP_sum,bin_size)[0], smooth(uni_rhot,bin_size)[0])
#plt.savefig('RIP_gauss_r_smooth.png')
plt.show()

plt.clf()
plt.title('fluctuation Ratio Based On RIP and BES(gamma weighted)')
plt.xlabel(r'$\frac{\int n_e \delta B_r dR}{B_0\int n_e dR} \frac{n_e}{\delta n_e}$',fontsize=10)
plt.ylabel('Height (m)',fontsize=13)
plt.plot(smooth(RIP_sum_gamma_weight,bin_size)[0], smooth(uni_rhot,bin_size)[0])
#plt.savefig('RIP_gauss_r_smooth.png')
plt.show()


