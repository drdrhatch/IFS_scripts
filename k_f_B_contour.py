import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

bin=200

n1_doppler=0.

f_min=100
f_max=600

n_min=4
path='2d_list/36th'


sim_scale=1.25


f_csv=pd.read_csv(path+'/0f_list.csv')  
f=f_csv['f(kHz)']

ky_csv=pd.read_csv(path+'/0ky_list.csv')  
ky=ky_csv['ky']

f1=open(path+"/0B1_matrix_f_ky.csv","r")
lines=f1.readlines() 

def smooth(avg_list,bin_size):
    l_list=len(avg_list)
    list_avg=np.zeros(l_list-bin_size+1)
    list_dev=np.zeros(l_list-bin_size+1)
    avg=0
    for i in range(0,len(list_avg)):
      list_avg[i]=np.mean(avg_list[i:i+bin_size]) 
      list_dev[i]=np.std(avg_list[i:i+bin_size]) 
    return list_avg, list_dev

def smooth_f_x(f,x,dx_size):
    x_min=np.min(x)
    f_avg_list=[]
    f_std_list=[]
    x_avg_list=[]
    x_std_list=[]
    while(x_min<np.max(x)):
        f_temp=[]
        for i in range(len(x)):
            if x_min<=x[i] and x[i]<x_min+dx_size:
                f_temp.append(f[i])
        f_avg_list.append(np.mean(f_temp))
        f_std_list.append(np.std(f_temp))
        x_avg_list.append(x_min+dx_size/2.)
        x_std_list.append(dx_size/2.)
        x_min=x_min+dx_size
        print('x_min='+str(x_min))
    return f_avg_list,f_std_list,x_avg_list,x_std_list


n_list=np.arange(n_min,n_min+len(ky))

print(np.min(f))
print(np.max(f))

uni_freq=np.linspace(np.min(f), np.max(f)+np.max(n_list)*n1_doppler,num=len(f)*10)

frequency_kHZ_uni=np.zeros( (len(ky), len(uni_freq)) )
amplitude_frequency_uni=np.zeros( (len(ky), len(uni_freq)) )
ky_plot=np.zeros( (len(ky), len(uni_freq)) )

i_ky=0

for line in lines:
    line=line[:-1]
    n_temp=n_list[i_ky]
    f_temp=f+n_temp*n1_doppler
    print('**********')
    #print(f_temp)

    #print(line.split(',')) 
    B_f_temp=[]
    for ele in line.split(','):
        B_f_temp.append(float(ele))
    print(len(f_temp))
    print(len(B_f_temp))

    frequency_kHZ_uni[i_ky,:]=uni_freq
    amplitude_frequency_uni[i_ky,:]=np.interp(uni_freq,f_temp,B_f_temp)
    ky_plot[i_ky,:]=[ky[i_ky]]*len(uni_freq)

    i_ky=i_ky+1
    
plt.clf()
plt.ylabel(r'$k_y \rho_s$',fontsize=10)
plt.xlabel(r'$f(kHz)$',fontsize=10)
plt.contourf(frequency_kHZ_uni,ky_plot,amplitude_frequency_uni)#,level=[50,50,50])#,cmap='RdGy')
#plt.contourf(-frequency_kHZ_uni,ky_plot,amplitude_frequency_uni)#,level=[50,50,50])#,cmap='RdGy')
for iky in ky:
    plt.axhline(iky,color='red',alpha=0.5)#alpha control the transparency, alpha=0 transparency, alpha=1 solid
plt.axhline(ky[0],color='red',alpha=0.5,label='n starts from'+str(n_list[0]) )#alpha control the transparency, alpha=0 transparency, alpha=1 solid
plt.xlim(30,600)
plt.ylim(0.03,0.15)
plt.legend()
plt.colorbar()
plt.title(r'$B_r$ contour plot',fontsize=10)
#plt.title(r'$B_r$ contour plot',fontsize=10)
plt.show()

plt.clf()
plt.ylabel(r'$k_y \rho_s$',fontsize=10)
plt.xlabel(r'$f(kHz)$',fontsize=10)
plt.contourf(frequency_kHZ_uni,ky_plot,amplitude_frequency_uni)#,level=[50,50,50])#,cmap='RdGy')
#plt.contourf(-frequency_kHZ_uni,ky_plot,amplitude_frequency_uni)#,level=[50,50,50])#,cmap='RdGy')
plt.xlim(30,600)
plt.ylim(0.03,0.15)
plt.colorbar()
plt.title(r'$B_r$ contour plot',fontsize=10)
#plt.title(r'$B_r$ contour plot',fontsize=10)
plt.show()


plt.clf()
plt.ylabel(r'$k_y \rho_s$',fontsize=10)
plt.xlabel(r'$f(kHz)$',fontsize=10)
plt.contourf(frequency_kHZ_uni,ky_plot,np.log(amplitude_frequency_uni))#,level=[50,50,50])#,cmap='RdGy')
#plt.contourf(-frequency_kHZ_uni,ky_plot,amplitude_frequency_uni)#,level=[50,50,50])#,cmap='RdGy')
for iky in ky:
    plt.axhline(iky,color='red',alpha=0.5)#alpha control the transparency, alpha=0 transparency, alpha=1 solid
plt.axhline(ky[0],color='red',alpha=0.5,label='n starts from'+str(n_list[0]) )#alpha control the transparency, alpha=0 transparency, alpha=1 solid
plt.xlim(0,600)
plt.legend()
plt.colorbar()
plt.title(r'log($B_r$) contour plot',fontsize=10)
#plt.title(r'$B_r$ contour plot',fontsize=10)
plt.show()


Exp_data=pd.read_csv('Experiment/delta_br_vs_f.csv')  

Exp_data['delta_br_vs_f']=Exp_data['delta_br_vs_f']*10000.

row1=uni_freq
row2=np.sum(amplitude_frequency_uni,axis=0)

plt.clf()
plt.plot(row1,row2)
plt.show()

index_min=np.argmin(abs(row1-f_min))
index_max=np.argmin(abs(row1-f_max))

row1=row1[index_min-bin:index_max+bin]
row2=row2[index_min-bin:index_max+bin]

plt.clf()
plt.plot(row1,row2)
plt.show()


f_avg0, f_dev0=smooth(row1,bin)
B1_avg0, B1_dev0=smooth(row2,bin)



f_avg=[]
f_dev=[]

B1_avg=[]
B1_dev=[]

for i in range(len(f_avg0)):
    if i%(bin/2)==0:
        f_avg.append(f_avg0[i])
        f_dev.append(f_dev0[i])
        B1_avg.append(B1_avg0[i])
        B1_dev.append(B1_dev0[i])
print(B1_avg)
sim_max=np.max(B1_avg)*sim_scale
print("sim_max="+str(sim_max))


f_avg_list,f_std_list,x_avg_list,x_std_list=smooth_f_x(Exp_data['delta_br_vs_f'],Exp_data['x'],np.mean(f_dev))

d = {'f(kHz)':x_avg_list,'f_err(kHz)':x_std_list,'B_R(Gauss)':f_avg_list,'B_R_err(Gauss)':f_std_list}
df_exp=pd.DataFrame(d, columns=['f(kHz)','f_err(kHz)','B_R(Gauss)','B_R_err(Gauss)'])
df_exp.to_csv('0_exp_B_r_f_smooth.csv',index=False)

d = {'f(kHz)':f_avg,'f_err(kHz)':f_dev,'B_R(Gauss)':B1_avg,'B_R_err(Gauss)':B1_dev}
df=pd.DataFrame(d, columns=['f(kHz)','f_err(kHz)','B_R(Gauss)','B_R_err(Gauss)'])
df.to_csv('0B_r_f_smooth.csv',index=False)


df_exp['B_R_err(Gauss)']=df_exp['B_R(Gauss)']*0.25

plt.clf()

ax=df_exp.plot(kind='scatter',x='f(kHz)',xerr='f_err(kHz)',y='B_R(Gauss)',yerr='B_R_err(Gauss)',\
	   #ylim=(0.00,sim_max),xlim=(f_min,f_max),
	   color='red',label='Experiment',alpha=1)
df.plot(kind='scatter',x='f(kHz)',xerr='f_err(kHz)',y='B_R(Gauss)',yerr='B_R_err(Gauss)',\
	  #ylim=(0.00,sim_max),xlim=(f_min,f_max),
	   secondary_y=True,color='green',label='GENE simulation',ax=ax)

ax.set_xlabel(r'$frequency(kHz)$',fontsize=15)
ax.set_ylabel(r'$\bar{B}_r(Gauss/\sqrt{Hz})$',fontsize=15)

plt.title(r'$\bar{B}_r$ spectrogram',fontsize=20)
plt.show()
plt.savefig('B1_f.png')


#********RIP***********

