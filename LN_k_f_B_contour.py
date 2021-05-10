import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

bin=10

n1_doppler=-0

contour_plot_ky_max=0.25

f_min=0
f_max=600

n_min=0
kymin=0.03
#n_step=3.77
ky_1=0.00775536103953903
n_step=kymin/ky_1
print(n_step)
path='2d_list/36_midplan_hann_no_square'

Frequency_flip=True  #change to True if one wants to look at electron direction

if Frequency_flip==True:
    label0='electron direction f(kHz)'
else:
    label0='electron direction f(kHz)'

#path='2d_list/36_midplan_hann_V2'

#ky rho_s for n=1
#ky_1=0.008


plot_scale=1.5 #How much blank space one wants in y axis
sim_scale=1.    #the y_plot_sim=y_sim*sim_scale
norm_to_max=True    #Ture if one wants to noralized sim and exp to their respected maximum

plot_linear=False  #plot linear contour
plot_log=False    #plot log contour
plot_subplot=True      #Plot subplot for the different ky
total_row=3			#Total rows for the subplot

#omstar=(1.2*8.77304279894492+abs(n1_doppler))
eta=2.16596001839912
omn=2.77105293432639

omegastar_nor=omn*(1.+eta)
omstar=(omegastar_nor+abs(n1_doppler))
#omegastar_nor_120=omn*(1.+1.2*eta)
#omstar120=(omegastar_nor_120+abs(n1_doppler))

f_csv=pd.read_csv(path+'/0f_list.csv')  
f=f_csv['f(kHz)']


LL_csv=pd.read_csv('./LL_MTM.csv')  
LL_f_plasma=LL_csv['omega(kHz)']

LL_gamma=LL_csv['gamma(cs/a)']
LL_ky=LL_csv['kymin']
LL_f=LL_f_plasma+LL_ky/ky_1*n1_doppler

ky_csv=pd.read_csv(path+'/0ky_list.csv')  
ky=ky_csv['ky']

f1=open(path+"/0B1_matrix_f_ky.csv","r")
lines=f1.readlines() 

if Frequency_flip==True:
    Frequency_flip_constant=-1.
else:
    Frequency_flip_constant=1.

def smooth(avg_list,bin_size):
    if bin_size==1:
        list_avg=avg_list
        dev=np.mean(abs(list_avg[1:]-list_avg[:-1]))*0.5
        list_dev=np.array([dev]*len(list_avg))
    else:
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


n_list=np.arange(n_min,n_min+len(ky)*n_step,n_step)

print('n_list'+str(n_list))

#print(np.min(f))
#print(np.max(f))

uni_freq=np.linspace(np.min(f), np.max(f)+np.max(n_list)*n1_doppler,num=len(f)*10)

frequency_kHZ_uni=np.zeros( (len(ky), len(uni_freq)) )
amplitude_frequency_uni=np.zeros( (len(ky), len(uni_freq)) )
ky_plot=np.zeros( (len(ky), len(uni_freq)) )

i_ky=0

for line in lines:
    line=line[:-1]
    n_temp=n_list[i_ky]
    print('n_temp='+str(n_temp))
    f_temp=f+n_temp*n1_doppler
    #print('**********')
    #print(f_temp)

    #print(line.split(',')) 
    B_f_temp=[]
    for ele in line.split(','):
        B_f_temp.append(float(ele))
    #print(len(f_temp))
    #print(len(B_f_temp))

    frequency_kHZ_uni[i_ky,:]=uni_freq
    amplitude_frequency_uni[i_ky,:]=np.interp(uni_freq,f_temp,B_f_temp)
    ky_plot[i_ky,:]=[ky[i_ky]]*len(uni_freq)

    i_ky=i_ky+1

uni_freq=Frequency_flip_constant*uni_freq
frequency_kHZ_uni=Frequency_flip_constant*frequency_kHZ_uni
ky_list=np.arange(np.min(ky_plot),np.max(ky_plot),0.0000001)

print('ky_list'+str(ky_list))

if plot_linear==True:
    plt.clf()
    plt.ylabel(r'$k_y \rho_s$',fontsize=10)
    plt.xlabel(r'$f(kHz)$',fontsize=10)
    plt.plot(omstar*ky_list/ky_1,ky_list,label=r'line of $\omega_{*e}$')
    plt.contourf(frequency_kHZ_uni,ky_plot,amplitude_frequency_uni,levels=10000,extend='both')#[100,50,50])#,cmap='RdGy')
    #plt.contourf(-frequency_kHZ_uni,ky_plot,amplitude_frequency_uni)#,level=[50,50,50])#,cmap='RdGy')

    for i_n in range(int(np.max(n_list))):
        if i_n%5==0:
            plt.axhline(i_n*ky_1,color='red',alpha=0.5)#alpha control the transparency, alpha=0 transparency, alpha=1 solid
        else:
            plt.axhline(i_n*ky_1,color='red',alpha=0.1)#alpha control the transparency, alpha=0 transparency, alpha=1 solid
    plt.axhline(ky[0],color='red',alpha=0.5,label='n starts from '+str(int(n_list[0])) )#alpha control the transparency, alpha=0 transparency, alpha=1 solid
    plt.xlim(0,600)
    plt.ylim(0,contour_plot_ky_max)
    plt.legend()
    plt.colorbar()
    plt.title(r'$B_r$ contour plot',fontsize=10)
    plt.savefig('contour.png')
    plt.show()

if plot_log==True: 
    plt.clf()
    plt.ylabel(r'$k_y \rho_s$',fontsize=10)
    plt.xlabel(r'$f(kHz)$',fontsize=10)
    #plt.plot(omstar120*ky_list/ky_1,ky_list,label=r'line of $\omega_{*e}$ with 1.2 $\omega_{*Te}$')
    plt.plot(omstar*ky_list/ky_1,ky_list,label=r'line of $\omega_{*e}$')
    for i in range(len(LL_f)):
        if i == 1:
            plt.plot(-LL_f[i],LL_ky[i],'o',color='orange',markersize=10.*LL_gamma[i],alpha=0.5,label=r'Local linear simulations')
        else: 
            plt.plot(-LL_f[i],LL_ky[i],'o',color='orange',markersize=10.*LL_gamma[i],alpha=0.5)
    plt.contourf(frequency_kHZ_uni,ky_plot,np.log10(amplitude_frequency_uni),levels=1000,extend='both')#,cmap='RdGy')
    #plt.contourf(-frequency_kHZ_uni,ky_plot,amplitude_frequency_uni)#,level=[50,50,50])#,cmap='RdGy')
    for i_n in range(int(np.max(n_list))):
        if i_n%5==0:
            plt.axhline(i_n*ky_1,color='red',alpha=0.5)#alpha control the transparency, alpha=0 transparency, alpha=1 solid
        else:
            plt.axhline(i_n*ky_1,color='red',alpha=0.1)#alpha control the transparency, alpha=0 transparency, alpha=1 solid
    plt.axhline(ky[0],color='red',alpha=0.5,label='n starts from '+str(int(n_list[0])) )#alpha control the transparency, alpha=0 transparency, alpha=1 solid

    plt.xlim(0,1000)
    plt.ylim(0,contour_plot_ky_max*2)
    plt.colorbar()
    plt.legend(loc='upper left')
    plt.title(r'log($B_r$) contour plot',fontsize=10)
    #plt.title(r'$B_r$ contour plot',fontsize=10)
    plt.savefig('contour_log.png')
    plt.show()


if plot_subplot==True:
    for i in range(total_row):
        for j in range(int(np.ceil(len(n_list)/float(total_row)))):
            plt.clf()
            x_TEMP=frequency_kHZ_uni[3*i+j,:]
            y_TEMP=amplitude_frequency_uni[3*i+j,:]
            print(np.shape(frequency_kHZ_uni))
            print(np.shape(amplitude_frequency_uni))

            ky_TEMP=ky[3*i+j]
            print(ky_TEMP)
            print(round(ky_TEMP,4))
            plt.plot(x_TEMP,y_TEMP)
            plt.grid()
            plt.title('ky='+str(round(ky_TEMP,4)))
            plt.xlim([0, abs(f_min)])
            plt.savefig('ky='+str(round(ky_TEMP,4))+'.png')

    fig, ax=plt.subplots(nrows=3,ncols=int(np.ceil(len(n_list)/3.)),sharex=True,sharey=True) 
    		#nrows is the total rows
    		#ncols is the total columns
    		#sharex true means the xaxies will be shared
    i_TEMP=0
    for i in range(total_row):
        for j in range(int(np.ceil(len(n_list)/float(total_row)))):     
            x_TEMP=frequency_kHZ_uni[i_TEMP,:]
            y_TEMP=amplitude_frequency_uni[i_TEMP,:]

            ky_TEMP=ky[i_TEMP]
            print(ky_TEMP)
            print(round(ky_TEMP,4))
            ax[i,j].plot(x_TEMP,y_TEMP)
            ax[i,j].set_title('ky='+str(round(ky_TEMP,4)))
            ax[i,j].set_xlim([abs(f_max), abs(f_min)])
            ax[i,j].grid()
            i_TEMP=i_TEMP+1
            #i_TEMP=3*i+j
            if i==total_row-1:
                ax[i,j].set_xlabel('frequency(kHz)')
            #ax[i,j].set_ylabel(r'$B_r(Gauss/\sqrt{kHz})$')

    #ax[0,0].plot(x,np.sin(x),label='sin(x)')
    #ax[0,0].set_xlabel('x')
    #ax[0,0].set_ylabel('sin(x)',fontsize=fontsize0)
    #ax1.set_title()		#for the set the title name
    #ax[1,1].set_xlim([abs(f_max), abs(f_min)])

    plt.tight_layout()
    plt.show()


Exp_data=pd.read_csv('Experiment/delta_br_vs_f.csv')  

Exp_data['delta_br_vs_f']=Exp_data['delta_br_vs_f']*10000.

row1=uni_freq
row2=(np.sum(amplitude_frequency_uni**2.,axis=0)*2.)**0.5

plt.clf()
plt.plot(row1,row2)
plt.title('plot_all_range '+label0)
plt.xlabel('frequency(kHz)')
plt.show()

def sort_x_f(x_unsort,f_unsort): 
    arr_unsort=[x_unsort,f_unsort]
    f_x_unsort=tuple(map(tuple, np.transpose(arr_unsort)))
    f_x_sort=sorted(f_x_unsort, key=lambda f_x_unsort: f_x_unsort[0])
    f_x_sort=np.array(f_x_sort)
    f_x_sort=np.transpose(f_x_sort)
    x_sort=f_x_sort[0,:]
    f_sort=f_x_sort[1,:]
    return x_sort,f_sort

row1,row2=sort_x_f(row1,row2)

index_min=np.argmin(abs(row1-f_min))
index_max=np.argmin(abs(row1-f_max))

row1=row1[index_min-bin:index_max+bin]
row2=row2[index_min-bin:index_max+bin]

plt.clf()
plt.plot(row1,row2)
plt.title('plot_zoom '+label0)
plt.xlabel('frequency(kHz)')
plt.show()


f_avg0, f_dev0=smooth(row1,bin)
B1_avg0, B1_dev0=smooth(row2,bin)



f_avg=[]
f_dev=[]
B1_avg=[]
B1_dev=[]

if bin==1:
    f_avg=f_avg0
    f_dev=f_dev0
    B1_avg=B1_avg0
    B1_dev=B1_dev0
else: 
    for i in range(len(f_avg0)):
        if i%(bin/2)==0:
            f_avg.append(f_avg0[i])
            f_dev.append(f_dev0[i])
            B1_avg.append(B1_avg0[i])
            B1_dev.append(B1_dev0[i])
#print(B1_avg)

#print("sim_max="+str(sim_max))

print(f_dev)

f_avg_list,f_std_list,x_avg_list,x_std_list=smooth_f_x(Exp_data['delta_br_vs_f'],Exp_data['x'],np.mean(f_dev))

d = {'f(kHz)':x_avg_list,'f_err(kHz)':x_std_list,'B_R(Gauss)':f_avg_list,'B_R_err(Gauss)':f_std_list}
df_exp=pd.DataFrame(d, columns=['f(kHz)','f_err(kHz)','B_R(Gauss)','B_R_err(Gauss)'])
df_exp.to_csv('0_exp_B_r_f_smooth.csv',index=False)

d = {'f(kHz)':np.array(f_avg),'f_err(kHz)':f_dev,'B_R(Gauss)':B1_avg,'B_R_err(Gauss)':B1_dev}
df=pd.DataFrame(d, columns=['f(kHz)','f_err(kHz)','B_R(Gauss)','B_R_err(Gauss)'])
df.to_csv('0B_r_f_smooth.csv',index=False)

print('-1.*np.array(f_avg)'+str(-1.*np.array(f_avg)))


if norm_to_max==True:
    print('norm_to_max==True')
    #**********for normalized to max = 1***********
    df_exp['B_R(Gauss)']=df_exp['B_R(Gauss)']/np.max(df_exp['B_R(Gauss)'])
    df['B_R(Gauss)']=df['B_R(Gauss)']/np.max(df['B_R(Gauss)'])*sim_scale
    print('sim_scale'+str(sim_scale))
    print('df[B_R(Gauss)]'+str(df['B_R(Gauss)']))
    #25% uncertainty. 
    df_exp['B_R_err(Gauss)']=df_exp['B_R(Gauss)']*0.25/np.max(df_exp['B_R(Gauss)'])
    df['B_R_err(Gauss)']=df['B_R_err(Gauss)']/np.max(df['B_R(Gauss)'])*sim_scale
    #********for normalized to  max = 1************
else: 
    #*******the y_plot_sim=y_sim*sim_scale*********
    df['B_R(Gauss)']=df['B_R(Gauss)']*sim_scale
    df['B_R_err(Gauss)']=df['B_R_err(Gauss)']*sim_scale
    #*******the y_plot_sim=y_sim*sim_scale*********

sim_max=np.max(df_exp['B_R(Gauss)'])*plot_scale



plt.clf()

ax=df_exp.plot(kind='scatter',x='f(kHz)',xerr='f_err(kHz)',y='B_R(Gauss)',yerr='B_R_err(Gauss)',\
       ylim=(0.00,sim_max),xlim=(f_min,f_max),\
       color='red',label='Experiment',alpha=0.5)
df.plot(kind='scatter',x='f(kHz)',y='B_R(Gauss)',\
      ylim=(0.00,sim_max),xlim=(f_min,f_max),\
      secondary_y=True,color='green',label='GENE simulation',ax=ax)

#df.plot(kind='scatter',x='f(kHz)',xerr='f_err(kHz)',y='B_R(Gauss)',yerr='B_R_err(Gauss)',\
#     ylim=(0.00,sim_max),xlim=(-f_max,-f_min),\
#      secondary_y=True,color='green',label='GENE simulation',ax=ax)


ax.set_xlabel(r'$frequency(kHz)$',fontsize=15)
ax.set_ylabel(r'$\bar{B}_r(Gauss/\sqrt{Hz})$',fontsize=15)

plt.title(r'$\bar{B}_r$ spectrogram',fontsize=20)
plt.savefig('B1_f.png')
plt.show()



#********RIP***********

