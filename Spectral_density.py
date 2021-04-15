import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

def sort_x_f(x_unsort,f_unsort): 
   
    arr_unsort=[x_unsort,f_unsort]
    f_x_unsort=tuple(map(tuple, np.transpose(arr_unsort)))
      
    f_x_sort=sorted(f_x_unsort, key=lambda f_x_unsort: f_x_unsort[0])
    f_x_sort=np.array(f_x_sort)
    f_x_sort=np.transpose(f_x_sort)
    x_sort=f_x_sort[0,:]
    f_sort=f_x_sort[1,:]

    return x_sort,f_sort

#https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.periodogram.html
def spectral_density(function,time,plot=False):
	time=np.array(time)
	dt=time[1:]-time[:-1]
	dt_min=np.min(dt)

	if abs(np.std(dt))>=np.min(dt)*0.01:
		print('time step is NOT uniform. interperlating')
		uni_time = np.linspace(min(time),max(time),int(abs((max(time)-min(time))/dt_min)*1.5)) #uniform time
		uni_function = np.interp(uni_time,time,function)
	else:
		uni_time=time
		uni_function=function


	fs=1./np.mean(abs(uni_time[1:]-uni_time[:-1]))
	print('avg_dt='+str(np.mean(abs(uni_time[1:]-uni_time[:-1]))))
	print('std_dt='+str(np.std(abs(uni_time[1:]-uni_time[:-1]))))
	f, Pxx_den = signal.periodogram(uni_function, fs) #, scaling='spectrum')

	#Sort frequency to monotonic increase
	f, Pxx_den=sort_x_f(f, Pxx_den)

	if plot==True:
		plt.plot(f, Pxx_den,label='Pxx_den')
		plt.plot(f, np.sqrt(Pxx_den),label='sqrt(Pxx_den)')
		#plt.semilogy(f, Pxx_den)
		#plt.ylim([1e-7, 1e2])
		plt.xlabel('frequency [Hz]')
		plt.ylabel('PSD [V**2/Hz]')
		plt.legend()
		plt.show()
	return f, Pxx_den

'''
#*********Demo function****************
timestep0=0.002
time1 = np.arange(0.,2.,timestep0)
time2 = np.arange(2.00001,3.,timestep0*0.1)

#time = time1
time = np.append(time1, time2)
print(str(time))
#time.append(2.001)
#time.extend(time)

frequency=20.
omega = 2.*np.pi*frequency
#function=1+np.exp(-1.j * omega * time -3.j)+0.5*np.exp(+1.j * 2*omega * time-1.j)+0.5*np.exp(+1.j * 3.*omega * time-1.j)

function=np.exp(-1.j * omega * time )+0.5*np.exp(+1.j * 2*omega * time)+2.*np.exp(+1.j * 3.*omega * time)


f, Pxx_den=spectral_density(function,time,plot=True)

Sum_over_t=np.mean(abs(function)**2.)**0.5
Sum_over_f=(np.sum(Pxx_den)*np.mean(abs(f[:-1]-f[1:])))**0.5


print('f='+str(f))
print('Pxx_den='+str(Pxx_den))
print('np.sum(Pxx_den)='+str(np.sum(Pxx_den)))
print('np.mean(abs(f[:-1]-f[1:]))='+str(np.mean(abs(f[:-1]-f[1:]))))


print('Sum_over_t='+str(Sum_over_t))
print('Sum_over_f='+str(Sum_over_f))

#*********Demo function****************



'''