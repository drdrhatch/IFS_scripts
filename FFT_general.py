import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

#Created by Max Curie: 05/15/2021
#GitHub: https://github.com/maxcurie1996/Python_Demo/tree/main/FFT
#Un-uniform FFT: https://scicomp.stackexchange.com/questions/593/how-do-i-take-the-fft-of-unevenly-spaced-data

def sort_x_f(x_unsort,f_unsort): 
   
    arr_unsort=[x_unsort,f_unsort]
    f_x_unsort=tuple(map(tuple, np.transpose(arr_unsort)))
      
    f_x_sort=sorted(f_x_unsort, key=lambda f_x_unsort: f_x_unsort[0])
    f_x_sort=np.array(f_x_sort)
    f_x_sort=np.transpose(f_x_sort)

    x_sort=f_x_sort[0,:]
    f_sort=f_x_sort[1,:]
    x_sort=x_sort.astype(type(x_unsort[0]))
    f_sort=f_sort.astype(type(f_unsort[0]))

    #print('2')
    #print('type(x_sort[0])'+str(type(x_sort[0])))
    #print('type(f_sort[0])'+str(type(f_sort[0])))

    return x_sort,f_sort


def FFT_function_time(function,time,plot=False): 
    time=np.array(time)
    time,function=sort_x_f(time,function)

    dt=time[1:]-time[:-1]
    dt_min=np.mean(dt)
    
    if abs(np.std(dt))>=np.min(dt)*0.01:
        print('time step is NOT uniform. interperlating')
        uni_time = np.linspace(min(time),max(time),int(abs((max(time)-min(time))/dt_min)*1.5)) #uniform time
        uni_function = np.interp(uni_time,time,function)
    else:
        uni_time=time
        uni_function=function


    timestep=np.mean(abs(uni_time[1:]-uni_time[:-1]))
    print('avg_dt='+str(np.mean(abs(uni_time[1:]-uni_time[:-1]))))
    print('std_dt='+str(np.std(abs(uni_time[1:]-uni_time[:-1]))))
    norm=1./float(len(uni_time))  #normalizing factor
    amplitude_complex = np.fft.fft(uni_function)
    #print(str(time.shape[-1]))
    #output_x = np.fft.fftfreq(t.shape[-1])
    frequency = np.fft.fftfreq(uni_time.shape[-1], d=timestep)
    amplitude_frequency=abs(norm*amplitude_complex)
    #amplitude_frequency=norm*amplitude_complex.real
    phase_frequency=np.angle(amplitude_complex)


    #Sort frequency to monotonic increase
    frequency_sort,amplitude_frequency_sort=sort_x_f(frequency,amplitude_frequency)
    frequency_sort,phase_frequency_sort=sort_x_f(frequency,phase_frequency)

    if plot==True:
        #plt.plot(frequency,amplitude_frequency)
        plt.plot(frequency_sort,amplitude_frequency_sort)
        #plt.semilogy(f, Pxx_den)
        #plt.ylim([1e-7, 1e2])
        plt.xlabel('frequency [Hz]')
        plt.ylabel('amplitude')
        plt.grid()
        plt.legend()
        plt.show()
    return frequency_sort,amplitude_frequency_sort,phase_frequency_sort

def FFT_sum(f,amp_f,frequency_min,frequency_max,frequency_all):
    sum0_TEMP=0.
    if frequency_all==True:
        frequency_min=np.min(f)
        frequency_max=np.max(f)

    for i_f in range(len(f)):
        if frequency_min<=f[i_f] and f[i_f]<=frequency_max and i_f-1>=0:
            sum0_TEMP=sum0_TEMP+abs(amp_f[i_f])
    sum0=sum0_TEMP
    sum0_error=0.
    return sum0,sum0_error


def FFT_interp(frequency,amplitude_frequency,total_len_scale=1.5):
    frequency_sort,amplitude_frequency_sort=sort_x_f(frequency,amplitude_frequency)
    uni_frequency_sort=np.linspace(np.min(frequency_sort),np.max(frequency_sort),int(len(frequency_sort)*total_len_scale))
    uni_amplitude_frequency_sort = np.interp(uni_frequency_sort,frequency_sort,amplitude_frequency_sort)*float(len(frequency_sort))/float(len(uni_frequency_sort))
    return uni_frequency_sort,uni_amplitude_frequency_sort

#About welch method: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.welch.html
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.periodogram.html
#'boxcar' Also known as a rectangular window or Dirichlet window, this is equivalent to no window at all.
#window types: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.get_window.html#scipy.signal.get_window
#intruction video of Welch's method: https://youtu.be/YK1F0-3VvQI
def spectral_density(function,time,percent=0.5,window_for_FFT='hann',plot=False):
    time=np.array(time)
    
    time,function=sort_x_f(time,function)
    dt=time[1:]-time[:-1]
    dt_min=np.mean(dt)

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
    #f, Pxx_den = signal.welch(uni_function, fs, nperseg=len(uni_function), window=window_for_FFT) #, scaling='spectrum')

    f, Pxx_den = signal.welch(uni_function, fs, nperseg=int(percent*len(uni_function)), window=window_for_FFT,return_onesided=False, scaling='density')
    #f, Pxx_den = signal.periodogram(uni_function, fs)

    #Sort frequency to monotonic increase
    f, Pxx_den=sort_x_f(f, Pxx_den)

    if plot==True:
        plt.plot(f, Pxx_den,label='Pxx_den')
        plt.plot(f, np.sqrt(Pxx_den),label='sqrt(Pxx_den)')
        #plt.semilogy(f, Pxx_den)
        #plt.ylim([1e-7, 1e2])
        plt.xlabel('frequency [Hz]')
        plt.ylabel('PSD [V**2/Hz]')
        plt.grid()
        plt.legend()
        plt.show()
    return f, Pxx_den


def spectral_density_sum(f,amp_f,frequency_min,frequency_max,frequency_all):
    sum0_TEMP=0.
    sum0_list=[]
    if frequency_all==True:
        frequency_min=np.min(f)
        frequency_max=np.max(f)

    for i_f in range(len(f)):
        if frequency_min<=f[i_f] and f[i_f]<=frequency_max and (i_f-1)>=0:
            sum0_TEMP=sum0_TEMP+abs(amp_f[i_f])**2.*abs(f[i_f]-f[i_f-1])
            sum0_list.append(sum0_TEMP)

    sum0=(sum0_TEMP)**0.5
    sum0_error=0.
    return sum0,sum0_error


def spectral_density_interp(frequency,amplitude_frequency,total_len_scale=1.5):
    frequency_sort,amplitude_frequency_sort=sort_x_f(frequency,amplitude_frequency)
    uni_frequency_sort=np.linspace(np.min(frequency_sort),np.max(frequency_sort),int(len(frequency_sort)*total_len_scale))
    uni_amplitude_frequency_sort = np.interp(uni_frequency_sort,frequency_sort,amplitude_frequency_sort)
    return uni_frequency_sort,uni_amplitude_frequency_sort

def gaussian_max(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean) / (1.*stddev))**2.)
    
def test_functions(function_num=1):
    fs = 10e3*100
    N = 1e5
    
    time = np.arange(N) / fs
    if function_num==1:
        amp = 3
        freq = 1000.0
        np.random.seed(1234)
        noise_power = 0.001 * fs / 2
        function = amp*np.sin(2*np.pi*freq*time)
        function += np.random.normal(scale=np.sqrt(noise_power), size=time.shape)

    elif function_num==2: 
        omega=2.*np.pi*200
        function=1+np.exp(-1.j * omega * time -3.j)+0.5*np.exp(+1.j * 2*omega * time-1.j)

    elif function_num==3:
        omega=2.*np.pi*200
        omega_list=np.arange(omega*0.2,omega*1.9,0.1)
        x_list=np.linspace(-2,2,len(omega_list))
        mode_list=gaussian_max(x_list, 1, 0, 0.5)
        for i in range(len(omega_list)):
            omega=omega_list[i]
            mode = mode_list[i]
            if i==0:
                function = mode*np.exp(-1.j * omega * time )
            else:
                function += mode*np.exp(-1.j * omega * time )
    if 0==1:
        plt.clf()
        plt.plot(time,function)
        plt.show()

    return function, time


#*********Demo function****************
timestep0=0.00002
time1 = np.arange(0.,2.,timestep0)
time2 = np.arange(2.00001,3.,timestep0*0.1)

#time = time1
time = np.append(time1, time2)
print(str(time))
#time.append(2.001)
#time.extend(time)

time=np.array(time)

frequency=20.
omega = 2.*np.pi*frequency
function=1.+np.exp(-1.j * omega * time -3.j)+0.5*np.exp(+1.j * 2*omega * time-1.j)
#function=np.exp(-1.j * omega * time -3.j)+0.2*np.exp(-1.j *3* omega * time -3.j)+0.5*np.exp(+1.j * 2*omega * time-1.j)

#function, time=test_functions(2)

#uni_time = np.linspace(min(time),max(time),int(len(time)*1.5)) #uniform time
#uni_function = np.interp(uni_time,time,function)


'''
#*********Demo function****************

frequency,amplitude_frequency,amplitude_growth = FFT_function_time(function,time,plot=True)
sum0,sum0_error=FFT_sum(frequency,amplitude_frequency,0,2,True)
print('*********************')
print('sum0,sum0_error='+str(sum0)+', '+str(sum0_error))
print('*********************')

frequency,amplitude_frequency = spectral_density(function,time,percent=1.,window_for_FFT='hann',plot=True)
amplitude_frequency=np.sqrt(amplitude_frequency)
amplitude_frequency
sum0,sum0_error=spectral_density_sum(frequency,amplitude_frequency,0,2,True)
print('*********************')
print('sum0,sum0_error='+str(sum0)+', '+str(sum0_error))
print('*********************')

'''