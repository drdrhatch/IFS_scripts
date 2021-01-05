import pandas as pd 


RIP_file_path='RIP_csv'        #path one stores 0RIP_from_all_freq.csv, default: 'RIP_csv'
RIP_file_name='0RIP_from 150.0kHz to 500.0kHz.csv'	#the name of the RIP csv file, default: '0RIP_from_all_freq.csv'

BES_file_path='.'        	#path one stores 0BES.txt, default: '.'  (current folder)
BES_file_name='0BES.txt'	#the name of the BES text file, default: '0BES.txt'


RIP_info=pd.read_csv(RIP_file_path+'/'+RIP_file_name) 

file=open(BES_file_path+'/'+BES_file_name,"r")
line0=file.readline() 
print(line0)
n1_BES=float(line0[len('n1_BES='):-len('/m^3')-1])
print(n1_BES)
line1=file.readline() 
print(line1)
n0_BES=float(line1[len('n0='):-len('*10^19/m^3')-1])*10.**19.
print(n0_BES)

n1_BES_list=[n1_BES]*len(RIP_info['Z(cm)'])
n0_BES_list=[n0_BES]*len(RIP_info['Z(cm)'])
Ratio_list=(RIP_info['B_R(Gauss)']*n0_BES_list) / (RIP_info['B0(Gauss)']*n1_BES_list)

d = {'Z(cm)':RIP_info['Z(cm)'],'Z_err(cm)':RIP_info['Z_err(cm)'],\
	'B_R(Gauss)':RIP_info['B_R(Gauss)'],'B_R_err(Gauss)':RIP_info['B_R_err(Gauss)'],\
	'B0(Gauss)':RIP_info['B0(Gauss)'],'n1(/m^3)':n1_BES_list,'n0(/m^3)':n0_BES_list,\
	'B1*n0/B0*n1':Ratio_list}
df=pd.DataFrame(d, columns=['Z(cm)','Z_err(cm)','B_R(Gauss)','B_R_err(Gauss)','B0(Gauss)','n1(/m^3)','n0(/m^3)','B1*n0/B0*n1'])
df.to_csv('0B_n_Ratio.csv',index=False)
