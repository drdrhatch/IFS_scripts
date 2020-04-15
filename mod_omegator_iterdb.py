import numpy as np
import matplotlib.pyplot as plt
from read_iterdb import *
from write_iterdb import *


filename = 'AUG_B1.iterdb'
shotnum = '82585'
VROT_factor = 0.0
sign = 1
omtor_str = 'Omtor_factor_'+str(VROT_factor)
reverse_sign = 1 #-1 for reverse, 1 otherwise

rhot,profs,units = read_iterdb(filename)

include_impurity = False
if 'NM2' in list(profs.keys()):
   include_impurity = True

maxVROT = np.max(abs(profs['VROT']))
print("maxVROT",maxVROT)
offset = np.empty(len(profs['VROT']))
offset[:] = 1.0
VROT_out = reverse_sign*profs['VROT'] + offset*maxVROT*VROT_factor*sign

print(list(profs.keys()))
for x in list(rhot.keys()):
  diff =  np.sum(np.abs(rhot[x]-rhot['TE']))
  print("diff",diff)
  if diff > len(x)*1.0e-14:
     stop

plt.plot(rhot['VROT'],profs['VROT'],label = 'Original VROT')
plt.plot(rhot['VROT'],VROT_out,label = 'New VROT')
plt.legend()
plt.show()

if include_impurity:
   output_iterdb(rhot['TE'],rhot['TE'],profs['NE']/1e19,profs['TE']/1000.0,profs['NM1']/1e19,profs['TI']/1000.0,filename[:-7]+'_'+omtor_str,shotnum,'999',vrot=VROT_out,nimp=profs['NM2']/1.0e19)
else:
   output_iterdb(rhot['TE'],rhot['TE'],profs['NE']/1e19,profs['TE']/1000.0,profs['NM1']/1e19,profs['TI']/1000.0,filename[:-7]+'_'+omtor_str,shotnum,'999',vrot=VROT_out)


