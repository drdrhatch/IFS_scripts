import numpy as np

#-rw-r----- 1 drhatch drhatch 13589545006 Aug  9 20:27 1_checkpoint
#-rw-r----- 1 drhatch drhatch 13589545006 Aug 10 10:36 43_s_checkpoint
#-rw-r----- 1 drhatch drhatch 13589545006 Aug 10 13:06 47_s_checkpoint
#-rw-rw---- 1 drhatch drhatch 13589545006 Aug 10 19:22 checkpoint
#-rw-rw---- 1 drhatch drhatch 13589545006 Aug 10 16:59 s_checkpoint

file_name = 's_checkpoint'   
f = open(file_name,'rb')
f.seek(6)
time=np.fromfile(f,dtype='float64',count=1)
f.close()
print(time)

