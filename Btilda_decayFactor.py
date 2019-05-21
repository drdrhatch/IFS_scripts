import sys
from parIOHelper import *
from geomHelper import *
import numpy as np

suffix = sys.argv[1]
if not suffix == '.dat':
    suffix = '_' + suffix

efit_file_name = sys.argv[2]

ktheta_cm = 0.65
calcRZ = False # calculate Mirnov coil pos. based on exp.

pars = init_read_parameters(suffix)

geom_type, geom_pars, geom_coeff = init_read_geometry(suffix, pars)

R, Z = local_grid_points(geom_coeff, False)

ky_fluxsurface = ky(pars, geom_coeff, ktheta_cm, False)

if calcRZ:    # if this is the local linear run at rhot = 0.999...
    topInd = np.argmax(Z)
    botInd = np.argmin(Z)
    Z1 = Z[botInd : topInd + 1]
    mirnovZ = 0.1
    zInd = np.argmin(abs(Z1 - mirnovZ)) + botInd
    mirnovZ = Z[zInd]
    mirnovR = R[zInd]
else:    # if the probe's position at LCFS is found
    mirnovR = 0.86669 + 0.02
    mirnovZ = 0.1015
print('Mirnov probe is at ({:.4f}, {:.4f})'.format(mirnovR, mirnovZ))

# go along the surface to find the point closest to Mirnov Probe
minDist = 99.
for i in range(len(R)):
    if calcRZ:
        dist = np.sqrt((R[i]-mirnovR)**2 + (Z[i]-mirnovZ)**2)+ 0.02
    else:
        dist = np.sqrt((R[i]-mirnovR)**2 + (Z[i]-mirnovZ)**2)
    if dist < minDist:
        minDist = dist
        minDR = R[i]
        minDZ = Z[i]
        minDInd = i
print 'Point on flux surface closest to the Mirnov Probe is:\n({:.4f}, {:.4f})'.format(minDR, minDZ)
print 'Its ky * d is {:.4f} '.format(minDist * ky_fluxsurface[minDInd])
#print 'index =', minDInd

# go along the surface to find the point with minimum ky*distance
minArg = 99.
for i in range(len(R)):
    if calcRZ:
        dist = np.sqrt((R[i]-mirnovR)**2 + (Z[i]-mirnovZ)**2)+ 0.02
    else:
        dist = np.sqrt((R[i]-mirnovR)**2 + (Z[i]-mirnovZ)**2)
    kyD = ky_fluxsurface[i] * dist
    if kyD < minArg:
        minArg = kyD
        minR = R[i]
        minZ = Z[i]
        minInd = i
print 'Point on flux surface minimizes ky * d to the Mirnov Probe is:\n({:.4f}, {:.4f})'.format(minR, minZ)
decayFactor = np.exp(-minArg)
print 'Its ky * d is {:.4f}, decay factor is {:.4f}'.format(minArg, decayFactor)
#print 'index =', minInd
