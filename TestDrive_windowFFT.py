#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from windowFFT import *
import sys

show_plots = True
plot_format = 'display'
zi = complex(0,1)
nt = 500
nf = 200
lf = 60.

tgrid = np.linspace(0., 2., num = nt)
#testField1 = np.exp(zi * np.pi * tgrid) + np.exp(zi * np.pi * tgrid * 2)
testField1 = np.exp(zi * tgrid * 9.) * 2. + np.exp(zi * tgrid * 25.) * 3
fgrid, dens_f = windowFFT(tgrid, testField1, nf, lf, 'test', show_plots, plot_format)
