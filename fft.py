#!/usr/bin/env python
# -*- coding: utf-8 -*-

# FFT routines for nonstandard uses

import numpy as np
from nfft import nfft_adjoint
import matplotlib.pyplot as plt

def fft_freq(f, t):
    ntimes = t.size//2*2
    timestep = (t[-1]-t[0])/(ntimes-1)
    omegas = 2*np.pi*np.fft.fftfreq(ntimes, d=timestep)
    dom_omega = np.empty(f.shape[1])
    for i, vec in enumerate(f.T):
        vec_hat = nfft_adjoint(t, vec, ntimes)
        dom_omega[i] = omegas[np.argmax(abs(vec_hat[:]))]  # skip freq zero as dominant
    print(omegas)
    plt.plot(dom_omega, marker="o")
    plt.show()

