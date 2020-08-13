#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


class ky_mode(object):
    """Class for organizing ballooning structure for each ky mode"""

    def __init__(self, ky, field, pars):
        self.ky = ky
        self.field = field
        self.nx = field.nx
        self.nz = field.nz
        self.N = pars["nexc"]
        self.construct_ranges()
        self.define_phase(pars)

    def construct_ranges(self):
        self.kxrange()
        self.zrange()

    def kxrange(self):
        hmodes = np.arange(0, self.nx / 2, self.N * self.ky, dtype=np.intc)
        lmodes = np.arange(0, -self.nx / 2, -self.N * self.ky, dtype=np.intc)
        self.kx_modes = np.union1d(lmodes, hmodes)

    def zrange(self):
        nxmodes = self.kx_modes.size
        self.zgrid = np.pi * np.linspace(
            -nxmodes, nxmodes, nxmodes * self.nz, endpoint=False
        )
        self.zero_ind = self.zgrid.size // 2

    def define_phase(self, pars):
        if "n0_global" in pars:
            phase = np.e ** (-2 * np.pi * 1j * pars["n0_global"] * pars["q0"])
        else:
            phase = -1
        self.phase = phase ** (self.kx_modes / max(self.kx_modes))

    def read_phi(self, stime, etime):
        """ Read phi for a given time window, returning array
        """
        if hasattr(self, "phi"):
            pass
        else:
            times = self.get_times(stime, etime)
            print(times)
            self.phi = np.empty([times.size, self.zgrid.size], dtype=np.complex128)
            for i in range(times.size):
                self.field.set_time(times[i])
                self.phi[i, :] = (
                    self.field.phi()[:, self.ky, self.kx_modes] * self.phase
                ).ravel(order="F")
                complex_phase = (
                    abs(self.phi[i, self.zero_ind]) / self.phi[i, self.zero_ind]
                )
                self.phi[i, :] *= complex_phase

    def get_times(self, stime, etime):
        tarray = np.array(self.field.tfld)
        tind = (stime < tarray) * (tarray < etime)
        return tarray[tind]
