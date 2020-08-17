#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt


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
        self.zgrid = np.linspace(-nxmodes, nxmodes, nxmodes * self.nz, endpoint=False)
        self.zero_ind = self.zgrid.size // 2

    def define_phase(self, pars):
        if "n0_global" in pars:
            phase = np.e ** (-2 * np.pi * 1j * pars["n0_global"] * pars["q0"])
        else:
            phase = -1
        step = max(1, max(self.kx_modes))
        self.phase = phase ** (self.kx_modes / step)

    def read_phi(self):
        """ Read phi for a given time window, returning array
        """
        tmp = (self.field.phi()[:, self.ky, self.kx_modes] * self.phase).ravel(
            order="F"
        )
        if hasattr(self, "phi"):
            self.phi = np.vstack([self.phi, tmp])
        else:
            self.phi = [tmp]

    def pod(self):
        self.u, self.s, self.vh = la.svd(self.phi)

    def plot_modes(self, times):
        for phi, time in zip(self.phi, times):
            plt.title(
                r"$\phi$, $k_y=$" + str(self.ky) + " t = " + str("{:6.3f}").format(time)
            )
            norm = phi[self.zero_ind]
            phi /= norm
            self.plot(phi)

    def plot_pod(self, pods):
        for pod in pods:
            plt.title(
                r"$\phi$, $k_y=$" + str(self.ky) + ", POD mode # = " + str(pod + 1)
            )
            phi = np.conjugate(self.vh[pod])
            norm = phi[self.zero_ind]
            phi /= norm
            self.plot(phi)

    def plot(self, var):
        plt.plot(self.zgrid, np.real(var), color="red", label=r"$Re[\phi]$")
        plt.plot(self.zgrid, np.imag(var), color="blue", label=r"$Im[\phi]$")
        plt.plot(self.zgrid, np.abs(var), color="black", label=r"$|\phi|$")
        plt.legend()
        plt.xlabel(r"$z/\pi$", size=18)
        plt.show()


def get_times(field, stime, etime):
    tarray = np.array(field.tfld)
    tind = (stime < tarray) * (tarray < etime)
    return tarray[tind]

