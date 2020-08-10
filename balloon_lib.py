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
        self.define_phase()

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

    def define_phase(self):
        self.phase = (-1) ** abs(self.kx_modes)

    def read_phi(self):
        self.phi = (self.field.phi()[:, self.ky, self.kx_modes] * self.phase).ravel(
            order="F"
        )
        complex_phase = abs(self.phi[self.zero_ind]) / self.phi[self.zero_ind]
        self.phi *= complex_phase
