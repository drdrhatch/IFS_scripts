#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


class ky_mode(object):
    """Class for organizing ballooning structure for each ky mode"""

    def __init__(self, ky, field, pars):
        self.ky = ky
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
        ncon = (self.kx_modes.size - 1) // 2
        ncon1 = ncon + 1
        self.zgrid = np.pi * np.linspace(
            -(ncon + 1), ncon + 1, (ncon + 1) * self.nz, endpoint=False
        )

    def define_phase(self):
        self.phase = (-1) ** (self.ky * self.N)
