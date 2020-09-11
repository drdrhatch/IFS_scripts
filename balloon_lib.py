#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt


class ky_mode(object):
    """Class for organizing ballooning structure for each ky mode"""

    def __init__(self, ky, pars, field, mom=None):
        self.ky = ky
        self.field = field
        self.mom = mom
        self.nx = field.nx
        self.nz = field.nz
        self.N = pars["nexc"]
        self.construct_ranges()
        self.define_phase(pars)
        self.define_dictionary()

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

    def define_dictionary(self):
        self.field_vars = {
            "phi": self.field.phi,
            "apar": self.field.apar,
            "bpar": self.field.bpar,
        }
        if self.mom:
            self.field_vars.update(
                {
                    "tpar": self.mom.tpar,
                    "tperp": self.mom.tperp,
                }
            )
        fields = ["phi", "apar", "bpar", "tpar", "tperp"]
        self.fields = dict.fromkeys(fields, None)

    def read_field(self, var):
        """ Read field for a given time window, returning array"""
        tmp = (self.field_vars[var]()[:, self.ky, self.kx_modes] * self.phase).ravel(
            order="F"
        )
        return tmp

    def read_fields(self, times, fields):
        """Read given fields data for the given times"""
        tmp = np.empty((len(fields), times.size, self.zgrid.size), dtype=np.cdouble)
        for i, var in enumerate(fields):
            self.fields[var] = tmp[i, :, :]
        for i, time in enumerate(times):
            print("Reading fields at time t = " + str("{:6.3f}").format(time))
            self.field.set_time(time)
            for j, var in enumerate(fields):
                tmp[j, i, :] = self.read_field(var)
        self.define_variables()

    def define_variables(self):
        self.phi = self.fields["phi"]
        self.apar = self.fields["apar"]
        self.bpar = self.fields["bpar"]
        self.tpar = self.fields["tpar"]
        self.tperp = self.fields["tperp"]

    def pod(self, var):
        u, sv, vh = la.svd(var)
        return u, sv, vh

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

    def plot_singular_values(self):
        pods = range(1, self.sv.size + 1)

        fig, ax1 = plt.subplots()

        color = "red"
        ax1.set_ylabel("value", color=color)
        ax1.tick_params(axis="y", labelcolor=color)
        ax1.plot(pods, self.sv, marker="o", color=color)

        ax2 = ax1.twinx()

        cs = np.cumsum(self.sv) / self.sv.sum()
        color = "blue"
        ax2.plot(pods, cs, color=color)
        ax2.set_ylim(0, 1.0)
        ax2.set_ylabel("cumulative", color=color)
        ax2.tick_params(axis="y", labelcolor=color)
        ax2.grid()

        plt.title(r"Singular values for mode $k_y = $" + str(self.ky))
        plt.xlabel("POD #", size=18)
        plt.xticks(pods)
        plt.grid(True)
        plt.show()

    def plot_time_dependence(self, times, pods):
        plt.title(r"Time dependece of POD modes")
        plt.xlabel("Time")
        plt.ylabel(r"$|\Phi_s|$")
        # plt.xticks(pods)

        for pod in pods:
            plt.plot(times, np.abs(self.u[:, pod]), label=r"$s_" + str(pod + 1) + "$")
        plt.grid(True)
        plt.legend()
        plt.show()

    def plot(self, var):
        plt.plot(self.zgrid, np.real(var), color="red", label=r"$Re[\phi]$")
        plt.plot(self.zgrid, np.imag(var), color="blue", label=r"$Im[\phi]$")
        plt.plot(self.zgrid, np.abs(var), color="black", label=r"$|\phi|$")
        plt.legend()
        plt.xlabel(r"$z/\pi$", size=18)
        plt.show()

    def output(self, pods, times, norm):
        """Output various POD data"""
        self.output_sv()
        self.output_pod_modes(pods, norm)
        self.output_time_modes(pods, times)

    def output_sv(self):
        """Output singular values"""
        filename = "./sv_ky" + str("{:03d}").format(self.ky) + ".dat"
        header = "Singular values"
        np.savetxt(filename, self.sv, fmt="%g", header=header, encoding="UTF-8")

    def output_pod_modes(self, pods, norm):
        """Output right pod modes (spatial variation)"""
        if norm:
            filename = "./pod_ky" + str("{:03d}").format(self.ky) + "_norm.dat"
        else:
            filename = "./pod_ky" + str("{:03d}").format(self.ky) + ".dat"
        fp = open(filename, "w")
        for pod in range(pods):
            header = str(pod + 1)
            phi = self.vh[pod]
            if norm:
                phi /= phi[self.zero_ind]
            data = np.vstack((self.zgrid, np.real(phi), np.imag(phi))).T
            np.savetxt(
                fp,
                data,
                fmt="%g",
                header=header,
                encoding="UTF-8",
            )
            fp.write("\n\n")
        fp.close()

    def output_time_modes(self, pods, times):
        """Output left pod modes (time variation)"""
        filename = "./pod_time_ky" + str("{:03d}").format(self.ky) + ".dat"
        head = ["time"]
        for pod in range(pods):
            head.append(str(pod + 1))
        header = " ".join(head)
        data = np.hstack((times.reshape(-1, 1), np.abs(self.u[:, :pods])))
        np.savetxt(
            filename,
            data,
            fmt="%g",
            header=header,
            encoding="UTF-8",
        )

    def write_data(self, filename, data, indices):
        """Write data in text format for later plotting and analysis"""
        pass


def get_times(field, stime, etime):
    tarray = np.array(field.tfld)
    tind = (stime < tarray) * (tarray < etime)
    return tarray[tind]
