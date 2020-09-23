#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from operator import attrgetter


class ky_mode(object):
    """Class for organizing ballooning structure for each ky mode"""

    varnames = {
        "phi": "$\Phi$",
        "apar": "$A_\parallel$",
        "bpar": "$B_\parallel$",
        "tperp": "$T_\perp$",
        "tpar": "T_\parallel$",
        "dens": "$n$",
        "q": "$Q$",
    }

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
        if self.ky == 0:
            step = 1
        else:
            step = self.N * self.ky
        hmodes = np.arange(0, self.nx / 2, step, dtype=np.intc)
        lmodes = np.arange(0, -self.nx / 2, -step, dtype=np.intc)
        self.kx_modes = np.union1d(lmodes, hmodes)

    def zrange(self):
        nxmodes = self.kx_modes.size
        self.zgrid = np.linspace(-1, 1, self.nz, endpoint=False)
        self.zgrid_ext = np.linspace(
            -nxmodes, nxmodes, nxmodes * self.nz, endpoint=False
        )
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
                    "dens": self.mom.dens,
                    "tpar": self.mom.tpar,
                    "tperp": self.mom.tperp,
                }
            )
        fields = ["phi", "apar", "bpar", "dens", "tpar", "tperp", "q"]
        self.fields = dict.fromkeys(fields, None)

    def read_field(self, varname):
        """ Read field for a given time window, returning array"""
        var = self.field_vars[varname]()
        if var.shape[1] == 1:  # for linear scan data with single ky
            indy = 0
        else:
            indy = self.ky
        tmp = var[:, indy, :]
        return tmp

    def read_fields(self, times, fields):
        """Read given fields data for the given times"""
        self.fields_read = fields
        tmp = np.empty((len(fields), times.size, self.nz, self.nx), dtype=np.cdouble)
        for i, var in enumerate(fields):
            self.fields[var] = tmp[i, :, :]
        for j, time in enumerate(times):
            print("Reading fields at time t = " + str("{:6.3f}").format(time))
            self.field.set_time(time)
            self.mom.set_time(time)
            for i, var in enumerate(fields):
                tmp[i, j, :, :] = self.read_field(var)
        self.define_variables()

    def define_variables(self):
        self.phi = self.fields["phi"]
        self.fields["phi2"] = np.square(self.phi)
        self.phi2 = self.fields["phi2"]
        self.apar = self.fields["apar"]
        self.bpar = self.fields["bpar"]
        self.dens = self.fields["dens"]
        self.tpar = self.fields["tpar"]
        self.tperp = self.fields["tperp"]

    def pod(self, var):
        u, sv, vh = la.svd(var)
        return u, sv, vh

    def construct_q(self, times):
        req_fields = {"phi", "tpar", "tperp", "dens"}
        fields_toread = req_fields.difference(self.fields_read)
        if fields_toread:
            self.read_fields(times, fields_toread)
        phi = self.phi
        tpar = self.tpar
        tperp = self.tperp
        dens = self.dens
        tmp = -1j * self.ky * phi * np.conj(0.5 * tpar + tperp + 1.5 * dens)
        self.q = tmp + np.conj(tmp)
        self.fields["q"] = self.q

    def plot_modes(self, varname, times, extend=True):
        varlabel = ky_mode.get_varname(varname)
        for var, time in zip(self.fields[varname], times):
            plt.title(r"$k_y=$" + str(self.ky) + " t = " + str("{:6.3f}").format(time))
            if extend:
                pvar = (var[:, self.kx_modes] * self.phase).ravel(order="F")
                norm = pvar[self.zero_ind]
                zgrid = self.zgrid_ext
            else:
                pvar = var.mean(axis=-1)
                norm = var[0, 0]
                zgrid = self.zgrid
            if norm == 0:
                norm = 1
            pvar *= 1 / norm
            self.plot(zgrid, pvar, varlabel)

    def plot_pod(self, var, pods, varn):
        varname = ky_mode.get_varname(varn)
        for pod in pods:
            plt.title("$k_y=$" + str(self.ky) + ", POD mode # = " + str(pod + 1))
            pvar = np.conjugate(var[pod])
            norm = pvar[self.zero_ind]
            pvar /= norm
            self.plot(pvar, varname)

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

    def plot(self, zgrid, var, varname):
        plt.plot(zgrid, np.real(var), color="red", label=r"$Re[$" + varname + "$]$")
        plt.plot(zgrid, np.imag(var), color="blue", label=r"$Im[$" + varname + "$]$")
        plt.plot(zgrid, np.abs(var), color="black", label=r"$|$" + varname + "$|$")
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

    @classmethod
    def get_varname(cls, var):
        """returns formatted label for plots corresponding to input variable"""
        try:
            varname = ky_mode.varnames[var]
        except:
            varname = ""
        return varname


def get_times(field, stime, etime):
    """Get times between two extremes from either field or mom file"""
    try:
        tarray = np.array(field.tfld)
    except AttributeError:
        tarray = np.array(field.tmom)
    tind = (stime <= tarray) * (tarray <= etime)
    return tarray[tind]


def sum_modes(modes, varname):
    """Average variable var over modes (x & y)"""
    ntimes = modes[0].fields[varname].shape[0]
    tmp = np.empty(
        (len(modes), ntimes, modes[0].nz), dtype=modes[0].fields[varname].dtype
    )
    for i, mode in enumerate(modes):
        tmp[i] = sum_x(mode, varname)
    ysum = tmp.sum(axis=0, keepdims=False)
    return ysum


def sum_x(mode, varname):
    """Average variable over x dimension"""
    var = mode.fields[varname]
    xsum = np.sum(var, axis=-1, keepdims=False)
    return xsum
