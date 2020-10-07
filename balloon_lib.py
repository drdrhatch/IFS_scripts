#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from operator import attrgetter
from matplotlib.backends.backend_pdf import PdfPages


VARNAMES = {
    "phi": r"$\Phi$",
    "apar": r"$A_\parallel$",
    "bpar": r"$B_\parallel$",
    "tperp": r"$T_\perp$",
    "tpar": r"$T_\parallel$",
    "dens": "$n$",
    "q": "$Q$",
}


class KyMode:
    """Class for organizing ballooning structure for each ky mode"""

    def __init__(
        self, ky, pars, times, fields, field_file, mom_file=None, geom_file=None
    ):
        self.ky = ky
        self.times = times
        self.field = field_file
        self.mom = mom_file
        self.nx = field_file.nx
        self.nz = field_file.nz
        self.N = pars["nexc"]
        self.construct_ranges()
        self.define_phase(pars)
        self.define_dictionary()
        self.geometry = geom_file
        self.read_fields(times, fields)

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
        self.zero_ind = self.zgrid_ext.size // 2

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
        self.fields_read = set(fields)
        tmp = np.empty((len(fields), times.size, self.nz, self.nx), dtype=np.cdouble)
        for j, time in enumerate(times):
            self.field.set_time(time)
            if self.mom:
                self.mom.set_time(time)
            for i, var in enumerate(fields):
                tmp[i, j, :, :] = self.read_field(var)
        for i, var in enumerate(fields):
            self.fields[var] = tmp[i]
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

    def plot_pod(self, var, pods, varn):
        varname = get_varname(varn)
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


def plot(zgrid, var, varname, title):
    """Base plotting function for complex variables
    returns plot object"""
    fig = plt.figure()
    plt.title(title)
    plt.plot(zgrid, np.real(var), color="red", label=r"$Re[$" + varname + "$]$")
    plt.plot(zgrid, np.imag(var), color="blue", label=r"$Im[$" + varname + "$]$")
    plt.plot(zgrid, np.abs(var), color="black", label=r"$|$" + varname + "$|$")
    plt.legend()
    plt.xlabel(r"$z/\pi$", size=18)
    return fig


def plot_var(mode, varname, times, extend=True, show=True, output=False):
    """plot variable for mode with formatted key returns plot object"""
    varlabel = get_varname(varname)
    for var, time in zip(mode.fields[varname], times):
        title = r"$k_y=" + str(mode.ky) + ", t = " + str("{:6.3f}").format(time) + "$"
        if extend:
            pvar = (var[:, mode.kx_modes] * mode.phase).ravel(order="F")
            norm = pvar[mode.zero_ind]
            zgrid = mode.zgrid_ext
        else:
            pvar = var.mean(axis=-1)
            norm = var[0, 0]
            zgrid = mode.zgrid
        if norm == 0:
            norm = 1
        pvar *= 1 / norm
        fig = plot(zgrid, pvar, varlabel, title)
        if show:
            plt.show()
        if output:
            output.savefig(fig)
        plt.close()


def plot_vars(mode, varnames, times, extend=True, show=True, save=False):
    if save:
        pdf_figs = PdfPages("mode_" + str(mode.ky) + ".pdf")
        output = pdf_figs
    else:
        output = False
    for varname in varnames:
        plot_var(mode, varname, times, extend, show, output)
    if save:
        pdf_figs.close()


def get_varname(var):
    """returns formatted label for plots corresponding to input variable"""
    try:
        varname = VARNAMES[var]
    except KeyError:
        print("ERROR: Variable not found in dictionary")
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


def pod(mode, varname):
    var = mode.fields[varname]
    ntimes = var.shape[0]
    pvar = var[:, :, mode.kx_modes].reshape(ntimes, -1, order="F")
    u, sv, vtmp = la.svd(pvar, full_matrices=False)
    vh = vtmp.reshape(-1, mode.nz, mode.kx_modes.size, order="F")
    return u, sv, vh


# collective is (slightly, usually) different because it includes all kx modes
def collective_pod(mode, fields):
    ntimes = mode.fields[fields[0]].shape[0]
    nxnz = mode.nx * mode.nz
    all_fields = np.concatenate(
        ([mode.fields[field].reshape(ntimes, -1) for field in fields]), axis=1
    )
    u, sv, vh = la.svd(all_fields, full_matrices=False)
    VH = {}
    for i, field in enumerate(fields):
        VH[field] = vh[:, i * nxnz : (i + 1) * nxnz].reshape((-1, mode.nz, mode.nx))
    return u, sv, VH


def calc_heat_flux(ky, fields):
    phi = fields["phi"]
    tpar = fields["tpar"]
    tperp = fields["tperp"]
    dens = fields["dens"]
    tmp = -1j * ky * phi * np.conj(0.5 * tpar + tperp + 1.5 * dens)
    heat_flux = tmp + np.conj(tmp)
    return heat_flux
