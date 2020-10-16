#!/usr/bin/env python
# -*- coding: utf-8 -*-

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

try:
    import scipy.linalg as la
except ImportError:
    import numpy.linalg as la
import numpy as np

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
        self.iky = ky
        self.ky = ky * pars["kymin"]
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
            step = self.N * self.iky
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
            indy = self.iky
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


def plot_pod(mode, var, pods, varn, extend=True):
    varname = get_varname(varn)
    for ipod in pods:
        title = "$k_y=$" + str(mode.ky) + ", POD mode # " + str(ipod + 1)
        pvar, zgrid = get_plot_variable(mode, var[ipod], extend)
        plot(zgrid, np.conj(pvar), varname, title)
        plt.show()


def plot_time_dependence(mode, times, pods):
    plt.title(r"Time dependece of POD modes")
    plt.xlabel("Time")
    plt.ylabel(r"$|\Phi_s|$")
    for ipod in pods:
        plt.plot(times, np.abs(mode.u[:, ipod]), label=r"$s_" + str(ipod + 1) + "$")
    plt.grid(True)
    plt.legend()
    plt.show()


def output_pod(mode, u, sv, vh, fields, pods, times):
    """Output various POD data"""
    output_sv(mode, sv)
    output_pod_modes(mode, vh, fields, pods, norm=True)
    output_time_modes(mode, u, pods, times)


def output_sv(mode, sv):
    """Output singular values"""
    filename = "./sv_ky" + str("{:03d}").format(int(mode.ky)) + ".dat"
    header = "Singular values"
    np.savetxt(filename, sv, fmt="%g", header=header, encoding="UTF-8")


def output_pod_modes(mode, r_vec, fields, pods, norm):
    """Output right pod modes (spatial variation)"""
    if norm:
        filename = "./pod_ky" + str("{:03d}").format(int(mode.ky)) + "_norm.dat"
    else:
        filename = "./pod_ky" + str("{:03d}").format(int(mode.ky)) + ".dat"
    fp = open(filename, "w")
    fp.write("# theta Re Im\n")
    for ipod in pods:
        for field in fields:
            header = field + " POD " + str(ipod + 1)
            pvar, zgrid = get_plot_variable(mode, r_vec[field][ipod], extend=True)
            if norm:
                pvar /= pvar[mode.zero_ind]
            data = np.vstack((mode.zgrid_ext, np.real(pvar), np.imag(pvar))).T
            np.savetxt(
                fp,
                data,
                fmt="% E",
                header=header,
                encoding="UTF-8",
            )
            fp.write("\n\n")
    fp.close()


def output_time_modes(mode, l_vec, pods, times):
    """Output left pod modes (time variation)"""
    filename = "./pod_time_ky" + str("{:03d}").format(int(mode.ky)) + ".dat"
    head = ["time"]
    for ipod in pods:
        head.append(str(ipod + 1))
    header = " ".join(head)
    data = np.hstack((times.reshape(-1, 1), np.abs(l_vec[:, :pods.stop])))
    np.savetxt(
        filename,
        data,
        fmt="% E",
        header=header,
        encoding="UTF-8",
    )


def plot(zgrid, var, varname, title):
    """Base plotting function for complex variables
    returns plot object"""
    fig = plt.figure()
    plt.title(title)
    plt.plot(zgrid, np.real(var), color="red", label=r"$\Re[$" + varname + "$]$")
    plt.plot(zgrid, np.imag(var), color="blue", label=r"$\Im[$" + varname + "$]$")
    plt.plot(zgrid, np.abs(var), color="black", label=r"$|$" + varname + "$|$")
    plt.legend()
    plt.xlabel(r"$z/\pi$", size=18)
    return fig


def plot_var(mode, var, varlabel, title, extend=True, show=True, output=False):
    """plot variable for mode with formatted key returns plot object"""
    pvar, zgrid = get_plot_variable(mode, var, extend)
    fig = plot(zgrid, pvar, varlabel, title)
    if show:
        plt.show()
    if output:
        output.savefig(fig)
    plt.close()


def plot_vars(mode, varnames, times, extend=True, show=True, save=False):
    """Plot a given variable from mode for given times
    By default:
    plots extended ballooning structure
    shows plot
    Can also save plot"""
    if save:
        pdf_figs = PdfPages("mode_" + str(mode.ky) + ".pdf")
        output = pdf_figs
    else:
        output = False
    for varname in varnames:
        varlabel = get_varname(varname)
        for var, time in zip(mode.fields[varname], times):
            title = (
                r"$k_y=" + str(mode.ky) + ", t = " + str("{:6.3f}").format(time) + "$"
            )
            plot_var(mode, var, varlabel, title, extend, show, output)
    if save:
        pdf_figs.close()


def plot_singular_values(mode, sv, save=False):
    pods = range(1, sv.size + 1)

    fig, ax1 = plt.subplots()

    color = "red"
    ax1.set_ylabel("value", color=color)
    ax1.tick_params(axis="y", labelcolor=color)
    ax1.plot(pods, sv, marker="o", color=color)

    ax2 = ax1.twinx()

    sv_sum = np.cumsum(sv) / sv.sum()
    color = "blue"
    ax2.plot(pods, sv_sum, color=color)
    ax2.set_ylim(0, 1.0)
    ax2.set_ylabel("cumulative", color=color)
    ax2.tick_params(axis="y", labelcolor=color)
    ax2.grid()

    plt.title(r"Singular values for mode $k_y = $" + str(mode.ky))
    plt.xlabel("POD #", size=18)
    plt.xticks(pods)
    plt.grid(True)
    plt.show()
    if save:
        pdf_figs = PdfPages("mode_" + str(mode.ky) + "_sv.pdf")
        output = pdf_figs
        output.savefig(fig)
        pdf_figs.close()
    plt.close()


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
def collective_pod(mode, fields, extend=True):
    ntimes = mode.fields[fields[0]].shape[0]
    if extend:
        nx = len(mode.kx_modes)
        all_fields = np.concatenate(
            (
                [
                    mode.fields[field][:, :, mode.kx_modes].reshape(ntimes, -1)
                    for field in fields
                ]
            ),
            axis=1,
        )
    else:
        nx = mode.nx
        all_fields = np.concatenate(
            ([mode.fields[field].reshape(ntimes, -1) for field in fields]), axis=1
        )
    nxnz = nx * mode.nz
    u, sv, vh = la.svd(all_fields, full_matrices=False)
    VH = {}
    for i, field in enumerate(fields):
        VH[field] = vh[:, i * nxnz : (i + 1) * nxnz].reshape((-1, mode.nz, nx))
    return u, sv, VH


def calc_heat_flux(ky, fields):
    phi = fields["phi"]
    tpar = fields["tpar"]
    tperp = fields["tperp"]
    dens = fields["dens"]
    tmp = -1j * ky * phi * np.conj(0.5 * tpar + tperp + 1.5 * dens)
    heat_flux = tmp + np.conj(tmp)
    return heat_flux


def get_plot_variable(mode, var, extend):
    """Returns plot variable and zgrid formatted for extended balloning structure, or not"""
    if extend:
        if var.shape[-1] == mode.nx:
            pvar = (var[:, mode.kx_modes] * mode.phase).ravel(order="F")
        else:
            pvar = (var * mode.phase).ravel(order="F")
        norm = pvar[mode.zero_ind]
        zgrid = mode.zgrid_ext
    else:
        pvar = var[:, 0]
        mid = mode.nz // 2
        norm = pvar[mid]
        zgrid = mode.zgrid
    if norm == 0:
        norm = 1
    pvar *= 1 / norm
    return pvar, zgrid
