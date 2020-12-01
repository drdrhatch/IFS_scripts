#!/usr/bin/env python
# -*- coding: utf-8 -*-

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

try:
    import scipy.linalg as la
except ImportError:
    import numpy.linalg as la
import numpy as np
import ParIO as pario
import fieldlib
import momlib
import read_write_geometry as rwg

VARNAMES = {
    "phi": r"$\Phi$",
    "apar": r"$A_\parallel$",
    "bpar": r"$B_\parallel$",
    "tperp": r"$T_\perp$",
    "tpar": r"$T_\parallel$",
    "dens": "$n$",
    "q": "$Q$",
}

HEADER_NAMES = {
    "sv": "Singular values",
    "q": "Heat flux",
}


class KyMode:
    """Class for organizing ballooning structure for each ky mode"""

    def __init__(self, ky, times, fields, gene_files):
        pars = gene_files["pars"]
        field_file = gene_files["field"]
        mom_file = gene_files["mom"]
        geom_file = gene_files["geometry"]
        self.iky = ky
        self.ky = ky * pars["kymin"]
        self.nx = field_file.nx
        self.nz = field_file.nz
        self.construct_ranges(pars)
        self.define_phase(pars)
        self.define_dictionary(field_file, mom_file)
        self.geometry = geom_file
        self.read_fields(times, fields, field_file, mom_file)

    def construct_ranges(self, pars):
        self.kxrange(pars)
        self.zrange()

    def kxrange(self, pars):
        if self.ky == 0:
            step = 1
        else:
            step = pars["nexc"] * self.iky
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

    def define_dictionary(self, field_file, mom_file=None):
        self.field_vars = {
            "phi": field_file.phi,
            "apar": field_file.apar,
            "bpar": field_file.bpar,
        }
        if mom_file:
            self.field_vars.update(
                {
                    "dens": mom_file.dens,
                    "tpar": mom_file.tpar,
                    "tperp": mom_file.tperp,
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

    def read_fields(self, times, fields, field_file, mom_file):
        """Read given fields data for the given times"""
        self.fields_read = set(fields)
        tmp = np.empty((len(fields), times.size, self.nz, self.nx), dtype=np.cdouble)
        for j, time in enumerate(times):
            field_file.set_time(time)
            if mom_file:
                mom_file.set_time(time)
            for i, var in enumerate(fields):
                tmp[i, j, :, :] = self.read_field(var)
        for i, var in enumerate(fields):
            self.fields[var] = tmp[i]


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
    output_cum_sum(mode, sv, "sv")
    output_pod_modes(mode, vh, fields, pods, norm=True)
    output_time_modes(mode, u, pods, times)


def output_cum_sum(mode, var, varname):
    """Output variable and its cumulative sum"""
    filename = "./" + varname + "_ky" + str("{:03d}").format(int(mode.ky)) + ".dat"
    header = HEADER_NAMES[varname]
    var_sum = np.cumsum(var) / var.sum()
    data = np.vstack((var, var_sum)).T
    np.savetxt(filename, data, fmt="%g", header=header, encoding="UTF-8")


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
    data = np.hstack((times.reshape(-1, 1), np.abs(l_vec[:, : pods.stop])))
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


def plot_cumulative_array(mode, var, varname, show=True, fname=None):
    pods = np.arange(1, var.size + 1)

    fig, ax1 = plt.subplots()

    color = "red"
    ax1.set_ylabel("value", color=color)
    ax1.tick_params(axis="y", labelcolor=color)
    # ax1.plot(pods, var, marker="o", color=color)
    ax1.scatter(pods, var, marker="o", c=color)
    ax1.set_xlim(1, pods[-1])
    ax1.set_xlabel("POD #")
    ax1.set_xticks(np.arange(5, pods[-1] + 1, 5))
    ax1.set_xticks(pods, minor=True)

    ax2 = ax1.twinx()

    var_sum = np.cumsum(var) / var.sum()
    color = "blue"
    ax2.plot(pods, var_sum, color=color)
    # ax2.set_xlim(1, pods.stop)
    ax2.set_ylim(0, 1.0)
    ax2.set_ylabel("cumulative", color=color)
    ax2.tick_params(axis="y", labelcolor=color)
    ax2.grid()

    plt.title(varname + r" for mode $k_y = $" + str(mode.ky))
    plt.grid(True)
    if show:
        plt.show()
    if fname:
        pdf_figs = PdfPages("mode_" + str(int(mode.ky)) + "_" + fname + ".pdf")
        output = pdf_figs
        output.savefig(fig)
        pdf_figs.close()
    plt.close()


def plot_singular_values(mode, sv, show=True, save=False):
    if save:
        fname = "sv"
    else:
        fname = None
    plot_cumulative_array(mode, sv, "Singular values", show, fname)


def plot_heat_flux(mode, Q, show=True, save=False):
    heat = np.real(Q.sum(axis=(1, 2)))
    if save:
        fname = "qsum"
        output_cum_sum(mode, heat, "q")
    else:
        fname = None
    plot_cumulative_array(mode, heat, "Heat flux", show, fname)


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
        pvar = var.sum(axis=1)
        mid = mode.nz // 2
        norm = pvar[mid]
        zgrid = mode.zgrid
    if norm == 0:
        norm = 1
    pvar *= 1 / norm
    return pvar, zgrid


def get_input_params(directory, suffix, geom=None):
    par = pario.Parameters()
    par.Read_Pars(directory + "/parameters" + suffix)
    pars = par.pardict

    field = fieldlib.fieldfile(directory + "/field" + suffix, pars)
    mom_e = momlib.momfile(directory + "/mom_e" + suffix, pars)
    if geom:
        parameters, geometry = rwg.read_geometry_local(geom)
    else:
        geometry = None

    # min_time, max_time = field.get_minmaxtime()
    # stime = max(args.stime, min_time)
    # etime = min(args.etime, max_time)

    # ftimes = bl.get_times(field, stime, etime)
    # mtimes = bl.get_times(mom_e, stime, etime)
    # times = np.intersect1d(ftimes, mtimes)
    times = field.tfld
    gene_files = {"pars": pars, "field": field, "mom": mom_e, "geometry": geometry}
    return times, gene_files


def fft_freq(times, f, samplerate=2, axis=0):
    """Calculates fft of nonuniform data by first interpolating to uniform grid"""
    ntimes = times.size
    samples = samplerate * ntimes
    times_lin = np.linspace(times[0], times[-1], samples)
    if axis == 0:
        f_lin = np.empty((samples, f.shape[1]), dtype=np.cdouble)
        for i, row in enumerate(f.T):
            f_int = np.interp(times_lin, times, row)
            f_lin[:, i] = f_int.T
    else:
        f_lin = np.empty((f.shape[0], samples), dtype=np.cdouble)
        for i, row in enumerate(f):
            f_lin[i] = np.interp(times_lin, times, row)
    f_hat = np.fft.fft(f_lin, axis=axis)
    return f_hat, times_lin


def dom_freq(times, f, samplerate=2):
    """Returns the dominant frequency from field"""
    ntimes = times.size
    samples = samplerate * ntimes
    timestep = (times[-1] - times[0]) / samples
    omegas = np.fft.fftfreq(samples, d=timestep)
    f_hat, times_lin = fft_freq(times, f)
    dom_omega = omegas[
        np.argmax(abs(f_hat[1:]), axis=0) + 1
    ]  # skip freq zero as dominant
    return dom_omega
