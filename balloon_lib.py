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
import finite_differences as fd
import re

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
        self.T0 = pars["temp1"]
        self.n0 = pars["dens1"]
        self.construct_ranges(pars)
        self.define_phase(pars)
        self.define_dictionary(field_file, mom_file)
        self.geometry = geom_file
        self.read_fields(times, fields, field_file, mom_file, pars)

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
            fields = ("phi", "apar", "bpar", "dens", "tpar", "tperp", "q")
        else:
            fields = ("phi", "apar", "bpar")
        self.fields = dict.fromkeys(fields, None)

    def read_field(self, varname):
        """Read field for a given time window, returning array"""
        var = self.field_vars[varname]()
        if var.shape[1] == 1:  # for linear scan data with single ky
            indy = 0
        else:
            indy = self.iky
        tmp = var[:, indy, :]
        return tmp

    def read_fields(self, times, fields, field_file, mom_file, pars):
        """Read given fields data for the given times"""
        self.fields_read = set(fields)
        if pars["PRECISION"] == "DOUBLE":
            tmp = np.empty(
                (len(fields), times.size, self.nz, self.nx), dtype=np.cdouble
            )
        else:
            tmp = np.empty(
                (len(fields), times.size, self.nz, self.nx), dtype=np.csingle
            )
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


def plot_time_dependence(mode, u, times, pods):
    plt.title(r"Time dependece of POD modes")
    plt.xlabel("Time")
    plt.ylabel(r"$|\Phi_s|$")
    for ipod in pods:
        plt.plot(times, np.abs(u[:, ipod]), label=r"$s_" + str(ipod + 1) + "$")
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
    Q_x = np.sum(Q, axis=2)
    Q_xz = np.average(Q_x, weights=mode.geometry["gjacobian"], axis=1)
    if save:
        fname = "qsum"
        output_cum_sum(mode, Q_xz, "q")
    else:
        fname = None
    plot_cumulative_array(mode, Q_xz, "Heat flux", show, fname)


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


def pod(mode, var):
    ntimes = var.shape[0]
    pvar = var.reshape(ntimes, -1, order="F")
    u, sv, vtmp = la.svd(pvar, full_matrices=False)
    vh = vtmp.reshape(-1, mode.nz, mode.nx, order="F")
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


def calc_heat_flux(mode, fields, weights=None):
    phi = fields["phi"]
    tpar = fields["tpar"]
    tperp = fields["tperp"]
    dens = fields["dens"]
    ky = mode.ky
    n0 = mode.n0
    T0 = mode.T0
    if "C_xy" in mode.geometry:
        Cxy = mode.geometry["C_xy"]
    else:
        Cxy = 1
    temp1 = -1j * n0 * T0 * ky * phi / Cxy * np.conj(0.5 * tpar + tperp + 1.5 * dens)
    # \/ not divided by 2 because we only have half the ky modes
    temp2 = np.real_if_close(temp1 + np.conj(temp1))
    if np.any(weights):
        heat_flux = weights[:, np.newaxis, np.newaxis] * temp2
    else:
        heat_flux = temp2
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
    pvar_norm = pvar / norm
    return pvar_norm, zgrid


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


def fft_nonuniform(times, f, axis=0, samplerate=2):
    """Calculates fft of nonuniform data by first interpolating to uniform grid"""
    datatype = f.dtype
    ntimes = times.size
    samples = samplerate * ntimes
    times_lin = np.linspace(times[0], times[-1], samples)
    if f.ndim > 1:
        if axis == 0:
            f_lin = np.empty((samples, f.shape[1]), dtype=datatype)
            for i, row in enumerate(f.T):
                f_int = np.interp(times_lin, times, row)
                f_lin[:, i] = f_int.T
        else:
            f_lin = np.empty((f.shape[0], samples), dtype=datatype)
            for i, row in enumerate(f):
                f_lin[i] = np.interp(times_lin, times, row)
    else:
        f_lin = np.interp(times_lin, times, f)
        axis = 0
    f_hat = np.fft.fft(f_lin, axis=axis)
    return f_hat, times_lin


def avg_freq(times, f, axis=0, samplerate=2):
    """Returns the dominant frequency from field"""
    ntimes = times.size
    dt = np.diff(times)
    even_dt = np.all(dt == dt[0])
    if not even_dt:
        samples = samplerate * ntimes
        f_hat, times_lin = fft_nonuniform(times, f)
    else:
        samples = ntimes
        f_hat = np.fft.fft(f, axis=axis)
    timestep = (times[-1] - times[0]) / samples
    omegas = np.fft.fftfreq(samples, d=timestep)
    if f.ndim > 1:
        if axis == 0:
            weights = np.sum(np.expand_dims(omegas, -1) ** 2 * abs(f_hat) ** 2, axis=0)
        elif axis == 1:
            weights = np.sum(np.expand_dims(omegas, 0) ** 2 * abs(f_hat) ** 2, axis=1)
    else:
        weights = np.sum(omegas ** 2 * abs(f_hat) ** 2)
    norm = np.sum(abs(f_hat) ** 2, axis=axis)
    freq = np.sqrt(weights / norm)
    return freq


def get_extended_var(mode, var):
    """Flattens array over last two dimensions to return z-extended variable"""
    phase = np.expand_dims(mode.phase, axis=0)
    newshape = (var.shape[0], -1)
    ext_var = np.reshape(var * phase, newshape, order="F")
    return ext_var


def avg_kz(mode, var, outspect=False):
    """Calculate the average kz mode weighted by given field"""
    jacxBpi = mode.geometry["gjacobian"] * mode.geometry["gBfield"] * np.pi
    jacxBpi_ext = np.expand_dims(np.tile(jacxBpi, mode.kx_modes.size), -1)
    if var.ndim > 2:
        var_ext = get_extended_var(mode, var)
    else:
        var_ext = var
    if var.ndim > 1:
        field = var_ext.T
    else:
        field = np.expand_dims(var, axis=-1)

    zgrid = mode.zgrid_ext
    dfielddz = fd.fd_d1_o4(field, zgrid) / jacxBpi_ext

    # Select range, cutting off extreme ends of z domain
    zstart, zend = 5, len(zgrid) - 5
    dz = np.expand_dims(np.diff(zgrid)[zstart:zend], -1)
    dfdz1 = dfielddz[zstart:zend]
    dfdz2 = dfielddz[zstart + 1 : zend + 1]
    jac = jacxBpi_ext[zstart:zend]
    f1 = field[zstart:zend]
    f2 = field[zstart + 1 : zend + 1]

    ddz = 0.5 * (abs(dfdz1) ** 2 + abs(dfdz2) ** 2) / dz * jac
    sum_ddz = np.sum(ddz, axis=0)
    denom = np.sum(0.5 * (abs(f1) ** 2 + abs(f2) ** 2) / dz * jac, axis=0)
    akz = np.sqrt(sum_ddz / denom).T
    if outspect:
        return akz, ddz
    return akz


def output_scales(modes, scales, varname, intype="POD"):
    """Output a list of scales for a mode, e.g. frequencies or correlation lengths"""
    if intype == "POD":
        ky = str("{:03d}").format(int(modes.ky)) + "_pod"
        header = "POD # " + varname
    else:
        ky = "_all"
        header = "ky avg_omega avg_kz corr_time corr_len"
    filename = "./" + varname + "_ky" + ky + ".dat"
    if scales.ndim == 1:
        pods = np.arange(1, scales.size + 1)
        data = np.vstack((pods, scales)).T
    else:
        kylist = np.expand_dims(np.array([mode.ky for mode in modes]).T, -1)
        data = np.hstack((kylist, scales))
    np.savetxt(
        filename,
        data,
        fmt="% E",
        header=header,
        encoding="UTF-8",
    )


def autocorrelate(mode, var, domain, axis=-1, samplerate=2, tol=1e-6):
    """Calculate correlation length/time for given input field"""
    datatype = var.dtype
    if var.ndim > 2:
        fvar = get_extended_var(mode, var)
    else:
        fvar = var

    dt = np.diff(domain)
    test_dt = np.floor(dt / tol)
    even_dt = np.all(test_dt == test_dt[0])
    if not even_dt:
        npts = domain.size
        samples = samplerate * npts
        dom_lin = np.linspace(domain[0], domain[-1], samples)
        if axis == 0:
            f_lin = np.empty((fvar.shape[1], samples), dtype=datatype)
            for i, row in enumerate(fvar.T):
                f_int = np.interp(dom_lin, domain, row).T
                f_lin[i] = f_int.T
        else:
            f_lin = np.empty((fvar.shape[0], samples), dtype=datatype)
            if fvar.ndim > 1:
                for i, row in enumerate(fvar):
                    f_lin[i] = np.interp(dom_lin, domain, row)
            else:
                f_lin = np.interp(dom_lin, domain, fvar)
        dom = dom_lin
        f = f_lin
    else:
        dom = domain
        if axis == 0:
            f = fvar.T
        else:
            f = fvar

    N = f.shape[-1]
    N2 = N // 2
    norm = N - np.arange(0, N2)
    if f.ndim > 1:
        corr = np.empty((f.shape[0], N2), dtype=datatype)
        for i, row in enumerate(f):
            f1 = row
            corr[i] = np.correlate(f1, f1, mode="same")[N2:] / norm
            corr[i] /= corr[i, 0]
    else:
        f1 = f
        corr = np.correlate(f1, f1, mode="same")[N2:] / norm
        corr /= corr[0]
    r = np.linspace(0, (dom[-1] - dom[0]) / 2, N2)
    scale = r[1] - r[0]
    corr_len = scale * np.real(np.sum(corr, axis=-1))
    return r, corr, corr_len


def norm_z_field(mode, var):
    fvar = var[:, :, mode.kx_modes]
    evar = get_extended_var(mode, fvar)
    norm_var = la.norm(evar, axis=-1)
    return norm_var


def avg_t_field(mode, var):
    fvar = var[:, :, mode.kx_modes]
    evar = get_extended_var(mode, fvar)
    avg_var = np.mean(evar, axis=0)
    return avg_var


def avg_kz_tz(mode, var):
    fvar = var[:, :, mode.kx_modes]
    evar = get_extended_var(mode, fvar)
    kz = avg_kz(mode, evar)
    mean_kz = np.sqrt(np.mean(kz ** 2))
    return mean_kz


def avg_freq_tz(mode, times, var):
    fvar = var[:, :, mode.kx_modes]
    evar = get_extended_var(mode, fvar)
    omega = avg_freq(times, evar)
    avg_omega = np.sqrt(np.mean(omega ** 2))
    return avg_omega


def freq_spec(mode, times, varname, axis=0, samplerate=2, output=False):
    var = mode.fields[varname]
    fvar = var[:, :, mode.kx_modes]
    evar = get_extended_var(mode, fvar)
    f = evar
    ntimes = times.size
    dt = np.diff(times)
    even_dt = np.all(dt == dt[0])
    if not even_dt:
        samples = samplerate * ntimes
        f_hat, times_lin = fft_nonuniform(times, f)
    else:
        samples = ntimes
        f_hat = np.fft.fft(f, axis=axis)
    timestep = (times[-1] - times[0]) / samples
    omegas = np.fft.fftfreq(samples, d=timestep)

    dims = tuple(range(f_hat.ndim))
    saxes = dims[dims != axis]
    spectrum = np.sum(np.abs(f_hat) ** 2, axis=saxes)

    om = np.array(omegas[0 : samples // 2])
    spec = np.union1d(
        spectrum[0], spectrum[1 : samples // 2] + spectrum[-1 : -samples // 2 : -1]
    )

    if output:
        output_spec(mode, om, spec, varname)
    return omegas, spec


def output_spec(mode, omegas, spec, varname):
    """Output a list of scales for a mode, e.g. frequencies or correlation lengths"""
    header = "omega " + varname + "^2"
    ky = str("{:03d}").format(int(mode.ky))
    filename = "./" + varname + "_ky" + ky + "_spec.dat"
    data = np.vstack((omegas, spec)).T
    np.savetxt(
        filename,
        data,
        fmt="% E",
        header=header,
        encoding="UTF-8",
    )


def check_suffix(run_number):
    if re.search("dat$", run_number):
        suffix = ".dat"
    elif re.search("[0-9]{1,4}$", run_number):
        match = re.search(r"([0-9]{1,4})$", run_number)
        suffix = "_" + match.group(0).zfill(2)
    else:
        print("Please enter a valid run number, e.g. .dat or 0123")
        return None
    return suffix
def test_pod(mode, u, sv, vh, fields):
    """testing that pod behaved in the expected way"""
    npods = u.shape[0]
    nx = mode.nx
    nz = mode.nz
    print("Shapes")
    print("------------")
    print("u : ", u.shape)
    print("sv : ", sv.shape)
    print("vh['phi'] : ", vh["phi"].shape)
    for field in fields:
        original = mode.fields[field]
        new = np.empty(original.shape, dtype=original.dtype)
        temp2 = np.zeros((npods, nx * nz), dtype=original.dtype)
        for i in range(npods):
            temp1 = vh[field].reshape((npods, -1))
            temp2 += sv[i] * np.outer(u[:, i], temp1[i])
        new = temp2.reshape((-1, nz, nx))
        print("Are the arrays close?....", np.allclose(original, new, atol=1e-6))
    Q = calc_heat_flux(mode, vh, sv ** 2)
    Q_sum = np.mean(
        np.average(Q.sum(axis=2), axis=1, weights=mode.geometry["gjacobian"]), axis=0
    )
    print("Q_sum(ky = %d) = ", mode.ky, Q_sum)
