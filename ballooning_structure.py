#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import genelib as glib
import ParIO as pario
import fieldlib
import momlib
import matplotlib.pyplot as plt
import balloon_lib as bl
import read_write_geometry as rwg

parser = argparse.ArgumentParser()
parser.add_argument("suffix", help="run number or .dat suffix of output data")
parser.add_argument(
    "--geom", "-g", action="store", metavar="GEOMETRY_FILE", help="geometry file"
)
parser.add_argument(
    "--stime", "-s", action="store", type=float, default=0, help="start time window"
)
parser.add_argument(
    "--etime", "-e", action="store", type=float, default=999, help="end time window"
)
parser.add_argument("--plot", "-p", action="store_true", help="plot individual modes")
parser.add_argument("--output", "-o", action="store_true", help="output plots to pdfs")
parser.add_argument("--noshow", "-n", action="store_false", help="suppress popup plots")
parser.add_argument("--heat", "-q", action="store_true", help="calculate heat flux")
parser.add_argument("--debug", "-d", action="store_true", help="debug switch")
parser.add_argument(
    "--avgs", "-a", action="store_true", help="find avg omega and kz for each ky mode"
)
parser.add_argument(
    "--eigen",
    "-E",
    action="store",
    nargs="+",
    metavar="EIGENVALUE_DIRECTORY RUN_NUMBER_0001 [RUN_NUMBER_0002 ...]",
    help="eigenvalue directory, followed by a list of run numbers",
)
parser.add_argument(
    "--pod",
    "-P",
    action="store",
    type=int,
    metavar="N",
    help="POD analysis, plotting N modes",
)
parser.add_argument(
    "--kylist",
    "-k",
    dest="ky_list",
    action="store",
    default=0,
    nargs="+",
    type=int,
    help="list of ky modes",
)
args = parser.parse_args()

suffix = bl.check_suffix(args.suffix)
if args.eigen:
    esuffix = [glib.check_suffix(eigen) for eigen in args.eigen[1:]]

save_figs = args.output
show_figs = args.noshow

par = pario.Parameters()
par.Read_Pars("parameters" + suffix)
pars = par.pardict

field = fieldlib.fieldfile("field" + suffix, pars)
mom_e = momlib.momfile("mom_e" + suffix, pars)
parameters, geometry = rwg.read_geometry_local(args.geom)

min_time, max_time = field.get_minmaxtime()
stime = max(args.stime, min_time)
etime = min(args.etime, max_time)

ftimes = bl.get_times(field, stime, etime)
mtimes = bl.get_times(mom_e, stime, etime)
if args.heat:  # moment values needed for heat flux calc
    times = np.intersect1d(ftimes, mtimes)
    fields = ("phi", "tpar", "tperp", "dens")
else:  # otherwise, default to phi
    times = ftimes
    mom_e = None
    fields = ("phi",)
print("Analyzing for times: ", times)

if args.ky_list == 0:
    ky_list = [range(0, field.ny)]
else:
    ky_list = args.ky_list
print("ky modes to analyze: ", ky_list)

if args.pod:
    pods = range(min(args.pod, times.size))
else:
    pods = None

gene_files = {"pars": pars, "field": field, "mom": mom_e, "geometry": geometry}

ky_modes = [bl.KyMode(ky, times, fields, gene_files) for ky in ky_list]

if args.debug:
    for mode in ky_modes:
        bl.plot_vars(mode, fields, times, show=show_figs, save=save_figs)

if args.avgs and not pods:
    scales = np.empty((len(ky_list), 4))

for i, mode in enumerate(ky_modes):
    ky = mode.ky
    if pods:
        u, sv, VH = bl.collective_pod(mode, fields)
        bl.plot_singular_values(mode, sv, show_figs, save_figs)
        if save_figs:
            bl.output_pod(mode, u, sv, VH, fields, pods, times)
        if args.heat:
            Q = bl.calc_heat_flux(ky, VH)
            bl.plot_heat_flux(mode, Q, show_figs, save_figs)
        if args.plot:
            bl.plot_time_dependence(mode, u, times, pods)
            if args.heat:
                bl.plot_pod(mode, Q, pods, "q", extend=False)
            for var in fields:
                bl.plot_pod(mode, VH[var], pods, var)
        if args.avgs:
            t, t_corr, corr_time = bl.autocorrelate(mode, u, times, axis=0)
            r, r_corr, corr_len = bl.autocorrelate(
                mode, VH["phi"], mode.zgrid_ext * np.pi, axis=-1
            )
            avg_freq = bl.avg_freq(times, u)
            avg_kz = bl.avg_kz(mode, VH["phi"])
            bl.output_scales(mode, avg_freq, "avg_freq")
            bl.output_scales(mode, avg_kz, "avg_kz")
            bl.output_scales(mode, corr_time, "corr_time")
            bl.output_scales(mode, corr_len, "corr_len")
    else:
        if args.avgs:
            phi = mode.fields["phi"]
            norm_phi = bl.norm_z_field(mode, phi)
            avg_phi = bl.avg_t_field(mode, phi)
            t, t_corr, corr_time = bl.autocorrelate(mode, norm_phi, times, axis=-1)
            r, r_corr, corr_len = bl.autocorrelate(
                mode, avg_phi, mode.zgrid_ext * np.pi, axis=-1
            )
            avg_freq = bl.avg_freq_tz(mode, times, phi)
            avg_kz = bl.avg_kz_tz(mode, phi)
            scales[i] = [avg_freq, avg_kz, corr_time, corr_len]
            omegas, spec = bl.freq_spec(mode, times, "phi", output=True)

if args.avgs and not pods:
    bl.output_scales(ky_modes, scales, "avgs", "avgs")
