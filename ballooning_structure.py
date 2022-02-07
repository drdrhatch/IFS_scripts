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
import time

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
parser.add_argument("--step", action="store", type=int, default=1, help="time step")
parser.add_argument("--plot", "-p", action="store_true", help="plot individual modes")
parser.add_argument("--output", "-o", action="store_true", help="output plots to pdfs")
parser.add_argument("--noshow", "-n", action="store_false", help="suppress popup plots")
parser.add_argument("--heat", "-q", action="store_true", help="calculate heat flux")
parser.add_argument("--debug", "-d", action="store_true", help="debug switch")
parser.add_argument("--notimeavg", "-N", action="store_false", help="Time average")
parser.add_argument(
    "--avgs",
    "-a",
    action="store",
    type=int,
    metavar="METHOD",
    help="compute averages of kz and omega",
)
parser.add_argument(
    "--corr",
    "-c",
    action="store_true",
    help="compute parallel correlation length and correlation time",
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
parser.add_argument(
    "--kx",
    dest="kx_cent",
    action="store",
    type=int,
    default=0,
    help="kx center mode",
)
args = parser.parse_args()

suffix = bl.check_suffix(args.suffix)

save_figs = args.output
show_figs = args.noshow
time_avg = args.notimeavg

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
times = times[:: args.step]
print("Analyzing for times: ", times)

kx_cent = args.kx_cent
print("kx_cent = ", kx_cent)

if args.ky_list == 0:
    ky_list = list(range(0, field.ny))
else:
    ky_list = args.ky_list
print("ky modes to analyze: ", ky_list)

if args.pod:
    pods = np.arange(1, min(args.pod, times.size) + 1)
else:
    pods = None

gene_files = {"pars": pars, "field": field, "mom": mom_e, "geometry": geometry}

print("Loading data for fields" + str(fields) + "...", end="")
start = time.time()
ky_modes = [bl.KyMode(ky, kx_cent, times, fields, gene_files) for ky in ky_list]
print("Done: ", str("{:6.3f}").format(time.time() - start), "s")

if args.debug:
    for mode in ky_modes:
        bl.plot_vars(mode, fields, times, show=show_figs, save=save_figs)

if not np.any(pods):
    spec = np.empty((len(ky_list), times.size))
    if args.avgs and not args.corr:
        scales = np.empty((len(ky_list), 2))
    if not args.avgs and args.corr:
        scales = np.empty((len(ky_list), 2))
    if args.avgs and args.corr:
        scales = np.empty((len(ky_list), 4))

print("Working on :")
for i, mode in enumerate(ky_modes):
    ky = mode.ky
    kx = mode.kx_cent
    print("ky = ", ky, "...", end="")
    print("kx modes = ", mode.kx_modes)
    start = time.time()
    if np.any(pods):
        u, sv, VH = bl.collective_pod(mode, fields, extend=True)
        bl.plot_singular_values(mode, sv, show_figs, save_figs)
        if args.debug:
            bl.pod_orthog_test(mode, u, VH)
            bl.test_pod(mode, u, sv, VH, fields)
        if save_figs:
            bl.output_pod(mode, u, sv, VH, fields, pods, times)
        if args.heat:
            weights = sv ** 2 / times.size
            Q_pod = bl.calc_heat_flux(mode, VH, weights)
            bl.plot_heat_flux(mode, Q_pod, show_figs, save_figs)
        if args.plot:
            bl.plot_time_dependence(mode, u, times, pods)
            if args.heat:
                bl.plot_pod(mode, Q_pod, pods, "q", extend=False)
            for var in fields:
                bl.plot_pod(mode, VH[var], pods, var)
        if args.corr:
            t, t_corr, corr_time = bl.autocorrelate(mode, u, times, axis=0)
            r, r_corr, corr_len = bl.autocorrelate(
                mode,
                VH["phi"],
                mode.zgrid_ext * np.pi,
                weights=mode.geometry["gjacobian"],
                axis=-1,
            )
            bl.output_scales(mode, corr_time, "corr_time")
            bl.output_scales(mode, corr_len, "corr_len")
        if args.avgs:
            if args.avgs == 1:
                avg_freq = bl.avg_freq(times, u)
                avg_kz = bl.avg_kz(mode, VH["phi"])
            elif args.avgs == 2:
                avg_freq, spec, omegas = bl.avg_freq2(
                    times, u, samplerate=2, spec_out=True
                )
                # bl.pod_kz_test(mode, u, sv, VH)
                avg_kz = bl.avg_kz2(mode, VH["phi"])
                avg_kz = bl.avg_kz2_pod(mode, VH["phi"])
                if args.plot:
                    bl.freq_spec_pod_plot(mode, omegas, spec, pods, output=True)
                varname = (
                    "pod_ky" + str(int(ky)).zfill(3) + "_kx" + str(int(kx)).zfill(3)
                )
                bl.output_spec_all_pod(pods, omegas, np.abs(spec), varname)
            bl.output_scales(mode, avg_freq, "avg_freq")
            bl.output_scales(mode, avg_kz, "avg_kz")
            end = time.time()
    else:
        phi = mode.fields["phi"]
        scale_list = []
        if args.avgs:
            if args.avgs == 1:
                if time_avg:
                    avg_freq = bl.avg_freq_tz(mode, times, phi)
                    avg_kz = bl.avg_kz_tz(mode, phi)
                else:
                    avg_freq = bl.avg_freq(times, phi)
                    avg_kz = bl.avg_kz(mode, phi)
            elif args.avgs == 2:
                if time_avg:
                    avg_freq = bl.avg_freq2_tz(mode, times, phi)
                    avg_kz = bl.avg_kz2_tz(mode, phi)
                else:
                    avg_freq = bl.avg_freq2(times, phi)
                    avg_kz = bl.avg_kz2(mode, phi)
            print("avg_kz = ", avg_kz)
            scale_list.append(avg_freq)
            scale_list.append(avg_kz)
        if args.corr:
            dphi = phi[:, :, 0]  # average over kx
            w1 = np.expand_dims(mode.geometry["gjacobian"], 0)
            w2 = mode.geometry["gjacobian"]
            doms, corr = bl.autocorrelate_tz(dphi, (times, mode.zgrid), w1)
            corr_time = bl.corr_len(doms[0], corr, axis=0)
            corr_len = bl.corr_len(doms[1], corr, 1, w2)
            if args.debug:
                bl.test_corr(mode, doms, corr)
            scale_list.append(corr_time)
            scale_list.append(corr_len)
        if args.avgs or args.corr:
            if time_avg:
                scales[i] = np.array(scale_list)
            else:
                scales = np.array(scale_list)
                bl.output_scales(ky_modes, scales, "phi" + suffix, "ev")
        omegas, spec[i] = bl.freq_spec(mode, times, phi, "phi", output=False)
    print(str("{:6.3f}").format(time.time() - start), "s")

if args.avgs and not np.any(pods):
    if time_avg:
        bl.output_scales(ky_modes, scales, "avgs", "avgs")
        varname = "phi2_kx" + str(int(kx_cent)).zfill(3)
        bl.output_spec_all_ky(ky_list, omegas, spec, varname)
