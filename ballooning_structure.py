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
    "species_list", nargs="+", metavar="SPECIES", help="list of species to load"
)
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
    "--avg",
    "-a",
    action="store_true",
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
    dest="iky_list",
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
mom_list = []
for species in args.species_list:
    momfile = "mom_" + species + suffix
    mom_list.append(momlib.momfile(momfile, pars))

parameters, geometry = rwg.read_geometry_local(args.geom)

min_time, max_time = field.get_minmaxtime()
stime = max(args.stime, min_time)
etime = min(args.etime, max_time)

ftimes = bl.get_times(field, stime, etime)
mtimes = bl.get_times(mom_list[0], stime, etime)  # momtimes from first species
if args.heat:  # moment values needed for heat flux calc
    times = np.intersect1d(ftimes, mtimes)
    fields = ("phi", "tpar", "tperp", "dens")
else:  # otherwise, default to phi
    times = ftimes
    mom_list = None
    fields = ("phi",)
times = times[:: args.step]
print("Analyzing for times: ", times)

kx_cent = args.kx_cent
print("kx_cent = ", kx_cent)

if args.iky_list == 0:
    iky_list = list(range(0, field.ny))
else:
    iky_list = args.iky_list
print("ky modes to analyze: ", iky_list)

if args.pod:
    pods = np.arange(1, min(args.pod, times.size) + 1)
else:
    pods = None

gene_files = {"pars": pars, "field": field, "mom_list": mom_list, "geometry": geometry}

# print("Loading data for fields" + str(fields) + "...", end="")
# ky_modes = [bl.KyMode(ky, kx_cent, times, fields, gene_files) for ky in iky_list]
# print("Done: ", str("{:6.3f}").format(time.time() - start), "s")

# if args.debug:
#     for mode in ky_modes:
#         bl.plot_vars(mode, fields, times, show=show_figs, save=save_figs)

if not np.any(pods):
    spec = np.empty((len(iky_list), times.size))

scale_dict = {}
ky_list = []

for i, iky in enumerate(iky_list):
    print(
        "Loading data for ky=" + str(iky) + " and fields" + str(fields) + "...", end=""
    )
    start = time.time()
    mode = bl.KyMode(iky, kx_cent, times, fields, gene_files)
    print("Done: ", str("{:6.3f}").format(time.time() - start), "s")
    ky = mode.ky
    kx = mode.kx_cent
    ky_list.append(ky)
    print("Working on :")
    print("ky = ", ky, "...", end="")
    print("connected kx modes = ", mode.kx_modes)
    start = time.time()
    if np.any(pods):
        ltimes, ldata = bl.resample_time(mode, fields, times)
        u, sv, VH = bl.collective_pod(mode, ldata, fields, extend=True)
        bl.plot_singular_values(mode, sv, show_figs, save_figs)
        if args.debug:
            bl.pod_orthog_test(mode, u, VH)
            bl.test_pod(mode, u, sv, VH, fields)
        if save_figs:
            bl.output_pod(mode, u, sv, VH, fields, pods, times)
        if args.heat:
            weights = sv**2 / times.size
            Q_pod = bl.calc_heat_flux(mode, VH, weights)
            bl.plot_heat_flux(mode, Q_pod, show_figs, save_figs)
        if args.plot:
            bl.plot_time_dependence(mode, u, times, pods)
            if args.heat:
                bl.plot_pod(mode, Q_pod, pods, "q", extend=False)
            for var in fields:
                bl.plot_pod(mode, VH[var], pods, var)
        if args.corr:
            t, t_corr, corr_time = bl.autocorrelate(mode, u, ltimes, axis=0)
            r, r_corr, corr_len = bl.autocorrelate(
                mode,
                VH["phi"],
                mode.zgrid_ext * np.pi,
                weights=mode.geometry["gjacobian"],
                axis=-1,
            )
            scale_dict["corr_len"] = corr_len
            scale_dict["corr_time"] = corr_time
        if args.avg:
            avg_freq = bl.avg_freq(ltimes, u)
            avg_kz = bl.avg_kz_pod(mode, VH["phi"], sv)
            avg_freq2, spec, omegas = bl.avg_freq2(
                ltimes, u, samplerate=2, spec_out=True
            )
            wspec = sv**2 / np.sum(sv**2) * spec
            avg_kz2 = bl.avg_kz2_pod(mode, VH["phi"], sv)
            if save_figs:
                bl.freq_spec_pod_plot(mode, omegas, wspec, pods, output=True)
            varname = "pod_ky" + str(int(ky)).zfill(3) + "_kx" + str(int(kx)).zfill(3)
            bl.output_spec_all_pod(pods, omegas, wspec, varname)
            scale_dict["avg_freq"] = avg_freq
            scale_dict["avg_freq_rms"] = avg_freq2
            scale_dict["avg_kz"] = avg_kz
            scale_dict["avg_kz_rms"] = avg_kz2
        bl.output_scales(ky, kx_cent, scale_dict, "phi", "POD")
        end = time.time()
    else:
        phi = mode.fields["phi"]
        w1 = np.expand_dims(mode.geometry["gjacobian"], (0, 2))
        w2 = mode.geometry["gjacobian"]
        if args.corr:
            # doms, corr = bl.autocorrelate_tz(phi, (times, mode.zgrid), w1)
            doms, corr = bl.autocorrelate_tz(phi, (times, mode.zgrid))
            corr_time = bl.corr_len(doms[0], corr, axis=0)
            corr_len = bl.corr_len(doms[1], corr, 1, w2)
            if args.debug:
                bl.test_corr(mode, doms, corr)
            if i == 0:
                scale_dict.update({"corr_len": [], "corr_time": []})
            scale_dict["corr_len"].append(corr_len)
            scale_dict["corr_time"].append(corr_time)
        if args.avg:
            if time_avg:
                avg_freq = bl.avg_freq_tz(mode, times, phi)
                avg_kz = bl.avg_kz_tz(mode, phi)
                avg_freq2 = bl.avg_freq2_tz(mode, times, phi)
                avg_kz2 = bl.avg_kz2_tz(mode, phi)
            else:
                avg_freq = bl.avg_freq(times, phi)
                avg_kz = bl.avg_kz(mode, phi)
                avg_freq2 = bl.avg_freq2(times, phi)
                avg_kz2 = bl.avg_kz2(mode, phi)
            if i == 0:
                scale_dict.update(
                    {"avg_freq": [], "avg_freq_rms": [], "avg_kz": [], "avg_kz_rms": []}
                )
            scale_dict["avg_freq"].append(avg_freq)
            scale_dict["avg_freq_rms"].append(avg_freq2)
            scale_dict["avg_kz"].append(avg_kz)
            scale_dict["avg_kz_rms"].append(avg_kz2)
        omegas, spec[i] = bl.freq_spec(
            mode, times, phi, "phi", weights=w2, output=False
        )
    print(str("{:6.3f}").format(time.time() - start), "s")

if args.avg and not np.any(pods):
    if time_avg:
        bl.output_scales(ky_list, kx_cent, scale_dict, "avgs_phi", "avgs")
if not np.any(pods):
    varname = "phi2_kx" + str(int(kx_cent)).zfill(3)
    bl.output_spec_all_ky(ky_list, omegas, spec, varname)
