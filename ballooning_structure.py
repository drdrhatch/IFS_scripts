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
import fft

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

suffix = glib.check_suffix(args.suffix)
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
times = np.intersect1d(ftimes, mtimes)
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
fields = ("phi", "tpar", "tperp", "dens")

ky_modes = [bl.KyMode(ky, times, fields, gene_files) for ky in ky_list]

if args.debug:
    for mode in ky_modes:
        bl.plot_vars(mode, fields, times, show=show_figs, save=save_figs)

if pods:
    for mode in ky_modes:
        ky = mode.ky
        u, sv, VH = bl.collective_pod(mode, fields)
        bl.plot_singular_values(mode, sv, show_figs, save_figs)
        if save_figs:
            bl.output_pod(mode, u, sv, VH, fields, pods, times)
        if args.heat:
            Q = bl.calc_heat_flux(ky, VH)
            bl.plot_heat_flux(mode, Q, show_figs, save_figs)
        if args.plot:
            if args.heat:
                bl.plot_pod(mode, Q, pods, "q", extend=False)
            for var in fields:
                bl.plot_pod(mode, VH[var], pods, var)
        # bl.kz_pod_modes(mode, VH, pods, "phi", pars)
        dom_omega = bl.dom_freq(times,u)
