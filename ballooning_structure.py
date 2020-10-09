#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import re
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

if re.search("dat$", args.suffix):
    suffix = ".dat"
elif re.search("[0-9]{1,4}$", args.suffix):
    match = re.search(r"([0-9]{1,4})$", args.suffix)
    suffix = "_" + match.group(0).zfill(4)
    print(suffix)
else:
    print("Please enter a valid run number, e.g. .dat or 0123")

save_figs = args.output
show_figs = args.noshow

par = pario.Parameters()
par.Read_Pars("parameters" + suffix)
pars = par.pardict

field = fieldlib.fieldfile("field" + suffix, pars)
mom_e = momlib.momfile("mom_e" + suffix, pars)
geometry = rwg.read_geometry_local(args.geom)

min_time, max_time = field.get_minmaxtime()
stime = max(args.stime, min_time)
etime = min(args.etime, max_time)

ftimes = bl.get_times(field, stime, etime)
mtimes = bl.get_times(mom_e, stime, etime)
times = np.intersect1d(ftimes, mtimes)
print("Analyzing for times: ", times)

if args.ky_list == 0:
    ky_list = [ky for ky in range(0, field.ny)]
else:
    ky_list = args.ky_list
print("ky modes to analyze: ", ky_list)

if args.pod:
    pods = range(args.pod)
else:
    pods = None

fields = ("phi", "tpar", "tperp", "dens")
ky_modes = [
    bl.KyMode(ky, pars, times, fields, field, mom_e, geometry) for ky in ky_list
]

for mode in ky_modes:
    if args.plot:
        bl.plot_vars(mode, fields, times, show=show_figs, save=save_figs)

# sumq = bl.sum_modes(ky_modes, "q")
# zgrid = ky_modes[-1].zgrid
# if args.plot:
#     for i, time in enumerate(times):
#         title = r"$k_y=" + str(ky_list) + ",t = " + str("{:6.3f}").format(time) + "$"
#         varname = r"$\sum_k Q_k$"
#         bl.plot(zgrid, sumq[i], varname, title)
#         plt.show()
#     pass

if pods:
    for mode in ky_modes:
        ky = mode.ky
        up, svp, VHp = bl.pod(mode, "phi")
        u, sv, VH = bl.collective_pod(mode, fields)
        bl.plot_singular_values(mode, svp, save_figs)
        bl.plot_singular_values(mode, sv, save_figs)
        Q = bl.calc_heat_flux(ky, VH)
        bl.plot_pod(mode, VHp, pods, "phi")
        bl.plot_pod(mode, Q, pods, "q")
        for var in fields:
            bl.plot_pod(mode, VH[var], pods, var)
