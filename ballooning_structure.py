#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import re
import ParIO as pario
import fieldlib
import matplotlib.pyplot as plt
import balloon_lib as bl

parser = argparse.ArgumentParser()
parser.add_argument("suffix", help="run number or .dat suffix of output data")
parser.add_argument(
    "--stime", "-s", action="store", type=float, default=0, help="start time window"
)
parser.add_argument(
    "--etime", "-e", action="store", type=float, default=999, help="end time window"
)
parser.add_argument("--plot", "-p", action="store_true", help="plot individual modes")
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

par = pario.Parameters()
par.Read_Pars("parameters" + suffix)
pars = par.pardict

field = fieldlib.fieldfile("field" + suffix, pars)

min_time, max_time = field.get_minmaxtime()
stime = max(args.stime, min_time)
etime = min(args.etime, max_time)

ky_list = args.ky_list
ky_modes = [bl.ky_mode(ky, field, pars) for ky in ky_list]

for time in bl.get_times(field, stime, etime):
    field.set_time(time)
    for mode in ky_modes:
        mode.read_phi()

if args.plot:
    for mode in ky_modes:
        mode.plot_modes(bl.get_times(field, stime, etime))

if args.pod:
    for mode in ky_modes:
        mode.pod()
        print("singular values = ", mode.s)
        mode.plot_singular_values()
        mode.plot_time_dependence(bl.get_times(field, stime, etime), range(args.pod))
        mode.plot_pod(range(args.pod))
