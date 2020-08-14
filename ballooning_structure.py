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
parser.add_argument("--stime", "-s", action="store", help="start time window")
parser.add_argument("--etime", "-e", action="store", help="end time window")
parser.add_argument("--plot", "-p", action="store_true", help="plot modes")
parser.add_argument("--pod", "-P", action="store_true", help="perform POD analysis")
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
stime = args.stime

if args.stime == None:
    stime = min_time
else:
    stime = float(args.stime)

if args.etime == None:
    etime = max_time
else:
    etime = float(args.etime)

ky_list = args.ky_list
ky_modes = [bl.ky_mode(ky, field, pars) for ky in ky_list]

for time in bl.get_times(field, stime, etime):
    field.set_time(time)
    for mode in ky_modes:
        mode.read_phi()

if args.plot:
    for mode in ky_modes:
        zgrid = mode.zgrid / np.pi
        times = bl.get_times(field, stime, etime)
        for phi, time in zip(mode.phi, bl.get_times(field, stime, etime)):
            plt.title(
                r"$\phi$, $k_y=$" + str(mode.ky) + " t = " + str("{:6.3f}").format(time)
            )
            plt.plot(zgrid, np.real(phi), color="red", label=r"$Re[\phi]$")
            plt.plot(zgrid, np.imag(phi), color="blue", label=r"$Im[\phi]$")
            plt.plot(zgrid, np.abs(phi), color="black", label=r"$|\phi|$")
            ax = plt.axis()
            plt.legend()
            plt.xlabel(r"$z/\pi$", size=18)
            plt.show()
