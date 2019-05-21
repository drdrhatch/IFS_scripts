#!/usr/bin/env python
# -*- coding: utf-8 -*-

from omega_tool import omega_calc
import os

cwd = os.getcwd()
filelist = []
for filename in os.listdir(cwd):
    if filename.startswith("omega"):
        filelist.append(filename[-5:])
filelist.sort()

for suffix in filelist:
    omega_calc(suffix)
