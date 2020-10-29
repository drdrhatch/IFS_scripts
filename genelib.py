#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re


def check_suffix(run_number):
    if re.search("dat$", run_number):
        suffix = ".dat"
    elif re.search("[0-9]{1,4}$", run_number):
        match = re.search(r"([0-9]{1,4})$", run_number)
        suffix = "_" + match.group(0).zfill(4)
    else:
        print("Please enter a valid run number, e.g. .dat or 0123")
        return None
    return suffix
