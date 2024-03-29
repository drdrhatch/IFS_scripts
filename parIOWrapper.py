#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ParIO import *
import optparse as op


def init_read_parameters_file(suffix):
    #    print 'Reading', 'parameters'+suffix
    par = Parameters()
    par.Read_Pars("parameters" + suffix)
    pars = par.pardict

    return pars


def read_ref_values(suffix, pars):
    if "Bref" in pars:
        Bref = pars["Bref"]
    else:
        Bref = 1.0
        print("Bref not in parameters" + suffix + ". Bref = 1")

    if "Tref" in pars:
        Tref = pars["Tref"]
    else:
        Tref = 1.0
        print("Tref not in parameters" + suffix + ". Tref = 1")

    if "nref" in pars:
        nref = pars["nref"]
    else:
        nref = 1.0
        print("nref not in parameters" + suffix + ". nref = 1")

    if "Lref" in pars:
        Lref = pars["Lref"]
    else:
        Lref = 1.0
        print("Lref not in parameters" + suffix + ". Lref = 1")

    if "mref" in pars:
        mref = pars["mref"]
    else:
        mref = 1.0
        print("mref not in parameters" + suffix + ". mref = 1")

    return Bref, Tref, nref, Lref, mref


def read_species_gradients(q_charge, pars):
    if "x_local" in pars:
        if pars["x_local"]:
            x_local = True
        else:
            x_local = False
    else:
        x_local = True

    if x_local:
        if "omn1" in pars:
            if pars["charge1"] == q_charge:
                return pars["omn1"], pars["omt1"]
            elif "omn2" in pars:
                if pars["charge2"] == q_charge:
                    return pars["omn2"], pars["omt2"]
                elif "omn3" in pars:
                    if pars["charge3"] == q_charge:
                        return pars["omn3"], pars["omt3"]
                    else:
                        return 0, 0
                        print(
                            "No species with charge = " + str(q_charge) + " is found."
                            "omn = 0, omt = 0"
                        )
    else:
        return 1.0, 1.0


def read_species_tempdens(q_charge, pars):
    if "temp1" in pars:
        if pars["charge1"] == q_charge:
            return pars["temp1"], pars["dens1"]
        elif "temp2" in pars:
            if pars["charge2"] == q_charge:
                return pars["temp2"], pars["dens2"]
            elif "temp3" in pars:
                if pars["charge3"] == q_charge:
                    return pars["temp3"], pars["dens3"]
                else:
                    return 0, 0
                    print(
                        "No species with charge = " + str(q_charge) + " iis found."
                        "temp = 0, dens = 0"
                    )


def create_parameters_dict(parameters_filepath: str):
    # Create a parameter dictionary using the Parameters class
    par = Parameters()
    par.Read_Pars(parameters_filepath)  # Read the parameter file
    parameter_dict = par.pardict

    # Add the filepath and parameter keys to the parameter dictionary
    parameter_dict["filepath"] = parameters_filepath
    parameter_dict["key_list"] = list(parameter_dict.keys())

    return parameter_dict


def check_suffix(run_number):
    """Check that given suffix matches typical GENE output"""
    if re.search("dat$", run_number):
        suffix = ".dat"
    elif re.search("[0-9]{1,4}$", run_number):
        match = re.search(r"([0-9]{1,4})$", run_number)
        suffix = "_" + match.group(0).zfill(2)
    else:
        print("Please enter a valid run number, e.g. .dat or 0123")
        return None
    return suffix


def read_parameters_files(runlist):
    """Given a list of GENE runs as numbers,
    return a generator containing the parameters files for each run in list"""
    return (init_read_parameters_file(check_suffix(suffix)) for suffix in runlist)
