# -*- coding: utf-8 -*-
"""
Module handling nrg files

"""
# pylint: disable=E1101
import numpy as np
import sys
import csv
import errors as err
from baseplot import plt, Plotting
from bisect import bisect_left
from ParIO import Parameters


class Nrgdata(object):
    """ Nrgdata: Class to read a nrg file from GENE """

    def __init__(self, fext, starttime, endtime):
        self.isdatapresent = False
        self.fext = fext
        self.filename = 'nrg{}'.format(self.fext)
        self.starttime = starttime
        self.endtime = endtime
        self.timefld = []
        self.isdatapresent = False
        self.n_col = 10  # Number of columns in nrg files
        # Names of the columns in LaTeX readable form
        self.colnames = [r"|n_1|^2", r"|u_{1,\parallel}|^2", r"|T_{1,\parallel}|^2",
                         r"|T_{1,\perp}|^2", r"\Gamma_{es}", r"\Gamma_{em}", r"Q_{es}", r"Q_{em}",
                         r"\Pi_{es}", r"\Pi_{em}"]
        self.units = {'nref': 1.0, 'Tref': 1.0, 'Bref': 1.0, 'Lref': 1.0,
                      'mref': 1.0}
        self.nrgcols = []
        self.pnt = 0
        self.specname = [[], []]
        self.readpar()

    def readpar(self):
        """ Get some essential parameters from par file"""
        par = Parameters()
        par.Read_Pars('parameters{}'.format(self.fext))
        self.pnt = par.asnamedtuple()
        for nsp in range(self.pnt.n_spec):
            self.specname[nsp] = str(par.pardict['name{:1d}'.format(nsp+1)])
            self.specname[nsp] = self.specname[nsp][1:-1]
        for k in self.units:
            if k in par.pardict:
                self.units[k] = par.pardict[k]

    def get_minmaxtime(self):
        return self.timefld[0], self.timefld[-1]

    def readnrg(self):
        """ Fill the Nrgdata object with data """
        n_col = self.n_col

        try:
            with open(self.filename) as nrgfile:
                csvnrg = csv.reader(nrgfile, delimiter=' ', skipinitialspace=True)
                for line in csvnrg:
                    if len(line) == 0:
                        continue
                    if len(line) == 1:
                        self.timefld.append(float(line[0]))
                        self.nrgcols.append([[] for _ in range(self.pnt.n_spec)])
                        ispec = 0  # Reset index for the species at the current time step
                    elif len(line) == n_col:
                        self.nrgcols[-1][ispec] = line
                        ispec += 1
                    else:
                        raise IOError("Incorrect number of columns")
        except IOError:
            sys.exit("IOError: nrg file does not exist or has"
                     " wrong number of columns: {}".format(self.filename))
        self.timefld = np.array(self.timefld)
        first_time, last_time = self.get_minmaxtime()
        if len(self.timefld) != len(set(self.timefld)):
            raise RuntimeError("Error: {} contains 2 blocks with identical"
                               " timestamp".format(self.filename))
        if self.starttime == -1 or (0 < self.starttime < first_time):
            print("Using first time present in nrg data for starttime")
            self.starttime = first_time
        if self.endtime == -1:
            print("Using first time present in nrg data for endtime")
            self.endtime = first_time
        if self.starttime == -2:
            print("Using last time present in nrg data for starttime")
            self.starttime = last_time
        if (self.endtime == -2) or (self.endtime > last_time):
            print("Using last time present in nrg data for endtime")
            self.endtime = last_time
        if (self.endtime < first_time) or (self.starttime > last_time):
            print("Time window not contained in nrg data")
            return
        print(("starttime={}, endtime={}, "
               "first_time={}, last_time={}".format(self.starttime, self.endtime,
                                                    first_time, last_time)))
        #  For a single time, find element closest to given input time
        if self.starttime == self.endtime:
            pos = np.array([bisect_left(self.timefld, self.starttime)])
        else:
            pos = np.where((self.timefld >= float(self.starttime)) &
                           (self.timefld <= self.endtime))[0]

        self.timefld = self.timefld[pos]
        self.nrgcols = np.array(self.nrgcols).astype(float, copy=False)
        # Reduce the nrgcols to only the required time frame
        self.nrgcols = self.nrgcols[pos, ...]
        self.isdatapresent = True


def gluenrgdata(nrglist):
    """ Function to combine a list of nrg files from a continuation run.

    :param nrglist: the list of Nrgdata objects to combine
    :returns: the combined Nrgdata object
    """
    # TODO: More sanity checks for agreeing parameters in the list
    if len(nrglist) == 1:
        return nrglist[0]
    result = nrglist.pop(0)
    for nrg in nrglist:
        # When times overlap, give following Nrgdata preference
        while result.timefld[-1] > nrg.timefld[0]:
            del result.timefld[-1]
        result.nrgcols = result.nrgcols[:len(result.timefld), ...]
        result.timefld = np.concatenate((result.timefld, nrg.timefld))
        result.nrgcols = np.concatenate((result.nrgcols, nrg.nrgcols))
    result.endtime = nrglist[-1].endtime
    del nrglist
    return result


class PlotNrgdata(Plotting):
    """ PlotNrgdata: Class to generate plots from nrgdata objects

    So far only implements a routine to calculate windowed mean and errors
    """
    def __init__(self, nrgdata):
        super().__init__()
        self.nrgdata = nrgdata
        self.n_col1 = 4
        self.n_col2 = 6
        if self.n_col1+self.n_col2 != nrgdata[0].n_col:
            raise RuntimeError("The number of nrg columsn in the Nrgdata class does not agree "
                               "with the current implementation in the plotting routine")
        for nrgd in nrgdata:
            if nrgd.isdatapresent is False:
                print("No nrg data present in object ", nrgdata.fext)

    def calcmean_err(self):
        """ Calculate and print the mean and error estimate for each nrg column """
        for nrgd in self.nrgdata:
            mean_nrg = err.mytrapz(nrgd.nrgcols, nrgd.timefld)
            err_nrg, ctimes_nrg = err.windowerr(nrgd.nrgcols, nrgd.timefld)
            print("Mean and errors for ", nrgd.filename)
            for i_sp in range(nrgd.pnt.n_spec):
                print("Species: ", i_sp, nrgd.specname[i_sp])
                print("============================")
                for col in range(nrgd.n_col):
                    print(nrgd.colnames[col]+": ", mean_nrg[i_sp, col], "+-", err_nrg[i_sp, col])
                print("corr times: ", ctimes_nrg[i_sp, :])
                print("============================")

    def plot_ttrace(self, poutput=False):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        for nrgd in self.nrgdata:
            mean_nrg = err.mytrapz(nrgd.nrgcols, nrgd.timefld)
            err_nrg, ctimes_nrg = err.windowerr(nrgd.nrgcols, nrgd.timefld)
            for i_sp in range(nrgd.pnt.n_spec):
                # Plot first group: n
                for col in range(self.n_col1):
                    texlabel = r"$" + nrgd.colnames[col] + r"$"
                    base_line, = ax1.plot(nrgd.timefld, nrgd.nrgcols[:, i_sp, col], label=texlabel)
                    if err_nrg[i_sp, col] != 0.0:
                        ax1.axhline(y=mean_nrg[i_sp, col], linestyle='--',
                                    color=base_line.get_color())
                        ax1.axhspan(mean_nrg[i_sp, col]-err_nrg[i_sp, col],
                                    mean_nrg[i_sp, col]+err_nrg[i_sp, col],
                                    alpha=0.3, facecolor=base_line.get_color())
                ax1.set_xlabel(r'$t$', fontsize=self.xyfs)
                ax1.legend()
                fig1.tight_layout()
                for col in range(self.n_col1, self.n_col1 + self.n_col2):
                    texlabel = r"$" + nrgd.colnames[col] + r"$"
                    base_line, = ax2.plot(nrgd.timefld, nrgd.nrgcols[:, i_sp, col], label=texlabel)
                    ax2.axhline(y=mean_nrg[i_sp, col], linestyle='--', color=base_line.get_color())
                    ax2.axhspan(mean_nrg[i_sp, col]-err_nrg[i_sp, col],
                                mean_nrg[i_sp, col]+err_nrg[i_sp, col],
                                alpha=0.3, facecolor=base_line.get_color())
                ax2.set_xlabel(r'$t$', fontsize=self.xyfs)
                ax2.legend()
                fig2.tight_layout()
                if poutput == 1:
                    fig1.frameon = False
                    fig1.savefig("nrg_ttrace1{}{}.pdf".format(nrgd.specname[i_sp], nrgd.fext))
                    fig2.frameon = False
                    fig2.savefig("nrg_ttrace2{}{}.pdf".format(nrgd.specname[i_sp], nrgd.fext))

    @classmethod
    def show(cls):
        plt.show()
