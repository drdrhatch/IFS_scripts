#!/usr/bin/env python3

import numpy as np

class Test:
    def __init__(self,testno,timeloop,parall,pv,init_time):
        self.number=testno
        self.timeloop=timeloop
        self.parallelization=parall
        self.perf_vec=pv
        self.init_time=init_time

    def output(self):
        print(("\nTest no. %d:"%self.number))
        print(("\ttime for time loop was %8.3f"%self.timeloop))
        print(("\tparallelization used: %s"%self.parallelization))
        print(("\tperf_vec = ",self.perf_vec))
        print(("\tinit_time= ",self.init_time))

    def get_timeloop(self):
        return self.timeloop

    def get_init_time(self):
        return self.init_time

class Testset:
    def __init__(self,filename,label,linebuffer=[]):
        if (not filename):
            self.filename=None
            self.label=''
            self.linebuffer=linebuffer
        else:
            self.filename=filename
            self.label=label
            print(("Open perf file %s\n"%self.filename))
            fh=open(self.filename,"r")
            self.linebuffer=fh.readlines()
            fh.close()
        self.tests=[]

    def read(self):
        line=self.linebuffer[0]
        tokens=line.split()
        self.testset=tokens[0]
        self.machine=tokens[3]
        line=self.linebuffer[1]
        mo=re.search(r"^started ([0-9.]+), ([0-9:]+) with GIT master ([0-9a-zA-Z]+), GIT branch ([0-9a-zA-Z]+), (.*)$",line)
        if mo:
            self.date=mo.group(1)
            self.time=mo.group(2)
            self.master=mo.group(3)
            self.branch=mo.group(4)
            self.flags=mo.group(5)
        else:
            print(("header of file %s is malformed.\n"%self.filename))

        for line in self.linebuffer[2:]:
            mo=re.search(r"test (\d+): time loop: ([0-9.]+) s \(parallelization: (\d \d \d \d \d \d), perf_vec: ([0-9 ]+)\), init_time ([0-9.]+) s",line)
            if mo:
                testno=int(mo.group(1))
                timeloop=float(mo.group(2))
                parallelization=mo.group(3)
                perf_vec=mo.group(4)
                init_time=float(mo.group(5))
                self.tests.append(Test(testno,timeloop,parallelization,perf_vec,init_time))

    def output(self):
        print("\n\n=============================================")
        print(("File %s has %d tests:\n"%(self.filename,len(self.tests))))
        for t in self.tests:
            t.output()

    def get_timeloop(self):
        return [x.get_timeloop() for x in self.tests]

    def get_init_time(self):
        return [x.get_init_time() for x in self.tests]

    def get_ntests(self):
        return len(self.tests)

    def get_machine(self):
        return self.machine

    def get_testset(self):
        return self.testset

    def get_label(self):
        return self.label

# -----------------------------------------------------------------------------------------------
# ------------------ Main program ---------------
# -----------------------------------------------------------------------------------------------

#from optparse import OptionParser
import re
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Compare two or more output files from the testsuites.')
parser.add_argument('-f',nargs=2,action='append',required=True)
parser.add_argument('--tm',nargs=3,action='append',help='path to perf_data, testset, machine')

myopt=parser.parse_args()

#print("args = ",myopt.filenames)
print(("myopt=",myopt))
files=[]

if myopt.tm:
    for tm in myopt.tm:
        header_line=tm[1]+' testsuite on '+tm[2]
        print(("Look for header_line(",header_line,") in ",tm[0]+"/perf_data"))
        fh=open(tm[0]+"/perf_data","r")
        store_lines=False
        linebuffer=[]
        for line in fh:
            if line.find(header_line)==0:
                linebuffer.append(line)
                store_lines=True
            else:
                if store_lines:
                    if len(line.strip())==0:
                        store_lines=False
                    else:
                        linebuffer.append(line)
        fh.close()
        files.append(Testset(None,label=None,linebuffer=linebuffer))
        files[-1].read()

for (perf_file,label) in myopt.f:
    files.append(Testset(perf_file,label))
    files[-1].read()


# now everything is read an parsed, next is to plot it
# therefore we have to create the right numpy arrays from the parsed data

# first plot is a bar plot with x axis the test number
# and y axis the time (either timeloop or init_time)

nfiles=len(files)
width=0.8/nfiles
colarr=[(0.8,0.0,0.0),(0.1,0.8,0.1),(0.1,0.1,0.8),(0.5,0.5,0.1),(0.1,0.6,0.6)]
fig=plt.figure(1)
labels=[]
counter=0
for tf in files:
    timeloop_arr=np.array(tf.get_timeloop())
    xpos=list(range(tf.get_ntests()))
    print(timeloop_arr)
    #print([x+counter*width for x in xpos])
    labels.append(tf.get_machine()+' '+tf.get_label())
    plt.bar([x+counter*width for x in xpos],height=timeloop_arr,width=width,color=colarr[counter])
    counter = counter +1

plt.legend(labels,loc='upper center')

plt.title("Testsuite runtime comparison for testset "+files[0].get_testset())
plt.show()
