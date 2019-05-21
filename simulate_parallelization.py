#!/usr/bin/env python3

class Parallelization(object):
    def __init__(self):
        self.npx=1
        self.npy=1
        self.npz=1
        self.npv=1
        self.npw=1
        self.nps=1

    def total(self):
        return self.npx*self.npy*self.npz*self.npv*self.npw*self.nps

    def output(self):
        print("Parallelization: %u %u %u %u %u %u"%(self.npx,self.npy,self.npz,self.npv,self.npw,self.nps))

class ParallelizationRange(object):
    def __init__(self):
        self.npx={'min':None,'max':None,'fix':None}
        self.npy={'min':None,'max':None,'fix':None}
        self.npz={'min':None,'max':None,'fix':None}
        self.npv={'min':None,'max':None,'fix':None}
        self.npw={'min':None,'max':None,'fix':None}
        self.nps={'min':None,'max':None,'fix':None}

    def setFixValues(self,npx=None,npy=None,npz=None,npv=None,npw=None,nps=None):
        if npx:
            self.npx['fix']=npx
            self.npx['min']=npx
            self.npx['max']=npx

        if npy:
            self.npy['fix']=npy
            self.npy['min']=npy
            self.npy['max']=npy

        if npz:
            self.npz['fix']=npz
            self.npz['min']=npz
            self.npz['max']=npz

        if npv:
            self.npv['fix']=npv
            self.npv['min']=npv
            self.npv['max']=npv

        if npw:
            self.npw['fix']=npw
            self.npw['min']=npw
            self.npw['max']=npw

        if nps:
            self.nps['fix']=nps
            self.nps['min']=nps
            self.nps['max']=nps

    def fillRemainingValues(self,g):
        if not self.npx['fix']:
            self.npx['min']=1
            self.npx['max']=g.nx

        if not self.npy['fix']:
            self.npy['min']=1
            self.npy['max']=g.ny

        if not self.npz['fix']:
            self.npz['min']=1
            self.npz['max']=g.nz

        if not self.npv['fix']:
            self.npv['min']=1
            self.npv['max']=g.nv

        if not self.npw['fix']:
            self.npw['min']=1
            self.npw['max']=g.nw

        if not self.nps['fix']:
            self.nps['min']=1
            self.nps['max']=g.ns

    def getNextParallelization(self,n_mpi):
        p=Parallelization()
        if (self.npx['fix']):
            p.npx=self.npx['fix']
        if (self.npy['fix']):
            p.npy=self.npy['fix']
        if (self.npz['fix']):
            p.npz=self.npz['fix']
        if (self.npv['fix']):
            p.npv=self.npv['fix']
        if (self.npw['fix']):
            p.npw=self.npw['fix']
        if (self.nps['fix']):
            p.nps=self.nps['fix']

        for p.nps in range(self.nps['min'],self.nps['max']+1):
            for p.npw in range(self.npw['min'],self.npw['max']+1):
                for p.npz in range(self.npz['min'],self.npz['max']+1):
                    for p.npv in range(self.npv['min'],self.npv['max']+1):
                        for p.npx in range(self.npx['min'],self.npx['max']+1):
                            for p.npy in range(self.npy['min'],self.npy['max']+1):
                                if p.total() == n_mpi:
                                    yield p

class Grid(object):
    def __init__(self,nx,ny,nz,nv,nw,ns):
        self.nx=nx
        self.ny=ny
        self.nz=nz
        self.nv=nv
        self.nw=nw
        self.ns=ns
        self.nxb=2
        self.evenx= (self.nx+1) % 2

        self.x_local=False
        self.y_local=True
        self.xy_local=self.x_local and self.y_local
        self.nonlinear = True
        self.arakawa = True

    def checkParallelization(self,p):
        if self.ns % p.nps != 0:
            valid_parall=False
            print("species-parallelization failed: n_spec has to be a multiple of n_procs_s!")
        elif (self.nv % p.npv != 0) or (self.nv/p.npv < 2):
            valid_parall=False
            print("v-parallelization failed: nv0 has to be a multiple of n_procs_v and nv0/n_procs_v has to be greater 1!")
        elif self.nw % p.npw != 0:
            valid_parall=False
            print("w-parallelization failed: nw0 has to be a multiple of n_procs_w!")
        elif ((p.npx != 1) and self.x_local) :
            valid_parall = False
            print("parallelization failed: p.npx must be 1 in self.x_local mode!")
        elif ((p.npy != 1) and (not self.y_local)) :
            valid_parall = False
            print("parallelization failed: p.npy must be 1 for self.y_local=False!")
        elif (self.ny % p.npy != 0 ):
            valid_parall=False
            print("y-parallelization failed: nky0 has to be a multiple of p.npy!")
        elif (self.nx % p.npx != 0 ):
            valid_parall=False
            print("x-parallelization failed: nx0 has to be a multiple of p.npx!")
        elif (p.npx*p.npy > self.nx) :
            valid_parall=False
            print("parallelization failed: p.npx*p.npy must not exceed nx0!")
        elif (self.nxb > (self.nx/p.npx)):
            valid_parall = False
            print("x-parallelization failed: number of x grid points is smaller than required number of ghost cells!")
        elif ((self.nz % p.npz != 0 ) or ((self.nz/p.npz<2) and (self.nz > 1))) :
            valid_parall=False
            print("z-parallelization failed: nz0 has to be a multiple of p.npz and nz0/p.npz has to be greater 1!")
        elif (self.nonlinear and ((self.nx+1-self.evenx) % (2*p.npy)!=0)) :
            valid_parall=False
            print("nonlinear y-parallelization failed: mod(nx0,2*p.npy) has to be 0 (nonlinear run)!")
        elif (self.arakawa and (self.nx<2*self.nxb*p.npx*p.npy)):
            valid_parall=False
            print("Arakawa and x,y parallelization failed: nx0 has to be >= 2*nxb*p.npx*p.npy!")
        elif (not self.xy_local and self.nx/p.npx % p.npy !=0) :
            valid_parall=False
            print("x,y parallelization failed: nx0/p.npx has to be a multiple of p.npy!")
        else:
            valid_parall=True

        return valid_parall

def read_parallelization_nml(filename):
    par=ParallelizationRange()
    par.setFixValues(nps=2,npv=1,npy=4)
    return par


def read_grid(filename):
    g=Grid(288,32,24,96,32,2)
    return g




n_mpi=2304*4

pr=read_parallelization_nml('bla')
g=read_grid('bla')
pr.fillRemainingValues(g)

skip=False
for p in pr.getNextParallelization(n_mpi):
    #p.output()
    if g.checkParallelization(p):
        print("successful test of: ",end='')
        p.output()
    else:
        print("\tFailed test of: ",end='')
        p.output()
