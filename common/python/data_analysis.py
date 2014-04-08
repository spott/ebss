#!/usr/bin/env python3

import numpy as np
import os
import pandas as pd
import math
import functools as ft


''' import a binary file of type 'd': '''
def get_file( filename, dtype='d' ):
    with open(os.path.join(filename), 'rb') as f:
        npy = np.fromfile(f, dtype)
    return npy

''' import a PETSc vec: '''
def import_petsc_vec( filename ):
    with open(filename, "rb") as f:
        byte = f.read(8)
        size = struct.unpack('>ii',byte)[1]
        npy = np.fromfile(f, '>D', size)
        return npy

class atomic(object):
    ''' class that contains atomic unit conversions and constants '''
    pressure = 3.44386e-9
    kT = 0.000942716
    intensity = 3.5094452e16
    c = 137.035999074
    e0 = 0.0795775

    @staticmethod
    def from_wavelength( l ):
        return 45.56335/l

    @staticmethod
    def averaged_intensity( i, n ):
        return (math.gamma(.5 + 2.*n) / (np.sqrt(np.pi) * math.gamma(1.+2.*n))) * (i/atomic.intensity)**n

class perturbative_set(object):
    ''' holds the dataframe that has all the perturbative calculations. '''

    chi_labels = ["1","-1,1,1","-1,-1,1,1,1","-1,-1,-1,1,1,1,1","-1,-1,-1,-1,1,1,1,1,1","-1,-1,-1,-1,-1,1,1,1,1,1,1"]
    def __init__(self, folder = None, df = None ):
        ''' parent folder is the input, containing at least a "rmax/npoints/nmax_nonlinear_DNWM" structure '''
        if (folder is None and df is None):
            raise Exception("no arguments passed!")
        elif folder is not None:
            self.data = None

            self.rmax = []
            self.nmax = []
            self.npoints = []
            for rmax in filter( lambda x: os.path.isdir(os.path.join(folder,x)), os.listdir(folder)):
                for npoints in filter( 
                        lambda x: os.path.isdir(os.path.join(folder,rmax,x)), 
                        os.listdir(os.path.join(folder,rmax))):
                    for nmax_folder in filter( 
                            lambda x: x.count("nonlinear") == 1 
                                and os.path.isdir(os.path.join(folder,rmax,npoints,x)), 
                            os.listdir(os.path.join(folder,rmax,npoints))):
                        #try:
                        self.rmax += [int(rmax)]
                        self.nmax += [int(nmax_folder.split("_")[0])]
                        self.npoints += [int(npoints)]
                        try:
                            d = perturbative_set.perturbative( os.path.join(folder,rmax,npoints,nmax_folder), 
                                    rmax, npoints, nmax_folder.split("_")[0])
                            if self.data is None:
                                self.data = d
                            else:
                                self.data = self.data.join(d)
                        except:
                            pass
                        #except Exception as inst:
                            #print(inst)
                            #print(os.path.join(folder,rmax,npoints,nmax_folder))
                            #print("==== failed ====")
                            #continue
        else:
            self.data = df

    @staticmethod
    def perturbative(folder, rmax, npoints, nmax):
        ''' create the object, populate all the information on it. '''
        #import the frequencies:
        freqs = [ str(i)[:5] for i in get_file( os.path.join(folder, "frequencies.dat") )]
        rmax = int(rmax)
        npoints = int(npoints)
        nmax = int(nmax)

        #import the imaginary component:
        imgs = []
        with open(os.path.join(folder,"que"), 'r') as f:
            terms = []
            for line in f:
                if (line.strip().startswith("mpiexec")):
                    terms = line.split(" ")
                    break

            for i in range(len(terms)):
                if (terms[i] =="-nonlinear_img"):
                    imgs = [ float(i) for i in terms[i+1].split(',') ]
                    break

        chis = {}
        for i in filter( lambda x: x.count("chi") == 1, os.listdir(folder) ):
            s = i.split("_")
            if (len(s) == 2):
                s = s[1].split(".")[0]
            else:
                s = s[1:-1]
            chis[",".join(s)] = get_file( os.path.join(folder,i), 'D')
        data = pd.DataFrame( chis)
        data.index = pd.MultiIndex.from_product([imgs,freqs], names=["epsilon","frequency"])
        data.columns = pd.MultiIndex.from_product([[rmax], [npoints], [nmax], data.columns], names=["rmax","npoints","nmax","chi"])
        return data
        #self.data.columns.set_levels = 

    def chis(self, freq="0.056"):
        return self.data.xs(freq, level="frequency").xs(max(self.nmax), level="nmax", axis=1).xs(max(self.npoints), level="npoints", axis=1)[max(self.rmax)].iloc[-1].T

    def dnwm(self, freq="0.056"):
        return self.data.xs(freq, level="frequency").xs(max(self.nmax), level="nmax", axis=1).xs(max(self.npoints), level="npoints", axis=1)[max(self.rmax)][["1","-1,1,1","-1,-1,1,1,1","-1,-1,-1,1,1,1,1","-1,-1,-1,-1,1,1,1,1,1","-1,-1,-1,-1,-1,1,1,1,1,1,1"]].iloc[-1]

    def nlchi_vs_intensity(self, intensities, freq="0.056"):
        l = {}
        dnwm_data = self.dnwm(freq)
        for m in range(1,len(perturbative_set.chi_labels)):
            #print(int(m/2))
            ch = lambda i : sum([ dnwm_data[perturbative_set.chi_labels[j]] * atomic.averaged_intensity(i,j) for j in range(1,m+1) ])
            l["chi" + str(m*2+1)] = [ ch(i) for i in intensities ]

        return pd.DataFrame(l, index=intensities)

    def fractional_chis(self, intensities, compared_to=3, freq="0.056"):
        l = {}
        dnwm_data = self.dnwm(freq)
        for m in range(int(compared_to/2)+1,len(perturbative_set.chi_labels)):
            #print(int(m/2))
            l["chi" + str(m*2+1)] = [ (dnwm_data[perturbative_set.chi_labels[m]] * atomic.averaged_intensity(i,m)) / (dnwm_data[int(compared_to/2)] * atomic.averaged_intensity(i,int(compared_to/2))) for i in intensities ]

        return pd.DataFrame(l, index=intensities)

    def rmax_convergence(self, nmax=None, npoints=None, epsilon=None):
        if nmax is None:
            nmax = max(self.nmax)
        if npoints is None:
            npoints = max(self.npoints)
        if epsilon is None:
            epsilon = min(self.data.index.get_level_values("epsilon"))
        return self.data.xs(epsilon,level="epsilon").xs(nmax,level="nmax",axis=1).xs(npoints, level="npoints",axis=1)

    def nmax_convergence(self, rmax=None, npoints=None, epsilon=None):
        if rmax is None:
            rmax = max(self.rmax)
        if npoints is None:
            npoints = max(self.npoints)
        if epsilon is None:
            epsilon = min(self.data.index.get_level_values("epsilon"))
        return self.data.xs(epsilon,level="epsilon").xs(npoints, level="npoints",axis=1)[rmax]

    def npoints_convergence(self, rmax=None, nmax=None, epsilon=None):
        if rmax is None:
            rmax = max(self.rmax)
        if nmax is None:
            nmax = max(self.nmax)
        if epsilon is None:
            epsilon = min(self.data.index.get_level_values("epsilon"))
        return self.data.xs(epsilon,level="epsilon").xs(nmax, level="nmax",axis=1)[rmax]

class nonperturbative_set(object):
    ''' the perturbative set of values '''

    def __init__(self, folder):
        self.data = None
        self.folders = []

        os.path.walk(folder, nonperturbative_set.visit, self)
        self.data = self.data.groupby(self.data.index).sum()

        self.data.sort_index(inplace=True)

    def visit(self, dirname, names):
        if not "wf_final.dat" in names:
            return
        try:
            if (self.data is None):
                self.data = pd.DataFrame(nonperturbative_set.nonperturbative( dirname ).data)
                self.folders.append([dirname])
            else:
                self.data = pd.concat([self.data,nonperturbative_set.nonperturbative( dirname ).data])
                self.folders.append([dirname])
        except:
            return

    def nl_chis(self):
        new_data = self.data.copy()
        for c in new_data.columns:
            m = (new_data[c].iloc[1] - new_data[c].iloc[0]) / (new_data[c].index[1] - new_data[c].index[0])
            y0 = - (new_data[c].index[0] * m - new_data[c].iloc[0])
            new_data[c] = new_data[c] - y0
        return new_data

    class nonperturbative(object):

        def __init__(self, folder):
            if ( not os.path.exists(os.path.join(folder,"wf_final.dat"))):
                raise Exception("folder doesn't have wf_final.dat", folder)
            with open(os.path.join(folder,"que"),'r') as q:
                self.intensity = float()
                line_of_interest = str()
                for l in q:
                    if (l.startswith("mpirun")):
                        line_of_interest = l
                line_of_interest = line_of_interest.replace('\t','   ')
                for c in line_of_interest.split('    '):
                    if (c.startswith("-laser_intensity")):
                        self.intensity = float(c.split(" ")[1])

            with open(os.path.join(folder,"Laser.config"),'r') as lc:
                self.cycles = int()
                self.dt = float()
                self.wavelength = float()
                for l in lc:
                    if (l.startswith("-laser_cycles")):
                        self.cycles = int(l.split(" ")[1])
                    if (l.startswith("-laser_lambda")):
                        self.wavelength = float(l.split(" ")[1])
                    if (l.startswith("-laser_dt ")):
                        self.dt = float(l.split(" ")[1])

            with open(os.path.join(folder, "dipole.dat"), 'rb') as f:
                temp = np.fromfile(f, 'd', -1)
                ft = (2. * np.abs(np.fft.rfft(temp)) / len(temp))
                time = self.cycles * 2. * np.pi / atomic.from_wavelength(self.wavelength)
                df = 2. * np.pi / ( time )
                i = int(round( atomic.from_wavelength(self.wavelength) / df))
                self.chi = 2. * ft[i] / np.sqrt(self.intensity / atomic.intensity)

            mi = pd.MultiIndex.from_arrays([[self.cycles], [self.wavelength]], names=["cycles","wavelength"])
            self.data = pd.DataFrame( self.chi, columns = mi, index = [self.intensity] )
            self.data.index.name = "intensity"


    @staticmethod
    #@ft.lru_cache(maxsize=10)
    def get_prototype( filename = "prototype.csv" ):
        prototype_f = file( filename )
        n = []
        l = []
        j = []
        m = []
        e = []
        for l in prototype_f:
            i = l.split(',')
            n.append(int(i[0]))
            l.append(int(i[1]))
            j.append(float( int( i[2] ) ) /2.)
            m.append(int(i[3]))
            e.append(float(i[4]))
        return pd.MultiIndex.from_arrays([n,l,j,m,e], names=["n,l,j,m,e"])


class basis(object):
    def __init__(self, folder="/lustre/janus_scratch/ansp6066/basis/hydrogen_new_160kp_n1100_l100_r500"):
        self.folder = folder
        self.grid = get_file( os.path.join(folder, "grid.dat"))
        self.points = len(self.grid)
    
    def wf(self, n, l):
        with open(os.path.join(self.folder, "l_" + str(l) + ".dat"), 'rb') as f:
            npy = np.fromfile(f, 'd', (n-(l))*self.points)[-self.points:]
        return npy
    
    def prototype(self):
        with open(os.path.join(self.folder, "prototype.dat"), 'rb') as f:
            dt = numpy.dtype([('n', numpy.int32), ('l', numpy.int32), ('j', numpy.int32), 
                          ('m', numpy.int32), ('e', numpy.complex128)])
            npy = numpy.fromfile(f, dt)
        return npy
