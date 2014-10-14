"""
nonperturbative.py:
    Mostly just contains the nonperturbative_set and nonperturbative classes.

"""


import numpy as np
import os
import pandas as pd
import scipy.signal
from common import import_petsc_vec
from units import atomic

import fourier


class NonperturbativeSet(object):

    '''A set of nonperturbative (full TDSE) runs,
    specifically their chi values.  '''

    def __init__(self, folder, n=1, wanted="susceptibility", window=scipy.signal.boxcar):
        self.data = None
        self.folders = []
        self.n = n
        self.window = window

        if wanted not in ["susceptibility","peak_dipole","efield","integrated_harmonic_power","peak_harmonic_power"]:
            raise Exception("wanted value <" + str(wanted) + "> not known")

        self.wanted = wanted

        os.path.walk(folder, NonperturbativeSet.__visit, self)
        self.data = self.data.groupby(self.data.index).sum()

        self.data.sort_index(inplace=True)

    def __visit(self, dirname, names):
        """ a helper function for the directory walk """
        if not "wf_final.dat" in names:
            return
        # try:
        if self.data is None:
            self.data = pd.DataFrame(
                NonperturbativeSet.Nonperturbative(dirname, self.n, self.wanted, self.window).data)
            self.folders.append([dirname])
        else:
            new_data = NonperturbativeSet.Nonperturbative(dirname, self.n, self.wanted, self.window).data
            self.data = pd.concat(
                [self.data, new_data])
            self.folders.append([dirname])
        # except Exception as inst:
            #print("Failed",dirname, inst)
            # return

    def nl_chis(self):
        """finds the nonlinear part of the susceptibility for each
        individual column by subtracting out the y-intercept of a
        linear fit of the first two points in intensity.

        Returns:
            a new dataframe with just the nonlinear part of the susceptibility
        """
        nl_susc = self.data.copy()
        for column in nl_susc.columns:
            slope = (nl_susc[column].iloc[1] - nl_susc[column].iloc[0]) / \
                (nl_susc[column].index[1] - nl_susc[column].index[0])
            _y0 = - \
                (nl_susc[column].index[0] * slope - nl_susc[column].iloc[0])
            nl_susc[column] = nl_susc[column] - _y0
        return nl_susc

    def ionizations(self, zero=-1):
        """ Find the ionization for the set of runs.
        Calls Nonperturbative.ionization(n) for each run
        in the set, and puts them all in a DataFrame.

        Args:
            n: the zero of the field to find the ionization in.

        Returns:
            a pandas single column DataFrame with a index
            parameterized on intensity, wavelength and cycles.
        """
        ionization = {}
        for i in self.folders:
            run = self.Nonperturbative(i[0])
            ionization[(run.intensity, run.wavelength, run.cycles)] = float(
                run.ionization(zero))
        index = pd.MultiIndex.from_tuples(
            ionization.keys(), names=["intensity", "wavelength", "cycles"])
        return pd.DataFrame(ionization.values(), index=index).sort_index()

    class Nonperturbative(object):

        """ nonperturbative object representing a single nonperturbative run
        allows access to properties of the run, and to calculations done on
        the run
        """

        def __init__(self, folder, n=1, wanted="susceptibility" ,window=scipy.signal.boxcar):
            self.window=window
            if not os.path.exists(os.path.join(folder, "wf_final.dat")):
                raise Exception("folder doesn't have wf_final.dat", folder)
            self.folder = folder
            with open(os.path.join(folder, "que"), 'r') as que_file:
                self.intensity = float()
                line_of_interest = str()
                for line in que_file:
                    if line.startswith("mpirun"):
                        line_of_interest = line
                line_of_interest = line_of_interest.replace('\t', '   ')
                for element in line_of_interest.split('    '):
                    if element.startswith("-laser_intensity"):
                        self.intensity = float(element.split(" ")[1])

            with open(os.path.join(folder, "Laser.config"), 'r') as laser_conf:
                self.cycles = int()
                self.dt = float()
                self.wavelength = float()
                for line in laser_conf:
                    if line.startswith("-laser_cycles"):
                        self.cycles = int(line.split(" ")[1])
                    if line.startswith("-laser_lambda"):
                        self.wavelength = float(line.split(" ")[1])
                    if line.startswith("-laser_dt "):
                        self.dt = float(line.split(" ")[1])

            if wanted == "susceptibility":
                self.chi = self.harmonic(n, self.window)
            elif wanted == "peak_dipole":
                self.chi = self.dipole(self.window)(n * atomic.from_wavelength(self.wavelength))
            elif wanted == "peak_harmonic_power":
                self.chi = np.square(np.abs(self.dipole(self.window)(n * atomic.from_wavelength(self.wavelength))))
            elif wanted == "efield":
                self.chi = self.efield(self.window)(n * atomic.from_wavelength(self.wavelength))
            elif wanted == "integrated_harmonic_power":
                self.chi = self.dipole(self.window).integrated_freq(lambda x: np.abs(x)**2, ((n-1) * atomic.from_wavelength(self.wavelength)),((n+1) * atomic.from_wavelength(self.wavelength)) )
            else:
                raise Exception("wanted value <" + str(wanted) + "> not known")


            mindex = pd.MultiIndex.from_arrays(
                [[self.cycles], [self.wavelength]], 
                names=["cycles", "wavelength"])
            self.data = pd.DataFrame(
                self.chi, columns=mindex, index=[self.intensity])
            self.data.index.name = "intensity"

        # def dipole_t(self):
            # with open(os.path.join(self.folder, "dipole.dat"), 'rb') as f:
            #dp = np.fromfile(f, 'd', -1)
            # with open(os.path.join(self.folder, "time.dat"), 'rb') as f:
            #time = np.fromfile(f, 'd', -1)

            # return (time, dp)

        def dipole(self, window=None, t=None):
            """
            returns a Fourier object of the dipole moment selected of the run.
            t = None : the regular, total, dipole moment.
            t = "ab","ba","bb", etc.  the section of the total dipole moment.  if t[0] > t[1], flip and take the complex conjugate.
            t = [(x, "ab"), (y, "bc")], etc.  x(dipole moment of "ab") + y(dipole moment of "bc"), etc.
            """

            with open(os.path.join(self.folder, "time.dat"), 'rb') as time_f:
                time_f.seek(0, os.SEEK_END)
                timesize = time_f.tell()
                time_f.seek(0)
                time = np.fromfile(time_f, 'd', -1)

            files = []
            if t is None:
                files = [(lambda x: x, "dipole.dat")]
            if t is not None:
                if not isinstance(t, (list, tuple)):
                    files = [(lambda x: x, "dipole_" + t + ".dat")]
                else:
                    for fun, f in t:
                        if (f is None):
                            files.append( (fun, "dipole.dat") )
                        else:
                            if (f[0] > f[1]):
                                files.append( (lambda x: fun(np.conj(x)), "dipole_" + f[1] + f[0] + ".dat") )
                            else:
                                files.append( (fun, "dipole_" + f + ".dat"))
            print files
            for func, f in files:
                dp = np.zeros(timesize/8, dtype='D')
                with open(os.path.join(self.folder, f), 'rb') as dipolef:
                    dipolef.seek(0, os.SEEK_END)
                    dipolesize = dipolef.tell()
                    dipolef.seek(0)
                    if dipolesize == 2 * timesize:
                        dp = np.add(dp,func(np.fromfile(dipolef, 'D', -1)))
                    elif dipolesize == timesize:
                        dp = np.add(dp,func(np.fromfile(dipolef, 'd', -1)))
                    else:
                        raise Exception("dipole: dipole file does not have a filesize equal to, or double that of the time file")
            dp *= -1
            if window:
                return fourier.Fourier(time, dp, window)
            else:
                return fourier.Fourier(time, dp, self.window)

        def efield(self, window=None):
            """
            returns a Fourier object of the efield of the run.
            """
            with open(os.path.join(self.folder, "efield.dat"), 'rb') as ef_file:
                ef = np.fromfile(ef_file, 'd', -1)
            with open(os.path.join(self.folder, "time.dat"), 'rb') as time_f:
                time = np.fromfile(time_f, 'd', -1)
            if window:
                return fourier.Fourier(time, ef, window)
            else:
                return fourier.Fourier(time, ef, self.window)

        # def efield_t(self):
            # with open(os.path.join(self.folder, "efield.dat"), 'rb') as f:
            #ef = np.fromfile(f, 'd', -1)
            # with open(os.path.join(self.folder, "time.dat"), 'rb') as f:
            #time = np.fromfile(f, 'd', -1)
            # return (time,ef)

        def harmonic(self, order=1, window=scipy.signal.boxcar):
            """
            return the susceptibility for a specific harmonic `order`.
            """
            dipole = self.dipole(window)
            dipole = dipole(order * atomic.from_wavelength(self.wavelength))
            efield = self.efield(window)
            # We want the susceptibility at order omega from omega
            efield = efield(atomic.from_wavelength(self.wavelength))
            return dipole / efield

        def get_prototype(self):
            """get the prototype from the csv file in the run's
            directory.  returns a pandas MultiIndex"""
            n = []
            l = []
            j = []
            m = []
            e = []
            with open(os.path.join(self.folder, "prototype.csv"), 'r') as prototype_f:
                for line in prototype_f:
                    i = line.split(',')
                    n.append(int(i[0]))
                    l.append(int(i[1]))
                    j.append(float(int(i[2])) / 2.)
                    m.append(int(i[3]))
                    e.append(float(i[4]))
            return pd.MultiIndex.from_arrays([n, l, j, m, e], names=["n", "l", "j", "m", "e"])

        def wf(self, zero=-1):
            """ The wavefunction at the zero of the field number `zero`.
            """
            if zero == -1:
                wavefn = import_petsc_vec(
                    os.path.join(self.folder, "wf_final.dat"))
                return pd.DataFrame({"wf": wavefn}, index=self.get_prototype())
            elif zero < self.cycles * 2:
                wavefn = import_petsc_vec(
                    os.path.join(self.folder, "wf_" + str(zero) + ".dat"))
                return pd.DataFrame({"wf": wavefn}, index=self.get_prototype())
            else:
                raise Exception(zero, "n not less than ", self.cycles)

        def gs_population(self, zero=-1):
            """get the ground state population at the zero of the
            field represented by `zero`."""
            wavefn = self.wf(zero)
            return wavefn.apply(lambda x: abs(x) ** 2).query('n == 1').iloc[0]

        def bound_population(self, zero=-1):
            """get the bound state population (including the gs) at
            the zero of the field represented by `zero`."""
            wavefn = self.wf(zero)
            return wavefn.apply(lambda x: abs(x) ** 2).query('e < 0').sum()

        def ionization(self, zero=-1):
            """get the ionized population at the zero of the
            field represented by `zero`."""
            wavefn = self.wf(zero)
            absorbed = 1 - wavefn.apply(lambda x: abs(x) ** 2).sum()
            return absorbed + wavefn.apply(lambda x: abs(x) ** 2).query('e > 0').sum()

    @staticmethod
    def get_prototype(filename="prototype.csv"):
        """get the prototype from a csv file
        returns a pandas MultiIndex"""
        with open(filename, 'r') as prototype_f:
            n = []
            l = []
            j = []
            m = []
            e = []
            for line in prototype_f:
                i = line.split(',')
                n.append(int(i[0]))
                l.append(int(i[1]))
                j.append(float(int(i[2])) / 2.)
                m.append(int(i[3]))
                e.append(float(i[4]))
            return pd.MultiIndex.from_arrays([n, l, j, m, e], names=["n,l,j,m,e"])
