"""
fourier module. contains the Fourier class.
"""
from __future__ import print_function
import numpy as np
import scipy.signal
import scipy.interpolate

class Fourier(object):
    """
    A Fourier object takes a time series and a data series and takes the fourier
    transform of them and stores both in the object for later information 
    retrieval.

    Attributes:
        time: the list of times.
        tdata: the data corresponding to the list of times.
        freq: the energies of the fourier transform.
        fdata: the fourier transform data
        df: the spacing between the frequency points.
    """

    def __init__(self, time, data, window=scipy.signal.boxcar):

        self.time = np.array(time)
        self.tdata = window(len(data)) * np.array(data)
        assert len(self.time) == len(self.tdata), "Fourier.__init__: time and tdata must have same length"
        #self.fdata = None
        #self.freq = None

    @property
    def fdata(self):
        if self.fdata_ is None:
            if tdata.dtype == np.dtype('d'):
                # real!
                self.fdata_ = 2 * np.fft.rfft(self.tdata) / len(self.tdata)
                self.freq_ = np.fft.rfftfreq(len(self.tdata), (self.time[1]-self.time[0]) / (2. * np.pi))

            elif data.dtype == np.dtype('D'):
                self.fdata_ = 2 * np.fft.fft(self.tdata) / len(self.tdata)
                self.freq_ = np.fft.fftfreq(len(self.tdata), (self.time[1]-self.time[0]) / (2. * np.pi))
            else:
                raise Exception("Fourier.__init__: Not using double complex, or double... this is wrong")
            self.df = self.freq[2] - self.freq[1]
        return self.fdata_

    @property
    def freq(self):
        if self.fdata_ is None:
            self.fdata()
        return self.freq_

    def __call__(self, freq):
        point = freq/self.df
        i = int(np.round(point))
        return self.fdata[i]

    def interpolated_freq(self, freq):
        """ returns the frequency asked for, interpolated instead of
        rounded"""
        i = int(np.round(freq/self.df))
        iinitial = i - 10 if i > 10 else 0
        ifinal = i + 10 if len(self.freq) - i < 10 else len(self.freq)
        return scipy.interpolate.InterpolatedUnivariateSpline(self.freq[iinitial:ifinal], np.abs(self.fdata[iinitial:ifinal]))(freq)

    def integrated_freq(self, fn, a, b):
        """ integrates the frequency over the range. Uses an interpolant"""
        return scipy.interpolate.InterpolatedUnivariateSpline(self.freq, fn(self.fdata), k=1).integral(a,b)

    def stft(self, window_size, overlap, window_fn=scipy.signal.flattop):
        hop = window_size / overlap
        # better reconstruction with this trick +1)[:-1]
        w = window_fn(window_size)#[:-1]
        df = 2. * np.pi / self.time[window_size]
        time = np.array([self.time[j+(window_size)/(2)] for j in range(0, len(self.time)-window_size, hop)])
        #remove the DC average
        chi = np.array([np.fft.rfft(w*(self.tdata[i:i+window_size] 
                                    - np.average(self.tdata[i:i+window_size])))
                        * 2. / window_size for i in range(0, len(self.tdata)-window_size, hop)])
        assert len(time) == len(chi), "stft: time and chi have different lengths"
        return (time, chi)


"""
A number of helper functions. for dealing with fourier data
"""

def get_stft_data_from_folder(folder, t=None, name="All", window_fn=scipy.signal.flattop):
    """
    given a folder, find the STFT of the data in that folder.
    t = which dipole to get
    name = what to name that column
    window_fn = what window_fn to use for the stft.
    """
    if t is not None:
        print("calculating for " + folder + " " + str(t) + "...", end="")
    elif t is None:
        print("calculating for " + folder + " All...", end="")
    run = data_analysis.NonperturbativeSet.Nonperturbative(folder)
    # number of points of time dependent susceptibility to calc:
    points = 50

    # number of cycles to average over
    cycles = 4

    # 1 period
    period = 2. * pi / energy(run.wavelength)

    dipole = run.dipole(t=t)
    efield = run.efield()
    window = np.argmin( np.abs(dipole.time - period) ) * cycles
    jump = (len(dipole.time) - window) / points

    time, td_dipole = dipole.stft(window, jump, window_fn=window_fn)
    _, td_ef = efield.stft(window, jump, window_fn=window_fn)

    run_df = pd.DataFrame.from_dict({"dipole": td_dipole.T[cycles], "efield": td_ef.T[cycles], "susceptibility": td_dipole.T[cycles] / td_ef.T[cycles]})
    run_df.index = pd.Index(time, name="time")
    run_df.columns = pd.MultiIndex.from_tuples([(run.intensity, name , "dipole"),(run.intensity, name , "efield"),(run.intensity, name , "susceptibility")], names=["intensity", "decomp" ,"value"])
    print("done")
    return run_df

def get_raw_data(folder, t=None, name="All"):
    """
    get the dipole moment and efield as a function of time for:
    t = the dipole to get
    name = what to call that column
    """
    if t is not None:
        print("calculating for " + folder + " " + str(t) + "...", end="")
    elif t is None:
        print("calculating for " + folder + " All...", end="")
    run = data_analysis.NonperturbativeSet.Nonperturbative(folder)

    dipole = run.dipole(t=t)
    efield = run.efield()

    run_df = pd.DataFrame.from_dict({"dipole": dipole.tdata, "efield": efield.tdata})
    run_df.index = pd.Index(dipole.time, name="time")
    run_df.columns = pd.MultiIndex.from_tuples([(run.intensity, name , "dipole"),(run.intensity, name , "efield")], names=["intensity", "decomp" ,"value"])

    print("done")
    return run_df


