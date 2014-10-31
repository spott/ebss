"""
fourier module. contains the Fourier class.
"""
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
        if data.dtype == np.dtype('d'):
            # real!
            self.fdata = 2 * np.fft.rfft(self.tdata) / len(self.tdata)
            self.freq = np.fft.rfftfreq(len(self.tdata), (self.time[1]-self.time[0]) / (2. * np.pi))

        elif data.dtype == np.dtype('D'):
            self.fdata = 2 * np.fft.fft(self.tdata) / len(self.tdata)
            self.freq = np.fft.fftfreq(len(self.tdata), (self.time[1]-self.time[0]) / (2. * np.pi))
        else:
            raise Exception("Fourier.__init__: Not using double complex, or double... this is wrong")

        self.df = self.freq[2] - self.freq[1]

    def __call__(self, freq):
        point = freq/self.df
        i = int(np.round(point))
        #print "error in frequency: ", point - i
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
        return scipy.interpolate.InterpolatedUnivariateSpline(self.freq, fn(self.fdata),k=1).integral(a,b)

    def stft(self, window_size, overlap, window_fn=scipy.signal.flattop):
        hop = window_size / overlap
        # better reconstruction with this trick +1)[:-1]
        w = window_fn(window_size)#[:-1]
        df = 2. * np.pi / self.time[window_size]
        time = np.array([self.time[j+(window_size)/(2)] for j in range(0, len(self.time)-window_size, hop)])
        chi = np.array([np.fft.rfft(w*self.tdata[i:i+window_size]) * 2. / window_size for i in range(0, len(self.tdata)-window_size, hop)])
        assert len(time) == len(chi), "stft: time and chi have different lengths"
        return (time, chi)

