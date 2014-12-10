"""
window module.  Contains a number of windows for fourier transforming things.
"""

from __future__ import print_function, division
import numpy as np


class CosWindowFunction(object):
    """
    A window function represented as a sum of cosines.

    required attributes:

    c : a list of the coefficients for a sum of cosines.
    rov : recommended overlap
    """

    def __init__(self, c, rov):
        self.c = c
        self.rov = rov

    def __call__(self, M, sym=True):
        if M < 1:
            return np.array([])
        if M == 1:
            return np.ones(1, 'd')
        odd = M % 2
        if not sym and not odd:
            M = M + 1
        n = np.arange(0, M)
        fac = n * 2 * np.pi / (M - 1.0)

        w = np.zeros_like(len(n))
        for i, c in enumerate(self.c):
            if i == 0:
                w += c[i]
            else:
                w += c[i] * np.cos(i * fac)

        if not sym and not odd:
            w = w[:-1]
        return w


sft3f = CosWindowFunction([0.26526, -.5, .23474], .667)
sft4f = CosWindowFunction([0.21706, -.42103, .28294, -.07897], .655)
sft5f = CosWindowFunction([0.1881, -.36923, .28702, 0.13077, .02488], .75)
sft3m = CosWindowFunction([0.28235, -.52105, .19659], .721)
sft4m = CosWindowFunction([.241906, -.460841, .255381, -.041872], .785)
sft5m = CosWindowFunction([.209671, -.407331, .281225, -.092669, .0091036], .76)
hft90d = CosWindowFunction([1, -1.942604, 1.340318, -.440811, .043097], .76)
hft95 = CosWindowFunction([1, -1.9383379, 1.3045202, -.4028270, .0350665], .756)
hft116d = CosWindowFunction([1, -1.9575375, 1.4780705, -0.6367431, .1228389,
                             -.0066288], .782)
hft144d = CosWindowFunction([1, -1.96760033, 1.57983607, -0.81123644, .22583558,
                             -.02773848, .00090360], .799)
hft169d = CosWindowFunction([1, -1.97441842, 1.65409888, -0.95788186, 0.33673420,
                             -0.06364621, 0.00521942, -0.00010599], .812)
hft248d = CosWindowFunction([1, -1.985844164102, 1.791176438506, -1.282075284005,
                             0.667777530266, -0.240160796576, 0.056656381764,
                             -0.008134974479, 0.000624544650, -0.000019808998,
                             0.000000132974], .841)
nutall4 = CosWindowFunction([.3125, -.46875, .1875, -.03125], .704)

