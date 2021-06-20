from __future__ import division, print_function
from numpy import pi
import numpy as np
from scipy import interpolate


def intensity_calc(p, radi):
    '''
    calculate the beam intensity
    :param p: beam power [mW]
    :param radi: beam radius 1/e^2 [mm]
    :return: peak intensity mW/cm^2
    '''
    ret_val = 0.1 * (2 *p * 1e-3)/(pi * (radi * 1e-3)**2)
    return ret_val

def half_max_roots(x, y):
    half_max = np.max(y) / 2
    spline = interpolate.UnivariateSpline(x, y - half_max, s=0)
    r1, r2 = spline.roots()
    return half_max, r1, r2


def full_width_at_half_max(x, y):
    half_max, r1, r2 = half_max_roots(x, y)
    return r2 - r1