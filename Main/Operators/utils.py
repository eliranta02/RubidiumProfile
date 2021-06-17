from __future__ import division, print_function
from numpy import pi

def intensity_calc(p, radi):
    '''
    calculate the beam intensity
    :param p: beam power [mW]
    :param radi: beam radius 1/e^2 [mm]
    :return: peak intensity mW/cm^2
    '''
    ret_val = 0.1 * (2 *p * 1e-3)/(pi * (radi * 1e-3)**2)
    return ret_val
