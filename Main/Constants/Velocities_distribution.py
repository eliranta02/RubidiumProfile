from __future__ import division, print_function
from numpy import sqrt, exp, pi

def maxwell(v,vp):
    '''
    Maxwell-Boltzmann Velocity Distribution
    :param v: velocity [m/s]
    :param vp: most probable velocity [m/s] (you can use the temp2velocity() function to calculate the most probable velocity)
    :return: velocity distribution
    '''
    return (1/(sqrt(pi)*vp))*exp(-v**2 / vp**2)