from __future__ import division, print_function
from Main.Constants.Rb_constants import *
from numpy import *
from scipy.constants import k

def p(T):
    '''
    calculate the pressure in the atomic media
    :param T: temperature [k]
    :return: the pressure of atomic media
    '''
    ret_val = 10**(pa - pb/T + pc*T + pd*log10(T))
    return ret_val

def temp2den(T):
    '''
    calculate the atomic density of rubidium atoms
    :param T: temperature (K)
    :return: atomic density (1/m^3)
    '''
    ret_val = p(T)*133.3/(k * T)
    return ret_val

def temp2den_cm(T):
    '''
    calculate the atomic density of rubidium atoms
    :param T: temperature (K)
    :return: atomic density (1/cm^3)
    '''
    ret_val = temp2den(T)/(100**3)
    return ret_val

def temp2velocity(T):
    '''
    calculate the most probable velocity for a specific temperature
    :param T: temperature (K)
    :return: most probable velocity (m/s)
    '''
    ret_val = sqrt(2 * k * T / m)
    return ret_val

'''
from Global_unit_converter import *
T = celsius2kelvin(25)
print(temp2velocity(T))
'''

