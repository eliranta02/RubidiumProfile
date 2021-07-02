from __future__ import division, print_function
from numpy import pi,sqrt,log10
from scipy.constants import k,N_A
from Buffer_gas_constants import *

pa = -94.04826
pb =1961.258
pc = -0.03771687
pd =42.57526

def p(T):
    ret_val = 10**(pa - pb/T + pc*T + pd*log10(T))
    return ret_val

def celsius2kelvin(Tc):
    ret_val = Tc + 273.15
    return ret_val

def n(T):
    ret_val = p(T)*133.3/(k * T)
    return ret_val

def n_cm(T):
    ret_val = n(T)/(100**3)
    return ret_val

def n_zeltser(T):
    A = 4.312
    B = 4040
    ret_val = (1/T) * 10**(21.866 + A -B/T)
    return ret_val * (1/((0.01)**3))

def wall_relaxation(T,P,cell_shape,param,D0):
    '''
    :param T: cell temperature [K]
    :param P: pressure [Torr]
    :param cell_shape: 'rect' / 'circ' / 'cylind'
    :param param: param of the cell dimensions in cm
    :return:
    '''
    k = 0
    if cell_shape == 'rect':
        k = pi**2 * (1/(param[0])**2 + 1/(param[1])**2 + 1/(param[2])**2)
    if cell_shape == 'circ':
        # param[0] cell radius
        k = (pi / param[0]) **2
    if cell_shape == 'cylind':
        #param[0] cell length
        #param[1] cell radius
        k = (pi/param[0])**2 + (2.405/param[1])**2

    P0 = 760 # [Torr]
    T0 = celsius2kelvin(100) #[C]

    D = D0 * (P0/P) * (T/T0)**(3/2)
    ret_val = k * D
    return ret_val

def spin_exchange_relaxation(rubidium_num,T):
    if rubidium_num == 85:
        I = 5/2
    elif rubidium_num == 87:
        I = 3/2

    qse = (3 * (2 * I + 1)**2)/(2 * I * (2 * I - 1))

    sigma_se = 1.9e-14 #[cm^2]
    sigma_se = sigma_se*(0.01*0.01) #[m^2]
    mass = 86.909184
    if rubidium_num == 85:
        mass =  84.911
    elif rubidium_num == 87:
        mass =  86.909184
    m = mass/N_A #[gram]
    m = 0.001 * m #convert to Kg
    M = m/2
    den = n(T)#n_zeltser(T)
    vrel = sqrt(8 * k * T/ (pi * M))
    ret_val = den * sigma_se * vrel
    return (1/qse)*ret_val

def spin_destruction_relaxation(rubidium_num,T,buffer_gas_list):
    P = 0.5
    if rubidium_num == 85:
        q = (38 + 52 * P ** 2 + 6 * P**4) / (3 + 10 * P ** 2 + 3 * P**4)
    elif rubidium_num == 87:
        q = (6 + 2 * P**2)/(1 + P**2)

    sigma_d = 1.6e-17 #[cm^2]
    sigma_d = sigma_d * (0.01 * 0.01) #[m^2]
    mass = 86.909184
    if rubidium_num == 85:
        mass = 85.4678
    elif rubidium_num == 87:
        mass = 86.909184
    m = mass / N_A  # [gram]
    m = 0.001 * m  # -> convert to Kg
    M = m / 2
    den = n(T)
    vrel = sqrt(8 * k * T / (pi * M))
    rb_sd = den * sigma_d * vrel
    bg_sd = 0
    for bg in buffer_gas_list:
        sigma_bg_d = bg.get_sigma()
        vrel_bg = bg.get_relative_velocity(T,m) * 100 #
        den_bg = bg.get_density() #cm^-3
        bg_sd += den_bg * sigma_bg_d * vrel_bg

    return (1/q) * (rb_sd + bg_sd)

