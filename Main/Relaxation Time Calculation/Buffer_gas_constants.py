from __future__ import division,print_function
from scipy.constants import k,N_A
from numpy import sqrt, pi

class GAS:

    def set_pressure(self,pressure):
        '''

        :param pressure: set pressure in torr
        :return:
        '''
        self.__pressure = pressure

    def get_pressure(self):
        '''

        :return: pressure [Torr]
        '''
        return self.__pressure

    def set_D0(self,D0):
        self.__Do = D0

    def get_D0(self):
        return self.__Do

    def set_sigma(self, sigma):
        '''

        :param sigma: collisional cross section [cm^2]
        :return:
        '''
        self.__sigma_sd = sigma

    def get_sigma(self):
        return self.__sigma_sd

    def set_mass(self,mass):
        '''
        :param mas: mas in atomic mass unit[u] (basically taken from the periodic table)
        '''
        retVal = mass/N_A #[gram]
        retVal = 0.001 * retVal #convert to Kg
        self.__mass = retVal

    def get_mass(self):
        '''
        :return:Mass in Kg
        '''
        return self.__mass

    def set_density(self,temp,pressure=0):
        '''
        :param temp: cell temperatuer [K]
        :param pressure: buffer gas pressure [Torr]
        :return: density of the buffer gas [cm^-3]
        '''
        pressure = pressure/760 #convert Torr to atmosphere
        R = 82.06 #[cm^3 * atm * mol^-1 * K^-1]
        self.__density = pressure * N_A / (R * temp)

    def get_density(self):
        return self.__density

    def get_relative_velocity(self,T,m):
        # k = 1.380649e-23 [m^2 * kg * s^-2 * K^-1]
        M = (self.__mass * m) / (self.__mass + m)
        ret_Val = sqrt(8 * k * T/ (pi * M))
        return ret_Val



class Ar(GAS):
    def __init__(self):
        GAS.set_D0(self,0.22)
        GAS.set_sigma(self,3.7e-22)
        GAS.set_mass(self,39.948)

class N2(GAS):
    def __init__(self):
        GAS.set_D0(self,0.23)
        GAS.set_sigma(self,1e-22)
        GAS.set_mass(self,2*14.007)
