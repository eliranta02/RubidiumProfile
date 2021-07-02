from __future__ import division
from numpy import pi, sqrt
from FundamentalConstants import *

'''
"""Constants relating to the rubidium-85 atom"""
I  = 2.5         #Nuclear spin
As = 1011.910813 #Ground state hyperfine constant in units of MHz
gI = -0.00029364 #nuclear spin g-factor
mass = 84.911789732*amu
'''

"""Constants relating to the rubidium-87 atom"""
I  = 3/2         #Nuclear spin
As =  3417.341305452145 #Ground state hyperfine constant in units of MHz
gI = -0.0009951414 #nuclear spin g-factor
mass = 86.909180520*amu


"""Constants relating to the rubidium D1 transition"""
wavelength=794.978969380e-9 #The weighted linecentre of the rubidium D1 line in m
wavevectorMagnitude=2.0*pi/wavelength #Magnitude of the wavevector
NatGamma=5.746 #Rubidium D1 natural linewidth in MHz
dipoleStrength=3.0*sqrt(e0*hbar*(2.0*NatGamma*(10.0**6))*(wavelength**3)/(8.0*pi))
v0=377107407.299e6 #The weighted linecentre of the rubidium D1 line in Hz