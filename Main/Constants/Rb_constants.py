from __future__ import division, print_function
from scipy.constants import atomic_mass, c
#from Global_constants import S
from numpy import pi

S = 1/2

D_line = 'D2'

# pressure atomic vapor constants
pa = -94.04826
pb = 1961.258
pc = -0.03771687
pd = 42.57526

#mass of the rubidium
#Rubidium 85
mRb85 = 85 * atomic_mass
I_85 = 5/2

#Rubidium 87
mRb85 = 87 * atomic_mass
I_87 = 3/2

#mass uses for the entire simulation
m = mRb85
I = I_85
#ground state angular momentum
L = 0
jg = L + S


# wavelength and wavevector
# D2 line

if D_line == 'D2':
    #excited state angular momentum
    L = 1
    wavelen = 780e-9
    gamma = 2 * pi * 6.06e6 # [MHz]
    je = L + S
elif D_line == 'D1':
    #excited state angular momentum
    L = 1
    wavelen = 795e-9
    gamma = 6  # [MHz]
    je = L - S

k_num = 2 * pi / wavelen
frequency = k_num/ c
