from __future__ import division, print_function
from Main.Operators.Wigner_Calc import *
from numpy import sqrt, pi
from Main.Constants.Rb_constants import *
from scipy.constants import c, hbar, h
from q_state import State
from Main.Constants.Rb_constants import I
#from Main.Operators.utils import intensity_calc

def i_sat(wavelen, gamma):
    '''
    calculate the saturation intensity
    :param wavelen: wave length [m]
    :param gamma: relaxation rate [1/sec]
    :return: saturation intensity [mW/cm^2]
    '''
    tau_lt = 1/gamma
    ret_val = (pi * h * c)/(3 * wavelen**3 * tau_lt)
    return 0.1 * ret_val

#calculate the Rabi frequency of define transition. the calculation is based on aritcle "On the consistency of Rabi frequency calculations"
def rabi_freq(inten,state1,state2):
    '''
    rabi frequency calculation of the transition from s1 --> s2
    :param I:
    :param state1:
    :param state2:
    :return: rabi frequency [Hz]
    '''
    s1 = state1
    f1 = s1.get_F()
    m_f1 = s1.get_mF()
    j_f1= s1.get_J()
    l_f1 = s1.get_L()

    s2 = state2
    f2 = s2.get_F()
    m_f2 = s2.get_mF()
    j_f2 = s2.get_J()
    l_f2 = s2.get_L()

    q = m_f1 - m_f2

    a1 = -1**(1/2 * q * (q+1) + 1 + 2 * f1 - m_f1 + I + j_f2 + S + l_f2 + j_f1)
    a2 = sqrt((2 * f1 + 1) * (2 * f2 + 1) * (2 * j_f1 + 1)* (2 * j_f2 + 1)* (2 * l_f1 + 1))
    a3 = Wigner3j(f2, 1, f1, m_f2, q, -m_f1) \
         * Wigner6j(j_f1, j_f2, 1, f2, f1, I) *\
         Wigner6j(l_f1, l_f2, 1, j_f2, j_f1, S)


    #ret_val = a1 * a2 * a3 * sqrt(3 * wavelen**3 * gamma * inten / (4 * pi**2 * c * hbar))
    ret_val = a1 * a2 * a3 * sqrt((gamma**2 * I)/(2 * i_sat(wavelen, gamma)))

    return ret_val

'''
Example
f1 = State(2,1,3,3,0,0,True)
f1.set_S(S)
f2 = State(2,1,4,4,1,0,False)
f2.set_S(S)
inten = intensity_calc(2,2)
print('Rabi frequency : {} MHz'.format(rabi_freq(inten,f1,f2) * 1e-6))
'''