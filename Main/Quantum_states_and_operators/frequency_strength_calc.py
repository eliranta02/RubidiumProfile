from __future__ import division, print_function
from Main.Operators.Wigner_Calc import *
from numpy import sqrt, pi
from Main.Constants.Rb_constants import *
from scipy.constants import c, hbar
from q_state import State
from Main.Constants.Rb_constants import I

#calculate the Rabi frequency of define transition. the calculation is based on aritcle "On the consistency of Rabi frequency calculations"
def rabi_freq(inten,state1,state2):
    '''
    rabi frequency calculation of the transition from s1 --> s2
    :param I:
    :param state1:
    :param state2:
    :return:
    '''
    s1 = State(state1)
    f1 = s1.get_F()
    m_f1 = s1.get_mF()
    j_f1= s1.get_gF()
    l_f1 = s1.get_L()

    s2 = State(state2)
    f2 = s2.get_F()
    m_f2 = s2.get_mF()
    j_f2 = s2.get_gF()
    l_f2 = s2.get_L()

    q = m_f2 - m_f1

    a1 = -1**(1/2 * q * (q+1) + 1 + 2 * f1 - m_f1 + I + j_f2 + S + l_f2 + j_f1)
    a2 = sqrt((2 * f1 + 1) * (2 * f2 + 1) * (2 * j_f1 + 1)* (2 * j_f2 + 1)* (2 * l_f1 + 1))
    a3 = Wigner3j(f1, 1, f2, -m_f1, q, m_f2 ) * Wigner6j(j_f1, f1, I, f2, j_f2, 1) * Wigner6j(l_f1, j_f1, S, j_f2, l_f2, 1)
    ret_val = a1 * a2 * a3 * sqrt(3 * wavelen**3 * gamma * inten / (4 * pi**2 * c * hbar)) * sqrt()

    return ret_val


